#include "task_system.h"
#include <TaskScheduler.h>
#include <core/common.h>
#include <core/log.h>
#include <core/sync.h>

template <typename T, uint32_t N>
struct LIFO {
    T data[N] = {};
    uint32_t count = 0;

    bool enqueue(const T& item) {
        if (count == N) return false;
        data[count] = item;
        count++;
        return true;
    }

    bool dequeue(T& item) {
        if (count == 0) return false;
        item = data[count - 1];
        count--;
        return true;
    }

    bool pop() {
        if (count == 0) return false;
        --count;
        return true;
    }

    const T& back() const { return data[count - 1]; }
    T& back() { return data[count - 1]; }

    void clear() { count = 0; }
    bool empty() const { return count == 0; }
    uint32_t size() const { return count; }

    const T* begin() const { return data; }
    const T* end() const { return data + count; }

    T* begin() { return data; }
    T* end() { return data + count; }

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }

    operator Array<T>() const { return {data, count}; }
};

namespace task_system {

constexpr uint32_t MAX_TASKS = 256;
constexpr uint32_t LABEL_SIZE = 64;

class Task : public enki::ITaskSet {
public:
    Task() = default;
    Task(TaskSetFunction func_, const char* lbl) : m_function(func_), m_set_completed(0), m_interrupt(false) { strncpy(m_label, lbl, LABEL_SIZE); }
    Task(uint32_t set_size_, TaskSetFunction func_, const char* lbl)
        : ITaskSet(set_size_), m_function(func_), m_set_completed(0), m_interrupt(false) {
        strncpy(m_label, lbl, LABEL_SIZE);
    }

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        if (!m_interrupt) {
            m_function({range.start, range.end}, {threadnum, m_interrupt});
        }
        m_set_completed += (range.end - range.start);
    }

    TaskSetFunction m_function = nullptr;
    enki::Dependency m_dependency = {};
    std::atomic<uint32_t> m_set_completed = 0;
    std::atomic<bool> m_interrupt = false;
    char m_label[LABEL_SIZE] = {};
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(MainFunction func, const char* lbl) : m_function(func) { strncpy(m_label, lbl, LABEL_SIZE); }
    virtual void Execute() final {
        m_function();
    }
    MainFunction m_function = nullptr;
    char m_label[LABEL_SIZE] = {};
};

constexpr ID base_id = {0xdeadb00b};
static enki::TaskScheduler ts{};

namespace main {
static std::mutex mutex;
static MainTask task_data[MAX_TASKS]{};
static std::atomic_uint32_t task_count = 0;
}

namespace pool {
static std::mutex mutex;
static LIFO<ID, MAX_TASKS> free {};
static LIFO<ID, MAX_TASKS> used {};
static Task task_data[MAX_TASKS] {};
void clear_completed_tasks();
}

void initialize() {
    ts.Initialize();
    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        pool::free.enqueue({base_id.id + i});
    }
}
void shutdown() { ts.WaitforAllAndShutdown(); }

namespace main {

bool enqueue(const char* label, MainFunction func) {
    uint32_t idx = task_count++;
    if (idx >= MAX_TASKS) return false;

    mutex.lock();
    MainTask* task = &task_data[idx];
    PLACEMENT_NEW(task) MainTask(func, label);
    ts.AddPinnedTask(task);
    mutex.unlock();
}

void run_tasks() {
    mutex.lock();
    ts.RunPinnedTasks();
    task_count = 0;
    mutex.unlock();
    pool::clear_completed_tasks();
}

}

namespace pool {

static Task* get_task(ID id) {
    const uint32_t idx = id.id - base_id.id;
    if (idx < 0 || MAX_TASKS <= idx) return nullptr;
    return &task_data[idx];
}

uint32_t get_num_threads() { return (uint32_t)ts.GetNumTaskThreads(); }
uint32_t get_num_tasks() { return (uint32_t)used.size(); }
ID* get_tasks() { return used.data; }

void clear_completed_tasks() {
    mutex.lock();
    for (uint32_t i = 0; i < used.size();) {
        const ID id = used[i];
        const uint32_t idx = id.id - base_id.id;
        if (task_data[idx].GetIsComplete()) {
            free.enqueue(id);
            used[i] = used.back();
            used.pop();
        } else {
            i++;
        }
    }
    mutex.unlock();
}

ID enqueue(const char* label, uint32_t size, TaskSetFunction func) {
    if (free.empty()) {
        LOG_ERROR("Task queue is full! Cannot create task");
        return INVALID_ID;
    }

    mutex.lock();
    const ID id = free.back();
    free.pop();
    used.enqueue(id);
    mutex.unlock();

    Task* task = get_task(id);
    ASSERT(task);
    PLACEMENT_NEW(task) Task(size, func, label);

    ts.AddTaskSetToPipe(task);

    return id;
}

const char* get_task_label(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return NULL;
    }

    return task->m_label;
}

bool get_task_complete(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return false;
    }

    return task->GetIsComplete();
}

float get_task_fraction_complete(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return 0.0f;
    }
    return (float)task->m_set_completed / (float)task->m_SetSize;
}

void wait_for_task(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    ts.WaitforTask(task);
}

void interrupt_task(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = true;
}

void interrupt_and_wait(ID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = true;
    ts.WaitforTask(task);
}
}


};  // namespace task_system
