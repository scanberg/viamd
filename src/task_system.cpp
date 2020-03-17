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
        if (count == N) {
            return false;
        }
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

class PoolTask : public enki::ITaskSet {
public:
    PoolTask() = default;
    PoolTask(uint32_t set_size_, TaskSet func_, const char* lbl_ = "")
        : ITaskSet(set_size_), m_set_func(func_), m_set_completed(0), m_interrupt(false) {
        strncpy(m_label, lbl_, LABEL_SIZE);
    }

    virtual ~PoolTask() {}

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        (void)threadnum;
        if (!m_interrupt) {
            m_set_func({range.start, range.end});
        }
        m_set_completed += (range.end - range.start);
    }

    TaskSet m_set_func = nullptr;
    enki::Dependency m_dependency = {};
    std::atomic_uint32_t m_set_completed = 0;
    std::atomic_bool m_interrupt = false;
    char m_label[LABEL_SIZE] = {};
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(Task func, const char* lbl) : m_function(func) { strncpy(m_label, lbl, LABEL_SIZE); }
    virtual void Execute() final {
        m_function();
    }
    Task m_function = nullptr;
    char m_label[LABEL_SIZE] = {};
};

constexpr ID base_id = {0xdeadb00b};
static enki::TaskScheduler ts{};

namespace main {
static std::mutex mutex;
static MainTask task_data[MAX_TASKS];
static std::atomic_uint32_t task_count = 0;
}

namespace pool {
static std::mutex mutex;
static LIFO<ID, MAX_TASKS> free {};
static LIFO<ID, MAX_TASKS> used {};
static PoolTask task_data[MAX_TASKS];

}
void clear_completed_tasks();

void initialize() {
    ts.Initialize();
    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        pool::free.enqueue({base_id.id + i});
    }
}
void shutdown() { ts.WaitforAllAndShutdown(); }

bool enqueue_main(const char* label, Task func) {
    using namespace main;
    uint32_t idx = task_count++;
    if (idx >= MAX_TASKS) return false;

    mutex.lock();
    MainTask* task = &task_data[idx];
    PLACEMENT_NEW(task) MainTask(func, label);
    ts.AddPinnedTask(task);
    mutex.unlock();
    return true;
}

void run_main_tasks() {
    using namespace main;
    mutex.lock();
    ts.RunPinnedTasks();
    task_count = 0;
    mutex.unlock();
    clear_completed_tasks();
}

static PoolTask* get_task(ID id) {
    using namespace pool;
    const uint32_t idx = id.id - base_id.id;
    if (idx < 0 || MAX_TASKS <= idx) return nullptr;
    return &task_data[idx];
}

uint32_t get_num_threads() { return (uint32_t)ts.GetNumTaskThreads(); }
uint32_t get_num_tasks() { return (uint32_t)pool::used.size(); }
ID* get_tasks() { return pool::used.data; }

void clear_completed_tasks() {
    using namespace pool;
    mutex.lock();
    for (uint32_t i = 0; i < used.size();) {
        const ID id = used[i];
        const uint32_t idx = id.id - base_id.id;
        if (task_data[idx].GetIsComplete()) {
            free.enqueue(id);
            used.dequeue(used[i]);
        } else {
            i++;
        }
    }
    mutex.unlock();
}

ID enqueue_pool(const char* label, uint32_t size, TaskSet func, TaskSet on_complete) {
    using namespace pool;
    const uint32_t num_tasks = (on_complete != nullptr) ? 2 : 1;

    mutex.lock();
    if (free.size() < num_tasks) {
        LOG_ERROR("Task queue is full! Cannot create task");
        mutex.unlock();
        return INVALID_ID;
    }

    ID set_id = INVALID_ID;
    free.dequeue(set_id);
    used.enqueue(set_id);

    ID com_id = INVALID_ID;
    if (num_tasks > 1) {
        free.dequeue(com_id);
        used.enqueue(com_id);
    }
    mutex.unlock();

    PoolTask* set_task = get_task(set_id);
    ASSERT(set_task);
    PLACEMENT_NEW(set_task) PoolTask(size, func, label);

    if (num_tasks > 1) {
        PoolTask* com_task = get_task(com_id);
        ASSERT(com_task);
        PLACEMENT_NEW(com_task) PoolTask(1, on_complete, "");
        set_task->m_dependency.SetDependency(set_task, com_task);
    }

    ts.AddTaskSetToPipe(set_task);
    return set_id;
}

const char* get_task_label(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return NULL;
    }

    return task->m_label;
}

bool get_task_complete(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return false;
    }

    return task->GetIsComplete();
}

float get_task_fraction_complete(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return 0.0f;
    }
    return (float)task->m_set_completed / (float)task->m_SetSize;
}

void wait_for_task(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    ts.WaitforTask(task);
}

void interrupt_task(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = true;
}

void interrupt_and_wait(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = true;
    ts.WaitforTask(task);
}

};  // namespace task_system
