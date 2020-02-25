#include "task_system.h"
#include <TaskScheduler.h>
#include <core/common.h>
#include <core/log.h>

namespace task_system {

constexpr uint32_t MAX_TASKS = 32;
constexpr uint32_t LABEL_SIZE = 64;

class Task : public enki::ITaskSet {
public:
    Task() = default;
    Task(TaskSetFunction func_, const char* lbl) : m_function(func_), m_set_completed(0), m_interrupt(false) {
        strncpy(m_label, lbl, LABEL_SIZE);
    }
    Task(uint32_t set_size_, TaskSetFunction func_, const char* lbl) : ITaskSet(set_size_), m_function(func_), m_set_completed(0), m_interrupt(false) {
        strncpy(m_label, lbl, LABEL_SIZE);
    }

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        if (!m_interrupt) {
            m_function({range.start, range.end}, {threadnum, m_interrupt});
        }
        m_set_completed += (range.end - range.start);
    }

    TaskSetFunction m_function;
    std::atomic<uint32_t> m_set_completed;
    std::atomic<bool> m_interrupt;
    char m_label[LABEL_SIZE] = {};
};


template <typename T, uint32_t N>
struct LIFO {
    T data[N] = {};
    uint32_t count = 0;

    bool push(const T& item) {
        if (count == N) return false;
        data[count] = item;
        count++;
        return true;
    }

    void pop() { if (count > 0) count--; }
    void clear() { count = 0; }
    bool empty() { return count == 0; }
    uint32_t size() { return count; }

    T& back() { ASSERT(size() > 0); return data[count - 1]; }
    const T& back() const { ASSERT(size() > 0); return data[count - 1]; }

    T* begin() { return data; }
    T* end() { return data + count; }

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }
};

constexpr TaskID base_id = 0xdeadb00b;
static enki::TaskScheduler ts{};

static LIFO<TaskID, MAX_TASKS> free = {};
static LIFO<TaskID, MAX_TASKS> used = {};
static Task task_data[MAX_TASKS] {};

static Task* get_task(TaskID id) {
    const uint32_t idx = id - base_id;
    if (idx < 0 || MAX_TASKS <= idx) return nullptr;
    return &task_data[idx];
}

void initialize() {
    ts.Initialize();

    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        free.push(base_id + i);
    }
}

void shutdown() {
    ts.WaitforAllAndShutdown();
}

uint32_t get_num_tasks() { return (uint32_t)used.size(); }
TaskID* get_tasks() { return used.data; }

void clear_completed_tasks() {
    LIFO<uint32_t, MAX_TASKS> remove = {};
    for (uint32_t i = 0; i < used.size();) {
        const TaskID id = used[i];
        const uint32_t idx = id - base_id;
        if (task_data[idx].GetIsComplete()) {
            free.push(id);
            used[i] = used.back();
            used.pop();
        } else {
            i++;
        }
    }
}

TaskID create_task(const char* label, TaskSetFunction func) {
    if (free.empty()) {
        LOG_ERROR("Task queue is full! Cannot create task");
        return 0;
    }

    const TaskID id = free.back();
    free.pop();
    used.push(id);

    Task* task = get_task(id);
    ASSERT(task);
    PLACEMENT_NEW(task) Task(1, func, label);

    ts.AddTaskSetToPipe(task);

    return id;
}

TaskID create_task(const char* label, uint32_t size, TaskSetFunction func) {
    if (free.empty()) {
        LOG_ERROR("Task queue is full! Cannot create task");
        return 0;
    }

    const TaskID id = free.back();
    free.pop();
    used.push(id);

    Task* task = get_task(id);
    ASSERT(task);
    PLACEMENT_NEW(task) Task(size, func, label);

    ts.AddTaskSetToPipe(task);

    return id;
}

const char* get_task_label(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return false;
    }

    return task->m_label;
}

bool get_task_complete(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return false;
    }

    return task->GetIsComplete();
}

float get_task_fraction_complete(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return 0.0f;
    }
    return (float)task->m_set_completed / (float) task->m_SetSize;
}

void wait_for_task(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    ts.WaitforTask(task);
}

void interrupt_task(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = 1;
}

void interrupt_and_wait(TaskID id) {
    Task* task = get_task(id);
    if (!task) {
        LOG_ERROR("Invalid id");
        return;
    }

    task->m_interrupt = 1;
    ts.WaitforTask(task);
}

};  // namespace task_system