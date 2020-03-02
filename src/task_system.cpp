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

    TaskSetFunction m_function = nullptr;
    std::atomic<uint32_t> m_set_completed = 0;
    std::atomic<bool> m_interrupt = false;
    char m_label[LABEL_SIZE] = {};
};

constexpr ID base_id = 0xdeadb00b;
static enki::TaskScheduler ts{};

static LIFO<ID, MAX_TASKS> free = {};
static LIFO<ID, MAX_TASKS> used = {};
static Task task_data[MAX_TASKS] {};

static Task* get_task(ID id) {
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

uint32_t get_num_threads() { return (uint32_t)ts.GetNumTaskThreads(); }

uint32_t get_num_tasks() { return (uint32_t)used.size(); }

ID* get_tasks() { return used.data; }

void clear_completed_tasks() {
    for (uint32_t i = 0; i < used.size();) {
        const ID id = used[i];
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

ID create_task(const char* label, TaskSetFunction func) {
    if (free.empty()) {
        LOG_ERROR("Task queue is full! Cannot create task");
        return 0;
    }

    const ID id = free.back();
    free.pop();
    used.push(id);

    Task* task_system = get_task(id);
    ASSERT(task_system);
    PLACEMENT_NEW(task_system) Task(1, func, label);

    ts.AddTaskSetToPipe(task_system);

    return id;
}

ID create_task(const char* label, uint32_t size, TaskSetFunction func) {
    if (free.empty()) {
        LOG_ERROR("Task queue is full! Cannot create task");
        return 0;
    }

    const ID id = free.back();
    free.pop();
    used.push(id);

    Task* task_system = get_task(id);
    ASSERT(task_system);
    PLACEMENT_NEW(task_system) Task(size, func, label);

    ts.AddTaskSetToPipe(task_system);

    return id;
}

const char* get_task_label(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return NULL;
    }

    return task_system->m_label;
}

bool get_task_complete(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return false;
    }

    return task_system->GetIsComplete();
}

float get_task_fraction_complete(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return 0.0f;
    }
    return (float)task_system->m_set_completed / (float) task_system->m_SetSize;
}

void wait_for_task(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return;
    }

    ts.WaitforTask(task_system);
}

void interrupt_task(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return;
    }

    task_system->m_interrupt = true;
}

void interrupt_and_wait(ID id) {
    Task* task_system = get_task(id);
    if (!task_system) {
        LOG_ERROR("Invalid id");
        return;
    }

    task_system->m_interrupt = true;
    ts.WaitforTask(task_system);
}

};  // namespace task_system
