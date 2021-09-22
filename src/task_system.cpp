#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "task_system.h"
#include <TaskScheduler.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <string.h>

// Blatantly stolen from ImGui (thanks Omar!)
struct NewDummy {};
inline void* operator new(size_t, NewDummy, void* ptr) { return ptr; }
inline void  operator delete(void*, NewDummy, void*)   {} // This is only required so we can use the symetrical new()
#define PLACEMENT_NEW(_PTR)					new(NewDummy(), _PTR)

namespace task_system {

constexpr uint32_t MAX_TASKS = 256;
constexpr uint32_t LABEL_SIZE = 64;

namespace main {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> used_slots;
}

class PoolTask;

namespace pool {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
    uint32_t get_slot_idx(PoolTask* task);
}

class PoolTask : public enki::ITaskSet {
public:
    PoolTask() = default;
    PoolTask(uint32_t set_size_, TaskSet func_, const char* lbl_ = "")
        : ITaskSet(set_size_), m_set_func(func_), m_set_completed(0), m_interrupt(false), m_running(true) {
        strncpy(m_label, lbl_, LABEL_SIZE);
    }

    virtual ~PoolTask() {}

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        (void)threadnum;
        if (!m_interrupt) {
            m_set_func({range.start, range.end});
        }
        m_set_completed += (range.end - range.start);
        if (m_set_completed == m_SetSize) {
            uint32_t slot_idx = pool::get_slot_idx(this);
            ASSERT(0 <= slot_idx && slot_idx < MAX_TASKS);
            m_running = false;
            pool::free_slots.push(slot_idx);
        }
    }

    bool Running() {
        return m_running;
    }

    TaskSet m_set_func = nullptr;
    enki::Dependency m_dependency = {};
    std::atomic_uint32_t m_set_completed = 0;
    std::atomic_bool m_interrupt = false;
    char m_label[LABEL_SIZE] = {};
    bool m_running = false;
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(Task func, const char* lbl) : m_function(func), m_running(true) { strncpy(m_label, lbl, LABEL_SIZE); }
    virtual void Execute() final {
        m_function();
        m_running = false;
    }

    bool Running() {
        return m_running;
    }

    Task m_function = nullptr;
    char m_label[LABEL_SIZE] = {};
    bool m_running = false;
};

namespace main {
    static MainTask task_data[MAX_TASKS];
}

namespace pool {
    static PoolTask task_data[MAX_TASKS];

    // These are not kept in exact sync, as that would be very costly.
    // This list is filled in upon request
    static uint32_t num_running_tasks;
    static ID running_tasks[MAX_TASKS];

    uint32_t get_slot_idx(PoolTask* task) {
        return (uint32_t)(task - task_data);
    }
}

constexpr ID base_id = {0xdeadb00b};
static enki::TaskScheduler ts{};

void initialize() {
    ts.Initialize();
    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        pool::free_slots.push(i);
    }

    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        main::free_slots.push(i);
    }
}

void shutdown() { ts.WaitforAllAndShutdown(); }

bool enqueue_main(const char* label, Task func) {
    using namespace main;
    uint32_t idx = main::free_slots.pop();
    main::used_slots.push(idx);
    MainTask* task = &task_data[idx];
    PLACEMENT_NEW(task) MainTask(func, label);
    ts.AddPinnedTask(task);
    return true;
}

void run_main_tasks() {
    ts.RunPinnedTasks();
    while (!main::used_slots.was_empty()) {
        uint32_t idx = main::used_slots.pop();
        main::free_slots.push(idx);
    }
}

static inline PoolTask* get_task(ID id) {
    using namespace pool;
    const uint32_t idx = id.id - base_id.id;
    if (idx < 0 || MAX_TASKS <= idx) return nullptr;
    return &task_data[idx];
}

uint32_t get_num_threads() { return (uint32_t)ts.GetNumTaskThreads(); }
uint32_t get_num_tasks() { return (uint32_t)(pool::free_slots.capacity() - pool::free_slots.was_size()); }

uint32_t get_tasks(ID** ids) {
    pool::num_running_tasks = 0;
    for (uint32_t i = 0; i < MAX_TASKS; ++i) {
        if (pool::task_data[i].Running()) {
            ID id = {base_id.id + i};
            pool::running_tasks[pool::num_running_tasks++] = id;
        }
    }
    *ids = pool::running_tasks;
    return pool::num_running_tasks;
}

ID enqueue_pool(const char* label, uint32_t size, TaskSet func, TaskSet on_complete) {
    using namespace pool;
    const uint32_t num_tasks = (on_complete != nullptr) ? 2 : 1;

    uint32_t set_idx = free_slots.pop();
    ID set_id = {base_id.id + set_idx};

    ID com_id = INVALID_ID;
    if (num_tasks > 1) {
        uint32_t com_idx = free_slots.pop();
        com_id = {base_id.id + com_idx};
    }

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
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return NULL;
    }

    return task->m_label;
}

bool get_task_complete(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return false;
    }

    return task->GetIsComplete();
}

float get_task_fraction_complete(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return 0.0f;
    }
    return (float)task->m_set_completed / (float)task->m_SetSize;
}

void wait_for_task(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return;
    }

    ts.WaitforTask(task);
}

void interrupt_task(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return;
    }

    task->m_interrupt = true;
}

void interrupt_and_wait(ID id) {
    using namespace pool;

    PoolTask* task = get_task(id);
    if (!task) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid id");
        return;
    }

    task->m_interrupt = true;
    ts.WaitforTask(task);
}

};  // namespace task_system
