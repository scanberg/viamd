#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "task_system.h"
#include <TaskScheduler.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_os.h>

#include <string.h>
#include <atomic_queue.h>

// Blatantly stolen from ImGui (thanks Omar!)
struct NewDummy {};
inline void* operator new(size_t, NewDummy, void* ptr) { return ptr; }
inline void  operator delete(void*, NewDummy, void*)   {} // This is only required so we can use the symetrical new()
#define PLACEMENT_NEW(_PTR)					new(NewDummy(), _PTR)

namespace task_system {

constexpr uint32_t MAX_TASKS = 256;
constexpr uint32_t LABEL_SIZE = 64;

static inline ID generate_id(uint32_t slot_idx) {
    return (md_os_time_current() << 8) | (slot_idx & (MAX_TASKS - 1));
}

static inline uint32_t get_slot_idx(ID id) {
    return (uint32_t)(id & (MAX_TASKS - 1));
}

namespace main {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> queued_slots;
}

namespace pool {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> queued_slots;
}

class PoolTask : public enki::ITaskSet {
public:
    PoolTask() = default;
    PoolTask(uint32_t set_size_, execute_range_fn set_func_, const char* lbl_ = "", ID id = INVALID_ID, enki::ICompletable* dependency = 0)
        : ITaskSet(set_size_), m_set_func(set_func_), m_set_completed(0), m_interrupt(false), m_id(id) {
        strncpy(m_label, lbl_, LABEL_SIZE);
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }

    PoolTask(execute_fn func_, const char* lbl_ = "", ID id = INVALID_ID, enki::ICompletable* dependency = 0)
        : ITaskSet(1), m_func(func_), m_set_completed(0), m_interrupt(false), m_id(id) {
        strncpy(m_label, lbl_, LABEL_SIZE);
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }

    virtual ~PoolTask() {}

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        (void)threadnum;
        if (!m_interrupt) {
            if (m_set_func)
                m_set_func(range.start, range.end);
            else if (m_func)
                m_func();
        }
       
        uint32_t range_ext = (range.end - range.start);

        // This should be protected with a mutex or something to ensure that only a single thread does this
        if (m_set_completed + range_ext == m_SetSize) {
            pool::free_slots.push(get_slot_idx(m_id));
        }

        m_set_completed += range_ext;
    }

    bool Running() {
        return !GetIsComplete();
    }

    execute_range_fn m_set_func = nullptr;  // either of these two are executed
    execute_fn       m_func     = nullptr;
    std::atomic_uint32_t m_set_completed = 0;
    std::atomic_bool m_interrupt = false;
    enki::Dependency m_dependency;
    char m_label[LABEL_SIZE] = {};
    ID m_id = INVALID_ID;
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(execute_fn func, const char* lbl = "", ID id = INVALID_ID, enki::ICompletable* dependency = 0) :
        IPinnedTask(0), m_function(func), m_id(id) {
        strncpy(m_label, lbl, LABEL_SIZE);
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }
    virtual void Execute() final {
        m_function();
        main::free_slots.push(get_slot_idx(m_id));
    }

    execute_fn m_function = nullptr;
    enki::Dependency m_dependency;
    char m_label[LABEL_SIZE] = {};
    ID m_id = INVALID_ID;
};

namespace main {
    static MainTask task_data[MAX_TASKS];
}

namespace pool {
    static PoolTask task_data[MAX_TASKS];
}

static inline enki::ICompletable* get_task(ID id) {
    if (id != INVALID_ID) {
        uint32_t slot_idx = get_slot_idx(id);
        PoolTask* ptask = &pool::task_data[slot_idx];
        MainTask* mtask = &main::task_data[slot_idx];
        if (ptask->m_id == id)       return ptask;
        else if (mtask->m_id == id)  return mtask;
    }
    return NULL;
}

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

void execute_tasks() {
    while (!pool::queued_slots.was_empty()) {
        uint32_t idx = pool::queued_slots.pop();
        ts.AddTaskSetToPipe(&pool::task_data[idx]);
    }
    while (!main::queued_slots.was_empty()) {
        uint32_t idx = main::queued_slots.pop();
        ts.AddPinnedTask(&main::task_data[idx]);
    }
    ts.RunPinnedTasks();
}

ID main_enqueue(const char* label, execute_fn func, ID dependency) {
    using namespace main;
    uint32_t idx = free_slots.pop();

    ID id = generate_id(idx);
    MainTask* task = &task_data[idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(task) MainTask(func, label, id, dep_task);

    if (!dep_task) {
        queued_slots.push(idx);
    }
    
    return id;
}

uint32_t pool_num_threads() { return ts.GetNumTaskThreads(); }

ID* pool_running_tasks(md_allocator_i* alloc) {
    ASSERT(alloc);
    ID* arr = 0;
    for (uint32_t i = 0; i < MAX_TASKS; ++i) {
        if (pool::task_data[i].Running()) {
            md_array_push(arr, pool::task_data[i].m_id, alloc);
        }
    }
    return arr;
}

void pool_interrupt_running_tasks() {
    for (uint32_t i = 0; i < MAX_TASKS; ++i) {
        if (pool::task_data[i].Running()) {
            pool::task_data[i].m_interrupt = true;
        }
    }
}

ID pool_enqueue(const char* label, execute_fn func, ID dependency) {
    using namespace pool;

    uint32_t slot_idx = free_slots.pop();

    ID id = generate_id(slot_idx);
    PoolTask* task = &pool::task_data[slot_idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(task) PoolTask(func, label, id, dep_task);

    if (!dep_task) {
        queued_slots.push(slot_idx);   
    }

    return id;
}

ID pool_enqueue(const char* label, uint32_t range_size, execute_range_fn range_func, ID dependency) {
    using namespace pool;

    uint32_t slot_idx = free_slots.pop();

    ID id = generate_id(slot_idx);
    PoolTask* task = &pool::task_data[slot_idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(task) PoolTask(range_size, range_func, label, id, dep_task);

    if (!dep_task) {
        queued_slots.push(slot_idx);
    }

    return id;
}

bool task_is_running(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    return task->m_id == id ? task->Running() : false;
}

const char* task_label(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    return task->m_id == id ? task->m_label : "";
}

float task_fraction_complete(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    return task->m_id == id ? (float)task->m_set_completed / (float)task->m_SetSize : 0.f;
}

void task_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    if (task->m_id == id && task->Running()) {
        ts.WaitforTask(task);
    }
}

void task_interrupt(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    if (task->m_id == id) {
        task->m_interrupt = true;
    }
}

void task_interrupt_and_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* task = &pool::task_data[slot_idx];
    if (task->m_id == id && task->Running()) {
        task->m_interrupt = true;
        ts.WaitforTask(task);
    }
}

};  // namespace task_system
