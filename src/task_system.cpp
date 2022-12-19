#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "task_system.h"
#include <TaskScheduler.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
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
    return (md_time_current() << 8) | (slot_idx & (MAX_TASKS - 1));
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
    PoolTask(uint32_t set_beg_, uint32_t set_end_, RangeTask set_func_, void* user_data_, str_t lbl_ = {}, ID id = INVALID_ID, enki::ICompletable* dependency = 0)
        : ITaskSet(set_end_-set_beg_), m_set_func(set_func_), m_user_data(user_data_), m_range_offset(set_beg_), m_set_completed(0), m_interrupt(false), m_id(id) {
        int64_t len = MIN(lbl_.len, LABEL_SIZE-1);
        m_label = {strncpy(m_buf, lbl_.ptr, len), len};
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }

    PoolTask(Task func_, void* user_data_, str_t lbl_ = {}, ID id = INVALID_ID, enki::ICompletable* dependency = 0)
        : ITaskSet(1), m_func(func_), m_user_data(user_data_), m_set_completed(0), m_interrupt(false), m_id(id) {
        int64_t len = MIN(lbl_.len, LABEL_SIZE-1);
        m_label = {strncpy(m_buf, lbl_.ptr, len), len};
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }

    virtual ~PoolTask() {}

    virtual void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        (void)threadnum;
        if (!m_interrupt) {
            if (m_set_func)
                m_set_func(m_range_offset + range.start, m_range_offset + range.end, m_user_data);
            else if (m_func)
                m_func(m_user_data);
        }
       
        uint32_t range_ext = (range.end - range.start);

        // This should be protected with a mutex or something to ensure that only a single thread does this
        uint32_t set_size = m_set_completed += range_ext;
        if (set_size == m_SetSize) {
            pool::free_slots.push(get_slot_idx(m_id));
        }
    }

    bool Running() {
        return !GetIsComplete();
    }

    RangeTask  m_set_func = nullptr;  // either of these two are executed
    Task       m_func     = nullptr;
    void*      m_user_data = nullptr;
    uint32_t   m_range_offset = 0;
    std::atomic_uint32_t m_set_completed = 0;
    std::atomic_bool m_interrupt = false;
    enki::Dependency m_dependency;
    char m_buf[LABEL_SIZE];
    str_t m_label = {};
    ID m_id = INVALID_ID;
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(Task func, void* user_data, str_t lbl = {}, ID id = INVALID_ID, enki::ICompletable* dependency = 0) :
        IPinnedTask(0), m_function(func), m_user_data(user_data), m_id(id) {
        int64_t len = MIN(lbl.len, LABEL_SIZE-1);
        m_label = {strncpy(m_buf, lbl.ptr, len), len};
        if (dependency) {
            SetDependency(m_dependency, dependency);
        }
    }
    virtual void Execute() final {
        m_function(m_user_data);
        main::free_slots.push(get_slot_idx(m_id));
    }

    Task m_function = nullptr;
    void* m_user_data = nullptr;
    enki::Dependency m_dependency;
    char m_buf[LABEL_SIZE];
    str_t m_label = {};
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

void initialize(uint32_t num_threads = 0) {
    ts.Initialize(num_threads);
    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        pool::free_slots.push(i);
        main::free_slots.push(i);
    }
}

void shutdown() { ts.WaitforAllAndShutdown(); }

void execute_queued_tasks() {
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

ID main_enqueue(str_t label, Task func, void* user_data, ID dependency) {
    using namespace main;
    uint32_t idx = free_slots.pop();

    ID id = generate_id(idx);
    MainTask* Task = &task_data[idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(Task) MainTask(func, user_data, label, id, dep_task);

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

ID pool_enqueue(str_t label, Task func, void* user_data, ID dependency) {
    using namespace pool;

    uint32_t slot_idx = free_slots.pop();

    ID id = generate_id(slot_idx);
    PoolTask* Task = &pool::task_data[slot_idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(Task) PoolTask(func, user_data, label, id, dep_task);

    if (!dep_task) {
        queued_slots.push(slot_idx);
    }

    return id;
}

ID pool_enqueue(str_t label, uint32_t range_beg, uint32_t range_end, RangeTask range_func, void* user_data, ID dependency) {
    using namespace pool;

    uint32_t slot_idx = free_slots.pop();

    ID id = generate_id(slot_idx);
    PoolTask* Task = &pool::task_data[slot_idx];
    enki::ICompletable* dep_task = get_task(dependency);
    PLACEMENT_NEW(Task) PoolTask(range_beg, range_end, range_func, user_data, label, id, dep_task);

    if (!dep_task) {
        queued_slots.push(slot_idx);
    }

    return id;
}

bool task_is_running(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? Task->Running() : false;
}

str_t task_label(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? Task->m_label : str_t{};
}

float task_fraction_complete(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? (float)Task->m_set_completed / (float)Task->m_SetSize : 0.f;
}

void task_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id && Task->Running()) {
        ts.WaitforTask(Task);
    }
}

void task_interrupt(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id) {
        Task->m_interrupt = true;
    }
}

void task_interrupt_and_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    PoolTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id && Task->Running()) {
        Task->m_interrupt = true;
        ts.WaitforTask(Task);
    }
}

};  // namespace task_system
