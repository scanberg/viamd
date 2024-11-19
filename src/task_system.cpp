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

#define MAX_TASKS 256
#define LABEL_SIZE 64

static inline ID generate_id(uint32_t slot_idx) {
    return (md_time_current() << 8) | (slot_idx & (MAX_TASKS - 1));
}

static inline uint32_t get_slot_idx(ID id) {
    return (uint32_t)(id & (MAX_TASKS - 1));
}

namespace main {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
}

namespace pool {
    static atomic_queue::AtomicQueue<uint32_t, MAX_TASKS, 0xFFFFFFFF> free_slots;
}

struct CompletionActionFreePoolSlot : public enki::ICompletable {
public:
    CompletionActionFreePoolSlot() = default;
    void OnDependenciesComplete(enki::TaskScheduler* taskscheduler, uint32_t threadnum ) final {
        ICompletable::OnDependenciesComplete(taskscheduler, threadnum);
        pool::free_slots.push(m_slot_idx);
    }
    uint32_t m_slot_idx = UINT32_MAX;
    enki::Dependency m_dependency = {};
};

class AsyncTask : public enki::ITaskSet {
public:
    AsyncTask() = default;
    AsyncTask(uint32_t set_size, RangeTask set_func, str_t lbl = {}, ID id = INVALID_ID, uint32_t grain_size = 1)
        : ITaskSet(set_size, grain_size), m_set_func(set_func), m_set_complete(0), m_interrupt(false), m_dependency(), m_id(id), m_grain_size(grain_size) {
        size_t len = str_copy_to_char_buf(m_buf, sizeof(m_buf), lbl);
        m_label = {m_buf, len};
        m_completion_action.m_slot_idx = get_slot_idx(id);
        m_completion_action.SetDependency(m_completion_action.m_dependency, this);
    }

    AsyncTask(Task func, str_t lbl = {}, ID id = INVALID_ID)
        : ITaskSet(1), m_func(func), m_set_complete(0), m_interrupt(false), m_dependency(), m_id(id) {
        size_t len = str_copy_to_char_buf(m_buf, sizeof(m_buf), lbl);
        m_label = {m_buf, len};
        m_completion_action.m_slot_idx = get_slot_idx(id);
        m_completion_action.SetDependency(m_completion_action.m_dependency, this);
    }

    void ExecuteRange(enki::TaskSetPartition range, uint32_t threadnum) final {
        (void)threadnum;
        if (!m_interrupt) {
            if (m_set_func) {
                m_set_func(range.start, range.end, threadnum);
                m_set_complete += (range.end - range.start);
            }
            else if (m_func) {
                m_func();
                m_set_complete += 1;
            }
        }
    }

    inline bool Running() const {
        return !GetIsComplete();
    }

    RangeTask  m_set_func = nullptr;  // either of these two are executed
    Task       m_func     = nullptr;
    uint32_t   m_range_offset = 0;
    uint32_t   m_grain_size = 1;
    std::atomic_uint32_t m_set_complete = 0;
    std::atomic_bool m_interrupt = false;
    enki::Dependency m_dependency = {};
    CompletionActionFreePoolSlot m_completion_action = {};
    char m_buf[LABEL_SIZE] = "";
    str_t m_label = {};
    ID m_id = INVALID_ID;
};

class MainTask : public enki::IPinnedTask {
public:
    MainTask() = default;
    MainTask(Task func, str_t lbl = {}, ID id = INVALID_ID) :
        IPinnedTask(0), m_function(func), m_id(id) {
        size_t len = MIN(lbl.len, LABEL_SIZE-1);
        m_label = {strncpy(m_buf, lbl.ptr, len), len};
    }

    void Execute() final {
        m_function();
        main::free_slots.push(get_slot_idx(m_id));
    }

    Task m_function = nullptr;
    enki::Dependency m_dependency = {};
    char m_buf[LABEL_SIZE] = "";
    str_t m_label = {};
    ID m_id = INVALID_ID;
};

namespace main {
    static MainTask task_data[MAX_TASKS];
}

namespace pool {
    static AsyncTask task_data[MAX_TASKS];
}

static inline enki::ICompletable* get_task(ID id) {
    if (id != INVALID_ID) {
        uint32_t slot_idx = get_slot_idx(id);
        AsyncTask* ptask = &pool::task_data[slot_idx];
        MainTask*  mtask = &main::task_data[slot_idx];
        if (ptask->m_id == id)       return ptask;
        else if (mtask->m_id == id)  return mtask;
    }
    return NULL;
}

static enki::TaskScheduler ts{};

void initialize(size_t num_threads = 0) {
    ts.Initialize((uint32_t)num_threads);
    for (uint32_t i = 0; i < MAX_TASKS; i++) {
        pool::free_slots.push(i);
        main::free_slots.push(i);
    }
}

void shutdown() { ts.WaitforAllAndShutdown(); }

ID create_main_task(str_t label, Task func) {
    const uint32_t idx = main::free_slots.pop();
    ID id = generate_id(idx);
    MainTask* task = &main::task_data[idx];
    PLACEMENT_NEW(task) MainTask(func, label, id);
    return id;
}

ID create_pool_task(str_t label, Task func) {
    const uint32_t idx = pool::free_slots.pop();
    ID id = generate_id(idx);
    AsyncTask* task = &pool::task_data[idx];
    PLACEMENT_NEW(task) AsyncTask(func, label, id);
    return id;
}

ID create_pool_task(str_t label, uint32_t range_size, RangeTask func) {
    const uint32_t idx = pool::free_slots.pop();
    ID id = generate_id(idx);
    AsyncTask* task = &pool::task_data[idx];
    PLACEMENT_NEW(task) AsyncTask(range_size, func, label, id);
    return id;
}

void enqueue_task(ID id) {
    const uint32_t slot_idx = get_slot_idx(id);
    {
        AsyncTask* task = &pool::task_data[slot_idx];
        if (task->m_id == id) {
            if (task->m_dependency.GetDependencyTask() != NULL) goto dep_error;
            ts.AddTaskSetToPipe(&pool::task_data[slot_idx]);
            return;
        }
    }
    {
        MainTask* task = &main::task_data[slot_idx];
        if (task->m_id == id) {
            if (task->m_dependency.GetDependencyTask() != NULL) goto dep_error;
            ts.AddPinnedTask(task);
            return;
        }
    }
    MD_LOG_DEBUG("Invalid Operation: Attempting to enquque invalid task id");
    ASSERT(false);
    return;
dep_error:
    MD_LOG_DEBUG("Invalid Operation: Attempting to enquque task which has dependency set.");
    ASSERT(false);
    return;
}

void execute_main_task_queue() {
    ts.RunPinnedTasks();
}

size_t pool_num_threads() { return ts.GetNumTaskThreads(); }

size_t pool_running_tasks(ID* out_id_arr, size_t id_arr_cap) {
    size_t num_tasks = 0;
    for (size_t i = 0; i < MAX_TASKS; ++i) {
        if (pool::task_data[i].Running()) {
            out_id_arr[num_tasks++] = pool::task_data[i].m_id;
            if (num_tasks == id_arr_cap) break;
        }
    }
    return num_tasks;
}

void pool_interrupt_running_tasks() {
    for (uint32_t i = 0; i < MAX_TASKS; ++i) {
        if (pool::task_data[i].Running()) {
            pool::task_data[i].m_interrupt = true;
        }
    }
}

void pool_wait_for_completion() {
    ts.WaitforAll();
}

void set_task_dependency(ID task_id, ID dep_id) {
    enki::ICompletable* dep = get_task(dep_id);
    if (dep == NULL) return;

    uint32_t task_idx = get_slot_idx(task_id);
    if (pool::task_data[task_idx].m_id == task_id) {
        pool::task_data[task_idx].SetDependency(pool::task_data[task_idx].m_dependency, dep);
        return;
    }
    if (main::task_data[task_idx].m_id == task_id) {
        main::task_data[task_idx].SetDependency(main::task_data[task_idx].m_dependency, dep);
        return;
    }
}

bool task_is_running(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? Task->Running() : false;
}

str_t task_label(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? Task->m_label : str_t{};
}

float task_fraction_complete(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    return Task->m_id == id ? (float)Task->m_set_complete / (float)Task->m_SetSize : 0.f;
}

void task_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id && Task->Running()) {
        ts.WaitforTask(Task);
    }
}

void task_interrupt(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id) {
        Task->m_interrupt = true;
    }
}

void task_interrupt_and_wait_for(ID id) {
    uint32_t slot_idx = get_slot_idx(id);
    AsyncTask* Task = &pool::task_data[slot_idx];
    if (Task->m_id == id && Task->Running()) {
        Task->m_interrupt = true;
        ts.WaitforTask(Task);
    }
}

};  // namespace task_system
