#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>
#include <functional>

namespace task_system {

typedef uint64_t ID;
constexpr ID INVALID_ID = 0;

using Task = std::function<void()>;
using RangeTask = std::function<void(uint32_t range_beg, uint32_t range_end, uint32_t thread_num)>;

/*
using Task      = void (*)(void* user_data);
using RangeTask = void (*)(uint32_t range_beg, uint32_t range_end, void* user_data, uint32_t thread_num);
*/

void initialize(size_t num_threads);
void shutdown();

// Call once per frame at some approriate time, if there are items in the main queue, the main thread will be stalled.
// Pool tasks will not stall the main thread.
//void execute_queued_tasks();

ID create_main_task(str_t label, Task task);
ID create_pool_task(str_t label, Task task);
ID create_pool_task(str_t label, uint32_t range_size, RangeTask task, uint32_t grain_size = 1);

// Sets a dependency for a task such that the task will only be executed upon the completion of 'dependency'
void set_task_dependency(ID task, ID dependency);

// Once the task is created and potential dependencies has been declared, the task is enqueued through this procedure.
void enqueue_task(ID task);

// Call once per frame at the appropriate time to execute items queued up for the main thread. The main thread will be stalled and the tasks will be performed
void execute_main_task_queue();

// This signals interruption for all running tasks
void pool_interrupt_running_tasks();

// This halts the calling thread until all running tasks have completed
void pool_wait_for_completion();

size_t pool_num_threads();

// These do not really reflect the 'current' state since that is illdefined. But rather what the state was at the time of the function call.
size_t pool_running_tasks(ID* out_id_arr, size_t id_arr_cap);

// These are safe to call with an invalid id, in such case, they will just return some default value
bool  task_is_running(ID);
str_t task_label(ID);
float task_fraction_complete(ID);

// These are safe to call with an invalid id, and in such case, they do nothing
void task_wait_for(ID);
void task_interrupt(ID);
void task_interrupt_and_wait_for(ID);


}  // namespace task_system