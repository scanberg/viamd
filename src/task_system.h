#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <functional>

struct md_allocator_i;

namespace task_system {

typedef uint64_t ID;
constexpr ID INVALID_ID = 0;

using Task = std::function<void()>;
using RangeTask = std::function<void(uint32_t range_beg, uint32_t range_end)>;

//using Task = void (*)(void* user_data);
//using RangeTask = void (*)(uint32_t range_beg, uint32_t range_end, void* user_data);

/*
typedef void (*Task) (void* user_data);
typedef void (*RangeTask)(uint32_t range_beg, uint32_t range_end, void *user_data);
*/

void initialize(uint32_t num_threads);
void shutdown();

// Call once per frame at some approriate time, if there are items in the main queue, the main thread will be stalled.
// Pool tasks will not stall the main thread.
void execute_queued_tasks();

// Execute the task immediately.
// If the task is queued for the main thread, it will stall the thread and wait for completion.
// If the task us queued for the thread-pool, it will continue execution immediately.
void execute_task(ID);

// This is to generate tasks for the main thread ("render" thread)
ID main_enqueue(str_t label, Task task, ID dependency = 0);

// This is to generate tasks for the thread-pool (async operations)
ID pool_enqueue(str_t label, Task task, ID dependency = 0);
ID pool_enqueue(str_t label, uint32_t range_beg, uint32_t range_end, RangeTask task, ID dependency = 0);

uint32_t pool_num_threads();

// These do not really reflect the 'current' state since that is illdefined. But rather what the state was at the time of the function call.
void pool_interrupt_running_tasks();
ID*  pool_running_tasks(md_allocator_i* alloc);

// These are safe to call with an invalid id, in such case, they will just return some 'zero' default value
bool  task_is_running(ID);
str_t task_label(ID);
float task_fraction_complete(ID);

// These are safe to call with an invalid id, and in such case, they do nothing
void task_wait_for(ID);
void task_interrupt(ID);
void task_interrupt_and_wait_for(ID);

/*
ID task_create(str_t label, Task Task);
ID task_create(str_t label, uint32_t range_size, RangeTask RangeTask);
ID task_create(str_t label, uint32_t range_beg, uint32_t range_end, RangeTask RangeTask);

void task_set_dependency(ID, ID dependency);

void pool_enqueue(ID);
void main_enqueue(ID);
*/

}  // namespace task_system