#pragma once

#include <stdint.h>

struct md_allocator_i;

namespace task_system {

typedef uint64_t ID;
constexpr ID INVALID_ID = 0;

typedef void (*execute_fn)(void* user_data);
typedef void (*execute_range_fn)(uint32_t range_beg, uint32_t range_end, void* user_data);

void initialize();
void shutdown();

// This is to generate tasks for the main thread ("render" thread)
bool main_enqueue(const char* label, void* user_data, execute_fn task);

// Call once per frame at some approriate time
void main_execute_tasks();

// This is to generate tasks for the thread-pool (async operations)
ID pool_enqueue(const char* label, void* user_data, execute_fn task, execute_fn on_complete = NULL);
ID pool_enqueue(const char* label, void* user_data, uint32_t range_size, execute_range_fn range_task, execute_fn on_complete = NULL);

uint32_t pool_num_threads();

// These do not really reflect the 'current' state since that is illdefined. But rather what the state was at the time of the function call.
void pool_interrupt_running_tasks();
ID*  pool_running_tasks(md_allocator_i* alloc);

// These are safe to call with an invalid id, in such case, they will just return some 'zero' default value
bool        task_is_running(ID);
const char* task_label(ID);
float       task_fraction_complete(ID);

// These are safe to call with an invalid id, and in such case, they do nothing
void task_wait_for(ID);
void task_interrupt(ID);
void task_interrupt_and_wait_for(ID);

}  // namespace task_system