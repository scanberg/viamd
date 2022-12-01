#pragma once

#include <stdint.h>
#include <functional>

struct md_allocator_i;

namespace task_system {

typedef uint64_t ID;
constexpr ID INVALID_ID = 0;

using execute_fn = std::function<void()>;
using execute_range_fn = std::function<void(uint32_t range_beg, uint32_t range_end)>;

void initialize();
void shutdown();

// Call once per frame at some approriate time, if there are items in the main queue, the main thread will be stalled.
// Pool tasks will not stall the main thread.
void execute_tasks();

// This is to generate tasks for the main thread ("render" thread)
ID main_enqueue(const char* label, execute_fn task, ID dependency = 0);

// This is to generate tasks for the thread-pool (async operations)
ID pool_enqueue(const char* label, execute_fn task, ID dependency = 0);
ID pool_enqueue(const char* label, uint32_t range_size, execute_range_fn range_task, ID dependency = 0);

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