#pragma once

#include <stdint.h>
#include <functional>
#include <atomic>

namespace task_system {

struct TaskSetRange {
    uint32_t beg;
    uint32_t end;
};

typedef std::function<void()> Task;
typedef std::function<void(TaskSetRange range)> TaskSet;

struct ID {
    uint32_t id;
};

constexpr ID INVALID_ID = {0};

void initialize();
void shutdown();

// This is to generate tasks for the main thread (render thread)
bool enqueue_main(const char* label, Task task);
void run_main_tasks();

// This is to generate tasks for the thread-pool (async operations)
ID enqueue_pool(const char* label, uint32_t size, TaskSet task_set, TaskSet on_complete = nullptr);

uint32_t get_num_threads();
uint32_t get_num_tasks();
ID* get_tasks();

bool get_task_complete(ID);
const char* get_task_label(ID);
float get_task_fraction_complete(ID);

void wait_for_task(ID);
void interrupt_task(ID);
void interrupt_and_wait(ID);

}  // namespace task_system