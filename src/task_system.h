#pragma once

#include <stdint.h>
#include <functional>
#include <atomic>

namespace task_system {

struct TaskSetRange {
    uint32_t beg;
    uint32_t end;
};

struct TaskData {
    uint32_t threadnum;
    std::atomic_bool& interrupt;
};

typedef std::function<void()> MainFunction;
typedef std::function<void(TaskSetRange range, TaskData data)> TaskSetFunction;

struct ID {
    uint32_t id;
};

constexpr ID INVALID_ID = {0};

void initialize();
void shutdown();

namespace main {
// This is to generate tasks for the main thread (render thread)
bool enqueue(const char* label, MainFunction func);
void run_tasks();
}

namespace pool {

// This is to generate tasks for the thread-pool (async operations)
ID enqueue(const char* label, uint32_t size, TaskSetFunction func);

uint32_t get_num_threads();
uint32_t get_num_tasks();
ID* get_tasks();

bool get_task_complete(ID);
const char* get_task_label(ID);
float get_task_fraction_complete(ID);

void wait_for_task(ID);
void interrupt_task(ID);
void interrupt_and_wait(ID);
}



}  // namespace task_system