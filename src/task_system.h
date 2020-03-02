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

typedef std::function<void(TaskSetRange range, TaskData data)> TaskSetFunction;

using ID = uint32_t;

void initialize();
void shutdown();

uint32_t get_num_threads();
uint32_t get_num_tasks();
ID* get_tasks();

void clear_completed_tasks();

ID create_task(const char* label, TaskSetFunction func);
ID create_task(const char* label, uint32_t size, TaskSetFunction func);

bool get_task_complete(ID);
const char* get_task_label(ID);
float get_task_fraction_complete(ID);
void wait_for_task(ID);
void interrupt_task(ID);
void interrupt_and_wait(ID);


}  // namespace task_system