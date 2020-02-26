#pragma once

#include <stdint.h>
#include <functional>
#include <atomic>

namespace task {

struct TaskSetRange {
    uint32_t beg;
    uint32_t end;
};

struct TaskData {
    uint32_t threadnum;
    std::atomic_bool& interrupt;
};

typedef std::function<void(TaskSetRange range, TaskData data)> TaskSetFunction;

using TaskID = uint32_t;

void initialize();
void shutdown();

uint32_t get_num_tasks();
TaskID* get_tasks();

void clear_completed_tasks();

TaskID create_task(const char* label, TaskSetFunction func);
TaskID create_task(const char* label, uint32_t size, TaskSetFunction func);

bool get_task_complete(TaskID);
const char* get_task_label(TaskID);
float get_task_fraction_complete(TaskID);
void wait_for_task(TaskID);
void interrupt_task(TaskID);
void interrupt_and_wait(TaskID);


}  // namespace task