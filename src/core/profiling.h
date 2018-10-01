#pragma once

#include <core/types.h>
#include <core/string_utils.h>

namespace cpu_profiling {
void initialize();
void shutdown();

void push_section(CString label);
void pop_section();

void clear();
void finish();
void draw_window(bool* window_open);
}  // namespace cpu_profiling

namespace gpu_profiling {
void initialize();
void shutdown();

void push_section(CString label);
void pop_section();

void begin_frame();
void end_frame();

void draw_window(bool* window_open);
void print();

}  // namespace gpu_profiling
