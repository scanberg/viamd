#pragma once

#include <core/types.h>
#include <core/string_utils.h>

namespace cpu_profiling {
void initialize();
void shutdown();

void push_section(CString sec);
void pop_section(CString sec);

void clear();
void finish();
void draw_window();

struct ScopedSection {
    ScopedSection(CString section) : sec(section) { push_section(sec); }
    ~ScopedSection() { pop_section(sec); }
    CString sec;
};
}  // namespace cpu_profiling

namespace gpu_profiling {
void initialize();
void shutdown();

void push_section(CString sec);
void pop_section(CString sec);

void clear();
void finish();
void draw_window();
void print();

struct ScopedSection {
    ScopedSection(CString section) : sec(section) { push_section(sec); }
    ~ScopedSection() { pop_section(sec); }
    CString sec;
};

}  // namespace gpu_profiling
