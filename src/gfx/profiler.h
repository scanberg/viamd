#pragma once

// PROFILER
// A frame in nested sections: PUSH_GPU_SECTION("name") ... POP_GPU_SECTION() around a part of it, or
// PROFILE_SCOPE("name") for the rest of a block. Each section is timed on the CPU, from its push to its pop, in
// milliseconds per frame, averaged over a number of frames, shown in Windows > Profiler. The sections of one name under
// one parent are one: their times add up (a pass per representation, per view).
//
// The GPU's time is not measured: the frame is drawn with md_gpu, whose work a CPU section only records. Sections are
// for the main thread. Outside frame_begin() / frame_end(), and nested deeper than the profiler keeps, they are not timed.

struct md_writer_i;

namespace profiler {

// The frame: its root section, around all that a frame does. frame_begin() also publishes the averages of the frames
// before it, once there are enough.
void frame_begin();
void frame_end();
void shutdown();

void push(const char* name);
void pop();

// For the rest of the block
struct Scope {
    explicit Scope(const char* name) { push(name); }
    ~Scope() { pop(); }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
};

// The published times per frame, in ms, 0 until there are any: of the sections named 'name' (all of them, added up),
// or only of those within a section named 'within' (anywhere above them)
float cpu_ms(const char* name, const char* within = nullptr);

// Of the whole frame
float frame_cpu_ms();

// The sections as a tree, a row each with its times: what Windows > Profiler shows (within a window)
void draw_table();
// The frame and the sections at its top, a few lines: what the tooltip of the frame time shows
void draw_summary();

// The published times as CSV, what Windows > Profiler exports: '#' lines on what they are of (when, the build, the
// frames averaged, then the caller's 'info', a line each, NULL for none), the header and a row per section shown, in the
// order of the table: its path ("Frame/Render/..."), its CPU time in ms per frame, and its calls per frame. Before
// anything is published, the header alone.
void write_csv(md_writer_i* w, const char* info);

}  // namespace profiler

#define PROFILER_CONCAT_(a, b) a##b
#define PROFILER_CONCAT(a, b)  PROFILER_CONCAT_(a, b)

// Statements of their own, as they always were: with or without a semicolon after them
#define PUSH_GPU_SECTION(name) { profiler::push(name); }
#define POP_GPU_SECTION()      { profiler::pop(); }
#define PROFILE_SCOPE(name)    profiler::Scope PROFILER_CONCAT(profile_scope_, __LINE__)(name)
