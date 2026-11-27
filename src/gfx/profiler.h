#pragma once

// PROFILER
// A frame in nested sections: PUSH_GPU_SECTION("name") ... POP_GPU_SECTION() around a part of it, or
// PROFILE_SCOPE("name") for the rest of a block. Each section is timed on the CPU, from its push to its pop, in
// milliseconds per frame, averaged over a number of frames, shown in Windows > Profiler. The sections of one name under
// one parent are one: their times add up (a pass per representation, per view).
//
// With GPU timing on (set_gpu), the work recorded within a section on the streams given (set_gpu_streams) is timed on
// the GPU as well: each launch, render pass and copy (md_gpu_stream_profile), by the name of its kernel or pass, under
// the section it was recorded in. A dotted name is a path ("md_compose.ao_estimate": "md_compose", within it
// "ao_estimate"), and a name a section of the same name under the same parent goes into it. A section's GPU time is
// that of all within it. The GPU's times arrive a few frames after their CPU's, averaged the same way.
//
// Sections are for the main thread. Outside frame_begin() / frame_end(), and nested deeper than the profiler keeps,
// they are not timed.

struct md_writer_i;
struct md_gpu_stream;

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

// The streams whose work is timed on the GPU (at most 4), and whether it is: a timestamp per operation while on
void set_gpu_streams(md_gpu_stream* const* streams, int count);
void set_gpu(bool on);

// The GPU's time per frame, in ms, 0 while it is not timed: of the sections named 'name', or of the whole frame
float gpu_ms(const char* name, const char* within = nullptr);
float frame_gpu_ms();

// The sections as a tree, a row each with its times: what Windows > Profiler shows (within a window)
void draw_table();
// The frame and the sections at its top, a few lines: what the tooltip of the frame time shows
void draw_summary();

// The published times as CSV, what Windows > Profiler exports: '#' lines on what they are of (when, the build, the
// frames averaged, then the caller's 'info', a line each, NULL for none), the header and a row per section shown, in the
// order of the table: its path ("Frame/Render/..."), its CPU and GPU times in ms per frame (empty where not timed), and
// its calls per frame (of a GPU operation: how often it ran). Before anything is published, the header alone.
void write_csv(md_writer_i* w, const char* info);

}  // namespace profiler

#define PROFILER_CONCAT_(a, b) a##b
#define PROFILER_CONCAT(a, b)  PROFILER_CONCAT_(a, b)

// Statements of their own, as they always were: with or without a semicolon after them
#define PUSH_GPU_SECTION(name) { profiler::push(name); }
#define POP_GPU_SECTION()      { profiler::pop(); }
#define PROFILE_SCOPE(name)    profiler::Scope PROFILER_CONCAT(profile_scope_, __LINE__)(name)
