#pragma once

// PROFILER
// A frame in nested sections: PUSH_GPU_SECTION("name") ... POP_GPU_SECTION() around a part of it, or
// PROFILE_SCOPE("name") for the rest of a block. Each section is
//   - a debug group (KHR_debug): what Nsight, RenderDoc and apitrace show a captured frame by, where the context has
//     it (macOS does not: OpenGL 4.1 has no debug groups);
//   - timed: on the GPU, the work issued from its push to its pop, and on the CPU, the time from its push to its pop;
//     in milliseconds per frame, averaged over a number of frames, shown in Windows > Profiler.
// The sections of one name under one parent are one: their times add up (a pass per representation, per view).
//
// The GPU time is measured with GL_TIME_ELAPSED queries, read back a few frames late so that reading them never
// stalls. They cannot nest, one runs at a time: a push ends the parent's query and starts the section's, a pop ends
// the section's and starts the parent's again. A section's own time is what ran while it was the innermost, its time
// that and its sections'. So no other GL_TIME_ELAPSED query may run inside a frame: one that has to (in a library)
// runs within a profiler::Pause, which stops the GPU timing until it ends.
//
// Sections are for the thread of the GL context. Outside frame_begin() / frame_end() they are debug groups only, and
// so are sections nested deeper than the profiler keeps.

struct md_writer_i;

namespace profiler {

// The frame: its root section, around all that a frame does. frame_begin() also collects the times of the frames whose
// results are in, and publishes their averages. The queries are made as they are needed (with the GL context current);
// shutdown() deletes them.
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

// Stops the GPU timing until the end of the block, for a GL_TIME_ELAPSED query of someone else's within it. What runs
// meanwhile is timed by no section; the sections pushed are still debug groups and timed on the CPU.
struct Pause {
    Pause();
    ~Pause();
    Pause(const Pause&) = delete;
    Pause& operator=(const Pause&) = delete;
};

// The published times per frame, in ms, 0 until there are any: of the sections named 'name' (all of them, added up),
// or only of those within a section named 'within' (anywhere above them)
float gpu_ms(const char* name, const char* within = nullptr);
float cpu_ms(const char* name, const char* within = nullptr);

// Of the whole frame
float frame_gpu_ms();
float frame_cpu_ms();

// The sections as a tree, a row each with its times: what Windows > Profiler shows (within a window)
void draw_table();
// The frame and the sections at its top, a few lines: what the tooltip of the frame time shows
void draw_summary();

// The published times as CSV, what Windows > Profiler exports: '#' lines on what they are of (when, the build, the GL
// context's renderer, the frames averaged, then the caller's 'info', a line each, NULL for none), the header and a row
// per section shown, in the order of the table: its path ("Frame/Render/..."), its GPU time with its sections' and
// without, its CPU time, in ms per frame, and its calls per frame. Before anything is published, the header alone.
void write_csv(md_writer_i* w, const char* info);

}  // namespace profiler

#define PROFILER_CONCAT_(a, b) a##b
#define PROFILER_CONCAT(a, b)  PROFILER_CONCAT_(a, b)

// Statements of their own, as they always were: with or without a semicolon after them
#define PUSH_GPU_SECTION(name) { profiler::push(name); }
#define POP_GPU_SECTION()      { profiler::pop(); }
#define PROFILE_SCOPE(name)    profiler::Scope PROFILER_CONCAT(profile_scope_, __LINE__)(name)
