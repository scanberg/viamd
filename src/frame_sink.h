#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>

#include <string>
#include <vector>

// Takes finished movie frames off the render thread. Frames are handed over as RGBA pixels, bottom row
// first (the order OpenGL reads them), and are written by worker threads, either as a numbered PNG
// sequence or piped straight into an ffmpeg process that encodes an mp4. A fixed pool of frame buffers
// bounds the memory in use: when the writers fall behind, acquire() has nothing to give and the caller
// has to wait, rather than the queue growing without limit.
namespace frame_sink {

enum class Kind {
    PngSequence,  // <dir>/<prefix>_00000.png, ...
    Ffmpeg,       // <dir>/<prefix>.mp4
};

// What ffmpeg encodes to
enum class Codec {
    H264,  // libx264 in an mp4
    H265,  // libx265 in an mp4, tagged hvc1 so that QuickTime plays it
    Vp9,   // libvpx-vp9 in a webm
};

// The extension of the file ffmpeg writes, without the dot
const char* file_extension(Codec codec);

struct Desc {
    Kind  kind = Kind::PngSequence;
    str_t dir = {};
    str_t prefix = {};
    int   width = 0;
    int   height = 0;
    float fps = 24.0f;
    int   max_in_flight = 0;  // Frame buffers, 0 picks a count that keeps the pool near 256 MB

    // Ffmpeg only
    str_t ffmpeg = {};        // Executable, "ffmpeg" (found on PATH) when empty
    Codec codec = Codec::H264;
    int   crf = 18;           // Quality, lower is better and larger
};

struct Status {
    int  submitted = 0;
    int  written = 0;
    int  failed = 0;
    int  queued = 0;      // Submitted and not yet written
    bool done = false;    // Closed and every writer has finished (for ffmpeg: the process has exited)
    bool ok = true;       // Nothing failed so far
    char message[256] = {};
};

struct Sink;

// Returns null with a message in err if the sink could not be set up, e.g. ffmpeg was not found
Sink* create(const Desc& desc, char* err, size_t err_cap);

// A buffer of width * height * 4 bytes to fill and hand to submit(), or null if every buffer is in
// use. With wait set it blocks until one is free.
uint8_t* acquire(Sink* sink, bool wait);

// Takes ownership of a buffer from acquire(). Frames are numbered by index in the file names.
// With ffmpeg they must be submitted in order.
void submit(Sink* sink, uint8_t* frame, int index);

// Gives back a buffer from acquire() that will not be used
void release(Sink* sink, uint8_t* frame);

// No more frames are coming. What is queued is still written, in the background.
void close(Sink* sink);

// Drops what has not been written yet and closes
void cancel(Sink* sink);

Status status(Sink* sink);

// Cancels if still open, then waits for the writers to finish what they have
void destroy(Sink* sink);

// The arguments handed to ffmpeg, the executable first. Exposed for testing.
void ffmpeg_arguments(const Desc& desc, std::vector<std::string>& out);

}  // namespace frame_sink
