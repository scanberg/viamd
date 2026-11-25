#pragma once

// THE CLOCK: the time shown, which every instance that follows the time shows its run at (viamd.h, ApplicationState::
// animation.time). One for the application, in the unit of the reference's run: the run whose frames stepping, the
// frame slider and the speed (frames per second) count. A run is placed on it by its own times (frame_axis.h), taken
// into the clock's unit: a run sampled every 800 ps beside one every 200 ps shows the frame of the same time, not the
// frame of the same number. A run whose times cannot be taken into the clock's unit (frames numbered beside times, or
// the other way around) is placed frame by frame as the reference: its frame is the reference's at the time. A run is
// moved along the clock by an offset of its own: its start aligned with another's, a continuation whose times read from
// 0 again placed after the run it continues.
//
// This unit is free of the application so that it can be tested; viamd.h has the clock of the application.

#include <frame_axis.h>

#include <math.h>
#include <stddef.h>

// A run's frame axis as placed on the clock: its times, in its own unit and never decreasing, the offset it is moved by
// along the clock, in the same unit, and the factor which takes them into the clock's unit. A scale of 0 for a run
// which is not placed by its times.
struct ClockAxis {
    const double* times = nullptr;
    size_t count = 0;
    double scale = 1.0;     // clock time = scale * (the run's time + offset)
    double offset = 0.0;
};

static inline bool clock_axis_timed(const ClockAxis& a) { return a.scale > 0.0 && a.times && a.count > 0; }

// The time on the clock at a fractional frame of a timed axis, clamped to its first and last frame
static inline double clock_axis_time(const ClockAxis& a, double frame) {
    return (time_at_frame(a.times, a.count, frame) + a.offset) * a.scale;
}

// The fractional frame of a timed axis at a time on the clock, clamped to its frames
static inline double clock_axis_frame(const ClockAxis& a, double time) {
    return frame_at_time(a.times, a.count, time / a.scale - a.offset);
}

// The offset which puts the first frame of a timed axis at a time on the clock: its start aligned with another's
static inline double clock_axis_offset_to_start_at(const ClockAxis& a, double time) {
    return time / a.scale - a.times[0];
}

// The mean time from frame to frame of a timed axis, on the clock: what playback advances the clock by per frame of the
// reference. 1 for fewer than two frames, or none between its first and its last.
static inline double clock_axis_frame_duration(const ClockAxis& a) {
    if (!clock_axis_timed(a) || a.count < 2) return 1.0;
    const double d = (a.times[a.count - 1] - a.times[0]) * a.scale / (double)(a.count - 1);
    return d > 0.0 ? d : 1.0;
}

// The time on the clock moved by a number of frames of the reference (fractional, negative backwards), within
// [beg, end]: to the reference's own frame while the time is within the reference's frames and so is the one moved to,
// else by its mean time per frame (another run goes on past the reference's last frame, or begins before its first).
static inline double clock_step_time(const ClockAxis& ref, double time, double frames, double beg, double end) {
    double t = time + frames * clock_axis_frame_duration(ref);
    if (clock_axis_timed(ref)) {
        const double first = clock_axis_time(ref, 0.0);
        const double last  = clock_axis_time(ref, (double)(ref.count - 1));
        const double to    = clock_axis_frame(ref, time) + frames;
        if (time >= first && time <= last && to >= 0.0 && to <= (double)(ref.count - 1)) t = clock_axis_time(ref, to);
    }
    return t < beg ? beg : t > end ? end : t;
}

// What playback does at the ends of the span it goes over (playback_range): stops there, goes on from the other end, or
// turns back. The values are those a workspace stores ([Animation] Loop), so a new one goes last.
enum PlaybackLoop : int {
    PlaybackLoop_Once = 0,
    PlaybackLoop_Loop,
    PlaybackLoop_Bounce,
    PlaybackLoop_Count
};

// The time moved by delta (negative backwards) in the direction it is played in (+1, -1), within [beg, end] as the loop
// has it at the ends: once, to the end it reaches, which ends the playing (false); looping, on from the other end; and
// bouncing, back from it, the direction turned. A span of nothing ends it at once.
static inline bool clock_advance(double* time, double delta, double beg, double end, PlaybackLoop loop, int* direction) {
    const double len = end - beg;
    const double t = *time + delta * (double)*direction;
    if (!(len > 0.0)) {
        *time = beg;
        return false;
    }
    if (t >= beg && t <= end) {
        *time = t;
        return true;
    }
    switch (loop) {
    case PlaybackLoop_Loop: {
        double u = fmod(t - beg, len);
        if (u < 0.0) u += len;
        *time = beg + u;
        return true;
    }
    case PlaybackLoop_Bounce: {
        // Reflected off the ends as many times as it went past one: an odd number turns it
        double u = fmod(t - beg, 2.0 * len);
        if (u < 0.0) u += 2.0 * len;
        *time = u <= len ? beg + u : beg + 2.0 * len - u;
        const double k = floor((t - beg) / len);
        if (fmod(fabs(k), 2.0) == 1.0) *direction = -*direction;
        return true;
    }
    default:
        *time = t < beg ? beg : end;
        return false;
    }
}
