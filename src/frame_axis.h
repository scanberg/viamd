#pragma once

#include <stddef.h>

// Frame <-> time along a frame axis: the time of every frame, in double, never decreasing. The run's
// own axis is one ("<run>/time", which md_run_publish keeps from decreasing), and so is any other axis
// of a run (an energy file's, a section of a trajectory written at its own rate). A fractional frame
// lies on the straight line between the two frames around it, so the two functions are each other's
// inverse wherever the time strictly increases.
//
// These are the only conversions between the two: everything that relates a frame to a time, or a
// point on the timeline to a frame, goes through them, in the axis' own unit. A display unit is
// applied where a time is shown, never to what is converted here.

// The time at a fractional frame, clamped to the first and the last frame
static inline double time_at_frame(const double* times, size_t count, double frame) {
    if (count == 0) return 0.0;
    if (!(frame > 0.0)) return times[0];
    if (frame >= (double)(count - 1)) return times[count - 1];
    const size_t f = (size_t)frame;
    return times[f] + (times[f + 1] - times[f]) * (frame - (double)f);
}

// The fractional frame at a time, clamped to the first and the last frame. Of frames which share a
// time, the first: there is no time between them to place anything at.
static inline double frame_at_time(const double* times, size_t count, double time) {
    if (count == 0 || !(time > times[0])) return 0.0;
    if (time > times[count - 1]) return (double)(count - 1);
    // times[lo] < time <= times[hi], until hi is the first frame at or after the time
    size_t lo = 0;
    size_t hi = count - 1;
    while (hi - lo > 1) {
        const size_t mid = lo + (hi - lo) / 2;
        if (times[mid] < time) lo = mid;
        else hi = mid;
    }
    return (double)lo + (time - times[lo]) / (times[hi] - times[lo]);
}
