#pragma once

// Named spans of time on the timeline, which the user marks out, names and annotates: a period covers the time [beg, end]
// on the clock (clock.h), in its unit, a marker is a period of a single time. The timeline owns them and they are saved
// with the workspace. What is shown over a period refers to it by id: a script property over its frames ("r
// (Production)"), a period of a distribution. Each run takes the frames of its own the period covers on the clock
// (timeline_period_frames): a run sampled every 800 ps the same stretch of time as one every 200 ps, not the same frames.
//
// A period is either fixed, its ends are times, or it follows the current time, its ends are offsets from it: -60..0
// ns is the last 60 ns, -25..25 ns is centered on the current time. The two ends are kept in whatever order they were
// dragged into, so an end dragged past the other keeps going; everything that reads a period takes the lower one as the
// beginning.
//
// This unit is free of ImGui so that it can be tested; timeline_periods_ui.h draws and edits them.

#include <clock.h>
#include <core/md_vec_math.h>

#include <stdint.h>
#include <stddef.h>

#define TIMELINE_MAX_PERIODS     64
#define TIMELINE_PERIOD_NAME_CAP 64
#define TIMELINE_PERIOD_NOTE_CAP 512 // 64 of them live in ApplicationState, which is on the stack of main()

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

struct TimelinePeriod {
    uint32_t id = 0;                            // Unique among the periods for as long as they live, never 0
    char     name[TIMELINE_PERIOD_NAME_CAP] = "";
    char     note[TIMELINE_PERIOD_NOTE_CAP] = "";  // The annotation shown when it is hovered
    md_float4_t   color = md_float4_set(1, 1, 1, 1);
    bool     follow = false;                    // The ends are offsets from the current time
    double   beg = 0.0;                         // On the clock, in its unit
    double   end = 0.0;
};

enum TimelinePeriodLayout : int {
    TimelinePeriodLayout_OneLane = 0,    // Every period in one lane, the shorter ones on top of the longer
    TimelinePeriodLayout_LanePerPeriod,   // A lane of its own for each period
    TimelinePeriodLayout_Count
};

struct TimelinePeriods {
    TimelinePeriod items[TIMELINE_MAX_PERIODS];
    int      count = 0;
    uint32_t next_id = 1;
    // The period made from the timeline filter of a workspace from before there were periods, while it is read: what was
    // shown over the filter is shown over it. 0 for none. Not saved.
    uint32_t legacy_filter_id = 0;

    TimelinePeriodLayout layout = TimelinePeriodLayout_OneLane;
    bool     shade_all = false;         // Shade every period in the plots, not only the one under the mouse
};

static inline bool timeline_period_is_marker(const TimelinePeriod& r) { return r.beg == r.end; }

// The time [beg, end] on the clock the period covers with the current time where it is, the lower end first
void timeline_period_span(double* out_beg, double* out_end, const TimelinePeriod& r, double current_time);

// The frames [beg, end] of a run the period covers, the run placed on the clock by its axis: those whose times lie
// within the period, and of a period which lies between two of them (a marker), the one nearest its middle. False when
// it covers none of the run: it is before its first frame or after its last, or the axis has no frames.
bool timeline_period_frames(int* out_beg, int* out_end, const TimelinePeriod& r, double current_time, const ClockAxis& axis);

// The lower end first
void timeline_period_normalize(TimelinePeriod* r);

// Fixed or following. Made to follow, it keeps its extent and is centered on the current time, on a whole number of
// steps (frames of the reference) to each side where it can be, the extra one on the side of the past; made fixed, it
// covers the time it covers now, within the clock's span [span_beg, span_end].
void timeline_period_set_follow(TimelinePeriod* r, bool follow, double current_time, double step, double span_beg, double span_end);

// A period read from a workspace from before the clock, whose ends are frames of the reference's run (offsets in frames
// of it, when it follows), taken into time on the clock by the reference's axis
void timeline_period_frames_to_time(TimelinePeriod* r, const ClockAxis& reference);

// Adds a period and returns it, NULL when there is no room for another. Without a name it gets one
// of its own ("Period 2", "Marker 3"); its colour is the one of the palette the fewest periods have.
TimelinePeriod* timeline_periods_add(TimelinePeriods* rs, const char* name, double beg, double end, bool follow);
void           timeline_periods_remove(TimelinePeriods* rs, uint32_t id);
void           timeline_periods_clear(TimelinePeriods* rs);
int            timeline_periods_index(const TimelinePeriods* rs, uint32_t id);
TimelinePeriod* timeline_periods_find(TimelinePeriods* rs, uint32_t id);
const TimelinePeriod* timeline_periods_find(const TimelinePeriods* rs, uint32_t id);

// "Period 1", "Marker 1": the lowest number no period is named with yet
void timeline_periods_default_name(char* buf, size_t cap, const TimelinePeriods* rs, bool marker);

// ## Lanes
//
// The order periods are drawn in, bottom first: the ones covering the most time, and of two that cover equally much the
// older one. Shorter periods, and markers, end up on top, where they can be grabbed. Periods that cover none of the
// clock's span [span_beg, span_end] right now are left out. Returns how many were written.
int timeline_periods_draw_order(int* out_indices, const TimelinePeriods* rs, double current_time, double span_beg, double span_end);

// How deep each period in the draw order is nested in one lane: how many periods drawn before it
// overlap it. Written per position in the order.
void timeline_periods_nesting(int* out_depth, const int* order, int count, const TimelinePeriods* rs, double current_time, double span_beg, double span_end);

enum TimelinePeriodPart : int {
    TimelinePeriodPart_None = 0,
    TimelinePeriodPart_Beg,      // The handle at its lower end
    TimelinePeriodPart_End,      // The handle at its upper end
    TimelinePeriodPart_Body,     // The bar itself, or a marker
};

// Where a period is drawn, in pixels. A marker is drawn around x0, which equals x1.
struct TimelinePeriodBar {
    float x0 = 0, x1 = 0;
    float y0 = 0, y1 = 0;
    bool  marker = false;
    bool  visible = false;
};

struct TimelinePeriodHit {
    int index = -1;
    TimelinePeriodPart part = TimelinePeriodPart_None;
};

// What the pointer would grab. bars is indexed like the periods, top_first lists the visible ones in
// the order they are on top of each other. Markers first, as they are the smallest targets. Then
// the handles: the nearest within reach, and of two at the same place (a period ending where the next
// begins) the one whose period lies on the pointer's side. Then the bars, the top one first.
TimelinePeriodHit timeline_periods_hit_test(const TimelinePeriodBar* bars, const int* top_first, int count, float px, float py,
                                          float handle_reach, float marker_reach);

// ## Workspace
//
// One [TimelinePeriod] section per period, in the order of the list, its ends times in the unit given (Time, TimeUnit).
// A period keeps its id across a save and a load, so what refers to it by id (a plotted series) finds it again.
void timeline_periods_serialize(viamd::serialization_state_t& state, const TimelinePeriods* rs, const char* time_unit);
// Reads the section the state is at as one more period, returned. NULL when there is no room for it. Its ends as they
// were written: times in the unit written to time_unit (cap bytes, empty for none), or, from a workspace from before
// the clock, frames (Frames), which *in_frames says (timeline_period_frames_to_time).
TimelinePeriod* timeline_periods_deserialize_period(viamd::deserialization_state_t& state, TimelinePeriods* rs, char* time_unit, size_t cap, bool* in_frames);
