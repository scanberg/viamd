#pragma once

// Named spans of frames on the timeline, which the user marks out, names and annotates: a period
// covers frames [beg, end], a marker is a period of a single frame. The timeline owns them and they
// are saved with the workspace. What is shown over a period refers to it by id: a script property
// over its frames ("r (Production)"), a period of a distribution.
//
// A period is either fixed, its ends are frames, or it follows the current frame, its ends are
// offsets from it: -60..0 is the last 60 frames, -25..25 is centered on the current frame. The two
// ends are kept in whatever order they were dragged into, so an end dragged past the other keeps
// going; everything that reads a period takes the lower one as the beginning.
//
// This unit is free of ImGui so that it can be tested; timeline_periods_ui.h draws and edits them.

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
    vec4_t   color = vec4_set(1, 1, 1, 1);
    bool     follow = false;                    // The ends are offsets from the current frame
    int32_t  beg = 0;
    int32_t  end = 0;
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

// The frames [beg, end] the period covers with the current frame where it is, clamped to the
// trajectory. False when it covers none of its frames: a window ahead of the last frame, or a
// fixed period beyond the end of a shorter trajectory.
bool timeline_period_frames(int* out_beg, int* out_end, const TimelinePeriod& r, double current_frame, int num_frames);

// The lower end first
void timeline_period_normalize(TimelinePeriod* r);

// Fixed or following, covering the frames it covers now (clamped to the trajectory when it stops following)
void timeline_period_set_follow(TimelinePeriod* r, bool follow, double current_frame, int num_frames);

// Adds a period and returns it, NULL when there is no room for another. Without a name it gets one
// of its own ("Period 2", "Marker 3"); its colour is the one of the palette the fewest periods have.
TimelinePeriod* timeline_periods_add(TimelinePeriods* rs, const char* name, int32_t beg, int32_t end, bool follow);
void           timeline_periods_remove(TimelinePeriods* rs, uint32_t id);
void           timeline_periods_clear(TimelinePeriods* rs);
int            timeline_periods_index(const TimelinePeriods* rs, uint32_t id);
TimelinePeriod* timeline_periods_find(TimelinePeriods* rs, uint32_t id);
const TimelinePeriod* timeline_periods_find(const TimelinePeriods* rs, uint32_t id);

// "Period 1", "Marker 1": the lowest number no period is named with yet
void timeline_periods_default_name(char* buf, size_t cap, const TimelinePeriods* rs, bool marker);

// ## Lanes
//
// The order periods are drawn in, bottom first: the ones covering the most frames, and of two that
// cover equally many the older one. Shorter periods, and markers, end up on top, where they can be
// grabbed. Periods that cover no frames right now are left out. Returns how many were written.
int timeline_periods_draw_order(int* out_indices, const TimelinePeriods* rs, double current_frame, int num_frames);

// How deep each period in the draw order is nested in one lane: how many periods drawn before it
// overlap it. Written per position in the order.
void timeline_periods_nesting(int* out_depth, const int* order, int count, const TimelinePeriods* rs, double current_frame, int num_frames);

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
// One [TimelinePeriod] section per period, in the order of the list. A period keeps its id across a save and a load, so
// what refers to it by id (a plotted series) finds it again.
void timeline_periods_serialize(viamd::serialization_state_t& state, const TimelinePeriods* rs);
// Reads the section the state is at as one more period. False when there is no room for it.
bool timeline_periods_deserialize_period(viamd::deserialization_state_t& state, TimelinePeriods* rs);
