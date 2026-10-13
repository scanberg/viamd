#pragma once

// Named spans of frames on the timeline, which the user marks out, names and annotates: a range
// covers frames [beg, end], a marker is a range of a single frame. The timeline owns them and they
// are saved with the workspace. One of them can be marked as the timeline filter, which then
// follows it (see timeline_ranges_filter_frames).
//
// A range is either fixed, its ends are frames, or it follows the current frame, its ends are
// offsets from it: -60..0 is the last 60 frames, -25..25 is centered on the current frame. The two
// ends are kept in whatever order they were dragged into, so an end dragged past the other keeps
// going; everything that reads a range takes the lower one as the beginning.
//
// This unit is free of ImGui so that it can be tested; timeline_ranges_ui.h draws and edits them.

#include <core/md_vec_math.h>

#include <stdint.h>
#include <stddef.h>

#define TIMELINE_MAX_RANGES     64
#define TIMELINE_RANGE_NAME_CAP 64
#define TIMELINE_RANGE_NOTE_CAP 512 // 64 of them live in ApplicationState, which is on the stack of main()

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

struct TimelineRange {
    uint32_t id = 0;                            // Unique among the ranges for as long as they live, never 0
    char     name[TIMELINE_RANGE_NAME_CAP] = "";
    char     note[TIMELINE_RANGE_NOTE_CAP] = "";  // The annotation shown when it is hovered
    vec4_t   color = vec4_set(1, 1, 1, 1);
    bool     follow = false;                    // The ends are offsets from the current frame
    int32_t  beg = 0;
    int32_t  end = 0;
};

enum TimelineRangeLayout : int {
    TimelineRangeLayout_OneLane = 0,    // Every range in one lane, the shorter ones on top of the longer
    TimelineRangeLayout_LanePerRange,   // A lane of its own for each range
    TimelineRangeLayout_Count
};

struct TimelineRanges {
    TimelineRange items[TIMELINE_MAX_RANGES];
    int      count = 0;
    uint32_t next_id = 1;
    uint32_t filter_id = 0;             // The range the timeline filter follows, 0 for none

    TimelineRangeLayout layout = TimelineRangeLayout_OneLane;
    bool     shade_all = false;         // Shade every range in the plots, not only the one under the mouse
};

static inline bool timeline_range_is_marker(const TimelineRange& r) { return r.beg == r.end; }

// The frames [beg, end] the range covers with the current frame where it is, clamped to the
// trajectory. False when it covers none of its frames: a window ahead of the last frame, or a
// fixed range beyond the end of a shorter trajectory.
bool timeline_range_frames(int* out_beg, int* out_end, const TimelineRange& r, double current_frame, int num_frames);

// The lower end first
void timeline_range_normalize(TimelineRange* r);

// Fixed or following, covering the frames it covers now (clamped to the trajectory when it stops following)
void timeline_range_set_follow(TimelineRange* r, bool follow, double current_frame, int num_frames);

// Adds a range and returns it, NULL when there is no room for another. Without a name it gets one
// of its own ("Range 2", "Marker 3"); its colour is the one of the palette the fewest ranges have.
TimelineRange* timeline_ranges_add(TimelineRanges* rs, const char* name, int32_t beg, int32_t end, bool follow);
void           timeline_ranges_remove(TimelineRanges* rs, uint32_t id);
void           timeline_ranges_clear(TimelineRanges* rs);
int            timeline_ranges_index(const TimelineRanges* rs, uint32_t id);
TimelineRange* timeline_ranges_find(TimelineRanges* rs, uint32_t id);
const TimelineRange* timeline_ranges_find(const TimelineRanges* rs, uint32_t id);

// "Range 1", "Marker 1": the lowest number no range is named with yet
void timeline_ranges_default_name(char* buf, size_t cap, const TimelineRanges* rs, bool marker);

// The frames of the range the filter follows. False when no range is the filter, or it covers no frames.
bool timeline_ranges_filter_frames(int* out_beg, int* out_end, const TimelineRanges* rs, double current_frame, int num_frames);

// ## Lanes
//
// The order ranges are drawn in, bottom first: the ones covering the most frames, and of two that
// cover equally many the older one. Shorter ranges, and markers, end up on top, where they can be
// grabbed. Ranges that cover no frames right now are left out. Returns how many were written.
int timeline_ranges_draw_order(int* out_indices, const TimelineRanges* rs, double current_frame, int num_frames);

// How deep each range in the draw order is nested in one lane: how many ranges drawn before it
// overlap it. Written per position in the order.
void timeline_ranges_nesting(int* out_depth, const int* order, int count, const TimelineRanges* rs, double current_frame, int num_frames);

enum TimelineRangePart : int {
    TimelineRangePart_None = 0,
    TimelineRangePart_Beg,      // The handle at its lower end
    TimelineRangePart_End,      // The handle at its upper end
    TimelineRangePart_Body,     // The bar itself, or a marker
};

// Where a range is drawn, in pixels. A marker is drawn around x0, which equals x1.
struct TimelineRangeBar {
    float x0 = 0, x1 = 0;
    float y0 = 0, y1 = 0;
    bool  marker = false;
    bool  visible = false;
};

struct TimelineRangeHit {
    int index = -1;
    TimelineRangePart part = TimelineRangePart_None;
};

// What the pointer would grab. bars is indexed like the ranges, top_first lists the visible ones in
// the order they are on top of each other. Markers first, as they are the smallest targets. Then
// the handles: the nearest within reach, and of two at the same place (a range ending where the next
// begins) the one whose range lies on the pointer's side. Then the bars, the top one first.
TimelineRangeHit timeline_ranges_hit_test(const TimelineRangeBar* bars, const int* top_first, int count, float px, float py,
                                          float handle_reach, float marker_reach);

// ## Workspace
//
// One [TimelineRange] section per range, in the order of the list.
void timeline_ranges_serialize(viamd::serialization_state_t& state, const TimelineRanges* rs);
// Reads the section the state is at as one more range. False when there is no room for it.
bool timeline_ranges_deserialize_range(viamd::deserialization_state_t& state, TimelineRanges* rs);
