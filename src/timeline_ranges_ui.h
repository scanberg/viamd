#pragma once

// Drawing and editing the timeline's ranges (timeline_ranges.h): a lane above the plots, or a lane
// per range, where they are created, moved, resized, named and removed directly, and their shading
// in the plots below.
//
//   shift-drag     marks a new range, in the lane or in any plot; shift-click marks a single frame
//   drag           a handle resizes a range, the bar moves it; an end dragged past the other keeps going
//   double-click   renames
//   right-click    the annotation, the colour, whether it follows the current frame, the filter, delete
//   delete         removes the range under the mouse
//
// Empty lane space pans and zooms the timeline like the plots do.

#include <timeline_ranges.h>

// How frames map to the screen, shared with the timeline's plots
struct TimelineRangeAxis {
    const float* frame_times = nullptr;     // per frame, in the timeline's unit (timeline.x_values)
    int          num_frames = 0;
    double*      view_beg = nullptr;        // the visible span of time, which the lane pans and zooms
    double*      view_end = nullptr;
    const char*  time_unit = "";            // "ns"; empty when the frames have no time of their own
};

// Call where the lane goes, before the plots of the timeline. Ctrl-drag in it moves current_frame.
void timeline_ranges_lane(TimelineRanges* rs, const TimelineRangeAxis& axis, double* current_frame);

// Call inside each plot of the timeline, after its setup: shades the ranges in the plot, and starts
// a new range on shift-press. The first plot each frame tells the lane where the plots' area is.
void timeline_ranges_plot(TimelineRanges* rs, const TimelineRangeAxis& axis, double current_frame);

// The contents of the timeline's Ranges menu: the layout of the lanes, the shading, the filter
void timeline_ranges_menu(TimelineRanges* rs, double current_frame, int num_frames);

// While a range is created, moved or resized, or the lane pans: the plots should not pan meanwhile
bool timeline_ranges_dragging();
