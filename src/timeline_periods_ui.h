#pragma once

// Drawing and editing the timeline's periods (timeline_periods.h): a lane above the plots, or a lane
// per period, where they are created, moved, resized, named and removed directly, and their shading
// in the plots below.
//
//   shift-drag     marks a new period, in the lane or in any plot; shift-click marks a single frame
//   drag           a handle resizes a period, the bar moves it; an end dragged past the other keeps going
//   double-click   renames
//   right-click    the annotation, the colour, whether it follows the current frame, delete
//   delete         removes the period under the mouse
//   pulled out     a period dragged well out of the lane is taken elsewhere (TIMELINE_PERIOD_DND), and stays put:
//                  dropped on a distribution, it is shown over the period too
//
// Empty lane space pans and zooms the timeline like the plots do.

#include <timeline_periods.h>

// How frames map to the screen, shared with the timeline's plots
struct TimelinePeriodAxis {
    const float* frame_times = nullptr;     // per frame, in the timeline's unit (timeline.x_values)
    int          num_frames = 0;
    double*      view_beg = nullptr;        // the visible span of time, which the lane pans and zooms
    double*      view_end = nullptr;
    const char*  time_unit = "";            // "ns"; empty when the frames have no time of their own
};

// Call where the lane goes, before the plots of the timeline. Ctrl-drag in it moves current_frame.
void timeline_periods_lane(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double* current_frame);

// Call inside each plot of the timeline, after its setup: shades the periods in the plot, and starts
// a new period on shift-press. The first plot each frame tells the lane where the plots' area is.
void timeline_periods_plot(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double current_frame);

// The contents of the timeline's Periods menu: the layout of the lanes, the shading
void timeline_periods_menu(TimelinePeriods* rs);

// While a period is created, moved or resized, or the lane pans: the plots should not pan meanwhile
bool timeline_periods_dragging();

// ImGui drag and drop payload of a period pulled out of the lane: its id, a uint32_t
#define TIMELINE_PERIOD_DND "TIMELINE_PERIOD"

// The period the pointer is on in the lane, or holds (moving it, resizing it, or taking it elsewhere): what other views
// may show a preview of. 0 for none.
uint32_t timeline_periods_hovered();

// Lights the period up in the lane and the plots, as if it were hovered, for this frame and the next: a view showing
// something over it says which period it is
void timeline_periods_highlight(uint32_t id);
