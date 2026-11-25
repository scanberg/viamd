#pragma once

// Drawing and editing the timeline's periods (timeline_periods.h): a lane above the plots, or a lane
// per period, where they are created, moved, resized, named and removed directly, and their shading
// in the plots below.
//
//   shift-drag     marks a new period, in the lane or in any plot; shift-click marks a single time
//   drag           a handle resizes a period, the bar moves it; an end dragged past the other keeps going
//                  (the ends go to the reference's frames, and beyond them by its mean step)
//   double-click   renames
//   right-click    the annotation, the colour, whether it follows the current time, delete
//   delete         removes the period under the mouse
//   pulled out     a period dragged well out of the lane is taken elsewhere (TIMELINE_PERIOD_DND), and stays put:
//                  dropped on a distribution, it is shown over the period too
//
// Empty lane space pans and zooms the timeline like the plots do.

#include <timeline_periods.h>
#include <core/md_bitset.h>

// The timeline's time, shared with its plots: the clock (clock.h), in its unit, and the frames of its reference's run on it
struct TimelinePeriodAxis {
    const double* frame_times = nullptr;    // the reference's frames' times on the clock, never decreasing
    int           num_frames = 0;
    double        span_beg = 0.0;           // the clock's span: where a fixed period stays within
    double        span_end = 0.0;
    double*       view_beg = nullptr;       // the visible span of time, which the lane pans and zooms
    double*       view_end = nullptr;
    const char*   reference_name = "";      // of the reference's run, "catalyst"
    const char*   time_unit = "";           // the display unit, "ns"; empty when the frames have no time of their own
    double        display_scale = 1.0;      // from the clock's unit into time_unit
};

// Call where the lane goes, before the plots of the timeline. Ctrl-drag in it moves current_time, to the reference's
// frames.
void timeline_periods_lane(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double* current_time);

// Call inside each plot of the timeline, after its setup: shades the periods in the plot, and starts
// a new period on shift-press. The first plot each frame tells the lane where the plots' area is.
void timeline_periods_plot(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double current_time);

// Call inside each plot of the timeline, after timeline_periods_plot: the frames of a selection over frames (frames, as
// md_bitmask_t members) shaded faintly over the plot, with a band of their colour along its bottom edge, and those of a
// highlight over that, brighter. Runs of frames closer than a pixel are drawn as one. Either may be NULL.
void timeline_frames_plot(const md_bitmask_t* selected, const md_bitmask_t* highlighted, const TimelinePeriodAxis& axis, uint32_t selected_color, uint32_t highlight_color);

// The contents of the timeline's Periods menu: the layout of the lanes, the shading
void timeline_periods_menu(TimelinePeriods* rs);

// While a period is created, moved or resized, or the lane pans: the plots should not pan meanwhile
bool timeline_periods_dragging();

// Where the plots' area was across the screen, as the first plot of the last frame drew it (timeline_periods_plot):
// what another lane above the plots aligns with. False before the plots have been drawn.
bool timeline_plot_extent(float* x0, float* x1);

// ImGui drag and drop payload of a period pulled out of the lane: its id, a uint32_t
#define TIMELINE_PERIOD_DND "TIMELINE_PERIOD"

// The period the pointer is on in the lane, or holds (moving it, resizing it, or taking it elsewhere): what other views
// may show a preview of. 0 for none.
uint32_t timeline_periods_hovered();

// Lights the period up in the lane and the plots, as if it were hovered, for this frame and the next: a view showing
// something over it says which period it is
void timeline_periods_highlight(uint32_t id);
