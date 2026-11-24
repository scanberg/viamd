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
#include <core/md_bitset.h>

// How frames map to the screen, shared with the timeline's plots: in the run's time, in its own unit (frame_axis.h)
struct TimelinePeriodAxis {
    const double* frame_times = nullptr;    // per frame, never decreasing
    int           num_frames = 0;
    double*       view_beg = nullptr;       // the visible span of time, which the lane pans and zooms
    double*       view_end = nullptr;
    const char*   time_unit = "";           // the display unit, "ns"; empty when the frames have no time of their own
    double        display_scale = 1.0;      // from the unit of frame_times into time_unit
};

// Call where the lane goes, before the plots of the timeline. Ctrl-drag in it moves current_frame.
void timeline_periods_lane(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double* current_frame);

// Call inside each plot of the timeline, after its setup: shades the periods in the plot, and starts
// a new period on shift-press. The first plot each frame tells the lane where the plots' area is.
void timeline_periods_plot(TimelinePeriods* rs, const TimelinePeriodAxis& axis, double current_frame);

// Call inside each plot of the timeline, after timeline_periods_plot: the frames of a selection over frames (frames, as
// md_bitmask_t members) shaded faintly over the plot, with a band of their colour along its bottom edge, and those of a
// highlight over that, brighter. Runs of frames closer than a pixel are drawn as one. Either may be NULL.
void timeline_frames_plot(const md_bitmask_t* selected, const md_bitmask_t* highlighted, const TimelinePeriodAxis& axis, uint32_t selected_color, uint32_t highlight_color);

// Call inside a plot of the timeline, between ImPlot's BeginItem and EndItem of the series: a series of truths over the
// frames (values[frame * stride + k], true above 0.5: a script's bools) as the frames it is true in, shaded over the plot
// as timeline_frames_plot shades a selection, in color, with a band of it along the top edge: the lane-th one down, so
// that the truths of a plot are each seen apart. Brighter when hovered.
void timeline_truth_plot(const float* values, int stride, int k, const TimelinePeriodAxis& axis, uint32_t color, int lane, bool hovered);
// Where the band of lane is, in pixels [*y0, *y1): false when the plot has no room for it
bool timeline_truth_band(int lane, float* y0, float* y1);

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
