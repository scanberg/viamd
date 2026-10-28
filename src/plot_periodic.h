#pragma once

// PERIODIC AXES
//
// Periodic values (an angle about an axis, see SeriesPeriodic in plot_series.h) lie on a circle, and a plot lays that circle
// out flat by cutting it open somewhere. Where it is cut is a choice: a dihedral may be shown in [-180, 180[, in [0, 360[,
// or cut wherever its values are least likely to be, so that they cross the cut as seldom as possible. Every value which
// crosses it is drawn leaving the plot at one edge and coming back in at the other.
//
// An axis of periodic values shows at most one period, and panning it turns the circle under the cut: the values wrap
// around as it goes, and the cut is wherever the axis begins. Its ticks are labelled as the values are given (200 degrees
// of a dihedral as -160), whatever the cut.
//
// One PeriodicView is the cut of every axis of one period: all the subplots, in the Timelines and in the Distributions
// alike, which show values of the same period in the same unit are linked to it, and so cut alike and read together.
// Values of another period, or of the same length given in another one (a pseudorotation phase in [0, 360[ beside a
// dihedral in [-180, 180[), have a view of their own, cut where suits them.

#include <plot_series.h>
#include <plot_periodic_core.h>

#ifndef IMGUI_DEFINE_MATH_OPERATORS
#define IMGUI_DEFINE_MATH_OPERATORS
#endif
#include <imgui.h>
#include <implot.h>

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

struct PeriodicView {
    SeriesPeriodic period = {};     // what it is the view of, and the seam of the series it was made for: display unit
    char   unit[32] = "";           // the display unit
    double min = 0, max = 0;        // what the axes show: the cut is at min, and max - min is at most a period
    double pending_min = 0;         // a cut asked for, applied before the axes are next set up
    bool   has_pending = false;
    bool   used = false;
};

// The view of values of period p in unit, made the first time it is asked for: cut at the seam of p, a period wide.
// Never NULL for a periodic p; the pointer stays valid for as long as the application runs.
PeriodicView* periodic_view(const SeriesPeriodic& p, const char* unit);

// Every view forgotten: what a workspace is read over
void periodic_views_clear();

static inline double periodic_view_length(const PeriodicView& v) {
    return series_periodic_length(v.period);
}

// Cuts the circle at cut: the axes show [cut, cut + their span[ from the next frame on
void periodic_view_cut(PeriodicView* view, double cut);

// Cuts the circle where it is least populated (periodic_least_populated): density is num_bins bins over the period of the
// view, as series_resolve_histogram bins periodic values
void periodic_view_cut_least_populated(PeriodicView* view, const float* density, int num_bins);

// Sets up axis of the plot being set up as an axis of the view: linked to it, at most a period wide, its ticks labelled in the
// period the values are given in, and without ImPlot's own menu (periodic_axis_popup takes its place). Between BeginPlot
// and SetupFinish, in place of ImPlot::SetupAxis.
void periodic_axis_setup(ImAxis axis, const char* label, PeriodicView* view, ImPlotAxisFlags flags);

// The menu of a periodic axis, opened by right-clicking it: where to cut the circle. Between SetupFinish and EndPlot.
// True when the cut is asked to go where the values are least populated, which only the caller can find out: it is
// expected to call periodic_view_cut_least_populated.
bool periodic_axis_popup(ImAxis axis, PeriodicView* view);

// Workspace: the cut of every view, so that a workspace opens cut as it was saved
void periodic_views_serialize(viamd::serialization_state_t& state);
// Reads the section the state is at if it is a view's. False otherwise.
bool periodic_views_deserialize(viamd::deserialization_state_t& state);
