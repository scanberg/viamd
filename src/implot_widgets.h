#include <implot.h>

typedef int ImPlotDragRangeFlags;          // -> enum ImPlotBin_

// Options for plots (see BeginPlot).
enum ImPlotDragRangeFlags_ {
    ImPlotDragRangeFlags_None     = 0,  // default
    ImPlotDragRangeFlags_NoLabel  = 16, // No displayed label
    ImPlotDragRangeFlags_NoBar    = 32, // No draggable bar at bottom
};

// INSET AXES
// The tick labels and the names of axes drawn within the plot, along its edges, rather than in margins around it: the
// plot keeps all of its area. For plots stacked or side by side (subplots), where the margins take a large share.
//
// Axis by axis: an axis set up without ImPlot's tick labels (ImPlotAxisFlags_NoTickLabels) has them drawn within, and
// one set up without its label (ImPlotAxisFlags_NoLabel) has its name drawn within. InsetAxisFlags is both. An axis
// set up as usual is left to ImPlot, margin and all, so the two mix: only the y axis inset is
//
//   ImPlot::SetupAxes("time", "value", 0, ImPlot::InsetAxisFlags);    // ImPlot reserves no room for y's labels
//   ImPlot::SetupFinish();
//   ...                                                                // what is plotted
//   ImPlot::DrawInsetAxes();                                           // over it, before EndPlot
//
// An axis whose labels are to be hidden rather than drawn within (a plot under others sharing its x axis) is left out
// with ImPlotInsetFlags_NoX / _NoY.
//
// The labels are the ones ImPlot made of the ticks, default or custom, through the axis's formatter. An x axis has its
// row along the bottom (the top with ImPlotAxisFlags_Opposite) with its name at the right end; a y axis its column along
// the left (the right) with its name at the top. More axes on a side stack inwards. A label is left out where it would
// run into the one before it, a corner, the name of the axis, the legend (where it is within the plot), or the region
// clear_top_left (in pixels from the plot's top left corner, such as a title): zoomed in, they drop out rather than
// pile up.
//
// The row or column of an inset axis takes the place of the margin ImPlot would hover it in: IsInsetAxisHovered says
// whether the mouse is on it (as of the last DrawInsetAxes of the plot), for what an axis does when hovered (a menu).
typedef int ImPlotInsetFlags;

enum ImPlotInsetFlags_ {
    ImPlotInsetFlags_None         = 0,
    ImPlotInsetFlags_NoTickLabels = 1 << 0, // No tick labels within, of any axis
    ImPlotInsetFlags_NoAxisLabels = 1 << 1, // No names of axes within
    ImPlotInsetFlags_NoX          = 1 << 2, // Nothing of the x axes: their labels are hidden
    ImPlotInsetFlags_NoY          = 1 << 3, // Nothing of the y axes: their labels are hidden
};

struct ImPlotDragRangeStyle {
    ImVec4 line_col = IMPLOT_AUTO_COL;
    ImVec4 range_col = IMPLOT_AUTO_COL;
    ImVec4 scrollbar_col = IMPLOT_AUTO_COL;
    float line_thickness = 1;
};

namespace ImPlot {

IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, double min_value, double max_value, ImPlotDragRangeFlags flags = 0, const ImPlotDragRangeStyle& style = ImPlotDragRangeStyle());

// An axis whose labels DrawInsetAxes draws within the plot: it draws none of its own, and so takes no room for them
constexpr ImPlotAxisFlags InsetAxisFlags = ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels;

IMPLOT_API void DrawInsetAxes(ImPlotInsetFlags flags = 0, ImVec2 clear_top_left = ImVec2(0, 0), float font_scale = 0.8f);

// Whether the mouse is on the row or column of the inset axis of the current plot (see DrawInsetAxes)
IMPLOT_API bool IsInsetAxisHovered(ImAxis axis);

IMPLOT_API bool ColormapSelection(const char* id, ImPlotColormap* idx, ImVec2 size = ImVec2(0, 0));

IMPLOT_API bool ColormapQualitative(ImPlotColormap idx);

IMPLOT_API void SyncAxesY(double padding_fraction = 0.1);

IMPLOT_API void SyncAxesWithPadding(int master_axis, int aux_axis, double aux_to_master, double padding_fraction = 0.1);

IMPLOT_API void PlotLineSpline(const char* label_id, const double* xs, const double* ys, int count);

}  // namespace ImPlot
