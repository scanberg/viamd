#include <implot.h>

typedef int ImPlotDragRangeFlags;          // -> enum ImPlotBin_

// Options for plots (see BeginPlot).
enum ImPlotDragRangeFlags_ {
    ImPlotDragRangeFlags_None     = 0,  // default
    ImPlotDragRangeFlags_NoLabel  = 16, // No displayed label
    ImPlotDragRangeFlags_NoBar    = 32, // No draggable bar at bottom
};

// INSET AXES
// The tick labels and the names of the axes drawn within the plot, along its edges, rather than in margins around it:
// the plot keeps all of its area. For small plots side by side (subplots), where the margins take a large share.
//
//   ImPlot::SetupAxes("x", "y", ImPlot::InsetAxisFlags, ImPlot::InsetAxisFlags);   // ImPlot reserves no room for labels
//   ImPlot::SetupFinish();
//   ...                                                                         // what is plotted
//   ImPlot::DrawInsetAxes();                                                     // over it, before EndPlot
//
// The labels are the ones ImPlot made of the ticks, default or custom, through the axis's formatter. An x axis has its
// row along the bottom (the top with ImPlotAxisFlags_Opposite) with its name at the right end; a y axis its column
// along the left (the right) with its name at the top. More axes on a side stack inwards. Along an edge a label is left
// out where it would run into the one before it, a corner, or the name of the axis: zoomed in, they drop out rather
// than pile up. clear_top_left is a rectangle, in pixels from the plot's top left corner, they keep out of (a title).
typedef int ImPlotInsetFlags;

enum ImPlotInsetFlags_ {
    ImPlotInsetFlags_None         = 0,
    ImPlotInsetFlags_NoTickLabels = 1 << 0, // Only the names of the axes
    ImPlotInsetFlags_NoAxisLabels = 1 << 1, // Only the tick labels
};

struct ImPlotDragRangeStyle {
    ImVec4 line_col = IMPLOT_AUTO_COL;
    ImVec4 range_col = IMPLOT_AUTO_COL;
    ImVec4 scrollbar_col = IMPLOT_AUTO_COL;
    float line_thickness = 1;
};

namespace ImPlot {

IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, double min_value, double max_value, ImPlotDragRangeFlags flags = 0, const ImPlotDragRangeStyle& style = ImPlotDragRangeStyle());

// The axis flags for DrawInsetAxes: the axis draws no labels of its own, and so takes no room for them
constexpr ImPlotAxisFlags InsetAxisFlags = ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels;

IMPLOT_API void DrawInsetAxes(ImPlotInsetFlags flags = 0, ImVec2 clear_top_left = ImVec2(0, 0), float font_scale = 0.8f);

IMPLOT_API bool ColormapSelection(const char* id, ImPlotColormap* idx, ImVec2 size = ImVec2(0, 0));

IMPLOT_API bool ColormapQualitative(ImPlotColormap idx);

IMPLOT_API void SyncAxesY(double padding_fraction = 0.1);

IMPLOT_API void SyncAxesWithPadding(int master_axis, int aux_axis, double aux_to_master, double padding_fraction = 0.1);

IMPLOT_API void PlotLineSpline(const char* label_id, const double* xs, const double* ys, int count);

}  // namespace ImPlot
