#include <implot.h>

namespace ImPlot {

typedef int ImPlotDragRangeFlags;          // -> enum ImPlotBin_

                                // Options for plots (see BeginPlot).
enum ImPlotDragRangeFlags_ {
    ImPlotDragRangeFlags_None = 0,          // default
    ImPlotDragRangeFlags_Periodic = 1 << 0, // the range is periodic and should wrap around the min and max values.
    ImPlotDragRangeFlags_NoLabel  = 1 << 1, // do not display any label
};

struct ImPlotDragRangeStyle {
    ImVec4 line_col = IMPLOT_AUTO_COL;
    ImVec4 range_col = IMPLOT_AUTO_COL;
    ImVec4 scrollbar_col = IMPLOT_AUTO_COL;
    float line_thickness = 1;
};

IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, double min_value, double max_value, ImPlotDragRangeFlags flags = 0, const ImPlotDragRangeStyle& style = ImPlotDragRangeStyle());

IMPLOT_API bool ColorMapSelection(const char* id, ImPlotColormap* idx, float* cur_range_min, float* cur_range_max, float min = 0, float max = 0, ImVec2 size = ImVec2(150,0));

}  // namespace ImPlot
