#include <implot.h>

namespace ImPlot {
// Shows a draggable vertical guide range between two x-values. #col defaults to ImGuiCol_Text.
IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, bool show_label = true, const ImVec4& line_col = IMPLOT_AUTO_COL, const ImVec4& range_col = IMPLOT_AUTO_COL, float line_thickness = 1);

}  // namespace ImPlot
