#include <imgui.h>

namespace ImGui {

enum LinePlotFlags_ { LinePlotFlags_AxisX = 1 << 0, LinePlotFlags_AxisY = 1 << 1, LinePlotFlags_ShowXVal = 1 << 2, LinePlotFlats_Grid = 1 << 3 };

typedef int LinePlotFlags;

IMGUI_API void BeginPlot(const char* label, ImVec2 frame_size, ImVec2 x_range, ImVec2 y_range, LinePlotFlags flags = LinePlotFlags_ShowXVal);
IMGUI_API void EndPlot();

IMGUI_API void PlotVerticalBars(const float* bar_opacity, int count, ImU32 color);
IMGUI_API void PlotVariance(const float* avg, const float* var, int count, float var_scl = 1.f, ImU32 line_color = ImColor(1.f, 1.f, 0.2f, 0.3f),
                            ImU32 fill_color = ImColor(1.f, 1.f, 0.2f, 0.1f));
IMGUI_API void PlotValues(const char* line_label, const float* value, int count, ImU32 line_color = 0xffffffff);


IMGUI_API void DrawHistogram(ImVec2 frame_min, ImVec2 frame_max, const float* value, const int count, float max_val = 0.f,
                             ImU32 color = GetColorU32(ImGuiCol_PlotHistogram));
IMGUI_API void DrawFilledLine(ImVec2 frame_min, ImVec2 frame_max, const float* value, const int count, float max_val = 0.f,
                              ImU32 line_color = GetColorU32(ImGuiCol_PlotLines), ImU32 fill_color = GetColorU32(ImGuiCol_PlotHistogram));

IMGUI_API bool PlotHistogram(const char* label, ImVec2 frame_size, const float* value, int count, bool periodic = false,
                             ImVec2 value_range = ImVec2(0, 1), ImVec2* selection_range = nullptr);
IMGUI_API bool PlotPeriodic(const char* label, float outer_radius, float inner_radius_ratio, const float* value, int count, ImVec2 value_range,
                            ImU32 line_color = 0xffffffff);

}  // namespace ImGui
