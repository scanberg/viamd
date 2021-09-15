#include <imgui.h>

namespace ImGui {

IMGUI_API bool RangeSliderFloat(const char* label, float* v1, float* v2, float v_min, float v_max, const char* display_format = "(%.3f, %.3f)", ImGuiSliderFlags flags = 0);

}  // namespace ImGui
