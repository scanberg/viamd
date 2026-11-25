#include <imgui.h>

#include <core/md_vec.h>

struct ImRect;
struct ImGuiWindow;

namespace ImGui {

IMGUI_API bool RangeSliderFloat(const char* label, float* v1, float* v2, float v_min, float v_max, const char* display_format = "(%.3f, %.3f)", ImGuiSliderFlags flags = 0);

// custom ImGui procedures
bool DeleteButton(const char* label, const ImVec2& size = ImVec2(0, 0));
// The dockspace over the main viewport, which windows dock into; its central node, which they leave free, is returned
// (screen coordinates): where the main view is (BeginCanvas)
ImRect CreateDockspace();
// An invisible window over rect, behind every other, which takes the mouse where no window is over it
void BeginCanvas(const char* id, const ImRect& rect, bool allow_inputs = false);
void EndCanvas();

inline void PushInvalid() {
	ImGui::PushStyleColor(ImGuiCol_FrameBg, 0xAA222299);
}

inline void PopInvalid() {
	ImGui::PopStyleColor();
}


inline bool InputQuery(const char* label, char* buf, size_t buf_size, bool is_valid, const char* err_text = 0, ImGuiInputTextFlags flags = 0) {
	if (!is_valid) ImGui::PushInvalid();
	bool result = ImGui::InputText(label, buf, buf_size, flags);
	if (!is_valid) ImGui::PopInvalid();
    size_t err_len = err_text ? strlen(err_text) : 0;
	if (ImGui::IsItemHovered() && !is_valid && err_len > 0) {
		ImGui::SetTooltip("%s", err_text);
	}
	return result;
}



void init_theme();

void PushDisabled();
void PopDisabled();

bool ColorEdit3Minimal(const char* label, float color[3]);
bool ColorEdit4Minimal(const char* label, float color[4]);

void DrawCheckerboard(ImDrawList* draw_list, ImVec2 p_min, ImVec2 p_max, ImU32 col1, ImU32 col2, float grid_step, ImVec2 grid_off);

bool CoordinateSystemWidget(md_quat_t* out_orientation, const md_quat_t& in_orientation, const ImVec2& size);

}  // namespace ImGui
