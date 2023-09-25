#include <imgui_widgets.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

namespace ImGui {

extern template IMGUI_API float RoundScalarWithFormatT<float>(const char* format, ImGuiDataType data_type, float v);
extern template IMGUI_API float ScaleRatioFromValueT<float, float, float>(ImGuiDataType data_type, float t, float v_min, float v_max, bool is_logarithmic, float logarithmic_zero_epsilon, float zero_deadzone_halfsize);
extern template IMGUI_API float ScaleValueFromRatioT<float, float, float>(ImGuiDataType data_type, float t, float v_min, float v_max, bool is_logarithmic, float logarithmic_zero_epsilon, float zero_deadzone_halfsize);

// https://github.com/ocornut/imgui/issues/76

bool RangeSliderBehavior(const ImRect& frame_bb, const char* str_id, float* v1, float* v2, float v_min, float v_max, const char* format, ImGuiSliderFlags flags) {
    ImGuiContext& g = *GImGui;
    ImGuiWindow* window = GetCurrentWindow();
    const ImGuiStyle& style = g.Style;

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    const bool is_logarithmic = (flags & ImGuiSliderFlags_Logarithmic);
    const bool is_horizontal = (flags & ImGuiSliderFlags_Vertical) == 0;

    const float grab_padding = 2.0f;
    const float slider_sz = is_horizontal ? (frame_bb.GetWidth() - grab_padding * 2.0f) : (frame_bb.GetHeight() - grab_padding * 2.0f);
    const float grab_sz = ImMin(style.GrabMinSize, slider_sz);
    const float slider_usable_sz = slider_sz - grab_sz;
    const float slider_usable_pos_min = (is_horizontal ? frame_bb.Min.x : frame_bb.Min.y) + grab_padding + grab_sz * 0.5f;
    const float slider_usable_pos_max = (is_horizontal ? frame_bb.Max.x : frame_bb.Max.y) - grab_padding - grab_sz * 0.5f;

    bool changed = false;

    float* values[2] = {v1, v2};

    for (int i = 0; i < 2; ++i) {
        float t = ImClamp(ScaleRatioFromValueT<float, float, float>(ImGuiDataType_Float, *values[i], v_min, v_max, is_logarithmic, 0.0f, 0.0f), 0.0f, 1.0f);
        if (!is_horizontal) t = 1.0f - t;
        float grab_pos = ImLerp(slider_usable_pos_min, slider_usable_pos_max, t);
        ImRect grab_bb;
        if (is_horizontal)
            grab_bb = ImRect(ImVec2(grab_pos - grab_sz * 0.5f, frame_bb.Min.y + grab_padding),
                ImVec2(grab_pos + grab_sz * 0.5f, frame_bb.Max.y - grab_padding));
        else
            grab_bb = ImRect(ImVec2(frame_bb.Min.x + grab_padding, grab_pos - grab_sz * 0.5f),
                ImVec2(frame_bb.Max.x - grab_padding, grab_pos + grab_sz * 0.5f));


        SetCursorScreenPos(grab_bb.Min);

        ImGui::PushID(i);
        ImGui::InvisibleButton(str_id, grab_bb.GetSize());
        bool active = IsItemActive();
        bool hovered = IsItemHovered();

        if (active || hovered)
        {
            ImGui::SetTooltip("%f", float(*values[i]));
            hovered = true;
        }
        if (active && IsMouseDragging(0))
        {
            const float mouse_abs_pos = is_horizontal ? g.IO.MousePos.x : g.IO.MousePos.y;
            float mouse_t = (slider_usable_sz > 0.0f) ? ImClamp((mouse_abs_pos - slider_usable_pos_min) / slider_usable_sz, 0.0f, 1.0f) : 0.0f;
            float new_value = ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, mouse_t, v_min, v_max, is_logarithmic, 0.0f, 0.0f);
            *values[i] = RoundScalarWithFormatT<float>(format, ImGuiDataType_Float, new_value);
            changed = true;
        }

        ImGui::PopID();

        // Draw
        ImGuiCol col = active ? ImGuiCol_ButtonActive : (hovered ? ImGuiCol_ButtonHovered : ImGuiCol_Button);
        window->DrawList->AddRectFilled(grab_bb.Min, grab_bb.Max, GetColorU32(col), style.GrabRounding);
    }


    float grab_pos[2];
    for (int i = 0; i < 2; ++i) {
        float t = ImClamp(ScaleRatioFromValueT<float, float, float>(ImGuiDataType_Float, *values[i], v_min, v_max, is_logarithmic, 0.0f, 0.0f), 0.0f, 1.0f);
        if (!is_horizontal) t = 1.0f - t;
        grab_pos[i] = ImLerp(slider_usable_pos_min, slider_usable_pos_max, t);
    }

    ImRect grab_bb;
    if (is_horizontal)
        grab_bb = ImRect(ImVec2(grab_pos[0] + grab_sz * 0.5f, frame_bb.Min.y + grab_padding),
                         ImVec2(grab_pos[1] - grab_sz * 0.5f, frame_bb.Max.y - grab_padding));
    else
        grab_bb = ImRect(ImVec2(frame_bb.Min.x + grab_padding, grab_pos[0] + grab_sz * 0.5f),
                         ImVec2(frame_bb.Max.x - grab_padding, grab_pos[1] - grab_sz * 0.5f));

    SetCursorScreenPos(grab_bb.Min);

    ImGui::PushID(-1);
    ImGui::InvisibleButton(str_id, grab_bb.GetSize());
    bool active = IsItemActive();
    bool hovered = IsItemHovered();

    if (hovered) {
        const float wheel_delta = GetIO().MouseWheel;
        if (wheel_delta) {
            float t[2] = {
                ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, *values[0], v_min, v_max, is_logarithmic, 0.0f, 0.0f),
                ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, *values[1], v_min, v_max, is_logarithmic, 0.0f, 0.0f),
            };

            const float delta = wheel_delta * (v_max - v_min) * 0.025f;

            t[0] -= delta;
            t[1] += delta;

            for (int i = 0; i < 2; ++i) {
                float new_value = ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, t[i], v_min, v_max, is_logarithmic, 0.0f, 0.0f);
                *values[i] = RoundScalarWithFormatT<float>(format, ImGuiDataType_Float, new_value);
            }
            changed = true;
        }
    }

    if (active && IsMouseDragging(0))
    {
        // This is a bit cumbersome, but we want to make sure we don't squash the range when we go past min or max with the range

        float t[2] = {
            ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, *values[0], v_min, v_max, is_logarithmic, 0.0f, 0.0f),
            ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, *values[1], v_min, v_max, is_logarithmic, 0.0f, 0.0f),
        };

        const float t_diff = t[1] - t[0];

        const float mouse_delta_t = (is_horizontal ? GetIO().MouseDelta.x  : GetIO().MouseDelta.y) / slider_sz;
        t[0] += mouse_delta_t;
        t[1] += mouse_delta_t;

        if (t[0] < 0.0f) {
            t[0] = 0.0f;
            t[1] = t[0] + t_diff;
        }
        if (t[1] > 1.0f) {
            t[1] = 1.0f;
            t[0] = t[1] - t_diff;
        }

        for (int i = 0; i < 2; ++i) {
            float new_value = ScaleValueFromRatioT<float, float, float>(ImGuiDataType_Float, t[i], v_min, v_max, is_logarithmic, 0.0f, 0.0f);
            *values[i] = RoundScalarWithFormatT<float>(format, ImGuiDataType_Float, new_value);
        }
        
        changed = true;
    }

    ImGui::PopID();

    // Draw
    ImGuiCol col = active ? ImGuiCol_SliderGrabActive : ImGuiCol_SliderGrab;
    window->DrawList->AddRectFilled(grab_bb.Min, grab_bb.Max, GetColorU32(col), style.GrabRounding);

    return changed;
}

// ~95% common code with ImGui::SliderFloat
bool RangeSliderFloat(const char* label, float* v1, float* v2, float v_min, float v_max, const char* display_format, ImGuiSliderFlags flags) {
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems) return false;

    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const float w = CalcItemWidth();

    const ImVec2 label_size = CalcTextSize(label, NULL, true);
    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(w, label_size.y + style.FramePadding.y * 2.0f));
    const ImRect total_bb(frame_bb.Min, frame_bb.Max + ImVec2(label_size.x > 0.0f ? style.ItemInnerSpacing.x + label_size.x : 0.0f, 0.0f));

    ImGui::BeginChild(label, total_bb.GetSize());

    // Actual slider behavior + render grab
    bool value_changed = RangeSliderBehavior(frame_bb, label, v1, v2, v_min, v_max, display_format, flags);

    // Display value using user-provided display format so user can add prefix/suffix/decorations to
    // the value.
    char value_buf[64];
    const char* value_buf_end = value_buf + ImFormatString(value_buf, IM_ARRAYSIZE(value_buf), display_format, *v1, *v2);
    RenderTextClipped(frame_bb.Min, frame_bb.Max, value_buf, value_buf_end, NULL, ImVec2(0.5f, 0.5f));

    if (label_size.x > 0.0f) RenderText(ImVec2(frame_bb.Max.x + style.ItemInnerSpacing.x, frame_bb.Min.y + style.FramePadding.y), label);

    ImGui::EndChild();

    return value_changed;
}

void init_theme() {
    ImVec4* colors = ImGui::GetStyle().Colors;
    colors[ImGuiCol_Text]                   = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
    colors[ImGuiCol_WindowBg]               = ImVec4(0.00f, 0.00f, 0.00f, 0.67f);
    colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_PopupBg]                = ImVec4(0.11f, 0.11f, 0.14f, 0.92f);
    colors[ImGuiCol_Border]                 = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
    colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_FrameBg]                = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
    colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.70f, 0.70f, 0.70f, 0.40f);
    colors[ImGuiCol_FrameBgActive]          = ImVec4(0.65f, 0.65f, 0.65f, 0.69f);
    colors[ImGuiCol_TitleBg]                = ImVec4(0.07f, 0.07f, 0.07f, 0.83f);
    colors[ImGuiCol_TitleBgActive]          = ImVec4(0.00f, 0.00f, 0.00f, 0.87f);
    colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.00f, 0.00f, 0.00f, 0.20f);
    colors[ImGuiCol_MenuBarBg]              = ImVec4(0.55f, 0.55f, 0.55f, 0.80f);
    colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.21f, 0.21f, 0.21f, 0.60f);
    colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.83f, 0.83f, 0.83f, 0.30f);
    colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.83f, 0.83f, 0.83f, 0.40f);
    colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.81f, 0.81f, 0.81f, 0.60f);
    colors[ImGuiCol_CheckMark]              = ImVec4(0.90f, 0.90f, 0.90f, 0.50f);
    colors[ImGuiCol_SliderGrab]             = ImVec4(1.00f, 1.00f, 1.00f, 0.30f);
    colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.82f, 0.82f, 0.82f, 0.60f);
    colors[ImGuiCol_Button]                 = ImVec4(0.64f, 0.64f, 0.64f, 0.62f);
    colors[ImGuiCol_ButtonHovered]          = ImVec4(0.72f, 0.72f, 0.72f, 0.79f);
    colors[ImGuiCol_ButtonActive]           = ImVec4(0.80f, 0.80f, 0.81f, 0.85f);
    colors[ImGuiCol_Header]                 = ImVec4(0.64f, 0.64f, 0.64f, 0.45f);
    colors[ImGuiCol_HeaderHovered]          = ImVec4(0.65f, 0.65f, 0.65f, 0.80f);
    colors[ImGuiCol_HeaderActive]           = ImVec4(0.85f, 0.85f, 0.85f, 0.80f);
    colors[ImGuiCol_Separator]              = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
    colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.71f, 0.71f, 0.71f, 1.00f);
    colors[ImGuiCol_SeparatorActive]        = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    colors[ImGuiCol_ResizeGrip]             = ImVec4(1.00f, 1.00f, 1.00f, 0.16f);
    colors[ImGuiCol_ResizeGripHovered]      = ImVec4(1.00f, 1.00f, 1.00f, 0.47f);
    colors[ImGuiCol_ResizeGripActive]       = ImVec4(1.00f, 0.96f, 1.00f, 0.63f);
    colors[ImGuiCol_Tab]                    = ImVec4(0.56f, 0.56f, 0.56f, 0.78f);
    colors[ImGuiCol_TabHovered]             = ImVec4(0.87f, 0.87f, 0.87f, 0.80f);
    colors[ImGuiCol_TabActive]              = ImVec4(0.73f, 0.73f, 0.73f, 0.84f);
    colors[ImGuiCol_TabUnfocused]           = ImVec4(0.57f, 0.57f, 0.57f, 0.82f);
    colors[ImGuiCol_TabUnfocusedActive]     = ImVec4(0.65f, 0.65f, 0.65f, 0.84f);
    colors[ImGuiCol_DockingPreview]         = ImVec4(0.90f, 0.90f, 0.90f, 0.31f);
    colors[ImGuiCol_DockingEmptyBg]         = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    colors[ImGuiCol_PlotLines]              = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_PlotLinesHovered]       = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.00f, 0.00f, 1.00f, 0.35f);
    colors[ImGuiCol_DragDropTarget]         = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
    colors[ImGuiCol_NavHighlight]           = ImVec4(0.45f, 0.45f, 0.90f, 0.80f);
    colors[ImGuiCol_NavWindowingHighlight]  = ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    colors[ImGuiCol_NavWindowingDimBg]      = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);
    colors[ImGuiCol_ModalWindowDimBg]       = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
}

void CreateDockspace() {
    // Invisible dockspace
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos);
    ImGui::SetNextWindowSize(viewport->Size);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::SetNextWindowBgAlpha(0.0f);

    const ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_NoNavFocus;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin("DockspaceWindow", NULL, window_flags);
    ImGui::PopStyleVar(3);

    // ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_PassthruDockspace;
    const ImGuiID id = ImGui::GetID("Dockspace");
    const ImGuiDockNodeFlags flags = ImGuiDockNodeFlags_PassthruCentralNode;
    ImGui::DockSpace(id, ImVec2(0.0f, 0.0f), flags);

    ImGui::End();
}

void BeginCanvas(const char* id, bool allow_inputs) {
    // Invisible Canvas
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos);
    ImGui::SetNextWindowSize(viewport->Size);
    ImGui::SetNextWindowBgAlpha(0.0f);
    ImGui::SetNextWindowViewport(viewport->ID);

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoNavInputs | ImGuiWindowFlags_NoDocking;

    if (!allow_inputs) {
        window_flags |= ImGuiWindowFlags_NoInputs;
    }

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin(id, NULL, window_flags);
    ImGui::PopStyleVar(3);
}

void EndCanvas() { ImGui::End(); }

constexpr unsigned int DEL_BTN_COLOR = 0xFF1111CC;
constexpr unsigned int DEL_BTN_HOVER_COLOR = 0xFF3333DD;
constexpr unsigned int DEL_BTN_ACTIVE_COLOR = 0xFF5555FF;

bool DeleteButton(const char* label, const ImVec2& size) {
    PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
    PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
    PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
    bool res = ImGui::Button(label, size);
    PopStyleColor(3);
    return res;
}

void PushDisabled() {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
}
void PopDisabled() {
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
}

bool ColorEdit3Minimal(const char* label, float color[3]) {
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
    bool result = ImGui::ColorEdit3(label, color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_Float);
    ImGui::PopStyleVar();
    return result;
}

bool ColorEdit4Minimal(const char* label, float color[4]) {
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
    bool result = ImGui::ColorEdit4(label, color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_Float);
    ImGui::PopStyleVar();
    return result;
}

void DrawCheckerboard(ImDrawList* draw_list, ImVec2 p_min, ImVec2 p_max, ImU32 col1, ImU32 col2, float grid_step, ImVec2 grid_off) {
    draw_list->AddRectFilled(p_min, p_max, col1);

    int yi = 0;
    for (float y = p_min.y + grid_off.y; y < p_max.y; y += grid_step, yi++)
    {
        float y1 = ImClamp(y, p_min.y, p_max.y), y2 = ImMin(y + grid_step, p_max.y);
        if (y2 <= y1)
            continue;
        for (float x = p_min.x + grid_off.x + (yi & 1) * grid_step; x < p_max.x; x += grid_step * 2.0f)
        {
            float x1 = ImClamp(x, p_min.x, p_max.x), x2 = ImMin(x + grid_step, p_max.x);
            if (x2 <= x1)
                continue;
            draw_list->AddRectFilled(ImVec2(x1, y1), ImVec2(x2, y2), col2);
        }
    }
}

bool SliderDouble(const char* label, double* v, double v_min, double v_max, const char* format, ImGuiSliderFlags flags) {
    return SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max, format, flags);
}

// NB: You likely want to specify the ImGuiSliderFlags_AlwaysClamp when using this.
bool ImGui::DragDoubleRange2(const char* label, double* v_current_min, double* v_current_max, double v_speed, double v_min, double v_max, const char* format, const char* format_max, ImGuiSliderFlags flags)
{
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems)
        return false;

    ImGuiContext& g = *GImGui;
    PushID(label);
    BeginGroup();
    PushMultiItemsWidths(2, CalcItemWidth());

    double min_min = (v_min >= v_max) ? -DBL_MAX : v_min;
    double min_max = (v_min >= v_max) ? *v_current_max : ImMin(v_max, *v_current_max);
    ImGuiSliderFlags min_flags = flags | ((min_min == min_max) ? ImGuiSliderFlags_ReadOnly : 0);
    bool value_changed = DragScalar("##min", ImGuiDataType_Double, v_current_min, (float)v_speed, &min_min, &min_max, format, min_flags);
    PopItemWidth();
    SameLine(0, g.Style.ItemInnerSpacing.x);

    double max_min = (v_min >= v_max) ? *v_current_min : ImMax(v_min, *v_current_min);
    double max_max = (v_min >= v_max) ? DBL_MAX : v_max;
    ImGuiSliderFlags max_flags = flags | ((max_min == max_max) ? ImGuiSliderFlags_ReadOnly : 0);
    value_changed |= DragScalar("##max", ImGuiDataType_Double, v_current_max, (float)v_speed, &max_min, &max_max, format_max ? format_max : format, max_flags);
    PopItemWidth();
    SameLine(0, g.Style.ItemInnerSpacing.x);

    TextEx(label, FindRenderedTextEnd(label));
    EndGroup();
    PopID();

    return value_changed;
}

}  // namespace ImGui
