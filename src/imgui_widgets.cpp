#include "imgui_widgets.h"
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

namespace ImGui {

extern template IMGUI_API float RoundScalarWithFormatT<float, float>(const char* format, ImGuiDataType data_type, float v);

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
            *values[i] = RoundScalarWithFormatT<float, float>(format, ImGuiDataType_Float, new_value);
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
            *values[i] = RoundScalarWithFormatT<float, float>(format, ImGuiDataType_Float, new_value);
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

}  // namespace ImGui
