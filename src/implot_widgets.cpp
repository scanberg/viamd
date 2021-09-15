#include "implot_widgets.h"
#include <implot_internal.h>

namespace ImPlot {

IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, bool show_label, const ImVec4& line_col, const ImVec4& fill_col, float line_thickness) {
    ImPlotContext& gp = *GImPlot;
    IM_ASSERT_USER_ERROR(gp.CurrentPlot != NULL, "DragRangeX() needs to be called between BeginPlot() and EndPlot()!");
    const float grab_size = ImMax(5.0f, line_thickness);
    float yt = gp.CurrentPlot->PlotRect.Min.y;
    float yb = gp.CurrentPlot->PlotRect.Max.y;
    float x_min = IM_ROUND(PlotToPixels(*x_range_min,0).x);
    float x_max = IM_ROUND(PlotToPixels(*x_range_max,0).x);
    const bool outside = x_max < (gp.CurrentPlot->PlotRect.Min.x - grab_size / 2) || x_min > (gp.CurrentPlot->PlotRect.Max.x + grab_size / 2);
    if (outside)
        return false;

    int dragging = 0;
    double* x_range[2] = {x_range_min, x_range_max};
    for (int i = 0; i < 2; ++i) {
        ImGui::PushID(i);
        char buff[64];
        snprintf(buff, sizeof(buff), "%s %s", id, i == 0 ? "min" : "max");
        if (DragLineX(buff, x_range[i], show_label, line_col, line_thickness)) dragging = i + 1;
        ImGui::PopID();
    }

    //float len = gp.Style.MajorTickLen.x;
    const float drag_bar_height = 15;

    bool bar_active = false;
    bool bar_hovered = false;
    if (!gp.CurrentPlot->Selecting && !gp.CurrentPlot->Querying) {
        if (dragging) {
            if (dragging == 1) {
                if (*x_range_min > *x_range_max - 1) {
                    *x_range_min = *x_range_max - 1;
                }
            } else {
                if (*x_range_max < *x_range_min + 1) {
                    *x_range_max = *x_range_min + 1;
                }
            }
        }
        else {
            ImVec2 old_cursor_pos = ImGui::GetCursorScreenPos();
            ImVec2 new_cursor_pos = ImVec2(x_min, yb - drag_bar_height);
            ImGui::GetCurrentWindow()->DC.CursorPos = new_cursor_pos;
            ImGui::InvisibleButton(id, ImVec2(x_max - x_min, drag_bar_height));
            ImGui::GetCurrentWindow()->DC.CursorPos = old_cursor_pos;
            bar_hovered = ImGui::IsItemHovered();
            bar_active = ImGui::IsItemActive();
            if (ImGui::IsItemActive() && ImGui::IsMouseDragging(0)) {
                const float delta = ImGui::GetIO().MouseDelta.x;
                double x_min_temp = PlotToPixels(*x_range_min, 0).x;
                double x_max_temp = PlotToPixels(*x_range_max, 0).x;

                x_min_temp += delta;
                x_max_temp += delta;

                *x_range_min = PixelsToPlot((float)x_min_temp, 0).x;
                *x_range_max = PixelsToPlot((float)x_max_temp, 0).x;

                double x_diff = *x_range_max - *x_range_min;

                if (*x_range_min < gp.CurrentPlot->XAxis.Range.Min) {
                    *x_range_min = gp.CurrentPlot->XAxis.Range.Min;
                    *x_range_max = *x_range_min + x_diff;
                }
                if (*x_range_max > gp.CurrentPlot->XAxis.Range.Max) {
                    *x_range_max = gp.CurrentPlot->XAxis.Range.Max;
                    *x_range_min = *x_range_max - x_diff;
                }

                dragging = 3;
            }
        }
    }

    ImVec4 range_color = IsColorAuto(fill_col) ? ImGui::GetStyleColorVec4(ImGuiCol_ScrollbarBg) : fill_col;
    ImVec4 bar_color = ImGui::GetStyleColorVec4(ImGuiCol_ScrollbarGrab);

    ImDrawList& DrawList = *GetPlotDrawList();
    PushPlotClipRect();
    DrawList.AddRectFilled(ImVec2(x_min+1, yt), ImVec2(x_max-1, yb), ImGui::ColorConvertFloat4ToU32(range_color));
    DrawList.AddRectFilled(ImVec2(x_min+1, yb), ImVec2(x_max-1, yb-drag_bar_height), ImGui::ColorConvertFloat4ToU32(bar_color));
    PopPlotClipRect();

    if (bar_hovered || bar_active) {
        gp.CurrentPlot->PlotHovered = false;
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        if (show_label) {
            char min_buff[32];
            char max_buff[32];
            char buff[64];
            float x  = IM_ROUND(PlotToPixels((*x_range_max + *x_range_min) * 0.5, 0).x);
            LabelAxisValue(gp.CurrentPlot->XAxis, gp.XTicks, *x_range_min, min_buff, sizeof(min_buff));
            LabelAxisValue(gp.CurrentPlot->XAxis, gp.XTicks, *x_range_max, max_buff, sizeof(max_buff));
            snprintf(buff, sizeof(buff), "[%s, %s]", min_buff, max_buff);
            ImVec4 color = IsColorAuto(line_col) ? ImGui::GetStyleColorVec4(ImGuiCol_Text) : line_col;
            ImU32  col32 = ImGui::ColorConvertFloat4ToU32(color);
            gp.Annotations.Append(ImVec2(x,yb),ImVec2(0,0),col32,CalcTextColor(color),true,"%s = %s", id, buff);
        }
    }

    return dragging;
}

}  // namespace ImGui
