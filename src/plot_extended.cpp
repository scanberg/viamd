#include "plot_extended.h"
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

namespace ImGui {

struct PlotState {
    ImGuiID id = 0;
    ImRect inner_bb;
    ImRect coord_view;
    // ImVec2* selection_range;
    // float selection_start;
    // bool is_selecting = false;
};

static PlotState ps;

const ImU32 HOVER_LINE_COLOR = 0xaaffffff;
const ImU32 CURRENT_LINE_COLOR = 0xaa33ffff;
const ImU32 AXIS_COLOR = 0xaaffffff;
const ImU32 SELECTION_RANGE_COLOR = 0x33bbbbbb;

IMGUI_API void BeginPlot(const char* label, ImVec2 frame_size, ImVec2 x_range, ImVec2 y_range, LinePlotFlags flags) {
    ImGuiWindow* window = GetCurrentWindow();
    const ImGuiID id = window->GetID(label);
    ps.id = id;

    if (window->SkipItems) return;

    ImGuiContext& ctx = *GImGui;
    const ImGuiStyle& style = ctx.Style;

    if (frame_size.x == 0.0f) frame_size.x = CalcItemWidth();
    if (frame_size.y == 0.0f) frame_size.y = (style.FramePadding.y * 2);

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, 0)) return;

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    if (flags & LinePlotFlags_AxisY) {
        if (x_range.x < 0 && 0 < x_range.y) {
            const float t = (0 - x_range.x) / (x_range.y - x_range.x);
            const float x = ImLerp(inner_bb.Min.x, inner_bb.Max.x, t);
            window->DrawList->AddLine(ImVec2(x, inner_bb.Min.y), ImVec2(x, inner_bb.Max.y), AXIS_COLOR);
        }
    }

    if (flags & LinePlotFlags_AxisX) {
        if (y_range.x < 0 && 0 < y_range.y) {
            const float t = (0 - y_range.x) / (y_range.y - y_range.x);
            const float y = roundf(ImLerp(inner_bb.Min.y, inner_bb.Max.y, 1.f - t));
            window->DrawList->AddLine(ImVec2(roundf(inner_bb.Min.x), y), ImVec2(roundf(inner_bb.Max.x), y), AXIS_COLOR);
        }
    }

    /*
if (selection_range) {
    selection_range->x = ImClamp(selection_range->x, x_range.x, x_range.y);
    selection_range->y = ImClamp(selection_range->y, x_range.x, x_range.y);

    // Draw selection range
    const float t0 = (selection_range->x - x_range.x) / (x_range.y - x_range.x);
    const float t1 = (selection_range->y - x_range.x) / (x_range.y - x_range.x);
    ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
    ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
    window->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
}
    */

    ps.inner_bb = inner_bb;
    ps.coord_view = ImRect(ImVec2(x_range.x, y_range.x), ImVec2(x_range.y, y_range.y));
    // ps.selection_range = selection_range;

    RenderTextClipped(ImVec2(frame_bb.Min.x + style.FramePadding.x, frame_bb.Min.y + style.FramePadding.y), frame_bb.Max, label, NULL, NULL,
                      ImVec2(0.0f, 0.0f));

    // return interacting_x_val;
}

static inline ImVec2 compute_frame_coord(ImRect frame, ImVec2 coord) {
    float cx = ImClamp((coord.x - frame.Min.x) / (frame.Max.x - frame.Min.x), 0.f, 1.f);
    float cy = ImClamp((coord.y - frame.Min.y) / (frame.Max.y - frame.Min.y), 0.f, 1.f);
    return {cx, cy};
}

IMGUI_API void PlotVerticalBars(const float* bar_opacity, int count, ImU32 color) {
    if (!ps.id) return;
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems) return;
    if (count < 2) return;

    ImVec4 base_color = ImColor(color);

    // ImVec2 prev_c = ImVec2(0, values[0]);
    // bool prev_v = valid_range.x <= values[0] && values[0] <= valid_range.y;
    for (int i = 0; i < count; i++) {
        if (bar_opacity[i] > 0) {
            ImVec2 prev_c = ImVec2((float)i, ps.coord_view.Min.y);
            float prev_o = bar_opacity[i];
            i++;
            while (i < count) {
                if (bar_opacity[i] != prev_o) break;
                i++;
            }
            ImVec2 next_c = ImVec2((float)i, ps.coord_view.Max.y);

            if (next_c.x < ps.coord_view.Min.x) continue;
            if (prev_c.x > ps.coord_view.Max.x) break;

            ImVec2 p = compute_frame_coord(ps.coord_view, prev_c);
            ImVec2 n = compute_frame_coord(ps.coord_view, next_c);

            ImVec2 pos0 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(p.x, 1.f - p.y));
            ImVec2 pos1 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(n.x, 1.f - n.y));
            window->DrawList->AddRectFilled(pos0, pos1, ImColor(base_color.x, base_color.y, base_color.z, base_color.w * prev_o));
        }
    }
}

IMGUI_API void PlotVariance(const float* avg, const float* var, int count, float var_scl, ImU32 line_color, ImU32 fill_color) {
    if (!ps.id) return;
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems) return;
    if (count < 2) return;

    ImVec2 prev_c = ImVec2(0, avg[0]);
    float pw = var[0] * var_scl / (ps.coord_view.Max.y - ps.coord_view.Min.y);
    ImVec2 p = compute_frame_coord(ps.coord_view, prev_c);
    ImVec2 prev_screen_pos_btm = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(p.x, ImClamp(1.f - p.y + pw, 0.f, 1.f)));
    ImVec2 prev_screen_pos_top = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(p.x, ImClamp(1.f - p.y - pw, 0.f, 1.f)));
    for (int i = 1; i < count; i++) {
        ImVec2 next_c = ImVec2((float)i, avg[i]);

        if (prev_c.x < ps.coord_view.Min.x && next_c.x < ps.coord_view.Min.x) continue;
        if (prev_c.x > ps.coord_view.Max.x && next_c.x > ps.coord_view.Max.x) break;

        ImVec2 n = compute_frame_coord(ps.coord_view, next_c);
        float nw = var[i] * var_scl / (ps.coord_view.Max.y - ps.coord_view.Min.y);

        ImVec2 next_screen_pos_btm = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(n.x, ImClamp(1.f - n.y + nw, 0.f, 1.f)));
        ImVec2 next_screen_pos_top = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(n.x, ImClamp(1.f - n.y - nw, 0.f, 1.f)));

        if ((fill_color & IM_COL32_A_MASK) > 0) {
            int flags = GetCurrentWindow()->DrawList->Flags;
            GetCurrentWindow()->DrawList->Flags &= ~ImDrawListFlags_AntiAliasedFill;
            GetCurrentWindow()->DrawList->AddQuadFilled(prev_screen_pos_top, prev_screen_pos_btm, next_screen_pos_btm, next_screen_pos_top,
                                                        fill_color);
            GetCurrentWindow()->DrawList->Flags = flags;
        }
        if ((line_color & IM_COL32_A_MASK) > 0) {
            GetCurrentWindow()->DrawList->AddLine(prev_screen_pos_btm, next_screen_pos_btm, line_color);
            GetCurrentWindow()->DrawList->AddLine(prev_screen_pos_top, next_screen_pos_top, line_color);
        }

        prev_c = next_c;
        prev_screen_pos_btm = {(next_screen_pos_btm.x), next_screen_pos_btm.y};
        prev_screen_pos_top = {(next_screen_pos_top.x), next_screen_pos_top.y};
    }
}

IMGUI_API void PlotValues(const char* line_label, const float* values, int count, ImU32 line_color) {
    (void)line_label;
    if (!ps.id) return;
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems) return;
    if (count < 2) return;

    ImVec2 prev_c = ImVec2(0, values[0]);
    for (int i = 1; i < count; i++) {
        ImVec2 next_c = ImVec2((float)i, values[i]);

        if (prev_c.x < ps.coord_view.Min.x && next_c.x < ps.coord_view.Min.x) continue;
        if (prev_c.y < ps.coord_view.Min.y && next_c.y < ps.coord_view.Min.y) continue;
        if (prev_c.x > ps.coord_view.Max.x && next_c.x > ps.coord_view.Max.x) break;
        if (prev_c.y > ps.coord_view.Max.y && next_c.y > ps.coord_view.Max.y) continue;

        ImVec2 p = compute_frame_coord(ps.coord_view, prev_c);
        ImVec2 n = compute_frame_coord(ps.coord_view, next_c);

        ImVec2 pos0 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(p.x, 1.f - p.y));
        ImVec2 pos1 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(n.x, 1.f - n.y));
        window->DrawList->AddLine(pos0, pos1, line_color);

        prev_c = next_c;
    }
}

IMGUI_API void EndPlot() {
    IM_ASSERT(ps.id != 0);
    ps.id = 0;
}

IMGUI_API void DrawHistogram(ImVec2 frame_min, ImVec2 frame_max, const float* values, int count, float max_val, ImU32 color) {
    if (max_val == 0.f) {
        for (int i = 0; i < count; i++) {
            max_val = ImMax(max_val, values[i]);
        }
        if (max_val == 0.f) return;
    }

    PushClipRect(frame_min, frame_max, true);
    const float x_scl = 1.f / (float)(count - 1);
    const float y_scl = 1.f / max_val;
    for (int i = 1; i < count; i++) {
        ImVec2 pos0 = ImLerp(frame_min, frame_max, ImVec2((i - 1) * x_scl, 1.f));
        ImVec2 pos1 = ImLerp(frame_min, frame_max, ImVec2(i * x_scl, 1.f - values[i - 1] * y_scl));

        if (pos1.x >= pos0.x + 2.0f) pos1.x -= 1.0f;
        GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, color);
    }
    PopClipRect();
}

IMGUI_API void DrawFilledLine(ImVec2 frame_min, ImVec2 frame_max, const float* values, int count, float max_val, ImU32 line_color, ImU32 fill_color) {
    if (max_val == 0.f) {
        for (int i = 0; i < count; i++) {
            max_val = ImMax(max_val, values[i]);
        }
        if (max_val == 0.f) return;
    }

    PushClipRect(frame_min, frame_max, true);
    const float x_scl = 1.f / (float)(count - 1);
    const float y_scl = 1.f / max_val;
    for (int i = 1; i < count; i++) {
        ImVec2 pos0 = ImLerp(frame_min, frame_max, ImVec2((i - 1) * x_scl, 1.f - values[i - 1] * y_scl));
        ImVec2 pos1 = ImLerp(frame_min, frame_max, ImVec2(i * x_scl, 1.f - values[i] * y_scl));
        pos0 = ImVec2(roundf(pos0.x), roundf(pos0.y));
        pos1 = ImVec2(roundf(pos1.x), roundf(pos1.y));

        GetCurrentWindow()->DrawList->AddLine(pos0, pos1, line_color);
        const ImVec2 points[4]{{pos0.x, frame_max.y}, pos0, pos1, {pos1.x, frame_max.y}};
        GetCurrentWindow()->DrawList->AddConvexPolyFilled(points, 4, fill_color);
    }
    PopClipRect();
}

IMGUI_API bool PlotHistogram(const char* label, ImVec2 frame_size, const float* values, int count, bool periodic, ImVec2 value_range,
                             ImVec2* selection_range) {
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems) return false;
    const ImGuiID id = window->GetID(label);
    ImGuiContext& ctx = *GImGui;
    const ImGuiStyle& style = ctx.Style;

    if (frame_size.x == 0.0f) frame_size.x = CalcItemWidth();

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, 0)) return false;

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    static float* drag_target = nullptr;
    bool modifying = false;

    if (GetActiveID() == id) {
        if (!ctx.IO.MouseDown[0]) {
            ClearActiveID();
            drag_target = nullptr;
        }
    }

    if (IsItemHovered() && ctx.IO.MouseClicked[0]) {
        SetActiveID(id, window);
        FocusWindow(window);
    }

    if (GetActiveID() == id) {
        if (selection_range && ctx.IO.MouseDown[0] && ctx.IO.KeyCtrl && ctx.IO.MouseDelta.x != 0) {
            modifying = true;
            float t = ImClamp((ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.f, 1.f);
            static int periodic_counter = 0;
            if (ctx.IO.MouseClicked[0]) {
                selection_range->x = ImLerp(value_range.x, value_range.y, t);
                selection_range->y = selection_range->x;

                drag_target = &selection_range->y;
                periodic_counter = 0;
            }
            if (periodic) {
                if (ctx.IO.MousePos.x < inner_bb.Min.x) {
                    ctx.IO.WantSetMousePos = true;
                    ctx.IO.MousePos.x = inner_bb.Max.x - (inner_bb.Min.x - ctx.IO.MousePos.x);
                    periodic_counter--;
                }
                if (ctx.IO.MousePos.x > inner_bb.Max.x) {
                    ctx.IO.WantSetMousePos = true;
                    ctx.IO.MousePos.x = inner_bb.Min.x + (ctx.IO.MousePos.x - inner_bb.Max.x);
                    periodic_counter++;
                }
            }
            if (drag_target) {
                *drag_target = ImLerp(value_range.x, value_range.y, t);
                if (drag_target == &selection_range->y) {
                    if (selection_range->y < selection_range->x) {
                        if (periodic_counter == 0) {
                            ImSwap(selection_range->x, selection_range->y);
                            drag_target = &selection_range->x;
                        }
                    }
                } else if (drag_target == &selection_range->x) {
                    if (selection_range->x > selection_range->y) {
                        if (periodic_counter == 0) {
                            ImSwap(selection_range->x, selection_range->y);
                            drag_target = &selection_range->y;
                        }
                    }
                }
            }
        }
    }

    float max_val = 0.f;
    for (int i = 0; i < count; i++) {
        max_val = ImMax(max_val, values[i]);
    }

    if (max_val == 0.f) return false;

    for (int i = 1; i < count; i++) {
        ImVec2 tp0 = ImVec2((i - 1) / (float)(count - 1), 1.f);
        ImVec2 tp1 = ImVec2(i / (float)(count - 1), 1.f - values[i - 1] / max_val);

        ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, tp0);
        ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, tp1);

        if (pos1.x >= pos0.x + 2.0f) pos1.x -= 1.0f;
        window->DrawList->AddRectFilled(pos0, pos1, GetColorU32(ImGuiCol_PlotHistogram));
    }

    if (selection_range) {
        // Draw selection range
        const float t0 = (selection_range->x - value_range.x) / (value_range.y - value_range.x);
        const float t1 = (selection_range->y - value_range.x) / (value_range.y - value_range.x);

        if (t0 != t1) {
            ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
            ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));

            if (t1 < t0) {
                window->DrawList->AddRectFilled(pos0, inner_bb.Max, SELECTION_RANGE_COLOR);
                window->DrawList->AddRectFilled(inner_bb.Min, pos1, SELECTION_RANGE_COLOR);
            } else {
                window->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }
        }
    }

    if (IsItemHovered()) {
        window->DrawList->AddLine(ImVec2(ctx.IO.MousePos.x, inner_bb.Min.y), ImVec2(ctx.IO.MousePos.x, inner_bb.Max.y), 0xffffffff);
        float t = (ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
        float val = values[ImClamp((int)(t * (count - 1)), 0, count - 1)];
        BeginTooltip();
        Text("%.3f: %.3f\n", ImLerp(value_range.x, value_range.y, t), val);
        EndTooltip();
    }

    RenderTextClipped(ImVec2(frame_bb.Min.x + style.FramePadding.x, frame_bb.Min.y + style.FramePadding.y), frame_bb.Max, label, NULL, NULL,
                      ImVec2(0.0f, 0.0f));

    return modifying;
}

IMGUI_API bool PlotPeriodic(const char* label, float outer_radius, float inner_radius_ratio, const float* values, int count, ImVec2 value_range,
                            ImU32 line_color) {
    ImGuiWindow* window = GetCurrentWindow();
    const ImGuiID id = window->GetID(label);
    ps.id = id;

    if (window->SkipItems) return false;

    ImGuiContext& ctx = *GImGui;
    const ImGuiStyle& style = ctx.Style;

    if (outer_radius == 0.0f) outer_radius = CalcItemWidth();
    ImVec2 frame_size(2.f * outer_radius, 2.f * outer_radius);

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, 0)) return false;

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    if (GetActiveID() == id) {
        if (!IsItemHovered()) {
            ClearActiveID();
        }
    }

    if (IsItemHovered() && ctx.IO.MouseClicked[0]) {
        SetActiveID(id, window);
        FocusWindow(window);
    }

    float max_val = 0.f;
    for (int i = 0; i < count; i++) {
        max_val = ImMax(max_val, values[i]);
    }

    ImVec2 center = ImLerp(inner_bb.Min, inner_bb.Max, 0.5f);
    float inner_radius = outer_radius * inner_radius_ratio;
    window->DrawList->AddCircleFilled(center, outer_radius, GetColorU32(ImGuiCol_ChildBg), 32);
    window->DrawList->AddCircleFilled(center, inner_radius, GetColorU32(ImGuiCol_WindowBg), 32);
    // window->DrawList->AddCircle(center, inner_radius, 0xffffffff, 64);
    window->DrawList->AddCircle(center, outer_radius, 0xffffffff, 64);

    if (count == 0) return true;
    const ImVec2 zero_line(center.x, center.y - ImLerp(inner_radius, outer_radius, values[0] / max_val));
    const ImVec2 zero_inner(center.x, center.y - inner_radius);
    ImVec2 prev_line = zero_line;
    ImVec2 next_line;
    ImVec2 prev_inner = zero_inner;
    ImVec2 next_inner;
    const float PI_HALF = IM_PI * 0.5f;
    const float TWO_PI = IM_PI * 2.f;
    for (int i = 1; i < count; i++) {
        float angle = PI_HALF - TWO_PI * (i / (float)count);
        float radius = ImLerp(inner_radius, outer_radius, values[i] / max_val);
        next_line = ImVec2(center.x + cosf(angle) * radius, center.y - sinf(angle) * radius);
        next_inner = ImVec2(center.x + cosf(angle) * inner_radius, center.y - sinf(angle) * inner_radius);
        ImVec2 points[4]{prev_inner, prev_line, next_line, next_inner};
        window->DrawList->AddConvexPolyFilled(points, 4, 0x44aabb33);
        window->DrawList->AddLine(prev_line, next_line, line_color);
        prev_line = next_line;
        prev_inner = next_inner;
    }
    ImVec2 points[4]{prev_inner, prev_line, zero_line, zero_inner};
    window->DrawList->AddConvexPolyFilled(points, 4, 0x44aabb33);
    window->DrawList->AddLine(prev_line, zero_line, line_color);

    if (IsItemHovered()) {
        const ImVec2 delta(ctx.IO.MousePos.x - center.x, ctx.IO.MousePos.y - center.y);
        const float len_sqr = ImLengthSqr(delta);
        if (inner_radius * inner_radius < len_sqr && len_sqr < outer_radius * outer_radius) {
            float angle = atan2f(delta.y, delta.x);
            // @TODO: FIX THIS!
            float a2 = (angle + PI_HALF);
            if (a2 < 0) a2 += TWO_PI;
            if (a2 > TWO_PI) a2 -= TWO_PI;
            float t = a2 / TWO_PI;
            float val = values[ImClamp((int)(t * count), 0, count - 1)] / max_val;

            ImVec2 inner(center.x + inner_radius * cosf(angle), center.y + inner_radius * sinf(angle));
            ImVec2 outer(center.x + outer_radius * cosf(angle), center.y + outer_radius * sinf(angle));
            window->DrawList->AddLine(inner, outer, 0xffffffff);
            BeginTooltip();
            Text("%.3f: %.3f\n", ImLerp(value_range.x, value_range.y, t), val);
            EndTooltip();
        }
    }

    return true;
}

}  // namespace ImGui
