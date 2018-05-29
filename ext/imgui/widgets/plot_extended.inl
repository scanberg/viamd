namespace ImGui {

int PlotExtended(ImGuiPlotType plot_type, const char* label, const float* data, int count, int offset, int* selected_idx, const char* tooltip,
                 std::function<std::string(int, float)> caption_func, float scale_min, float scale_max, ImVec2 graph_size, int res,
                 const int* highlighted_indices, int highlight_count) {
    ImGuiWindow* window = GetCurrentWindow();

    int value = -1;
    if (window->SkipItems) return value;

    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;

    PushItemWidth(-1);

    // const ImVec2 label_size = CalcTextSize(label, NULL, true);
    if (graph_size.x == 0.0f) graph_size.x = CalcItemWidth();
    if (graph_size.y == 0.0f) graph_size.y = (style.FramePadding.y * 2);

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(graph_size.x, graph_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, NULL)) return value;

    if (res < 1) return value;

    // Determine scale from values if not specified
    if (scale_min == FLT_MAX || scale_max == FLT_MAX) {
        float v_min = FLT_MAX;
        float v_max = -FLT_MAX;
        for (int i = 0; i < count; i++) {
            const float v = data[i];
            v_min = ImMin(v_min, v);
            v_max = ImMax(v_max, v);
        }
        if (scale_min == FLT_MAX) scale_min = v_min;
        if (scale_max == FLT_MAX) scale_max = v_max;
    }

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    if (count > 0 && graph_size.x > 0 && graph_size.y > 0) {
        auto get_data = [data, offset](int idx) -> float { return data[idx + offset]; };

        int res_w = ImMin((int)graph_size.x / res, count) + ((plot_type == ImGuiPlotType_Lines) ? -1 : 0);
        int item_count = count + ((plot_type == ImGuiPlotType_Lines) ? -1 : 0);

        const ImU32 col_base = GetColorU32((plot_type == ImGuiPlotType_Lines) ? ImGuiCol_PlotLines : ImGuiCol_PlotHistogram);
        const ImU32 col_hovered = GetColorU32((plot_type == ImGuiPlotType_Lines) ? ImGuiCol_PlotLinesHovered : ImGuiCol_PlotHistogramHovered);

        const ImGuiID id = window->GetID(label);

        // Tooltip on hover
        int v_hovered = -1;
        if (IsItemHovered()) {
            const float t = ImClamp((g.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.0f, 1.0f);
            const int v_idx = (int)(t * (count - 1) + 0.5f);
            IM_ASSERT(v_idx >= 0 && v_idx < count);

            if (tooltip) {
                std::string formatted = tooltip;
                auto i_pos = formatted.find_first_of("[idx]");
                if (i_pos != std::string::npos) {
                    formatted = formatted.substr(0, i_pos) + std::to_string(v_idx) + formatted.substr(i_pos + 5);
                }

                auto v_pos = formatted.find_first_of("[val]");
                if (v_pos != std::string::npos) {
                    formatted = formatted.substr(0, v_pos) + std::to_string(get_data(v_idx)) + formatted.substr(v_pos + 5);
                }

                SetTooltip(formatted.c_str());
                // float val = get_data(v_idx);
            }

            // Draw vertical-line
            const ImVec2 btm = ImVec2(roundf(g.IO.MousePos.x + 0.5f), inner_bb.Min.y);
            const ImVec2 top = ImVec2(roundf(g.IO.MousePos.x + 0.5f), inner_bb.Max.y);
            window->DrawList->AddLine(btm, top, col_base);

            v_hovered = v_idx;
        }

        // Draw zero axis
        if (plot_type == ImGuiPlotType_Lines && scale_min < 0.f && scale_max > 0.f) {
            auto t = 1.0f - ImSaturate((0.f - scale_min) / (scale_max - scale_min));
            auto y_pos = roundf(ImLerp(inner_bb.Min.y, inner_bb.Max.y, t) + 0.5f);
            auto color = IM_COL32(255, 255, 255, 140);
            window->DrawList->AddLine(ImVec2(inner_bb.Min.x, y_pos), ImVec2(inner_bb.Max.x, y_pos), color);
        }

        const float t_step = 1.0f / (float)res_w;
        const int samples_per_step = count / res_w;

        // Draw highlighted indices
        if (highlight_count > 0) {
            const ImU32 highlight_col = IM_COL32(255, 255, 255, 30);
            float t0 = 0.f;
            int v0_idx = offset;
            int current_idx = 0;
            for (int n = 0; n < res_w; n++) {
                const float t1 = t0 + t_step;
                const int v1_idx = ImClamp((int)(t1 * item_count + 0.5f), 0, count - 1) + offset;

                bool color_background = false;
                while (current_idx < highlight_count && v0_idx <= highlighted_indices[current_idx] && highlighted_indices[current_idx] < v1_idx) {
                    current_idx++;
                    color_background = true;
                }

                if (color_background) {
                    ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                    ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                    window->DrawList->AddRectFilled(pos0, pos1, highlight_col);
                }

                t0 = t1;
                v0_idx = v1_idx;
            }
        }

        auto avg_samples = [get_data](const int beg, const int end) -> float {
            float val = 0.f;
            for (int i = beg; i < end; i++) {
                val += get_data(i);
            }
            return val /= (float)(end - beg);
        };

        int v0_idx = 0;
        float v0 = avg_samples(v0_idx, (int)(t_step * item_count + 0.5f));
        float t0 = 0.0f;
        ImVec2 tp0 = ImVec2(t0,
                            1.0f - ImSaturate((v0 - scale_min) / (scale_max - scale_min)));  // Point in the normalized space of our target rectangle

        for (int n = 0; n < res_w; n++) {
            const float t1 = t0 + t_step;
            const int v1_idx = ImClamp((int)(t1 * item_count + 0.5f), 0, count - 1);
            float v1 = avg_samples(v0_idx, v1_idx);
            const ImVec2 tp1 = ImVec2(t1, 1.0f - ImSaturate((v1 - scale_min) / (scale_max - scale_min)));

            // NB: Draw calls are merged together by the DrawList system. Still, we should
            // render our batch are lower level to save a bit of CPU.
            ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, tp0);
            ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, (plot_type == ImGuiPlotType_Lines) ? tp1 : ImVec2(tp1.x, 1.0f));

            auto col = v_hovered == v1_idx ? col_hovered : col_base;

            if (plot_type == ImGuiPlotType_Lines) {
                window->DrawList->AddLine(pos0, pos1, col);
            } else if (plot_type == ImGuiPlotType_Histogram) {
                if (pos1.x >= pos0.x + 2.0f) pos1.x -= 1.0f;
                window->DrawList->AddRectFilled(pos0, pos1, col);
            }

            t0 = t1;
            tp0 = tp1;
            v0_idx = v1_idx;
        }

        if (IsItemHovered() && g.IO.MouseClicked[0]) {
            SetActiveID(id, window);
            FocusWindow(window);
        }

        if (g.ActiveId == id) {
            if (g.IO.MouseDown[0]) {
                const float t = ImClamp((g.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.0f, 1.0f);
                const int v_idx = (int)(t * (count - 1) + 0.5f);
                IM_ASSERT(v_idx >= 0 && v_idx < count);
                value = offset + v_idx;
            } else {
                ClearActiveID();
            }
        }

        if (selected_idx) {
            // Draw line if in visible range
            if (offset <= *selected_idx && *selected_idx < offset + count) {
                const float t = ImClamp((*selected_idx - offset) / (float)(count - 1), 0.0f, 1.0f);
                // cast to int for sharp edges due to dpi-scaling
                const float x = roundf(ImLerp(inner_bb.Min.x, inner_bb.Max.x, t) + 0.5f);
                ImVec2 btm = ImVec2(x, inner_bb.Min.y);
                ImVec2 top = ImVec2(x, inner_bb.Max.y);
                window->DrawList->AddLine(btm, top, col_hovered);
            }
        }
    }

    // Text overlay

    // Caption
    if (caption_func && *selected_idx > -1) {
        std::string caption = caption_func(*selected_idx, data[*selected_idx]);
        RenderTextClipped(ImVec2(frame_bb.Min.x, frame_bb.Min.y + style.FramePadding.y), frame_bb.Max, caption.c_str(), NULL, NULL,
                          ImVec2(0.5f, 0.0f));
    }

    // Label
    RenderTextClipped(ImVec2(frame_bb.Min.x + style.FramePadding.x, frame_bb.Min.y + style.FramePadding.y), frame_bb.Max, label, NULL, NULL,
                      ImVec2(0.0f, 0.0f));

    return value;
}

int PlotLinesExtended(const char* label, const float* values, int count, int offset, int* selected_idx, const char* tooltip,
                      std::function<std::string(int, float)> caption_func, float scale_min, float scale_max, ImVec2 graph_size, int res,
                      const int* highlighted_indices, int highlight_count) {
    return PlotExtended(ImGuiPlotType_Lines, label, values, count, offset, selected_idx, tooltip, caption_func, scale_min, scale_max, graph_size, res,
                        highlighted_indices, highlight_count);
}

int PlotHistogramExtended(const char* label, const float* values, int count, int offset, int* selected_idx, const char* tooltip,
                          std::function<std::string(int, float)> caption_func, float scale_min, float scale_max, ImVec2 graph_size, int res,
                          const int* highlighted_indices, int highlight_count) {
    return PlotExtended(ImGuiPlotType_Histogram, label, values, count, offset, selected_idx, tooltip, caption_func, scale_min, scale_max, graph_size,
                        res, highlighted_indices, highlight_count);
}

PlotFrame BeginPlotFrame(const char* label, ImVec2 frame_size, int offset, int count, float scale_min, float scale_max,
                         std::function<std::string(int)> x_label_func, const int* highlight_indices, int highlight_count) {
    ImGuiWindow* window = GetCurrentWindow();

    if (window->SkipItems) return {};

    ImGuiContext& ctx = *GImGui;
    const ImGuiStyle& style = ctx.Style;

    PushItemWidth(-1);

    if (frame_size.x == 0.0f) frame_size.x = CalcItemWidth();
    if (frame_size.y == 0.0f) frame_size.y = (style.FramePadding.y * 2);

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, NULL)) return {};

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    const ImGuiID id = window->GetID(label);

    // Draw zero axis
    if (scale_min < 0.f && scale_max > 0.f) {
        float t = 1.0f - ImSaturate((0.f - scale_min) / (scale_max - scale_min));
        float y_pos = roundf(ImLerp(inner_bb.Min.y, inner_bb.Max.y, t) + 0.5f);
        ImU32 color = IM_COL32(255, 255, 255, 140);
        window->DrawList->AddLine(ImVec2(inner_bb.Min.x, y_pos), ImVec2(inner_bb.Max.x, y_pos), color);
    }

    // Draw highlighted indices
    if (highlight_count > 0) {
        const int res_w = ImMin((int)(inner_bb.Max.x - inner_bb.Min.x), count);
        const float t_step = 1.0f / (float)res_w;
        const ImU32 highlight_col = IM_COL32(255, 255, 255, 30);

        // float t0 = 0.f;
        // int l0_idx = 0;
        // int g0_idx = offset;
        // int v0_idx = offset;
        int current_idx = 0;

        for (int n = 0; n < res_w; n++) {
            // const float t1 = ImClamp(t0 + t_step, 0.f, 1.f);
            const int l0_idx = n * count / res_w;
            const int g0_idx = l0_idx + offset;

            const int l1_idx = (n + 1) * count / res_w;
            const int g1_idx = l1_idx + offset;

            // const int v1_idx = (n + 1) * count / res_w + offset; //ImClamp((int)(t1 * count), 0, count - 1) + offset;

            bool color_background = false;
            while (current_idx < highlight_count && g0_idx <= highlight_indices[current_idx] && highlight_indices[current_idx] < g1_idx) {
                current_idx++;
                color_background = true;
            }

            if (color_background) {
                const float t0 = l0_idx / (float)count;
                const float t1 = l1_idx / (float)count;
                ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                window->DrawList->AddRectFilled(pos0, pos1, highlight_col);
            }

            // t0 = t1;
            // v0_idx = v1_idx;
        }
    }

    std::string tooltip;

    int hovered_index = -1;
    if (IsItemHovered()) {
        const float t = ImSaturate((ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x));
        hovered_index = ImClamp((int)(t * count), 0, count - 1) + offset;

        if (x_label_func) {
            tooltip += x_label_func(hovered_index) + '\n';
        } else {
            tooltip += std::to_string(hovered_index) + ":\n";
        }
    }

    return {id, label, hovered_index, count, offset, scale_min, scale_max, frame_bb.Min, frame_bb.Max, inner_bb.Min, inner_bb.Max, tooltip};
}

int EndPlotFrame(const PlotFrame& frame, int selected_idx) {
    ImGuiContext& ctx = *GImGui;
    ImGuiWindow* window = GetCurrentWindow();
    bool manhattan_mode = false;

    const ImU32 col_base = GetColorU32(ImGuiCol_PlotLines);

    if (frame.hovered_index > -1) {
        // Draw vertical-line
        const ImVec2 btm = ImVec2(roundf(ctx.IO.MousePos.x + 0.5f), frame.inner_bb_min.y);
        const ImVec2 top = ImVec2(roundf(ctx.IO.MousePos.x + 0.5f), frame.inner_bb_max.y);
        window->DrawList->AddLine(btm, top, col_base);

        SetTooltip("%s", frame.tooltip.c_str());
    }

    if (frame.hovered_index > -1 && ctx.IO.MouseClicked[0]) {
        SetActiveID(frame.id, window);
        FocusWindow(window);
    }

    int value = -1;
    if (ctx.ActiveId == frame.id) {
        if (ctx.IO.MouseDown[0]) {
            const float t = ImClamp((ctx.IO.MousePos.x - frame.inner_bb_min.x) / (frame.inner_bb_max.x - frame.inner_bb_min.x), 0.0f, 1.0f);
            const int v_idx = (int)(t * frame.count);
            value = ImClamp(v_idx, 0, frame.count - 1) + frame.offset;
        } else {
            ClearActiveID();
        }
    }

    if (selected_idx > -1) {
        float t_offset = 0.f;
        if (manhattan_mode) {
            t_offset = 0.5f;
        }
        // Draw line if in visible range
        if (frame.offset <= selected_idx && selected_idx < frame.offset + frame.count) {
            const float t = ImClamp((selected_idx - frame.offset + t_offset) / (float)(frame.count), 0.0f, 1.0f);
            // cast to int for sharp edges due to dpi-scaling
            const float x = roundf(ImLerp(frame.inner_bb_min.x, frame.inner_bb_max.x, t) + 0.5f);
            ImVec2 btm = ImVec2(x, frame.inner_bb_min.y);
            ImVec2 top = ImVec2(x, frame.inner_bb_max.y);
            window->DrawList->AddLine(btm, top, GetColorU32(ImGuiCol_PlotLinesHovered));
        }
    }

    const ImGuiStyle& style = ctx.Style;

    // Caption
    if (frame.label) {
        RenderTextClipped(ImVec2(frame.frame_bb_min.x, frame.frame_bb_min.y + style.FramePadding.y), frame.frame_bb_max, frame.label, NULL, NULL,
                          ImVec2(0.0f, 0.0f));
    }

    return value;
}

void PlotFrameLine(PlotFrame& frame, const char* label, const float* values, FrameLineStyle style, int selected_idx) {

    if (values) {
        ImGuiWindow* window = GetCurrentWindow();

        float w = frame.inner_bb_max.x - frame.inner_bb_min.x;
        int res_w = ImMin((int)w, frame.count);
        // bool manhattan_mode = style.line_mode == LineMode_Histogram;
        bool fill = style.fill_mode > FillMode_None;

        if (style.fill_mode == FillMode_Automatic) {
            ImVec4 tmp = ColorConvertU32ToFloat4(style.line_color);
            tmp.w = 0.25f;
            style.custom_fill_color = ColorConvertFloat4ToU32(tmp);
        }

        const float t_step = 1.0f / (float)res_w;

        auto avg_samples = [values](const int beg, const int end) -> float {
            int ext = end - beg;
            if (ext == 0) return values[beg];

            float val = 0.f;
            for (int i = beg; i < end; i++) {
                val += values[i];
            }
            return val /= (float)ext;
        };

        int v0_idx = frame.offset;
        float v0 = avg_samples(v0_idx, v0_idx + (int)(t_step * frame.count));
        float t0 = 0.0f;
        ImVec2 tp0 = ImVec2(t0,
                            1.0f - ImSaturate((v0 - frame.scale_min) /
                                              (frame.scale_max - frame.scale_min)));  // Point in the normalized space of our target rectangle

        // Brighten the color in HSV and store result in hovered_color
        ImVec4 rgba = ImGui::ColorConvertU32ToFloat4(style.line_color);
        ImGui::ColorConvertRGBtoHSV(rgba.x, rgba.y, rgba.z, rgba.x, rgba.y, rgba.z);
        rgba.z = ImMin(rgba.z * 1.3f, 1.f);
        ImGui::ColorConvertHSVtoRGB(rgba.x, rgba.y, rgba.z, rgba.x, rgba.y, rgba.z);
        ImU32 hovered_color = ImGui::ColorConvertFloat4ToU32(rgba);

        float zero =
            ImLerp(frame.inner_bb_min.y, frame.inner_bb_max.y, 1.f - ImSaturate((0 - frame.scale_min) / (frame.scale_max - frame.scale_min)));

        for (int n = 0; n < res_w; n++) {
            const float t1 = t0 + t_step;
            const int v1_idx = ImClamp((int)(t1 * frame.count), 0, frame.count) + frame.offset;
            float v1 = avg_samples(v0_idx, v1_idx);
            const ImVec2 tp1 = ImVec2(t1, 1.0f - ImSaturate((v1 - frame.scale_min) / (frame.scale_max - frame.scale_min)));

            // NB: Draw calls are merged together by the DrawList system. Still, we should
            // render our batch are lower level to save a bit of CPU.
            ImVec2 pos0 = ImLerp(frame.inner_bb_min, frame.inner_bb_max, tp0);
            ImVec2 pos2 = ImLerp(frame.inner_bb_min, frame.inner_bb_max, tp1);

            ImU32 c = style.line_color;
            if (frame.hovered_index > -1 && frame.hovered_index == v0_idx) c = hovered_color;

            if (style.line_mode == LineMode_Histogram) {
                ImVec2 pos1(pos0.x, pos2.y);
                window->DrawList->AddLine(pos0, pos1, c);
                window->DrawList->AddLine(pos1, pos2, c);

                if (fill) {
                    ImVec2 bl(pos0.x, zero);
                    ImVec2 tl = pos1;
                    ImVec2 tr = pos2;
                    ImVec2 br(pos2.x, zero);
                    window->DrawList->AddQuadFilled(bl, tl, tr, br, style.custom_fill_color);
                }
            } else if (style.line_mode == LineMode_Line) {
                window->DrawList->AddLine(pos0, pos2, c);

                if (fill) {
                    ImVec2 points[4]{{pos0.x, zero}, pos0, pos2, {pos2.x, zero}};

                    window->DrawList->AddConvexPolyFilled(points, 4, style.custom_fill_color);
                }
            }

            t0 = t1;
            tp0 = tp1;
            v0_idx = v1_idx;
        }

        if (frame.hovered_index > -1) {
            // const char* c = reinterpret_cast<char*>(&style.line_color);
            // std::string color_modifier{ '\033', c[0], c[1], c[2], c[3] };
            // frame.tooltip += color_modifier + label + ": " + std::to_string(values[frame.hovered_index]) + "\n";
            frame.tooltip += std::string(label) + ": " + std::to_string(values[frame.hovered_index]) + "\n";
        }

        if (selected_idx > -1) {
            // Draw line if in visible range
            if (frame.offset <= selected_idx && selected_idx < frame.offset + frame.count) {
                const float t = ImClamp((selected_idx - frame.offset) / (float)(frame.count), 0.0f, 1.0f);
                // cast to int for sharp edges due to dpi-scaling
                const float x = roundf(ImLerp(frame.inner_bb_min.x, frame.inner_bb_max.x, t) + 0.5f);
                ImVec2 btm = ImVec2(x, frame.inner_bb_min.y);
                ImVec2 top = ImVec2(x, frame.inner_bb_max.y);
                window->DrawList->AddLine(btm, top, style.line_color);
            }
        }
    }
}

struct PlotState {
    ImGuiID id = 0;
    ImRect inner_bb;
    ImRect coord_view;
    ImVec2* selection_range;
    float selection_start;
    bool is_selecting = false;
};

static PlotState ps;

const ImU32 HOVER_LINE_COLOR = 0xaaffffff;
const ImU32 CURRENT_LINE_COLOR = 0xaa33ffff;
const ImU32 AXIS_COLOR = 0xaaffffff;
const ImU32 SELECTION_RANGE_COLOR = 0x33bbbbbb;

IMGUI_API bool BeginPlot(const char* label, ImVec2 frame_size, ImVec2 x_range, ImVec2 y_range, float* x_val, ImVec2* selection_range,
                         LinePlotFlags flags) {
    ImGuiWindow* window = GetCurrentWindow();
    const ImGuiID id = window->GetID(label);
    ps.id = id;

    if (window->SkipItems) return false;

    ImGuiContext& ctx = *GImGui;
    const ImGuiStyle& style = ctx.Style;

    if (frame_size.x == 0.0f) frame_size.x = CalcItemWidth();
    if (frame_size.y == 0.0f) frame_size.y = (style.FramePadding.y * 2);

    const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
    const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
    const ImRect total_bb(frame_bb.Min, frame_bb.Max);
    ItemSize(total_bb, style.FramePadding.y);
    if (!ItemAdd(total_bb, NULL)) return false;

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

    if (IsItemHovered() && ctx.IO.MouseClicked[0]) {
        SetActiveID(id, window);
        FocusWindow(window);
    }

    if (selection_range && IsItemHovered() && ctx.IO.MouseClicked[1] && ctx.IO.KeyCtrl) {
        *selection_range = x_range;
    }

    bool interacting_x_val = false;

    if (ctx.ActiveId == id) {
        if (selection_range && ctx.IO.MouseClicked[0] && ctx.IO.KeyCtrl) {
            float t = (ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
            ps.selection_start = ImLerp(x_range.x, x_range.y, t);
            selection_range->x = ps.selection_start;
            selection_range->y = ps.selection_start;
            ps.is_selecting = true;
        } else if (ps.is_selecting) {
            float t = (ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
            float v = ImLerp(x_range.x, x_range.y, t);
            if (v < selection_range->x) {
                selection_range->x = v;
            } else if (v > selection_range->x && v < selection_range->y) {
                if (ps.selection_start < v) {
                    selection_range->y = v;
                } else {
                    selection_range->x = v;
                }
            } else if (v > selection_range->y) {
                selection_range->y = v;
            }
        } else if (ctx.IO.MouseDown[0]) {
            if (x_val) {
                float t = ImClamp((ctx.IO.MousePos.x - ps.inner_bb.Min.x) / (ps.inner_bb.Max.x - ps.inner_bb.Min.x), 0.f, 1.f);
                *x_val = ImLerp(x_range.x, x_range.y, t);
                interacting_x_val = true;
            }
        }

        if (!ctx.IO.MouseDown[0] && !IsItemHovered()) {
            ClearActiveID();
            ps.is_selecting = false;
        }
    }

    if (flags & LinePlotFlags_AxisY) {
        if (x_range.x < 0 && 0 < x_range.y) {
            const float t = (0 - x_range.x) / (x_range.y - x_range.x);
            const float x = ImLerp(inner_bb.Min.x, inner_bb.Max.x, t);
            window->DrawList->AddLine(ImVec2(x, inner_bb.Min.y), ImVec2(x, inner_bb.Max.y), AXIS_COLOR);
        }
    }

    if (flags & LinePlotFlags_ShowXVal && x_val) {
        if (x_range.x < *x_val && *x_val < x_range.y) {
            const float t = (*x_val - x_range.x) / (x_range.y - x_range.x);
            const float x = roundf(ImLerp(inner_bb.Min.x, inner_bb.Max.x, t));
            window->DrawList->AddLine(ImVec2(x, inner_bb.Min.y), ImVec2(x, inner_bb.Max.y), CURRENT_LINE_COLOR);
        }
    }

    if (flags & LinePlotFlags_AxisX) {
        if (y_range.x < 0 && 0 < y_range.y) {
            const float t = (0 - y_range.x) / (y_range.y - y_range.x);
            const float y = roundf(ImLerp(inner_bb.Min.y, inner_bb.Max.y, 1.f - t));
            window->DrawList->AddLine(ImVec2(roundf(inner_bb.Min.x), y), ImVec2(roundf(inner_bb.Max.x), y), AXIS_COLOR);
        }
    }

    if (IsItemHovered()) {
        // Draw vertical line to show current position
        ImVec2 pos0(roundf(ctx.IO.MousePos.x), inner_bb.Min.y);
        ImVec2 pos1(roundf(ctx.IO.MousePos.x), inner_bb.Max.y);
        window->DrawList->AddLine(pos0, pos1, HOVER_LINE_COLOR);
    }

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

    ps.inner_bb = inner_bb;
    ps.coord_view = ImRect(ImVec2(x_range.x, y_range.x), ImVec2(x_range.y, y_range.y));
    ps.selection_range = selection_range;

    RenderTextClipped(ImVec2(frame_bb.Min.x + style.FramePadding.x, frame_bb.Min.y + style.FramePadding.y), frame_bb.Max, label, NULL, NULL,
                      ImVec2(0.0f, 0.0f));

    return interacting_x_val;
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

            // if (prev_c.x < ps.coord_view.Min.x && next_c.x < ps.coord_view.Min.x) continue;
            // if (prev_c.y < ps.coord_view.Min.y && next_c.y < ps.coord_view.Min.y) continue;
            if (prev_c.x > ps.coord_view.Max.x && next_c.x > ps.coord_view.Max.x) break;
            // if (prev_c.y > ps.coord_view.Max.y && next_c.y > ps.coord_view.Max.y) continue;

            float px = ImClamp((prev_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
            float py = ImClamp((prev_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);
            float nx = ImClamp((next_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
            float ny = ImClamp((next_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);

            ImVec2 pos0 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(px, 1.f - py));
            ImVec2 pos1 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(nx, 1.f - ny));
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
    float prev_w = var[0] * var_scl;
    for (int i = 1; i < count; i++) {
        ImVec2 next_c = ImVec2((float)i, avg[i]);
        float next_w = var[i] * var_scl;

        if (prev_c.x < ps.coord_view.Min.x && next_c.x < ps.coord_view.Min.x) continue;
        if (prev_c.y < ps.coord_view.Min.y && next_c.y < ps.coord_view.Min.y) continue;
        if (prev_c.x > ps.coord_view.Max.x && next_c.x > ps.coord_view.Max.x) break;
        if (prev_c.y > ps.coord_view.Max.y && next_c.y > ps.coord_view.Max.y) continue;

        float px = ImClamp((prev_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
        float py = ImClamp((prev_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);
        float nx = ImClamp((next_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
        float ny = ImClamp((next_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);
        float pw = prev_w / (ps.coord_view.Max.y - ps.coord_view.Min.y);
        float nw = next_w / (ps.coord_view.Max.y - ps.coord_view.Min.y);

        ImVec2 pos0 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(px, 1.f - py + pw));
        ImVec2 pos1 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(px, 1.f - py - pw));
        ImVec2 pos2 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(nx, 1.f - ny + nw));
        ImVec2 pos3 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(nx, 1.f - ny - nw));

        pos0.x = roundf(pos0.x);
        pos1.x = roundf(pos1.x);
        pos2.x = roundf(pos2.x);
        pos3.x = roundf(pos3.x);

        if (pos0.x == pos2.x) continue;

        // float scl_prev_w = prev_w / (ps.coord_view.Max.y - ps.coord_view.Min.y);
        // float scl_next_w = next_w / (ps.coord_view.Max.y - ps.coord_view.Min.y);
        const ImVec2 pos[] = {pos1, pos0, pos2, pos3};

        // GetCurrentWindow()->DrawList->AddConvexPolyFilled(pos, 4, fill_color);
        GetCurrentWindow()->DrawList->AddQuadFilled(pos1, pos0, pos2, pos3, fill_color);
        // window->DrawList->AddLine(pos0, pos1, line_color);
        // window->DrawList->AddLine(pos2, pos3, line_color);

        prev_c = next_c;
        prev_w = next_w;
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

        float px = ImClamp((prev_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
        float py = ImClamp((prev_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);
        float nx = ImClamp((next_c.x - ps.coord_view.Min.x) / (ps.coord_view.Max.x - ps.coord_view.Min.x), 0.f, 1.f);
        float ny = ImClamp((next_c.y - ps.coord_view.Min.y) / (ps.coord_view.Max.y - ps.coord_view.Min.y), 0.f, 1.f);

        ImVec2 pos0 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(px, 1.f - py));
        ImVec2 pos1 = ImLerp(ps.inner_bb.Min, ps.inner_bb.Max, ImVec2(nx, 1.f - ny));
        window->DrawList->AddLine(pos0, pos1, line_color);

        prev_c = next_c;
    }

    if (GetActiveID() == ps.id || IsItemHovered()) {
        ImGuiContext& ctx = *GImGui;
        float t = ImClamp((ctx.IO.MousePos.x - ps.inner_bb.Min.x) / (ps.inner_bb.Max.x - ps.inner_bb.Min.x), 0.f, 1.f);
        int i = ImClamp((int)ImLerp(ps.coord_view.Min.x, ps.coord_view.Max.x, t), 0, count - 1);
        float y = values[i];

        BeginTooltip();
        Text("%i: %g\n", i, y);
        EndTooltip();
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
    if (!ItemAdd(total_bb, NULL)) return false;

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
    if (!ItemAdd(total_bb, NULL)) return false;

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
