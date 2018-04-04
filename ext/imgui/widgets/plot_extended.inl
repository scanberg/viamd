namespace ImGui {

int PlotExtended(ImGuiPlotType plot_type, const char* label, const float* data, int count,
                  int offset, int* selected_idx, const char* tooltip, std::function<std::string(int,float)> caption_func,
                  float scale_min, float scale_max, ImVec2 graph_size, int res,
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

    const ImRect frame_bb(window->DC.CursorPos,
                          window->DC.CursorPos + ImVec2(graph_size.x, graph_size.y));
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

    RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true,
                style.FrameRounding);
    
    if (count > 0 && graph_size.x > 0 && graph_size.y > 0) {
		auto get_data = [data, offset](int idx) -> float {
			return data[idx + offset];
        };

        int res_w =
            ImMin((int)graph_size.x / res, count) + ((plot_type == ImGuiPlotType_Lines) ? -1 : 0);
        int item_count = count + ((plot_type == ImGuiPlotType_Lines) ? -1 : 0);

        const ImU32 col_base = GetColorU32(
            (plot_type == ImGuiPlotType_Lines) ? ImGuiCol_PlotLines : ImGuiCol_PlotHistogram);
        const ImU32 col_hovered =
            GetColorU32((plot_type == ImGuiPlotType_Lines) ? ImGuiCol_PlotLinesHovered
                                                           : ImGuiCol_PlotHistogramHovered);

        const ImGuiID id = window->GetID(label);

        // Tooltip on hover
        int v_hovered = -1;
        if (IsItemHovered()) {
            const float t = ImClamp(
                (g.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.0f, 1.0f);
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
				float val = get_data(v_idx);
			}

            // Draw vertical-line
            const ImVec2 btm = ImVec2((int)(g.IO.MousePos.x + 0.5f), inner_bb.Min.y);
            const ImVec2 top = ImVec2((int)(g.IO.MousePos.x + 0.5f), inner_bb.Max.y);
            window->DrawList->AddLine(btm, top, col_base);

            v_hovered = v_idx;
        }

        // Draw zero axis
        if (plot_type == ImGuiPlotType_Lines && scale_min < 0.f && scale_max > 0.f) {
            auto t = 1.0f - ImSaturate((0.f - scale_min) / (scale_max - scale_min));
            auto y_pos = (int)(ImLerp(inner_bb.Min.y, inner_bb.Max.y, t) + 0.5f);
            auto color = IM_COL32(255, 255, 255, 140);
            window->DrawList->AddLine(ImVec2(inner_bb.Min.x, y_pos), ImVec2(inner_bb.Max.x, y_pos),
                                      color);
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
                while (current_idx < highlight_count && v0_idx <= highlighted_indices[current_idx] &&
                       highlighted_indices[current_idx] < v1_idx) {
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
        ImVec2 tp0 = ImVec2(
            t0,
            1.0f - ImSaturate(
                        (v0 - scale_min) /
                        (scale_max -
                        scale_min)));  // Point in the normalized space of our target rectangle

        for (int n = 0; n < res_w; n++) {
            const float t1 = t0 + t_step;
            const int v1_idx = ImClamp((int)(t1 * item_count + 0.5f), 0, count - 1);
            float v1 = avg_samples(v0_idx, v1_idx);
            const ImVec2 tp1 =
                ImVec2(t1, 1.0f - ImSaturate((v1 - scale_min) / (scale_max - scale_min)));

            // NB: Draw calls are merged together by the DrawList system. Still, we should
            // render our batch are lower level to save a bit of CPU.
            ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, tp0);
            ImVec2 pos1 =
                ImLerp(inner_bb.Min, inner_bb.Max,
                        (plot_type == ImGuiPlotType_Lines) ? tp1 : ImVec2(tp1.x, 1.0f));

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
                const float t = ImClamp(
                    (g.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x),
                    0.0f, 1.0f);
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
                const float x = (int)(ImLerp(inner_bb.Min.x, inner_bb.Max.x, t) + 0.5f);
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
		RenderTextClipped(ImVec2(frame_bb.Min.x, frame_bb.Min.y + style.FramePadding.y),
			frame_bb.Max, caption.c_str(), NULL, NULL, ImVec2(0.5f, 0.0f));
	}

    // Label
    RenderTextClipped(
        ImVec2(frame_bb.Min.x + style.FramePadding.x, frame_bb.Min.y + style.FramePadding.y),
        frame_bb.Max, label, NULL, NULL, ImVec2(0.0f, 0.0f));

	return value;
}

int PlotLinesExtended(const char* label, const float* values, int count, int offset, int* selected_idx,
	const char* tooltip, std::function<std::string(int,float)> caption_func, float scale_min, float scale_max, ImVec2 graph_size,
                       int res, const int* highlighted_indices, int highlight_count) {
    return PlotExtended(ImGuiPlotType_Lines, label, values, count, offset, selected_idx,
                 tooltip, caption_func, scale_min, scale_max, graph_size, res,
                 highlighted_indices, highlight_count);
}

int PlotHistogramExtended(const char* label, const float* values, int count, int offset, int* selected_idx,
                           const char* tooltip, std::function<std::string(int, float)> caption_func, float scale_min, float scale_max, ImVec2 graph_size,
                           int res, const int* highlighted_indices, int highlight_count) {
    return PlotExtended(ImGuiPlotType_Histogram, label, values, count, offset, selected_idx,
                 tooltip, caption_func, scale_min, scale_max, graph_size, res,
                 highlighted_indices, highlight_count);
}

static int ComputeIndex(float val, float val_min, float val_max, int idx_min, int idx_max) {
	float t = ImSaturate((val - val_min) / (val_max - val_min));
	return 0;
}

PlotFrame BeginPlotFrame(const char* label, ImVec2 frame_size, int offset, int count, float scale_min, float scale_max, std::function<std::string(int)> x_label_func, const int* highlight_indices, int highlight_count) {
	ImGuiWindow* window = GetCurrentWindow();

	if (window->SkipItems) return{};

	ImGuiContext& ctx = *GImGui;
	const ImGuiStyle& style = ctx.Style;

	PushItemWidth(-1);

	if (frame_size.x == 0.0f) frame_size.x = CalcItemWidth();
	if (frame_size.y == 0.0f) frame_size.y = (style.FramePadding.y * 2);

	const ImRect frame_bb(window->DC.CursorPos,
		window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
	const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
	const ImRect total_bb(frame_bb.Min, frame_bb.Max);
	ItemSize(total_bb, style.FramePadding.y);
	if (!ItemAdd(total_bb, NULL)) return{};

	RenderFrame(frame_bb.Min, frame_bb.Max, GetColorU32(ImGuiCol_FrameBg), true,
		style.FrameRounding);

	const ImGuiID id = window->GetID(label);

	// Draw zero axis
	if (scale_min < 0.f && scale_max > 0.f) {
		auto t = 1.0f - ImSaturate((0.f - scale_min) / (scale_max - scale_min));
		auto y_pos = (int)(ImLerp(inner_bb.Min.y, inner_bb.Max.y, t) + 0.5f);
		auto color = IM_COL32(255, 255, 255, 140);
		window->DrawList->AddLine(ImVec2(inner_bb.Min.x, y_pos), ImVec2(inner_bb.Max.x, y_pos),
			color);
	}

	// Draw highlighted indices
	if (highlight_count > 0) {
		const int res_w = ImMin((int)(inner_bb.Max.x - inner_bb.Min.x), count);
		const float t_step = 1.0f / (float)res_w;
		const ImU32 highlight_col = IM_COL32(255, 255, 255, 30);

		//float t0 = 0.f;
		//int l0_idx = 0;
		//int g0_idx = offset;
		//int v0_idx = offset;
		int current_idx = 0;

		for (int n = 0; n < res_w; n++) {
			//const float t1 = ImClamp(t0 + t_step, 0.f, 1.f);
			const int l0_idx = n * count / res_w;
			const int g0_idx = l0_idx + offset;

			const int l1_idx = (n + 1) * count / res_w;
			const int g1_idx = l1_idx + offset;

			//const int v1_idx = (n + 1) * count / res_w + offset; //ImClamp((int)(t1 * count), 0, count - 1) + offset;

			bool color_background = false;
			while (current_idx < highlight_count && g0_idx <= highlight_indices[current_idx] &&
				highlight_indices[current_idx] < g1_idx) {
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

			//t0 = t1;
			//v0_idx = v1_idx;
		}
	}

	std::string tooltip;

	int hovered_index = -1;
	if (IsItemHovered()) {
		const float t = ImSaturate((ctx.IO.MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x));
		hovered_index = ImClamp((int)(t * count), 0, count - 1) + offset;

		if (x_label_func) {
			tooltip += x_label_func(hovered_index) + '\n';
		}
		else {
			tooltip += std::to_string(hovered_index) + ":\n";
		}
	}

	return { id, label, hovered_index, count, offset, scale_min, scale_max, frame_bb.Min, frame_bb.Max, inner_bb.Min, inner_bb.Max, tooltip };
}

int EndPlotFrame(const PlotFrame& frame, int selected_idx) {
	ImGuiContext& ctx = *GImGui;
	ImGuiWindow* window = GetCurrentWindow();
	bool manhattan_mode = false;

	const ImU32 col_base = GetColorU32(ImGuiCol_PlotLines);

	if (frame.hovered_index > -1) {
		// Draw vertical-line
		const ImVec2 btm = ImVec2((int)(ctx.IO.MousePos.x + 0.5f), frame.inner_bb_min.y);
		const ImVec2 top = ImVec2((int)(ctx.IO.MousePos.x + 0.5f), frame.inner_bb_max.y);
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
			const float t = ImClamp(
				(ctx.IO.MousePos.x - frame.inner_bb_min.x) / (frame.inner_bb_max.x - frame.inner_bb_min.x),
				0.0f, 1.0f);
			const int v_idx = (int)(t * frame.count);
			value = ImClamp(v_idx, 0, frame.count - 1) + frame.offset;
		}
		else {
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
			const float x = (int)(ImLerp(frame.inner_bb_min.x, frame.inner_bb_max.x, t) + 0.5f);
			ImVec2 btm = ImVec2(x, frame.inner_bb_min.y);
			ImVec2 top = ImVec2(x, frame.inner_bb_max.y);
			window->DrawList->AddLine(btm, top, GetColorU32(ImGuiCol_PlotLinesHovered));
		}
	}

	const ImGuiStyle& style = ctx.Style;

	// Caption
    if (frame.label) {
        RenderTextClipped(ImVec2(frame.frame_bb_min.x, frame.frame_bb_min.y + style.FramePadding.y),
            frame.frame_bb_max, frame.label, NULL, NULL, ImVec2(0.0f, 0.0f));
    }

	return value;
}

void PlotFrameLine(PlotFrame& frame, const char* label, const float* values, FrameLineStyle style, int selected_idx) {

	if (values) {
		ImGuiWindow* window = GetCurrentWindow();

		float w = frame.inner_bb_max.x - frame.inner_bb_min.x;
		int res_w = ImMin((int)w, frame.count);
		//bool manhattan_mode = style.line_mode == LineMode_Histogram;
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
		ImVec2 tp0 = ImVec2(
			t0,
			1.0f - ImSaturate( (v0 - frame.scale_min) / (frame.scale_max - frame.scale_min)));  // Point in the normalized space of our target rectangle

		// Brighten the color in HSV and store result in hovered_color
		ImVec4 rgba = ImGui::ColorConvertU32ToFloat4(style.line_color);
		ImGui::ColorConvertRGBtoHSV(rgba.x, rgba.y, rgba.z, rgba.x, rgba.y, rgba.z);
		rgba.z = ImMin(rgba.z * 1.3f, 1.f);
		ImGui::ColorConvertHSVtoRGB(rgba.x, rgba.y, rgba.z, rgba.x, rgba.y, rgba.z);
		ImU32 hovered_color = ImGui::ColorConvertFloat4ToU32(rgba);

		float zero = ImLerp(frame.inner_bb_min.y, frame.inner_bb_max.y, 1.f - ImSaturate((0 - frame.scale_min) / (frame.scale_max - frame.scale_min)));

		for (int n = 0; n < res_w; n++) {
			const float t1 = t0 + t_step;
			const int v1_idx = ImClamp((int)(t1 * frame.count), 0, frame.count) + frame.offset;
			float v1 = avg_samples(v0_idx, v1_idx);
			const ImVec2 tp1 =
				ImVec2(t1, 1.0f - ImSaturate((v1 - frame.scale_min) / (frame.scale_max - frame.scale_min)));

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
			}
			else if (style.line_mode == LineMode_Line) {
				window->DrawList->AddLine(pos0, pos2, c);

				if (fill) {
					ImVec2 points[4]{
						{pos0.x, zero},
						pos0,
						pos2,
						{pos2.x, zero}
					};

					window->DrawList->AddConvexPolyFilled(points, 4, style.custom_fill_color);
				}
			}

			t0 = t1;
			tp0 = tp1;
			v0_idx = v1_idx;
		}

		if (frame.hovered_index > -1) {
			//const char* c = reinterpret_cast<char*>(&style.line_color);
			//std::string color_modifier{ '\033', c[0], c[1], c[2], c[3] };
			//frame.tooltip += color_modifier + label + ": " + std::to_string(values[frame.hovered_index]) + "\n";
			frame.tooltip += std::string(label) + ": " + std::to_string(values[frame.hovered_index]) + "\n";
		}

		if (selected_idx > -1) {
			// Draw line if in visible range
			if (frame.offset <= selected_idx && selected_idx < frame.offset + frame.count) {
				const float t = ImClamp((selected_idx - frame.offset) / (float)(frame.count), 0.0f, 1.0f);
				// cast to int for sharp edges due to dpi-scaling
				const float x = (int)(ImLerp(frame.inner_bb_min.x, frame.inner_bb_max.x, t) + 0.5f);
				ImVec2 btm = ImVec2(x, frame.inner_bb_min.y);
				ImVec2 top = ImVec2(x, frame.inner_bb_max.y);
				window->DrawList->AddLine(btm, top, style.line_color);
			}
		}
	}
}

}  // namespace ImGui