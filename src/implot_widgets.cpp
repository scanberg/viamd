#include <implot_widgets.h>
#include <implot_internal.h>

namespace ImPlot {

IMPLOT_API bool DragRangeX(const char* id, double* x_range_min, double* x_range_max, double min_value, double max_value, ImPlotDragRangeFlags flags, const ImPlotDragRangeStyle& style) {
    ImPlotContext& gp = *GImPlot;
    IM_ASSERT_USER_ERROR(gp.CurrentPlot != NULL, "DragRangeX() needs to be called between BeginPlot() and EndPlot()!");
    const float grab_size = ImMax(5.0f, style.line_thickness);
    float yt = gp.CurrentPlot->PlotRect.Min.y;
    float yb = gp.CurrentPlot->PlotRect.Max.y;
    float x_min = IM_ROUND(PlotToPixels(*x_range_min,0).x);
    float x_max = IM_ROUND(PlotToPixels(*x_range_max,0).x);
    const bool outside = x_max < (gp.CurrentPlot->PlotRect.Min.x - grab_size / 2) || x_min > (gp.CurrentPlot->PlotRect.Max.x + grab_size / 2);
    if (outside)
        return false;

    int dragging = 0;
    bool active = false;
    bool hovered = false;

    ImGui::PushID(id);

    //double* x_range[2] = {x_range_max, x_range_min};
    if ((gp.CurrentPlot->PlotRect.Min.x - grab_size / 2) < x_min && x_min < (gp.CurrentPlot->PlotRect.Max.x + grab_size / 2))  {
        if (DragLineX(1, x_range_min, style.line_col, style.line_thickness, flags)) dragging = 1;
        *x_range_min = ImClamp(*x_range_min, min_value, max_value);
        active  |= ImGui::GetActiveID()  == 1;
        hovered |= ImGui::GetHoveredID() == 1;
    }
    if ((gp.CurrentPlot->PlotRect.Min.x - grab_size / 2) < x_max && x_max < (gp.CurrentPlot->PlotRect.Max.x + grab_size / 2))  {
        if (DragLineX(2, x_range_max, style.line_col, style.line_thickness, flags)) dragging = 2;
        *x_range_max = ImClamp(*x_range_max, min_value, max_value);
        active  |= ImGui::GetActiveID()  == 2;
        hovered |= ImGui::GetHoveredID() == 2;
    }

    //float len = gp.Style.MajorTickLen.x;
    const float drag_bar_height = 15;
    const bool bar_enabled = !(flags & ImPlotDragRangeFlags_NoBar);

    if (bar_enabled) {
        if (!(ImGui::GetItemFlags() & ImGuiItemFlags_Disabled)) {
            if (!gp.CurrentPlot->Selecting) {
                if (dragging) {
                    if (dragging == 2) {
                        if (*x_range_min > *x_range_max) {
                            *x_range_min = *x_range_max;
                        }
                    } else {
                        if (*x_range_max < *x_range_min) {
                            *x_range_max = *x_range_min;
                        }
                    }
                }
                else {
                    ImVec2 old_cursor_pos = ImGui::GetCursorScreenPos();
                    ImVec2 new_cursor_pos = ImVec2(x_min, yb - drag_bar_height);
                    ImGui::GetCurrentWindow()->DC.CursorPos = new_cursor_pos;
                    const float btn_width = ImMax(x_max - x_min, 1.0f);
                    ImGui::InvisibleButton(id, ImVec2(btn_width, drag_bar_height));
                    ImGui::GetCurrentWindow()->DC.CursorPos = old_cursor_pos;

                    hovered |= ImGui::IsItemHovered();
                    active  |= ImGui::IsItemActive();

                    if (ImGui::IsItemActive() && ImGui::IsMouseDragging(0)) {
                        const float delta = ImGui::GetIO().MouseDelta.x;
                        double x_min_temp = PlotToPixels(*x_range_min, 0).x;
                        double x_max_temp = PlotToPixels(*x_range_max, 0).x;

                        x_min_temp += delta;
                        x_max_temp += delta;

                        *x_range_min = PixelsToPlot((float)x_min_temp, 0).x;
                        *x_range_max = PixelsToPlot((float)x_max_temp, 0).x;

                        double x_diff = *x_range_max - *x_range_min;

                        if (*x_range_min < min_value) {
                            *x_range_min = min_value;
                            *x_range_max = *x_range_min + x_diff;
                        }
                        if (*x_range_max > max_value) {
                            *x_range_max = max_value;
                            *x_range_min = *x_range_max - x_diff;
                        }

                        dragging = 3;
                    }
                }
            }
        }
    }

    ImVec4 range_color = IsColorAuto(style.range_col) ? ImGui::GetStyleColorVec4(ImGuiCol_ScrollbarBg) : style.range_col;
    ImVec4 bar_color = ImGui::GetStyleColorVec4(ImGuiCol_ScrollbarGrab);

    range_color.w *= ImGui::GetStyle().Alpha;
    bar_color.w *= ImGui::GetStyle().Alpha;

    ImDrawList& DrawList = *GetPlotDrawList();
    PushPlotClipRect(0.0f);
    DrawList.AddRectFilled(ImVec2(x_min+1, yt), ImVec2(x_max-1, yb), ImGui::ColorConvertFloat4ToU32(range_color));
    if (bar_enabled) {
        DrawList.AddRectFilled(ImVec2(x_min+1, yb), ImVec2(x_max-1, yb-drag_bar_height), ImGui::ColorConvertFloat4ToU32(bar_color));
    }
    PopPlotClipRect();

    if (hovered || active) {
        gp.CurrentPlot->Hovered = false;
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        if (!ImHasFlag(flags, ImPlotDragRangeFlags_NoLabel)) {
            char min_buff[32];
            char max_buff[32];
            char buff[128];
            float x  = IM_ROUND(PlotToPixels((*x_range_max + *x_range_min) * 0.5, 0).x);
            LabelAxisValue(gp.CurrentPlot->XAxis(0), *x_range_min, min_buff, sizeof(min_buff));
            LabelAxisValue(gp.CurrentPlot->XAxis(0), *x_range_max, max_buff, sizeof(max_buff));
            snprintf(buff, sizeof(buff), "[%s, %s]", min_buff, max_buff);
            ImVec4 color = IsColorAuto(style.line_col) ? ImGui::GetStyleColorVec4(ImGuiCol_Text) : style.line_col;
            ImU32  col32 = ImGui::ColorConvertFloat4ToU32(color);
            gp.Annotations.Append(ImVec2(x,yb),ImVec2(0,0),col32,CalcTextColor(color),true,"%s = %s", id, buff);
        }
    }

    ImGui::PopID();

    return dragging;
}

// Text with a dark shadow, legible over whatever is plotted under it
static void inset_text(ImDrawList* dl, ImFont* font, float size, ImVec2 p, const char* text, ImU32 col) {
    dl->AddText(font, size, ImVec2(p.x + 1, p.y + 1), IM_COL32(0, 0, 0, 170), text);
    dl->AddText(font, size, p, col, text);
}

IMPLOT_API void DrawInsetAxes(ImPlotInsetFlags flags, ImVec2 clear_top_left, float font_scale) {
    ImPlotContext& gp = *GImPlot;
    IM_ASSERT_USER_ERROR(gp.CurrentPlot != nullptr, "DrawInsetAxes() needs to be called between BeginPlot() and EndPlot()!");
    SetupLock();
    ImPlotPlot& plot = *gp.CurrentPlot;
    const ImRect r = plot.PlotRect;

    ImFont* font = ImGui::GetFont();
    const float fs = ImGui::GetFontSize() * font_scale;
    const float name_fs = fs * 1.15f;
    const float pad = 3.0f;
    const float line_h = font->CalcTextSizeA(fs, FLT_MAX, 0.0f, "0").y;
    ImVec4 text = GetStyleColorVec4(ImPlotCol_AxisText);
    const ImU32 name_col = ImGui::ColorConvertFloat4ToU32(text);
    text.w *= 0.75f;
    const ImU32 tick_col = ImGui::ColorConvertFloat4ToU32(text);
    const bool ticks = !(flags & ImPlotInsetFlags_NoTickLabels);
    const bool names = !(flags & ImPlotInsetFlags_NoAxisLabels);

    auto has_ticks = [&](const ImPlotAxis& ax) {
        return ticks && ax.Enabled && ax.Ticker.TickCount() > 0;
    };
    auto name_of = [&](const ImPlotAxis& ax) -> const char* {
        return names && ax.Enabled && ax.LabelOffset != -1 ? plot.GetAxisLabel(ax) : nullptr;
    };
    auto label_width = [&](const ImPlotAxis& ax) {
        float w = 0.0f;
        for (int j = 0; j < ax.Ticker.TickCount(); ++j) {
            const ImPlotTick& tk = ax.Ticker.Ticks[j];
            if (tk.ShowLabel && tk.TextOffset != -1) w = ImMax(w, font->CalcTextSizeA(fs, FLT_MAX, 0.0f, ax.Ticker.GetText(j)).x);
        }
        return w;
    };

    // The rows and columns the axes take along the edges, from the outside in: what the others keep clear of
    float row_bottom = 0.0f, row_top = 0.0f, col_left = 0.0f, col_right = 0.0f;
    for (int i = 0; i < IMPLOT_NUM_X_AXES; ++i) {
        const ImPlotAxis& ax = plot.XAxis(i);
        if (!has_ticks(ax) && !name_of(ax)) continue;
        (ax.IsOpposite() ? row_top : row_bottom) += line_h + pad;
    }
    for (int i = 0; i < IMPLOT_NUM_Y_AXES; ++i) {
        const ImPlotAxis& ax = plot.YAxis(i);
        if (!has_ticks(ax)) continue;
        (ax.IsOpposite() ? col_right : col_left) += label_width(ax) + pad * 2;
    }
    const float top_clear = ImMax(clear_top_left.y, row_top);

    PushPlotClipRect();
    ImDrawList* dl = GetPlotDrawList();

    float off_bottom = 0.0f, off_top = 0.0f;
    for (int i = 0; i < IMPLOT_NUM_X_AXES; ++i) {
        const ImPlotAxis& ax = plot.XAxis(i);
        const char* name = name_of(ax);
        if (!has_ticks(ax) && !name) continue;
        const bool top = ax.IsOpposite();
        float& off = top ? off_top : off_bottom;
        const float y = top ? r.Min.y + off + pad : r.Max.y - off - line_h - pad;
        off += line_h + pad;
        // A top row starts right of what is kept clear
        const float x_min = r.Min.x + (top ? ImMax(col_left, clear_top_left.x) : col_left) + pad;
        float x_max = r.Max.x - col_right - pad;
        if (name) {
            const ImVec2 ns = font->CalcTextSizeA(name_fs, FLT_MAX, 0.0f, name);
            x_max -= ns.x;
            inset_text(dl, font, name_fs, ImVec2(x_max, y + line_h - ns.y), name, name_col);
            x_max -= pad * 2;
        }
        if (!has_ticks(ax)) continue;
        float last = -FLT_MAX;
        for (int j = 0; j < ax.Ticker.TickCount(); ++j) {
            const ImPlotTick& tk = ax.Ticker.Ticks[j];
            if (!tk.ShowLabel || tk.TextOffset == -1) continue;
            const char* text = ax.Ticker.GetText(j);
            const float w = font->CalcTextSizeA(fs, FLT_MAX, 0.0f, text).x;
            const float x = IM_ROUND(ax.PlotToPixels(tk.PlotPos) - w * 0.5f);
            if (x < x_min || x + w > x_max || x < last + pad * 2) continue;
            inset_text(dl, font, fs, ImVec2(x, y), text, tick_col);
            last = x + w;
        }
    }

    float off_left = 0.0f, off_right = 0.0f;
    for (int i = 0; i < IMPLOT_NUM_Y_AXES; ++i) {
        const ImPlotAxis& ax = plot.YAxis(i);
        const char* name = name_of(ax);
        if (!has_ticks(ax) && !name) continue;
        const bool right = ax.IsOpposite();
        const float w_col = has_ticks(ax) ? label_width(ax) : 0.0f;
        float& off = right ? off_right : off_left;
        const float x0 = right ? r.Max.x - off - pad - w_col : r.Min.x + off + pad;
        off += w_col + pad * 2;
        // Below what is kept clear at the top (on the left), above the rows of the x axes
        float y_min = r.Min.y + (right ? row_top : top_clear) + pad;
        const float y_max = r.Max.y - row_bottom - pad;
        if (name) {
            const ImVec2 ns = font->CalcTextSizeA(name_fs, FLT_MAX, 0.0f, name);
            inset_text(dl, font, name_fs, ImVec2(right ? x0 + w_col - ns.x : x0, y_min), name, name_col);
            y_min += ns.y + pad;
        }
        if (!has_ticks(ax)) continue;
        float last = FLT_MAX;   // The top of the label drawn before, going up
        for (int j = 0; j < ax.Ticker.TickCount(); ++j) {
            const ImPlotTick& tk = ax.Ticker.Ticks[j];
            if (!tk.ShowLabel || tk.TextOffset == -1) continue;
            const char* text = ax.Ticker.GetText(j);
            const float w = font->CalcTextSizeA(fs, FLT_MAX, 0.0f, text).x;
            const float y = IM_ROUND(ax.PlotToPixels(tk.PlotPos) - line_h * 0.5f);
            if (y < y_min || y + line_h > y_max) continue;
            if (last != FLT_MAX && y + line_h > last - pad && y < last + line_h + pad) continue;
            inset_text(dl, font, fs, ImVec2(right ? x0 + w_col - w : x0, y), text, tick_col);
            last = y;
        }
    }

    PopPlotClipRect();
}

IMPLOT_API bool ColormapSelection(const char* id, ImPlotColormap* idx, ImVec2 size) {
    IM_ASSERT(id);
    IM_ASSERT(idx);

    bool value = false;

    if (size.x == 0)
        size.x = ImGui::CalcItemWidth();

    ImGui::PushID(ImGui::GetID(id));
    if (ImPlot::ColormapButton(ImPlot::GetColormapName(*idx), size, *idx)) {
        ImGui::OpenPopup("Color Map Selector");
    }
    if (id && id[0] != '#' && id[1] != '#')
        ImGui::Text("%s", id);
    if (ImGui::BeginPopup("Color Map Selector")) {
        for (int map = 0; map < ImPlot::GetColormapCount(); ++map) {
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), size, map)) {
                *idx = map;
                value = true;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    ImGui::PopID();

    return value;
}

IMPLOT_API bool ColormapQualitative(ImPlotColormap idx) {
    return ImPlot::GetCurrentContext()->ColormapData.IsQual(idx);
}

IMPLOT_API void SyncAxesY(double padding_fraction) {
    ImPlotPlot& plot = *GImPlot->CurrentPlot;
    for (int s = ImAxis_Y1; s < ImAxis_COUNT; ++s)
    {
        ImPlotAxis& axis = plot.Axes[s];
        double v = ImMax(fabs(axis.FitExtents.Min), fabs(axis.FitExtents.Max));
        v += v * ImMax(0.0, padding_fraction);
        axis.FitExtents.Min = -v;
        axis.FitExtents.Max = v;
    }
}

IMPLOT_API void SyncAxesWithPadding(int master_axis, int aux_axis, double aux_to_master, double padding_fraction) {
    IM_ASSERT(GImPlot && GImPlot->CurrentPlot);
    ImPlotContext& gp = *GImPlot;
    ImPlotPlot& plot = *gp.CurrentPlot;

    // Guard
    if (aux_to_master == 0.0) return;
    padding_fraction = ImMax(0.0, padding_fraction);

    // Read current fit extents (ImPlot computed these while plotting)
    double mmin = plot.Axes[master_axis].FitExtents.Min;
    double mmax = plot.Axes[master_axis].FitExtents.Max;
    double amin = plot.Axes[aux_axis].FitExtents.Min;
    double amax = plot.Axes[aux_axis].FitExtents.Max;

    double mspan = mmax - mmin;
    double aspan = amax - amin;

    // If spans are degenerate, fall back to small epsilon to avoid collapse
    const double EPS = 1e-12;
    if (mspan <= EPS && aspan <= EPS) return;

    // Convert aux span into master units
    double aspan_master = aspan * aux_to_master;

    // Choose the larger span and apply padding
    double target_span_master = ImMax(mspan, aspan_master);
    target_span_master *= (1.0 + padding_fraction);

    // Compute centers (preserve original centers rather than forcing symmetric about zero)
    double mcenter = 0.5 * (mmax + mmin);
    double acenter = 0.5 * (amax + amin);

    // New master extents centered on original master center
    double new_mmin = mcenter - 0.5 * target_span_master;
    double new_mmax = mcenter + 0.5 * target_span_master;

    // Convert target span back to aux units and center around original aux center
    double target_span_aux = target_span_master / aux_to_master;
    double new_amin = acenter - 0.5 * target_span_aux;
    double new_amax = acenter + 0.5 * target_span_aux;

    // Apply to plot (ImPlot reads FitExtents for layout)
    plot.Axes[master_axis].FitExtents.Min = new_mmin;
    plot.Axes[master_axis].FitExtents.Max = new_mmax;
    plot.Axes[aux_axis].FitExtents.Min = new_amin;
    plot.Axes[aux_axis].FitExtents.Max = new_amax;
}

// Plot an implot line spline using cubic spline.
IMPLOT_API void PlotLineSpline(const char* label_id, const double* xs, const double* ys, int count) {
    if (ImPlot::BeginItem(label_id)) {
        ImPlot::PushPlotClipRect();
        ImDrawList& DrawList = *ImPlot::GetPlotDrawList();
        const ImPlotNextItemData& s = ImPlot::GetItemData();
        ImVec4 color = s.Colors[ImPlotCol_Line];

        for (int i = 0; i < count; i++) {
            if (ImPlot::FitThisFrame()) {
                ImPlot::FitPoint(ImPlotPoint(xs[i], ys[i]));
            }
        }

        // Get explicit rect to test points against for skipping some work
        ImRect rect = ImPlot::GetCurrentPlot()->PlotRect;

        if (count >= 2) {
            const float target_pixel_spacing = 4.0;
            const float line_weight   = s.LineWeight;
            const ImU32 col32         = ImGui::ColorConvertFloat4ToU32(color);

            const int n = count;

            // Record the visible segment span in x
            int seg_min = INT_MAX;
            int seg_max = INT_MIN;

            // Points in pixel space
            ImVector<ImVec2> p;
            p.resize(n);
            for (int i = 0; i < n; ++i) {
                p[i] = ImPlot::PlotToPixels(xs[i], ys[i]);
                if (p[i].x > rect.Min.x) {
                    seg_min = ImMin(seg_min, i);
                }
                if (p[i].x < rect.Max.x) {
                    seg_max = ImMax(seg_max, i);
                }
            }

            ImVector<float> M, diag;
            M.resize(n, 0.0f);
            diag.resize(n, 0.0f);

            // Natural cubic spline:
            // M[0] = M[n-1] = 0
            // Solve for interior second derivatives
            if (n > 2) {
                ImVector<float> lower, upper;
                lower.resize(n, 0.0f);
                upper.resize(n, 0.0f);

                for (int i = 1; i < n - 1; ++i) {
                    const float h0 = p[i].x - p[i-1].x;
                    const float h1 = p[i+1].x - p[i].x;

                    lower[i] = h0;
                    diag[i]  = 2.0f * (h0 + h1);
                    upper[i] = h1;

                    M[i] = 6.0f * ((p[i+1].y - p[i].y) / h1 -
                                   (p[i].y   - p[i-1].y) / h0);
                }

                // Thomas algorithm forward elimination
                for (int i = 2; i < n - 1; ++i) {
                    const float w = lower[i] / diag[i-1];
                    diag[i] -= w * upper[i-1];
                    M[i]    -= w * M[i-1];
                }

                // Back substitution
                M[n-2] /= diag[n-2];
                for (int i = n - 3; i > 0; --i)
                    M[i] = (M[i] - upper[i] * M[i+1]) / diag[i];
            }

            // Evaluate spline segments
            // We need to pad the min and max with 1 extra segment to ensure we draw the entire visible region
            int seg_beg = ImMax(seg_min - 1, 0);
            int seg_end = ImMin(seg_max + 1, n - 1);

            for (int seg = seg_beg; seg < seg_end; ++seg) {
                const float x0 = p[seg].x;
                const float x1 = p[seg+1].x;
                const float y0 = p[seg].y;
                const float y1 = p[seg+1].y;

                const float h = x1 - x0;

                const float m0 = M[seg];
                const float m1 = M[seg+1];

                // Cubic coefficients:
                // S(t)=y0 + b*t + c*t² + d*t³, t=[0,h]
                const float b = (y1-y0)/h - h*(2.0f*m0+m1)/6.0f;
                const float c = m0/2.0f;
                const float d = (m1-m0)/(6.0f*h);

                const int step_start = (seg == 0) ? 0 : 1;

                const int num_steps = ImClamp(int(h / target_pixel_spacing), 4, 128);

                for (int step = step_start; step <= num_steps; ++step) {
                    const float t = (float)step / num_steps * h;
                    DrawList.PathLineTo(ImVec2(x0 + t, y0 + t*(b + t*(c + t*d))));
                }
            }

            DrawList.PathStroke(col32, 0, line_weight);
        }

        ImPlot::PopPlotClipRect();
        ImPlot::EndItem();
    }
}

}  // namespace ImGui
