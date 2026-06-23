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

// Plot an implot line spline using a centripetal Catmull-Rom spline through the control points.
IMPLOT_API void PlotLineSpline(const char* label_id, const double* xs, const double* ys, int count) {
    if (ImPlot::BeginItem(label_id)) {
        ImPlot::PushPlotClipRect();
        ImDrawList& DrawList = *ImPlot::GetPlotDrawList();
        const ImPlotNextItemData& s = ImPlot::GetItemData();
        ImVec4 color = s.Colors[ImPlotCol_Line];

        for (int i = 0; i < count; i++) {
            if (ImPlot::FitThisFrame()) {
                ImPlot::FitPoint(ImPlotPoint(xs[i], 0.0));
                ImPlot::FitPoint(ImPlotPoint(xs[i], ys[i]));
            }
        }

        if (count >= 2) {
            const int   segs_per_span = 20;
            const float line_weight   = s.LineWeight;
            const ImU32 col32         = ImGui::ColorConvertFloat4ToU32(color);

            // Uniform-X natural cubic spline (C2, interpolates all points).
            // M[i] = second derivative at each knot; natural BCs fix M[0]=M[n-1]=0.
            // Uniform spacing reduces the tridiagonal to diagonals {1,4,1}.
            // M[] doubles as the RHS during the forward sweep, then is overwritten
            // in-place during back substitution — no separate rhs vector needed.
            const int    n = count;
            const double h = (xs[n-1] - xs[0]) / (double)(n - 1);

            ImVector<double> M, diag;
            M.resize(n, 0.0);

            if (n > 2) {
                diag.resize(n, 0.0);
                const double r_scale = 6.0 / (h * h);

                // Forward sweep: eliminate sub-diagonal, accumulate RHS into M[]
                diag[1] = 4.0;
                M[1]    = r_scale * (ys[2] - 2.0*ys[1] + ys[0]);
                for (int i = 2; i <= n - 2; ++i) {
                    double f = 1.0 / diag[i-1];
                    diag[i] = 4.0 - f;
                    M[i]    = r_scale * (ys[i+1] - 2.0*ys[i] + ys[i-1]) - f * M[i-1];
                }

                // Back substitution: M[] becomes the final second derivatives
                M[n-2] /= diag[n-2];
                for (int i = n - 3; i >= 1; --i)
                    M[i] = (M[i] - M[i+1]) / diag[i];
            }

            // Evaluate each segment: S(t) = y0 + b*t + c*t² + d*t³,  t ∈ [0, h]
            for (int seg = 0; seg < n - 1; ++seg) {
                double m0 = M[seg], m1 = M[seg+1];
                double x0 = xs[seg], y0 = ys[seg], y1 = ys[seg+1];
                double b  = (y1 - y0)/h - h/6.0 * (2.0*m0 + m1);
                double c  = m0 / 2.0;
                double d  = (m1 - m0) / (6.0 * h);

                const int step_start = (seg == 0) ? 0 : 1;
                for (int step = step_start; step <= segs_per_span; ++step) {
                    double t = (double)step / (double)segs_per_span * h;
                    DrawList.PathLineTo(ImPlot::PlotToPixels(x0 + t, y0 + t*(b + t*(c + t*d))));
                }
            }
            DrawList.PathStroke(col32, 0, line_weight);
        }

        ImPlot::PopPlotClipRect();
        ImPlot::EndItem();
    }
}

}  // namespace ImGui
