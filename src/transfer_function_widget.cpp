#include <transfer_function_widget.h>

#include <imgui.h>
#include <implot.h>

#include <math.h>
#include <stdio.h>

void transfer_function_bake_rgba(const TransferFunction& tf, float* rgba, int n) {
    if (!rgba || n <= 0) return;
    for (int i = 0; i < n; ++i) {
        const float t = n > 1 ? (float)i / (float)(n - 1) : 0.0f;
        const ImVec4 c = ImPlot::SampleColormap(t, tf.colormap);
        rgba[4 * i + 0] = c.x;
        rgba[4 * i + 1] = c.y;
        rgba[4 * i + 2] = c.z;
    }
    transfer_function_bake_alpha(tf, rgba + 3, n, 4);
}

static inline ImU32 with_alpha(ImVec4 c, float a) {
    c.w = a;
    return ImGui::ColorConvertFloat4ToU32(c);
}

static void fit_view(TransferFunctionView* view, const TransferFunction& tf, float data_min, float data_max) {
    double lo = fmin(tf.range_beg, tf.range_end), hi = fmax(tf.range_beg, tf.range_end);
    if (view->fit == 2 && data_max > data_min) {
        lo = data_min;
        hi = data_max;
    }
    // and a little either side: the range's lines off the plot's edges, where they could not be taken hold of
    const double m = 0.1 * (hi - lo);
    lo -= m;
    hi += m;
    if (!(hi > lo)) hi = lo + 1.0;
    view->x_min = lo;
    view->x_max = hi;
}

bool transfer_function_editor(const char* id, TransferFunction* tf, TransferFunctionView* view,
                              const TransferFunctionHistogram* hist, float data_min, float data_max,
                              bool colors_only, float width, float height) {
    if (!tf || !view) return false;
    bool changed = false;
    ImGui::PushID(id);
    if (width <= 0.0f) width = ImGui::GetContentRegionAvail().x;

    const bool refit = view->fit != 0;
    if (refit) {
        fit_view(view, *tf, data_min, data_max);
        view->fit = 0;
    }

    // the opacity as a vertical slider beside the plot, its scale the curve's: the y axis is labelled in opacity
    const float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
    const float vslider_w = colors_only ? 0.0f : ImGui::GetFrameHeight();
    const float plot_w = colors_only ? width : width - vslider_w - spacing;

    const ImPlotFlags plot_flags = ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText;
    if (ImPlot::BeginPlot("##tf", ImVec2(plot_w, height), plot_flags)) {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoHighlight,
                          ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoHighlight | (colors_only ? ImPlotAxisFlags_NoTickLabels : 0));
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.06, ImPlotCond_Always);
        if (!colors_only) {
            // the curve's 0, 1/2 and 1 as the opacities they are
            char l0[16], l1[16], l2[16];
            snprintf(l0, sizeof(l0), "0");
            snprintf(l1, sizeof(l1), "%.3g", 0.5f * tf->opacity);
            snprintf(l2, sizeof(l2), "%.3g", tf->opacity);
            const double ticks[3] = {0.0, 0.5, 1.0};
            const char* labels[3] = {l0, l1, l2};
            ImPlot::SetupAxisTicks(ImAxis_Y1, ticks, 3, labels);
        }
        ImPlot::SetupAxisLimits(ImAxis_X1, view->x_min, view->x_max, refit ? ImPlotCond_Always : ImPlotCond_Once);
        ImPlot::SetupFinish();

        const ImPlotRect lim = ImPlot::GetPlotLimits();
        view->x_min = lim.X.Min;
        view->x_max = lim.X.Max;
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();

        // the histogram: log of the counts, so what few values there are away from the bulk still shows
        if (hist && hist->counts && hist->bins > 0 && hist->hi > hist->lo) {
            float cmax = 0.0f;
            for (int i = 0; i < hist->bins; ++i) cmax = fmaxf(cmax, hist->counts[i]);
            const float norm = cmax > 0.0f ? 1.0f / logf(1.0f + cmax) : 0.0f;
            const double w = ((double)hist->hi - hist->lo) / hist->bins;
            const ImU32 col = ImGui::GetColorU32(ImGuiCol_Text, 0.22f);   // neutral: the colour map is the colour
            for (int i = 0; i < hist->bins; ++i) {
                if (hist->counts[i] <= 0.0f) continue;
                const double h = logf(1.0f + hist->counts[i]) * norm;
                const ImVec2 a = ImPlot::PlotToPixels(hist->lo + i * w, 0.0);
                const ImVec2 b = ImPlot::PlotToPixels(hist->lo + (i + 1) * w, h);
                dl->AddRectFilled(ImVec2(a.x, b.y), ImVec2(fmaxf(b.x, a.x + 1.0f), a.y), col);
            }
        }

        // the colours: under the curve, and as a strip along the bottom (they are there where the curve is 0 too)
        {
            const int slices = 160;
            for (int i = 0; i < slices; ++i) {
                const double x0 = lim.X.Min + (lim.X.Max - lim.X.Min) * i / slices;
                const double x1 = lim.X.Min + (lim.X.Max - lim.X.Min) * (i + 1) / slices;
                const float t = transfer_function_param(*tf, (float)(0.5 * (x0 + x1)));
                const ImVec4 c = ImPlot::SampleColormap(t, tf->colormap);
                const float a = colors_only ? 1.0f : transfer_function_curve(*tf, t);
                const ImVec2 p0 = ImPlot::PlotToPixels(x0, 0.0);
                const ImVec2 p1 = ImPlot::PlotToPixels(x1, (double)a);
                const ImVec2 s1 = ImPlot::PlotToPixels(x1, 0.045);
                if (a > 0.0f) dl->AddRectFilled(ImVec2(p0.x, p1.y), ImVec2(p1.x + 0.5f, p0.y), with_alpha(c, colors_only ? 0.55f : 0.75f));
                dl->AddRectFilled(ImVec2(p0.x, s1.y), ImVec2(p1.x + 0.5f, p0.y), with_alpha(c, 1.0f));
            }
        }

        // the curve, flat beyond its first and last points to the edges of the view
        if (!colors_only) {
            ImVec2 pts[TRANSFER_FUNCTION_MAX_POINTS + 2];
            int n = 0;
            const float first = tf->num_points > 0 ? tf->points[0].alpha : 0.0f;
            const float last = tf->num_points > 0 ? tf->points[tf->num_points - 1].alpha : 0.0f;
            const bool forward = tf->range_end >= tf->range_beg;
            pts[n++] = ImPlot::PlotToPixels(lim.X.Min, forward ? first : last);
            for (int k = 0; k < tf->num_points; ++k) {
                const int i = forward ? k : tf->num_points - 1 - k;
                pts[n++] = ImPlot::PlotToPixels(transfer_function_value(*tf, tf->points[i].t), tf->points[i].alpha);
            }
            pts[n++] = ImPlot::PlotToPixels(lim.X.Max, forward ? last : first);
            dl->AddPolyline(pts, n, IM_COL32(0, 0, 0, 140), 0, 3.0f);
            dl->AddPolyline(pts, n, IM_COL32(255, 255, 255, 230), 0, 1.5f);
        }
        ImPlot::PopPlotClipRect();

        // the points: dragged, removed by a right click; one added by ctrl + click elsewhere. Before the range's lines:
        // the first item under the mouse is the one it takes, and a point on a line is to be taken, not the line
        bool point_hovered = false;
        int remove = -1;
        const float ext = tf->range_end - tf->range_beg;
        for (int i = 0; !colors_only && i < tf->num_points; ++i) {
            double x = transfer_function_value(*tf, tf->points[i].t);
            double y = tf->points[i].alpha;
            bool hovered = false;
            if (ImPlot::DragPoint(100 + i, &x, &y, ImVec4(1, 1, 1, 1), 5.0f, ImPlotDragToolFlags_None, nullptr, &hovered)) {
                const float t = fabsf(ext) > 0.0f ? (float)((x - tf->range_beg) / ext) : tf->points[i].t;
                transfer_function_move(tf, i, t, (float)y);
                changed = true;
            }
            // under the mouse by its pixels: in the frame of a right click the plot holds the mouse, and the point's
            // own hover is false then
            const ImVec2 px = ImPlot::PlotToPixels(x, y);
            const ImVec2 mouse = ImGui::GetMousePos();
            const bool under = (px.x - mouse.x) * (px.x - mouse.x) + (px.y - mouse.y) * (px.y - mouse.y) <= 7.0f * 7.0f;
            if (hovered || under) {
                point_hovered = true;
                ImGui::SetTooltip("%.4g: %.0f%%  (right click to remove)", transfer_function_value(*tf, tf->points[i].t), 100.0f * tf->points[i].alpha);
                if (under && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) remove = i;
            }
        }
        // the range: dragged as lines; symmetric, the other follows mirrored
        {
            double beg = tf->range_beg, end = tf->range_end;
            const ImVec4 col = ImVec4(1.0f, 1.0f, 1.0f, 0.8f);
            const bool db = ImPlot::DragLineX(1, &beg, col, 1.5f);
            const bool de = ImPlot::DragLineX(2, &end, col, 1.5f);
            if (db || de) {
                if (tf->symmetric) {
                    const float r = (float)fabs(db ? beg : end);
                    tf->range_beg = -r;
                    tf->range_end = r;
                } else {
                    tf->range_beg = (float)beg;
                    tf->range_end = (float)end;
                }
                changed = true;
            }
        }

        // the points over the lines, drawn again
        if (!colors_only) {
            ImPlot::PushPlotClipRect();
            for (int i = 0; i < tf->num_points; ++i) {
                const ImVec2 px = ImPlot::PlotToPixels(transfer_function_value(*tf, tf->points[i].t), tf->points[i].alpha);
                dl->AddCircleFilled(px, 6.5f, IM_COL32(0, 0, 0, 120), 16);
                dl->AddCircleFilled(px, 5.0f, IM_COL32(255, 255, 255, 255), 16);
            }
            ImPlot::PopPlotClipRect();
        }
        if (remove >= 0) changed |= transfer_function_remove(tf, remove);
        // (ctrl: ImPlot leaves the mouse alone while it is held, where a double click is its fitting of the axes)
        if (!colors_only && !point_hovered && ImPlot::IsPlotHovered() && ImGui::GetIO().KeyCtrl && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && fabsf(ext) > 0.0f) {
            const ImPlotPoint m = ImPlot::GetPlotMousePos();
            const float t = (float)((m.x - tf->range_beg) / ext);
            if (t >= 0.0f && t <= 1.0f) changed |= transfer_function_insert(tf, t, (float)m.y) >= 0;
        }
        // How it is used, once the mouse rests on the plot (not while it is being worked: it would be in the way)
        static float rest = 0.0f;
        const ImGuiIO& io = ImGui::GetIO();
        const bool resting = ImPlot::IsPlotHovered() && io.MouseDelta.x == 0.0f && io.MouseDelta.y == 0.0f && !ImGui::IsAnyMouseDown();
        rest = resting ? rest + io.DeltaTime : 0.0f;
        if (rest > 1.0f) {
            ImGui::SetTooltip(colors_only ? "The colour map over its range, the values' histogram behind.\n"
                                            "Drag the two lines of the range. Scroll to zoom, drag to pan."
                                          : "The opacity over the values, the colour map's colours under it, the values' histogram behind.\n"
                                            "Drag the points and the two lines of the range; the opacity by the slider beside it.\n"
                                            "Ctrl + click to add a point, right click one to remove it.\n"
                                            "Scroll to zoom, drag to pan.");
        }
        ImPlot::EndPlot();
    }
    if (!colors_only) {
        ImGui::SameLine(0.0f, spacing);
        float opacity = tf->opacity;
        if (ImGui::VSliderFloat("##opacity", ImVec2(vslider_w, height), &opacity, 0.001f, 1.0f, "", ImGuiSliderFlags_Logarithmic)) {
            tf->opacity = opacity;
            changed = true;
        }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip) || ImGui::IsItemActive()) {
            ImGui::SetTooltip("opacity %.3f: at the top of the curve (its 1), of a slab 1/150 of the volume thick", tf->opacity);
        }
    }

    // under the plot: the colour map, the range, symmetric, fitting the view
    const float cmap_w = width;
    if (ImPlot::ColormapButton(ImPlot::GetColormapName(tf->colormap), ImVec2(cmap_w, 0), tf->colormap)) {
        ImGui::OpenPopup("##colormaps");
    }
    if (ImGui::BeginPopup("##colormaps")) {
        for (int i = 0; i < ImPlot::GetColormapCount(); ++i) {
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(i), ImVec2(fminf(cmap_w, 240.0f), 0), i)) {
                tf->colormap = i;
                changed = true;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    const float sym_w = ImGui::CalcTextSize("symmetric").x + ImGui::GetFrameHeight() + spacing;
    const float fit_w = ImGui::CalcTextSize("fit").x + 2.0f * ImGui::GetStyle().FramePadding.x;
    ImGui::SetNextItemWidth(width - sym_w - fit_w - 2.0f * spacing);
    float range[2] = {tf->range_beg, tf->range_end};
    if (ImGui::InputFloat2("##range", range, "%.4g")) {
        // any range but none: one of no extent maps every value to the middle
        if (range[0] == range[1]) range[1] = range[0] + 1.0e-6f * fmaxf(1.0f, fabsf(range[0]));
        transfer_function_set_range(tf, range[0], range[1]);
        changed = true;
    }
    ImGui::SetItemTooltip("The values the colour map runs over, from its first colour to its last");
    ImGui::SameLine(0.0f, spacing);
    if (ImGui::Checkbox("symmetric", &tf->symmetric)) {
        if (tf->symmetric) transfer_function_set_range(tf, tf->range_beg, tf->range_end);
        changed = true;
    }
    ImGui::SetItemTooltip("The range as -r to r: for a signed field, the middle of the colour map at 0");
    ImGui::SameLine(0.0f, spacing);
    if (ImGui::Button("fit")) ImGui::OpenPopup("##fit");
    ImGui::SetItemTooltip("Fit the view or the range");
    if (ImGui::BeginPopup("##fit")) {
        if (ImGui::Selectable("View to the range")) view->fit = 1;
        if (ImGui::Selectable("View to the data", false, data_max > data_min ? 0 : ImGuiSelectableFlags_Disabled)) view->fit = 2;
        if (ImGui::Selectable("Range to the data", false, data_max > data_min ? 0 : ImGuiSelectableFlags_Disabled)) {
            transfer_function_set_range(tf, data_min, data_max);
            view->fit = 1;
            changed = true;
        }
        if (ImGui::Selectable("Range to the view")) {
            transfer_function_set_range(tf, (float)view->x_min, (float)view->x_max);
            changed = true;
        }
        ImGui::EndPopup();
    }

    ImGui::PopID();
    return changed;
}
