#include <color_scale.h>

#include <display_units.h>
#include <imgui_widgets.h>
#include <app/IconsFontAwesome6.h>

#include <core/md_common.h>

#include <imgui.h>
#include <implot.h>
#include <implot_internal.h>   // SampleColormapU32, the ticker of ColormapScale

#include <math.h>
#include <stdio.h>

uint32_t color_scale_color_u32(const ColorScale& scale, float value) {
    return ImPlot::SampleColormapU32(color_scale_param(scale, value), scale.colormap);
}

// The unit goes into a format string, so a '%' in it would be read as a conversion
static void escape_unit_for_format(char* dst, size_t cap, const char* unit_str) {
    if (cap == 0) return;
    dst[0] = '\0';
    for (size_t i = 0, j = 0; unit_str[i] && j + 2 < cap; ++i) {
        if (unit_str[i] == '%') dst[j++] = '%';
        dst[j++] = unit_str[i];
        dst[j] = '\0';
    }
}

// ---------------------------------------------------------------------------
// Isolines
// ---------------------------------------------------------------------------

int color_scale_isolines_add(ColorScaleIsolines* lines, float value) {
    ASSERT(lines);
    if (lines->count < 0 || lines->count >= ISOLINES_MAX) return -1;
    int i = lines->count;
    while (i > 0 && lines->values[i - 1] > value) {
        lines->values[i] = lines->values[i - 1];
        --i;
    }
    lines->values[i] = value;
    lines->count += 1;
    if (lines->drag >= i) lines->drag += 1;
    return i;
}

void color_scale_isolines_remove(ColorScaleIsolines* lines, int idx) {
    if (!lines || idx < 0 || idx >= lines->count) return;
    for (int i = idx; i + 1 < lines->count; ++i) {
        lines->values[i] = lines->values[i + 1];
    }
    lines->count -= 1;
    if (lines->drag == idx)     lines->drag = -1;
    else if (lines->drag > idx) lines->drag -= 1;
}

// In order of value again, after one was moved; returns where line 'idx' went
static int isolines_sort(ColorScaleIsolines* lines, int idx) {
    for (int i = 1; i < lines->count; ++i) {
        const float v = lines->values[i];
        const bool moving = idx == i;
        int j = i;
        while (j > 0 && lines->values[j - 1] > v) {
            lines->values[j] = lines->values[j - 1];
            if (!moving && idx == j - 1) idx = j;
            --j;
        }
        lines->values[j] = v;
        if (moving) idx = j;
    }
    return idx;
}

IsolineDesc color_scale_isoline_desc(const ColorScaleIsolines& lines, float px_per_point) {
    IsolineDesc desc;
    if (lines.show) {
        desc.count = (uint32_t)CLAMP(lines.count, 0, ISOLINES_MAX);
        for (uint32_t i = 0; i < desc.count; ++i) desc.values[i] = lines.values[i];
    }
    desc.hover       = lines.hover;
    desc.hover_value = lines.hover_value;
    desc.width       = MAX(lines.width, 0.25f) * MAX(px_per_point, 0.25f);
    return desc;
}

bool color_scale_draw_isoline_controls(ColorScaleIsolines* lines, const ColorScale& scale, md_unit_t unit) {
    ASSERT(lines);
    bool changed = false;
    ImGui::PushID("isolines");

    changed |= ImGui::Checkbox("isolines", &lines->show);
    ImGui::SetItemTooltip("Lines on the surface where the field takes given values.\n"
                          "On the legend: hover a value to see its line, click to add a line there,\n"
                          "drag a line to move it, right click to edit or remove one.");
    if (!lines->show) {
        ImGui::PopID();
        return changed;
    }

    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    const float fscl = (float)scl;

    // Under the checkbox and narrower by the indent, so that the labels stay in their column
    const ImGuiStyle& style = ImGui::GetStyle();
    const float indent = style.IndentSpacing;
    const float item_w = MAX(ImGui::CalcItemWidth() - indent, 60.0f);
    const float btn_w  = ImGui::GetFrameHeight();
    const float gap    = style.ItemInnerSpacing.x;
    ImGui::Indent(indent);

    ImGui::SetNextItemWidth(item_w);
    changed |= ImGui::SliderFloat("width", &lines->width, 0.5f, 6.0f, "%.1f pt");

    // The values in the display unit, which is said once, as the label of the list
    char label[64];
    if (unit_str[0]) snprintf(label, sizeof(label), "values (%s)", unit_str);
    else             snprintf(label, sizeof(label), "values");

    // A drag of the whole range over about 300 points
    const float ext   = fabsf(scale.range_end - scale.range_beg) * fscl;
    const float speed = ext > 0.0f ? ext / 300.0f : 0.001f;
    int remove = -1;
    for (int i = 0; i < lines->count; ++i) {
        ImGui::PushID(i);
        float v = lines->values[i] * fscl;
        ImGui::SetNextItemWidth(item_w - btn_w - gap);
        if (ImGui::DragFloat("##value", &v, speed, 0.0f, 0.0f, "%.4g")) {
            lines->values[i] = v / fscl;
            changed = true;
        }
        // Back in order once the value is set, not while it is being dragged past another
        if (ImGui::IsItemDeactivatedAfterEdit()) {
            isolines_sort(lines, i);
        }
        ImGui::SameLine(0.0f, gap);
        if (ImGui::DeleteButton(ICON_FA_XMARK, ImVec2(btn_w, 0))) {
            remove = i;
        }
        ImGui::SetItemTooltip("Remove the line");
        if (i == 0) {
            ImGui::SameLine(0.0f, gap);
            ImGui::TextUnformatted(label);
        }
        ImGui::PopID();
    }
    if (remove >= 0) {
        color_scale_isolines_remove(lines, remove);
        changed = true;
    }

    const float row_x = ImGui::GetCursorPosX();
    ImGui::BeginDisabled(lines->count >= ISOLINES_MAX);
    if (ImGui::Button(ICON_FA_PLUS, ImVec2(btn_w, 0))) {
        // Zero where the range holds it (the neutral contour of a potential), else the middle of the range
        const float lo = MIN(scale.range_beg, scale.range_end), hi = MAX(scale.range_beg, scale.range_end);
        const float v = (lo <= 0.0f && 0.0f <= hi) ? 0.0f : 0.5f * (lo + hi);
        changed |= color_scale_isolines_add(lines, v) >= 0;
    }
    ImGui::EndDisabled();
    if (lines->count >= ISOLINES_MAX) {
        ImGui::SetItemTooltip("%d lines at most", ISOLINES_MAX);
    } else {
        ImGui::SetItemTooltip("Add a line (or click the legend's colour bar where it should be)");
    }
    if (lines->count == 0) {
        // Nothing to list: the label goes with the button, in the column of the labels
        ImGui::SameLine(0.0f, gap);
        ImGui::TextDisabled("none yet");
        ImGui::SameLine(row_x + item_w + gap);
        ImGui::TextUnformatted(label);
    }

    ImGui::Unindent(indent);
    ImGui::PopID();
    return changed;
}

// Where the isolines of a legend take the mouse, by the legend's id: the frame of its scale as last drawn. Read
// before the legend's window is begun, to decide whether it takes input at all.
static ImGuiStorage legend_hit_storage;

static ImGuiID legend_hit_key(int id, int k) {
    const int key[2] = {id, k};
    return ImHashData(key, sizeof(key), 0x15011E5u);
}

static ImRect legend_hit_get(int id) {
    return ImRect(legend_hit_storage.GetFloat(legend_hit_key(id, 0), 0.0f), legend_hit_storage.GetFloat(legend_hit_key(id, 1), 0.0f),
                  legend_hit_storage.GetFloat(legend_hit_key(id, 2), 0.0f), legend_hit_storage.GetFloat(legend_hit_key(id, 3), 0.0f));
}

static void legend_hit_set(int id, const ImRect& r) {
    legend_hit_storage.SetFloat(legend_hit_key(id, 0), r.Min.x);
    legend_hit_storage.SetFloat(legend_hit_key(id, 1), r.Min.y);
    legend_hit_storage.SetFloat(legend_hit_key(id, 2), r.Max.x);
    legend_hit_storage.SetFloat(legend_hit_key(id, 3), r.Max.y);
}

// The line a popup of the legend edits (-1: none, the popup adds one at popup_value). One popup is open at a time.
static int   legend_popup_line  = -1;
static float legend_popup_value = 0.0f;

// The drag of a line on a legend (one at a time): how far the line was from the cursor when it was taken, so that
// one taken off-centre does not jump to the cursor, and whether it has moved since, so that a click on a line
// leaves it where it was
static float legend_drag_offset = 0.0f;
static float legend_drag_press_y = 0.0f;
static bool  legend_drag_moved = false;

// The isolines over the scale ImPlot::ColormapScale just drew at frame_min, frame_size: their markers, and with
// 'interactive' the hover, adding, dragging and editing
static void legend_isolines(ColorScaleIsolines* lines, const ColorScale& scale, double scl, const char* unit_str, int id,
                            ImVec2 frame_min, ImVec2 frame_size, bool interactive) {
    // The colour bar as ColormapScale lays it out (implot.cpp): inset by the plot padding, and as wide as its tick
    // labels leave, whose widths are in the ticker it has just filled
    const ImPlotStyle& ps = ImPlot::GetStyle();
    const float pad   = ps.LabelPadding.x + GImPlot->CTicker.MaxSize.x;
    const float bar_w = MAX(frame_size.x - (pad + 2.0f * ps.PlotPadding.x), ps.MajorTickLen.y);
    const ImRect frame(frame_min, ImVec2(frame_min.x + frame_size.x, frame_min.y + frame_size.y));
    const ImRect bar(ImVec2(frame_min.x + ps.PlotPadding.x, frame_min.y + ps.PlotPadding.y),
                     ImVec2(frame_min.x + ps.PlotPadding.x + bar_w, frame_min.y + frame_size.y - ps.PlotPadding.y));
    legend_hit_set(id, frame);

    const float beg = scale.range_beg, end = scale.range_end;
    const float bar_h = bar.GetHeight();
    if (!(bar_h > 1.0f) || beg == end) return;
    const float lo = MIN(beg, end), hi = MAX(beg, end);

    // range_beg at the bottom of the bar and range_end at the top, as ColormapScale draws them, also when the
    // range runs backwards
    auto y_of = [&](float v) { return bar.Max.y - (v - beg) / (end - beg) * bar_h; };
    auto v_of = [&](float y) {
        const float t = CLAMP((bar.Max.y - y) / bar_h, 0.0f, 1.0f);
        return beg + t * (end - beg);
    };

    // A value taken from the bar is rounded to the coarsest 1, 2 or 5 times a power of ten (in the display unit)
    // that is still finer than a pixel of the bar: what a click or a drag gives reads as a round number, and
    // nothing the bar can resolve is lost to the rounding
    double step = 0.0;
    {
        const double per_px = fabs((double)(end - beg) * scl) / (double)bar_h;
        if (per_px > 0.0) {
            const double p = pow(10.0, floor(log10(per_px)));
            const double m = per_px / p;
            step = (m >= 5.0 ? 5.0 : m >= 2.0 ? 2.0 : 1.0) * p;
        }
    }
    auto pick = [&](float y) {
        float v = v_of(y);
        if (step > 0.0) {
            v = (float)(round((double)v * scl / step) * step / scl);
            v = CLAMP(v, lo, hi);
        }
        return v;
    };

    char unit_fmt[64];
    escape_unit_for_format(unit_fmt, sizeof(unit_fmt), unit_str);
    char fmt[96];
    snprintf(fmt, sizeof(fmt), "%%.4g %s", unit_fmt);

    const ImVec2 mouse = ImGui::GetIO().MousePos;
    bool hovered = false;
    int  nearest = -1;

    if (interactive) {
        ImGui::SetCursorScreenPos(frame_min);
        ImGui::InvisibleButton("##isolines", frame_size, ImGuiButtonFlags_MouseButtonLeft);
        hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        const bool active = ImGui::IsItemActive();

        // The line whose marker is within reach of the cursor, the nearest if several are
        float best = 5.0f;
        for (int i = 0; i < lines->count; ++i) {
            const float v = lines->values[i];
            if (v < lo || v > hi) continue;
            const float d = fabsf(y_of(v) - mouse.y);
            if (d <= best) {
                best = d;
                nearest = i;
            }
        }

        // A press takes the line within reach, or adds one where it is, to be dragged into place
        if (ImGui::IsItemActivated()) {
            if (nearest >= 0) {
                lines->drag = nearest;
                legend_drag_offset = y_of(lines->values[nearest]) - mouse.y;
            } else {
                lines->drag = color_scale_isolines_add(lines, pick(mouse.y));
                legend_drag_offset = 0.0f;
            }
            legend_drag_press_y = mouse.y;
            legend_drag_moved = false;
            if (lines->drag >= 0) lines->show = true;
        }
        if (lines->drag >= 0) {
            if (active && lines->drag < lines->count) {
                legend_drag_moved |= mouse.y != legend_drag_press_y;
                if (legend_drag_moved) {
                    lines->values[lines->drag] = pick(mouse.y + legend_drag_offset);
                }
            } else {
                isolines_sort(lines, lines->drag);
                lines->drag = -1;
            }
        }

        // The cursor says what a press does: move the line within reach (or being moved), add one, or nothing when
        // there is no room for another
        if (lines->drag >= 0 || (hovered && nearest >= 0)) {
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeNS);
        } else if (hovered && lines->count >= ISOLINES_MAX) {
            ImGui::SetMouseCursor(ImGuiMouseCursor_NotAllowed);
        }

        if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            legend_popup_line  = nearest;
            legend_popup_value = pick(mouse.y);
            ImGui::OpenPopup("##isoline");
        }

        // What the cursor points at is drawn as the transient line, on the surface and here: the line being moved,
        // or the one within reach, or else the value under the cursor
        if (lines->drag >= 0) {
            lines->hover = true;
            lines->hover_value = lines->values[lines->drag];
        } else if (hovered) {
            lines->hover = true;
            lines->hover_value = nearest >= 0 ? lines->values[nearest] : pick(mouse.y);
        }
        if (lines->hover && !ImGui::IsPopupOpen("##isoline")) {
            ImGui::BeginTooltip();
            ImGui::Text(fmt, lines->hover_value * scl);
            if (lines->drag < 0) {
                if (nearest >= 0) {
                    ImGui::TextDisabled("drag to move, right click to edit");
                } else if (lines->count >= ISOLINES_MAX) {
                    ImGui::TextDisabled("%d lines at most", ISOLINES_MAX);
                } else {
                    ImGui::TextDisabled("click to add a line, right click for more");
                }
            }
            ImGui::EndTooltip();
        }
    } else if (lines->drag >= 0) {
        isolines_sort(lines, lines->drag);
        lines->drag = -1;
    }

    if (ImGui::BeginPopup("##isoline")) {
        const float fscl = (float)scl;
        if (legend_popup_line >= 0 && legend_popup_line < lines->count) {
            int& i = legend_popup_line;
            float v = lines->values[i] * fscl;
            const float speed = fabsf(end - beg) * fscl / 300.0f;
            ImGui::SetNextItemWidth(160.0f);
            if (ImGui::DragFloat("##value", &v, speed > 0.0f ? speed : 0.001f, 0.0f, 0.0f, fmt)) {
                lines->values[i] = v / fscl;
                lines->show = true;
            }
            if (ImGui::IsItemDeactivatedAfterEdit()) {
                i = isolines_sort(lines, i);
            }
            if (ImGui::MenuItem(ICON_FA_XMARK " Remove line")) {
                color_scale_isolines_remove(lines, i);
                i = -1;
            }
        } else {
            char item[128];
            snprintf(item, sizeof(item), ICON_FA_PLUS " Add line at %.4g %s", legend_popup_value * scl, unit_str);
            if (ImGui::MenuItem(item, nullptr, false, lines->count < ISOLINES_MAX)) {
                legend_popup_line = color_scale_isolines_add(lines, legend_popup_value);
                lines->show = true;
            }
        }
        ImGui::MenuItem("Show lines", nullptr, &lines->show);
        if (ImGui::MenuItem("Remove all lines", nullptr, false, lines->count > 0)) {
            lines->count = 0;
            lines->drag  = -1;
            legend_popup_line = -1;
        }
        ImGui::EndPopup();
    }

    // The markers: across the bar and as a handle to the left of it, drawn as the lines are on the surface - the
    // persistent ones dark with a light rim, the transient one the other way around. A line beyond the range keeps
    // its handle at the end it is beyond; hidden lines are faint.
    ImDrawList* dl = ImGui::GetWindowDrawList();
    auto marker = [&](float y, bool persistent, bool across, bool emphasis, float alpha) {
        const ImU32 dark  = IM_COL32(15, 15, 15, (int)(255 * alpha));
        const ImU32 light = IM_COL32(245, 245, 245, (int)(255 * alpha));
        const ImU32 core = persistent ? dark : light;
        const ImU32 rim  = persistent ? light : dark;
        const float x0 = bar.Min.x, x1 = bar.Max.x;
        if (across) {
            dl->AddLine(ImVec2(x0, y), ImVec2(x1, y), rim, emphasis ? 5.0f : 4.0f);
            dl->AddLine(ImVec2(x0, y), ImVec2(x1, y), core, 2.0f);
        }
        const float s = emphasis ? 5.0f : 4.0f;
        const float xh = x0 - 1.0f;
        dl->AddTriangleFilled(ImVec2(xh, y), ImVec2(xh - 1.5f * s, y - s), ImVec2(xh - 1.5f * s, y + s), rim);
        dl->AddTriangle(ImVec2(xh, y), ImVec2(xh - 1.5f * s, y - s), ImVec2(xh - 1.5f * s, y + s), core, 1.0f);
    };
    for (int i = 0; i < lines->count; ++i) {
        const float v = lines->values[i];
        const bool inside = lo <= v && v <= hi;
        const float y = inside ? y_of(v) : (v < lo ? (y_of(lo)) : y_of(hi));
        const bool emphasis = i == lines->drag || (lines->drag < 0 && i == nearest && hovered);
        marker(y, true, inside, emphasis, (lines->show ? 1.0f : 0.35f) * (inside ? 1.0f : 0.5f));
    }
    if (lines->hover) {
        marker(y_of(lines->hover_value), false, true, lines->drag >= 0 || nearest >= 0, 1.0f);
    }
}

bool color_scale_draw_controls(ColorScale* scale, const ColorScaleSpan& span, md_unit_t unit, const char* span_label, const char* span_tooltip) {
    ASSERT(scale);
    bool changed = false;

    ImGui::PushID("color_scale");

    const ImVec2 button_size = {ImGui::CalcItemWidth(), 0};
    if (ImPlot::ColormapButton(ImPlot::GetColormapName(scale->colormap), button_size, scale->colormap)) {
        ImGui::OpenPopup("Colormap");
    }
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
    changed |= ImGui::Checkbox("legend", &scale->show_legend);
    ImGui::SetItemTooltip("Show the colour map with its range and unit in the view.\nHold Alt to move or resize it.");
    if (ImGui::BeginPopup("Colormap")) {
        // The qualitative maps are left out: their colours do not order, so they cannot carry a value
        for (int m = COLOR_SCALE_FIRST_CONTINUOUS_COLORMAP; m < ImPlot::GetColormapCount(); ++m) {
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(m), button_size, m)) {
                scale->colormap = m;
                changed = true;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }

    // Shown in the display unit, kept in the values' own
    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    const float fscl = (float)scl;

    bool mode_changed = false;
    mode_changed |= ImGui::Checkbox("auto range", &scale->auto_range);
    if (span_tooltip && span_tooltip[0]) {
        ImGui::SetItemTooltip("Follow the span of the values: %s", span_tooltip);
    } else {
        ImGui::SetItemTooltip("Follow the span of the values");
    }
    ImGui::SameLine();
    mode_changed |= ImGui::Checkbox("symmetric", &scale->symmetric);
    ImGui::SetItemTooltip("A range centred on zero");
    if (mode_changed) {
        if (scale->symmetric && !scale->auto_range) {
            const float r = MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
            scale->range_beg = -r;
            scale->range_end =  r;
        }
        color_scale_update_range(scale, span);
        changed = true;
    }

    char unit_fmt[64];
    escape_unit_for_format(unit_fmt, sizeof(unit_fmt), unit_str);
    char fmt[96];
    if (scale->symmetric) {
        snprintf(fmt, sizeof(fmt), "%%.4g %s", unit_fmt);
    } else {
        snprintf(fmt, sizeof(fmt), "%%.4g to %%.4g %s", unit_fmt);
    }
    if (scale->symmetric) {
        // From zero to half as far again as the largest magnitude, and never short of the range itself
        float r = MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
        float limit = span.valid ? 1.5f * MAX(fabsf(span.min), fabsf(span.max)) : 0.0f;
        limit = MAX(limit, r);
        if (!(limit > 0.0f)) limit = 1.0f;
        float r_disp = r * fscl;
        if (ImGui::SliderFloat("range ±", &r_disp, 0.0f, limit * fscl, fmt)) {
            scale->range_beg = -r_disp / fscl;
            scale->range_end =  r_disp / fscl;
            scale->auto_range = false;
            changed = true;
        }
    } else {
        // The span padded by half its width on either side, and never short of the range itself
        float lo = scale->range_beg, hi = scale->range_end;
        if (span.valid) {
            float pad = 0.5f * (span.max - span.min);
            if (!(pad > 0.0f)) pad = MAX(0.5f * fabsf(span.max), 1.0f);
            lo = MIN(lo, span.min - pad);
            hi = MAX(hi, span.max + pad);
        }
        if (!(hi > lo)) { lo -= 1.0f; hi += 1.0f; }
        float beg = scale->range_beg * fscl;
        float end = scale->range_end * fscl;
        if (ImGui::RangeSliderFloat("range", &beg, &end, lo * fscl, hi * fscl, fmt)) {
            scale->range_beg = beg / fscl;
            scale->range_end = end / fscl;
            scale->auto_range = false;
            changed = true;
        }
    }

    if (span.valid) {
        ImGui::TextDisabled("%s: %.4g to %.4g %s", span_label ? span_label : "values", span.lo * scl, span.hi * scl, unit_str);
        if (span.lo != span.min || span.hi != span.max) {
            ImGui::SetItemTooltip("%s%sextremes %.4g to %.4g %s", span_tooltip ? span_tooltip : "", span_tooltip ? "; " : "", span.min * scl, span.max * scl, unit_str);
        } else if (span_tooltip && span_tooltip[0]) {
            ImGui::SetItemTooltip("%s", span_tooltip);
        }
    }

    ImGui::PopID();
    return changed;
}

void color_scale_draw_legend(const ColorScale& scale, md_unit_t unit, const char* label, const char* owner, int id, int slot, ColorScaleIsolines* lines) {
    if (lines) lines->hover = false;
    if (scale.range_end == scale.range_beg) return;

    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    char header[128];
    if (unit_str[0]) snprintf(header, sizeof(header), "%s (%s)", label ? label : "", unit_str);
    else             snprintf(header, sizeof(header), "%s", label ? label : "");

    const ImGuiStyle& style = ImGui::GetStyle();
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float width = MAX(ImGui::CalcTextSize(header).x + 2.0f * style.WindowPadding.x, 120.0f);
    const float height = 300.0f;
    // Where a legend starts out, the first time: down the right edge of the view in the order the
    // legends are drawn, a new column to the left of the last when the view is full
    const float top = 60.0f, gap = 10.0f;
    const int per_column = MAX(1, (int)((viewport->WorkSize.y - top) / (height + gap)));
    const int column = MAX(slot, 0) / per_column, row = MAX(slot, 0) % per_column;
    ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - (float)(column + 1) * (width + 20.0f),
                                   viewport->WorkPos.y + top + (float)row * (height + gap)), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(width, height), ImGuiCond_FirstUseEver);
    // Never narrower than its heading: the size is the user's and kept in the .ini, and what it heads
    // changes - another attribute, another unit - without the user having resized anything
    ImGui::SetNextWindowSizeConstraints(ImVec2(MAX(width, 80.0f), 120), ImVec2(MAX(width, 1000.0f), 2000));

    const bool editable = ImGui::IsKeyDown(ImGuiMod_Alt);
    // With isolines the scale takes the mouse where the cursor is over it (as it was drawn the frame before) and
    // while a line is dragged; everywhere else the mouse goes through to the view, as without them
    const bool interactive = lines && !editable && (lines->drag >= 0 || legend_hit_get(id).Contains(ImGui::GetIO().MousePos));
    ImGuiWindowFlags flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse |
                             ImGuiWindowFlags_NoFocusOnAppearing;
    if (!editable) {
        flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoNavFocus;
        flags |= interactive ? ImGuiWindowFlags_NoNav : ImGuiWindowFlags_NoInputs;
    }

    // The look of the text the script draws into the view: light on a translucent dark backing,
    // which reads on any background and any colour map
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.0f, 0.0f, 0.0f, editable ? 0.65f : 0.5f));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 5.0f);
    ImPlot::PushStyleColor(ImPlotCol_FrameBg, ImVec4(0, 0, 0, 0));

    char title[160];
    snprintf(title, sizeof(title), "Legend: %s###legend_%d", owner ? owner : "", id);
    if (ImGui::Begin(title, nullptr, flags)) {
        // Kept within the view: it grows to the right to fit a longer heading, and the view can
        // shrink under a position that was kept from a larger one
        const ImVec2 pos = ImGui::GetWindowPos(), size = ImGui::GetWindowSize();
        const ImVec2 lo = viewport->WorkPos;
        const ImVec2 hi = ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - size.x, viewport->WorkPos.y + viewport->WorkSize.y - size.y);
        const ImVec2 clamped = ImVec2(MAX(lo.x, MIN(pos.x, hi.x)), MAX(lo.y, MIN(pos.y, hi.y)));
        if (clamped.x != pos.x || clamped.y != pos.y) {
            ImGui::SetWindowPos(clamped);
        }

        ImGui::TextUnformatted(header);
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const ImVec2 frame_min = ImGui::GetCursorScreenPos();
        char scale_id[32];
        snprintf(scale_id, sizeof(scale_id), "##scale%d", id);
        ImPlot::ColormapScale(scale_id, scale.range_beg * scl, scale.range_end * scl, ImVec2(avail.x, avail.y), "%g", ImPlotColormapScaleFlags_NoLabel, scale.colormap);
        if (lines) {
            legend_isolines(lines, scale, scl, unit_str, id, frame_min, avail, interactive);
        }
    }
    ImGui::End();

    ImPlot::PopStyleColor();
    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor(2);
}
