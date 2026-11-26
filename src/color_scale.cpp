#include <color_scale.h>

#include <display_units.h>
#include <imgui_widgets.h>
#include <property_grid.h>
#include <app/IconsFontAwesome6.h>

#include <core/md_common.h>

#include <imgui.h>
#include <implot.h>
#include <implot_internal.h>   // SampleColormapU32, the ticker of ColormapScale

#include <math.h>
#include <stdio.h>

const char* color_scale_colormap_name(int colormap) {
    return (0 <= colormap && colormap < ImPlot::GetColormapCount()) ? ImPlot::GetColormapName(colormap) : nullptr;
}

int color_scale_colormap_index(md_str_t name) {
    char buf[64];
    if (md_str_empty(name) || name.len >= sizeof(buf)) return -1;
    md_str_copy_to_char_buf(buf, sizeof(buf), name);
    return ImPlot::GetColormapIndex(buf);
}

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
    MD_ASSERT(lines);
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
        desc.count = (uint32_t)MD_CLAMP(lines.count, 0, ISOLINES_MAX);
        for (uint32_t i = 0; i < desc.count; ++i) desc.values[i] = lines.values[i];
    }
    desc.hover       = lines.hover;
    desc.hover_value = lines.hover_value;
    desc.width       = MD_MAX(lines.width, 0.25f) * MD_MAX(px_per_point, 0.25f);
    return desc;
}

bool color_scale_draw_isoline_controls(ColorScaleIsolines* lines, const ColorScale& scale, md_unit_t unit) {
    MD_ASSERT(lines);
    bool changed = false;
    ImGui::PushID("isolines");
    const ImGuiStyle& style = ImGui::GetStyle();

    // Shown, and how wide, on one row: the width means nothing without the lines
    prop::row("isolines", "Lines on the surface where the field takes the given values.\n"
                          "On the legend: hover a value to see its line, click to add a line there,\n"
                          "drag a line to move it, right click to edit or remove one.");
    changed |= ImGui::Checkbox("##show", &lines->show);
    if (!lines->show) {
        ImGui::PopID();
        return changed;
    }
    ImGui::SameLine(0.0f, style.ItemInnerSpacing.x);
    ImGui::SetNextItemWidth(-FLT_MIN);
    changed |= ImGui::SliderFloat("##width", &lines->width, 0.5f, 6.0f, "width %.1f pt");

    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    const float fscl = (float)scl;

    // The values in the display unit, which is said once, beside the button that adds one
    const float ext   = fabsf(scale.range_end - scale.range_beg) * fscl;
    const float speed = ext > 0.0f ? ext / 300.0f : 0.001f;     // the whole range over about 300 points
    const float btn_w = ImGui::GetFrameHeight();
    const float gap   = style.ItemInnerSpacing.x;
    int remove = -1;
    for (int i = 0; i < lines->count; ++i) {
        ImGui::PushID(i);
        prop::row(i == 0 ? "values" : "");
        float v = lines->values[i] * fscl;
        ImGui::SetNextItemWidth(MD_MAX(prop::control_width() - btn_w - gap, 1.0f));
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
        ImGui::PopID();
    }
    if (remove >= 0) {
        color_scale_isolines_remove(lines, remove);
        changed = true;
    }

    prop::row(lines->count == 0 ? "values" : "");
    ImGui::BeginDisabled(lines->count >= ISOLINES_MAX);
    if (ImGui::Button(ICON_FA_PLUS, ImVec2(btn_w, 0))) {
        // Zero where the range holds it (the neutral contour of a potential), else the middle of the range
        const float lo = MD_MIN(scale.range_beg, scale.range_end), hi = MD_MAX(scale.range_beg, scale.range_end);
        const float v = (lo <= 0.0f && 0.0f <= hi) ? 0.0f : 0.5f * (lo + hi);
        changed |= color_scale_isolines_add(lines, v) >= 0;
    }
    ImGui::EndDisabled();
    if (lines->count >= ISOLINES_MAX) {
        ImGui::SetItemTooltip("%d lines at most", ISOLINES_MAX);
    } else {
        ImGui::SetItemTooltip("Add a line (or click the legend's colour bar where it should be)");
    }
    ImGui::SameLine(0.0f, 2.0f * gap);
    ImGui::AlignTextToFramePadding();
    if (lines->count == 0) {
        ImGui::TextDisabled("none yet");
    } else if (unit_str[0]) {
        ImGui::TextDisabled("%s", unit_str);
    }

    ImGui::PopID();
    return changed;
}

// Where a legend takes the mouse, by the legend's id, as it was last drawn: the frame of its scale (with isolines)
// and its heading (for its layout). Read before the legend's window is begun, to decide whether it takes input at all.
enum LegendHit { LegendHit_Scale, LegendHit_Heading };
static ImGuiStorage legend_hit_storage;

static ImGuiID legend_hit_key(int id, int k) {
    const int key[2] = {id, k};
    return ImHashData(key, sizeof(key), 0x15011E5u);
}

static ImRect legend_hit_get(int id, LegendHit hit) {
    const int k = 4 * (int)hit;
    return ImRect(legend_hit_storage.GetFloat(legend_hit_key(id, k + 0), 0.0f), legend_hit_storage.GetFloat(legend_hit_key(id, k + 1), 0.0f),
                  legend_hit_storage.GetFloat(legend_hit_key(id, k + 2), 0.0f), legend_hit_storage.GetFloat(legend_hit_key(id, k + 3), 0.0f));
}

static void legend_hit_set(int id, LegendHit hit, const ImRect& r) {
    const int k = 4 * (int)hit;
    legend_hit_storage.SetFloat(legend_hit_key(id, k + 0), r.Min.x);
    legend_hit_storage.SetFloat(legend_hit_key(id, k + 1), r.Min.y);
    legend_hit_storage.SetFloat(legend_hit_key(id, k + 2), r.Max.x);
    legend_hit_storage.SetFloat(legend_hit_key(id, k + 3), r.Max.y);
}

// The items of a legend's menu that lay it out; returns true when the layout changed
static bool legend_layout_items(bool* horizontal) {
    bool h = *horizontal;
    if (ImGui::MenuItem("Vertical", nullptr, !h))  h = false;
    if (ImGui::MenuItem("Horizontal", nullptr, h)) h = true;
    if (h == *horizontal) return false;
    *horizontal = h;
    return true;
}

// The colour bar of a scale, its values in the display unit (scl of them to one of its own), filling 'size' at the
// cursor: range_beg at the bottom of a vertical bar and the left of a horizontal one, also when the range runs
// backwards, with the ticks of ImPlot's axes. The vertical one is ImPlot::ColormapScale without its label; the
// horizontal one has its tick labels under the bar, which is inset from the ends by half the widest so that those
// at the ends are whole. Returns the bar, empty when there is no room.
static ImRect legend_bar(const ColorScale& scale, double scl, ImVec2 size, bool horizontal) {
    ImPlotContext& gp = *GImPlot;
    const ImPlotStyle& ps = gp.Style;
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    const ImRect frame(window->DC.CursorPos, ImVec2(window->DC.CursorPos.x + size.x, window->DC.CursorPos.y + size.y));
    ImGui::ItemSize(frame);
    if (!ImGui::ItemAdd(frame, ImGui::GetID("##legend_bar"))) return ImRect();

    const double a = scale.range_beg * scl, b = scale.range_end * scl;
    const ImPlotRange range(MD_MIN(a, b), MD_MAX(a, b));
    static char format[] = "%g";
    gp.CTicker.Reset();
    ImPlot::Locator_Default(gp.CTicker, range, horizontal ? size.x : size.y, !horizontal, ImPlot::Formatter_Default, format);

    ImRect bar;
    if (horizontal) {
        const float pad   = ps.LabelPadding.y + gp.CTicker.MaxSize.y;
        const float inset = MD_MAX(ps.PlotPadding.x, 0.5f * gp.CTicker.MaxSize.x);
        const float bar_h = MD_MAX(size.y - (pad + 2.0f * ps.PlotPadding.y), ps.MajorTickLen.x);
        bar = ImRect(ImVec2(frame.Min.x + inset, frame.Min.y + ps.PlotPadding.y), ImVec2(frame.Max.x - inset, frame.Min.y + ps.PlotPadding.y + bar_h));
    } else {
        const float pad   = ps.LabelPadding.x + gp.CTicker.MaxSize.x;
        const float bar_w = MD_MAX(size.x - (pad + 2.0f * ps.PlotPadding.x), ps.MajorTickLen.y);
        bar = ImRect(ImVec2(frame.Min.x + ps.PlotPadding.x, frame.Min.y + ps.PlotPadding.y), ImVec2(frame.Min.x + ps.PlotPadding.x + bar_w, frame.Max.y - ps.PlotPadding.y));
    }
    if (!(bar.GetWidth() > 1.0f && bar.GetHeight() > 1.0f) || a == b) return ImRect();

    // The labels that fit, evenly: every k-th along the bar counted from a major tick, k the least that leaves none
    // closer to the next than half a line. (The ticker thins them by their sum only, so labels could run into each
    // other along a horizontal bar.)
    {
        const float gap = 0.5f * ImGui::GetFontSize();
        const float len = horizontal ? bar.GetWidth() : bar.GetHeight();
        int order[64];
        int n = 0, anchor = 0;
        for (int i = 0; i < gp.CTicker.TickCount() && n < (int)MD_COUNTOF(order); ++i) {
            if (gp.CTicker.Ticks[i].ShowLabel) order[n++] = i;
        }
        auto at = [&](int i) { return (float)((gp.CTicker.Ticks[i].PlotPos - a) / (b - a)) * len; };
        auto half = [&](int i) { return 0.5f * (horizontal ? gp.CTicker.Ticks[i].LabelSize.x : gp.CTicker.Ticks[i].LabelSize.y); };
        for (int j = 1; j < n; ++j) {      // by position along the bar
            const int t = order[j];
            int k = j;
            for (; k > 0 && at(order[k - 1]) > at(t); --k) order[k] = order[k - 1];
            order[k] = t;
        }
        for (int j = 0; j < n; ++j) if (gp.CTicker.Ticks[order[j]].Major) { anchor = j; break; }
        int stride = 1;
        for (; stride < n; ++stride) {
            bool fits = true;
            int prev = -1;
            for (int j = anchor % stride; j < n && fits; j += stride) {
                if (prev >= 0) fits = at(order[prev]) + half(order[prev]) + gap <= at(order[j]) - half(order[j]);
                prev = j;
            }
            if (fits) break;
        }
        for (int j = 0; j < n; ++j) gp.CTicker.Ticks[order[j]].ShowLabel = (j - anchor) % stride == 0;
    }

    ImDrawList& dl = *window->DrawList;
    ImGui::RenderFrame(frame.Min, frame.Max, ImPlot::GetStyleColorU32(ImPlotCol_FrameBg), true, ImGui::GetStyle().FrameRounding);
    ImGui::PushClipRect(frame.Min, frame.Max, true);
    const ImU32 col_text = ImGui::GetColorU32(ImGuiCol_Text);
    const int cmap = scale.colormap;
    // Vertical: reversed, the first colour of the map at the bottom
    ImPlot::RenderColorBar(gp.ColormapData.GetKeys(cmap), gp.ColormapData.GetKeyCount(cmap), dl, bar, !horizontal, !horizontal, !gp.ColormapData.IsQual(cmap));
    for (int i = 0; i < gp.CTicker.TickCount(); ++i) {
        const ImPlotTick& tick = gp.CTicker.Ticks[i];
        const float t = (float)((tick.PlotPos - a) / (b - a));
        const ImU32 tick_col = ImPlot::CalcTextColor(gp.ColormapData.LerpTable(cmap, t));
        if (horizontal) {
            const float x = bar.Min.x + t * bar.GetWidth();
            const float len   = tick.Major ? ps.MajorTickLen.x : ps.MinorTickLen.x;
            const float thick = tick.Major ? ps.MajorTickSize.x : ps.MinorTickSize.x;
            if (x > bar.Min.x + 2.0f && x < bar.Max.x - 2.0f) {
                dl.AddLine(ImVec2(x, bar.Max.y - 1.0f), ImVec2(x, bar.Max.y - len), tick_col, thick);
            }
            if (tick.ShowLabel) dl.AddText(ImVec2(x - 0.5f * tick.LabelSize.x, bar.Max.y + ps.LabelPadding.y), col_text, gp.CTicker.GetText(i));
        } else {
            const float y = bar.Max.y - t * bar.GetHeight();
            const float len   = tick.Major ? ps.MajorTickLen.y : ps.MinorTickLen.y;
            const float thick = tick.Major ? ps.MajorTickSize.y : ps.MinorTickSize.y;
            if (y < bar.Max.y - 2.0f && y > bar.Min.y + 2.0f) {
                dl.AddLine(ImVec2(bar.Max.x - 1.0f, y), ImVec2(bar.Max.x - len, y), tick_col, thick);
            }
            if (tick.ShowLabel) dl.AddText(ImVec2(bar.Max.x + ps.LabelPadding.x, y - 0.5f * tick.LabelSize.y), col_text, gp.CTicker.GetText(i));
        }
    }
    dl.AddRect(bar.Min, bar.Max, ImPlot::GetStyleColorU32(ImPlotCol_PlotBorder));
    ImGui::PopClipRect();
    return bar;
}

// The line a popup of the legend edits (-1: none, the popup adds one at popup_value). One popup is open at a time.
static int   legend_popup_line  = -1;
static float legend_popup_value = 0.0f;

// The drag of a line on a legend (one at a time): how far the line was from the cursor when it was taken, so that
// one taken off-centre does not jump to the cursor, and whether it has moved since, so that a click on a line
// leaves it where it was
static float legend_drag_offset = 0.0f;
static float legend_drag_press = 0.0f;
static bool  legend_drag_moved = false;

// The isolines over the colour bar legend_bar just drew in 'frame': their markers, and with 'interactive' the
// hover, adding, dragging and editing. Positions are along the bar: y of a vertical one, x of a horizontal one. Its
// menu also lays the legend out; returns true when that changed 'horizontal'.
static bool legend_isolines(ColorScaleIsolines* lines, const ColorScale& scale, double scl, const char* unit_str,
                            ImRect frame, ImRect bar, bool* horizontal_layout, bool interactive) {
    const bool horizontal = *horizontal_layout;
    bool relaid = false;
    const float beg = scale.range_beg, end = scale.range_end;
    const float bar_h = horizontal ? bar.GetWidth() : bar.GetHeight();     // its length
    if (!(bar_h > 1.0f) || beg == end) return false;
    const float lo = MD_MIN(beg, end), hi = MD_MAX(beg, end);

    // range_beg at the bottom (left) of the bar and range_end at the top (right), as legend_bar draws them, also when
    // the range runs backwards
    auto y_of = [&](float v) {
        const float t = (v - beg) / (end - beg);
        return horizontal ? bar.Min.x + t * bar_h : bar.Max.y - t * bar_h;
    };
    auto v_of = [&](float y) {
        const float t = MD_CLAMP(horizontal ? (y - bar.Min.x) / bar_h : (bar.Max.y - y) / bar_h, 0.0f, 1.0f);
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
            v = MD_CLAMP(v, lo, hi);
        }
        return v;
    };

    char unit_fmt[64];
    escape_unit_for_format(unit_fmt, sizeof(unit_fmt), unit_str);
    char fmt[96];
    snprintf(fmt, sizeof(fmt), "%%.4g %s", unit_fmt);

    const ImVec2 mouse_pos = ImGui::GetIO().MousePos;
    const float mouse = horizontal ? mouse_pos.x : mouse_pos.y;     // along the bar
    bool hovered = false;
    int  nearest = -1;

    if (interactive) {
        ImGui::SetCursorScreenPos(frame.Min);
        ImGui::InvisibleButton("##isolines", frame.GetSize(), ImGuiButtonFlags_MouseButtonLeft);
        hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        const bool active = ImGui::IsItemActive();

        // The line whose marker is within reach of the cursor, the nearest if several are
        float best = 5.0f;
        for (int i = 0; i < lines->count; ++i) {
            const float v = lines->values[i];
            if (v < lo || v > hi) continue;
            const float d = fabsf(y_of(v) - mouse);
            if (d <= best) {
                best = d;
                nearest = i;
            }
        }

        // A press takes the line within reach, or adds one where it is, to be dragged into place
        if (ImGui::IsItemActivated()) {
            if (nearest >= 0) {
                lines->drag = nearest;
                legend_drag_offset = y_of(lines->values[nearest]) - mouse;
            } else {
                lines->drag = color_scale_isolines_add(lines, pick(mouse));
                legend_drag_offset = 0.0f;
            }
            legend_drag_press = mouse;
            legend_drag_moved = false;
            if (lines->drag >= 0) lines->show = true;
        }
        if (lines->drag >= 0) {
            if (active && lines->drag < lines->count) {
                legend_drag_moved |= mouse != legend_drag_press;
                if (legend_drag_moved) {
                    lines->values[lines->drag] = pick(mouse + legend_drag_offset);
                }
            } else {
                isolines_sort(lines, lines->drag);
                lines->drag = -1;
            }
        }

        // The cursor says what a press does: move the line within reach (or being moved), add one, or nothing when
        // there is no room for another
        if (lines->drag >= 0 || (hovered && nearest >= 0)) {
            ImGui::SetMouseCursor(horizontal ? ImGuiMouseCursor_ResizeEW : ImGuiMouseCursor_ResizeNS);
        } else if (hovered && lines->count >= ISOLINES_MAX) {
            ImGui::SetMouseCursor(ImGuiMouseCursor_NotAllowed);
        }

        if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            legend_popup_line  = nearest;
            legend_popup_value = pick(mouse);
            ImGui::OpenPopup("##isoline");
        }

        // What the cursor points at is drawn as the transient line, on the surface and here: the line being moved,
        // or the one within reach, or else the value under the cursor
        if (lines->drag >= 0) {
            lines->hover = true;
            lines->hover_value = lines->values[lines->drag];
        } else if (hovered) {
            lines->hover = true;
            lines->hover_value = nearest >= 0 ? lines->values[nearest] : pick(mouse);
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
        ImGui::Separator();
        relaid = legend_layout_items(horizontal_layout);
        ImGui::EndPopup();
    }

    // The markers: across the bar and as a handle to the left of (above) it, drawn as the lines are on the surface - the
    // persistent ones dark with a light rim, the transient one the other way around. A line beyond the range keeps
    // its handle at the end it is beyond; hidden lines are faint.
    ImDrawList* dl = ImGui::GetWindowDrawList();
    auto marker = [&](float y, bool persistent, bool across, bool emphasis, float alpha) {
        const ImU32 dark  = IM_COL32(15, 15, 15, (int)(255 * alpha));
        const ImU32 light = IM_COL32(245, 245, 245, (int)(255 * alpha));
        const ImU32 core = persistent ? dark : light;
        const ImU32 rim  = persistent ? light : dark;
        // In the bar's own terms: y along it, x across it from the side of the handle (flipped into the view's
        // axes for a horizontal bar)
        auto at = [&](float along, float across_) { return horizontal ? ImVec2(along, across_) : ImVec2(across_, along); };
        const float x0 = horizontal ? bar.Min.y : bar.Min.x, x1 = horizontal ? bar.Max.y : bar.Max.x;
        if (across) {
            dl->AddLine(at(y, x0), at(y, x1), rim, emphasis ? 5.0f : 4.0f);
            dl->AddLine(at(y, x0), at(y, x1), core, 2.0f);
        }
        const float s = emphasis ? 5.0f : 4.0f;
        const float xh = x0 - 1.0f;
        dl->AddTriangleFilled(at(y, xh), at(y - s, xh - 1.5f * s), at(y + s, xh - 1.5f * s), rim);
        dl->AddTriangle(at(y, xh), at(y - s, xh - 1.5f * s), at(y + s, xh - 1.5f * s), core, 1.0f);
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
    return relaid;
}

bool color_scale_draw_controls(ColorScale* scale, const ColorScaleSpan& span, md_unit_t unit, const char* span_label, const char* span_tooltip) {
    MD_ASSERT(scale);
    bool changed = false;

    ImGui::PushID("color_scale");

    prop::row("colormap");
    const ImVec2 button_size = {ImGui::CalcItemWidth(), 0};
    if (ImPlot::ColormapButton(ImPlot::GetColormapName(scale->colormap), button_size, scale->colormap)) {
        ImGui::OpenPopup("Colormap");
    }
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

    prop::row("legend", "The colour map with its range and unit, in the view. Right click its heading to lay it out\n"
                        "vertically or horizontally; hold Alt to move or resize it.");
    changed |= ImGui::Checkbox("##legend", &scale->show_legend);

    // Shown in the display unit, kept in the values' own
    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    const float fscl = (float)scl;

    char unit_fmt[64];
    escape_unit_for_format(unit_fmt, sizeof(unit_fmt), unit_str);
    char fmt[96];
    if (scale->symmetric) {
        snprintf(fmt, sizeof(fmt), (const char*)u8"± %%.4g %s", unit_fmt);
    } else {
        snprintf(fmt, sizeof(fmt), "%%.4g to %%.4g %s", unit_fmt);
    }

    // The range, and whether it follows the values: on the same row, as setting the range by hand stops it following
    bool mode_changed = false;
    prop::row("range");
    mode_changed |= ImGui::Checkbox("auto", &scale->auto_range);
    if (span_tooltip && span_tooltip[0]) {
        ImGui::SetItemTooltip("Follow the span of the values: %s", span_tooltip);
    } else {
        ImGui::SetItemTooltip("Follow the span of the values");
    }
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (scale->symmetric) {
        // From zero to half as far again as the largest magnitude, and never short of the range itself
        float r = MD_MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
        float limit = span.valid ? 1.5f * MD_MAX(fabsf(span.min), fabsf(span.max)) : 0.0f;
        limit = MD_MAX(limit, r);
        if (!(limit > 0.0f)) limit = 1.0f;
        float r_disp = r * fscl;
        if (ImGui::SliderFloat("##range", &r_disp, 0.0f, limit * fscl, fmt)) {
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
            if (!(pad > 0.0f)) pad = MD_MAX(0.5f * fabsf(span.max), 1.0f);
            lo = MD_MIN(lo, span.min - pad);
            hi = MD_MAX(hi, span.max + pad);
        }
        if (!(hi > lo)) { lo -= 1.0f; hi += 1.0f; }
        float beg = scale->range_beg * fscl;
        float end = scale->range_end * fscl;
        if (ImGui::RangeSliderFloat("##range", &beg, &end, lo * fscl, hi * fscl, fmt)) {
            scale->range_beg = beg / fscl;
            scale->range_end = end / fscl;
            scale->auto_range = false;
            changed = true;
        }
    }

    prop::row("symmetric", "A range centred on zero");
    mode_changed |= ImGui::Checkbox("##symmetric", &scale->symmetric);
    if (mode_changed) {
        if (scale->symmetric && !scale->auto_range) {
            const float r = MD_MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
            scale->range_beg = -r;
            scale->range_end =  r;
        }
        color_scale_update_range(scale, span);
        changed = true;
    }

    if (span.valid) {
        prop::hint("%s: %.4g to %.4g %s", span_label ? span_label : "values", span.lo * scl, span.hi * scl, unit_str);
        if (span.lo != span.min || span.hi != span.max) {
            ImGui::SetItemTooltip("%s%sextremes %.4g to %.4g %s", span_tooltip ? span_tooltip : "", span_tooltip ? "; " : "", span.min * scl, span.max * scl, unit_str);
        } else if (span_tooltip && span_tooltip[0]) {
            ImGui::SetItemTooltip("%s", span_tooltip);
        }
    }

    ImGui::PopID();
    return changed;
}

// The legend whose layout changed, to start the other layout's window where it was (one at a time)
static int    legend_relaid_id = -1;
static ImVec2 legend_relaid_pos = {};

void color_scale_draw_legend(ColorScale* scale_ptr, md_unit_t unit, const char* label, const char* owner, int id, int slot, ColorScaleIsolines* lines) {
    MD_ASSERT(scale_ptr);
    const ColorScale& scale = *scale_ptr;
    if (lines) lines->hover = false;
    if (scale.range_end == scale.range_beg) return;
    const bool horizontal = scale.legend_horizontal;

    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    char header[128];
    if (unit_str[0]) snprintf(header, sizeof(header), "%s (%s)", label ? label : "", unit_str);
    else             snprintf(header, sizeof(header), "%s", label ? label : "");

    const ImGuiStyle& style = ImGui::GetStyle();
    const ImPlotStyle& ps = ImPlot::GetStyle();
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float header_w = ImGui::CalcTextSize(header).x + 2.0f * style.WindowPadding.x;
    // What it starts out as: a vertical bar 300 high, or a horizontal one 18 high and 300 wide, its values under it.
    // The least: its heading, and a bar of some length (or, horizontal, some thickness)
    const float h_height = 2.0f * style.WindowPadding.y + ImGui::GetTextLineHeightWithSpacing() + 2.0f * ps.PlotPadding.y + ps.LabelPadding.y + ImGui::GetTextLineHeight();
    const float width  = horizontal ? MD_MAX(header_w, 300.0f) : MD_MAX(header_w, 120.0f);
    const float height = horizontal ? h_height + 18.0f : 300.0f;
    const ImVec2 min_size = horizontal ? ImVec2(MD_MAX(header_w, 160.0f), h_height + 6.0f) : ImVec2(MD_MAX(header_w, 80.0f), 120.0f);
    const ImVec2 max_size = horizontal ? ImVec2(MD_MAX(header_w, 2000.0f), 600.0f) : ImVec2(MD_MAX(header_w, 1000.0f), 2000.0f);

    // Where a legend starts out, the first time: down the right edge of the view in the order the legends are drawn,
    // a new column to the left of the last when the view is full, in places as large as a vertical one's whichever
    // way each is laid out, so that none starts out over another. One laid out anew starts where it was.
    if (legend_relaid_id == id) {
        ImGui::SetNextWindowPos(legend_relaid_pos, ImGuiCond_Always);
        legend_relaid_id = -1;
    } else {
        const float top = 60.0f, gap = 10.0f, place_h = 300.0f;
        const int per_column = MD_MAX(1, (int)((viewport->WorkSize.y - top) / (place_h + gap)));
        const int column = MD_MAX(slot, 0) / per_column, row = MD_MAX(slot, 0) % per_column;
        ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - (float)(column + 1) * (width + 20.0f),
                                       viewport->WorkPos.y + top + (float)row * (place_h + gap)), ImGuiCond_FirstUseEver);
    }
    ImGui::SetNextWindowSize(ImVec2(width, height), ImGuiCond_FirstUseEver);
    // Never narrower than its heading: the size is the user's and kept in the .ini, and what it heads
    // changes - another attribute, another unit - without the user having resized anything
    ImGui::SetNextWindowSizeConstraints(min_size, max_size);

    const bool editable = ImGui::IsKeyDown(ImGuiMod_Alt);
    // With isolines the scale takes the mouse where the cursor is over it (as it was drawn the frame before) and
    // while a line is dragged, and the heading takes it for the menu of the layout; everywhere else the mouse goes
    // through to the view
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    const bool lines_interactive = lines && !editable && (lines->drag >= 0 || legend_hit_get(id, LegendHit_Scale).Contains(mouse));
    const bool interactive = !editable && (lines_interactive || legend_hit_get(id, LegendHit_Heading).Contains(mouse));
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

    // A window of each layout, each with its own place and size in the .ini
    char title[160];
    snprintf(title, sizeof(title), "Legend: %s###legend_%d%s", owner ? owner : "", id, horizontal ? "_h" : "");
    bool relaid = false;
    if (ImGui::Begin(title, nullptr, flags)) {
        // Kept within the view: it grows to the right to fit a longer heading, and the view can
        // shrink under a position that was kept from a larger one
        const ImVec2 pos = ImGui::GetWindowPos(), size = ImGui::GetWindowSize();
        const ImVec2 lo = viewport->WorkPos;
        const ImVec2 hi = ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - size.x, viewport->WorkPos.y + viewport->WorkSize.y - size.y);
        const ImVec2 clamped = ImVec2(MD_MAX(lo.x, MD_MIN(pos.x, hi.x)), MD_MAX(lo.y, MD_MIN(pos.y, hi.y)));
        if (clamped.x != pos.x || clamped.y != pos.y) {
            ImGui::SetWindowPos(clamped);
        }

        ImGui::TextUnformatted(header);
        // The heading: the top of the window, across it. A right click there (anywhere while Alt is held) lays it out.
        const ImRect heading(pos, ImVec2(pos.x + size.x, ImGui::GetItemRectMax().y + 0.5f * style.ItemSpacing.y));
        legend_hit_set(id, LegendHit_Heading, heading);
        const bool hovered = ImGui::IsWindowHovered();
        if (hovered && (editable || heading.Contains(mouse)) && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            ImGui::OpenPopup("##legend_layout");
        }

        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const ImVec2 frame_min = ImGui::GetCursorScreenPos();
        const ImRect frame(frame_min, ImVec2(frame_min.x + avail.x, frame_min.y + avail.y));
        const ImRect bar = legend_bar(scale, scl, avail, horizontal);
        if (lines) {
            legend_hit_set(id, LegendHit_Scale, frame);
            if (bar.GetWidth() > 0.0f) {
                relaid |= legend_isolines(lines, scale, scl, unit_str, frame, bar, &scale_ptr->legend_horizontal, lines_interactive);
            }
        }

        if (ImGui::BeginPopup("##legend_layout")) {
            relaid |= legend_layout_items(&scale_ptr->legend_horizontal);
            ImGui::EndPopup();
        }
        if (relaid) {
            legend_relaid_id = id;
            legend_relaid_pos = pos;
        }
    }
    ImGui::End();

    ImPlot::PopStyleColor();
    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor(2);
}
