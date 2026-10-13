#define IMGUI_DEFINE_MATH_OPERATORS
#include <timeline_ranges_ui.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>
#include <imgui_widgets.h>
#include <app/IconsFontAwesome6.h>

#include <core/md_common.h>
#include <core/md_str.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

namespace {

enum DragKind {
    Drag_None = 0,
    Drag_Create,    // A new range, from where it was pressed to the pointer
    Drag_Move,
    Drag_Resize,
    Drag_Pan,       // Empty lane space, which pans the timeline
    Drag_Scrub,     // Ctrl-drag, which moves the current frame
    Drag_Out,       // A range pulled out of the lane, taken elsewhere as a drag and drop payload
};

struct Drag {
    DragKind kind = Drag_None;
    uint32_t id = 0;
    bool     hold_beg = false;          // Resize: the end held is the one stored in beg, wherever it has been dragged to
    int32_t  orig_beg = 0;
    int32_t  orig_end = 0;
    int      anchor = 0;                // The frame under the pointer when it was pressed
    int      current = 0;               // Create: the frame under it now
    float    press_x = 0;
    double   view_beg = 0;
    double   view_end = 0;
};

// One timeline, so one of these. Nothing here is saved: it is what the pointer is doing.
struct State {
    Drag drag;
    TimelineRangeHit hover;
    uint32_t hover_id = 0;

    uint32_t edit_id = 0;               // The range whose popup is open
    uint32_t rename_id = 0;             // The range being named
    bool     rename_request = false;    // Open the name field when the lane is drawn next
    bool     rename_focus = false;
    char     rename_prev[TIMELINE_RANGE_NAME_CAP] = "";

    // Where the plots' area was: the lane is drawn above the plots and aligned with them
    float plot_x0 = 0;
    float plot_x1 = 0;
    int   plot_frame = -1;

    // Where the lane was last drawn: a range dragged far enough from it is pulled out
    float lane_y0 = 0;
    float lane_y1 = 0;

    // Lit up for another view (timeline_ranges_highlight), as of a frame
    uint32_t highlight_id = 0;
    int      highlight_frame = -1;
};

State ui;

ImU32 range_color(const TimelineRange& r, float alpha) {
    return ImGui::ColorConvertFloat4ToU32(ImVec4(r.color.x, r.color.y, r.color.z, r.color.w * alpha));
}

bool is_held(uint32_t id) {
    return id != 0 && (ui.drag.kind == Drag_Move || ui.drag.kind == Drag_Resize) && ui.drag.id == id;
}

// Under the pointer, held, being edited, or lit up for another view: lit up in the lane and shaded in the plots
bool is_hot(uint32_t id) {
    const bool lit = ui.highlight_id == id && ImGui::GetFrameCount() - ui.highlight_frame <= 1;
    return id != 0 && (ui.hover_id == id || ui.edit_id == id || ui.rename_id == id || is_held(id) || lit);
}

double frame_time(const TimelineRangeAxis& a, double frame) {
    const int last = a.num_frames - 1;
    const double f = CLAMP(frame, 0.0, (double)last);
    const int f0 = (int)f;
    const int f1 = MIN(f0 + 1, last);
    const double t = f - f0;
    return a.frame_times[f0] * (1.0 - t) + a.frame_times[f1] * t;
}

int nearest_frame(const TimelineRangeAxis& a, double time) {
    const float* ft = a.frame_times;
    int lo = 0;
    int hi = a.num_frames - 1;
    if (hi <= 0 || time <= ft[0]) return 0;
    if (time >= ft[hi]) return hi;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (ft[mid] <= time) lo = mid;
        else hi = mid;
    }
    return (time - ft[lo] <= ft[hi] - time) ? lo : hi;
}

float plot_width() { return MAX(ui.plot_x1 - ui.plot_x0, 1.0f); }

double x_to_time(const TimelineRangeAxis& a, float x) {
    return *a.view_beg + (double)(x - ui.plot_x0) / plot_width() * (*a.view_end - *a.view_beg);
}

float time_to_x(const TimelineRangeAxis& a, double t) {
    const double span = *a.view_end - *a.view_beg;
    return span > 0 ? ui.plot_x0 + (float)((t - *a.view_beg) / span) * plot_width() : ui.plot_x0;
}

float frame_x(const TimelineRangeAxis& a, double frame) { return time_to_x(a, frame_time(a, frame)); }

void add_dashed_rect(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 col) {
    const float dash = 4.0f;
    const float gap  = 3.0f;
    auto line = [&](ImVec2 p, ImVec2 q) {
        const ImVec2 d = q - p;
        const float len = sqrtf(d.x * d.x + d.y * d.y);
        if (len <= 0.0f) return;
        const ImVec2 u = d / len;
        for (float t = 0.0f; t < len; t += dash + gap) {
            dl->AddLine(p + u * t, p + u * MIN(t + dash, len), col);
        }
    };
    line(a, ImVec2(b.x, a.y));
    line(ImVec2(b.x, a.y), b);
    line(b, ImVec2(a.x, b.y));
    line(ImVec2(a.x, b.y), a);
}

// The end of a range which is lit: the one under the pointer, or the one held, wherever it has been dragged to
TimelineRangePart lit_part(const TimelineRange& r) {
    if (ui.drag.kind == Drag_Resize && ui.drag.id == r.id) {
        const int32_t held  = ui.drag.hold_beg ? r.beg : r.end;
        const int32_t other = ui.drag.hold_beg ? r.end : r.beg;
        return held <= other ? TimelineRangePart_Beg : TimelineRangePart_End;
    }
    if (ui.hover_id == r.id) return ui.hover.part;
    return TimelineRangePart_None;
}

void label_of(char* buf, size_t cap, const TimelineRanges* rs, const TimelineRange& r) {
    snprintf(buf, cap, "%s%s", rs->filter_id == r.id ? ICON_FA_FILTER " " : "", r.name);
}

void tooltip_frames(const TimelineRangeAxis& a, bool marker, int beg, int end) {
    if (marker) {
        if (a.time_unit[0]) ImGui::Text("Frame %d (%.4g %s)", beg, a.frame_times[beg], a.time_unit);
        else                ImGui::Text("Frame %d", beg);
    } else {
        if (a.time_unit[0]) ImGui::Text("Frames %d to %d, %d frames (%.4g to %.4g %s)", beg, end, end - beg + 1, a.frame_times[beg], a.frame_times[end], a.time_unit);
        else                ImGui::Text("Frames %d to %d, %d frames", beg, end, end - beg + 1);
    }
}

void range_tooltip(const TimelineRanges* rs, const TimelineRange& r, const TimelineRangeAxis& a, double current_frame) {
    int beg, end;
    if (!timeline_range_frames(&beg, &end, r, current_frame, a.num_frames)) return;
    ImGui::BeginTooltip();
    char label[TIMELINE_RANGE_NAME_CAP + 8];
    label_of(label, sizeof(label), rs, r);
    ImGui::TextUnformatted(label);
    tooltip_frames(a, timeline_range_is_marker(r), beg, end);
    if (r.follow) {
        ImGui::TextDisabled("Follows the current frame, %+d to %+d", MIN(r.beg, r.end), MAX(r.beg, r.end));
    }
    if (rs->filter_id == r.id) {
        ImGui::TextDisabled("The timeline filter follows it");
    }
    if (r.note[0]) {
        ImGui::Separator();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(r.note);
        ImGui::PopTextWrapPos();
    }
    ImGui::EndTooltip();
}

void request_rename(const TimelineRange& r) {
    ui.rename_id = r.id;
    ui.rename_request = true;
    str_copy_to_char_buf(ui.rename_prev, sizeof(ui.rename_prev), str_from_cstr(r.name));
}

void end_drag(TimelineRanges* rs, const TimelineRangeAxis& axis, bool commit) {
    const Drag d = ui.drag;
    ui.drag = {};
    switch (d.kind) {
    case Drag_Create:
        if (commit) {
            // A click marks the frame, a drag the frames it went over
            const ImGuiIO& io = ImGui::GetIO();
            const bool click = fabsf(io.MousePos.x - d.press_x) < io.MouseDragThreshold;
            const int beg = click ? d.anchor : MIN(d.anchor, d.current);
            const int end = click ? d.anchor : MAX(d.anchor, d.current);
            if (TimelineRange* r = timeline_ranges_add(rs, nullptr, beg, end, false)) {
                request_rename(*r);
            }
        }
        break;
    case Drag_Move:
    case Drag_Resize:
        if (TimelineRange* r = timeline_ranges_find(rs, d.id)) {
            if (!commit) {
                r->beg = d.orig_beg;
                r->end = d.orig_end;
            }
            timeline_range_normalize(r);
        }
        break;
    case Drag_Pan:
        if (!commit) {
            *axis.view_beg = d.view_beg;
            *axis.view_end = d.view_end;
        }
        break;
    default:
        break;
    }
}

// Follows the pointer with whatever is held, wherever the pointer is: it may leave the lane while it drags
void update_drag(TimelineRanges* rs, const TimelineRangeAxis& axis, double* current_frame) {
    Drag& d = ui.drag;
    if (d.kind == Drag_None) return;
    if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
        end_drag(rs, axis, false);
        return;
    }

    const float mx = ImGui::GetIO().MousePos.x;
    const int   f  = nearest_frame(axis, x_to_time(axis, mx));
    const int   last = axis.num_frames - 1;

    switch (d.kind) {
    case Drag_Create:
        d.current = f;
        break;
    case Drag_Move:
        if (TimelineRange* r = timeline_ranges_find(rs, d.id)) {
            // Pulled well out of the lane, up or down: taken elsewhere, and put back where it was
            const float my = ImGui::GetIO().MousePos.y;
            const float reach = ImGui::GetFontSize() * 4.0f;
            if (ui.lane_y1 > ui.lane_y0 && (my < ui.lane_y0 - reach || my > ui.lane_y1 + reach)) {
                r->beg = d.orig_beg;
                r->end = d.orig_end;
                d.kind = Drag_Out;
                break;
            }
            int delta = f - d.anchor;
            if (!r->follow) {
                // A fixed range stops at the ends of the trajectory, as a whole
                const int lo = MIN(d.orig_beg, d.orig_end);
                const int hi = MAX(d.orig_beg, d.orig_end);
                delta = CLAMP(delta, MIN(-lo, 0), MAX(last - hi, 0));
            }
            r->beg = d.orig_beg + delta;
            r->end = d.orig_end + delta;
        }
        break;
    case Drag_Resize:
        if (TimelineRange* r = timeline_ranges_find(rs, d.id)) {
            const int32_t v = r->follow ? f - (int32_t)llround(*current_frame) : f;
            if (d.hold_beg) r->beg = v;
            else            r->end = v;
        }
        break;
    case Drag_Pan: {
        const double dt = (double)(mx - d.press_x) / plot_width() * (d.view_end - d.view_beg);
        *axis.view_beg = d.view_beg - dt;
        *axis.view_end = d.view_end - dt;
        break;
    }
    case Drag_Scrub:
        *current_frame = f;
        break;
    default:
        break;
    }

    if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        end_drag(rs, axis, true);
    }
}

void start_create(int frame) {
    ui.drag = {};
    ui.drag.kind    = Drag_Create;
    ui.drag.anchor  = frame;
    ui.drag.current = frame;
    ui.drag.press_x = ImGui::GetIO().MousePos.x;
}

}  // namespace

bool timeline_ranges_dragging() {
    return ui.drag.kind != Drag_None;
}

uint32_t timeline_ranges_hovered() {
    if (ui.drag.kind == Drag_Move || ui.drag.kind == Drag_Resize || ui.drag.kind == Drag_Out) return ui.drag.id;
    return ui.drag.kind == Drag_None ? ui.hover_id : 0;
}

void timeline_ranges_highlight(uint32_t id) {
    ui.highlight_id = id;
    ui.highlight_frame = ImGui::GetFrameCount();
}

void timeline_ranges_lane(TimelineRanges* rs, const TimelineRangeAxis& axis, double* current_frame) {
    if (axis.num_frames <= 0 || !axis.frame_times) return;

    // First, so that the bars are drawn where the pointer has taken them
    update_drag(rs, axis, current_frame);

    const ImGuiIO& io = ImGui::GetIO();
    const double cur  = *current_frame;
    const int    nf   = axis.num_frames;
    const float  font = ImGui::GetFontSize();

    const bool  one_lane = rs->layout != TimelineRangeLayout_LanePerRange;
    const float row_h    = ImGui::GetTextLineHeight() + 4.0f;
    const float row_gap  = 2.0f;
    const int   rows     = one_lane ? 1 : MAX(rs->count, 1);
    const float height   = one_lane ? ImGui::GetFrameHeight() : rows * row_h + (rows - 1) * row_gap;

    const ImVec2 pos   = ImGui::GetCursorScreenPos();
    const float  width = MAX(ImGui::GetContentRegionAvail().x, 1.0f);
    if (ui.plot_x1 <= ui.plot_x0) {
        // Before the plots have been drawn once
        ui.plot_x0 = pos.x;
        ui.plot_x1 = pos.x + width;
    }

    ImGui::InvisibleButton("##timeline_ranges", ImVec2(width, height), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
    const bool hovered = ImGui::IsItemHovered();
    ImGui::SetItemKeyOwner(ImGuiKey_MouseWheelY);
    ui.lane_y0 = pos.y;
    ui.lane_y1 = pos.y + height;

    // A range pulled out of the lane: carried along as long as the button is held
    if (ui.drag.kind == Drag_Out) {
        const TimelineRange* r = timeline_ranges_find(rs, ui.drag.id);
        if (r && ImGui::BeginDragDropSource(ImGuiDragDropFlags_SourceNoDisableHover)) {
            const uint32_t id = r->id;
            ImGui::SetDragDropPayload(TIMELINE_RANGE_DND, &id, sizeof(id));
            ImGui::ColorButton("##color", ImVec4(r->color.x, r->color.y, r->color.z, 1.0f), ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(ImGui::GetTextLineHeight(), ImGui::GetTextLineHeight()));
            ImGui::SameLine();
            ImGui::TextUnformatted(r->name);
            ImGui::TextDisabled("Drop on a distribution to show it over these frames too");
            ImGui::EndDragDropSource();
        }
    }

    const float x0 = ui.plot_x0;
    const float x1 = ui.plot_x1;
    const float y0 = pos.y;
    const float y1 = pos.y + height;

    // Where each range is
    int order[TIMELINE_MAX_RANGES];
    int depth[TIMELINE_MAX_RANGES];
    int top_first[TIMELINE_MAX_RANGES];
    TimelineRangeBar bars[TIMELINE_MAX_RANGES];
    const int n = timeline_ranges_draw_order(order, rs, cur, nf);
    timeline_ranges_nesting(depth, order, n, rs, cur, nf);
    for (int k = 0; k < n; ++k) {
        const int i = order[k];
        const TimelineRange& r = rs->items[i];
        int beg, end;
        timeline_range_frames(&beg, &end, r, cur, nf);
        TimelineRangeBar& b = bars[i];
        b.visible = true;
        b.marker  = timeline_range_is_marker(r);
        b.x0 = frame_x(axis, beg);
        b.x1 = frame_x(axis, end);
        if (!b.marker && b.x1 - b.x0 < 4.0f) {
            const float mid = (b.x0 + b.x1) * 0.5f;
            b.x0 = mid - 2.0f;
            b.x1 = mid + 2.0f;
        }
        if (one_lane) {
            // Nested ranges are drawn a little thinner than the ones they lie on, so that both can be seen
            const float inset = b.marker ? 0.0f : MIN(depth[k], 2) * 3.0f;
            b.y0 = y0 + inset;
            b.y1 = y1 - inset;
        } else {
            b.y0 = y0 + i * (row_h + row_gap);
            b.y1 = b.y0 + row_h;
        }
        top_first[n - 1 - k] = i;
    }

    // What the pointer is over, and what it does there
    ui.hover = {};
    ui.hover_id = 0;
    uint32_t remove_id = 0;
    if (hovered && ui.drag.kind == Drag_None) {
        ui.hover = timeline_ranges_hit_test(bars, top_first, n, io.MousePos.x, io.MousePos.y, font * 0.5f, font * 0.45f);
        const TimelineRange* hit = ui.hover.index != -1 ? &rs->items[ui.hover.index] : nullptr;
        ui.hover_id = hit ? hit->id : 0;

        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const int f = nearest_frame(axis, x_to_time(axis, io.MousePos.x));
            if (io.KeyShift) {
                start_create(f);
            } else if (io.KeyCtrl) {
                ui.drag.kind = Drag_Scrub;
                *current_frame = f;
            } else if (hit) {
                Drag d = {};
                d.id       = hit->id;
                d.orig_beg = hit->beg;
                d.orig_end = hit->end;
                d.anchor   = f;
                d.press_x  = io.MousePos.x;
                if (ui.hover.part == TimelineRangePart_Body) {
                    d.kind = Drag_Move;
                } else {
                    // The handle at the lower end holds whichever of beg and end is the lower one
                    d.kind = Drag_Resize;
                    const bool beg_is_lower = hit->beg <= hit->end;
                    d.hold_beg = (ui.hover.part == TimelineRangePart_Beg) == beg_is_lower;
                }
                ui.drag = d;
            } else {
                ui.drag = {};
                ui.drag.kind     = Drag_Pan;
                ui.drag.press_x  = io.MousePos.x;
                ui.drag.view_beg = *axis.view_beg;
                ui.drag.view_end = *axis.view_end;
            }
        }
        if (hit && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            request_rename(*hit);
        }
        if (hit && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            ui.edit_id = hit->id;
            ImGui::OpenPopup("##timeline_range_edit");
        }
        if (hit && !io.WantTextInput && (ImGui::IsKeyPressed(ImGuiKey_Delete, false) || ImGui::IsKeyPressed(ImGuiKey_Backspace, false))) {
            remove_id = hit->id;
        }
        if (io.MouseWheel != 0.0f) {
            const double t = x_to_time(axis, io.MousePos.x);
            const double scale = pow(0.9, (double)io.MouseWheel);
            *axis.view_beg = t - (t - *axis.view_beg) * scale;
            *axis.view_end = t + (*axis.view_end - t) * scale;
        }
    }

    if (hovered || ui.drag.kind != Drag_None) {
        const DragKind k = ui.drag.kind;
        const TimelineRangePart part = ui.hover.part;
        if (k == Drag_Resize || k == Drag_Create || (k == Drag_None && (part == TimelineRangePart_Beg || part == TimelineRangePart_End))) {
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        } else if (k == Drag_Move || (k == Drag_None && part == TimelineRangePart_Body)) {
            ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
        }
    }

    // Drawing: the track, the ranges bottom first, the range being created, the current frame
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(ImVec2(x0, y0 - 1.0f), ImVec2(x1, y1 + 1.0f), true);

    const ImU32 track_col = ImGui::GetColorU32(ImGuiCol_FrameBg);
    const ImU32 text_col  = ImGui::GetColorU32(ImGuiCol_Text);
    if (one_lane) {
        dl->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, y1), track_col, 3.0f);
    } else {
        for (int r = 0; r < rows; ++r) {
            const float ry = y0 + r * (row_h + row_gap);
            dl->AddRectFilled(ImVec2(x0, ry), ImVec2(x1, ry + row_h), track_col, 3.0f);
        }
    }

    for (int k = 0; k < n; ++k) {
        const int i = order[k];
        const TimelineRange& r = rs->items[i];
        const TimelineRangeBar& b = bars[i];
        const bool hot = is_hot(r.id);
        char label[TIMELINE_RANGE_NAME_CAP + 8];
        label_of(label, sizeof(label), rs, r);
        const float cy = (b.y0 + b.y1) * 0.5f;

        if (b.marker) {
            const float h = (b.y1 - b.y0) * 0.42f;
            const ImVec2 p0(b.x0, cy - h), p1(b.x0 + h, cy), p2(b.x0, cy + h), p3(b.x0 - h, cy);
            dl->AddQuadFilled(p0, p1, p2, p3, range_color(r, hot ? 1.0f : 0.85f));
            if (hot) dl->AddQuad(p0, p1, p2, p3, text_col, 1.5f);

            // The name to the right, as far as the next marker in the same lane
            float limit = x1;
            for (int j = 0; j < rs->count; ++j) {
                const TimelineRangeBar& o = bars[j];
                if (j == i || !o.visible || !o.marker || o.y0 != b.y0 || o.x0 <= b.x0) continue;
                limit = MIN(limit, o.x0 - h - 2.0f);
            }
            const float tx = b.x0 + h + 3.0f;
            if (limit - tx > font) {
                dl->PushClipRect(ImVec2(tx, b.y0), ImVec2(limit, b.y1), true);
                dl->AddText(ImVec2(tx, cy - font * 0.5f), text_col, label);
                dl->PopClipRect();
            }
        } else {
            dl->AddRectFilled(ImVec2(b.x0, b.y0), ImVec2(b.x1, b.y1), range_color(r, hot ? 0.5f : 0.28f));
            const ImU32 edge = range_color(r, hot ? 1.0f : 0.75f);
            if (r.follow) {
                // Dashed: it moves with the current frame
                add_dashed_rect(dl, ImVec2(b.x0 + 0.5f, b.y0 + 0.5f), ImVec2(b.x1 - 0.5f, b.y1 - 0.5f), edge);
            } else {
                dl->AddRect(ImVec2(b.x0, b.y0), ImVec2(b.x1, b.y1), edge);
            }

            if (b.x1 - b.x0 > 10.0f) {
                const TimelineRangePart lit = lit_part(r);
                const float gh = MAX(b.y1 - b.y0 - 8.0f, 4.0f);
                const float gy = b.y0 + (b.y1 - b.y0 - gh) * 0.5f;
                dl->AddRectFilled(ImVec2(b.x0 + 2.0f, gy), ImVec2(b.x0 + 4.0f, gy + gh), lit == TimelineRangePart_Beg ? text_col : range_color(r, 0.95f));
                dl->AddRectFilled(ImVec2(b.x1 - 4.0f, gy), ImVec2(b.x1 - 2.0f, gy + gh), lit == TimelineRangePart_End ? text_col : range_color(r, 0.95f));
            }

            float lx0 = MAX(b.x0, x0) + 7.0f;
            float lx1 = MIN(b.x1, x1) - 7.0f;
            if (one_lane) {
                // The stretches the ranges drawn on top of it leave free: its name goes in the first wide enough for
                // all of it, or else the widest
                float cover[TIMELINE_MAX_RANGES][2];
                int num_cover = 0;
                for (int kk = k + 1; kk < n; ++kk) {
                    const TimelineRangeBar& o = bars[order[kk]];
                    const float h = o.marker ? (o.y1 - o.y0) * 0.42f : 0.0f;
                    const float c0 = o.x0 - h - 4.0f;
                    const float c1 = (o.marker ? o.x0 + h : o.x1) + 7.0f;
                    if (c1 <= lx0 || c0 >= lx1) continue;
                    int at = num_cover++;
                    while (at > 0 && cover[at - 1][0] > c0) {
                        cover[at][0] = cover[at - 1][0];
                        cover[at][1] = cover[at - 1][1];
                        --at;
                    }
                    cover[at][0] = c0;
                    cover[at][1] = c1;
                }
                const float need = ImGui::CalcTextSize(label).x;
                float best0 = lx0, best1 = lx0, from = lx0;
                bool fits = false;
                for (int c = 0; c <= num_cover && !fits; ++c) {
                    const float to = c < num_cover ? MIN(cover[c][0], lx1) : lx1;
                    if (to - from > best1 - best0) {
                        best0 = from;
                        best1 = to;
                    }
                    fits = to - from >= need;
                    if (fits) {
                        best0 = from;
                        best1 = to;
                    }
                    if (c < num_cover) from = MAX(from, cover[c][1]);
                }
                lx0 = best0;
                lx1 = best1;
            }
            if (lx1 - lx0 > font && b.y1 - b.y0 >= font - 2.0f) {
                dl->PushClipRect(ImVec2(lx0, b.y0), ImVec2(lx1, b.y1), true);
                dl->AddText(ImVec2(lx0, cy - font * 0.5f), text_col, label);
                dl->PopClipRect();
            }
        }
    }

    if (ui.drag.kind == Drag_Create) {
        const float px0 = frame_x(axis, MIN(ui.drag.anchor, ui.drag.current));
        const float px1 = MAX(frame_x(axis, MAX(ui.drag.anchor, ui.drag.current)), px0 + 2.0f);
        dl->AddRectFilled(ImVec2(px0, y0), ImVec2(px1, y1), ImGui::GetColorU32(ImGuiCol_Text, 0.1f));
        add_dashed_rect(dl, ImVec2(px0 + 0.5f, y0 + 0.5f), ImVec2(px1 - 0.5f, y1 - 0.5f), ImGui::GetColorU32(ImGuiCol_Text, 0.8f));
    }

    {
        const float px = frame_x(axis, cur);
        dl->AddLine(ImVec2(px, y0 - 1.0f), ImVec2(px, y1 + 1.0f), IM_COL32(255, 255, 0, 128));
    }

    if (rs->count == 0 && ui.drag.kind != Drag_Create) {
        const float h = one_lane ? height : row_h;
        dl->AddText(ImVec2(x0 + 6.0f, y0 + (h - font) * 0.5f), ImGui::GetColorU32(ImGuiCol_TextDisabled),
            "Shift-drag to mark a range, shift-click to mark a frame");
    }

    dl->PopClipRect();

    // Tooltips: what is under the pointer, or what is being dragged
    const bool editing = ui.edit_id != 0 || ui.rename_id != 0;
    if (ui.drag.kind == Drag_Create) {
        const int beg = MIN(ui.drag.anchor, ui.drag.current);
        const int end = MAX(ui.drag.anchor, ui.drag.current);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(beg == end ? "New marker" : "New range");
        tooltip_frames(axis, beg == end, beg, end);
        ImGui::EndTooltip();
    } else if (ui.drag.kind == Drag_Move || ui.drag.kind == Drag_Resize) {
        if (const TimelineRange* r = timeline_ranges_find(rs, ui.drag.id)) range_tooltip(rs, *r, axis, cur);
    } else if (ui.drag.kind == Drag_None && ui.hover_id && !editing) {
        if (const TimelineRange* r = timeline_ranges_find(rs, ui.hover_id)) range_tooltip(rs, *r, axis, cur);
    }

    // Naming a new range, or renaming one: a field under the lane, at the range
    if (ui.rename_request) {
        ui.rename_request = false;
        if (const TimelineRange* r = timeline_ranges_find(rs, ui.rename_id)) {
            int beg = 0, end = 0;
            float x = x0;
            if (timeline_range_frames(&beg, &end, *r, cur, nf)) x = CLAMP(frame_x(axis, beg), x0, MAX(x0, x1 - font * 16.0f));
            ImGui::SetNextWindowPos(ImVec2(x, y1 + 2.0f));
            ImGui::OpenPopup("##timeline_range_name");
            ui.rename_focus = true;
        } else {
            ui.rename_id = 0;
        }
    }
    if (ImGui::BeginPopup("##timeline_range_name")) {
        TimelineRange* r = timeline_ranges_find(rs, ui.rename_id);
        if (!r) {
            ImGui::CloseCurrentPopup();
        } else {
            if (ui.rename_focus) {
                ImGui::SetKeyboardFocusHere();
                ui.rename_focus = false;
            }
            ImGui::SetNextItemWidth(font * 14.0f);
            if (ImGui::InputText("##name", r->name, sizeof(r->name), ImGuiInputTextFlags_AutoSelectAll | ImGuiInputTextFlags_EnterReturnsTrue)) {
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
                str_copy_to_char_buf(r->name, sizeof(r->name), str_from_cstr(ui.rename_prev));
                ImGui::CloseCurrentPopup();
            }
            ImGui::TextDisabled("Enter keeps the name");
        }
        ImGui::EndPopup();
    } else if (ui.rename_id) {
        // Closed, by Enter or a click elsewhere: a name is never left empty
        if (TimelineRange* r = timeline_ranges_find(rs, ui.rename_id)) {
            if (!r->name[0]) str_copy_to_char_buf(r->name, sizeof(r->name), str_from_cstr(ui.rename_prev));
        }
        ui.rename_id = 0;
    }

    // Everything about one range, on right-click
    if (ImGui::BeginPopup("##timeline_range_edit")) {
        TimelineRange* r = timeline_ranges_find(rs, ui.edit_id);
        if (!r) {
            ImGui::CloseCurrentPopup();
        } else {
            ImGui::PushItemWidth(font * 14.0f);
            ImGui::InputText("Name", r->name, sizeof(r->name));
            ImGui::InputTextMultiline("Note", r->note, sizeof(r->note), ImVec2(font * 14.0f, ImGui::GetTextLineHeight() * 5.0f));
            ImGui::SetItemTooltip("Shown when the range is hovered");
            ImGui::ColorEdit3("Color", r->color.elem, ImGuiColorEditFlags_NoInputs);

            int follow = r->follow ? 1 : 0;
            bool mode_changed = ImGui::RadioButton("Fixed", &follow, 0);
            ImGui::SameLine();
            mode_changed |= ImGui::RadioButton("Follows the current frame", &follow, 1);
            if (mode_changed) timeline_range_set_follow(r, follow == 1, cur, nf);

            int ends[2] = { MIN(r->beg, r->end), MAX(r->beg, r->end) };
            if (ImGui::InputInt2(r->follow ? "Offsets" : "Frames", ends)) {
                if (!r->follow) {
                    ends[0] = CLAMP(ends[0], 0, nf - 1);
                    ends[1] = CLAMP(ends[1], 0, nf - 1);
                }
                r->beg = ends[0];
                r->end = ends[1];
            }
            if (r->follow) {
                ImGui::SetItemTooltip("From and to, in frames from the current frame: -60 0 is the last 60 frames, -25 25 is centered on it");
            } else {
                ImGui::SetItemTooltip("The first and the last frame. Equal, it marks a single frame");
            }

            bool is_filter = rs->filter_id == r->id;
            if (ImGui::Checkbox(ICON_FA_FILTER " Use as filter", &is_filter)) {
                rs->filter_id = is_filter ? r->id : 0;
            }
            ImGui::SetItemTooltip("The timeline filter follows this range: what is evaluated over the filter uses its frames");
            ImGui::PopItemWidth();

            ImGui::Separator();
            if (ImGui::DeleteButton("Delete")) {
                remove_id = r->id;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    } else {
        ui.edit_id = 0;
    }

    if (remove_id) {
        timeline_ranges_remove(rs, remove_id);
        if (ui.hover_id == remove_id) {
            ui.hover = {};
            ui.hover_id = 0;
        }
    }
}

void timeline_ranges_plot(TimelineRanges* rs, const TimelineRangeAxis& axis, double current_frame) {
    if (axis.num_frames <= 0 || !axis.frame_times) return;

    const ImVec2 plot_pos  = ImPlot::GetPlotPos();
    const ImVec2 plot_size = ImPlot::GetPlotSize();
    const int frame = ImGui::GetFrameCount();
    if (ui.plot_frame != frame) {
        ui.plot_frame = frame;
        ui.plot_x0 = plot_pos.x;
        ui.plot_x1 = plot_pos.x + plot_size.x;
    }

    const int nf = axis.num_frames;
    const float y0 = plot_pos.y;
    const float y1 = plot_pos.y + plot_size.y;
    auto px = [&](double f) { return ImPlot::PlotToPixels(frame_time(axis, f), 0.0).x; };

    int order[TIMELINE_MAX_RANGES];
    const int n = timeline_ranges_draw_order(order, rs, current_frame, nf);

    ImPlot::PushPlotClipRect();
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    for (int k = 0; k < n; ++k) {
        const TimelineRange& r = rs->items[order[k]];
        const bool hot    = is_hot(r.id);
        const bool filter = rs->filter_id == r.id;
        if (!(hot || filter || rs->shade_all)) continue;
        int beg, end;
        timeline_range_frames(&beg, &end, r, current_frame, nf);
        if (timeline_range_is_marker(r)) {
            const float x = px(beg);
            dl->AddLine(ImVec2(x, y0), ImVec2(x, y1), range_color(r, hot ? 1.0f : 0.6f), hot ? 2.0f : 1.0f);
        } else {
            const float xa = px(beg);
            const float xb = MAX(px(end), xa + 1.0f);
            dl->AddRectFilled(ImVec2(xa, y0), ImVec2(xb, y1), range_color(r, hot ? 0.16f : (filter ? 0.10f : 0.07f)));
            const ImU32 edge = range_color(r, hot ? 0.9f : 0.5f);
            dl->AddLine(ImVec2(xa, y0), ImVec2(xa, y1), edge);
            dl->AddLine(ImVec2(xb, y0), ImVec2(xb, y1), edge);
        }
    }
    if (ui.drag.kind == Drag_Create) {
        const float xa = px(MIN(ui.drag.anchor, ui.drag.current));
        const float xb = MAX(px(MAX(ui.drag.anchor, ui.drag.current)), xa + 1.0f);
        dl->AddRectFilled(ImVec2(xa, y0), ImVec2(xb, y1), ImGui::GetColorU32(ImGuiCol_Text, 0.07f));
        dl->AddLine(ImVec2(xa, y0), ImVec2(xa, y1), ImGui::GetColorU32(ImGuiCol_Text, 0.5f));
        dl->AddLine(ImVec2(xb, y0), ImVec2(xb, y1), ImGui::GetColorU32(ImGuiCol_Text, 0.5f));
    }
    ImPlot::PopPlotClipRect();

    // Shift-press in a plot starts a range there; the lane follows it from the next frame on
    if (ui.drag.kind == Drag_None && ImPlot::IsPlotHovered() && ImGui::GetIO().KeyShift && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        start_create(nearest_frame(axis, ImPlot::GetPlotMousePos().x));
    }
}

void timeline_ranges_menu(TimelineRanges* rs, double current_frame, int num_frames) {
    ImGui::SeparatorText("Lanes");
    int layout = (int)rs->layout;
    ImGui::RadioButton("One lane", &layout, (int)TimelineRangeLayout_OneLane);
    ImGui::RadioButton("Lane per range", &layout, (int)TimelineRangeLayout_LanePerRange);
    rs->layout = (TimelineRangeLayout)layout;
    ImGui::Checkbox("Shade every range in the plots", &rs->shade_all);

    ImGui::SeparatorText("Filter");
    const TimelineRange* filter = timeline_ranges_find(rs, rs->filter_id);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
    if (ImGui::BeginCombo("##filter", filter ? filter->name : "None")) {
        if (ImGui::Selectable("None", filter == nullptr)) rs->filter_id = 0;
        for (int i = 0; i < rs->count; ++i) {
            const TimelineRange& r = rs->items[i];
            int beg, end;
            char label[TIMELINE_RANGE_NAME_CAP + 48];
            if (timeline_range_frames(&beg, &end, r, current_frame, num_frames)) {
                snprintf(label, sizeof(label), "%s (%d to %d)###%u", r.name, beg, end, r.id);
            } else {
                snprintf(label, sizeof(label), "%s###%u", r.name, r.id);
            }
            if (ImGui::Selectable(label, r.id == rs->filter_id)) rs->filter_id = r.id;
        }
        ImGui::EndCombo();
    }
    ImGui::SetItemTooltip("The range the timeline filter follows");

    ImGui::SeparatorText("Marking");
    if (ImGui::BeginTable("##marking", 2, ImGuiTableFlags_SizingFixedFit)) {
        auto row = [](const char* input, const char* action) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(input);
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", action);
        };
        row("Shift-drag",   "mark a range");
        row("Shift-click",  "mark a frame");
        row("Drag",         "move, or resize by a handle");
        row("Double-click", "rename");
        row("Right-click",  "annotate, follow the current frame, filter, delete");
        row("Delete",       "remove the range under the mouse");
        row("Pull out",     "drag a range out of the lane onto a distribution");
        ImGui::EndTable();
    }
}
