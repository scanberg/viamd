#include <timeline_periods.h>

#include <serialization_utils.h>

#include <core/md_common.h>
#include <core/md_str.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// Distinct on the dark background of the plots, and none of them the yellow of the current frame
static vec4_t palette_color(int i) {
    static const float palette[][3] = {
        {0.353f, 0.682f, 0.941f},   // blue
        {0.941f, 0.639f, 0.353f},   // orange
        {0.561f, 0.820f, 0.486f},   // green
        {0.835f, 0.545f, 0.878f},   // purple
        {0.941f, 0.486f, 0.486f},   // red
        {0.373f, 0.827f, 0.776f},   // teal
    };
    const int n = (int)ARRAY_SIZE(palette);
    const float* c = palette[((i % n) + n) % n];
    return vec4_set(c[0], c[1], c[2], 1.0f);
}

static bool same_color(vec4_t a, vec4_t b) {
    return fabsf(a.x - b.x) < 1e-3f && fabsf(a.y - b.y) < 1e-3f && fabsf(a.z - b.z) < 1e-3f;
}

bool timeline_period_frames(int* out_beg, int* out_end, const TimelinePeriod& r, double current_frame, int num_frames) {
    if (num_frames <= 0) return false;
    int64_t beg = MIN(r.beg, r.end);
    int64_t end = MAX(r.beg, r.end);
    if (r.follow) {
        const int64_t cur = (int64_t)llround(current_frame);
        beg += cur;
        end += cur;
    }
    const int64_t last = num_frames - 1;
    if (end < 0 || beg > last) return false;
    if (out_beg) *out_beg = (int)CLAMP(beg, (int64_t)0, last);
    if (out_end) *out_end = (int)CLAMP(end, (int64_t)0, last);
    return true;
}

void timeline_period_normalize(TimelinePeriod* r) {
    if (r->beg > r->end) {
        const int32_t tmp = r->beg;
        r->beg = r->end;
        r->end = tmp;
    }
}

void timeline_period_set_follow(TimelinePeriod* r, bool follow, double current_frame, int num_frames) {
    if (r->follow == follow) return;
    timeline_period_normalize(r);
    if (follow) {
        // Centered on the current frame (an offset of 0), as long as it was. When it cannot be exactly centered,
        // the extra frame is on the side of the past.
        const int32_t extent = r->end - r->beg;
        r->beg = -(extent - extent / 2);
        r->end = extent / 2;
    } else {
        const int32_t cur  = (int32_t)llround(current_frame);
        const int32_t last = MAX(num_frames - 1, 0);
        r->beg = CLAMP(r->beg + cur, 0, last);
        r->end = CLAMP(r->end + cur, 0, last);
    }
    r->follow = follow;
}

int timeline_periods_index(const TimelinePeriods* rs, uint32_t id) {
    for (int i = 0; id != 0 && i < rs->count; ++i) {
        if (rs->items[i].id == id) return i;
    }
    return -1;
}

TimelinePeriod* timeline_periods_find(TimelinePeriods* rs, uint32_t id) {
    const int i = timeline_periods_index(rs, id);
    return i != -1 ? &rs->items[i] : nullptr;
}

const TimelinePeriod* timeline_periods_find(const TimelinePeriods* rs, uint32_t id) {
    const int i = timeline_periods_index(rs, id);
    return i != -1 ? &rs->items[i] : nullptr;
}

void timeline_periods_default_name(char* buf, size_t cap, const TimelinePeriods* rs, bool marker) {
    const char* stem = marker ? "Marker" : "Period";
    for (int k = 1; ; ++k) {
        snprintf(buf, cap, "%s %d", stem, k);
        bool taken = false;
        for (int i = 0; i < rs->count && !taken; ++i) {
            taken = strcmp(rs->items[i].name, buf) == 0;
        }
        if (!taken) return;
    }
}

TimelinePeriod* timeline_periods_add(TimelinePeriods* rs, const char* name, int32_t beg, int32_t end, bool follow) {
    if (rs->count >= TIMELINE_MAX_PERIODS) return nullptr;

    // The colour of the palette the fewest periods have, the first of those on a tie
    int best = 0;
    int best_uses = INT32_MAX;
    for (int c = 0; c < 6; ++c) {
        int uses = 0;
        for (int i = 0; i < rs->count; ++i) {
            uses += same_color(rs->items[i].color, palette_color(c)) ? 1 : 0;
        }
        if (uses < best_uses) {
            best_uses = uses;
            best = c;
        }
    }

    TimelinePeriod r = {};
    r.id     = rs->next_id++;
    r.color  = palette_color(best);
    r.follow = follow;
    r.beg    = MIN(beg, end);
    r.end    = MAX(beg, end);
    if (name && name[0]) {
        str_copy_to_char_buf(r.name, sizeof(r.name), str_from_cstr(name));
    } else {
        timeline_periods_default_name(r.name, sizeof(r.name), rs, r.beg == r.end);
    }

    rs->items[rs->count] = r;
    return &rs->items[rs->count++];
}

void timeline_periods_remove(TimelinePeriods* rs, uint32_t id) {
    const int i = timeline_periods_index(rs, id);
    if (i == -1) return;
    for (int j = i; j + 1 < rs->count; ++j) {
        rs->items[j] = rs->items[j + 1];
    }
    rs->count -= 1;
    if (rs->legacy_filter_id == id) rs->legacy_filter_id = 0;
}

void timeline_periods_clear(TimelinePeriods* rs) {
    rs->count = 0;
    rs->legacy_filter_id = 0;
}

int timeline_periods_draw_order(int* out, const TimelinePeriods* rs, double current_frame, int num_frames) {
    int len[TIMELINE_MAX_PERIODS];
    int n = 0;
    for (int i = 0; i < rs->count; ++i) {
        int beg, end;
        if (!timeline_period_frames(&beg, &end, rs->items[i], current_frame, num_frames)) continue;
        len[i] = end - beg;
        out[n++] = i;
    }
    // Insertion sort, stable: the longest first, and of equal ones the older, which comes first in the list
    for (int a = 1; a < n; ++a) {
        const int idx = out[a];
        int b = a;
        while (b > 0 && len[out[b - 1]] < len[idx]) {
            out[b] = out[b - 1];
            --b;
        }
        out[b] = idx;
    }
    return n;
}

void timeline_periods_nesting(int* out_depth, const int* order, int count, const TimelinePeriods* rs, double current_frame, int num_frames) {
    int beg[TIMELINE_MAX_PERIODS], end[TIMELINE_MAX_PERIODS];
    for (int k = 0; k < count; ++k) {
        timeline_period_frames(&beg[k], &end[k], rs->items[order[k]], current_frame, num_frames);
    }
    for (int k = 0; k < count; ++k) {
        int depth = 0;
        for (int j = 0; j < k; ++j) {
            // Periods which only touch, one ending where the next begins, are side by side and not nested
            if (beg[j] < end[k] && end[j] > beg[k]) depth += 1;
        }
        out_depth[k] = depth;
    }
}

TimelinePeriodHit timeline_periods_hit_test(const TimelinePeriodBar* bars, const int* top_first, int count, float px, float py, float handle_reach, float marker_reach) {
    const float slack_y = 2.0f;
    auto within_y = [&](const TimelinePeriodBar& b) { return b.visible && b.y0 - slack_y <= py && py <= b.y1 + slack_y; };

    TimelinePeriodHit hit;

    // Markers: the nearest one in reach
    float best_d = FLT_MAX;
    for (int k = 0; k < count; ++k) {
        const int i = top_first[k];
        const TimelinePeriodBar& b = bars[i];
        if (!b.marker || !within_y(b)) continue;
        const float d = fabsf(px - b.x0);
        if (d <= marker_reach && d < best_d) {
            best_d = d;
            hit = { i, TimelinePeriodPart_Body };
        }
    }
    if (hit.index != -1) return hit;

    // Handles. A handle reaches into its bar by at most a third of the bar, so that the middle of a
    // narrow bar can still be grabbed to move it. Pointing at a handle from outside its period costs
    // a full reach, which is what lets two handles at the same place be told apart by the side.
    float best_score = FLT_MAX;
    for (int k = 0; k < count; ++k) {
        const int i = top_first[k];
        const TimelinePeriodBar& b = bars[i];
        if (b.marker || !within_y(b)) continue;
        const float reach = CLAMP((b.x1 - b.x0) / 3.0f, 3.0f, handle_reach);
        const struct { TimelinePeriodPart part; float x; } ends[2] = { { TimelinePeriodPart_Beg, b.x0 }, { TimelinePeriodPart_End, b.x1 } };
        for (const auto& e : ends) {
            const float d = fabsf(px - e.x);
            if (d > reach) continue;
            const bool inside = e.part == TimelinePeriodPart_Beg ? px >= e.x : px <= e.x;
            const float score = d + (inside ? 0.0f : reach) + (float)k * 0.01f;
            if (score < best_score) {
                best_score = score;
                hit = { i, e.part };
            }
        }
    }
    if (hit.index != -1) return hit;

    for (int k = 0; k < count; ++k) {
        const int i = top_first[k];
        const TimelinePeriodBar& b = bars[i];
        if (b.marker || !within_y(b)) continue;
        if (b.x0 <= px && px <= b.x1) return { i, TimelinePeriodPart_Body };
    }
    return hit;
}

void timeline_periods_serialize(viamd::serialization_state_t& state, const TimelinePeriods* rs) {
    for (int i = 0; i < rs->count; ++i) {
        const TimelinePeriod& r = rs->items[i];
        viamd::write_section_header(state, STR_LIT("TimelinePeriod"));
        viamd::write_int(state, STR_LIT("Id"), r.id);
        viamd::write_str(state, STR_LIT("Name"), str_from_cstr(r.name));
        if (r.note[0]) {
            viamd::write_str(state, STR_LIT("Note"), str_from_cstr(r.note));
        }
        viamd::write_vec4(state, STR_LIT("Color"), r.color);
        viamd::write_bool(state, STR_LIT("Follow"), r.follow);
        const int frames[2] = { r.beg, r.end };
        viamd::write_int_vec(state, STR_LIT("Frames"), frames, 2);
    }
}

bool timeline_periods_deserialize_period(viamd::deserialization_state_t& state, TimelinePeriods* rs) {
    TimelinePeriod r = {};
    bool has_color = false;
    str_t ident, arg;
    int id = 0;
    while (viamd::next_entry(ident, arg, state)) {
        if (str_eq(ident, STR_LIT("Id"))) {
            viamd::extract_int(id, arg);
        } else if (str_eq(ident, STR_LIT("Name"))) {
            viamd::extract_to_char_buf(r.name, sizeof(r.name), arg);
        } else if (str_eq(ident, STR_LIT("Note"))) {
            viamd::extract_to_char_buf(r.note, sizeof(r.note), arg);
        } else if (str_eq(ident, STR_LIT("Color"))) {
            has_color = viamd::extract_vec4(r.color, arg);
        } else if (str_eq(ident, STR_LIT("Follow"))) {
            viamd::extract_bool(r.follow, arg);
        } else if (str_eq(ident, STR_LIT("Frames"))) {
            int frames[2];
            if (viamd::extract_int_vec(frames, 2, arg)) {
                r.beg = frames[0];
                r.end = frames[1];
            }
        }
    }

    TimelinePeriod* dst = timeline_periods_add(rs, r.name, r.beg, r.end, r.follow);
    if (!dst) return false;
    // The id it was saved with, which what refers to it was saved with too, unless another period has it already
    if (id > 0 && timeline_periods_index(rs, (uint32_t)id) == -1) {
        dst->id = (uint32_t)id;
        rs->next_id = MAX(rs->next_id, (uint32_t)id + 1);
    }
    str_copy_to_char_buf(dst->note, sizeof(dst->note), str_from_cstr(r.note));
    if (has_color) dst->color = r.color;
    return true;
}
