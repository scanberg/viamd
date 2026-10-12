#include "movie_overlay.h"

#include <algorithm>
#include <math.h>

bool movie_marker_matches_subplot(const MovieMarker& marker, uint32_t subplot) {
    return marker.subplot == 0 || marker.subplot == subplot;
}

void movie_marker_color(const MovieMarker* markers, size_t n, size_t i, float out[4]) {
    // Light colors that read on the dark plate of an overlay and in the editor's lanes
    static const float palette[][3] = {
        {1.00f, 0.86f, 0.35f}, // yellow
        {0.40f, 0.85f, 1.00f}, // sky blue
        {1.00f, 0.55f, 0.30f}, // orange
        {0.55f, 0.95f, 0.50f}, // green
        {1.00f, 0.55f, 0.80f}, // pink
        {0.75f, 0.65f, 1.00f}, // lavender
        {0.35f, 0.95f, 0.85f}, // turquoise
        {1.00f, 0.45f, 0.45f}, // red
    };
    const size_t num = sizeof(palette) / sizeof(palette[0]);
    if (i >= n) {
        out[0] = out[1] = out[2] = out[3] = 1.0f;
        return;
    }
    const MovieMarker& mk = markers[i];
    if (mk.color[3] > 0.0f) {
        for (int c = 0; c < 4; ++c) out[c] = mk.color[c];
        return;
    }
    size_t rank = 0;
    for (size_t j = 0; j < n; ++j) {
        if (markers[j].time < mk.time || (markers[j].time == mk.time && j < i)) rank += 1;
    }
    const float* p = palette[rank % num];
    out[0] = p[0]; out[1] = p[1]; out[2] = p[2]; out[3] = 1.0f;
}

int movie_distribution_bins(const MovieOverlay& overlay, int source_bins) {
    return overlay.num_bins == 0 ? std::max(source_bins, 8)
        : std::clamp(overlay.num_bins, MOVIE_DISTRIBUTION_MIN_BINS, MOVIE_DISTRIBUTION_MAX_BINS);
}

MovieOverlay movie_overlay_default_logo() {
    MovieOverlay o;
    o.type = MovieOverlayType::Logo;
    o.begin = 0.0;
    o.end = 3600.0;   // The longest a movie can be, so it follows the length of the movie
    o.fade_in = 0.0f;
    o.fade_out = 0.0f;
    o.anchor = MovieOverlayAnchor::TopLeft;
    o.size = 0.08f;
    return o;
}

float movie_overlay_size_px(const MovieOverlay& o, float frame_height_px) {
    return o.size_unit == MovieOverlaySizeUnit::Points ? o.size * frame_height_px / MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT
                                                       : o.size * frame_height_px;
}

float movie_overlay_convert_size(float size, MovieOverlaySizeUnit from, MovieOverlaySizeUnit to) {
    if (from == to) return size;
    return to == MovieOverlaySizeUnit::Points ? size * MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT : size / MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT;
}

void movie_overlay_size_range(MovieOverlaySizeUnit unit, float* lo, float* hi) {
    if (unit == MovieOverlaySizeUnit::Points) {
        *lo = 6.0f;
        *hi = 400.0f;
    } else {
        *lo = 0.01f;
        *hi = 0.6f;
    }
}

const char* movie_plot_palette_name(int palette) {
    static const char* names[MOVIE_PLOT_PALETTE_COUNT] = {"Colours of the plots", "Deep", "Dark", "Pastel", "Paired"};
    return names[palette < 0 || palette >= MOVIE_PLOT_PALETTE_COUNT ? 0 : palette];
}

bool movie_figure_x_labels(const MoviePlotView* views, size_t n, size_t i) {
    if (i + 1 >= n || views[i] == MoviePlotView::Distribution) return true;
    return views[i + 1] != MoviePlotView::Timeline;
}

void movie_overlay_plot_defaults(MovieOverlay* o) {
    o->size_unit = MovieOverlaySizeUnit::FrameHeight;
    if (o->type == MovieOverlayType::Timeline) {
        o->anchor = MovieOverlayAnchor::BottomCenter;
        o->width = 0.6f;
        o->size = 0.2f;
    } else if (o->type == MovieOverlayType::Distribution) {
        o->anchor = MovieOverlayAnchor::MiddleRight;
        o->width = 0.22f;
        o->size = 0.4f;
    }
}

void movie_overlays_migrate(std::vector<MovieOverlay>* overlays, const uint32_t* timeline_ids, int num_timeline, const uint32_t* distribution_ids, int num_distribution) {
    std::vector<MovieOverlay> out;
    for (MovieOverlay o : *overlays) {
        if ((o.type == MovieOverlayType::Timeline || o.type == MovieOverlayType::Distribution) && o.legacy_subplot_mask != 0) {
            const bool timeline = o.type == MovieOverlayType::Timeline;
            const uint32_t* ids = timeline ? timeline_ids : distribution_ids;
            const int n = timeline ? num_timeline : num_distribution;
            for (int i = 0; i < n && i < 31; ++i) {
                if ((o.legacy_subplot_mask >> i) & 1) o.panels.push_back({timeline ? MoviePlotView::Timeline : MoviePlotView::Distribution, ids[i]});
            }
            o.legacy_subplot_mask = 0;
        }
        if (o.type != MovieOverlayType::Figure) {
            out.push_back(o);
            continue;
        }
        MovieOverlay t = o, d = o;
        t.type = MovieOverlayType::Timeline;
        d.type = MovieOverlayType::Distribution;
        t.panels.clear();
        d.panels.clear();
        for (const MoviePlotPanel& p : o.panels) (p.view == MoviePlotView::Timeline ? t : d).panels.push_back(p);
        movie_overlay_plot_defaults(&t);
        movie_overlay_plot_defaults(&d);
        if (!t.panels.empty() || d.panels.empty()) out.push_back(t);
        if (!d.panels.empty()) out.push_back(d);
    }
    *overlays = out;
}

float movie_overlay_alpha(const MovieOverlay& o, double time) {
    if (!o.enabled || o.end <= o.begin || time < o.begin || time > o.end) return 0.0f;

    const double len = o.end - o.begin;
    double fi = o.fade_in  > 0.0f ? (double)o.fade_in  : 0.0;
    double fo = o.fade_out > 0.0f ? (double)o.fade_out : 0.0;
    if (fi + fo > len) {
        const double s = len / (fi + fo);
        fi *= s;
        fo *= s;
    }

    double a = 1.0;
    if (fi > 0.0 && time < o.begin + fi) a = (time - o.begin) / fi;
    if (fo > 0.0 && time > o.end - fo)   a = fmin(a, (o.end - time) / fo);
    return (float)(a < 0.0 ? 0.0 : (a > 1.0 ? 1.0 : a));
}

float movie_scale_bar_length(double units_per_pixel, double span_px, double target) {
    const double want = units_per_pixel * span_px * target;
    if (!(want > 0.0) || !isfinite(want)) return 0.0f;

    // The largest of 1, 2, 5 x 10^k that is not longer than what is wanted
    const double mag = pow(10.0, floor(log10(want)));
    const double m = want / mag;
    const double step = m >= 5.0 ? 5.0 : (m >= 2.0 ? 2.0 : 1.0);
    return (float)(step * mag);
}

void movie_frame_fit(float vw, float vh, float fw, float fh, float fill, float* x, float* y, float* w, float* h) {
    const float s = (fw > 0.0f && fh > 0.0f) ? fill * fminf(vw / fw, vh / fh) : 0.0f;
    *w = fw * s;
    *h = fh * s;
    *x = (vw - *w) * 0.5f;
    *y = (vh - *h) * 0.5f;
}

float movie_guide_fov_y(float fov_y, float view_h, float guide_h) {
    if (guide_h <= 0.0f) return fov_y;
    const float wide = 2.0f * atanf(tanf(fov_y * 0.5f) * view_h / guide_h);
    return fminf(wide, 170.0f * 3.14159265358979f / 180.0f);
}

double movie_units_per_pixel(float distance, float fov_y, float frame_height_px) {
    if (frame_height_px <= 0.0f) return 0.0;
    return 2.0 * (double)distance * tan((double)fov_y * 0.5) / (double)frame_height_px;
}

void movie_time_bar_profile(MovieTimeBarProfile* out, double duration, int samples, const std::function<double(double)>& quantity_at) {
    out->duration = duration;
    out->distance.assign((size_t)(samples < 1 ? 2 : samples + 1), 0.0);
    const int n = (int)out->distance.size() - 1;
    double prev = quantity_at(0.0);
    out->q.assign(out->distance.size(), prev);
    out->lo.assign(out->distance.size(), prev);
    out->hi.assign(out->distance.size(), prev);
    for (int i = 1; i <= n; ++i) {
        const double q = quantity_at(duration * (double)i / (double)n);
        out->distance[i] = out->distance[i - 1] + fabs(q - prev);
        out->q[i] = q;
        out->lo[i] = fmin(out->lo[i - 1], q);
        out->hi[i] = fmax(out->hi[i - 1], q);
        prev = q;
    }
}

namespace {
// Positions in the profile's samples, and what it holds between them
double profile_pos(const MovieTimeBarProfile& p, double t) {
    return fmin(fmax(t / p.duration, 0.0), 1.0) * (double)(p.q.size() - 1);
}
double profile_lerp(const std::vector<double>& v, double x) {
    const size_t i = (size_t)floor(x);
    if (i + 1 >= v.size()) return v.back();
    return v[i] + (v[i + 1] - v[i]) * (x - (double)i);
}
bool profile_valid(const MovieTimeBarProfile& p) {
    return p.q.size() >= 2 && p.q.size() == p.distance.size() && p.duration > 0.0;
}
}

void movie_elapsed_curve_between(std::vector<MovieCurvePoint>* out, const MovieTimeBarProfile& p, double t0, double t1, const float* xs, int num_samples,
    const std::function<double(int)>& value_of_sample, const std::function<double(double)>& value_at) {
    out->clear();
    if (!profile_valid(p) || num_samples < 1) return;

    // One stretch of the path from qa to qb that starts at s0 on the axis
    auto stretch = [&](double qa, double qb, double s0) {
        if (qb > qa) {
            for (const float* it = std::upper_bound(xs, xs + num_samples, (float)qa); it != xs + num_samples && (double)*it < qb; ++it) {
                out->push_back({s0 + ((double)*it - qa), value_of_sample((int)(it - xs))});
            }
        } else if (qb < qa) {
            const float* it = std::lower_bound(xs, xs + num_samples, (float)qa);
            while (it != xs) {
                --it;
                if ((double)*it <= qb) break;
                out->push_back({s0 + (qa - (double)*it), value_of_sample((int)(it - xs))});
            }
        }
        if (qb != qa) out->push_back({s0 + fabs(qb - qa), value_at(qb)});
    };

    const double x0 = profile_pos(p, t0);
    const double x1 = fmax(profile_pos(p, t1), x0);
    const double d0 = profile_lerp(p.distance, x0);
    out->push_back({0.0, value_at(profile_lerp(p.q, x0))});
    const size_t m = p.q.size() - 1;
    for (size_t i = (size_t)floor(x0); i < m && (double)i < x1; ++i) {
        const double a = fmax((double)i, x0), b = fmin((double)(i + 1), x1);
        if (b <= a) continue;
        stretch(profile_lerp(p.q, a), profile_lerp(p.q, b), profile_lerp(p.distance, a) - d0);
    }
}

void movie_elapsed_curve(std::vector<MovieCurvePoint>* out, const MovieTimeBarProfile& p, double time, const float* xs, int num_samples,
    const std::function<double(int)>& value_of_sample, const std::function<double(double)>& value_at) {
    movie_elapsed_curve_between(out, p, 0.0, time, xs, num_samples, value_of_sample, value_at);
}

bool movie_time_bar_visited_between(const MovieTimeBarProfile& p, double t0, double t1, double* lo, double* hi) {
    if (!profile_valid(p)) {
        if (p.q.empty()) return false;
        *lo = *hi = p.q[0];
        return true;
    }
    const double x0 = profile_pos(p, t0);
    const double x1 = fmax(profile_pos(p, t1), x0);
    *lo = *hi = profile_lerp(p.q, x0);
    for (size_t i = (size_t)ceil(x0); i < p.q.size() && (double)i <= x1; ++i) {
        *lo = fmin(*lo, p.q[i]);
        *hi = fmax(*hi, p.q[i]);
    }
    const double qe = profile_lerp(p.q, x1);
    *lo = fmin(*lo, qe);
    *hi = fmax(*hi, qe);
    return true;
}

void movie_panel_span(const MoviePlotPanel& p, double overlay_begin, double overlay_end, double* begin, double* end) {
    *begin = fmax(p.begin, overlay_begin);
    *end = p.end > p.begin ? fmin(p.end, overlay_end) : overlay_end;
}

float movie_panel_alpha(const MoviePlotPanel& p, const MovieOverlay& overlay, double time) {
    MovieOverlay span;
    movie_panel_span(p, overlay.begin, overlay.end, &span.begin, &span.end);
    span.fade_in = overlay.fade_in;
    span.fade_out = overlay.fade_out;
    return movie_overlay_alpha(span, time);
}

bool movie_time_bar_visited(const MovieTimeBarProfile& p, double time, double* lo, double* hi) {
    if (p.lo.size() < 2 || p.lo.size() != p.hi.size() || p.duration <= 0.0) {
        if (p.lo.empty()) return false;
        *lo = p.lo[0];
        *hi = p.hi[0];
        return true;
    }
    const double x = fmin(fmax(time / p.duration, 0.0), 1.0) * (double)(p.lo.size() - 1);
    const size_t i = (size_t)floor(x);
    if (i + 1 >= p.lo.size()) {
        *lo = p.lo.back();
        *hi = p.hi.back();
        return true;
    }
    const double f = x - (double)i;
    *lo = p.lo[i] + (p.lo[i + 1] - p.lo[i]) * f;
    *hi = p.hi[i] + (p.hi[i + 1] - p.hi[i]) * f;
    return true;
}

int movie_nice_ticks(double lo, double hi, int max_ticks, double* out, int cap, double* step_out) {
    if (step_out) *step_out = 0.0;
    if (!(hi > lo) || max_ticks < 1 || cap < 1) return 0;
    const double raw = (hi - lo) / (double)max_ticks;
    const double mag = pow(10.0, floor(log10(raw)));
    const double m = raw / mag;
    const double step = (m <= 1.0 ? 1.0 : (m <= 2.0 ? 2.0 : (m <= 5.0 ? 5.0 : 10.0))) * mag;
    if (step_out) *step_out = step;
    int n = 0;
    for (double t = ceil(lo / step - 1e-9) * step; t <= hi + step * 1e-9 && n < cap; t += step) {
        out[n++] = fabs(t) < step * 1e-9 ? 0.0 : t;
    }
    return n;
}

void movie_histogram_counts(std::vector<float>* counts, int num_bins, double v_min, double v_max, const float* x, const float* y,
    int stride, double y_scale, int num_samples, double x_lo, double x_hi) {
    counts->assign((size_t)(num_bins > 0 ? num_bins : 0), 0.0f);
    if (num_bins < 1 || !(v_max > v_min) || stride < 1) return;
    const double inv = (double)num_bins / (v_max - v_min);
    for (int i = 0; i < num_samples; ++i) {
        if ((double)x[i] < x_lo || (double)x[i] > x_hi) continue;
        const double v = (double)y[(size_t)i * (size_t)stride] * y_scale;
        if (v < v_min || v > v_max) continue;
        int b = (int)((v - v_min) * inv);
        if (b >= num_bins) b = num_bins - 1;
        (*counts)[b] += 1.0f;
    }
}

double movie_time_bar_moved(const MovieTimeBarProfile& p, double time) {
    if (p.distance.size() < 2 || p.duration <= 0.0) return 0.0;
    const double x = fmin(fmax(time / p.duration, 0.0), 1.0) * (double)(p.distance.size() - 1);
    const size_t i = (size_t)floor(x);
    if (i + 1 >= p.distance.size()) return p.distance.back();
    return p.distance[i] + (p.distance[i + 1] - p.distance[i]) * (x - (double)i);
}

double movie_time_bar_progress(const MovieTimeBarProfile& p, double time) {
    const double total = p.total();
    return total > 0.0 ? movie_time_bar_moved(p, time) / total : 0.0;
}

double movie_quantity_speed(const std::function<double(double)>& quantity_at, double time, double step) {
    if (step <= 0.0) return 0.0;
    return fabs(quantity_at(time + step) - quantity_at(time - step)) / (2.0 * step);
}

void movie_overlay_avoid(const MovieOverlayBox* boxes, size_t n, float gap, float frame_y0, float frame_y1, float* dy_out) {
    if (!boxes || !dy_out) return;
    auto rank = [&](const MovieOverlayBox& b) { return (b.anchor / 3 == 1 ? 2 : 0) + (b.anchor % 3 == 1 ? 0 : 1); };
    std::vector<size_t> order;
    order.reserve(n);
    for (int r = 0; r < 4; ++r) {
        for (size_t i = 0; i < n; ++i) {
            if (rank(boxes[i]) == r) order.push_back(i);
        }
    }
    std::vector<size_t> placed;
    placed.reserve(n);
    // Where the box ends up moving only up or only down, or false when it leaves the frame on the way
    auto push = [&](const MovieOverlayBox& b, bool up, float* dy_res) {
        float dy = 0.0f;
        // Each step clears at least one placed box, which it cannot hit again in the same direction
        for (size_t step = 0; step <= placed.size(); ++step) {
            bool moved = false;
            for (size_t j : placed) {
                const MovieOverlayBox& p = boxes[j];
                const float py0 = p.y0 + dy_out[j], py1 = p.y1 + dy_out[j];
                const bool same_time = b.begin < p.end && p.begin < b.end;
                const bool overlap = b.x0 < p.x1 && p.x0 < b.x1 && b.y0 + dy < py1 && py0 < b.y1 + dy;
                if (!same_time || !overlap) continue;
                dy = up ? py0 - gap - b.y1 : py1 + gap - b.y0;
                moved = true;
            }
            if (!moved) break;
        }
        *dy_res = dy;
        return b.y0 + dy >= frame_y0 - 0.5f && b.y1 + dy <= frame_y1 + 0.5f;
    };
    for (size_t i : order) {
        const MovieOverlayBox& b = boxes[i];
        const int row = b.anchor / 3;
        float dy_up = 0.0f, dy_down = 0.0f;
        const bool ok_up = push(b, true, &dy_up);
        const bool ok_down = push(b, false, &dy_down);
        float dy = 0.0f;
        if (row == 1) {
            if (ok_up && (!ok_down || -dy_up < dy_down)) dy = dy_up;
            else if (ok_down) dy = dy_down;
        } else {
            const bool prefer_up = row == 2;
            if (prefer_up ? ok_up : ok_down) dy = prefer_up ? dy_up : dy_down;
            else if (prefer_up ? ok_down : ok_up) dy = prefer_up ? dy_down : dy_up;
        }
        dy_out[i] = dy;
        placed.push_back(i);
    }
}

int movie_label_place(const MovieLabelSpan* placed, size_t n, float x, float w, float lo, float hi, float offset, float gap, int max_rows, float* x0_out) {
    max_rows = max_rows < 1 ? 1 : max_rows;
    const bool prefer_left = x > 0.5f * (lo + hi);
    float start[2];
    for (int side = 0; side < 2; ++side) {
        const bool left = (side == 0) == prefer_left;
        float lx = left ? x - offset - w : x + offset;
        if (lx + w > hi) lx = hi - w;
        if (lx < lo) lx = lo;
        start[side] = lx;
    }
    auto overlap = [&](float lx, int row) {
        float sum = 0.0f;
        for (size_t j = 0; j < n; ++j) {
            if (placed[j].row != row) continue;
            const float a = lx > placed[j].x0 - gap ? lx : placed[j].x0 - gap;
            const float b = lx + w < placed[j].x1 + gap ? lx + w : placed[j].x1 + gap;
            if (b > a) sum += b - a;
        }
        return sum;
    };
    int best_row = 0;
    float best_x = start[0], best = -1.0f;
    for (int side = 0; side < 2; ++side) {
        for (int row = 0; row < max_rows; ++row) {
            const float o = overlap(start[side], row);
            if (o <= 0.0f) {
                *x0_out = start[side];
                return row;
            }
            if (best < 0.0f || o < best) {
                best = o;
                best_row = row;
                best_x = start[side];
            }
        }
    }
    *x0_out = best_x;
    return best_row;
}
