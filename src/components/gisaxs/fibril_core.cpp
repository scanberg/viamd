#include "fibril_core.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <algorithm>

namespace {

struct v3 { double x, y, z; };
inline v3 operator+(v3 a, v3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline v3 operator-(v3 a, v3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline v3 operator*(v3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(v3 a, v3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline v3 cross(v3 a, v3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double length(v3 a) { return sqrt(dot(a, a)); }
inline v3 normalize(v3 a) { const double l = length(a); return l > 0.0 ? a * (1.0 / l) : v3{0, 0, 0}; }

inline double min_image_1d(double d, double L) {
    return L > 0.0 ? d - L * floor(d / L + 0.5) : d;
}
inline v3 min_image(v3 d, const double box[3]) {
    return {min_image_1d(d.x, box[0]), min_image_1d(d.y, box[1]), min_image_1d(d.z, box[2])};
}

inline v3 any_perpendicular(v3 t) {
    const v3 a = fabs(t.z) < 0.9 ? v3{0, 0, 1} : v3{1, 0, 0};
    return normalize(cross(a, t));
}

// Per slice geometry
struct SliceFrame {
    v3 c;           // centroid (unwrapped along the chain)
    v3 t;           // axis
    v3 e1, e2;      // transported in-plane basis
    double theta;   // rotation of the template in (e1, e2)
    double E;       // electrons
};

// Catmull-Rom position and derivative
inline v3 cr_pos(v3 p0, v3 p1, v3 p2, v3 p3, double u) {
    const double u2 = u * u, u3 = u2 * u;
    return (p1 * 2.0 + (p2 - p0) * u + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * u2 + (p1 * 3.0 - p0 - p2 * 3.0 + p3) * u3) * 0.5;
}
inline v3 cr_der(v3 p0, v3 p1, v3 p2, v3 p3, double u) {
    const double u2 = u * u;
    return ((p2 - p0) + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * (2.0 * u) + (p1 * 3.0 - p0 - p2 * 3.0 + p3) * (3.0 * u2)) * 0.5;
}

inline double wrap_angle_near(double a, double ref) {
    const double twopi = 6.283185307179586;
    return a - twopi * floor((a - ref) / twopi + 0.5);
}

// Computes the frames for all slices of one chain
void chain_frames(std::vector<SliceFrame>& fr, const FibrilInput& in, size_t s_beg, size_t s_end,
                  const float* au, const float* av, bool use_anchors, std::vector<v3>& scratch) {
    const int S = in.stride;
    const size_t n = s_end - s_beg;
    fr.resize(n);

    // Centroids (beads unwrapped relative to the first bead of the slice), then unwrapped along the chain
    for (size_t j = 0; j < n; ++j) {
        const size_t base = (s_beg + j) * S;
        const v3 b0 = {in.x[base], in.y[base], in.z[base]};
        v3 acc = {0, 0, 0};
        double E = 0.0;
        for (int b = 0; b < S; ++b) {
            const v3 p = {in.x[base + b], in.y[base + b], in.z[base + b]};
            const v3 d = min_image(p - b0, in.box);
            const double e = in.e ? in.e[base + b] : 1.0;
            acc = acc + (b0 + d) * e;
            E += e;
        }
        fr[j].E = E;
        fr[j].c = E != 0.0 ? acc * (1.0 / E) : b0;
        if (j > 0) {
            fr[j].c = fr[j - 1].c + min_image(fr[j].c - fr[j - 1].c, in.box);
        }
    }

    // Tangents
    for (size_t j = 0; j < n; ++j) {
        v3 d;
        if (n == 1) d = {1, 0, 0};
        else if (j == 0) d = fr[1].c - fr[0].c;
        else if (j == n - 1) d = fr[n - 1].c - fr[n - 2].c;
        else d = fr[j + 1].c - fr[j - 1].c;
        fr[j].t = normalize(d);
    }

    // Rotation minimizing frame (parallel transport)
    for (size_t j = 0; j < n; ++j) {
        v3 e1 = (j == 0) ? any_perpendicular(fr[0].t) : fr[j - 1].e1;
        e1 = normalize(e1 - fr[j].t * dot(e1, fr[j].t));
        if (length(e1) < 0.5) e1 = any_perpendicular(fr[j].t);
        fr[j].e1 = e1;
        fr[j].e2 = cross(fr[j].t, e1);
    }

    // Orientation
    scratch.resize(S);
    for (size_t j = 0; j < n; ++j) {
        const size_t base = (s_beg + j) * S;
        // Bead offsets in the slice plane
        double pu[64], pv[64];
        for (int b = 0; b < S && b < 64; ++b) {
            const v3 p = {in.x[base + b], in.y[base + b], in.z[base + b]};
            const v3 d = min_image(p - fr[j].c, in.box);
            pu[b] = dot(d, fr[j].e1);
            pv[b] = dot(d, fr[j].e2);
        }
        double theta;
        if (use_anchors) {
            // 2D Procrustes: rotation taking the anchors onto the observed offsets
            double sc = 0.0, ss = 0.0;
            for (int b = 0; b < S && b < 64; ++b) {
                const double w = in.e ? in.e[base + b] : 1.0;
                sc += w * (au[b] * pu[b] + av[b] * pv[b]);
                ss += w * (au[b] * pv[b] - av[b] * pu[b]);
            }
            theta = atan2(ss, sc);
        } else {
            const int r = (in.ref_bead >= 0 && in.ref_bead < S) ? in.ref_bead : 0;
            theta = atan2(pv[r], pu[r]);
        }
        fr[j].theta = (j > 0) ? wrap_angle_near(theta, fr[j - 1].theta) : theta;
    }
}

}  // namespace

bool fibril_mean_layout(std::vector<float>& out_u, std::vector<float>& out_v, std::vector<float>& out_frac, const FibrilInput& in) {
    const int S = in.stride;
    if (S <= 0 || S > 64 || in.num_slices == 0) return false;
    std::vector<double> su(S, 0.0), sv(S, 0.0), se(S, 0.0);
    size_t count = 0;
    std::vector<SliceFrame> fr;
    std::vector<v3> scratch;
    for (size_t c = 0; c < in.num_chains; ++c) {
        const size_t beg = in.chain_offset[c], end = in.chain_offset[c + 1];
        if (end - beg < 2) continue;
        chain_frames(fr, in, beg, end, nullptr, nullptr, false, scratch);
        for (size_t j = 0; j < fr.size(); ++j) {
            const size_t base = (beg + j) * S;
            const double ct = cos(-fr[j].theta), st = sin(-fr[j].theta);
            for (int b = 0; b < S; ++b) {
                const v3 p = {in.x[base + b], in.y[base + b], in.z[base + b]};
                const v3 d = min_image(p - fr[j].c, in.box);
                const double u = dot(d, fr[j].e1), v = dot(d, fr[j].e2);
                su[b] += ct * u - st * v;
                sv[b] += st * u + ct * v;
                se[b] += in.e ? in.e[base + b] : 1.0;
            }
            count += 1;
        }
    }
    if (!count) return false;
    out_u.resize(S); out_v.resize(S); out_frac.resize(S);
    double etot = 0.0;
    for (int b = 0; b < S; ++b) etot += se[b];
    for (int b = 0; b < S; ++b) {
        out_u[b] = (float)(su[b] / count);
        out_v[b] = (float)(sv[b] / count);
        out_frac[b] = (float)(etot != 0.0 ? se[b] / etot : 1.0 / S);
    }
    // Center on the weighted centroid
    double cu = 0.0, cv = 0.0;
    for (int b = 0; b < S; ++b) { cu += out_frac[b] * out_u[b]; cv += out_frac[b] * out_v[b]; }
    for (int b = 0; b < S; ++b) { out_u[b] -= (float)cu; out_v[b] -= (float)cv; }
    return true;
}

FibrilTemplate fibril_template_beads(const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& frac, float sigma) {
    FibrilTemplate t;
    t.anchor_u = u;
    t.anchor_v = v;
    double sum = 0.0;
    for (float f : frac) sum += f;
    for (size_t i = 0; i < u.size(); ++i) {
        t.gauss.push_back({u[i], v[i], (float)(sum > 0.0 ? frac[i] / sum : 1.0 / u.size()), sigma});
    }
    return t;
}

FibrilTemplate fibril_template_disk(const std::vector<float>& anchor_u, const std::vector<float>& anchor_v, float radius, float spacing) {
    FibrilTemplate t;
    t.anchor_u = anchor_u;
    t.anchor_v = anchor_v;
    spacing = std::max(spacing, 0.5f);
    const float sigma = 0.5f * spacing;
    const float dy = spacing * 0.8660254f;
    const int ny = (int)ceil(radius / dy) + 1;
    const int nx = (int)ceil(radius / spacing) + 1;
    for (int j = -ny; j <= ny; ++j) {
        const float off = (j & 1) ? 0.5f * spacing : 0.0f;
        for (int i = -nx - 1; i <= nx + 1; ++i) {
            const float u = i * spacing + off;
            const float v = j * dy;
            if (u * u + v * v <= radius * radius) {
                t.gauss.push_back({u, v, 1.0f, sigma});
            }
        }
    }
    if (t.gauss.empty()) t.gauss.push_back({0.0f, 0.0f, 1.0f, std::max(sigma, radius * 0.5f)});
    const float w = 1.0f / (float)t.gauss.size();
    for (auto& g : t.gauss) g.weight = w;
    return t;
}

bool fibril_template_load(FibrilTemplate* out, const char* path, char* err, size_t err_cap) {
    FILE* f = fopen(path, "r");
    if (!f) {
        snprintf(err, err_cap, "Could not open template file '%s'", path);
        return false;
    }
    FibrilTemplate t;
    char line[512];
    int line_no = 0;
    bool ok = true;
    while (fgets(line, sizeof(line), f)) {
        ++line_no;
        char* hash = strchr(line, '#');
        if (hash) *hash = '\0';
        char key[32] = "";
        float a, b, c, d;
        const int n = sscanf(line, "%31s %f %f %f %f", key, &a, &b, &c, &d);
        if (n <= 0) continue;
        if (strcmp(key, "anchor") == 0 && n >= 3) {
            t.anchor_u.push_back(a);
            t.anchor_v.push_back(b);
        } else if (strcmp(key, "gauss") == 0 && n >= 5) {
            if (d <= 0.0f) { snprintf(err, err_cap, "Line %i: sigma must be positive", line_no); ok = false; break; }
            t.gauss.push_back({a, b, c, d});
        } else {
            snprintf(err, err_cap, "Line %i: unrecognized entry", line_no);
            ok = false;
            break;
        }
    }
    fclose(f);
    if (!ok) return false;
    if (t.gauss.empty()) {
        snprintf(err, err_cap, "Template file contains no 'gauss' entries");
        return false;
    }
    double sum = 0.0;
    for (auto& g : t.gauss) sum += g.weight;
    if (sum <= 0.0) {
        snprintf(err, err_cap, "Template weights must sum to a positive value");
        return false;
    }
    for (auto& g : t.gauss) g.weight = (float)(g.weight / sum);
    *out = t;
    return true;
}

static double procrustes_rms(const std::vector<float>& au, const std::vector<float>& av, const std::vector<float>& bu, const std::vector<float>& bv, double mirror) {
    const size_t n = std::min(au.size(), bu.size());
    if (n == 0) return DBL_MAX;
    double sc = 0.0, ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double u = au[i], v = mirror * av[i];
        sc += u * bu[i] + v * bv[i];
        ss += u * bv[i] - v * bu[i];
    }
    const double th = atan2(ss, sc), c = cos(th), sn = sin(th);
    double r2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double u = au[i], v = mirror * av[i];
        const double du = c * u - sn * v - bu[i], dv = sn * u + c * v - bv[i];
        r2 += du * du + dv * dv;
    }
    return sqrt(r2 / n);
}

bool fibril_template_match_handedness(FibrilTemplate& tmpl, const std::vector<float>& ref_u, const std::vector<float>& ref_v, double* out_rms) {
    const double r_same = procrustes_rms(tmpl.anchor_u, tmpl.anchor_v, ref_u, ref_v, 1.0);
    const double r_mirr = procrustes_rms(tmpl.anchor_u, tmpl.anchor_v, ref_u, ref_v, -1.0);
    const bool mirror = r_mirr < r_same;
    if (mirror) {
        for (auto& v : tmpl.anchor_v) v = -v;
        for (auto& g : tmpl.gauss) g.v = -g.v;
    }
    if (out_rms) *out_rms = mirror ? r_mirr : r_same;
    return mirror;
}

bool fibril_sweep(FibrilOutput* out, const FibrilInput& in, const FibrilTemplate& tmpl, double ds, char* err, size_t err_cap) {
    const int S = in.stride;
    if (S <= 0 || S > 64) { snprintf(err, err_cap, "Invalid number of beads per slice"); return false; }
    if (tmpl.gauss.empty()) { snprintf(err, err_cap, "Empty cross-section template"); return false; }
    const bool use_anchors = (int)tmpl.anchor_u.size() == S && (int)tmpl.anchor_v.size() == S;

    double sigma_min = DBL_MAX;
    for (const auto& g : tmpl.gauss) sigma_min = std::min(sigma_min, (double)g.sigma);
    if (ds <= 0.0) ds = sigma_min;
    ds = std::max(ds, 0.1);

    out->x.clear(); out->y.clear(); out->z.clear(); out->w.clear(); out->s.clear();
    FibrilStats st = {};

    std::vector<SliceFrame> fr;
    std::vector<v3> scratch;
    double spacing_sum = 0.0; size_t spacing_cnt = 0;
    double r2_sum = 0.0, a2_sum = 0.0; size_t bead_cnt = 0;

    const size_t G = tmpl.gauss.size();
    auto emit = [&](v3 P, v3 e1, v3 e2, double theta, double weight) {
        const double ct = cos(theta), sn = sin(theta);
        for (size_t g = 0; g < G; ++g) {
            const FibrilGauss& gs = tmpl.gauss[g];
            const double u = ct * gs.u - sn * gs.v;
            const double v = sn * gs.u + ct * gs.v;
            const v3 p = P + e1 * u + e2 * v;
            out->x.push_back((float)p.x);
            out->y.push_back((float)p.y);
            out->z.push_back((float)p.z);
            out->w.push_back((float)(weight * gs.weight));
            out->s.push_back(gs.sigma);
            st.electrons_out += weight * gs.weight;
        }
    };

    for (size_t c = 0; c < in.num_chains; ++c) {
        const size_t beg = in.chain_offset[c], end = in.chain_offset[c + 1];
        if (end <= beg) continue;
        chain_frames(fr, in, beg, end, use_anchors ? tmpl.anchor_u.data() : nullptr, use_anchors ? tmpl.anchor_v.data() : nullptr, use_anchors, scratch);
        const size_t n = fr.size();

        for (size_t j = 0; j < n; ++j) {
            st.electrons_in += fr[j].E;
            const size_t base = (beg + j) * S;
            for (int b = 0; b < S; ++b) {
                const v3 p = {in.x[base + b], in.y[base + b], in.z[base + b]};
                const v3 d = min_image(p - fr[j].c, in.box);
                const double a = dot(d, fr[j].t);
                r2_sum += dot(d, d) - a * a;
                a2_sum += a * a;
                bead_cnt += 1;
            }
        }

        if (n == 1) {
            // Isolated slice: no axis information, place the template once
            emit(fr[0].c, fr[0].e1, fr[0].e2, fr[0].theta, fr[0].E);
            continue;
        }

        // Interior segments
        for (size_t j = 0; j + 1 < n; ++j) {
            const v3 p0 = (j == 0) ? fr[0].c * 2.0 - fr[1].c : fr[j - 1].c;
            const v3 p1 = fr[j].c;
            const v3 p2 = fr[j + 1].c;
            const v3 p3 = (j + 2 < n) ? fr[j + 2].c : fr[n - 1].c * 2.0 - fr[n - 2].c;
            const double L = length(p2 - p1);
            spacing_sum += L; spacing_cnt += 1;
            const int ns = std::max(1, (int)ceil(L / ds));
            const double wseg = 0.5 * (fr[j].E + fr[j + 1].E) / ns;
            for (int m = 0; m < ns; ++m) {
                const double u = (m + 0.5) / ns;
                const v3 P = cr_pos(p0, p1, p2, p3, u);
                const v3 T = normalize(cr_der(p0, p1, p2, p3, u));
                v3 e1 = fr[j].e1 * (1.0 - u) + fr[j + 1].e1 * u;
                e1 = normalize(e1 - T * dot(e1, T));
                const v3 e2 = cross(T, e1);
                const double theta = fr[j].theta * (1.0 - u) + fr[j + 1].theta * u;
                emit(P, e1, e2, theta, wseg);
            }
        }

        // End caps: half a spacing beyond the first and last slice, carrying the outer half of their electrons
        for (int side = 0; side < 2; ++side) {
            const size_t j = side == 0 ? 0 : n - 1;
            const size_t k = side == 0 ? 1 : n - 2;
            const v3 dir = normalize(fr[j].c - fr[k].c);
            const double Lh = 0.5 * length(fr[j].c - fr[k].c);
            const int ns = std::max(1, (int)ceil(Lh / ds));
            const double wseg = 0.5 * fr[j].E / ns;
            for (int m = 0; m < ns; ++m) {
                const double u = (m + 0.5) / ns;
                const v3 P = fr[j].c + dir * (u * Lh);
                emit(P, fr[j].e1, fr[j].e2, fr[j].theta, wseg);
            }
        }
    }

    st.num_points = out->x.size();
    st.mean_spacing = spacing_cnt ? spacing_sum / spacing_cnt : 0.0;
    st.mean_radius = bead_cnt ? sqrt(r2_sum / bead_cnt) : 0.0;
    st.mean_axial_offset = bead_cnt ? sqrt(a2_sum / bead_cnt) : 0.0;
    out->stats = st;
    return true;
}
