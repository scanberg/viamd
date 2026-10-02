#include "scattering_core.h"

#include <float.h>
#include <math.h>
#include <string.h>

namespace gisaxs {

namespace {
inline double clampd(double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); }
constexpr double FWHM_TO_SIGMA = 1.0 / 2.354820045030949;
}

// -------------------------------------------------------------------------------------------------
// Optics and geometry
// -------------------------------------------------------------------------------------------------

double critical_angle(double d_sld, double lambda) {
    if (!(d_sld > 0.0) || !(lambda > 0.0)) return 0.0;
    // k0 sin(alpha_c) = sqrt(4 pi d_sld)
    return asin(clampd(sqrt(4.0 * kPi * d_sld) * lambda / (2.0 * kPi), 0.0, 1.0));
}

double Beam::k0() const { return 2.0 * kPi / wavelength; }
double Beam::p() const  { return k0() * sin(alpha_i); }

double qz_from_alpha_f(const Beam& beam, double alpha_f) {
    return beam.k0() * (sin(alpha_f) + sin(beam.alpha_i));
}

double alpha_f_from_qz(const Beam& beam, double qz) {
    return asin(clampd(qz / beam.k0() - sin(beam.alpha_i), -1.0, 1.0));
}

double yoneda_qz(const Beam& beam, double d_sld) {
    return qz_from_alpha_f(beam, critical_angle(d_sld, beam.wavelength));
}

// -------------------------------------------------------------------------------------------------
// Setup -> model
// -------------------------------------------------------------------------------------------------

double contrast_amplitude(const Setup& s) {
    return s.rho_material > 0.0 ? 1.0 - s.rho_ambient / s.rho_material : 1.0;
}

double substrate_sld(const Setup& s) {
    return sld_from_electron_density(s.rho_substrate - s.rho_ambient);
}

double substrate_absorption_sld(const Setup& s) {
    return absorption_sld(s.beta_substrate, s.wavelength);
}

md_gisaxs_model_t make_model(const Setup& s) {
    md_gisaxs_model_t m = {};
    m.wavelength = s.wavelength;
    m.alpha_i = s.alpha_i;
    m.dwba = s.dwba;
    m.graded = s.graded;
    m.z_substrate = s.z_substrate;
    m.sld_ambient = sld_from_electron_density(s.rho_ambient);
    m.sld_substrate = sld_from_electron_density(s.rho_substrate);
    m.sld_substrate_abs = absorption_sld(s.beta_substrate, s.wavelength);
    m.substrate_roughness = s.roughness;
    const double c = contrast_amplitude(s);
    // Film: laterally averaged electron density, contrast corrected like the particles themselves
    m.profile_sld_scale = MD_GISAXS_R_E * c;
    m.profile_abs_scale = s.rho_material > 0.0 ? absorption_sld(s.beta_material, s.wavelength) / s.rho_material : 0.0;
    // md_gisaxs reads an intensity scale of exactly 0 as 'unset' (1.0). A contrast matched sample must give zero.
    const double scale = MD_GISAXS_R_E * MD_GISAXS_R_E * c * c;
    m.intensity_scale = scale > 0.0 ? scale : DBL_MIN;
    return m;
}

double film_sld(const double* profile, size_t n, double profile_sld_scale) {
    if (!profile || n == 0) return 0.0;
    double vmax = 0.0;
    for (size_t i = 0; i < n; ++i) vmax = profile[i] > vmax ? profile[i] : vmax;
    if (!(vmax > 0.0)) return 0.0;
    double sum = 0.0;
    size_t cnt = 0;
    for (size_t i = 0; i < n; ++i) {
        if (profile[i] >= 0.5 * vmax) { sum += profile[i]; cnt += 1; }
    }
    return cnt ? profile_sld_scale * sum / (double)cnt : 0.0;
}

// -------------------------------------------------------------------------------------------------
// q-map
// -------------------------------------------------------------------------------------------------

void qmap_init_rings(QMap* map, const md_gisaxs_t* ctx) {
    const size_t R = md_gisaxs_num_rings(ctx);
    const double* q = md_gisaxs_ring_q(ctx);
    const double* e = md_gisaxs_ring_edges(ctx);
    const unsigned* c = md_gisaxs_ring_count(ctx);
    map->ring_q.assign(q, q + R);
    map->ring_edge.assign(e, e + R + 1);
    map->ring_count.assign(c, c + R);
    map->rows = 0;
    map->qz0 = map->dqz = 0.0;
    map->I.clear();
}

long ring_at(const QMap& map, double qpar) {
    const size_t R = map.cols();
    if (R == 0 || map.ring_edge.size() != R + 1) return -1;
    const double* e = map.ring_edge.data();
    if (!(qpar >= e[0] && qpar <= e[R])) return -1;
    // Last edge <= q_par
    size_t lo = 0, hi = R;
    while (lo < hi) {
        const size_t mid = (lo + hi + 1) / 2;
        if (e[mid] <= qpar) lo = mid; else hi = mid - 1;
    }
    return (long)(lo < R ? lo : R - 1);
}

bool ring_interp(const QMap& map, double qpar, size_t* c0, size_t* c1, double* t) {
    const size_t R = map.cols();
    if (R == 0 || map.ring_edge.size() != R + 1) return false;
    if (!(qpar >= map.ring_edge[0] && qpar <= map.ring_edge[R])) return false;
    const double* rq = map.ring_q.data();
    if (qpar <= rq[0])     { *c0 = *c1 = 0;     *t = 0.0; return true; }
    if (qpar >= rq[R - 1]) { *c0 = *c1 = R - 1; *t = 0.0; return true; }
    // Last ring center <= q_par (R >= 2 here)
    size_t lo = 0, hi = R - 1;
    while (lo < hi) {
        const size_t mid = (lo + hi + 1) / 2;
        if (rq[mid] <= qpar) lo = mid; else hi = mid - 1;
    }
    *c0 = lo;
    *c1 = lo + 1 < R ? lo + 1 : R - 1;
    *t = *c1 > *c0 ? (qpar - rq[*c0]) / (rq[*c1] - rq[*c0]) : 0.0;
    return true;
}

bool qmap_sample(const QMap& map, double qpar, double qz, float* out) {
    if (map.empty()) return false;
    size_t c0, c1;
    double tc;
    if (!ring_interp(map, qpar, &c0, &c1, &tc)) return false;
    size_t r0 = 0, r1 = 0;
    double tr = 0.0;
    if (map.rows > 1 && map.dqz > 0.0) {
        const double rf = (qz - map.qz0) / map.dqz;
        if (rf < 0.0 || rf > (double)(map.rows - 1)) return false;
        r0 = (size_t)rf;
        if (r0 > map.rows - 1) r0 = map.rows - 1;
        r1 = r0 + 1 < map.rows ? r0 + 1 : r0;
        tr = rf - (double)r0;
    } else if (fabs(qz - map.qz0) > 1.0e-12) {
        return false;
    }
    const double v = (1 - tr) * ((1 - tc) * map.at(r0, c0) + tc * map.at(r0, c1)) +
                     tr       * ((1 - tc) * map.at(r1, c0) + tc * map.at(r1, c1));
    *out = (float)v;
    return true;
}

void qmap_apply_resolution(QMap* dst, const QMap& src, double fwhm_qpar, double fwhm_qz, double qz_valid_min) {
    if (dst != &src) *dst = src;
    const size_t R = dst->cols(), Q = dst->rows;
    if (R == 0 || Q == 0) return;
    float* I = dst->I.data();

    for (size_t r = 0; r < Q; ++r) {
        if (dst->qz(r) < qz_valid_min) memset(I + r * R, 0, sizeof(float) * R);
    }

    // q_par: Gaussian on the non-uniform rings, quadrature weight = ring width, each row of the kernel normalized
    const double sq = fwhm_qpar * FWHM_TO_SIGMA;
    if (sq > 0.0 && R > 1) {
        const double cutoff = 4.0 * sq;
        std::vector<size_t> kb(R), ke(R), ko(R + 1);
        std::vector<double> kern;
        ko[0] = 0;
        for (size_t c = 0; c < R; ++c) {
            size_t b = c, e = c + 1;
            while (b > 0 && dst->ring_q[c] - dst->ring_q[b - 1] <= cutoff) --b;
            while (e < R && dst->ring_q[e] - dst->ring_q[c] <= cutoff) ++e;
            kb[c] = b;
            ke[c] = e;
            const size_t off = kern.size();
            double wsum = 0.0;
            for (size_t cc = b; cc < e; ++cc) {
                const double d = (dst->ring_q[cc] - dst->ring_q[c]) / sq;
                const double w = exp(-0.5 * d * d) * (dst->ring_edge[cc + 1] - dst->ring_edge[cc]);
                kern.push_back(w);
                wsum += w;
            }
            for (size_t i = off; i < kern.size(); ++i) kern[i] = wsum > 0.0 ? kern[i] / wsum : 0.0;
            ko[c + 1] = kern.size();
        }
        std::vector<float> tmp(R);
        for (size_t r = 0; r < Q; ++r) {
            float* row = I + r * R;
            for (size_t c = 0; c < R; ++c) {
                const double* k = kern.data() + ko[c];
                double sum = 0.0;
                for (size_t cc = kb[c]; cc < ke[c]; ++cc) sum += k[cc - kb[c]] * row[cc];
                tmp[c] = (float)sum;
            }
            memcpy(row, tmp.data(), sizeof(float) * R);
        }
    }

    // q_z: Gaussian on the uniform rows, rows below qz_valid_min excluded (they carry no intensity)
    const double sr = dst->dqz > 0.0 ? fwhm_qz * FWHM_TO_SIGMA / dst->dqz : 0.0;
    if (sr > 0.05 && Q > 1) {
        const int rad = (int)ceil(3.0 * sr);
        std::vector<double> kern(2 * rad + 1);
        for (int i = -rad; i <= rad; ++i) kern[i + rad] = exp(-0.5 * (i * i) / (sr * sr));
        std::vector<float> tmp(Q);
        for (size_t c = 0; c < R; ++c) {
            for (size_t r = 0; r < Q; ++r) {
                if (dst->qz(r) < qz_valid_min) { tmp[r] = 0.0f; continue; }
                double sum = 0.0, wsum = 0.0;
                for (int i = -rad; i <= rad; ++i) {
                    const long rr = (long)r + i;
                    if (rr < 0 || rr >= (long)Q || dst->qz((size_t)rr) < qz_valid_min) continue;
                    sum += kern[i + rad] * I[rr * R + c];
                    wsum += kern[i + rad];
                }
                tmp[r] = (float)(wsum > 0.0 ? sum / wsum : 0.0);
            }
            for (size_t r = 0; r < Q; ++r) I[r * R + c] = tmp[r];
        }
    }
}

size_t qmap_cut_horizontal(const QMap& map, double qz, double width, double* out) {
    const size_t R = map.cols();
    for (size_t c = 0; c < R; ++c) out[c] = 0.0;
    const double half = width * 0.5 > map.dqz * 0.5 ? width * 0.5 : map.dqz * 0.5;
    size_t n = 0;
    for (size_t r = 0; r < map.rows; ++r) {
        if (fabs(map.qz(r) - qz) > half) continue;
        for (size_t c = 0; c < R; ++c) out[c] += map.at(r, c);
        n += 1;
    }
    if (n) for (size_t c = 0; c < R; ++c) out[c] /= (double)n;
    return n;
}

size_t qmap_cut_vertical(const QMap& map, double qpar, double width, double* out) {
    const size_t R = map.cols(), Q = map.rows;
    const double half = 0.5 * width;
    size_t n = 0, b = R, e = 0;
    for (size_t c = 0; c < R; ++c) {
        if (fabs(map.ring_q[c] - qpar) > half) continue;
        b = c < b ? c : b;
        e = c + 1 > e ? c + 1 : e;
        n += 1;
    }
    if (n) {
        for (size_t r = 0; r < Q; ++r) {
            double sum = 0.0;
            for (size_t c = b; c < e; ++c) sum += map.at(r, c);
            out[r] = sum / (double)n;
        }
        return n;
    }
    // No ring inside the band (a narrow band between non-uniform rings): interpolate between the neighbours
    size_t c0 = 0, c1 = 0;
    double t = 0.0;
    const bool ok = ring_interp(map, qpar, &c0, &c1, &t);
    for (size_t r = 0; r < Q; ++r) out[r] = ok ? (1.0 - t) * map.at(r, c0) + t * map.at(r, c1) : 0.0;
    return ok ? (c1 > c0 ? 2 : 1) : 0;
}

// -------------------------------------------------------------------------------------------------
// Detector
// -------------------------------------------------------------------------------------------------

bool detector_q(const Detector& det, const Beam& beam, double x_mm, double y_mm, DetectorQ* out) {
    const double k = beam.k0();
    const double ai = beam.alpha_i;
    // Direct beam direction and the vertical detector axis (sample frame, z normal to the surface)
    const double kin[3] = {cos(ai), 0.0, -sin(ai)};
    const double eup[3] = {sin(ai), 0.0, cos(ai)};
    double P[3] = {
        det.sdd_mm * kin[0] + y_mm * eup[0],
        x_mm,
        det.sdd_mm * kin[2] + y_mm * eup[2],
    };
    const double len = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
    P[0] /= len; P[1] /= len; P[2] /= len;
    const double qx = k * (P[0] - kin[0]);
    out->qy = k * P[1];
    out->qz = k * (P[2] - kin[2]);
    out->qpar = sqrt(qx * qx + out->qy * out->qy);
    out->alpha_f = asin(clampd(P[2], -1.0, 1.0));
    return P[2] >= 0.0;
}

bool detector_masked(const Detector& det, const Beam& beam, int px, int py) {
    if (det.gaps && det.module_w > 0 && det.module_h > 0) {
        const int pw = det.module_w + det.gap_w;
        const int ph = det.module_h + det.gap_h;
        if ((px % pw) >= det.module_w || (py % ph) >= det.module_h) return true;
    }
    const double dx = px + 0.5 - det.beam_x_px;
    const double dy = py + 0.5 - det.beam_y_px;
    const double bd = 0.5 * det.bs_direct_px;
    const double bs = 0.5 * det.bs_specular_px;
    const double spec_dy = tan(2.0 * beam.alpha_i) * det.sdd_mm / det.pixel_mm;   // specular spot above the direct beam
    if (fabs(dx) < bd && fabs(dy) < bd) return true;
    if (fabs(dx) < bs && fabs(dy - spec_dy) < bs) return true;
    return false;
}

void detector_bin_center_mm(const Detector& det, int bx, int by, double* x_mm, double* y_mm) {
    const int b = det.binning > 1 ? det.binning : 1;
    *x_mm = ((bx + 0.5) * b - det.beam_x_px) * det.pixel_mm;
    *y_mm = ((by + 0.5) * b - det.beam_y_px) * det.pixel_mm;
}

void detector_render(DetectorImage* img, const Detector& det, const Beam& beam, const QMap& map) {
    const int b = det.binning > 1 ? det.binning : 1;
    const int W = det.npx_h / b > 1 ? det.npx_h / b : 1;
    const int H = det.npx_v / b > 1 ? det.npx_v / b : 1;
    img->width = W;
    img->height = H;
    img->I.assign((size_t)W * H, -1.0f);

    for (int by = 0; by < H; ++by) {
        for (int bx = 0; bx < W; ++bx) {
            double sum = 0.0;
            int n = 0;
            for (int sy = 0; sy < b; ++sy) {
                for (int sx = 0; sx < b; ++sx) {
                    const int px = bx * b + sx;
                    const int py = by * b + sy;
                    if (detector_masked(det, beam, px, py)) continue;
                    DetectorQ q;
                    if (!detector_q(det, beam, (px + 0.5 - det.beam_x_px) * det.pixel_mm, (py + 0.5 - det.beam_y_px) * det.pixel_mm, &q)) continue;
                    float v;
                    if (!qmap_sample(map, q.qpar, q.qz, &v)) continue;
                    sum += v;
                    n += 1;
                }
            }
            if (n) img->I[(size_t)by * W + bx] = (float)(sum / n);
        }
    }

    // Approximate (linear) axis bounds along the central row / column, for display only
    DetectorQ q;
    const double x0 = (0.0 - det.beam_x_px) * det.pixel_mm, x1 = (W * b - det.beam_x_px) * det.pixel_mm;
    const double y0 = (0.0 - det.beam_y_px) * det.pixel_mm, y1 = (H * b - det.beam_y_px) * det.pixel_mm;
    detector_q(det, beam, x0, 0.0, &q); img->qy_min = q.qy;
    detector_q(det, beam, x1, 0.0, &q); img->qy_max = q.qy;
    detector_q(det, beam, 0.0, y0, &q); img->qz_min = q.qz;
    detector_q(det, beam, 0.0, y1, &q); img->qz_max = q.qz;
}

}  // namespace gisaxs
