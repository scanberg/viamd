#include "nci_surface_core.h"

#include "critical_point_distribution_core.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <thread>
#include <unordered_map>

namespace nci_surface {

namespace {

constexpr float s_out = 100.0f;   // what samples out of the analysis count as

struct Lattice {
    const float* rho;
    int dx, dy, dz;
    double h[3];
    size_t at(int i, int j, int k) const { return (size_t)i + (size_t)dx * ((size_t)j + (size_t)dy * (size_t)k); }
};

void derivatives(double* rho, double grad[3], double hess[6], const Lattice& l, int i, int j, int k) {
    auto r = [&](int a, int b, int c) { return (double)l.rho[l.at(i + a, j + b, k + c)]; };
    const double c = r(0, 0, 0);
    *rho = c;
    grad[0] = (r(1, 0, 0) - r(-1, 0, 0)) / (2.0 * l.h[0]);
    grad[1] = (r(0, 1, 0) - r(0, -1, 0)) / (2.0 * l.h[1]);
    grad[2] = (r(0, 0, 1) - r(0, 0, -1)) / (2.0 * l.h[2]);
    hess[0] = (r(1, 0, 0) - 2.0 * c + r(-1, 0, 0)) / (l.h[0] * l.h[0]);
    hess[3] = (r(0, 1, 0) - 2.0 * c + r(0, -1, 0)) / (l.h[1] * l.h[1]);
    hess[5] = (r(0, 0, 1) - 2.0 * c + r(0, 0, -1)) / (l.h[2] * l.h[2]);
    hess[1] = (r(1, 1, 0) - r(1, -1, 0) - r(-1, 1, 0) + r(-1, -1, 0)) / (4.0 * l.h[0] * l.h[1]);
    hess[2] = (r(1, 0, 1) - r(1, 0, -1) - r(-1, 0, 1) + r(-1, 0, -1)) / (4.0 * l.h[0] * l.h[2]);
    hess[4] = (r(0, 1, 1) - r(0, 1, -1) - r(0, -1, 1) + r(0, -1, -1)) / (4.0 * l.h[1] * l.h[2]);
}

double signed_density(const Lattice& l, int i, int j, int k) {
    i = std::clamp(i, 1, l.dx - 2);
    j = std::clamp(j, 1, l.dy - 2);
    k = std::clamp(k, 1, l.dz - 2);
    double rho, grad[3], hess[6];
    derivatives(&rho, grad, hess, l, i, j, k);
    const double l2 = critical_point_distribution::symmetric_eigenvalues(hess)[1];
    return l2 < 0.0 ? -rho : rho;
}

bool near_mask(std::vector<uint8_t>* mask, const Grid& grid, const Params& params, double offset,
               const std::atomic<int32_t>* cancel) {
    if (params.near_xyz.empty()) return true;
    mask->assign((size_t)grid.dim[0] * grid.dim[1] * grid.dim[2], 0);
    for (size_t p = 0; p + 2 < params.near_xyz.size(); p += 3) {
        if (cancel && cancel->load(std::memory_order_relaxed)) return false;
        int lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
            const double c = (params.near_xyz[p + a] - grid.origin[a]) / grid.spacing[a] - offset;
            const double r = params.near_radius / grid.spacing[a];
            lo[a] = std::max(0, (int)std::floor(c - r));
            hi[a] = std::min(grid.dim[a] - 2, (int)std::ceil(c + r));
        }
        for (int k = lo[2]; k <= hi[2]; ++k) {
            if (cancel && cancel->load(std::memory_order_relaxed)) return false;
            for (int j = lo[1]; j <= hi[1]; ++j) for (int i = lo[0]; i <= hi[0]; ++i) {
                double d2 = 0.0;
                const int ijk[3] = {i, j, k};
                for (int a = 0; a < 3; ++a) {
                    const double d = grid.origin[a] + (ijk[a] + offset) * grid.spacing[a] - params.near_xyz[p + a];
                    d2 += d * d;
                }
                if (d2 <= params.near_radius * params.near_radius) {
                    (*mask)[(size_t)i + (size_t)grid.dim[0] * ((size_t)j + (size_t)grid.dim[1] * k)] = 1;
                }
            }
        }
    }
    return true;
}

template <typename F>
void parallel_slabs(int n, F&& f) {
    const int threads = std::clamp((int)std::thread::hardware_concurrency(), 1, 32);
    std::atomic<int> next = 0;
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&]() {
            for (int k; (k = next.fetch_add(1)) < n;) f(k);
        });
    }
    for (auto& t : pool) t.join();
}

}  // namespace

void sample(double* s, double* signed_rho, const float* rho, const Grid& grid, int i, int j, int k) {
    const Lattice l = {rho, grid.dim[0], grid.dim[1], grid.dim[2], {grid.spacing[0], grid.spacing[1], grid.spacing[2]}};
    double r, grad[3], hess[6];
    derivatives(&r, grad, hess, l, i, j, k);
    *s = critical_point_distribution::reduced_density_gradient(r, grad);
    *signed_rho = critical_point_distribution::symmetric_eigenvalues(hess)[1] < 0.0 ? -r : r;
}

bool extract(Mesh* out, const float* rho, const Grid& grid, const Params& params, const std::atomic<int32_t>* cancel, std::atomic<float>* progress,
             Scatter* scatter) {
    *out = {};
    if (scatter) *scatter = {};
    const int dx = grid.dim[0], dy = grid.dim[1], dz = grid.dim[2];
    if (!rho || dx < 3 || dy < 3 || dz < 3 || !(params.s_iso > 0.0)) return false;
    for (double h : grid.spacing) if (!(h > 0.0)) return false;
    const Lattice l = {rho, dx, dy, dz, {grid.spacing[0], grid.spacing[1], grid.spacing[2]}};
    auto cancelled = [&]() { return cancel && cancel->load(std::memory_order_relaxed) != 0; };
    auto report = [&](float f) { if (progress) progress->store(f, std::memory_order_relaxed); };
    const size_t n = (size_t)dx * dy * dz;

    // s, out where the density is above rho_max, and on the outer layer
    std::vector<float> S(n, s_out);
    std::atomic<int> slabs_done = 0;
    parallel_slabs(dz - 2, [&](int kk) {
        if (cancelled()) return;
        const int k = kk + 1;
        for (int j = 1; j < dy - 1; ++j) {
            for (int i = 1; i < dx - 1; ++i) {
                const size_t v = l.at(i, j, k);
                if (!(rho[v] > 0.0f) || rho[v] > params.rho_max) continue;
                double r, grad[3], hess[6];
                derivatives(&r, grad, hess, l, i, j, k);
                const double s = critical_point_distribution::reduced_density_gradient(r, grad);
                S[v] = std::isfinite(s) ? (float)std::min(s, (double)s_out) : s_out;
            }
        }
        report(0.6f * (float)(slabs_done.fetch_add(1) + 1) / (float)(dz - 2));
    });
    if (cancelled()) return false;

    std::vector<uint8_t> allowed;
    if (scatter) {
        if (!near_mask(&allowed, grid, params, 0.0, cancel)) return false;
        std::vector<size_t> reservoir;
        reservoir.reserve(scatter_max_points);
        std::mt19937 rng(0);
        for (int k = 1; k < dz - 1; ++k) {
            if (cancelled()) return false;
            for (int j = 1; j < dy - 1; ++j) for (int i = 1; i < dx - 1; ++i) {
                const size_t v = l.at(i, j, k);
                if (!(S[v] <= scatter_s_max) || (!allowed.empty() && !allowed[v])) continue;
                ++scatter->eligible_count;
                if (reservoir.size() < scatter_max_points) {
                    reservoir.push_back(v);
                } else {
                    const size_t slot = std::uniform_int_distribution<size_t>(0, scatter->eligible_count - 1)(rng);
                    if (slot < scatter_max_points) reservoir[slot] = v;
                }
            }
        }
        scatter->signed_rho.reserve(reservoir.size());
        scatter->s.reserve(reservoir.size());
        for (size_t v : reservoir) {
            if (cancelled()) return false;
            const int i = (int)(v % dx), j = (int)((v / dx) % dy), k = (int)(v / ((size_t)dx * dy));
            const float r = (float)signed_density(l, i, j, k);
            if (!std::isfinite(r)) continue;
            scatter->signed_rho.push_back(r);
            scatter->s.push_back(S[v]);
        }
    }

    // Which cells may have a vertex
    if (!near_mask(&allowed, grid, params, 0.5, cancel)) return false;

    // Surface Nets: a vertex per cell the surface crosses, at the mean of its edges' crossings
    const float iso = (float)params.s_iso;
    std::unordered_map<size_t, uint32_t> vertex_of;
    std::unordered_map<size_t, float> signed_rho;   // per corner, as needed
    auto corner_value = [&](int i, int j, int k) {
        const size_t v = l.at(i, j, k);
        auto it = signed_rho.find(v);
        if (it != signed_rho.end()) return it->second;
        const float c = (float)signed_density(l, i, j, k);
        signed_rho.emplace(v, c);
        return c;
    };
    static constexpr int corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    static constexpr int edge[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (int k = 0; k < dz - 1; ++k) {
        if (cancelled()) return false;
        for (int j = 0; j < dy - 1; ++j) {
            for (int i = 0; i < dx - 1; ++i) {
                const size_t cell = l.at(i, j, k);
                if (!allowed.empty() && !allowed[cell]) continue;
                float s[8];
                int mask = 0;
                for (int c = 0; c < 8; ++c) {
                    s[c] = S[l.at(i + corner[c][0], j + corner[c][1], k + corner[c][2])];
                    mask |= (s[c] < iso) << c;
                }
                if (mask == 0 || mask == 0xFF) continue;
                double p[3] = {0, 0, 0};
                int crossings = 0;
                for (const auto& e : edge) {
                    if (((mask >> e[0]) & 1) == ((mask >> e[1]) & 1)) continue;
                    const double t = (iso - s[e[0]]) / (double)(s[e[1]] - s[e[0]]);
                    for (int a = 0; a < 3; ++a) p[a] += corner[e[0]][a] + t * (corner[e[1]][a] - corner[e[0]][a]);
                    ++crossings;
                }
                for (double& x : p) x /= crossings;
                // sign(lambda2) rho at the vertex, trilinear from the cell's corners
                double value = 0.0;
                for (int c = 0; c < 8; ++c) {
                    double w = 1.0;
                    for (int a = 0; a < 3; ++a) w *= corner[c][a] ? p[a] : 1.0 - p[a];
                    if (w > 0.0) value += w * corner_value(i + corner[c][0], j + corner[c][1], k + corner[c][2]);
                }
                vertex_of.emplace(cell, (uint32_t)out->value.size());
                const int ijk[3] = {i, j, k};
                for (int a = 0; a < 3; ++a) out->position.push_back((float)(grid.origin[a] + (ijk[a] + p[a]) * grid.spacing[a]));
                out->value.push_back((float)value);
            }
        }
        report(0.6f + 0.3f * (float)(k + 1) / (float)(dz - 1));
    }

    // A quad per crossed lattice edge, of the four cells around it, facing higher s
    auto vertex = [&](int i, int j, int k) -> int64_t {
        if (i < 0 || j < 0 || k < 0 || i >= dx - 1 || j >= dy - 1 || k >= dz - 1) return -1;
        auto it = vertex_of.find(l.at(i, j, k));
        return it == vertex_of.end() ? -1 : (int64_t)it->second;
    };
    for (int k = 1; k < dz - 1; ++k) {
        if (cancelled()) return false;
        for (int j = 1; j < dy - 1; ++j) {
            for (int i = 1; i < dx - 1; ++i) {
                const bool in = S[l.at(i, j, k)] < iso;
                for (int a = 0; a < 3; ++a) {
                    const int ni = i + (a == 0), nj = j + (a == 1), nk = k + (a == 2);
                    if (ni >= dx || nj >= dy || nk >= dz || (S[l.at(ni, nj, nk)] < iso) == in) continue;
                    // The cells around the edge, counter-clockwise about +axis a: in the plane of the two other
                    // axes b, c (a, b, c cyclic), cells offset by (0,0), (1,0), (1,1), (0,1) less one
                    const int b = (a + 1) % 3, c = (a + 2) % 3;
                    int64_t q[4];
                    static constexpr int ofs[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
                    bool ok = true;
                    for (int m = 0; m < 4 && ok; ++m) {
                        int at[3] = {i, j, k};
                        at[b] += ofs[m][0] - 1;
                        at[c] += ofs[m][1] - 1;
                        q[m] = vertex(at[0], at[1], at[2]);
                        ok = q[m] >= 0;
                    }
                    if (!ok) continue;
                    if (!in) std::swap(q[1], q[3]);   // s decreases along +a: the surface faces -a
                    const uint32_t t[6] = {(uint32_t)q[0], (uint32_t)q[1], (uint32_t)q[2], (uint32_t)q[0], (uint32_t)q[2], (uint32_t)q[3]};
                    out->index.insert(out->index.end(), t, t + 6);
                }
            }
        }
    }

    // Normals: of the triangles around each vertex, by area
    const size_t nv = out->value.size();
    out->normal.assign(3 * nv, 0.0f);
    for (size_t t = 0; t + 2 < out->index.size(); t += 3) {
        const float* p0 = &out->position[3 * out->index[t]];
        const float* p1 = &out->position[3 * out->index[t + 1]];
        const float* p2 = &out->position[3 * out->index[t + 2]];
        const float u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const float w[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        const float nrm[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
        for (int m = 0; m < 3; ++m) for (int a = 0; a < 3; ++a) out->normal[3 * out->index[t + m] + a] += nrm[a];
    }
    for (size_t v = 0; v < nv; ++v) {
        float* nrm = &out->normal[3 * v];
        const float len = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
        if (len > 0.0f) for (int a = 0; a < 3; ++a) nrm[a] /= len;
    }
    report(1.0f);
    return true;
}

}  // namespace nci_surface
