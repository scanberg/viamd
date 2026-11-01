#include "nci_surface_core.h"

#include "critical_point_distribution_core.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <thread>
#include <unordered_map>

namespace nci_surface {

static Family contact_family(const std::string& interaction) {
    if (interaction == "Hydrogen bond") return Family::HydrogenBond;
    if (interaction == "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80") return Family::PiStacking;
    if (interaction == "Halogen bond") return Family::HalogenBond;
    return Family::OtherContact;
}

Family family_of(const critical_point_distribution::ChemicalLabel& label) {
    using critical_point_distribution::Kind;
    switch (label.kind) {
    case Kind::Contact: return contact_family(label.interaction);
    case Kind::Ring: return Family::Ring;
    case Kind::Cage: return Family::Cage;
    case Kind::Bond: return Family::Bond;
    case Kind::BondPath: return Family::BondPath;
    case Kind::Atom: case Kind::NonNuclear: return Family::Atom;
    case Kind::Unassigned: return Family::Unassigned;
    }
    return Family::Unassigned;
}

const char* family_name(Family family) {
    switch (family) {
    case Family::HydrogenBond: return "Hydrogen bonds";
    case Family::Ring: return "Rings";
    case Family::PiStacking: return "\xCF\x80-stacking";
    case Family::HalogenBond: return "Halogen bonds";
    case Family::OtherContact: return "Other contacts";
    case Family::Bond: return "Bonds";
    case Family::BondPath: return "Bond paths (connectivity unknown)";
    case Family::Atom: return "Attractors";
    case Family::Cage: return "Cages";
    case Family::Unassigned: return "Unassigned";
    case Family::Count: break;
    }
    return "Unassigned";
}

Family parameter_family_of(const critical_point_distribution::ChemicalLabel& label) {
    using critical_point_distribution::Kind;
    // A ring or cage closed by contacts lies in the region of its closing interaction: it takes that interaction's
    // parameters, so the region is cut at one threshold
    if ((label.kind == Kind::Ring || label.kind == Kind::Cage) && !label.interaction.empty())
        return contact_family(label.interaction);
    return family_of(label);
}

bool has_nci(Family family) {
    return family != Family::Bond && family != Family::BondPath && family != Family::Atom;
}

bool single_family(Families families) {
    return families && !(families & (families - 1));
}

Defaults family_defaults(Families families) {
    if (families == family_bit(Family::HydrogenBond)) return {1.0, 0.07, 1.7, true};
    if (families == family_bit(Family::Ring)) return {1.0, 0.05, 2.0, true};
    if (families == family_bit(Family::PiStacking)) return {0.85, 0.07, 3.7, true};
    // Covalent cages: the 3D counterpart of ring centres
    if (families == family_bit(Family::Cage)) return {1.0, 0.05, 2.0, true};
    // Provisional, not calibrated: densities like medium hydrogen bonds over a longer contact
    if (families == family_bit(Family::HalogenBond)) return {0.9, 0.07, 2.2, false};
    return {};
}

UniformDefaults uniform_defaults() { return {}; }

ExpectedRegion expected_region(Family family) {
    switch (family) {
    case Family::HydrogenBond: return {-0.050f, -0.008f, 1.0f, true};
    case Family::HalogenBond:  return {-0.030f, -0.005f, 0.9f, true};
    case Family::PiStacking:   return {-0.015f, 0.008f, 0.85f, true};   // the contacts and the rings they close
    case Family::Ring:         return {0.004f, 0.060f, 1.0f, true};
    case Family::Cage:         return {0.002f, 0.030f, 1.0f, true};
    case Family::OtherContact: return {-0.012f, 0.010f, 0.5f, true};
    default:                   return {};
    }
}

bool in_expected_region(const ExpectedRegion& region, float signed_rho, float s) {
    return region.valid && s <= region.s_max && signed_rho >= region.rho_min && signed_rho <= region.rho_max;
}

RegionCounts expected_region_counts(const Scatter& scatter, Family family) {
    RegionCounts counts;
    const ExpectedRegion region = expected_region(family);
    if (!region.valid) return counts;
    for (size_t i = 0; i < scatter.s.size() && i < scatter.signed_rho.size(); ++i) {
        if (scatter.s[i] > region.s_max) continue;
        if (in_expected_region(region, scatter.signed_rho[i], scatter.s[i])) ++counts.inside;
        else ++counts.outside;
    }
    return counts;
}

CriticalPointGroups group_critical_points(std::span<const CriticalPointAnchor> points) {
    CriticalPointGroups groups;
    for (const CriticalPointAnchor& point : points) {
        const size_t family = (size_t)point.family;
        auto& all = groups.all_anchors[family];
        all.insert(all.end(), point.position.begin(), point.position.end());
        if (point.selected) {
            auto& selected = groups.selected_anchors[family];
            selected.insert(selected.end(), point.position.begin(), point.position.end());
        }
    }
    return groups;
}

namespace {

constexpr int outline_bins = 48, outline_dim = outline_bins + 2;
using OutlineField = std::array<float, outline_dim * outline_dim>;

// The smoothed (signed rho, s) histogram, on bins 1..outline_bins with an empty border, and its peak
bool smooth_histogram(OutlineField* smooth, float* peak, const Scatter& scatter, float rho_max) {
    constexpr int bins = outline_bins, dim = outline_dim;
    OutlineField count = {};
    smooth->fill(0.0f);
    *peak = 0.0f;
    for (size_t i = 0; i < scatter.s.size(); ++i) {
        const float r = scatter.signed_rho[i], s = scatter.s[i];
        if (!std::isfinite(r) || !std::isfinite(s)) return false;
        if (std::abs(r) > rho_max || s < 0 || s > scatter_s_max) continue;
        const int x = 1 + std::clamp((int)((r + rho_max) / (2 * rho_max) * bins), 0, bins - 1);
        const int y = 1 + std::clamp((int)(s / scatter_s_max * bins), 0, bins - 1);
        count[x + dim * y] += 1.0f;
    }
    for (int y = 1; y <= bins; ++y) for (int x = 1; x <= bins; ++x) {
        float sum = 0.0f;
        for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a)
            sum += count[x + a + dim * (y + b)] * (a == 0 ? 2.0f : 1.0f) * (b == 0 ? 2.0f : 1.0f);
        (*smooth)[x + dim * y] = sum;
        *peak = std::max(*peak, sum);
    }
    return true;
}

void contour(std::vector<std::array<float, 4>>* out, const OutlineField& field, float level, float rho_max) {
    constexpr int dim = outline_dim;
    constexpr int corner[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (int y = 0; y < dim - 1; ++y) for (int x = 0; x < dim - 1; ++x) {
        float value[4];
        for (int c = 0; c < 4; ++c) value[c] = field[x + corner[c][0] + dim * (y + corner[c][1])];
        std::array<std::array<float, 2>, 4> crossings;
        int num_crossings = 0;
        for (int c = 0; c < 4; ++c) {
            const int next = (c + 1) % 4;
            if ((value[c] >= level) == (value[next] >= level)) continue;
            const float t = (level - value[c]) / (value[next] - value[c]);
            const float px = x + corner[c][0] + t * (corner[next][0] - corner[c][0]);
            const float py = y + corner[c][1] + t * (corner[next][1] - corner[c][1]);
            crossings[num_crossings++] = {-rho_max + 2 * rho_max * px / (dim - 1), scatter_s_max * py / (dim - 1)};
        }
        const auto segment = [&](int a, int b) {
            out->push_back({crossings[a][0], crossings[a][1], crossings[b][0], crossings[b][1]});
        };
        if (num_crossings == 2) segment(0, 1);
        else if (num_crossings == 4) {
            // Resolve saddle cells by the bilinear value at their centre.
            const bool centre_inside = 0.25f * (value[0] + value[1] + value[2] + value[3]) >= level;
            if ((value[0] >= level) == centre_inside) {
                segment(0, 1);
                segment(2, 3);
            } else {
                segment(3, 0);
                segment(1, 2);
            }
        }
    }
}

bool valid_outline_input(const Scatter& scatter, float rho_max) {
    return rho_max > 0.0f && std::isfinite(rho_max) && scatter.s.size() == scatter.signed_rho.size();
}

}  // namespace

bool scatter_density_outline(std::vector<std::array<float, 4>>* out, const Scatter& scatter, float rho_max) {
    out->clear();
    if (!valid_outline_input(scatter, rho_max)) return false;
    OutlineField smooth;
    float peak = 0.0f;
    if (!smooth_histogram(&smooth, &peak, scatter, rho_max)) return false;
    if (peak > 0) contour(out, smooth, 0.1f * peak, rho_max);
    return true;
}

namespace {

// The union of the scatters' 10%-of-peak regions as a field whose level 1 is the boundary, grown by whole bins.
// Returns false on invalid input; *any tells whether a scatter had samples.
bool union_field(OutlineField* field, bool* any, std::span<const Scatter* const> scatters, float rho_max, int grow_bins) {
    constexpr int bins = outline_bins, dim = outline_dim;
    field->fill(0.0f);
    *any = false;
    for (const Scatter* scatter : scatters) {
        if (!scatter || !valid_outline_input(*scatter, rho_max)) return false;
        OutlineField smooth;
        float peak = 0.0f;
        if (!smooth_histogram(&smooth, &peak, *scatter, rho_max)) return false;
        if (!(peak > 0)) continue;
        *any = true;
        for (size_t i = 0; i < field->size(); ++i) (*field)[i] = std::max((*field)[i], smooth[i] / (0.1f * peak));
    }
    // Grown by whole bins, inside the empty border so the contour stays closed
    for (int g = 0; g < grow_bins; ++g) {
        OutlineField grown = *field;
        for (int y = 1; y <= bins; ++y) for (int x = 1; x <= bins; ++x)
            for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a)
                grown[x + dim * y] = std::max(grown[x + dim * y], (*field)[x + a + dim * (y + b)]);
        *field = grown;
    }
    return true;
}

}  // namespace

bool scatter_union_outline(std::vector<std::array<float, 4>>* out, std::span<const Scatter* const> scatters, float rho_max,
                           int grow_bins) {
    out->clear();
    OutlineField field;
    bool any = false;
    if (!union_field(&field, &any, scatters, rho_max, grow_bins)) return false;
    if (any) contour(out, field, 1.0f, rho_max);
    return true;
}

bool scatter_union_fill(std::vector<std::array<float, 2>>* triangles, std::span<const Scatter* const> scatters,
                        float rho_max, int grow_bins) {
    triangles->clear();
    OutlineField field;
    bool any = false;
    if (!union_field(&field, &any, scatters, rho_max, grow_bins)) return false;
    if (!any) return true;
    constexpr int dim = outline_dim;
    constexpr int corner[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    const auto point = [&](float px, float py) -> std::array<float, 2> {
        return {-rho_max + 2 * rho_max * px / (dim - 1), scatter_s_max * py / (dim - 1)};
    };
    // Per cell, the part inside the boundary: its inside corners and the crossings, around the cell, as a fan
    for (int y = 0; y < dim - 1; ++y) for (int x = 0; x < dim - 1; ++x) {
        float value[4];
        for (int c = 0; c < 4; ++c) value[c] = field[x + corner[c][0] + dim * (y + corner[c][1])];
        std::array<std::array<float, 2>, 8> poly;
        int n = 0;
        for (int c = 0; c < 4; ++c) {
            const int next = (c + 1) % 4;
            if (value[c] >= 1.0f) poly[n++] = point((float)(x + corner[c][0]), (float)(y + corner[c][1]));
            if ((value[c] >= 1.0f) == (value[next] >= 1.0f)) continue;
            const float t = (1.0f - value[c]) / (value[next] - value[c]);
            poly[n++] = point(x + corner[c][0] + t * (corner[next][0] - corner[c][0]),
                              y + corner[c][1] + t * (corner[next][1] - corner[c][1]));
        }
        for (int k = 1; k + 1 < n; ++k) {
            triangles->push_back(poly[0]);
            triangles->push_back(poly[k]);
            triangles->push_back(poly[k + 1]);
        }
    }
    return true;
}

bool scatter_heatmap(Heatmap* out, std::span<const Scatter* const> scatters, float rho_max) {
    constexpr int bins = outline_bins, dim = outline_dim;
    out->bins = bins;
    out->rho_max = rho_max;
    out->values.assign((size_t)bins * bins, 0.0f);
    OutlineField sum = {};
    for (const Scatter* scatter : scatters) {
        if (!scatter || !valid_outline_input(*scatter, rho_max)) return false;
        OutlineField smooth;
        float peak = 0.0f;
        if (!smooth_histogram(&smooth, &peak, *scatter, rho_max)) return false;
        for (size_t i = 0; i < sum.size(); ++i) sum[i] += smooth[i];
    }
    float peak = 0.0f;
    for (float v : sum) peak = std::max(peak, v);
    if (!(peak > 0)) return true;
    for (int y = 0; y < bins; ++y) for (int x = 0; x < bins; ++x)
        out->values[x + (size_t)bins * y] = sum[x + 1 + dim * (y + 1)] / peak;
    return true;
}

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
    mask->clear();
    if (params.near_xyz.empty() || !params.radius_cap) return true;
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
    if (!rho || dx < 3 || dy < 3 || dz < 3 || !(params.s_iso > 0.0) || !std::isfinite(params.s_iso) ||
        !(params.rho_max > 0.0) || !std::isfinite(params.rho_max) || params.near_xyz.size() % 3 != 0 ||
        (params.radius_cap && (!(params.near_radius > 0.0) || !std::isfinite(params.near_radius)))) return false;
    for (double h : grid.spacing) if (!(h > 0.0) || !std::isfinite(h)) return false;
    for (double x : grid.origin) if (!std::isfinite(x)) return false;
    for (double x : params.near_xyz) if (!std::isfinite(x)) return false;
    const Lattice l = {rho, dx, dy, dz, {grid.spacing[0], grid.spacing[1], grid.spacing[2]}};
    auto cancelled = [&]() { return cancel && cancel->load(std::memory_order_relaxed) != 0; };
    auto report = [&](float f) { if (progress) progress->store(f, std::memory_order_relaxed); };
    size_t n = 1;
    for (int d : grid.dim) {
        if (n > std::numeric_limits<uint32_t>::max() / (size_t)d) return false;
        n *= (size_t)d;
    }

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
    if (!near_mask(&allowed, grid, params, 0.0, cancel)) return false;
    std::vector<uint8_t> selected;
    std::vector<uint8_t> halo;
    if (params.connected && !params.near_xyz.empty()) {
        selected.assign(n, 0);
        std::vector<uint32_t> queue;
        // A CP rarely lies exactly on the lattice. Seed the nearest eligible sample within two voxels.
        for (size_t p = 0; p < params.near_xyz.size(); p += 3) {
            if (cancelled()) return false;
            int centre[3];
            bool inside = true;
            for (int a = 0; a < 3; ++a) {
                const double c = (params.near_xyz[p + a] - grid.origin[a]) / grid.spacing[a];
                inside &= c >= -2.0 && c <= grid.dim[a] + 1.0;
                centre[a] = (int)std::lround(std::clamp(c, -2.0, (double)grid.dim[a] + 1.0));
            }
            if (!inside) continue;
            size_t seed = n;
            double closest = INFINITY;
            for (int k = std::max(1, centre[2] - 2); k <= std::min(dz - 2, centre[2] + 2); ++k)
                for (int j = std::max(1, centre[1] - 2); j <= std::min(dy - 2, centre[1] + 2); ++j)
                    for (int i = std::max(1, centre[0] - 2); i <= std::min(dx - 2, centre[0] + 2); ++i) {
                        const size_t v = l.at(i, j, k);
                        if (!(S[v] < params.s_iso) || (!allowed.empty() && !allowed[v])) continue;
                        const int ijk[3] = {i, j, k};
                        double d2 = 0.0;
                        for (int a = 0; a < 3; ++a) {
                            const double d = grid.origin[a] + ijk[a] * grid.spacing[a] - params.near_xyz[p + a];
                            d2 += d * d;
                        }
                        if (d2 < closest) { closest = d2; seed = v; }
                    }
            if (seed != n && !selected[seed]) { selected[seed] = 1; queue.push_back((uint32_t)seed); }
        }
        for (size_t head = 0; head < queue.size(); ++head) {
            if ((head & 4095) == 0 && cancelled()) return false;
            const size_t v = queue[head];
            const int ijk[3] = {(int)(v % dx), (int)((v / dx) % dy), (int)(v / ((size_t)dx * dy))};
            for (int a = 0; a < 3; ++a) for (int d : {-1, 1}) {
                int q[3] = {ijk[0], ijk[1], ijk[2]};
                q[a] += d;
                if (q[a] < 1 || q[a] >= grid.dim[a] - 1) continue;
                const size_t u = l.at(q[0], q[1], q[2]);
                if (!selected[u] && S[u] < params.s_iso && (allowed.empty() || allowed[u])) {
                    selected[u] = 1;
                    queue.push_back((uint32_t)u);
                }
            }
        }
        halo = selected;
        for (size_t head = 0; head < queue.size(); ++head) {
            if ((head & 4095) == 0 && cancelled()) return false;
            const size_t v = queue[head];
            const int i = (int)(v % dx), j = (int)((v / dx) % dy), k = (int)(v / ((size_t)dx * dy));
            for (int c = -1; c <= 1; ++c) for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a)
                halo[l.at(i + a, j + b, k + c)] = 1;
        }
    }
    if (scatter) {
        std::vector<size_t> reservoir;
        reservoir.reserve(scatter_max_points);
        std::mt19937 rng(0);
        for (int k = 1; k < dz - 1; ++k) {
            if (cancelled()) return false;
            for (int j = 1; j < dy - 1; ++j) for (int i = 1; i < dx - 1; ++i) {
                const size_t v = l.at(i, j, k);
                if (!(S[v] <= scatter_s_max) || (!allowed.empty() && !allowed[v]) || (!halo.empty() && !halo[v])) continue;
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
        scatter->voxel.reserve(reservoir.size());
        for (size_t v : reservoir) {
            if (cancelled()) return false;
            const int i = (int)(v % dx), j = (int)((v / dx) % dy), k = (int)(v / ((size_t)dx * dy));
            const float r = (float)signed_density(l, i, j, k);
            if (!std::isfinite(r)) continue;
            scatter->signed_rho.push_back(r);
            scatter->s.push_back(S[v]);
            scatter->voxel.push_back((uint32_t)v);
        }
    }

    if (!params.surface) { report(1.0f); return true; }
    if (!selected.empty()) {
        for (size_t v = 0; v < n; ++v) {
            if ((v & 4095) == 0 && cancelled()) return false;
            if (S[v] < params.s_iso && !selected[v]) S[v] = s_out;
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
