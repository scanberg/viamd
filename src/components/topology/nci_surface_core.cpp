#include "nci_surface_core.h"

#include "critical_point_distribution_core.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <mutex>
#include <random>
#include <thread>
#include <unordered_map>

namespace nci_surface {

static constexpr ContactClass contact_class_table[] = {
    {"Hydrogen bond", Family::HydrogenBond},
    {"Ionic hydrogen bond", Family::HydrogenBond},
    {"Salt bridge", Family::HydrogenBond},
    {"\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"X-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"Cation\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"Anion\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"Lone pair\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", Family::PiStacking},
    {"\xCF\x80-hole bond", Family::PiStacking},
    // Sigma-hole bonds share one type
    {"Halogen bond", Family::HalogenBond},
    {"Chalcogen bond", Family::HalogenBond},
    {"Pnictogen bond", Family::HalogenBond},
    {"Tetrel bond", Family::HalogenBond},
    {"Triel bond", Family::HalogenBond},
    {"Aerogen bond", Family::HalogenBond},
    {"Halogen\xC2\xB7\xC2\xB7\xC2\xB7halogen", Family::OtherContact},
    {"Dihydrogen bond", Family::OtherContact},
    {"Hydride bond", Family::OtherContact},
    {"Agostic", Family::OtherContact},
    {"Coordination", Family::OtherContact},
    {"Ion pair", Family::OtherContact},
    {"Metal\xC2\xB7\xC2\xB7\xC2\xB7metal", Family::OtherContact},
    {"Metal contact", Family::OtherContact},
    {"C-H\xC2\xB7\xC2\xB7\xC2\xB7X", Family::OtherContact},
    {"H\xC2\xB7\xC2\xB7\xC2\xB7H", Family::OtherContact},
    {"van der Waals contact", Family::OtherContact},
};

std::span<const ContactClass> contact_classes() { return contact_class_table; }

Family contact_family(std::string_view interaction) {
    for (const ContactClass& c : contact_class_table)
        if (interaction == c.name) return c.family;
    return Family::OtherContact;
}

size_t reference_category(Family family) {
    switch (family) {
    case Family::HydrogenBond: return 0;
    case Family::Ring: return 1;
    case Family::PiStacking: return 2;
    case Family::HalogenBond: case Family::OtherContact: case Family::Cage: return 3;
    default: return reference_category_count;
    }
}

const char* reference_category_name(size_t category) {
    static constexpr const char* names[reference_category_count] = {"Hydrogen bonds", "Rings", "\xCF\x80 interactions",
                                                                     "Other"};
    return category < reference_category_count ? names[category] : "Unassigned";
}

Family reference_settings_family(size_t category) {
    static constexpr Family families[reference_category_count] = {Family::HydrogenBond, Family::Ring, Family::PiStacking,
                                                                  Family::OtherContact};
    return category < reference_category_count ? families[category] : Family::Unassigned;
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
    case Family::PiStacking: return "\xCF\x80 interactions";
    case Family::HalogenBond: return "\xCF\x83-hole bonds";
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
    const auto preset = [](double radius, bool established) {
        Defaults d;
        d.radius = radius;
        d.established = established;
        return d;
    };
    if (families == family_bit(Family::HydrogenBond)) return preset(1.7, true);
    if (families == family_bit(Family::Ring)) return preset(2.0, true);
    if (families == family_bit(Family::PiStacking)) return preset(3.7, true);
    // Covalent cages: the 3D counterpart of ring centres
    if (families == family_bit(Family::Cage)) return preset(2.0, true);
    // Provisional, not calibrated: a longer contact than a hydrogen bond
    if (families == family_bit(Family::HalogenBond)) return preset(2.2, false);
    return {};
}

UniformDefaults uniform_defaults() { return {}; }

Params reference_params(std::span<const double> near_xyz, std::span<const double> other_xyz) {
    Params p;
    p.s_iso = reference_scatter_s;
    p.scatter_s = reference_scatter_s;
    p.rho_max = reference_rho_max;
    p.near_xyz = near_xyz;
    p.other_xyz = other_xyz;
    p.near_radius = reference_radius;
    p.radius_cap = true;
    p.connected = true;
    p.surface = false;
    p.by_anchor = true;
    return p;
}

Scatter scatter_subset(const Scatter& scatter, std::span<const uint8_t> keep) {
    const auto kept = [&](uint32_t anchor) { return anchor < keep.size() && keep[anchor]; };
    Scatter out;
    for (size_t i = 0; i < scatter.anchor.size() && i < scatter.signed_rho.size(); ++i) {
        if (!kept(scatter.anchor[i])) continue;
        out.signed_rho.push_back(scatter.signed_rho[i]);
        out.s.push_back(scatter.s[i]);
        if (i < scatter.voxel.size()) out.voxel.push_back(scatter.voxel[i]);
        out.anchor.push_back(scatter.anchor[i]);
    }
    if (!scatter.signed_rho.empty())
        out.eligible_count = (size_t)std::llround((double)scatter.eligible_count * out.signed_rho.size() / scatter.signed_rho.size());
    const ScatterDensity& d = scatter.density;
    if (!d.volume.empty()) {
        out.density.rho_extent = d.rho_extent;
        out.density.volume.assign(d.volume.size(), 0.0f);
        double total = 0.0;
        for (const AnchorBin& b : d.by_anchor) {
            if (!kept(b.anchor) || b.bin >= out.density.volume.size()) continue;
            out.density.volume[b.bin] += b.volume;
            out.density.by_anchor.push_back(b);
            total += b.volume;
        }
        out.density.total_volume = total;
        const double share = d.total_volume > 0.0 ? total / d.total_volume : 0.0;
        out.density.cells = (size_t)std::llround(d.cells * share);
        out.density.subsamples = (size_t)std::llround(d.subsamples * share);
    }
    return out;
}

void scatter_append(Scatter* dst, const Scatter& src) {
    dst->signed_rho.insert(dst->signed_rho.end(), src.signed_rho.begin(), src.signed_rho.end());
    dst->s.insert(dst->s.end(), src.s.begin(), src.s.end());
    dst->voxel.insert(dst->voxel.end(), src.voxel.begin(), src.voxel.end());
    dst->anchor.insert(dst->anchor.end(), src.anchor.begin(), src.anchor.end());
    dst->eligible_count += src.eligible_count;
    const ScatterDensity& s = src.density;
    if (s.volume.empty()) return;
    ScatterDensity& d = dst->density;
    if (d.volume.empty()) {
        d.rho_extent = s.rho_extent;
        d.volume.assign(s.volume.size(), 0.0f);
    }
    if (d.volume.size() != s.volume.size()) return;
    for (size_t b = 0; b < s.volume.size(); ++b) d.volume[b] += s.volume[b];
    d.by_anchor.insert(d.by_anchor.end(), s.by_anchor.begin(), s.by_anchor.end());
    d.total_volume += s.total_volume;
    d.cells += s.cells;
    d.subsamples += s.subsamples;
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

constexpr int outline_bins = scatter_field_bins, outline_dim = outline_bins + 2;
using OutlineField = std::array<float, outline_dim * outline_dim>;
static_assert(scatter_density_bins % scatter_field_bins == 0);

void smooth_once(OutlineField* field) {
    constexpr int bins = outline_bins, dim = outline_dim;
    const OutlineField in = *field;
    for (int y = 1; y <= bins; ++y) for (int x = 1; x <= bins; ++x) {
        float sum = 0.0f;
        for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a)
            sum += in[x + a + dim * (y + b)] * (a == 0 ? 2.0f : 1.0f) * (b == 0 ? 2.0f : 1.0f);
        (*field)[x + dim * y] = sum / 16.0f;
    }
}

// The continuous density on the field's bins over signed rho [-rho_max, rho_max]: its bins split by overlap in
// signed rho; in s the bins nest.
void resample_density(OutlineField* field, const ScatterDensity& density, float rho_max) {
    constexpr int fine = scatter_density_bins, ratio = scatter_density_bins / outline_bins;
    const float fine_w = 2.0f * density.rho_extent / fine, field_w = 2.0f * rho_max / outline_bins;
    for (int fx = 0; fx < fine; ++fx) {
        const float lo = -density.rho_extent + fx * fine_w, hi = lo + fine_w;
        const float flo = (lo + rho_max) / field_w, fhi = (hi + rho_max) / field_w;
        if (fhi <= 0.0f || flo >= outline_bins) continue;
        for (int x = std::max(0, (int)std::floor(flo)); x < std::min(outline_bins, (int)std::ceil(fhi)); ++x) {
            const float overlap = (std::min(fhi, (float)(x + 1)) - std::max(flo, (float)x)) / (fhi - flo);
            if (!(overlap > 0.0f)) continue;
            for (int fy = 0; fy < fine; ++fy) {
                const float v = density.volume[fx + (size_t)fine * fy];
                if (v > 0.0f) (*field)[(x + 1) + outline_dim * (fy / ratio + 1)] += overlap * v;
            }
        }
    }
}

bool has_density(const Scatter& scatter) {
    return scatter.density.rho_extent > 0.0f &&
           scatter.density.volume.size() == (size_t)scatter_density_bins * scatter_density_bins;
}

// The smoothed (signed rho, s) density, on bins 1..outline_bins with an empty border, and its peak: the continuous
// density when there is one, smoothed lightly; the points' histogram otherwise, smoothed twice as its noise needs.
bool smooth_histogram(OutlineField* smooth, float* peak, const Scatter& scatter, float rho_max) {
    constexpr int bins = outline_bins, dim = outline_dim;
    smooth->fill(0.0f);
    *peak = 0.0f;
    for (size_t i = 0; i < scatter.s.size(); ++i)
        if (!std::isfinite(scatter.signed_rho[i]) || !std::isfinite(scatter.s[i])) return false;
    if (has_density(scatter)) {
        resample_density(smooth, scatter.density, rho_max);
        smooth_once(smooth);
    } else {
        for (size_t i = 0; i < scatter.s.size(); ++i) {
            const float r = scatter.signed_rho[i], s = scatter.s[i];
            if (std::abs(r) > rho_max || s < 0 || s > scatter_s_max) continue;
            const int x = 1 + std::clamp((int)((r + rho_max) / (2 * rho_max) * bins), 0, bins - 1);
            const int y = 1 + std::clamp((int)(s / scatter_s_max * bins), 0, bins - 1);
            (*smooth)[x + dim * y] += 1.0f;
        }
        smooth_once(smooth);
        smooth_once(smooth);
    }
    for (int y = 1; y <= bins; ++y) for (int x = 1; x <= bins; ++x) *peak = std::max(*peak, (*smooth)[x + dim * y]);
    return true;
}

// The highest-density level whose superlevel set holds the given fraction of the field's mass: the region edge.
float mass_level(const OutlineField& field, float fraction) {
    constexpr int bins = outline_bins, dim = outline_dim;
    std::vector<float> values;
    values.reserve((size_t)bins * bins);
    double total = 0.0;
    for (int y = 1; y <= bins; ++y) for (int x = 1; x <= bins; ++x) {
        const float v = field[x + dim * y];
        if (v > 0.0f) {
            values.push_back(v);
            total += v;
        }
    }
    if (!(total > 0.0)) return 0.0f;
    std::sort(values.begin(), values.end(), std::greater<float>());
    double sum = 0.0;
    for (float v : values) {
        sum += v;
        if (sum >= fraction * total) return v;
    }
    return values.back();
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
    if (peak > 0) contour(out, smooth, mass_level(smooth, region_mass_fraction), rho_max);
    return true;
}

namespace {

// The union of the scatters' regions holding mass_fraction of their volume, as a field whose level 1 is the boundary,
// grown by whole bins. Returns false on invalid input; *any tells whether a scatter had samples.
bool union_field(OutlineField* field, bool* any, std::span<const Scatter* const> scatters, float rho_max, int grow_bins,
                 float mass_fraction = region_mass_fraction) {
    constexpr int bins = outline_bins, dim = outline_dim;
    field->fill(0.0f);
    *any = false;
    for (const Scatter* scatter : scatters) {
        if (!scatter || !valid_outline_input(*scatter, rho_max)) return false;
        OutlineField smooth;
        float peak = 0.0f;
        if (!smooth_histogram(&smooth, &peak, *scatter, rho_max)) return false;
        if (!(peak > 0)) continue;
        const float level = mass_level(smooth, mass_fraction);
        if (!(level > 0)) continue;
        *any = true;
        for (size_t i = 0; i < field->size(); ++i) (*field)[i] = std::max((*field)[i], smooth[i] / level);
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
                           int grow_bins, float mass_fraction) {
    out->clear();
    OutlineField field;
    bool any = false;
    if (!union_field(&field, &any, scatters, rho_max, grow_bins, mass_fraction)) return false;
    if (any) contour(out, field, 1.0f, rho_max);
    return true;
}

ExpectedRegion expected_region(std::span<const Scatter* const> scatters, float rho_extent, int grow_bins) {
    ExpectedRegion region;
    if (!(rho_extent > 0.0f) || !std::isfinite(rho_extent) || grow_bins < 0) return region;

    OutlineField field;
    bool any = false;
    if (!union_field(&field, &any, scatters, rho_extent, grow_bins) || !any) return region;

    region.rho_extent = rho_extent;
    region.mask.resize((size_t)outline_bins * outline_bins);
    for (int y = 0; y < outline_bins; ++y) for (int x = 0; x < outline_bins; ++x)
        region.mask[x + (size_t)outline_bins * y] = field[(x + 1) + outline_dim * (y + 1)] >= 1.0f;
    contour(&region.outline, field, 1.0f, rho_extent);
    OutlineField core;
    if (union_field(&core, &any, scatters, rho_extent, 0, region_core_mass_fraction) && any)
        contour(&region.core_outline, core, 1.0f, rho_extent);
    region.density.assign((size_t)outline_bins * outline_bins, 0.0f);
    OutlineField density = {};
    for (const Scatter* scatter : scatters) {
        OutlineField smooth;
        float peak = 0.0f;
        if (!smooth_histogram(&smooth, &peak, *scatter, rho_extent) || !(peak > 0.0f)) continue;
        for (size_t i = 0; i < density.size(); ++i) density[i] += smooth[i] / peak;
    }
    float density_peak = 0.0f;
    for (int y = 1; y <= outline_bins; ++y) for (int x = 1; x <= outline_bins; ++x)
        density_peak = std::max(density_peak, density[x + outline_dim * y]);
    if (density_peak > 0.0f) {
        for (int y = 0; y < outline_bins; ++y) for (int x = 0; x < outline_bins; ++x)
            region.density[x + (size_t)outline_bins * y] =
                density[(x + 1) + outline_dim * (y + 1)] / density_peak;
    }
    region.valid = !region.outline.empty();
    return region;
}

void expected_region_fill(std::vector<std::array<float, 2>>* triangles, const ExpectedRegion& region) {
    constexpr int bins = outline_bins;
    triangles->clear();
    if (!region.valid || !(region.rho_extent > 0.0f) || region.mask.size() != (size_t)bins * bins) return;
    for (int y = 0; y < bins; ++y) for (int x = 0; x < bins; ++x) {
        if (!region.mask[x + (size_t)bins * y]) continue;
        const float left = -region.rho_extent + 2.0f * region.rho_extent * x / bins;
        const float right = -region.rho_extent + 2.0f * region.rho_extent * (x + 1) / bins;
        const float bottom = scatter_s_max * y / bins;
        const float top = scatter_s_max * (y + 1) / bins;
        triangles->push_back({left, bottom});
        triangles->push_back({right, bottom});
        triangles->push_back({right, top});
        triangles->push_back({left, bottom});
        triangles->push_back({right, top});
        triangles->push_back({left, top});
    }
}

bool in_expected_region(const ExpectedRegion& region, float signed_rho, float s) {
    if (!region.valid || !(region.rho_extent > 0.0f) || !std::isfinite(signed_rho) || !std::isfinite(s) ||
        signed_rho < -region.rho_extent || signed_rho > region.rho_extent || s < 0.0f || s > scatter_s_max ||
        region.mask.size() != (size_t)outline_bins * outline_bins) return false;
    const int x = std::clamp((int)((signed_rho + region.rho_extent) / (2.0f * region.rho_extent) * outline_bins), 0, outline_bins - 1);
    const int y = std::clamp((int)(s / scatter_s_max * outline_bins), 0, outline_bins - 1);
    return region.mask[x + (size_t)outline_bins * y] != 0;
}

RegionCounts expected_region_counts(const Scatter& scatter, const ExpectedRegion& region) {
    RegionCounts counts;
    if (!region.valid) return counts;
    for (size_t i = 0; i < scatter.s.size() && i < scatter.signed_rho.size(); ++i) {
        if (!std::isfinite(scatter.s[i]) || !std::isfinite(scatter.signed_rho[i]) ||
            scatter.s[i] < 0.0f || scatter.s[i] > scatter_s_max) continue;
        if (in_expected_region(region, scatter.signed_rho[i], scatter.s[i])) ++counts.inside;
        else ++counts.outside;
    }
    return counts;
}

namespace {

// The bordered field bilinearly at a point, on the contour's lattice: node 0 at the low edge, node outline_dim - 1
// at the high one, so the fill and the isolines agree. Negative outside the plot.
float reference_value(const std::vector<float>& field, float rho_extent, float signed_rho, float s) {
    constexpr int dim = outline_dim;
    if (field.size() != (size_t)dim * dim || !std::isfinite(signed_rho) || !std::isfinite(s) ||
        signed_rho < -rho_extent || signed_rho > rho_extent || s < 0.0f || s > scatter_s_max) return -1.0f;
    const float px = (signed_rho + rho_extent) / (2.0f * rho_extent) * (dim - 1);
    const float py = s / scatter_s_max * (dim - 1);
    const int x = std::min((int)px, dim - 2), y = std::min((int)py, dim - 2);
    const float tx = px - x, ty = py - y;
    const auto at = [&](int i, int j) { return field[(size_t)(x + i) + (size_t)dim * (y + j)]; };
    return (1 - ty) * ((1 - tx) * at(0, 0) + tx * at(1, 0)) + ty * ((1 - tx) * at(0, 1) + tx * at(1, 1));
}

int level_of(const ReferenceCategory& c, float value) {
    if (!c.valid || !(value > 0.0f)) return 0;
    int level = 0;
    for (float t : c.thresholds) level += value >= t;
    return level;
}

}  // namespace

bool build_reference_map(ReferenceMap* out, const ReferenceData& data, std::span<const float> fractions) {
    *out = {};
    if (data.bins != outline_bins || !(data.rho_extent > 0.0f) || !std::isfinite(data.rho_extent) ||
        std::abs(data.s_max - scatter_s_max) > 1e-6f || fractions.empty() || fractions.size() > reference_max_levels)
        return false;
    for (size_t k = 0; k < fractions.size(); ++k)
        if (!(fractions[k] > 0.0f && fractions[k] < 1.0f) || (k && !(fractions[k] > fractions[k - 1]))) return false;
    out->rho_extent = data.rho_extent;
    out->fractions.assign(fractions.begin(), fractions.end());
    const auto build = [&](ReferenceCategory* category, const float* values) {
        if (!values) return false;
        OutlineField field = {};
        bool finite = true;
        for (int y = 0; y < outline_bins; ++y) for (int x = 0; x < outline_bins; ++x) {
            const float v = values[x + (size_t)outline_bins * y];
            finite &= std::isfinite(v) && v >= 0.0f;
            field[(x + 1) + outline_dim * (y + 1)] = v;
        }
        if (!finite) return false;
        category->inside_threshold = mass_level(field, reference_inside_fraction);
        if (!(category->inside_threshold > 0.0f)) return false;
        category->thresholds.resize(fractions.size());
        category->isolines.resize(fractions.size());
        for (size_t k = 0; k < fractions.size(); ++k) {
            category->thresholds[k] = mass_level(field, fractions[k]);
            contour(&category->isolines[k], field, category->thresholds[k], data.rho_extent);
        }
        category->field.assign(field.begin(), field.end());
        category->valid = true;
        return true;
    };
    for (size_t c = 0; c < reference_category_count; ++c) {
        ReferenceCategory& category = out->categories[c];
        category.name = data.names[c] ? data.names[c] : "";
        category.systems = data.systems[c];
        category.category = c;
        out->valid |= build(&category, data.fields[c]);
    }
    for (const ReferenceSubfamilyData& sub : data.subfamilies) {
        if (sub.category >= reference_category_count || !out->categories[sub.category].valid) continue;
        ReferenceCategory subfamily;
        subfamily.name = sub.name ? sub.name : "";
        subfamily.systems = sub.systems;
        subfamily.category = sub.category;
        subfamily.overlap = sub.overlap;
        if (build(&subfamily, sub.field)) out->subfamilies.push_back(std::move(subfamily));
    }
    return out->valid;
}

const ReferenceMap& builtin_reference_map() {
    static const ReferenceMap map = [] {
        ReferenceMap m;
        build_reference_map(&m, reference_data(), reference_default_fractions);
        return m;
    }();
    return map;
}

int reference_level(const ReferenceMap& map, size_t category, float signed_rho, float s) {
    if (!map.valid || category >= reference_category_count) return 0;
    const ReferenceCategory& c = map.categories[category];
    return level_of(c, reference_value(c.field, map.rho_extent, signed_rho, s));
}

std::array<float, reference_category_count> reference_shares(const ReferenceMap& map, float signed_rho, float s) {
    std::array<float, reference_category_count> shares = {};
    if (!map.valid) return shares;
    float total = 0.0f;
    for (size_t c = 0; c < reference_category_count; ++c) {
        const ReferenceCategory& category = map.categories[c];
        if (!category.valid) continue;
        shares[c] = std::max(0.0f, reference_value(category.field, map.rho_extent, signed_rho, s));
        total += shares[c];
    }
    if (!(total > 0.0f)) return {};
    for (float& share : shares) share /= total;
    return shares;
}

int reference_subfamily_level(const ReferenceMap& map, size_t subfamily, float signed_rho, float s) {
    if (!map.valid || subfamily >= map.subfamilies.size()) return 0;
    const ReferenceCategory& c = map.subfamilies[subfamily];
    return level_of(c, reference_value(c.field, map.rho_extent, signed_rho, s));
}

std::vector<float> reference_subfamily_shares(const ReferenceMap& map, size_t category, float signed_rho, float s) {
    std::vector<float> shares(map.subfamilies.size(), 0.0f);
    float total = 0.0f;
    for (size_t i = 0; i < shares.size(); ++i) {
        const ReferenceCategory& sub = map.subfamilies[i];
        if (sub.category != category || !sub.valid) continue;
        shares[i] = std::max(0.0f, reference_value(sub.field, map.rho_extent, signed_rho, s));
        total += shares[i];
    }
    if (!(total > 0.0f)) return std::vector<float>(shares.size(), 0.0f);
    for (float& share : shares) share /= total;
    return shares;
}

void reference_image(ReferenceImage* out, const ReferenceMap& map, int size) {
    out->size = std::max(0, size);
    out->category.assign((size_t)out->size * out->size, -1);
    out->level.assign((size_t)out->size * out->size, 0);
    out->depth.assign((size_t)out->size * out->size * reference_category_count, 0.0f);
    if (!map.valid) return;
    for (int j = 0; j < out->size; ++j) for (int i = 0; i < out->size; ++i) {
        const float r = map.rho_extent * (2.0f * (i + 0.5f) / out->size - 1.0f);
        const float s = scatter_s_max * (j + 0.5f) / out->size;
        const size_t pixel = i + (size_t)out->size * j;
        float best = 0.0f;
        for (size_t c = 0; c < reference_category_count; ++c) {
            const ReferenceCategory& category = map.categories[c];
            const float v = reference_value(category.field, map.rho_extent, r, s);
            const int level = level_of(category, v);
            if (level == 0) continue;
            const float outer = std::log(category.thresholds.back()), inner = std::log(category.thresholds.front());
            const float t = inner > outer ? std::clamp((std::log(v) - outer) / (inner - outer), 0.0f, 1.0f) : 1.0f;
            out->depth[pixel * reference_category_count + c] = reference_edge_depth + (1.0f - reference_edge_depth) * t;
            if (v <= best) continue;
            best = v;
            out->category[i + (size_t)out->size * j] = (int8_t)c;
            out->level[i + (size_t)out->size * j] = (uint8_t)level;
        }
    }
}

RegionCounts reference_counts(const Scatter& scatter, const ReferenceMap& map, size_t category) {
    RegionCounts counts;
    if (!map.valid || category >= reference_category_count || !map.categories[category].valid) return counts;
    for (size_t i = 0; i < scatter.s.size() && i < scatter.signed_rho.size(); ++i) {
        const float r = scatter.signed_rho[i], s = scatter.s[i];
        if (!std::isfinite(s) || !std::isfinite(r) || s < 0.0f || s > scatter_s_max) continue;
        const ReferenceCategory& c = map.categories[category];
        if (reference_value(c.field, map.rho_extent, r, s) >= c.inside_threshold) ++counts.inside;
        else ++counts.outside;
    }
    return counts;
}

bool scatter_union_fill(std::vector<std::array<float, 2>>* triangles, std::span<const Scatter* const> scatters,
                        float rho_max, int grow_bins, float mass_fraction) {
    triangles->clear();
    OutlineField field;
    bool any = false;
    if (!union_field(&field, &any, scatters, rho_max, grow_bins, mass_fraction)) return false;
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

// Cheap deterministic jitter in [0, 1)
inline float hash_unit(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return (float)(x >> 40) * (1.0f / 16777216.0f);
}

// The continuous scatterplot of the eligible region. Each cell between interior samples, with a corner in scope and
// at or under s_max, is supersampled on a jittered m^3 lattice, m following how many density bins the cell's corner
// values span. At each subsample: trilinear rho, gradient and lambda2 from the corners; s from the interpolated rho
// and gradient; signed rho from the interpolated lambda2, so a cell crossing lambda2 = 0 splits between the two
// sides rather than smearing across signed rho = 0. A subsample counts when its rho is in (0, rho_max], its s at
// most s_cut, and its nearest sample is kept (in scope, nearest to one of the anchors); it adds its share of the
// cell's volume, with by_anchor also to that sample's anchor.
template <typename Keep>
bool continuous_density(ScatterDensity* out, const Lattice& l, const std::vector<float>& S, const Params& params,
                        double s_cut, Keep&& keep, const std::vector<uint32_t>& owner,
                        const std::atomic<int32_t>* cancel) {
    constexpr int bins = scatter_density_bins, max_m = 4;
    const int dx = l.dx, dy = l.dy, dz = l.dz;
    const float rho_extent = (float)params.rho_max;
    const double bin_w = 2.0 * rho_extent / bins, bin_h = scatter_s_max / bins;
    const double cell_volume = l.h[0] * l.h[1] * l.h[2];
    *out = {};
    out->rho_extent = rho_extent;
    std::vector<double> total((size_t)bins * bins, 0.0);
    std::unordered_map<uint64_t, double> total_by_anchor;
    const bool by_anchor = params.by_anchor && !owner.empty();
    std::mutex merge;
    std::atomic<bool> stop = false;
    // Cells (i, j, k) .. (i + 1, j + 1, k + 1), all corners interior: i in [1, dx - 3]
    parallel_slabs(dz - 3, [&](int kk) {
        if (stop.load(std::memory_order_relaxed)) return;
        if (cancel && cancel->load(std::memory_order_relaxed)) { stop = true; return; }
        const int k = kk + 1;
        // Per corner of the two planes: rho, gradient, lambda2, s; NaN rho until computed
        const size_t plane = (size_t)dx * dy;
        // Allocated with the slab's first candidate: scoped extractions leave most slabs empty
        std::vector<std::array<float, 6>> corner;
        auto corner_at = [&](int i, int j, int c) -> const std::array<float, 6>& {
            if (corner.empty()) corner.assign(2 * plane, {NAN, 0, 0, 0, 0, 0});
            std::array<float, 6>& q = corner[(size_t)c * plane + (size_t)i + (size_t)dx * j];
            if (std::isnan(q[0])) {
                double r, grad[3], hess[6];
                derivatives(&r, grad, hess, l, i, j, k + c);
                const double s = critical_point_distribution::reduced_density_gradient(r, grad);
                q = {(float)r, (float)grad[0], (float)grad[1], (float)grad[2],
                     (float)critical_point_distribution::symmetric_eigenvalues(hess)[1],
                     std::isfinite(s) ? (float)std::min(s, (double)s_out) : s_out};
            }
            return q;
        };
        std::vector<double> local;
        std::unordered_map<uint64_t, double> local_by_anchor;
        size_t cells = 0, subsamples = 0;
        double volume = 0.0;
        for (int j = 1; j < dy - 2; ++j) for (int i = 1; i < dx - 2; ++i) {
            size_t vertex[8];
            bool candidate = false;
            for (int c = 0; c < 8; ++c) {
                vertex[c] = l.at(i + (c & 1), j + ((c >> 1) & 1), k + (c >> 2));
                candidate |= S[vertex[c]] <= scatter_s_max && keep(vertex[c]);
            }
            if (!candidate) continue;
            if (local.empty()) local.assign((size_t)bins * bins, 0.0);
            const std::array<float, 6>* q[8];
            float rlo = INFINITY, rhi = -INFINITY, slo = INFINITY, shi = -INFINITY;
            bool both_signs = false, negative = false;
            for (int c = 0; c < 8; ++c) {
                q[c] = &corner_at(i + (c & 1), j + ((c >> 1) & 1), c >> 2);
                rlo = std::min(rlo, (*q[c])[0]);
                rhi = std::max(rhi, (*q[c])[0]);
                slo = std::min(slo, std::min((*q[c])[5], scatter_s_max));
                shi = std::max(shi, std::min((*q[c])[5], scatter_s_max));
                const bool neg = (*q[c])[4] < 0.0f;
                if (c && neg != negative) both_signs = true;
                negative = c ? negative : neg;
            }
            const double span = std::max({(rhi - rlo) / bin_w, (shi - slo) / bin_h, both_signs ? (double)max_m : 1.0});
            const int m = std::clamp((int)std::ceil(span), 1, max_m);
            const double weight = cell_volume / (m * m * m);
            ++cells;
            for (int c = 0, n = 0; c < m; ++c) for (int b = 0; b < m; ++b) for (int a = 0; a < m; ++a, ++n) {
                const uint64_t seed = ((uint64_t)l.at(i, j, k) << 9) + (uint64_t)n * 3;
                const float u[3] = {(a + hash_unit(seed)) / m, (b + hash_unit(seed + 1)) / m, (c + hash_unit(seed + 2)) / m};
                // The nearest sample decides the scope, as each point of the scatter stands for its voxel
                const int nearest = (u[0] >= 0.5f) | ((u[1] >= 0.5f) << 1) | ((u[2] >= 0.5f) << 2);
                if (!keep(vertex[nearest])) continue;
                float v[5] = {0, 0, 0, 0, 0};
                for (int corner_index = 0; corner_index < 8; ++corner_index) {
                    float w = 1.0f;
                    for (int axis = 0; axis < 3; ++axis)
                        w *= (corner_index >> axis) & 1 ? u[axis] : 1.0f - u[axis];
                    for (int f = 0; f < 5; ++f) v[f] += w * (*q[corner_index])[f];
                }
                ++subsamples;
                if (!(v[0] > 0.0f) || v[0] > rho_extent) continue;
                const double grad[3] = {v[1], v[2], v[3]};
                const double s = critical_point_distribution::reduced_density_gradient(v[0], grad);
                if (!(s >= 0.0) || s > s_cut) continue;
                const double x = v[4] < 0.0f ? -v[0] : v[0];
                const int bx = std::clamp((int)((x + rho_extent) / bin_w), 0, bins - 1);
                const int by = std::clamp((int)(s / bin_h), 0, bins - 1);
                local[bx + (size_t)bins * by] += weight;
                if (by_anchor)
                    local_by_anchor[((uint64_t)owner[vertex[nearest]] << 32) | (uint64_t)(bx + bins * by)] += weight;
                volume += weight;
            }
        }
        if (local.empty()) return;
        std::lock_guard<std::mutex> lock(merge);
        for (size_t b = 0; b < local.size(); ++b) total[b] += local[b];
        for (const auto& [key, part] : local_by_anchor) total_by_anchor[key] += part;
        out->cells += cells;
        out->subsamples += subsamples;
        out->total_volume += volume;
    });
    if (stop || (cancel && cancel->load(std::memory_order_relaxed))) {
        *out = {};
        return false;
    }
    out->volume.assign(total.begin(), total.end());
    out->by_anchor.reserve(total_by_anchor.size());
    for (const auto& [key, part] : total_by_anchor)
        out->by_anchor.push_back({(uint32_t)(key >> 32), (uint32_t)(key & 0xffffffffu), (float)part});
    std::sort(out->by_anchor.begin(), out->by_anchor.end(), [](const AnchorBin& a, const AnchorBin& b) {
        return a.anchor != b.anchor ? a.anchor < b.anchor : a.bin < b.bin;
    });
    return true;
}

// The nearest of the anchors, near_xyz's then other_xyz's, by a uniform bucket grid; ties go to the lower index, so
// to near_xyz. Index i < near count is near_xyz's point i.
class AnchorLocator {
public:
    AnchorLocator(std::span<const double> near_xyz, std::span<const double> other_xyz) : near_(near_xyz), other_(other_xyz) {
        count_ = (near_.size() + other_.size()) / 3;
        if (!count_) return;
        double lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
        for (size_t i = 0; i < count_; ++i) for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], at(i)[a]);
            hi[a] = std::max(hi[a], at(i)[a]);
        }
        cell_ = 2.0;
        for (int a = 0; a < 3; ++a) cell_ = std::max(cell_, (hi[a] - lo[a]) / 64.0);
        for (int a = 0; a < 3; ++a) {
            origin_[a] = lo[a];
            dim_[a] = (int)std::floor((hi[a] - lo[a]) / cell_) + 1;
        }
        start_.assign((size_t)dim_[0] * dim_[1] * dim_[2] + 1, 0);
        std::vector<uint32_t> bucket_of(count_);
        for (size_t i = 0; i < count_; ++i) {
            int b[3];
            for (int a = 0; a < 3; ++a) b[a] = std::clamp((int)((at(i)[a] - origin_[a]) / cell_), 0, dim_[a] - 1);
            bucket_of[i] = (uint32_t)(b[0] + dim_[0] * (b[1] + dim_[1] * b[2]));
            ++start_[bucket_of[i] + 1];
        }
        for (size_t b = 1; b < start_.size(); ++b) start_[b] += start_[b - 1];
        items_.resize(count_);
        std::vector<uint32_t> fill(start_.begin(), start_.end() - 1);
        for (size_t i = 0; i < count_; ++i) items_[fill[bucket_of[i]]++] = (uint32_t)i;
    }
    uint32_t nearest(const double p[3]) const {
        if (!count_) return UINT32_MAX;
        int c[3];
        for (int a = 0; a < 3; ++a) c[a] = std::clamp((int)std::floor((p[a] - origin_[a]) / cell_), 0, dim_[a] - 1);
        double best = INFINITY;
        uint32_t best_index = UINT32_MAX;
        for (int r = 0;; ++r) {
            int lo[3], hi[3];
            bool covers = true;
            for (int a = 0; a < 3; ++a) {
                lo[a] = std::max(0, c[a] - r);
                hi[a] = std::min(dim_[a] - 1, c[a] + r);
                covers &= lo[a] == 0 && hi[a] == dim_[a] - 1;
            }
            for (int z = lo[2]; z <= hi[2]; ++z) for (int y = lo[1]; y <= hi[1]; ++y) for (int x = lo[0]; x <= hi[0]; ++x) {
                // Only the shell r: the inner buckets were visited before
                if (std::max({std::abs(x - c[0]), std::abs(y - c[1]), std::abs(z - c[2])}) != r) continue;
                const size_t b = (size_t)x + (size_t)dim_[0] * ((size_t)y + (size_t)dim_[1] * z);
                for (uint32_t k = start_[b]; k < start_[b + 1]; ++k) {
                    const uint32_t i = items_[k];
                    const double* q = at(i);
                    const double d2 = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
                    if (d2 < best || (d2 == best && i < best_index)) { best = d2; best_index = i; }
                }
            }
            if (covers) break;
            // Any anchor outside the visited buckets is at least this far
            double bound = INFINITY;
            for (int a = 0; a < 3; ++a) {
                if (lo[a] > 0) bound = std::min(bound, p[a] - (origin_[a] + lo[a] * cell_));
                if (hi[a] < dim_[a] - 1) bound = std::min(bound, origin_[a] + (hi[a] + 1) * cell_ - p[a]);
            }
            if (bound > 0.0 && best < bound * bound) break;
        }
        return best_index;
    }
private:
    const double* at(size_t i) const { return i * 3 < near_.size() ? &near_[i * 3] : &other_[i * 3 - near_.size()]; }
    std::span<const double> near_, other_;
    size_t count_ = 0;
    double cell_ = 1.0, origin_[3] = {0, 0, 0};
    int dim_[3] = {0, 0, 0};
    std::vector<uint32_t> start_, items_;
};

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
        !(params.scatter_s >= 0.0) || !std::isfinite(params.scatter_s) ||
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
    // With connected scope, the surface keeps the components grown below s_iso; the scatter those grown below
    // scatter_s (s_iso when unset) and their one-voxel halo, for sub-voxel sampling at the components' edge
    const double scatter_s = params.scatter_s > 0.0 ? params.scatter_s : params.s_iso;
    std::vector<uint8_t> selected;
    std::vector<uint8_t> halo;
    if (params.connected && !params.near_xyz.empty()) {
        const auto grow = [&](std::vector<uint8_t>* grown, std::vector<uint8_t>* grown_halo, double s_max) -> bool {
            std::vector<uint8_t>& sel = *grown;
            sel.assign(n, 0);
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
                            if (!(S[v] < s_max) || (!allowed.empty() && !allowed[v])) continue;
                            const int ijk[3] = {i, j, k};
                            double d2 = 0.0;
                            for (int a = 0; a < 3; ++a) {
                                const double d = grid.origin[a] + ijk[a] * grid.spacing[a] - params.near_xyz[p + a];
                                d2 += d * d;
                            }
                            if (d2 < closest) { closest = d2; seed = v; }
                        }
                if (seed != n && !sel[seed]) { sel[seed] = 1; queue.push_back((uint32_t)seed); }
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
                    if (!sel[u] && S[u] < s_max && (allowed.empty() || allowed[u])) {
                        sel[u] = 1;
                        queue.push_back((uint32_t)u);
                    }
                }
            }
            if (!grown_halo) return true;
            std::vector<uint8_t>& h = *grown_halo;
            h = sel;
            for (size_t head = 0; head < queue.size(); ++head) {
                if ((head & 4095) == 0 && cancelled()) return false;
                const size_t v = queue[head];
                const int i = (int)(v % dx), j = (int)((v / dx) % dy), k = (int)(v / ((size_t)dx * dy));
                for (int c = -1; c <= 1; ++c) for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a)
                    h[l.at(i + a, j + b, k + c)] = 1;
            }
            return true;
        };
        const bool shared = scatter_s == params.s_iso;
        if ((params.surface || (shared && scatter)) && !grow(&selected, shared && scatter ? &halo : nullptr, params.s_iso)) return false;
        if (scatter && !shared) {
            std::vector<uint8_t> scatter_selected;
            if (!grow(&scatter_selected, &halo, scatter_s)) return false;
        }
    }
    if (scatter) {
        const auto in_scope = [&](size_t v) {
            return (allowed.empty() || allowed[v]) && (halo.empty() || halo[v]);
        };
        // Each sample belongs to its nearest anchor; those nearer to one of other_xyz are not kept
        const size_t near_count = params.near_xyz.size() / 3;
        const bool own = near_count && (params.by_anchor || !params.other_xyz.empty());
        std::vector<uint32_t> owner;
        if (own) {
            for (double x : params.other_xyz) if (!std::isfinite(x)) return false;
            if (params.other_xyz.size() % 3 != 0) return false;
            const AnchorLocator locator(params.near_xyz, params.other_xyz);
            owner.assign(n, UINT32_MAX);
            parallel_slabs(dz, [&](int k) {
                if (cancelled()) return;
                for (int j = 0; j < dy; ++j) for (int i = 0; i < dx; ++i) {
                    const size_t v = l.at(i, j, k);
                    if (!(S[v] <= scatter_s_max) || !in_scope(v)) continue;
                    const double p[3] = {grid.origin[0] + i * grid.spacing[0], grid.origin[1] + j * grid.spacing[1],
                                         grid.origin[2] + k * grid.spacing[2]};
                    owner[v] = locator.nearest(p);
                }
            });
            if (cancelled()) return false;
        }
        const auto keep = [&](size_t v) { return in_scope(v) && (!own || owner[v] < near_count); };
        std::vector<size_t> reservoir;
        reservoir.reserve(scatter_max_points);
        std::mt19937 rng(0);
        for (int k = 1; k < dz - 1; ++k) {
            if (cancelled()) return false;
            for (int j = 1; j < dy - 1; ++j) for (int i = 1; i < dx - 1; ++i) {
                const size_t v = l.at(i, j, k);
                if (!(S[v] <= scatter_s) || !keep(v)) continue;
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
        if (own) scatter->anchor.reserve(reservoir.size());
        for (size_t v : reservoir) {
            if (cancelled()) return false;
            const int i = (int)(v % dx), j = (int)((v / dx) % dy), k = (int)(v / ((size_t)dx * dy));
            const float r = (float)signed_density(l, i, j, k);
            if (!std::isfinite(r)) continue;
            scatter->signed_rho.push_back(r);
            scatter->s.push_back(S[v]);
            scatter->voxel.push_back((uint32_t)v);
            if (own) scatter->anchor.push_back(owner[v]);
        }
        if (!continuous_density(&scatter->density, l, S, params, scatter_s, keep, owner, cancel)) return false;
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
