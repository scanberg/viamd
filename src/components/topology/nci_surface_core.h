#pragma once

// UI-independent extraction and analysis of non-covalent interaction surfaces.
// NCI surfaces (Johnson et al., JACS 2010): the isosurface of the reduced density gradient
// s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)) where the density is low (rho <= rho_max), coloured by sign(lambda2) rho,
// lambda2 the second eigenvalue of the density's Hessian: negative attractive, near zero van der Waals, positive repulsive.
//
// From a density sampled on a regular, axis-aligned grid: the gradient and the Hessian by central finite differences,
// the surface by Surface Nets on s. Samples denser than rho_max, and the grid's outer layer, count as s = infinity (out),
// as NCIPLOT's density cutoff does, so the low gradient regions of covalent bonds and of the atoms are not shown.

#include <atomic>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace critical_point_distribution { struct ChemicalLabel; }

namespace nci_surface {

enum class View { Overview, Preview, Saved };
inline bool follow_selection(View* view, uint64_t* previous_key, uint64_t key, bool selected) {
    if (*previous_key == key) return false;
    *previous_key = key;
    *view = selected ? View::Preview : View::Overview;
    return true;
}

enum class Family : uint8_t {
    HydrogenBond, Ring, PiStacking, HalogenBond, OtherContact, Bond, BondPath, Atom, Cage, Unassigned, Count
};
constexpr size_t family_count = (size_t)Family::Count;
using Families = uint32_t;
constexpr Families family_bit(Family family) { return 1u << (unsigned)family; }
// The CP type, as the user sees it
Family family_of(const critical_point_distribution::ChemicalLabel& label);
// The type whose parameters the CP uses: its own, except for a ring or cage closed by contacts, which uses its
// closing interaction's. Internal: never shown.
Family parameter_family_of(const critical_point_distribution::ChemicalLabel& label);
const char* family_name(Family family);
bool single_family(Families families);
// Whether CPs of the family anchor NCI surfaces: not bonds, bond paths and attractors, whose density is far above
// any NCI cutoff.
bool has_nci(Family family);

struct Defaults {
    double s_iso = 0.5;
    double rho_max = 0.05;
    double radius = 2.7;
    bool established = false;
};
// Calibrated single-family H-bond, ring and pi-stacking presets; the other types (halogen bonds TBD) keep
// generic defaults.
Defaults family_defaults(Families families);

// Uniform mode: one isovalue and density cutoff for every CP type, and the colour range of the single NCI colour scale.
// Placeholder values until they are defined.
struct UniformDefaults {
    double s_iso = 0.7;
    double rho_max = 0.07;
    double color_range = 0.03;
};
UniformDefaults uniform_defaults();

// The reference of the scatter plot: where the spike of each CP type is expected, sign(lambda2)rho between rho_min and
// rho_max below s_max. Provisional boxes from the calibration systems (AT-AT, biphenyl, anthracene dimer), until a
// reference map replaces them; expected regions, not a validation.
struct ExpectedRegion {
    float rho_min = 0.0f;
    float rho_max = 0.0f;
    float s_max = 0.0f;
    bool valid = false;
};
ExpectedRegion expected_region(Family family);
bool in_expected_region(const ExpectedRegion& region, float signed_rho, float s);
// Of the samples below the region's s_max, how many fall inside the region and how many outside
struct RegionCounts {
    size_t inside = 0;
    size_t outside = 0;
};

// Sample (i, j, k) is at origin + (i, j, k) * spacing, Bohr; rho[i + dim[0] * (j + dim[1] * k)]
struct Grid {
    int    dim[3]     = {0, 0, 0};
    double origin[3]  = {0, 0, 0};
    double spacing[3] = {0, 0, 0};
};

// Supplies density values at Grid::origin + (i, j, k) * Grid::spacing, in the same order as extract().
struct DensityProvider {
    void* user_data = nullptr;
    bool (*evaluate)(void* user_data, const Grid& grid, std::vector<float>* density) = nullptr;
};

struct Params {
    double s_iso   = 0.5;
    double rho_max = 0.05;   // au
    // Anchor points in Bohr. Radius-only scope clips their union of spheres; connected scope seeds low-s components.
    std::span<const double> near_xyz = {};   // 3 per point
    double near_radius = 3.0;
    bool connected = false; // Six-connected components, seeded within two voxels of each anchor
    bool radius_cap = true;
    bool surface = true; // False for a scatter-only extraction
};

struct Mesh {
    std::vector<float>    position;   // xyz per vertex, Bohr
    std::vector<float>    normal;     // xyz per vertex, unit, towards higher s (out of the low gradient region)
    std::vector<float>    value;      // per vertex: sign(lambda2) rho, au
    std::vector<uint32_t> index;      // three per triangle, counter-clockwise seen from where the normals point
};

constexpr size_t scatter_max_points = 50000;
constexpr float scatter_s_max = 2.0f;

struct Scatter {
    std::vector<float> signed_rho;
    std::vector<float> s;
    std::vector<uint32_t> voxel;
    size_t eligible_count = 0;
};
RegionCounts expected_region_counts(const Scatter& scatter, Family family);

struct CriticalPointAnchor {
    std::array<double, 3> position = {};
    Family family = Family::Unassigned;
    bool selected = false;
};

struct CriticalPointGroups {
    std::array<std::vector<double>, family_count> all_anchors;
    std::array<std::vector<double>, family_count> selected_anchors;
};

// Groups CP coordinates for the overview and selection preview. Unclassified points use Unassigned.
CriticalPointGroups group_critical_points(std::span<const CriticalPointAnchor> points);

// Segments of the 10%-of-peak contour of the smoothed (signed rho, s) sample histogram.
// Returns false for invalid bounds or mismatched/nonfinite sample arrays; empty input yields an empty contour.
bool scatter_density_outline(std::vector<std::array<float, 4>>* out, const Scatter& scatter, float rho_max);
// The outline of the union of the scatters' 10%-of-peak regions (each at its own peak), grown by grow_bins histogram
// bins so it encloses each scatter's own outline with a margin. Same return values as scatter_density_outline.
bool scatter_union_outline(std::vector<std::array<float, 4>>* out, std::span<const Scatter* const> scatters, float rho_max,
                           int grow_bins);
// The inside of that union as triangles, three (signed rho, s) points each. Same return values.
bool scatter_union_fill(std::vector<std::array<float, 2>>* triangles, std::span<const Scatter* const> scatters,
                        float rho_max, int grow_bins);
// The scatters' smoothed sample histogram together, relative to its peak: values[x + bins * y] is the bin of signed rho
// -rho_max + 2 rho_max (x, x + 1) / bins and s scatter_s_max (y, y + 1) / bins. Returns false on invalid input.
struct Heatmap {
    int bins = 0;
    float rho_max = 0.0f;
    std::vector<float> values;
};
bool scatter_heatmap(Heatmap* out, std::span<const Scatter* const> scatters, float rho_max);

// Returns false when cancelled (cancel set non-zero) or the input is invalid; progress (optional) goes from 0 to 1.
// Optional scatter: a deterministic reservoir of grid samples with s <= 2, respecting rho_max and near_xyz.
// It samples the volume, not the surface. Connected scope keeps the selected low-s components
// and their one-voxel halo; its scatter therefore follows s_iso. voxel indexes into grid.
bool extract(Mesh* out, const float* rho, const Grid& grid, const Params& params,
             const std::atomic<int32_t>* cancel = nullptr, std::atomic<float>* progress = nullptr,
             Scatter* scatter = nullptr);

// The reduced density gradient and sign(lambda2) rho of sample (i, j, k), not on the grid's outer layer
void sample(double* s, double* signed_rho, const float* rho, const Grid& grid, int i, int j, int k);

}  // namespace nci_surface
