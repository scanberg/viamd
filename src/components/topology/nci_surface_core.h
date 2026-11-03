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
#include <string_view>
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
// The CP classification is the one source of the types: a CP's type is its branch of the CP tree (the contact class,
// ChemicalLabel::interaction, under Contacts), mapped by this table. Every pi contact (pi...pi, X-H...pi, cation,
// anion and lone pair...pi) shares the pi type, and so its radius; the sigma-hole bonds (halogen, chalcogen, pnictogen,
// tetrel) share the halogen-bond type; C-H...X, H...H, dihydrogen bonds, metal and type I halogen contacts are other
// contacts. A class not in the table is an other contact.
struct ContactClass {
    const char* name;   // ChemicalLabel::interaction
    Family family;
};
std::span<const ContactClass> contact_classes();
Family contact_family(std::string_view interaction);
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

// The scatter extraction, the same for every CP type, the reference map and the app's defaults: the regions grown
// from the CPs through s < reference_scatter_s within reference_radius, rho <= reference_rho_max, both signs of
// lambda2; only their samples with s <= reference_scatter_s; each sample belonging to its nearest CP.
constexpr double reference_scatter_s = 1.0;
constexpr double reference_rho_max = 0.07;
constexpr double reference_radius = 2.7;   // Bohr

// The reference extraction's parameters (scatter only, by anchor), around near_xyz, with other_xyz competing for the
// samples
struct Params;
Params reference_params(std::span<const double> near_xyz, std::span<const double> other_xyz = {});

struct Defaults {
    double radius = reference_radius;
    bool established = false;
};
// CP-adaptive mode's display radii: calibrated for H-bonds, rings and pi-stacking; the other types keep the shared
// radius. The scatter statistics and the reference never use them.
Defaults family_defaults(Families families);

// The displayed isovalue and density cutoff (every mode), the radius (one for every CP type in Uniform mode; per type
// in CP-adaptive mode), and the colour range of the single NCI colour scale. The reference statistics use their own
// extraction (reference_params), independent of these.
struct UniformDefaults {
    double s_iso = 0.7;
    // The scatter's s, rho and R are the reference's, so the scatter plot covers what the reference map does
    double scatter_s = reference_scatter_s;
    double rho_max = reference_rho_max;
    double radius = reference_radius;
    double color_range = 0.03;
};
UniformDefaults uniform_defaults();

struct Scatter;

// Region edges are highest-density contours: the smallest area of the smoothed (signed rho, s) density holding this
// fraction of its volume. The core line holds the smaller fraction.
constexpr float region_mass_fraction = 0.9f;
constexpr float region_core_mass_fraction = 0.5f;

// A data-derived region enclosing the union of each reference scatter's 90%-of-volume contour, with the union of their
// 50% cores. The reference contours are irregular; they are not fitted to rectangular or triangular bounds.
struct ExpectedRegion {
    float rho_extent = 0.0f;
    std::vector<uint8_t> mask;
    std::vector<float> density;
    std::vector<std::array<float, 4>> outline;
    std::vector<std::array<float, 4>> core_outline;
    bool valid = false;
};
ExpectedRegion expected_region(std::span<const Scatter* const> scatters, float rho_extent, int grow_bins = 1);
void expected_region_fill(std::vector<std::array<float, 2>>* triangles, const ExpectedRegion& region);
bool in_expected_region(const ExpectedRegion& region, float signed_rho, float s);
// Of the eligible scatter samples, how many fall inside the data-derived region and how many outside
struct RegionCounts {
    size_t inside = 0;
    size_t outside = 0;
};
RegionCounts expected_region_counts(const Scatter& scatter, const ExpectedRegion& region);

// The built-in reference of the NCI scatter plot, Ramachandran-like: per CP category (hydrogen bonds, rings,
// pi-stacking, other), the pooled (signed rho, s) density of reference extractions over many calculations.
// Generated by tools/nci_reference/bake_reference.py into nci_reference_data_core.cpp.
constexpr size_t reference_category_count = 4;
// The reference category of a CP type: hydrogen bonds, rings, pi-stacking, or other (halogen bonds, other contacts,
// cages); reference_category_count for the types without NCI or a category (bonds, attractors, unassigned)
size_t reference_category(Family family);
// The type whose settings (the radius) a reference category's extraction uses
Family reference_settings_family(size_t category);
// "Hydrogen bonds", "Rings", "\xCF\x80 interactions", "Other"; the baked data's and the batch tools' names
const char* reference_category_name(size_t category);
// A category's subfamily: a CP-tree branch under it (N-H···O, aromatic rings, C-H···pi, cages, ...), pooled from the
// reference systems whose CPs of that category all belong to the branch. overlap: of its 90 % area with its category's
// (intersection over union)
struct ReferenceSubfamilyData {
    size_t category = reference_category_count;
    const char* name = "";
    int systems = 0;
    float overlap = 0.0f;
    const float* field = nullptr;
};
struct ReferenceData {
    int bins = 0;
    float rho_extent = 0.0f;
    float s_max = 0.0f;
    std::array<const char*, reference_category_count> names = {};
    std::array<int, reference_category_count> systems = {};
    // field[x + bins * y]: x the signed-rho bin over [-rho_extent, rho_extent], y the s bin over [0, s_max]; unit sum
    std::array<const float*, reference_category_count> fields = {};
    const char* source = "";
    std::span<const ReferenceSubfamilyData> subfamilies = {};
};
const ReferenceData& reference_data();

// The default reference levels, as the Ramachandran layers': the highest-density areas holding these fractions of a
// category's volume, increasing (the core first)
constexpr std::array<float, 3> reference_default_fractions = {0.5f, 0.9f, 0.98f};
// A sample counts as inside the reference within its category's 90 % area, whatever levels are drawn
constexpr float reference_inside_fraction = 0.9f;
constexpr int reference_max_levels = 8;

struct ReferenceCategory {
    const char* name = "";
    int systems = 0;
    // The density on bins 1..scatter_field_bins with an empty border, its value at each level (the map's fractions,
    // so decreasing) and at the inside fraction, and an isoline per level
    std::vector<float> field;
    std::vector<float> thresholds;
    float inside_threshold = 0.0f;
    std::vector<std::vector<std::array<float, 4>>> isolines;
    // Subfamilies only: the category they belong to and their overlap with it
    size_t category = reference_category_count;
    float overlap = 0.0f;
    bool valid = false;
};
struct ReferenceMap {
    float rho_extent = 0.0f;
    std::vector<float> fractions;
    std::array<ReferenceCategory, reference_category_count> categories;
    // In the data's order, grouped by category
    std::vector<ReferenceCategory> subfamilies;
    bool valid = false;
};
// At fractions in (0, 1), sorted increasing, at most reference_max_levels. Returns false when the data does not match
// the scatter plot's bins and s range, has no density, or the fractions are invalid.
bool build_reference_map(ReferenceMap* out, const ReferenceData& data, std::span<const float> fractions);
// The built-in data at the default fractions
const ReferenceMap& builtin_reference_map();
// How deep a point lies in a category: the number of its levels the point is inside, 0 outside the outermost
int reference_level(const ReferenceMap& map, size_t category, float signed_rho, float s);
// Which category a point most likely belongs to: each category's share of the categories' reference densities there,
// every category's map weighing the same; all 0 where no category has density or the point is off the map
std::array<float, reference_category_count> reference_shares(const ReferenceMap& map, float signed_rho, float s);
// How deep a point lies in a subfamily, as reference_level
int reference_subfamily_level(const ReferenceMap& map, size_t subfamily, float signed_rho, float s);
// Within a category, each of its subfamilies' share of their reference densities at the point (each subfamily weighing
// the same), in map.subfamilies' order, 0 for the other categories' subfamilies; all 0 where none has density
std::vector<float> reference_subfamily_shares(const ReferenceMap& map, size_t category, float signed_rho, float s);
// A size x size image over [-rho_extent, rho_extent] x [0, scatter_s_max], row 0 at s = 0. category: of the categories
// a pixel is inside (level >= 1), the one with the highest density, -1 for none; level: its level. depth[pixel *
// reference_category_count + c]: 0 outside category c's outermost level, else rising in log density from
// reference_edge_depth on that isoline to 1 at its innermost one and above, for a combined, graded map.
constexpr float reference_edge_depth = 0.1f;
struct ReferenceImage {
    int size = 0;
    std::vector<int8_t> category;
    std::vector<uint8_t> level;
    std::vector<float> depth;
};
void reference_image(ReferenceImage* out, const ReferenceMap& map, int size);
// The scatter's eligible samples inside and outside the category's reference region (its 90 % area)
RegionCounts reference_counts(const Scatter& scatter, const ReferenceMap& map, size_t category);

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
    // The scatter's s (0: s_iso): connected scope grows the scatter's components below it, and the scatter keeps only
    // samples with s at most it, so the scatter plot can reach further than the surface
    double scatter_s = 0.0;
    double rho_max = 0.05;   // au
    // Anchor points in Bohr. Radius-only scope clips their union of spheres; connected scope seeds low-s components.
    std::span<const double> near_xyz = {};   // 3 per point
    // Competing anchor points (Bohr, 3 per point; may repeat near_xyz's): the scatter keeps only samples nearer to an
    // anchor of near_xyz than to any of them, so each sample belongs to one CP. The surface ignores them.
    std::span<const double> other_xyz = {};
    double near_radius = 3.0;
    bool connected = false; // Six-connected components, seeded within two voxels of each anchor
    bool radius_cap = true;
    bool surface = true; // False for a scatter-only extraction
    // Label each scatter sample with its nearest anchor (Scatter::anchor, ScatterDensity::by_anchor)
    bool by_anchor = false;
};

struct Mesh {
    std::vector<float>    position;   // xyz per vertex, Bohr
    std::vector<float>    normal;     // xyz per vertex, unit, towards higher s (out of the low gradient region)
    std::vector<float>    value;      // per vertex: sign(lambda2) rho, au
    std::vector<uint32_t> index;      // three per triangle, counter-clockwise seen from where the normals point
};

constexpr size_t scatter_max_points = 50000;
constexpr float scatter_s_max = 2.0f;
// Bins per axis of the outline, area and expected-region fields, and of the continuous density they are made from
constexpr int scatter_field_bins = 96;
constexpr int scatter_density_bins = 2 * scatter_field_bins;

// The continuous scatterplot: the volume, Bohr^3, of the eligible region per (signed rho, s) bin, from every grid cell
// supersampled (trilinear rho, gradient and lambda2; sign(lambda2) taken after interpolation), not from the capped
// points. volume[x + bins * y] covers signed rho -rho_extent + 2 rho_extent (x, x + 1) / bins and
// s scatter_s_max (y, y + 1) / bins.
struct AnchorBin {
    uint32_t anchor = 0;   // index of the point in Params::near_xyz
    uint32_t bin = 0;      // x + bins * y
    float volume = 0.0f;
};
struct ScatterDensity {
    float rho_extent = 0.0f;
    std::vector<float> volume;
    // With Params::by_anchor: the volume split by nearest anchor, sorted by anchor then bin
    std::vector<AnchorBin> by_anchor;
    double total_volume = 0.0;
    size_t cells = 0;
    size_t subsamples = 0;
};

struct Scatter {
    std::vector<float> signed_rho;
    std::vector<float> s;
    std::vector<uint32_t> voxel;
    // With Params::by_anchor: per point, the index of its nearest anchor in Params::near_xyz
    std::vector<uint32_t> anchor;
    size_t eligible_count = 0;
    // Empty for scatters built from points only; then the fields use the points' histogram.
    ScatterDensity density;
};

struct CriticalPointAnchor {
    std::array<double, 3> position = {};
    Family family = Family::Unassigned;
    bool selected = false;
};

struct CriticalPointGroups {
    std::array<std::vector<double>, family_count> all_anchors;
    std::array<std::vector<double>, family_count> selected_anchors;
};

// The part of a scatter extracted with by_anchor whose nearest anchors are kept (keep[anchor] != 0): its points
// and density. eligible_count scales with the kept share of the points.
Scatter scatter_subset(const Scatter& scatter, std::span<const uint8_t> keep);
// Appends src's points and density (same rho extent) to dst
void scatter_append(Scatter* dst, const Scatter& src);

// Groups CP coordinates for the overview and selection preview. Unclassified points use Unassigned.
CriticalPointGroups group_critical_points(std::span<const CriticalPointAnchor> points);

// Segments of the 90%-of-volume contour of the smoothed (signed rho, s) density: the continuous density when the
// scatter has one, the points' histogram otherwise.
// Returns false for invalid bounds or mismatched/nonfinite sample arrays; empty input yields an empty contour.
bool scatter_density_outline(std::vector<std::array<float, 4>>* out, const Scatter& scatter, float rho_max);
// The outline of the union of the scatters' 90%-of-volume regions (each at its own level), grown by grow_bins histogram
// bins so it encloses each scatter's own outline with a margin. Same return values as scatter_density_outline.
// With mass_fraction, the union of the regions holding that fraction instead.
bool scatter_union_outline(std::vector<std::array<float, 4>>* out, std::span<const Scatter* const> scatters, float rho_max,
                           int grow_bins, float mass_fraction = region_mass_fraction);
// The inside of that union as triangles, three (signed rho, s) points each. Same return values.
bool scatter_union_fill(std::vector<std::array<float, 2>>* triangles, std::span<const Scatter* const> scatters,
                        float rho_max, int grow_bins, float mass_fraction = region_mass_fraction);
// The scatters' smoothed sample histogram together, relative to its peak: values[x + bins * y] is the bin of signed rho
// -rho_max + 2 rho_max (x, x + 1) / bins and s scatter_s_max (y, y + 1) / bins. Returns false on invalid input.
struct Heatmap {
    int bins = 0;
    float rho_max = 0.0f;
    std::vector<float> values;
};
bool scatter_heatmap(Heatmap* out, std::span<const Scatter* const> scatters, float rho_max);

// Returns false when cancelled (cancel set non-zero) or the input is invalid; progress (optional) goes from 0 to 1.
// Optional scatter: a deterministic reservoir of grid samples with s <= scatter_s (s_iso when unset), respecting
// rho_max, near_xyz and other_xyz. It samples the volume, not the surface. Connected scope keeps the components grown
// below scatter_s and their one-voxel halo, then only the samples at most scatter_s. voxel indexes into grid.
// The scatter's continuous density covers the same region: each point inside a cell belongs to its nearest grid
// sample's scope and anchor.
bool extract(Mesh* out, const float* rho, const Grid& grid, const Params& params,
             const std::atomic<int32_t>* cancel = nullptr, std::atomic<float>* progress = nullptr,
             Scatter* scatter = nullptr);

// The reduced density gradient and sign(lambda2) rho of sample (i, j, k), not on the grid's outer layer
void sample(double* s, double* signed_rho, const float* rho, const Grid& grid, int i, int j, int k);

}  // namespace nci_surface
