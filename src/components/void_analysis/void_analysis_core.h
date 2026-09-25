#pragma once

// The computation behind the void analysis component, linkable without an application: no ImGui, no
// events, no ApplicationState. Everything here takes plain data and returns plain data, which is what
// lets it be checked against geometry with closed form answers rather than by loading a structure and
// looking at it. void_analysis.cpp owns the window, the parameters, the threading and the overlay.
//
// The sections follow the pipeline:
//
//   porosity and accessible volume   the profile the field is summarized into, and its reductions
//   the distance field               grid, tiles, and the pass that fills a profile
//   surface topography               the height a probe of radius R finds, from above and below
//   channels                         connectivity through z, from a materialized field

#include <core/md_array.h>
#include <core/md_vec_math.h>

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

struct md_allocator_i;
struct md_spatial_acc_t;
struct md_unitcell_t;
struct md_grid_t;

// =================================================================================================
// Porosity and accessible volume, as reductions of the distance field
// =================================================================================================
//
// The field stage produces, for every voxel, the additively weighted distance to the nearest bead
// surface, d = min_i(|p - c_i| - R_i): the radius of the largest probe sphere centred there. Every
// quantity here is the same reduction of it,
//
//     V(R, z) = Vol[ d > R ]   restricted to a slab of z,
//
// evaluated at different R. Porosity is the R = 0 case - the void fraction is the volume a point
// sized probe can occupy - so it is not a separate computation but the left endpoint of the
// accessible volume curve. Keeping them one function is what keeps them consistent: a porosity and
// an accessible volume derived from different passes can disagree, and then neither is trustworthy.
//
// The field is summarized while its tiles are still in cache, as a histogram of d resolved along z.
// A 2.4e9 voxel field becomes a few hundred kilobytes and every quantity below is a sum over bins.
// The price is a discretization in R of one bin width, which is below the voxel spacing and so below
// the error the voxelization itself carries; the price of not doing it is keeping 9.6 GB resident to
// answer a question whose answer is a scalar.
//
// Two things this deliberately does not claim:
//
// - This is the volume accessible to the probe CENTRE. The volume the probe body sweeps out is the
//   dilation of that set by R and is larger - that is the Gelb-Gubbins pore volume, and it needs the
//   field itself rather than a histogram of it, because a voxel's value then depends on geometry up
//   to 2R away.
// - Accessible does not mean reachable. Nothing here knows whether a cavity connects to the outside.
//   That is connectivity, and it lives in the channels section below. A closed pore is counted here
//   and is invisible to infiltration, so the two numbers are expected to differ and their difference
//   is the closed porosity.

// A histogram of the distance field, resolved along z.
//
// Bin b covers [b * bin_width, (b + 1) * bin_width). Only void voxels - those with d > 0 - are
// binned; a voxel inside the solid is counted in solid[] instead, and total[] counts every voxel the
// slab considered at all. Every voxel therefore lands in exactly one of the two, which is what makes
// porosity = (total - solid) / total and the R = 0 accessible fraction the same number.
//
// The arrays are owned by the caller and nothing here allocates.
typedef struct void_profile_t {
    const uint64_t* hist;       // [slab * num_bins + bin], void voxels by distance
    const uint64_t* solid;      // [slab], voxels with d <= 0
    const uint64_t* total;      // [slab], voxels considered

    uint32_t num_slabs;
    uint32_t num_bins;

    double bin_width;           // Distance bin width, world units
    double z_min;               // World z of the lower edge of slab 0
    double slab_height;         // World z height of one slab
    double voxel_volume;        // World volume of one voxel
} void_profile_t;

// The binning conventions, shared by whatever fills a profile and by everything that reads one.
// They live here rather than at the fill site so the two cannot drift: a histogram binned one way
// and read another is a class of bug that produces plausible numbers.

// Bin b covers [b * bin_width, (b + 1) * bin_width). A distance at or beyond the last edge lands in
// the last bin, which is where the field query's own max_dist clamp ends up. Only call this for
// d > 0 - a voxel inside the solid is not binned at all.
static inline uint32_t void_profile_bin_of(double d, double bin_width, uint32_t num_bins) {
    if (!(d > 0.0) || !(bin_width > 0.0) || num_bins == 0) return 0;
    const double u = d / bin_width;
    if (u >= (double)num_bins) return num_bins - 1;
    return (uint32_t)u;
}

// Voxel plane k of dim_z, assigned by its centre, so the slabs are the uniform division of the grid
// in z that z_lo and z_hi report whether or not num_slabs divides dim_z.
static inline uint32_t void_profile_slab_of(int k, int dim_z, uint32_t num_slabs) {
    if (dim_z <= 0 || num_slabs == 0) return 0;
    if (k < 0) return 0;
    const uint64_t s = ((uint64_t)(2 * k + 1) * (uint64_t)num_slabs) / ((uint64_t)2 * (uint64_t)dim_z);
    return (s >= (uint64_t)num_slabs) ? num_slabs - 1 : (uint32_t)s;
}

#ifdef __cplusplus
extern "C" {
#endif

bool void_profile_valid(const void_profile_t* prof);

// World z bounds of a slab, and the slab a world z falls in (clamped to the profile).
double   void_profile_z_lo(const void_profile_t* prof, uint32_t slab);
double   void_profile_z_hi(const void_profile_t* prof, uint32_t slab);
uint32_t void_profile_slab_at(const void_profile_t* prof, double z);

// Voxel counts over the half open slab range [slab_beg, slab_end).
uint64_t void_profile_num_total(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end);
uint64_t void_profile_num_solid(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end);

// Void voxels with d > r. The bin containing r is split linearly rather than counted whole, which
// makes this continuous and monotone in r instead of stepping by a whole bin at a time. Any r <= 0
// returns every void voxel, so porosity is this evaluated at zero.
double void_profile_count_above(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r);

// V(R) over the slab range in world volume units, and the same divided by the volume of the range.
double void_profile_accessible_volume  (const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r);
double void_profile_accessible_fraction(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r);

// Porosity: the void volume fraction, i.e. the accessible fraction at zero probe radius.
double void_profile_porosity(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end);

// Voxels per unit distance at r, averaged over a window of at least min_width. By the coarea formula
// and |grad d| = 1 almost everywhere, -dV/dr is the accessible surface area at r, so this is that
// derivative up to a factor of voxel_volume. A window narrower than the voxel spacing only splits
// the same samples into noisier buckets, which is what min_width is for.
double void_profile_density(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r, double min_width);

// Solid volume fraction of one slab, and the same smoothed over its immediate neighbours. The film
// extent below is read off the smoothed profile, and plotting both next to the answer is how the
// answer gets checked.
double void_profile_solid_fraction(const void_profile_t* prof, uint32_t slab);
double void_profile_solid_fraction_smooth(const void_profile_t* prof, uint32_t slab);

// The slab range a film occupies, by the half density convention: the outermost slabs at which the
// solid volume fraction is still at least frac of its interior value. For a symmetric surface
// profile frac = 0.5 places the boundary at the Gibbs dividing surface, which is what a thickness
// means when the surface is rough enough to have no single z.
//
// This exists because a film simulated with an open z axis sits in a box with vacuum above and below
// it, and a porosity averaged over that box measures the size of the box rather than the porosity of
// the film. Reporting one number over the whole grid is the mistake this is here to prevent.
//
// Returns false when the profile carries no solid at all, in which case the outputs are untouched.
bool void_profile_film_extent(const void_profile_t* prof, double frac, uint32_t* out_slab_beg, uint32_t* out_slab_end, double* out_interior_solid_fraction);

// Mass per slab, binned along the same uniform division of z the profile uses, so a density read off
// it and a porosity read off the histogram describe the same slab. A position below or above the
// profile's z range is wrapped into it when periodic_z is set and dropped otherwise - an atom outside
// an open axis is outside the volume the density is taken over, and folding it into the edge slab
// would invent a spike there. mass may be NULL, in which case every position counts as one, which
// is the number density. out_slab_mass is overwritten, [num_slabs]. Returns what was binned.
double void_profile_bin_mass(double* out_slab_mass, const void_profile_t* prof, const vec3_t* xyz, const float* mass, size_t count, bool periodic_z);

// Binned mass over the slab range divided by the volume the profile considered there, i.e. the same
// volume the porosity of that range is a fraction of. Units are whatever the mass and the profile
// were given in, per world volume.
double void_profile_mass_density(const void_profile_t* prof, const double* slab_mass, uint32_t slab_beg, uint32_t slab_end);

#ifdef __cplusplus
}
#endif

// =================================================================================================
// The distance field
// =================================================================================================
//
// For every voxel of a regular grid, the additively weighted distance to the nearest bead surface,
//
//     d(p) = min_i ( |p - c_i| - R_i ),
//
// which is the radius of the largest probe sphere that fits at p. The query itself is the weighted
// nearest query on md_spatial_acc, so the field is exact everywhere rather than banded, and
// periodicity is whatever the unit cell the structure was built with says it is.
//
// The grid is walked in tiles of VOID_FIELD_TILE_DIM^3 voxels, which is what the batched query is
// efficient at, and each tile is summarized into a profile accumulator while it is still in cache.
// Evaluation is serial. A range of tiles touches only its own accumulator and its own voxels of the
// field, so a caller parallelises by handing disjoint tile ranges to threads, each with its own
// accumulator, and merging afterwards - the result does not depend on how the tiles were split.
//
// Output contract: unsigned distance in the void, clamped to 0 inside the solid.

#define VOID_FIELD_TILE_DIM 8

typedef struct void_field_desc_t {
    const struct md_spatial_acc_t* acc;     // Built over the beads, one radius per bead
    const struct md_unitcell_t*    cell;    // The cell acc was built with, or NULL. Drives the in-cell mask.
    const struct md_grid_t*        grid;    // Axis aligned, see void_field_grid

    double   max_dist;      // Range of the query. A voxel with no bead within it reports max_dist.
    uint32_t num_slabs;     // z resolution of the profile, at most grid->dim[2]
    uint32_t num_bins;      // Distance bins over [0, max_dist]

    float*   field;         // Optional, md_grid_num_points(grid) entries, x fastest. Written for every voxel.
} void_field_desc_t;

// What the field is summarized into. The arrays are caller owned and sized as a void_profile_t's:
// hist [num_slabs * num_bins], solid [num_slabs], total [num_slabs].
typedef struct void_field_accum_t {
    uint64_t* hist;
    uint64_t* solid;
    uint64_t* total;
    uint64_t  num_clamped;  // Voxels with no bead within max_dist
    float     d_min;        // Over the voxels the statistics considered; d_min > d_max when there were none
    float     d_max;
} void_field_accum_t;

#ifdef __cplusplus
extern "C" {
#endif

// Axis aligned grid over the cartesian bounds of the unit cell, or over the points plus a margin of
// two voxels when there is no cell. The voxel count per axis is extent / spacing rounded, and the
// spacing is then re-derived from it, so the grid tiles the box exactly: a periodic field sampled
// with a spacing that does not divide the box would otherwise seam.
bool void_field_grid(struct md_grid_t* out_grid, const struct md_unitcell_t* cell, const vec3_t* xyz, size_t count, float spacing);

uint32_t void_field_num_tiles(const struct md_grid_t* grid);

void void_field_accum_reset(void_field_accum_t* accum, uint32_t num_slabs, uint32_t num_bins);
void void_field_accum_merge(void_field_accum_t* dst, const void_field_accum_t* src, uint32_t num_slabs, uint32_t num_bins);

// Evaluate the tiles [tile_beg, tile_end) into accum, and into desc->field when there is one.
void void_field_eval_tiles(void_field_accum_t* accum, const void_field_desc_t* desc, uint32_t tile_beg, uint32_t tile_end);

// A profile over an accumulator, with the geometry filled in from the description. It points into
// the accumulator's arrays and does not own them.
void_profile_t void_field_profile(const void_field_accum_t* accum, const void_field_desc_t* desc);

#ifdef __cplusplus
}
#endif

// =================================================================================================
// Surface topography
// =================================================================================================
//
// The height of the structure as a spherical probe of radius R, lowered straight down along z,
// finds it: for every (x, y) column, the lowest z the probe centre reaches from the top before
// d(p) = R, minus R, which is the height of the probe apex at contact. The same from below gives the
// lower face. This is exactly what an AFM tip of that radius records, and it is what "respects the
// probe radius" means here: the reported surface is the structure dilated by R and eroded back,
// so a gap narrower than 2R is bridged and a cavity the probe cannot enter is not part of the
// topography at all. At R = 0 it is the bead surface itself.
//
// No field is needed. d is a lower bound on the distance to the union of the beads, so it is
// 1-Lipschitz, and descending by d - R per step can never step past a contact - this is sphere
// tracing on the same weighted nearest query the field pass uses. Near a contact the step is held
// at no less than tol and the crossing interpolated, so a column grazing the flank of a bead lands
// on it rather than creeping towards it. A column costs a handful of
// queries rather than a voxel per plane, which is what lets this run on a grid the field could not
// be materialized for.
//
// Columns are walked in patches of VOID_FIELD_TILE_DIM^2 so that each batched query stays spatially
// compact. Evaluation is serial and a patch range touches only its own columns, so a caller
// parallelises by handing out disjoint patch ranges, as with the field.

typedef struct void_heightmap_desc_t {
    const struct md_spatial_acc_t* acc;     // Built over the beads, one radius per bead
    const struct md_grid_t*        grid;    // Columns at its xy voxel centres; the scan spans its z extent

    double probe_radius;    // R, world units, >= 0
    double max_dist;        // Range of the query, must exceed probe_radius
    double tol;             // Contact tolerance, world units. <= 0 picks 1% of the smallest spacing
} void_heightmap_desc_t;

// Summary of one height map over the columns which found a contact.
typedef struct void_heightmap_stats_t {
    uint64_t num_valid;     // Columns with a contact
    uint64_t num_open;      // Columns the probe fell straight through
    double   mean;
    double   min;
    double   max;
    double   rq;            // RMS deviation from the mean, the usual Rq roughness
    double   ra;            // Mean absolute deviation, Ra
} void_heightmap_stats_t;

#ifdef __cplusplus
extern "C" {
#endif

uint32_t void_heightmap_num_patches(const struct md_grid_t* grid);

// Evaluate the patches [patch_beg, patch_end). out_top and out_bot hold grid->dim[0] * grid->dim[1]
// heights, x fastest, in world z; either may be NULL. A column the probe passes all the way through
// is written as NaN in both - it is an open pore, not a surface at the bottom of the box.
void void_heightmap_eval_patches(float* out_top, float* out_bot, const void_heightmap_desc_t* desc, uint32_t patch_beg, uint32_t patch_end);

// Statistics over the finite entries of h. NaN entries are counted as open.
void void_heightmap_stats(void_heightmap_stats_t* out, const float* h, size_t count);

#ifdef __cplusplus
}
#endif

// =================================================================================================
// Channels through a structure along the sweep axis
// =================================================================================================
//
// Given a clearance field - for every voxel, the radius of the largest probe sphere centred there -
// a probe of radius r can occupy exactly the superlevel set {d >= r}. A channel through the
// structure is a connected component of that set which touches both faces of the sweep axis.
//
// Sweeping a plane from the far face towards the near one, the components of {d >= r} restricted to
// the part already swept can only appear or merge, never split, because the swept region only ever
// grows. That makes the sweep history a forest: a genuine tree rooted at the near face, whose
// branches are stretches of channel and whose leaves are where a branch first appears. It is the
// structure behind both the count and the diagram.
//
// The sweep axis is z (axis 2) and must be non periodic. x and y may be periodic.

#define CHANNEL_INVALID_INDEX 0xFFFFFFFFu

// A clearance field on a regular grid, indexed x fastest.
typedef struct channel_field_t {
    const float* data;
    int   dim[3];
    float spacing[3];
    float origin[3];
    bool  pbc[3];       // Periodicity per axis. pbc[2] must be false.
} channel_field_t;

// One branch of the sweep forest: a component over the z range in which it had its own identity,
// from where it appeared down to where it merged into its parent.
typedef struct channel_node_t {
    uint32_t parent;            // Merge node below this branch, CHANNEL_INVALID_INDEX for a root
    uint32_t first_child;
    uint32_t next_sibling;

    float    z_top;             // Where the branch appears
    float    z_bot;             // Where it merges into its parent, or the bottom of the grid
    float    max_clearance;     // Largest clearance anywhere in the branch
    float    min_clearance;     // Smallest clearance along its own centreline, i.e. its tightest point
    uint32_t num_voxels;

    bool     reaches_top;       // This branch or one below it touches the far face
    bool     reaches_bottom;    // Only a root can, and only if its component touches the near face

    // Representative centreline: one point per slab the branch spans, from its top down. xyz is the
    // world position of the widest voxel of the branch in that slab, w is the clearance there. It is
    // a representative route, not the exact widest geodesic - see channel_critical_radius for the
    // number which is exact.
    md_array(vec4_t) path;
} channel_node_t;

typedef struct channel_tree_t {
    md_array(channel_node_t) nodes;
    md_array(uint32_t) roots;

    uint32_t num_components;    // Connected components of {d >= r} in the whole grid
    uint32_t num_spanning;      // Of those, the ones touching both faces: the channels
    double   probe_radius;

    struct md_allocator_i* alloc;
} channel_tree_t;

// Percolation and reachability, from one pass in order of decreasing clearance.
//
// Insert voxels into a union-find from the widest outward. At the moment some component first
// touches both z faces, the clearance of the voxel just inserted IS the critical radius - it is the
// tightest point of the widest route, which is what a maximum-capacity path gives and what bisecting
// on "does anything span" only brackets. That voxel is also *where* the route is limited, which no
// amount of bisection recovers: a scalar r_c says a probe of that size gets through and nothing
// about what stops a larger one.
//
// Everything else falls out of the same order, at no extra cost:
//
// - A component touching a z face is reachable from outside, so open and closed porosity separate
//   without the second flood fill the design document still lists as owed. A cavity that fits the
//   probe but connects to nothing is counted by the accessible volume and not by this.
// - The lowest z any top-touching component reaches is how far a probe of that radius penetrates
//   from above, and the same from below. They approach each other as the probe shrinks and meet at
//   r_c, which is the percolation threshold arrived at from the other direction - and unlike a
//   single number, the pair of curves says how much of the film a probe too large to cross can
//   still get into.
//
// Counting connected components, which is what channel_sweep reports, is nearly useless on a real
// network: below r_c the void space is one component and the count is 1, above it 0. The curves
// here are what that step function was standing in for.
typedef struct channel_percolation_t {
    md_array(double)   radius;          // Ascending, from r_min to the widest clearance in the field
    md_array(double)   frac_void;       // Vol{d >= r} / V, the whole box
    md_array(double)   frac_open;       // Of the box: in a component which touches a z face
    md_array(double)   frac_spanning;   // Of the box: in a component which touches both
    md_array(double)   z_from_top;      // Lowest z reached from the top face; the top of the grid if none
    md_array(double)   z_from_bottom;   // Highest z reached from the bottom face
    md_array(uint32_t) num_components;
    md_array(uint32_t) num_spanning;

    bool   has_r_c;
    double r_c;                         // Clearance at which the two faces first connect
    float  throat[3];                   // Where: the voxel whose insertion connected them

    size_t num_active;                  // Voxels at or above r_min, which is what the run was sized by
    size_t bytes;                       // What it allocated, so the cost of a larger field is legible

    struct md_allocator_i* alloc;
} channel_percolation_t;

#ifdef __cplusplus
extern "C" {
#endif

void channel_tree_free(channel_tree_t* tree);
void channel_percolation_free(channel_percolation_t* perc);

// One pass. r_min bounds the memory - voxels below it never enter the structure - and is therefore
// also the lowest radius the curves reach. num_samples is the number of points on them.
//
// Costs 4 bytes per voxel of the field plus 17 per voxel at or above r_min, against 9 per active
// voxel for a single channel_sweep. It replaces a bisection of about a dozen of those plus one sweep
// per reported radius - measured at 6.5x faster than that pair on a 7.7 Mvoxel disordered network,
// and exact where they were bracketed.
bool channel_percolate(channel_percolation_t* out, const channel_field_t* field, double r_min, uint32_t num_samples, struct md_allocator_i* alloc);

// A route a probe of radius r can actually follow from the top face to the bottom: breadth first
// through {d >= r}, so it is the shortest such route in voxel steps, then straightened wherever the
// segment between two of its points stays inside the set.
//
// Unlike the representative centreline on a tree branch - the widest voxel per slab, which need not
// be connected to the widest voxel of the slab below - every point of this is reachable from the
// previous one at radius r. Trace it at r_c and it is the route the critical radius belongs to.
//
// out_path is xyz in world space with w the clearance there, wrapped into the box on a periodic
// axis. out_length is the length of the unwrapped route, so length / |z span| is the tortuosity.
// Returns false when nothing gets through at r.
bool channel_trace_path(md_array(vec4_t)* out_path, double* out_length, const channel_field_t* field, double r, struct md_allocator_i* alloc);

// Sweep the field at one probe radius. With want_tree false only the counts are produced, which is
// what the radius sweeps below use.
bool channel_sweep(channel_tree_t* out, const channel_field_t* field, double probe_radius, bool want_tree, struct md_allocator_i* alloc);

// Number of channels at each supplied radius.
void channel_spanning_counts(uint32_t* out_counts, const double* radii, size_t num_radii, const channel_field_t* field, struct md_allocator_i* alloc);

// Assign a dendrogram column to every branch belonging to a channel which gets all the way through.
// Leaves take consecutive columns in traversal order and a merge node sits at the mean of its
// children, so a branch is drawn above the span of what feeds it. Branches of any other component
// are left at -1: a blind pore would only crowd the diagram.
//
// out_slot must have room for md_array_size(tree->nodes) entries. Returns the number of columns used.
//
// The traversal carries its own stack. A merge tree has one node per pair of branches that join, so
// a real network produces a chain thousands of levels deep and recursing over it overflows.
float channel_tree_layout(float* out_slot, const channel_tree_t* tree, struct md_allocator_i* alloc);

// The largest probe radius which still gets through: bisected until the bracket is below tol.
//
// channel_percolate answers this exactly and in one pass, and is what the component uses. This stays
// because it asks a completely different question of the field - a threshold test repeated, rather
// than an ordering - so the two agreeing on r_c is evidence about the answer rather than about one
// implementation. The tests cross check them.
// This is a property of the tightest throat along the best route, not of any cavity along it, and it
// is the number to quote - the branch clearances above are read off a representative centreline.
// Returns a value below r_lo when nothing gets through even at r_lo.
double channel_critical_radius(const channel_field_t* field, double r_lo, double r_hi, double tol, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
