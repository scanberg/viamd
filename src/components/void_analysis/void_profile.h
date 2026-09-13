#pragma once

// Porosity and accessible volume, as reductions of the distance field.
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
//   That is connectivity, and it lives in channels.h. A closed pore is counted here and is invisible
//   to infiltration, so the two numbers are expected to differ and their difference is the closed
//   porosity.

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

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

#ifdef __cplusplus
}
#endif
