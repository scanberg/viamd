#include "void_profile.h"

#include <math.h>

namespace {

// Clamp a half open slab range to the profile. Returns false when nothing is left of it, which is
// the one case every reduction below has to short circuit rather than divide by.
bool clamp_range(const void_profile_t* prof, uint32_t* beg, uint32_t* end) {
    if (!void_profile_valid(prof)) return false;
    uint32_t b = *beg;
    uint32_t e = *end;
    if (b > prof->num_slabs) b = prof->num_slabs;
    if (e > prof->num_slabs) e = prof->num_slabs;
    if (b >= e) return false;
    *beg = b;
    *end = e;
    return true;
}

}  // namespace

bool void_profile_valid(const void_profile_t* prof) {
    return prof && prof->hist && prof->solid && prof->total &&
           prof->num_slabs > 0 && prof->num_bins > 0 && prof->bin_width > 0.0;
}

double void_profile_z_lo(const void_profile_t* prof, uint32_t slab) {
    if (!prof) return 0.0;
    return prof->z_min + (double)slab * prof->slab_height;
}

double void_profile_z_hi(const void_profile_t* prof, uint32_t slab) {
    if (!prof) return 0.0;
    return prof->z_min + (double)(slab + 1) * prof->slab_height;
}

uint32_t void_profile_slab_at(const void_profile_t* prof, double z) {
    if (!prof || prof->num_slabs == 0 || !(prof->slab_height > 0.0)) return 0;
    const double s = floor((z - prof->z_min) / prof->slab_height);
    if (s <= 0.0) return 0;
    if (s >= (double)(prof->num_slabs - 1)) return prof->num_slabs - 1;
    return (uint32_t)s;
}

uint64_t void_profile_num_total(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end) {
    if (!clamp_range(prof, &slab_beg, &slab_end)) return 0;
    uint64_t n = 0;
    for (uint32_t s = slab_beg; s < slab_end; ++s) n += prof->total[s];
    return n;
}

uint64_t void_profile_num_solid(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end) {
    if (!clamp_range(prof, &slab_beg, &slab_end)) return 0;
    uint64_t n = 0;
    for (uint32_t s = slab_beg; s < slab_end; ++s) n += prof->solid[s];
    return n;
}

double void_profile_count_above(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r) {
    if (!clamp_range(prof, &slab_beg, &slab_end)) return 0.0;

    // Where r falls in the binning. Everything at or below bin b0 except the tail of b0 itself is
    // excluded; the tail is what the linear split recovers, on the assumption that d is uniform
    // within a bin. At a bin width well below the voxel spacing that assumption costs less than the
    // voxelization already does.
    const double u = (r > 0.0) ? r / prof->bin_width : 0.0;
    if (u >= (double)prof->num_bins) return 0.0;

    const uint32_t b0   = (uint32_t)u;
    const double   frac = u - (double)b0;

    double count = 0.0;
    for (uint32_t s = slab_beg; s < slab_end; ++s) {
        const uint64_t* h = prof->hist + (size_t)s * prof->num_bins;
        uint64_t whole = 0;
        for (uint32_t b = b0 + 1; b < prof->num_bins; ++b) whole += h[b];
        count += (double)whole + (1.0 - frac) * (double)h[b0];
    }
    return count;
}

double void_profile_accessible_volume(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r) {
    return void_profile_count_above(prof, slab_beg, slab_end, r) * (prof ? prof->voxel_volume : 0.0);
}

double void_profile_accessible_fraction(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r) {
    const uint64_t total = void_profile_num_total(prof, slab_beg, slab_end);
    if (total == 0) return 0.0;
    return void_profile_count_above(prof, slab_beg, slab_end, r) / (double)total;
}

double void_profile_porosity(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end) {
    return void_profile_accessible_fraction(prof, slab_beg, slab_end, 0.0);
}

double void_profile_density(const void_profile_t* prof, uint32_t slab_beg, uint32_t slab_end, double r, double min_width) {
    if (!void_profile_valid(prof)) return 0.0;
    const double w  = (min_width > prof->bin_width) ? min_width : prof->bin_width;
    const double hi = r + 0.5 * w;
    const double lo = (r - 0.5 * w > 0.0) ? r - 0.5 * w : 0.0;
    if (!(hi > lo)) return 0.0;
    const double n_lo = void_profile_count_above(prof, slab_beg, slab_end, lo);
    const double n_hi = void_profile_count_above(prof, slab_beg, slab_end, hi);
    return (n_lo - n_hi) / (hi - lo);
}

double void_profile_solid_fraction(const void_profile_t* prof, uint32_t slab) {
    if (!void_profile_valid(prof) || slab >= prof->num_slabs) return 0.0;
    const uint64_t t = prof->total[slab];
    return t ? (double)prof->solid[slab] / (double)t : 0.0;
}

double void_profile_solid_fraction_smooth(const void_profile_t* prof, uint32_t slab) {
    if (!void_profile_valid(prof) || slab >= prof->num_slabs) return 0.0;
    // Edge clamped [1 2 1]/4. A single noisy slab should not be able to move a surface by its own
    // width, and at the free surface of a rough film exactly one slab is usually the noisy one.
    const uint32_t lo = (slab > 0) ? slab - 1 : 0;
    const uint32_t hi = (slab + 1 < prof->num_slabs) ? slab + 1 : prof->num_slabs - 1;
    return 0.25 * void_profile_solid_fraction(prof, lo)
         + 0.50 * void_profile_solid_fraction(prof, slab)
         + 0.25 * void_profile_solid_fraction(prof, hi);
}

bool void_profile_film_extent(const void_profile_t* prof, double frac, uint32_t* out_slab_beg, uint32_t* out_slab_end, double* out_interior_solid_fraction) {
    if (!void_profile_valid(prof)) return false;

    // The interior value is taken as the peak of the smoothed profile rather than an average over
    // some assumed interior, because which slabs are interior is precisely what is being decided.
    double peak = 0.0;
    for (uint32_t s = 0; s < prof->num_slabs; ++s) {
        const double v = void_profile_solid_fraction_smooth(prof, s);
        if (v > peak) peak = v;
    }
    if (out_interior_solid_fraction) *out_interior_solid_fraction = peak;
    if (!(peak > 0.0)) return false;

    const double thr = frac * peak;

    // Outermost crossings, not the first contiguous run: an internal void large enough to drop the
    // solid fraction below the threshold is part of the film, not a gap between two films.
    uint32_t beg = prof->num_slabs;
    uint32_t end = 0;
    for (uint32_t s = 0; s < prof->num_slabs; ++s) {
        if (void_profile_solid_fraction_smooth(prof, s) >= thr) {
            if (s < beg) beg = s;
            end = s + 1;
        }
    }
    if (beg >= end) return false;

    if (out_slab_beg) *out_slab_beg = beg;
    if (out_slab_end) *out_slab_end = end;
    return true;
}
