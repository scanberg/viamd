#include "void_analysis_core.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_grid.h>
#include <core/md_hash.h>
#include <core/md_simd.h>
#include <core/md_spatial_acc.h>
#include <md_types.h>
#include <md_unitcell.h>

#include <float.h>
#include <math.h>
#include <string.h>

#include <algorithm>

// =================================================================================================
// Porosity and accessible volume
// =================================================================================================

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
    const double z = prof->z_min + (double)(slab + 1) * prof->slab_height;
    return (z < prof->z_max) ? z : prof->z_max;
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

// =================================================================================================
// The weighted nearest query
// =================================================================================================

namespace {

// Points are searched in batches, and within a batch in groups of consecutive points, each with its own
// box and its own worst distance so far. A bead is tested against the points of a group only when a test
// against the group box says it could improve on one of them. Smaller groups cull more finely but pay the
// box test per group for every bead; on a fibril network 128 measured best, 1.3 to 1.6 times faster than 32.
#ifndef VOID_NEAREST_GROUP
#define VOID_NEAREST_GROUP 128
#endif
// The first box reaches this far past the largest radius. It grows from there by doubling, so this only
// trades the work of an oversized first box against the rounds a small one takes in the open.
#ifndef VOID_NEAREST_SEED
#define VOID_NEAREST_SEED 2.0
#endif
constexpr int NEAREST_BATCH      = 512;
constexpr int NEAREST_GROUP      = VOID_NEAREST_GROUP;
constexpr int NEAREST_MAX_GROUPS = NEAREST_BATCH / NEAREST_GROUP;
static_assert(NEAREST_GROUP % 8 == 0 && NEAREST_BATCH % NEAREST_GROUP == 0, "groups are whole vectors");

struct nearest_group_t {
    float lo[3], hi[3];     // Box of the group's points
    float worst;            // Largest best distance in the group
    int   beg, end;         // Lanes, end rounded up to whole vectors
};

struct nearest_batch_t {
    // Points relative to ref, with the best distance and bead found for each. A group's padding lanes
    // repeat its last point, so they never change the worst distance of the group.
    alignas(32) float    px[NEAREST_BATCH];
    alignas(32) float    py[NEAREST_BATCH];
    alignas(32) float    pz[NEAREST_BATCH];
    alignas(32) float    best[NEAREST_BATCH];
    alignas(32) uint32_t bead[NEAREST_BATCH];

    nearest_group_t group[NEAREST_MAX_GROUPS];
    int num_groups;

    const float* radii;
    double shift[3];        // Takes a position the query reports to the frame of the batch: image offset - ref

    // When the search box no longer fits in the period, a bead is moved to its image nearest ref and tried
    // there and at the neighbouring images as well
    bool   images;
    int    num_images;
    float  image[27][3];
    double A[3][3];         // [col][row], as md_spatial_acc_t has them
    double I[3][3];
    bool   pbc[3];
};

// Try one bead against every group it could improve
inline void nearest_try(nearest_batch_t* b, float ex, float ey, float ez, float r, uint32_t idx) {
    for (int g = 0; g < b->num_groups; ++g) {
        nearest_group_t* grp = &b->group[g];

        // The bead improves on a point exactly when |p - c| < best + r. Nothing in the group can be nearer than
        // the group box is.
        const float t = grp->worst + r;
        if (!(t > 0.0f)) continue;
        const float bx = MAX(0.0f, MAX(grp->lo[0] - ex, ex - grp->hi[0]));
        const float by = MAX(0.0f, MAX(grp->lo[1] - ey, ey - grp->hi[1]));
        const float bz = MAX(0.0f, MAX(grp->lo[2] - ez, ez - grp->hi[2]));
        if (bx * bx + by * by + bz * bz >= t * t) continue;

        const md_256  vx = md_mm256_set1_ps(ex);
        const md_256  vy = md_mm256_set1_ps(ey);
        const md_256  vz = md_mm256_set1_ps(ez);
        const md_256  vr = md_mm256_set1_ps(r);
        const md_256  vi = md_mm256_castsi256_ps(md_mm256_set1_epi32((int)idx));
        const md_256  zero = md_mm256_setzero_ps();
        md_256 worst = md_mm256_set1_ps(-FLT_MAX);
        bool hit = false;

        for (int i = grp->beg; i < grp->end; i += 8) {
            const md_256 dx = md_mm256_sub_ps(md_mm256_load_ps(b->px + i), vx);
            const md_256 dy = md_mm256_sub_ps(md_mm256_load_ps(b->py + i), vy);
            const md_256 dz = md_mm256_sub_ps(md_mm256_load_ps(b->pz + i), vz);
            const md_256 d2 = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(dx, dx), md_mm256_mul_ps(dy, dy)), md_mm256_mul_ps(dz, dz));

            // Squared, which keeps the square root off the path of the beads which improve nothing
            md_256 bd = md_mm256_load_ps(b->best + i);
            const md_256 tt   = md_mm256_add_ps(bd, vr);
            const md_256 mask = md_mm256_and_ps(md_mm256_cmplt_ps(d2, md_mm256_mul_ps(tt, tt)), md_mm256_cmpgt_ps(tt, zero));
            if (md_mm256_movemask_ps(mask)) {
                const md_256 d  = md_mm256_sub_ps(md_mm256_sqrt_ps(d2), vr);
                const md_256 bi = md_mm256_load_ps((const float*)(b->bead + i));
                bd = md_mm256_blendv_ps(bd, d, mask);
                md_mm256_store_ps(b->best + i, bd);
                md_mm256_store_ps((float*)(b->bead + i), md_mm256_blendv_ps(bi, vi, mask));
                hit = true;
            }
            worst = md_mm256_max_ps(worst, bd);
        }
        if (hit) grp->worst = md_mm256_reduce_max_ps(worst);
    }
}

void nearest_callback(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num, void* user) {
    nearest_batch_t* b = (nearest_batch_t*)user;
    for (size_t k = 0; k < num; ++k) {
        const float r = b->radii[idx[k]];
        double e[3] = {
            (double)x[k] + b->shift[0],
            (double)y[k] + b->shift[1],
            (double)z[k] + b->shift[2],
        };
        if (!b->images) {
            nearest_try(b, (float)e[0], (float)e[1], (float)e[2], r, idx[k]);
            continue;
        }

        // The image nearest ref. Rounding the fractional coordinate is that image for an orthorhombic cell and
        // near enough to it for a triclinic one that the neighbours tried next include it.
        double n[3];
        for (int a = 0; a < 3; ++a) {
            n[a] = b->pbc[a] ? round(b->I[0][a] * e[0] + b->I[1][a] * e[1] + b->I[2][a] * e[2]) : 0.0;
        }
        for (int a = 0; a < 3; ++a) {
            e[a] -= b->A[0][a] * n[0] + b->A[1][a] * n[1] + b->A[2][a] * n[2];
        }
        for (int m = 0; m < b->num_images; ++m) {
            nearest_try(b, (float)e[0] + b->image[m][0], (float)e[1] + b->image[m][1], (float)e[2] + b->image[m][2], r, idx[k]);
        }
    }
}

// Every bead in the box [lo, hi], relative to ref
void nearest_box(nearest_batch_t* b, const md_spatial_acc_t* acc, const double ref[3], const double lo[3], const double hi[3]) {
    double cen[3], rad[3];
    for (int a = 0; a < 3; ++a) {
        if (!(hi[a] > lo[a])) return;
        cen[a] = ref[a] + 0.5 * (lo[a] + hi[a]);
        rad[a] = 0.5 * (hi[a] - lo[a]);
    }
    // The query reports positions in the image of its own centre, which is cen folded into the cell
    double img[3];
    md_spatial_acc_aabb_query_center(img, acc, cen);
    for (int a = 0; a < 3; ++a) {
        b->shift[a] = (cen[a] - img[a]) - ref[a];
    }
    md_spatial_acc_for_each_point_in_aabb(acc, cen, rad, nearest_callback, b);
}

// Every bead in the box between two nested boxes, as the six slabs it decomposes into
void nearest_shell(nearest_batch_t* b, const md_spatial_acc_t* acc, const double ref[3], const double in_lo[3], const double in_hi[3], const double out_lo[3], const double out_hi[3]) {
    for (int a = 0; a < 3; ++a) {
        // Full extent along the axes before a, the inner one along the axes after it
        double lo[3], hi[3];
        for (int c = 0; c < 3; ++c) {
            lo[c] = (c < a) ? out_lo[c] : in_lo[c];
            hi[c] = (c < a) ? out_hi[c] : in_hi[c];
        }
        lo[a] = out_lo[a]; hi[a] = in_lo[a];
        nearest_box(b, acc, ref, lo, hi);
        lo[a] = in_hi[a];  hi[a] = out_hi[a];
        nearest_box(b, acc, ref, lo, hi);
    }
}

}  // namespace

void_beads_t void_beads(const md_spatial_acc_t* acc, const float* radii, size_t num_radii) {
    void_beads_t beads = {};
    beads.acc   = acc;
    beads.radii = radii;
    float max_radius = 0.0f;
    for (size_t i = 0; i < num_radii; ++i) {
        max_radius = MAX(max_radius, radii[i]);
    }
    beads.max_radius = max_radius;
    return beads;
}

void void_beads_nearest(const void_beads_t* beads, const float* x, const float* y, const float* z, size_t count,
                        double max_dist, uint32_t* out_idx, float* out_dist) {
    if (count == 0) return;
    const float d_max = (float)max_dist;

    if (!beads || !beads->acc || !beads->radii || beads->acc->num_elems == 0) {
        for (size_t i = 0; i < count; ++i) {
            if (out_idx)  out_idx[i]  = VOID_BEAD_NONE;
            if (out_dist) out_dist[i] = d_max;
        }
        return;
    }

    const md_spatial_acc_t* acc = beads->acc;
    const double R = MAX(0.0, (double)beads->max_radius);

    // Nothing further than this from a point can be below max_dist
    const double e_max = max_dist + R;

    nearest_batch_t b;
    b.radii  = beads->radii;
    b.images = false;

    // The periodic axes, and how large a search box may get before it holds a bead twice. Orthorhombic: less
    // than a period along each periodic axis. Triclinic: a box inside a ball of less than half the smallest
    // distance between lattice planes, which no two images fit in.
    const bool tri = (acc->flags & MD_UNITCELL_TRICLINIC) != 0;
    bool any_pbc = false;
    double max_half[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
    double max_ball = DBL_MAX;
    for (int a = 0; a < 3; ++a) {
        b.pbc[a] = (acc->flags & (MD_UNITCELL_PBC_X << a)) != 0;
        any_pbc |= b.pbc[a];
        for (int c = 0; c < 3; ++c) {
            b.A[a][c] = (double)acc->A[a][c];
            b.I[a][c] = (double)acc->I[a][c];
        }
    }
    if (any_pbc) {
        const double* va = b.A[0];
        const double* vb = b.A[1];
        const double* vc = b.A[2];
        const double bc[3] = { vb[1] * vc[2] - vb[2] * vc[1], vb[2] * vc[0] - vb[0] * vc[2], vb[0] * vc[1] - vb[1] * vc[0] };
        const double ca[3] = { vc[1] * va[2] - vc[2] * va[1], vc[2] * va[0] - vc[0] * va[2], vc[0] * va[1] - vc[1] * va[0] };
        const double ab[3] = { va[1] * vb[2] - va[2] * vb[1], va[2] * vb[0] - va[0] * vb[2], va[0] * vb[1] - va[1] * vb[0] };
        const double vol = fabs(va[0] * bc[0] + va[1] * bc[1] + va[2] * bc[2]);
        const double* nrm[3] = { bc, ca, ab };
        for (int a = 0; a < 3; ++a) {
            if (!b.pbc[a]) continue;
            const double len   = sqrt(nrm[a][0] * nrm[a][0] + nrm[a][1] * nrm[a][1] + nrm[a][2] * nrm[a][2]);
            const double width = (len > 0.0) ? vol / len : 0.0;
            if (tri) {
                max_ball = MIN(max_ball, 0.5 * width);
            } else {
                max_half[a] = 0.5 * width;
            }
        }

        b.num_images = 0;
        for (int k = -1; k <= 1; ++k) {
            if (k && !b.pbc[2]) continue;
            for (int j = -1; j <= 1; ++j) {
                if (j && !b.pbc[1]) continue;
                for (int i = -1; i <= 1; ++i) {
                    if (i && !b.pbc[0]) continue;
                    float* t = b.image[b.num_images++];
                    for (int c = 0; c < 3; ++c) {
                        t[c] = (float)(b.A[0][c] * i + b.A[1][c] * j + b.A[2][c] * k);
                    }
                }
            }
        }
    }

    for (size_t base = 0; base < count; base += NEAREST_BATCH) {
        const int n = (int)MIN((size_t)NEAREST_BATCH, count - base);

        // Points relative to the centre of the batch, which keeps the distances accurate however far from the
        // origin the batch is
        double ref[3];
        {
            double lo[3] = {  DBL_MAX,  DBL_MAX,  DBL_MAX };
            double hi[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
            for (int i = 0; i < n; ++i) {
                const double p[3] = { (double)x[base + i], (double)y[base + i], (double)z[base + i] };
                for (int a = 0; a < 3; ++a) {
                    lo[a] = MIN(lo[a], p[a]);
                    hi[a] = MAX(hi[a], p[a]);
                }
            }
            for (int a = 0; a < 3; ++a) ref[a] = 0.5 * (lo[a] + hi[a]);
        }

        double ulo[3] = {  DBL_MAX,  DBL_MAX,  DBL_MAX };
        double uhi[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
        b.num_groups = 0;
        for (int beg = 0; beg < n; beg += NEAREST_GROUP) {
            const int end = MIN(beg + NEAREST_GROUP, n);
            nearest_group_t* grp = &b.group[b.num_groups++];
            grp->beg   = beg;
            grp->end   = (int)ALIGN_TO(end, 8);
            grp->worst = d_max;
            for (int a = 0; a < 3; ++a) {
                grp->lo[a] =  FLT_MAX;
                grp->hi[a] = -FLT_MAX;
            }
            for (int i = beg; i < grp->end; ++i) {
                const size_t s = base + (size_t)MIN(i, end - 1);
                const float p[3] = {
                    (float)((double)x[s] - ref[0]),
                    (float)((double)y[s] - ref[1]),
                    (float)((double)z[s] - ref[2]),
                };
                b.px[i]   = p[0];
                b.py[i]   = p[1];
                b.pz[i]   = p[2];
                b.best[i] = d_max;
                b.bead[i] = VOID_BEAD_NONE;
                for (int a = 0; a < 3; ++a) {
                    grp->lo[a] = MIN(grp->lo[a], p[a]);
                    grp->hi[a] = MAX(grp->hi[a], p[a]);
                }
            }
            for (int a = 0; a < 3; ++a) {
                ulo[a] = MIN(ulo[a], (double)grp->lo[a]);
                uhi[a] = MAX(uhi[a], (double)grp->hi[a]);
            }
        }

        // Grow the box around the batch until no bead outside it can improve on a point. A bead outside the box
        // is at least as far from a point of a group as the group box is from the faces of the search box, and
        // its weighted distance is that less its radius.
        double e = MIN(e_max, R + VOID_NEAREST_SEED);
        double e_prev = -1.0;
        for (;;) {
            bool fits = true;
            {
                double ball = 0.0;
                for (int a = 0; a < 3; ++a) {
                    const double h = 0.5 * (uhi[a] - ulo[a]) + e;
                    if (h >= max_half[a]) fits = false;
                    ball += h * h;
                }
                if (sqrt(ball) >= max_ball) fits = false;
            }

            if (!fits) {
                // Everything, each bead at every image which can matter. Final, whatever was found before.
                double reach = 0.0;
                for (int a = 0; a < 3; ++a) {
                    reach += sqrt(b.A[a][0] * b.A[a][0] + b.A[a][1] * b.A[a][1] + b.A[a][2] * b.A[a][2]);
                    reach += fabs(ref[a] - (double)acc->origin[a]);
                }
                const double lo[3] = { -reach - 1.0, -reach - 1.0, -reach - 1.0 };
                const double hi[3] = {  reach + 1.0,  reach + 1.0,  reach + 1.0 };
                b.images = true;
                nearest_box(&b, acc, ref, lo, hi);
                b.images = false;
                break;
            }

            const double out_lo[3] = { ulo[0] - e, ulo[1] - e, ulo[2] - e };
            const double out_hi[3] = { uhi[0] + e, uhi[1] + e, uhi[2] + e };
            if (e_prev < 0.0) {
                nearest_box(&b, acc, ref, out_lo, out_hi);
            } else {
                const double in_lo[3] = { ulo[0] - e_prev, ulo[1] - e_prev, ulo[2] - e_prev };
                const double in_hi[3] = { uhi[0] + e_prev, uhi[1] + e_prev, uhi[2] + e_prev };
                nearest_shell(&b, acc, ref, in_lo, in_hi, out_lo, out_hi);
            }

            double need = -DBL_MAX;
            bool done = true;
            for (int g = 0; g < b.num_groups; ++g) {
                const nearest_group_t* grp = &b.group[g];
                double margin = DBL_MAX;
                for (int a = 0; a < 3; ++a) {
                    margin = MIN(margin, MIN(out_hi[a] - (double)grp->hi[a], (double)grp->lo[a] - out_lo[a]));
                }
                const double reach = (double)grp->worst + R;
                if (reach > margin) done = false;
                need = MAX(need, reach);
            }
            if (done || e >= e_max) break;

            e_prev = e;
            e = MIN(e_max, MIN(2.0 * e, need));
            if (!(e > e_prev)) break;
        }

        for (int i = 0; i < n; ++i) {
            if (out_idx)  out_idx[base + i]  = b.bead[i];
            if (out_dist) out_dist[base + i] = b.best[i];
        }
    }
}

// =================================================================================================
// The distance field
// =================================================================================================

double void_profile_bin_mass(double* out_slab_mass, const void_profile_t* prof, const vec3_t* xyz, const float* mass, size_t count, bool periodic_z) {
    if (!out_slab_mass || !prof || prof->num_slabs == 0 || !(prof->slab_height > 0.0)) return 0.0;
    memset(out_slab_mass, 0, prof->num_slabs * sizeof(double));
    if (!xyz) return 0.0;

    const double H = prof->z_max - prof->z_min;
    if (!(H > 0.0)) return 0.0;
    double total = 0.0;
    for (size_t i = 0; i < count; ++i) {
        double u = (double)xyz[i].z - prof->z_min;
        if (periodic_z) {
            u -= H * floor(u / H);
        } else if (u < 0.0 || u >= H) {
            continue;
        }
        uint32_t s = (uint32_t)(u / prof->slab_height);
        if (s >= prof->num_slabs) s = prof->num_slabs - 1;
        const double m = mass ? (double)mass[i] : 1.0;
        out_slab_mass[s] += m;
        total += m;
    }
    return total;
}

double void_profile_mass_density(const void_profile_t* prof, const double* slab_mass, uint32_t slab_beg, uint32_t slab_end) {
    if (!slab_mass || !clamp_range(prof, &slab_beg, &slab_end)) return 0.0;
    double m = 0.0;
    for (uint32_t s = slab_beg; s < slab_end; ++s) m += slab_mass[s];
    const double v = (double)void_profile_num_total(prof, slab_beg, slab_end) * prof->voxel_volume;
    return (v > 0.0) ? m / v : 0.0;
}

bool void_field_grid(md_grid_t* out_grid, const md_unitcell_t* cell, const vec3_t* xyz, size_t count, float spacing) {
    if (!out_grid || !(spacing > 0.0f)) return false;

    const uint32_t flags = cell ? md_unitcell_flags(cell) : (uint32_t)MD_UNITCELL_NONE;

    vec3_t origin = {0, 0, 0};
    vec3_t extent = {0, 0, 0};

    if (flags != MD_UNITCELL_NONE) {
        // Cartesian bounds of the cell parallelepiped: the sum of the positive components of the basis vectors
        double A[3][3];
        md_unitcell_A_extract_double(A, cell);
        for (int r = 0; r < 3; ++r) {
            double lo = 0.0, hi = 0.0;
            for (int c = 0; c < 3; ++c) {
                const double v = A[c][r];
                if (v < 0.0) lo += v; else hi += v;
            }
            origin.elem[r] = (float)lo;
            extent.elem[r] = (float)(hi - lo);
        }
    } else {
        if (count == 0 || !xyz) return false;
        vec3_t aabb_min = { FLT_MAX,  FLT_MAX,  FLT_MAX};
        vec3_t aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (size_t i = 0; i < count; ++i) {
            const vec3_t p = xyz[i];
            aabb_min = vec3_min(aabb_min, p);
            aabb_max = vec3_max(aabb_max, p);
        }
        const float margin = 2.0f * spacing;
        for (int a = 0; a < 3; ++a) {
            origin.elem[a] = aabb_min.elem[a] - margin;
            extent.elem[a] = (aabb_max.elem[a] - aabb_min.elem[a]) + 2.0f * margin;
        }
    }

    for (int a = 0; a < 3; ++a) {
        const int dim = (int)round((double)extent.elem[a] / (double)spacing);
        out_grid->dim[a] = MAX(1, dim);
        out_grid->spacing.elem[a] = extent.elem[a] / (float)out_grid->dim[a];
    }
    out_grid->origin = origin;
    out_grid->orientation = mat3_ident();

    return md_grid_num_points(out_grid) > 0;
}

uint32_t void_field_num_tiles(const md_grid_t* grid) {
    if (!grid) return 0;
    uint32_t n = 1;
    for (int a = 0; a < 3; ++a) {
        if (grid->dim[a] <= 0) return 0;
        n *= (uint32_t)((grid->dim[a] + VOID_FIELD_TILE_DIM - 1) / VOID_FIELD_TILE_DIM);
    }
    return n;
}

void void_field_accum_reset(void_field_accum_t* accum, uint32_t num_slabs, uint32_t num_bins) {
    if (!accum) return;
    if (accum->hist)  MEMSET(accum->hist,  0, (size_t)num_slabs * num_bins * sizeof(uint64_t));
    if (accum->solid) MEMSET(accum->solid, 0, (size_t)num_slabs * sizeof(uint64_t));
    if (accum->total) MEMSET(accum->total, 0, (size_t)num_slabs * sizeof(uint64_t));
    accum->num_clamped = 0;
    accum->d_min =  FLT_MAX;
    accum->d_max = -FLT_MAX;
}

void void_field_accum_merge(void_field_accum_t* dst, const void_field_accum_t* src, uint32_t num_slabs, uint32_t num_bins) {
    if (!dst || !src) return;
    const size_t stride = (size_t)num_slabs * num_bins;
    for (size_t i = 0; i < stride; ++i) dst->hist[i] += src->hist[i];
    for (uint32_t s = 0; s < num_slabs; ++s) {
        dst->solid[s] += src->solid[s];
        dst->total[s] += src->total[s];
    }
    dst->num_clamped += src->num_clamped;
    dst->d_min = MIN(dst->d_min, src->d_min);
    dst->d_max = MAX(dst->d_max, src->d_max);
}

void void_field_eval_tiles(void_field_accum_t* accum, const void_field_desc_t* desc, uint32_t tile_beg, uint32_t tile_end) {
    if (!accum || !desc || !desc->beads.acc || !desc->grid) return;
    if (desc->planes_per_slab == 0 || desc->num_bins == 0) return;

    const md_grid_t& grid = *desc->grid;
    const int TD = VOID_FIELD_TILE_DIM;
    const int tiles[3] = {
        (grid.dim[0] + TD - 1) / TD,
        (grid.dim[1] + TD - 1) / TD,
        (grid.dim[2] + TD - 1) / TD,
    };
    tile_end = MIN(tile_end, void_field_num_tiles(desc->grid));

    const uint32_t pps       = desc->planes_per_slab;
    const uint32_t num_bins  = desc->num_bins;
    const float    max_dist  = (float)desc->max_dist;

    // Bin edges over [0, max_dist], which is exactly the range the query resolves - a voxel with no
    // bead within max_dist comes back at max_dist and lands in the last bin.
    const double bin_width = desc->max_dist / (double)num_bins;

    // A triclinic cell is covered by a grid over its bounding box, and the query wraps the corners
    // which stick out back into the cell. Those voxels are periodic images of voxels already counted,
    // so including them would weight part of the cell twice and quietly bias every fraction read off
    // the profile. Test the fractional coordinate and leave them out of the statistics. The field is
    // still written for them, since the channel sweep wants a complete grid.
    double Icell[3][3] = {};
    uint32_t cell_flags = 0;
    if (desc->cell) {
        cell_flags = md_unitcell_flags(desc->cell);
        md_unitcell_I_extract_double(Icell, desc->cell);
    }
    const bool mask_cell = (cell_flags & MD_UNITCELL_TRICLINIC) != 0;
    const bool pbc_axis[3] = {
        (cell_flags & MD_UNITCELL_PBC_X) != 0,
        (cell_flags & MD_UNITCELL_PBC_Y) != 0,
        (cell_flags & MD_UNITCELL_PBC_Z) != 0,
    };

    enum { TILE_SIZE = VOID_FIELD_TILE_DIM * VOID_FIELD_TILE_DIM * VOID_FIELD_TILE_DIM };
    float qx[TILE_SIZE], qy[TILE_SIZE], qz[TILE_SIZE];
    float dist[TILE_SIZE];
    int   vi[TILE_SIZE], vj[TILE_SIZE], vk[TILE_SIZE];

    for (uint32_t t = tile_beg; t < tile_end; ++t) {
        const int tx = (int)(t % (uint32_t)tiles[0]);
        const int ty = (int)((t / (uint32_t)tiles[0]) % (uint32_t)tiles[1]);
        const int tz = (int)(t / ((uint32_t)tiles[0] * (uint32_t)tiles[1]));

        // A tile is a compact block of voxels, which is what the batched query is efficient at
        int n = 0;
        for (int k = 0; k < TD; ++k) {
            const int z = tz * TD + k;
            if (z >= grid.dim[2]) break;
            for (int j = 0; j < TD; ++j) {
                const int y = ty * TD + j;
                if (y >= grid.dim[1]) break;
                for (int i = 0; i < TD; ++i) {
                    const int x = tx * TD + i;
                    if (x >= grid.dim[0]) break;
                    qx[n] = grid.origin.x + ((float)x + 0.5f) * grid.spacing.x;
                    qy[n] = grid.origin.y + ((float)y + 0.5f) * grid.spacing.y;
                    qz[n] = grid.origin.z + ((float)z + 0.5f) * grid.spacing.z;
                    vi[n] = x;
                    vj[n] = y;
                    vk[n] = z;
                    n += 1;
                }
            }
        }
        if (n == 0) continue;

        void_beads_nearest(&desc->beads, qx, qy, qz, (size_t)n, desc->max_dist, NULL, dist);

        for (int p = 0; p < n; ++p) {
            const float d = dist[p];

            if (desc->field) {
                const size_t idx = ((size_t)vk[p] * (size_t)grid.dim[1] + (size_t)vj[p]) * (size_t)grid.dim[0] + (size_t)vi[p];
                desc->field[idx] = MAX(0.0f, d);
            }

            if (mask_cell) {
                bool outside = false;
                for (int a = 0; a < 3 && !outside; ++a) {
                    if (!pbc_axis[a]) continue;
                    const double f = Icell[0][a] * (double)qx[p]
                                   + Icell[1][a] * (double)qy[p]
                                   + Icell[2][a] * (double)qz[p];
                    outside = (f < 0.0 || f >= 1.0);
                }
                if (outside) continue;
            }

            const uint32_t slab = void_profile_slab_of(vk[p], pps);

            accum->total[slab] += 1;
            accum->d_min = MIN(accum->d_min, d);
            accum->d_max = MAX(accum->d_max, d);
            if (d >= max_dist) accum->num_clamped += 1;

            if (d <= 0.0f) {
                accum->solid[slab] += 1;
            } else {
                const uint32_t bin = void_profile_bin_of((double)d, bin_width, num_bins);
                accum->hist[(size_t)slab * num_bins + bin] += 1;
            }
        }
    }
}

void_profile_t void_field_profile(const void_field_accum_t* accum, const void_field_desc_t* desc) {
    void_profile_t p = {};
    if (!accum || !desc || !desc->grid || desc->planes_per_slab == 0 || desc->num_bins == 0) return p;
    const md_grid_t& grid = *desc->grid;
    p.hist         = accum->hist;
    p.solid        = accum->solid;
    p.total        = accum->total;
    p.num_slabs    = void_profile_num_slabs(grid.dim[2], desc->planes_per_slab);
    p.num_bins     = desc->num_bins;
    p.bin_width    = desc->max_dist / (double)desc->num_bins;
    p.z_min        = (double)grid.origin.z;
    p.z_max        = (double)grid.origin.z + (double)grid.spacing.z * (double)grid.dim[2];
    p.slab_height  = (double)grid.spacing.z * (double)desc->planes_per_slab;
    p.voxel_volume = (double)grid.spacing.x * (double)grid.spacing.y * (double)grid.spacing.z;
    return p;
}

// =================================================================================================
// Surface topography
// =================================================================================================

namespace {

// Past this many steps a column is taken to be in contact where it stands. Sphere tracing only
// slows down when the ray grazes a bead, and there the remaining error is along a wall that is
// nearly vertical, so this bounds the cost without moving an answer by anything measurable.
constexpr int HEIGHTMAP_MAX_ITER = 4096;

// March a patch of columns from z_start towards z_end (dir = -1 downwards, +1 upwards) and write the
// probe apex height at contact, or NaN for a column the probe clears entirely.
//
// Each step is the gap d - R, never less than tol. Stepping by the gap alone converges only
// asymptotically where the column grazes a bead, and stopping once the gap is below tol would then
// report a contact up to sqrt(2 R tol) too high on the flank of every bead. Taking at least tol
// instead lets the column cross into contact, and the crossing is located by interpolating the gap
// between the last step outside and the first inside, which is exact for a head on approach.
void heightmap_march(float* out, const size_t* col, const float* cx, const float* cy, int n,
                     double z_start, double z_end, double dir, const void_heightmap_desc_t* desc, double tol) {
    enum { N = VOID_FIELD_TILE_DIM * VOID_FIELD_TILE_DIM };
    float  qx[N], qy[N], qz[N], dist[N];
    double z[N], z_prev[N], g_prev[N];
    int    act[N];

    const double R = desc->probe_radius;
    int num_act = n;
    for (int i = 0; i < n; ++i) {
        z[i]      = z_start;
        z_prev[i] = z_start;
        g_prev[i] = -1.0;      // No step taken yet
        act[i]    = i;
    }

    for (int iter = 0; iter < HEIGHTMAP_MAX_ITER && num_act > 0; ++iter) {
        for (int a = 0; a < num_act; ++a) {
            const int i = act[a];
            qx[a] = cx[i];
            qy[a] = cy[i];
            qz[a] = (float)z[i];
        }
        void_beads_nearest(&desc->beads, qx, qy, qz, (size_t)num_act, desc->max_dist, NULL, dist);

        int keep = 0;
        for (int a = 0; a < num_act; ++a) {
            const int    i   = act[a];
            const double gap = (double)dist[a] - R;
            if (gap <= 0.0) {
                double zc = z[i];
                if (g_prev[i] > 0.0) {
                    // The gap changed sign between the last two samples
                    const double t = g_prev[i] / (g_prev[i] - gap);
                    zc = z_prev[i] + t * (z[i] - z_prev[i]);
                }
                // A column which starts in contact stays where it started
                out[col[i]] = (float)(zc + dir * R);
                continue;
            }
            z_prev[i] = z[i];
            g_prev[i] = gap;
            z[i] += dir * (gap > tol ? gap : tol);
            if ((dir < 0.0 && z[i] < z_end) || (dir > 0.0 && z[i] > z_end)) {
                out[col[i]] = NAN;
                continue;
            }
            act[keep++] = i;
        }
        num_act = keep;
    }

    for (int a = 0; a < num_act; ++a) {
        const int i = act[a];
        out[col[i]] = (float)(z[i] + dir * R);
    }
}

}  // namespace

uint32_t void_heightmap_num_patches(const md_grid_t* grid) {
    if (!grid || grid->dim[0] <= 0 || grid->dim[1] <= 0) return 0;
    const uint32_t TD = VOID_FIELD_TILE_DIM;
    return (((uint32_t)grid->dim[0] + TD - 1) / TD) * (((uint32_t)grid->dim[1] + TD - 1) / TD);
}

void void_heightmap_eval_patches(float* out_top, float* out_bot, const void_heightmap_desc_t* desc, uint32_t patch_beg, uint32_t patch_end) {
    if (!desc || !desc->beads.acc || !desc->grid) return;
    if (!out_top && !out_bot) return;
    const md_grid_t& grid = *desc->grid;
    if (grid.dim[0] <= 0 || grid.dim[1] <= 0 || grid.dim[2] <= 0) return;

    const int TD = VOID_FIELD_TILE_DIM;
    const int px = (grid.dim[0] + TD - 1) / TD;
    patch_end = MIN(patch_end, void_heightmap_num_patches(desc->grid));

    const double R = desc->probe_radius > 0.0 ? desc->probe_radius : 0.0;
    const double min_spacing = MIN((double)grid.spacing.x, MIN((double)grid.spacing.y, (double)grid.spacing.z));
    const double tol = (desc->tol > 0.0) ? desc->tol : 0.01 * min_spacing;

    const double z_lo = (double)grid.origin.z;
    const double z_hi = (double)grid.origin.z + (double)grid.spacing.z * (double)grid.dim[2];

    // A probe the query cannot resolve reports every column as touching at once. Say nothing rather
    // than something wrong.
    const bool resolvable = desc->max_dist > R;

    enum { N = VOID_FIELD_TILE_DIM * VOID_FIELD_TILE_DIM };
    float  cx[N], cy[N];
    size_t col[N];

    for (uint32_t t = patch_beg; t < patch_end; ++t) {
        const int tx = (int)(t % (uint32_t)px);
        const int ty = (int)(t / (uint32_t)px);

        int n = 0;
        for (int j = 0; j < TD; ++j) {
            const int y = ty * TD + j;
            if (y >= grid.dim[1]) break;
            for (int i = 0; i < TD; ++i) {
                const int x = tx * TD + i;
                if (x >= grid.dim[0]) break;
                cx[n]  = grid.origin.x + ((float)x + 0.5f) * grid.spacing.x;
                cy[n]  = grid.origin.y + ((float)y + 0.5f) * grid.spacing.y;
                col[n] = (size_t)y * (size_t)grid.dim[0] + (size_t)x;
                n += 1;
            }
        }
        if (n == 0) continue;

        if (!resolvable) {
            for (int i = 0; i < n; ++i) {
                if (out_top) out_top[col[i]] = NAN;
                if (out_bot) out_bot[col[i]] = NAN;
            }
            continue;
        }

        if (out_top) heightmap_march(out_top, col, cx, cy, n, z_hi, z_lo, -1.0, desc, tol);
        if (out_bot) heightmap_march(out_bot, col, cx, cy, n, z_lo, z_hi, +1.0, desc, tol);
    }
}

void void_heightmap_stats(void_heightmap_stats_t* out, const float* h, size_t count) {
    if (!out) return;
    memset(out, 0, sizeof(*out));
    if (!h) return;

    double sum = 0.0;
    double lo  = DBL_MAX;
    double hi  = -DBL_MAX;
    for (size_t i = 0; i < count; ++i) {
        const float v = h[i];
        if (!isfinite(v)) {
            out->num_open += 1;
            continue;
        }
        out->num_valid += 1;
        sum += (double)v;
        lo = MIN(lo, (double)v);
        hi = MAX(hi, (double)v);
    }
    if (out->num_valid == 0) return;

    // Two passes: the deviations are small against the heights, which sit at world z
    const double mean = sum / (double)out->num_valid;
    double s2 = 0.0, s1 = 0.0;
    for (size_t i = 0; i < count; ++i) {
        const float v = h[i];
        if (!isfinite(v)) continue;
        const double d = (double)v - mean;
        s2 += d * d;
        s1 += fabs(d);
    }
    out->mean = mean;
    out->min  = lo;
    out->max  = hi;
    out->rq   = sqrt(s2 / (double)out->num_valid);
    out->ra   = s1 / (double)out->num_valid;
}

// =================================================================================================
// Channels
// =================================================================================================

#define CHANNEL_FLAG_TOP    0x1
#define CHANNEL_FLAG_BOTTOM 0x2

typedef struct sweep_ctx_t {
    // Union find over the qualifying voxels, with ids allocated in sweep order so that an older
    // component always carries the smaller id.
    md_array(uint32_t) parent;
    md_array(uint32_t) size;
    md_array(uint8_t)  flags;     // Which faces the component touches
    md_array(uint32_t) node_at;   // Branch currently carrying this component; valid at the root only

    md_array(int32_t)  node_slab; // Scratch: last slab a branch was seen in

    channel_tree_t* tree;
    struct md_allocator_i* alloc;
    float z_now;
    bool  want_tree;
} sweep_ctx_t;

static uint32_t uf_make(sweep_ctx_t* c, uint8_t flags) {
    const uint32_t id = (uint32_t)md_array_size(c->parent);
    md_array_push(c->parent, id, c->alloc);
    md_array_push(c->size,   1u, c->alloc);
    md_array_push(c->flags,  flags, c->alloc);
    if (c->want_tree) {
        md_array_push(c->node_at, CHANNEL_INVALID_INDEX, c->alloc);
    }
    return id;
}

static uint32_t uf_find(sweep_ctx_t* c, uint32_t x) {
    while (c->parent[x] != x) {
        c->parent[x] = c->parent[c->parent[x]];   // path halving
        x = c->parent[x];
    }
    return x;
}

static uint32_t node_new(sweep_ctx_t* c, float z, bool at_top) {
    channel_node_t n;
    MEMSET(&n, 0, sizeof(n));
    n.parent        = CHANNEL_INVALID_INDEX;
    n.first_child   = CHANNEL_INVALID_INDEX;
    n.next_sibling  = CHANNEL_INVALID_INDEX;
    n.z_top         = z;
    n.z_bot         = z;
    n.max_clearance = 0.0f;
    n.min_clearance = FLT_MAX;
    n.reaches_top   = at_top;

    md_array_push(c->tree->nodes, n, c->alloc);
    md_array_push(c->node_slab, INT32_MIN, c->alloc);
    return (uint32_t)md_array_size(c->tree->nodes) - 1;
}

static void node_attach(sweep_ctx_t* c, uint32_t child, uint32_t parent, float z) {
    channel_node_t* nodes = c->tree->nodes;
    nodes[child].parent  = parent;
    nodes[child].z_bot   = MIN(nodes[child].z_bot, z);
    nodes[child].next_sibling = nodes[parent].first_child;
    nodes[parent].first_child = child;
}

static void uf_union(sweep_ctx_t* c, uint32_t a, uint32_t b) {
    uint32_t ra = uf_find(c, a);
    uint32_t rb = uf_find(c, b);
    if (ra == rb) return;

    const uint32_t na = c->want_tree ? c->node_at[ra] : CHANNEL_INVALID_INDEX;
    const uint32_t nb = c->want_tree ? c->node_at[rb] : CHANNEL_INVALID_INDEX;
    const uint8_t  fl = (uint8_t)(c->flags[ra] | c->flags[rb]);

    if (c->size[ra] < c->size[rb]) {
        const uint32_t t = ra; ra = rb; rb = t;
    }
    c->parent[rb] = ra;
    c->size[ra]  += c->size[rb];
    c->flags[ra]  = fl;

    if (!c->want_tree) return;

    uint32_t node;
    if (na == CHANNEL_INVALID_INDEX) {
        node = nb;
    } else if (nb == CHANNEL_INVALID_INDEX || na == nb) {
        node = na;
    } else {
        // Two branches which already have an identity meet here, so this is where they join.
        //
        // Three or more branches meeting in the same slab is one junction, not a chain of them. The
        // unions arrive pairwise, so a merge node already made at this z is folded into rather than
        // nested under, which keeps the junction n-ary. Nesting instead leaves the lower node with no
        // z extent and - because the per slab accumulation below only ever reaches the node the
        // component currently carries - with no voxels, no path and no clearance either, which is
        // what used to fill the diagram with zero rows.
        //
        // A node made at this z is necessarily one of those junctions: branch nodes for this slab are
        // not created until every union below has been applied.
        const channel_node_t* nodes = c->tree->nodes;
        const bool na_fresh = (nodes[na].z_top == c->z_now) && (nodes[na].num_voxels == 0);
        const bool nb_fresh = (nodes[nb].z_top == c->z_now) && (nodes[nb].num_voxels == 0);

        if (na_fresh && nb_fresh) {
            // Two junctions made at this z, reached from different sides: one adopts the other's
            // branches. Reading each next_sibling before reattaching keeps the list well formed. The
            // emptied node is left unreferenced and is never visited again, since every traversal
            // starts from a root.
            node = na;
            uint32_t child = c->tree->nodes[nb].first_child;
            c->tree->nodes[nb].first_child = CHANNEL_INVALID_INDEX;
            while (child != CHANNEL_INVALID_INDEX) {
                const uint32_t next = c->tree->nodes[child].next_sibling;
                node_attach(c, child, node, c->z_now);
                child = next;
            }
        } else if (na_fresh) {
            node = na;
            node_attach(c, nb, node, c->z_now);
        } else if (nb_fresh) {
            node = nb;
            node_attach(c, na, node, c->z_now);
        } else {
            node = node_new(c, c->z_now, false);
            node_attach(c, na, node, c->z_now);
            node_attach(c, nb, node, c->z_now);
        }
    }
    c->node_at[ra] = node;
}

void channel_tree_free(channel_tree_t* tree) {
    if (!tree || !tree->alloc) return;
    for (size_t i = 0; i < md_array_size(tree->nodes); ++i) {
        md_array_free(tree->nodes[i].path, tree->alloc);
    }
    md_array_free(tree->nodes, tree->alloc);
    md_array_free(tree->roots, tree->alloc);
    MEMSET(tree, 0, sizeof(channel_tree_t));
}

bool channel_sweep(channel_tree_t* out, const channel_field_t* field, double probe_radius, bool want_tree, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(field);
    ASSERT(alloc);

    MEMSET(out, 0, sizeof(channel_tree_t));
    out->alloc = alloc;
    out->probe_radius = probe_radius;

    if (!field->data) return false;
    if (field->pbc[2]) return false;                       // The sweep axis has to have two faces
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    const float r = (float)probe_radius;
    const size_t plane = (size_t)nx * (size_t)ny;

    sweep_ctx_t ctx;
    MEMSET(&ctx, 0, sizeof(ctx));
    ctx.tree      = out;
    ctx.alloc     = alloc;
    ctx.want_tree = want_tree;

    md_array(uint32_t) id_prev = 0;
    md_array(uint32_t) id_cur  = 0;
    md_array_resize(id_prev, plane, alloc);
    md_array_resize(id_cur,  plane, alloc);
    for (size_t i = 0; i < plane; ++i) id_prev[i] = CHANNEL_INVALID_INDEX;

    // Sweep from the far face towards the near one. Only the two slabs in flight are ever indexed by
    // position, so the per position storage is one plane rather than one volume.
    for (int k = nz - 1; k >= 0; --k) {
        const float z = field->origin[2] + ((float)k + 0.5f) * field->spacing[2];
        ctx.z_now = z;

        const float* slab = field->data + (size_t)k * plane;

        uint8_t face = 0;
        if (k == nz - 1) face |= CHANNEL_FLAG_TOP;
        if (k == 0)      face |= CHANNEL_FLAG_BOTTOM;

        for (size_t i = 0; i < plane; ++i) {
            id_cur[i] = (slab[i] >= r) ? uf_make(&ctx, face) : CHANNEL_INVALID_INDEX;
        }

        // In plane neighbours
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const size_t p = (size_t)y * nx + x;
                if (id_cur[p] == CHANNEL_INVALID_INDEX) continue;
                if (x > 0 && id_cur[p - 1] != CHANNEL_INVALID_INDEX)  uf_union(&ctx, id_cur[p], id_cur[p - 1]);
                if (y > 0 && id_cur[p - nx] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[p], id_cur[p - nx]);
            }
        }
        // Periodic seams, joined after the interior so the raster scan above stays a plain two neighbour pass
        if (field->pbc[0] && nx > 1) {
            for (int y = 0; y < ny; ++y) {
                const size_t a = (size_t)y * nx + (nx - 1), b = (size_t)y * nx;
                if (id_cur[a] != CHANNEL_INVALID_INDEX && id_cur[b] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[a], id_cur[b]);
            }
        }
        if (field->pbc[1] && ny > 1) {
            for (int x = 0; x < nx; ++x) {
                const size_t a = (size_t)(ny - 1) * nx + x, b = (size_t)x;
                if (id_cur[a] != CHANNEL_INVALID_INDEX && id_cur[b] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[a], id_cur[b]);
            }
        }
        // Onto the slab already swept
        for (size_t i = 0; i < plane; ++i) {
            if (id_cur[i] != CHANNEL_INVALID_INDEX && id_prev[i] != CHANNEL_INVALID_INDEX) {
                uf_union(&ctx, id_cur[i], id_prev[i]);
            }
        }

        if (want_tree) {
            // A component with no branch yet is one which appears at this slab
            for (size_t i = 0; i < plane; ++i) {
                if (id_cur[i] == CHANNEL_INVALID_INDEX) continue;
                const uint32_t root = uf_find(&ctx, id_cur[i]);
                if (ctx.node_at[root] == CHANNEL_INVALID_INDEX) {
                    ctx.node_at[root] = node_new(&ctx, z, k == nz - 1);
                }
            }

            // Representative point per branch per slab: its widest voxel in this slab
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const size_t i = (size_t)y * nx + x;
                    if (id_cur[i] == CHANNEL_INVALID_INDEX) continue;

                    const uint32_t node = ctx.node_at[uf_find(&ctx, id_cur[i])];
                    channel_node_t* n = ctx.tree->nodes + node;

                    n->num_voxels    += 1;
                    n->max_clearance  = MAX(n->max_clearance, slab[i]);
                    n->z_bot          = MIN(n->z_bot, z);

                    if (ctx.node_slab[node] != k) {
                        ctx.node_slab[node] = k;
                        vec4_t pt = {0, 0, 0, -FLT_MAX};
                        md_array_push(n->path, pt, alloc);
                        n = ctx.tree->nodes + node;   // push may have moved the array
                    }
                    vec4_t* pt = md_array_last(n->path);
                    if (slab[i] > pt->w) {
                        pt->x = field->origin[0] + ((float)x + 0.5f) * field->spacing[0];
                        pt->y = field->origin[1] + ((float)y + 0.5f) * field->spacing[1];
                        pt->z = z;
                        pt->w = slab[i];
                    }
                }
            }
        }

        md_array(uint32_t) tmp = id_prev;
        id_prev = id_cur;
        id_cur  = tmp;
    }

    // Components, and which of them get all the way through
    for (uint32_t i = 0; i < (uint32_t)md_array_size(ctx.parent); ++i) {
        if (uf_find(&ctx, i) != i) continue;
        out->num_components += 1;
        const bool spans = (ctx.flags[i] & CHANNEL_FLAG_TOP) && (ctx.flags[i] & CHANNEL_FLAG_BOTTOM);
        if (spans) out->num_spanning += 1;
        if (want_tree) {
            const uint32_t node = ctx.node_at[i];
            if (node != CHANNEL_INVALID_INDEX) {
                out->nodes[node].reaches_bottom = (ctx.flags[i] & CHANNEL_FLAG_BOTTOM) != 0;
                md_array_push(out->roots, node, alloc);
            }
        }
    }

    if (want_tree) {
        // A merge node is always created after the branches it joins, so one ascending pass carries
        // "something above me reaches the far face" all the way to the roots.
        for (size_t i = 0; i < md_array_size(out->nodes); ++i) {
            channel_node_t* n = out->nodes + i;
            for (size_t p = 0; p < md_array_size(n->path); ++p) {
                n->min_clearance = MIN(n->min_clearance, n->path[p].w);
            }
            if (n->min_clearance == FLT_MAX) n->min_clearance = 0.0f;
            if (n->reaches_top && n->parent != CHANNEL_INVALID_INDEX) {
                out->nodes[n->parent].reaches_top = true;
            }
        }
    }

    md_array_free(id_prev, alloc);
    md_array_free(id_cur,  alloc);
    md_array_free(ctx.parent, alloc);
    md_array_free(ctx.size, alloc);
    md_array_free(ctx.flags, alloc);
    md_array_free(ctx.node_at, alloc);
    md_array_free(ctx.node_slab, alloc);

    return true;
}

// ---------------------------------------------------------------------------------------------
// Percolation: one pass in order of decreasing clearance
// ---------------------------------------------------------------------------------------------

#define PERC_INVALID 0xFFFFFFFFu
#define CHANNEL_FLAG_BOTH (CHANNEL_FLAG_TOP | CHANNEL_FLAG_BOTTOM)

typedef struct perc_ctx_t {
    // Union find over the active voxels, indexed by position in the clearance-ordered list. That
    // ordering is what does the work: an id smaller than the one being inserted is, by construction,
    // a voxel that is already in - so "have I seen this neighbour" is a comparison, not a lookup.
    uint32_t* parent;
    uint32_t* size;
    uint8_t*  flags;
    uint16_t* min_k;
    uint16_t* max_k;

    // Running totals, maintained through the unions rather than recomputed per sample. A component
    // leaves the totals before it is merged and the merged one re-enters, so every sample is a read
    // of four integers.
    uint64_t open_voxels;
    uint64_t span_voxels;
    uint32_t num_roots;
    uint32_t num_spanning;

    // Monotone: a component's reach only grows and components only merge, so these never need to be
    // taken back out.
    int32_t deepest_top;
    int32_t highest_bot;
} perc_ctx_t;

static uint32_t perc_find(perc_ctx_t* c, uint32_t x) {
    while (c->parent[x] != x) {
        c->parent[x] = c->parent[c->parent[x]];   // path halving
        x = c->parent[x];
    }
    return x;
}

static void perc_enter(perc_ctx_t* c, uint32_t r) {
    const uint8_t f = c->flags[r];
    if (f) c->open_voxels += c->size[r];
    if (f == CHANNEL_FLAG_BOTH) {
        c->span_voxels += c->size[r];
        c->num_spanning += 1;
    }
    if (f & CHANNEL_FLAG_TOP)    c->deepest_top = MIN(c->deepest_top, (int32_t)c->min_k[r]);
    if (f & CHANNEL_FLAG_BOTTOM) c->highest_bot = MAX(c->highest_bot, (int32_t)c->max_k[r]);
}

static void perc_leave(perc_ctx_t* c, uint32_t r) {
    const uint8_t f = c->flags[r];
    if (f) c->open_voxels -= c->size[r];
    if (f == CHANNEL_FLAG_BOTH) {
        c->span_voxels -= c->size[r];
        c->num_spanning -= 1;
    }
}

// Returns true when this union is what made a component touch both faces for the first time.
static bool perc_union(perc_ctx_t* c, uint32_t a, uint32_t b) {
    uint32_t ra = perc_find(c, a);
    uint32_t rb = perc_find(c, b);
    if (ra == rb) return false;

    const bool was_spanning = (c->flags[ra] == CHANNEL_FLAG_BOTH) || (c->flags[rb] == CHANNEL_FLAG_BOTH);

    perc_leave(c, ra);
    perc_leave(c, rb);

    if (c->size[ra] < c->size[rb]) { const uint32_t t = ra; ra = rb; rb = t; }
    c->parent[rb] = ra;
    c->size[ra]  += c->size[rb];
    c->flags[ra]  = (uint8_t)(c->flags[ra] | c->flags[rb]);
    c->min_k[ra]  = MIN(c->min_k[ra], c->min_k[rb]);
    c->max_k[ra]  = MAX(c->max_k[ra], c->max_k[rb]);
    c->num_roots -= 1;

    perc_enter(c, ra);
    return !was_spanning && c->flags[ra] == CHANNEL_FLAG_BOTH;
}

void channel_percolation_free(channel_percolation_t* perc) {
    if (!perc || !perc->alloc) return;
    md_array_free(perc->radius, perc->alloc);
    md_array_free(perc->frac_void, perc->alloc);
    md_array_free(perc->frac_open, perc->alloc);
    md_array_free(perc->frac_spanning, perc->alloc);
    md_array_free(perc->z_from_top, perc->alloc);
    md_array_free(perc->z_from_bottom, perc->alloc);
    md_array_free(perc->num_components, perc->alloc);
    md_array_free(perc->num_spanning, perc->alloc);
    MEMSET(perc, 0, sizeof(channel_percolation_t));
}

bool channel_percolate(channel_percolation_t* out, const channel_field_t* field, double r_min, uint32_t num_samples, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(field);
    ASSERT(alloc);

    MEMSET(out, 0, sizeof(channel_percolation_t));
    out->alloc = alloc;

    if (!field->data) return false;
    if (field->pbc[2]) return false;                       // The sweep axis has to have two faces
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;
    if (nz > 65535) return false;                          // min_k / max_k are uint16

    const size_t plane = (size_t)nx * (size_t)ny;
    const size_t N     = plane * (size_t)nz;
    if (N > 0xFFFFFFF0u) return false;                     // Voxel indices are uint32 throughout

    const float rmin = (float)r_min;

    // What the run is sized by, and the widest clearance the curves have to reach.
    float  d_max = rmin;
    size_t M     = 0;
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d >= rmin) {
            M += 1;
            if (d > d_max) d_max = d;
        }
    }
    out->num_active = M;
    if (M == 0) return false;

    const uint32_t ns  = (uint32_t)CLAMP((int)num_samples, 2, 4096);
    // Internal buckets are far finer than the reported samples. The bucket width is what bounds the
    // error on r_c, and at these counts it lands orders of magnitude below the grid's own limit -
    // half a voxel of unresolved gap - so the grid stays the thing that limits the answer.
    const uint32_t per = MAX(1u, (4096u + ns - 1u) / ns);
    const uint32_t NB  = ns * per;

    const double span  = MAX(1.0e-6, (double)d_max - (double)rmin);
    const double width = span / (double)NB;

    const size_t bytes = N * sizeof(uint32_t)
                       + M * (sizeof(uint32_t) * 3 + sizeof(uint8_t) + sizeof(uint16_t) * 2)
                       + (size_t)NB * sizeof(uint32_t) * 2;
    out->bytes = bytes;

    uint32_t* cid = (uint32_t*)md_alloc(alloc, N * sizeof(uint32_t));
    uint32_t* vox = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    uint32_t* off = (uint32_t*)md_alloc(alloc, (size_t)NB * sizeof(uint32_t));
    uint32_t* cur = (uint32_t*)md_alloc(alloc, (size_t)NB * sizeof(uint32_t));

    perc_ctx_t c;
    MEMSET(&c, 0, sizeof(c));
    c.parent = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    c.size   = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    c.flags  = (uint8_t*) md_alloc(alloc, M * sizeof(uint8_t));
    c.min_k  = (uint16_t*)md_alloc(alloc, M * sizeof(uint16_t));
    c.max_k  = (uint16_t*)md_alloc(alloc, M * sizeof(uint16_t));
    c.deepest_top = INT32_MAX;
    c.highest_bot = INT32_MIN;

    if (!cid || !vox || !off || !cur || !c.parent || !c.size || !c.flags || !c.min_k || !c.max_k) {
        md_free(alloc, cid, N * sizeof(uint32_t));
        md_free(alloc, vox, M * sizeof(uint32_t));
        md_free(alloc, off, (size_t)NB * sizeof(uint32_t));
        md_free(alloc, cur, (size_t)NB * sizeof(uint32_t));
        md_free(alloc, c.parent, M * sizeof(uint32_t));
        md_free(alloc, c.size,   M * sizeof(uint32_t));
        md_free(alloc, c.flags,  M * sizeof(uint8_t));
        md_free(alloc, c.min_k,  M * sizeof(uint16_t));
        md_free(alloc, c.max_k,  M * sizeof(uint16_t));
        return false;
    }

    MEMSET(off, 0, (size_t)NB * sizeof(uint32_t));

    // Counting sort by clearance, descending. A comparison sort of a hundred million voxels would
    // cost more than the union find it feeds.
    for (size_t i = 0; i < N; ++i) {
        cid[i] = PERC_INVALID;
        const float d = field->data[i];
        if (d < rmin) continue;
        int b = (int)(((double)d - (double)rmin) / width);
        b = CLAMP(b, 0, (int)NB - 1);
        off[b] += 1;
    }
    {   // Bucket NB-1 is widest and goes first, so the offsets accumulate downward
        uint32_t acc = 0;
        for (int b = (int)NB - 1; b >= 0; --b) {
            const uint32_t n = off[b];
            off[b] = acc;
            cur[b] = acc;
            acc += n;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d < rmin) continue;
        int b = (int)(((double)d - (double)rmin) / width);
        b = CLAMP(b, 0, (int)NB - 1);
        const uint32_t p = cur[b]++;
        vox[p] = (uint32_t)i;
        cid[i] = p;
    }

    md_array_resize(out->radius,         ns, alloc);
    md_array_resize(out->frac_void,      ns, alloc);
    md_array_resize(out->frac_open,      ns, alloc);
    md_array_resize(out->frac_spanning,  ns, alloc);
    md_array_resize(out->z_from_top,     ns, alloc);
    md_array_resize(out->z_from_bottom,  ns, alloc);
    md_array_resize(out->num_components, ns, alloc);
    md_array_resize(out->num_spanning,   ns, alloc);

    const double z_lo   = (double)field->origin[2];
    const double z_hi   = (double)field->origin[2] + (double)field->spacing[2] * (double)nz;
    const double inv_N  = 1.0 / (double)N;

    uint32_t p = 0;
    for (int b = (int)NB - 1; b >= 0; --b) {
        const uint32_t end = (b > 0) ? off[b - 1] : (uint32_t)M;

        for (; p < end; ++p) {
            const uint32_t v = vox[p];
            const int k = (int)((size_t)v / plane);
            const int y = (int)(((size_t)v % plane) / (size_t)nx);
            const int x = (int)((size_t)v % (size_t)nx);

            c.parent[p] = p;
            c.size[p]   = 1;
            c.min_k[p]  = (uint16_t)k;
            c.max_k[p]  = (uint16_t)k;
            c.flags[p]  = (uint8_t)(((k == nz - 1) ? CHANNEL_FLAG_TOP : 0) | ((k == 0) ? CHANNEL_FLAG_BOTTOM : 0));
            c.num_roots += 1;
            perc_enter(&c, p);

            bool connected = (!out->has_r_c && c.flags[p] == CHANNEL_FLAG_BOTH);

            // Six neighbours. An id below p is a voxel already inserted, since the list is in
            // insertion order - no separate "visited" mark is needed or kept.
            //
            // Both wrap directions are probed even though either alone would do: a seam pair is the
            // same pair seen from both sides, and whichever voxel is inserted second makes the
            // union. Keeping the loop symmetric costs two comparisons and means the periodic case
            // reads the same as the interior one.
            uint32_t nb[6];
            int      nn = 0;
            if (x > 0)               nb[nn++] = cid[(size_t)v - 1];
            else if (field->pbc[0] && nx > 1) nb[nn++] = cid[(size_t)v + (nx - 1)];
            if (x < nx - 1)          nb[nn++] = cid[(size_t)v + 1];
            else if (field->pbc[0] && nx > 1) nb[nn++] = cid[(size_t)v - (nx - 1)];
            if (y > 0)               nb[nn++] = cid[(size_t)v - nx];
            else if (field->pbc[1] && ny > 1) nb[nn++] = cid[(size_t)v + (size_t)(ny - 1) * nx];
            if (y < ny - 1)          nb[nn++] = cid[(size_t)v + nx];
            else if (field->pbc[1] && ny > 1) nb[nn++] = cid[(size_t)v - (size_t)(ny - 1) * nx];
            if (k > 0)               nb[nn++] = cid[(size_t)v - plane];
            if (k < nz - 1)          nb[nn++] = cid[(size_t)v + plane];

            for (int q = 0; q < nn; ++q) {
                if (nb[q] == PERC_INVALID || nb[q] >= p) continue;
                if (perc_union(&c, p, nb[q])) connected = true;
            }

            if (connected && !out->has_r_c) {
                // The voxel whose insertion joined the two faces is the tightest point of the widest
                // route, and its clearance is the critical radius. Bisecting on "does anything span"
                // brackets this number and never tells you which voxel it was.
                out->has_r_c   = true;
                out->r_c       = (double)field->data[v];
                out->throat[0] = field->origin[0] + ((float)x + 0.5f) * field->spacing[0];
                out->throat[1] = field->origin[1] + ((float)y + 0.5f) * field->spacing[1];
                out->throat[2] = field->origin[2] + ((float)k + 0.5f) * field->spacing[2];
            }
        }

        if ((uint32_t)b % per == 0) {
            const uint32_t s = (uint32_t)b / per;
            out->radius[s]         = (double)rmin + (double)b * width;
            out->frac_void[s]      = (double)p * inv_N;
            out->frac_open[s]      = (double)c.open_voxels * inv_N;
            out->frac_spanning[s]  = (double)c.span_voxels * inv_N;
            out->num_components[s] = c.num_roots;
            out->num_spanning[s]   = c.num_spanning;
            out->z_from_top[s]     = (c.deepest_top == INT32_MAX) ? z_hi
                                   : (double)field->origin[2] + ((double)c.deepest_top + 0.5) * (double)field->spacing[2];
            out->z_from_bottom[s]  = (c.highest_bot == INT32_MIN) ? z_lo
                                   : (double)field->origin[2] + ((double)c.highest_bot + 0.5) * (double)field->spacing[2];
        }
    }

    md_free(alloc, cid, N * sizeof(uint32_t));
    md_free(alloc, vox, M * sizeof(uint32_t));
    md_free(alloc, off, (size_t)NB * sizeof(uint32_t));
    md_free(alloc, cur, (size_t)NB * sizeof(uint32_t));
    md_free(alloc, c.parent, M * sizeof(uint32_t));
    md_free(alloc, c.size,   M * sizeof(uint32_t));
    md_free(alloc, c.flags,  M * sizeof(uint8_t));
    md_free(alloc, c.min_k,  M * sizeof(uint16_t));
    md_free(alloc, c.max_k,  M * sizeof(uint16_t));

    return true;
}

void channel_spanning_counts(uint32_t* out_counts, const double* radii, size_t num_radii, const channel_field_t* field, struct md_allocator_i* alloc) {
    for (size_t i = 0; i < num_radii; ++i) {
        channel_tree_t t;
        out_counts[i] = channel_sweep(&t, field, radii[i], false, alloc) ? t.num_spanning : 0;
        channel_tree_free(&t);
    }
}

double channel_critical_radius(const channel_field_t* field, double r_lo, double r_hi, double tol, struct md_allocator_i* alloc) {
    channel_tree_t t;

    // Nothing gets through even at the loosest radius asked for
    if (!channel_sweep(&t, field, r_lo, false, alloc)) return r_lo - tol;
    const bool lo_spans = t.num_spanning > 0;
    channel_tree_free(&t);
    if (!lo_spans) return r_lo - tol;

    if (channel_sweep(&t, field, r_hi, false, alloc) && t.num_spanning > 0) {
        channel_tree_free(&t);
        return r_hi;
    }
    channel_tree_free(&t);

    // Spanning is monotone in the radius: the superlevel set only shrinks as the probe grows, so a
    // radius which gets through implies every smaller one does.
    while (r_hi - r_lo > tol) {
        const double mid = 0.5 * (r_lo + r_hi);
        const bool spans = channel_sweep(&t, field, mid, false, alloc) && t.num_spanning > 0;
        channel_tree_free(&t);
        if (spans) r_lo = mid; else r_hi = mid;
    }
    return r_lo;
}

// =================================================================================================
// Pore network
// =================================================================================================

namespace {

// One per basin the pass opens. A region stays a pore of its own until a persistence merge folds it
// into an elder one; conn_parent joins on every contact and is only there for r_c.
struct pore_region_t {
    uint32_t pore_parent;
    uint32_t conn_parent;
    uint32_t peak_vox;      // Valid at a pore root
    uint32_t count;
    float    peak;
    float    face_top;
    float    face_bot;
    uint8_t  conn_faces;    // Valid at a conn root
};

struct pore_contact_t {
    uint32_t a, b;          // Pore roots when recorded; resolved to their final roots afterwards
    uint32_t vox;
    float    d;
};

inline uint32_t pore_find(pore_region_t* r, uint32_t x) {
    while (r[x].pore_parent != x) {
        r[x].pore_parent = r[r[x].pore_parent].pore_parent;
        x = r[x].pore_parent;
    }
    return x;
}

inline uint32_t conn_find(pore_region_t* r, uint32_t x) {
    while (r[x].conn_parent != x) {
        r[x].conn_parent = r[r[x].conn_parent].conn_parent;
        x = r[x].conn_parent;
    }
    return x;
}

inline uint32_t vert_find(uint32_t* parent, uint32_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

// The hash map indexes by the low bits of the key, and a pair key built from two small ids would
// pile into few buckets. The mix is a bijection, and the two values the map reserves are moved out
// of the way.
inline uint64_t pore_pair_key(uint32_t a, uint32_t b) {
    uint64_t x = (a < b) ? (((uint64_t)a << 32) | b) : (((uint64_t)b << 32) | a);
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    if (x >= MD_HASH_TOMBSTONE) x ^= (1ull << 63);
    return x;
}

inline void voxel_centre(float out[3], const channel_field_t* f, uint32_t v) {
    const size_t plane = (size_t)f->dim[0] * (size_t)f->dim[1];
    const int k = (int)((size_t)v / plane);
    const int y = (int)(((size_t)v % plane) / (size_t)f->dim[0]);
    const int x = (int)((size_t)v % (size_t)f->dim[0]);
    out[0] = f->origin[0] + ((float)x + 0.5f) * f->spacing[0];
    out[1] = f->origin[1] + ((float)y + 0.5f) * f->spacing[1];
    out[2] = f->origin[2] + ((float)k + 0.5f) * f->spacing[2];
}

}  // namespace

void pore_network_free(pore_network_t* net) {
    if (!net || !net->alloc) return;
    md_array_free(net->vertices,   net->alloc);
    md_array_free(net->edges,      net->alloc);
    md_array_free(net->adj_offset, net->alloc);
    md_array_free(net->adj,        net->alloc);
    MEMSET(net, 0, sizeof(pore_network_t));
}

bool pore_network_build(pore_network_t* out, const channel_field_t* field, double r_min, double merge, struct md_allocator_i* alloc,
                        pore_network_progress_fn progress, void* user) {
    ASSERT(out);
    ASSERT(field);
    ASSERT(alloc);

    MEMSET(out, 0, sizeof(pore_network_t));
    out->alloc = alloc;
    out->r_min = r_min;
    out->merge = merge;

    if (!field->data) return false;
    if (field->pbc[2]) return false;
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    const size_t plane = (size_t)nx * (size_t)ny;
    const size_t N     = plane * (size_t)nz;
    if (N > 0xFFFFFFF0u) return false;

    for (int a = 0; a < 3; ++a) {
        out->box_min[a] = field->origin[a];
        out->box_ext[a] = field->spacing[a] * (float)field->dim[a];
        out->pbc[a]     = field->pbc[a];
    }
    out->voxel_volume = (double)field->spacing[0] * (double)field->spacing[1] * (double)field->spacing[2];

    const float rmin = (float)r_min;
    float  d_max = rmin;
    size_t M     = 0;
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d >= rmin) {
            M += 1;
            if (d > d_max) d_max = d;
        }
    }
    out->num_active = M;
    if (M == 0) return false;
    if (progress && !progress(0.05f, user)) return false;

    // Counting sort, descending, as in channel_percolate but finer: the order inside a bucket is
    // voxel order rather than clearance, and the only thing that sees it is the persistence test,
    // which is told to ignore anything shallower than two buckets.
    const uint32_t NB    = 1u << 16;
    const double   width = MAX(1.0e-6, (double)d_max - (double)rmin) / (double)NB;
    const float    h     = (float)MAX(merge, 2.0 * width);

    uint32_t* label = (uint32_t*)md_alloc(alloc, N * sizeof(uint32_t));
    uint32_t* vox   = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    uint32_t* off   = (uint32_t*)md_alloc(alloc, (size_t)NB * sizeof(uint32_t));
    if (!label || !vox || !off) {
        md_free(alloc, label, N * sizeof(uint32_t));
        md_free(alloc, vox,   M * sizeof(uint32_t));
        md_free(alloc, off,   (size_t)NB * sizeof(uint32_t));
        return false;
    }

    MEMSET(off, 0, (size_t)NB * sizeof(uint32_t));
    for (size_t i = 0; i < N; ++i) {
        label[i] = PORE_INVALID;
        const float d = field->data[i];
        if (d < rmin) continue;
        const int b = CLAMP((int)(((double)d - (double)rmin) / width), 0, (int)NB - 1);
        off[b] += 1;
    }
    {
        uint32_t acc = 0;
        for (int b = (int)NB - 1; b >= 0; --b) {
            const uint32_t n = off[b];
            off[b] = acc;
            acc += n;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d < rmin) continue;
        const int b = CLAMP((int)(((double)d - (double)rmin) / width), 0, (int)NB - 1);
        vox[off[b]++] = (uint32_t)i;
    }
    md_free(alloc, off, (size_t)NB * sizeof(uint32_t));

    md_array(pore_region_t)  reg = 0;
    md_array(pore_contact_t) con = 0;
    md_hashmap32_t seen = {};
    seen.allocator = alloc;

    // The sweep is nearly all of the time; the sort before it and the graph after share the rest
    const size_t report_every = (size_t)1 << 16;
    bool abandoned = progress && !progress(0.15f, user);

    for (size_t p = 0; p < M && !abandoned; ++p) {
        if (progress && (p % report_every) == 0 && p > 0) {
            abandoned = !progress(0.15f + 0.8f * (float)((double)p / (double)M), user);
            if (abandoned) break;
        }
        const uint32_t v = vox[p];
        const float    d = field->data[v];
        const int k = (int)((size_t)v / plane);
        const int y = (int)(((size_t)v % plane) / (size_t)nx);
        const int x = (int)((size_t)v % (size_t)nx);

        size_t nb[6];
        int    nn = 0;
        if (x > 0)                        nb[nn++] = (size_t)v - 1;
        else if (field->pbc[0] && nx > 1) nb[nn++] = (size_t)v + (nx - 1);
        if (x < nx - 1)                   nb[nn++] = (size_t)v + 1;
        else if (field->pbc[0] && nx > 1) nb[nn++] = (size_t)v - (nx - 1);
        if (y > 0)                        nb[nn++] = (size_t)v - nx;
        else if (field->pbc[1] && ny > 1) nb[nn++] = (size_t)v + (size_t)(ny - 1) * nx;
        if (y < ny - 1)                   nb[nn++] = (size_t)v + nx;
        else if (field->pbc[1] && ny > 1) nb[nn++] = (size_t)v - (size_t)(ny - 1) * nx;
        if (k > 0)                        nb[nn++] = (size_t)v - plane;
        if (k < nz - 1)                   nb[nn++] = (size_t)v + plane;

        // The distinct pores already inserted around this voxel
        uint32_t roots[6];
        int      nr = 0;
        for (int q = 0; q < nn; ++q) {
            const uint32_t l = label[nb[q]];
            if (l == PORE_INVALID) continue;
            const uint32_t r = pore_find(reg, l);
            bool dup = false;
            for (int s = 0; s < nr; ++s) dup |= (roots[s] == r);
            if (!dup) roots[nr++] = r;
        }

        uint32_t own;
        if (nr == 0) {
            own = (uint32_t)md_array_size(reg);
            pore_region_t r = {};
            r.pore_parent = own;
            r.conn_parent = own;
            r.peak_vox    = v;
            r.peak        = d;
            r.face_top    = -1.0f;
            r.face_bot    = -1.0f;
            md_array_push(reg, r, alloc);
        } else {
            // The neighbour with the highest maximum, which is the basin steepest ascent leads to
            own = roots[0];
            for (int s = 1; s < nr; ++s) {
                if (reg[roots[s]].peak > reg[own].peak) own = roots[s];
            }
        }

        label[v] = own;
        reg[own].count += 1;
        if (k == nz - 1) {
            reg[own].face_top = MAX(reg[own].face_top, d);
            reg[conn_find(reg, own)].conn_faces |= PORE_FACE_TOP;
        }
        if (k == 0) {
            reg[own].face_bot = MAX(reg[own].face_bot, d);
            reg[conn_find(reg, own)].conn_faces |= PORE_FACE_BOTTOM;
        }

        for (int s = 0; s < nr; ++s) {
            const uint32_t ro = pore_find(reg, own);
            const uint32_t rr = pore_find(reg, roots[s]);
            if (ro == rr) continue;

            const uint32_t co = conn_find(reg, ro);
            const uint32_t cr = conn_find(reg, rr);
            if (co != cr) {
                reg[cr].conn_parent = co;
                reg[co].conn_faces |= reg[cr].conn_faces;
            }

            // Nothing wider than this voxel joins the two, so its clearance is their throat
            if (MIN(reg[ro].peak, reg[rr].peak) - d < h) {
                const uint32_t elder   = (reg[ro].peak >= reg[rr].peak) ? ro : rr;
                const uint32_t younger = (elder == ro) ? rr : ro;
                reg[younger].pore_parent = elder;
                reg[elder].count   += reg[younger].count;
                reg[elder].face_top = MAX(reg[elder].face_top, reg[younger].face_top);
                reg[elder].face_bot = MAX(reg[elder].face_bot, reg[younger].face_bot);
            } else {
                const uint64_t key = pore_pair_key(ro, rr);
                if (!md_hashmap_get(&seen, key)) {
                    const pore_contact_t c = { ro, rr, v, d };
                    md_hashmap_add(&seen, key, (uint32_t)md_array_size(con));
                    md_array_push(con, c, alloc);
                }
            }
        }

        if (!out->has_r_c && reg[conn_find(reg, own)].conn_faces == (PORE_FACE_TOP | PORE_FACE_BOTTOM)) {
            out->has_r_c = true;
            out->r_c     = (double)d;
            voxel_centre(out->throat, field, v);
        }
    }

    out->bytes = N * sizeof(uint32_t) + M * sizeof(uint32_t) + (size_t)NB * sizeof(uint32_t)
               + md_array_bytes(reg) + md_array_bytes(con) + (size_t)seen.num_buckets * (sizeof(uint64_t) + sizeof(uint32_t));

    md_free(alloc, label, N * sizeof(uint32_t));
    md_free(alloc, vox,   M * sizeof(uint32_t));
    md_hashmap_free(&seen);

    if (abandoned) {
        md_array_free(reg, alloc);
        md_array_free(con, alloc);
        pore_network_free(out);
        return false;
    }

    // Pores are the roots that survived the merges
    const size_t num_reg = md_array_size(reg);
    uint32_t* vid = (uint32_t*)md_alloc(alloc, num_reg * sizeof(uint32_t));
    for (size_t i = 0; i < num_reg; ++i) {
        vid[i] = PORE_INVALID;
        if (pore_find(reg, (uint32_t)i) != (uint32_t)i) continue;
        const pore_region_t& r = reg[i];
        pore_vertex_t pv = {};
        voxel_centre(pv.pos, field, r.peak_vox);
        pv.radius      = r.peak;
        pv.face_top    = r.face_top;
        pv.face_bottom = r.face_bot;
        pv.num_voxels  = r.count;
        vid[i] = (uint32_t)md_array_size(out->vertices);
        md_array_push(out->vertices, pv, alloc);
    }

    // Contacts resolved to pores. Two recorded between pores that were merged later are inside one
    // pore now, and two recorded before a merge can name the same pair twice; the first of those is
    // the wider, since contacts were recorded in order of decreasing clearance.
    md_array(pore_edge_t) cand = 0;
    for (size_t i = 0; i < md_array_size(con); ++i) {
        uint32_t a = vid[pore_find(reg, con[i].a)];
        uint32_t b = vid[pore_find(reg, con[i].b)];
        if (a == b) continue;
        if (a > b) { const uint32_t t = a; a = b; b = t; }
        pore_edge_t e = {};
        e.a = a;
        e.b = b;
        e.radius = con[i].d;
        voxel_centre(e.pos, field, con[i].vox);
        md_array_push(cand, e, alloc);
    }
    md_free(alloc, vid, num_reg * sizeof(uint32_t));
    md_array_free(reg, alloc);
    md_array_free(con, alloc);

    std::stable_sort(cand, cand + md_array_size(cand), [](const pore_edge_t& x, const pore_edge_t& y) {
        return (x.a != y.a) ? (x.a < y.a) : (x.b != y.b) ? (x.b < y.b) : (x.radius > y.radius);
    });
    for (size_t i = 0; i < md_array_size(cand); ++i) {
        if (i > 0 && cand[i].a == cand[i - 1].a && cand[i].b == cand[i - 1].b) continue;
        md_array_push(out->edges, cand[i], alloc);
    }
    md_array_free(cand, alloc);

    const size_t V = md_array_size(out->vertices);
    const size_t E = md_array_size(out->edges);
    std::stable_sort(out->edges, out->edges + E, [](const pore_edge_t& x, const pore_edge_t& y) { return x.radius > y.radius; });

    md_array_resize(out->adj_offset, V + 1, alloc);
    md_array_resize(out->adj, 2 * E, alloc);
    MEMSET(out->adj_offset, 0, (V + 1) * sizeof(uint32_t));
    for (size_t e = 0; e < E; ++e) {
        out->adj_offset[out->edges[e].a + 1] += 1;
        out->adj_offset[out->edges[e].b + 1] += 1;
    }
    for (size_t i = 0; i < V; ++i) out->adj_offset[i + 1] += out->adj_offset[i];
    {
        md_array(uint32_t) fill = 0;
        md_array_resize(fill, V, alloc);
        for (size_t i = 0; i < V; ++i) fill[i] = out->adj_offset[i];
        for (size_t e = 0; e < E; ++e) {
            out->adj[fill[out->edges[e].a]++] = (uint32_t)e;
            out->adj[fill[out->edges[e].b]++] = (uint32_t)e;
        }
        md_array_free(fill, alloc);
    }

    return true;
}

void pore_network_classify(uint8_t* out_class, const pore_network_t* net, double r, struct md_allocator_i* temp) {
    ASSERT(out_class);
    ASSERT(net);
    ASSERT(temp);
    const size_t V = md_array_size(net->vertices);
    if (V == 0) return;

    uint32_t* parent = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    uint8_t*  faces  = (uint8_t*) md_alloc(temp, V * sizeof(uint8_t));
    for (size_t i = 0; i < V; ++i) {
        parent[i] = (uint32_t)i;
        faces[i]  = 0;
    }

    // Widest first, so the throats a probe of radius r passes are a prefix
    for (size_t e = 0; e < md_array_size(net->edges); ++e) {
        const pore_edge_t& edge = net->edges[e];
        if ((double)edge.radius < r) break;
        const uint32_t a = vert_find(parent, edge.a);
        const uint32_t b = vert_find(parent, edge.b);
        if (a != b) parent[b] = a;
    }

    for (size_t i = 0; i < V; ++i) {
        const pore_vertex_t& v = net->vertices[i];
        uint8_t f = 0;
        if ((double)v.face_top    >= r) f |= PORE_FACE_TOP;
        if ((double)v.face_bottom >= r) f |= PORE_FACE_BOTTOM;
        faces[vert_find(parent, (uint32_t)i)] |= f;
    }

    for (size_t i = 0; i < V; ++i) {
        if ((double)net->vertices[i].radius < r) {
            out_class[i] = PORE_CLASS_SMALL;
            continue;
        }
        switch (faces[vert_find(parent, (uint32_t)i)]) {
        case PORE_FACE_TOP | PORE_FACE_BOTTOM: out_class[i] = PORE_CLASS_SPANNING; break;
        case PORE_FACE_TOP:                    out_class[i] = PORE_CLASS_TOP;      break;
        case PORE_FACE_BOTTOM:                 out_class[i] = PORE_CLASS_BOTTOM;   break;
        default:                               out_class[i] = PORE_CLASS_CLOSED;   break;
        }
    }

    md_free(temp, parent, V * sizeof(uint32_t));
    md_free(temp, faces,  V * sizeof(uint8_t));
}

void pore_network_class_sweep(uint32_t* out_counts, const pore_network_t* net, const double* r, size_t num_r, struct md_allocator_i* temp) {
    ASSERT(out_counts);
    ASSERT(net);
    ASSERT(temp);
    if (num_r == 0) return;
    ASSERT(r);
    MEMSET(out_counts, 0, num_r * PORE_CLASS_COUNT * sizeof(uint32_t));

    const size_t V = md_array_size(net->vertices);
    const size_t E = md_array_size(net->edges);
    if (V == 0) return;
    const pore_vertex_t* vert = net->vertices;

    // A component's class is which faces it reaches
    static const uint8_t face_class[4] = { PORE_CLASS_CLOSED, PORE_CLASS_TOP, PORE_CLASS_BOTTOM, PORE_CLASS_SPANNING };
    STATIC_ASSERT(PORE_FACE_TOP == 1 && PORE_FACE_BOTTOM == 2, "face_class is indexed by the face bits");

    uint32_t* parent = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    uint32_t* size   = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    uint8_t*  faces  = (uint8_t*) md_alloc(temp, V * sizeof(uint8_t));
    uint8_t*  active = (uint8_t*) md_alloc(temp, V * sizeof(uint8_t));
    uint32_t* by_rad = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    uint32_t* by_top = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    uint32_t* by_bot = (uint32_t*)md_alloc(temp, V * sizeof(uint32_t));
    for (size_t i = 0; i < V; ++i) {
        parent[i] = (uint32_t)i;
        size[i]   = 1;
        faces[i]  = 0;
        active[i] = 0;
        by_rad[i] = by_top[i] = by_bot[i] = (uint32_t)i;
    }
    std::sort(by_rad, by_rad + V, [vert](uint32_t a, uint32_t b) { return vert[a].radius      > vert[b].radius; });
    std::sort(by_top, by_top + V, [vert](uint32_t a, uint32_t b) { return vert[a].face_top    > vert[b].face_top; });
    std::sort(by_bot, by_bot + V, [vert](uint32_t a, uint32_t b) { return vert[a].face_bottom > vert[b].face_bottom; });

    uint32_t count[PORE_CLASS_COUNT] = {};
    uint32_t num_active = 0;

    // A face entry or a throat is never wider than its pore, so the pore is normally in already;
    // taking it in here keeps the counts whole if that ever fails to hold.
    auto activate = [&](uint32_t i) {
        if (active[i]) return;
        active[i] = 1;
        num_active += 1;
        count[face_class[faces[vert_find(parent, i)]]] += 1;
    };
    auto add_face = [&](uint32_t i, uint8_t bit) {
        activate(i);
        const uint32_t c = vert_find(parent, i);
        const uint8_t  f = faces[c] | bit;
        if (f == faces[c]) return;
        count[face_class[faces[c]]] -= size[c];
        count[face_class[f]]        += size[c];
        faces[c] = f;
    };
    auto join = [&](uint32_t a, uint32_t b) {
        activate(a);
        activate(b);
        a = vert_find(parent, a);
        b = vert_find(parent, b);
        if (a == b) return;
        count[face_class[faces[a]]] -= size[a];
        count[face_class[faces[b]]] -= size[b];
        if (size[a] < size[b]) { const uint32_t tmp = a; a = b; b = tmp; }
        parent[b] = a;
        size[a]  += size[b];
        faces[a] |= faces[b];
        count[face_class[faces[a]]] += size[a];
    };

    // The same comparisons classify makes - a pore fits at radius <= r, a face entry and a throat
    // pass at >= r - so the counts at each radius are its counts, not an approximation of them.
    size_t iv = 0, it = 0, ib = 0, ie = 0;
    for (size_t k = num_r; k-- > 0;) {
        const double rk = r[k];
        while (iv < V && (double)vert[by_rad[iv]].radius      >= rk) activate(by_rad[iv++]);
        while (it < V && (double)vert[by_top[it]].face_top    >= rk) add_face(by_top[it++], PORE_FACE_TOP);
        while (ib < V && (double)vert[by_bot[ib]].face_bottom >= rk) add_face(by_bot[ib++], PORE_FACE_BOTTOM);
        while (ie < E && (double)net->edges[ie].radius        >= rk) {
            join(net->edges[ie].a, net->edges[ie].b);
            ++ie;
        }

        uint32_t* out = out_counts + k * PORE_CLASS_COUNT;
        MEMCPY(out, count, sizeof(count));
        out[PORE_CLASS_SMALL] = (uint32_t)V - num_active;
    }

    md_free(temp, parent, V * sizeof(uint32_t));
    md_free(temp, size,   V * sizeof(uint32_t));
    md_free(temp, faces,  V * sizeof(uint8_t));
    md_free(temp, active, V * sizeof(uint8_t));
    md_free(temp, by_rad, V * sizeof(uint32_t));
    md_free(temp, by_top, V * sizeof(uint32_t));
    md_free(temp, by_bot, V * sizeof(uint32_t));
}

bool pore_network_widest_route(md_array(uint32_t)* out_vertices, double* out_bottleneck, const pore_network_t* net, struct md_allocator_i* alloc) {
    ASSERT(out_vertices);
    ASSERT(net);
    ASSERT(alloc);
    md_array_shrink(*out_vertices, 0);

    const uint32_t V = (uint32_t)md_array_size(net->vertices);
    if (V == 0) return false;
    const uint32_t T = V, B = V + 1;

    // Kruskal over the throats and the face entries together, widest first, until the two faces
    // meet. The links taken form a spanning forest in which the path between the faces is the
    // widest route, and the last link taken is its narrowest point.
    struct link_t { uint32_t u, w; float radius; };
    md_array(link_t) links = 0;
    for (size_t e = 0; e < md_array_size(net->edges); ++e) {
        const link_t l = { net->edges[e].a, net->edges[e].b, net->edges[e].radius };
        md_array_push(links, l, alloc);
    }
    for (uint32_t i = 0; i < V; ++i) {
        if (net->vertices[i].face_top    >= 0.0f) { const link_t l = { i, T, net->vertices[i].face_top };    md_array_push(links, l, alloc); }
        if (net->vertices[i].face_bottom >= 0.0f) { const link_t l = { i, B, net->vertices[i].face_bottom }; md_array_push(links, l, alloc); }
    }
    const size_t L = md_array_size(links);
    std::stable_sort(links, links + L, [](const link_t& x, const link_t& y) { return x.radius > y.radius; });

    md_array(uint32_t) parent = 0;
    md_array_resize(parent, V + 2, alloc);
    for (uint32_t i = 0; i < V + 2; ++i) parent[i] = i;

    md_array(link_t) tree = 0;
    bool   joined = false;
    double bottleneck = 0.0;
    for (size_t i = 0; i < L && !joined; ++i) {
        const uint32_t a = vert_find(parent, links[i].u);
        const uint32_t b = vert_find(parent, links[i].w);
        if (a == b) continue;
        parent[b] = a;
        md_array_push(tree, links[i], alloc);
        if (vert_find(parent, T) == vert_find(parent, B)) {
            joined = true;
            bottleneck = (double)links[i].radius;
        }
    }
    md_array_free(links, alloc);

    if (joined) {
        // Breadth first through the forest from the top face; the path in a forest is unique
        const uint32_t NV = V + 2;
        md_array(uint32_t) off   = 0;
        md_array(uint32_t) nbr   = 0;
        md_array(uint32_t) prev  = 0;
        md_array(uint32_t) queue = 0;
        md_array_resize(off, NV + 1, alloc);
        MEMSET(off, 0, (NV + 1) * sizeof(uint32_t));
        for (size_t i = 0; i < md_array_size(tree); ++i) {
            off[tree[i].u + 1] += 1;
            off[tree[i].w + 1] += 1;
        }
        for (uint32_t i = 0; i < NV; ++i) off[i + 1] += off[i];
        md_array_resize(nbr, off[NV], alloc);
        md_array_resize(prev, NV, alloc);
        for (uint32_t i = 0; i < NV; ++i) prev[i] = off[i];     // Fill cursor first, parent pointer after
        for (size_t i = 0; i < md_array_size(tree); ++i) {
            nbr[prev[tree[i].u]++] = tree[i].w;
            nbr[prev[tree[i].w]++] = tree[i].u;
        }
        for (uint32_t i = 0; i < NV; ++i) prev[i] = PORE_INVALID;

        prev[T] = T;
        md_array_push(queue, T, alloc);
        for (size_t qi = 0; qi < md_array_size(queue) && prev[B] == PORE_INVALID; ++qi) {
            const uint32_t u = queue[qi];
            for (uint32_t j = off[u]; j < off[u + 1]; ++j) {
                const uint32_t w = nbr[j];
                if (prev[w] != PORE_INVALID) continue;
                prev[w] = u;
                md_array_push(queue, w, alloc);
            }
        }

        // Walked back from the bottom, so reversed into top first on the way out
        md_array(uint32_t) rev = 0;
        for (uint32_t u = prev[B]; u != T && u != PORE_INVALID; u = prev[u]) md_array_push(rev, u, alloc);
        for (size_t i = md_array_size(rev); i > 0; --i) md_array_push(*out_vertices, rev[i - 1], alloc);

        md_array_free(rev, alloc);
        md_array_free(off, alloc);
        md_array_free(nbr, alloc);
        md_array_free(prev, alloc);
        md_array_free(queue, alloc);
    }

    md_array_free(tree, alloc);
    md_array_free(parent, alloc);

    if (out_bottleneck) *out_bottleneck = joined ? bottleneck : 0.0;
    return joined && md_array_size(*out_vertices) > 0;
}

uint32_t pore_network_find_edge(const pore_network_t* net, uint32_t a, uint32_t b) {
    if (!net || a >= md_array_size(net->vertices)) return PORE_INVALID;
    for (uint32_t j = net->adj_offset[a]; j < net->adj_offset[a + 1]; ++j) {
        const pore_edge_t& e = net->edges[net->adj[j]];
        if ((e.a == a && e.b == b) || (e.a == b && e.b == a)) return net->adj[j];
    }
    return PORE_INVALID;
}

void pore_network_edge_points(vec3_t* out_a, vec3_t* out_throat, vec3_t* out_b, const pore_network_t* net, uint32_t edge) {
    ASSERT(net && edge < md_array_size(net->edges));
    const pore_edge_t& e = net->edges[edge];
    const vec3_t s = { e.pos[0], e.pos[1], e.pos[2] };
    vec3_t a = { net->vertices[e.a].pos[0], net->vertices[e.a].pos[1], net->vertices[e.a].pos[2] };
    vec3_t b = { net->vertices[e.b].pos[0], net->vertices[e.b].pos[1], net->vertices[e.b].pos[2] };
    for (int ax = 0; ax < 3; ++ax) {
        if (!net->pbc[ax] || !(net->box_ext[ax] > 0.0f)) continue;
        const float L = net->box_ext[ax];
        float da = a.elem[ax] - s.elem[ax];
        float db = b.elem[ax] - s.elem[ax];
        da -= L * roundf(da / L);
        db -= L * roundf(db / L);
        a.elem[ax] = s.elem[ax] + da;
        b.elem[ax] = s.elem[ax] + db;
    }
    if (out_a)      *out_a = a;
    if (out_throat) *out_throat = s;
    if (out_b)      *out_b = b;
}
