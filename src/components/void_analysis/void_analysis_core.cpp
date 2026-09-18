#include "void_analysis_core.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_grid.h>
#include <core/md_spatial_acc.h>
#include <core/md_coord_stream.h>
#include <md_types.h>

#include <float.h>
#include <math.h>
#include <string.h>

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

// =================================================================================================
// The distance field
// =================================================================================================

double void_profile_bin_mass(double* out_slab_mass, const void_profile_t* prof, const float* z, const float* mass, size_t count, bool periodic_z) {
    if (!out_slab_mass || !prof || prof->num_slabs == 0 || !(prof->slab_height > 0.0)) return 0.0;
    memset(out_slab_mass, 0, prof->num_slabs * sizeof(double));
    if (!z) return 0.0;

    const double H = prof->slab_height * (double)prof->num_slabs;
    double total = 0.0;
    for (size_t i = 0; i < count; ++i) {
        double u = (double)z[i] - prof->z_min;
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

bool void_field_grid(md_grid_t* out_grid, const md_unitcell_t* cell, const float* x, const float* y, const float* z, size_t count, float spacing) {
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
        if (count == 0 || !x || !y || !z) return false;
        vec3_t aabb_min = { FLT_MAX,  FLT_MAX,  FLT_MAX};
        vec3_t aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (size_t i = 0; i < count; ++i) {
            const vec3_t p = {x[i], y[i], z[i]};
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
    if (!accum || !desc || !desc->acc || !desc->grid) return;
    if (desc->num_slabs == 0 || desc->num_bins == 0) return;

    const md_grid_t& grid = *desc->grid;
    const int TD = VOID_FIELD_TILE_DIM;
    const int tiles[3] = {
        (grid.dim[0] + TD - 1) / TD,
        (grid.dim[1] + TD - 1) / TD,
        (grid.dim[2] + TD - 1) / TD,
    };
    tile_end = MIN(tile_end, void_field_num_tiles(desc->grid));

    const uint32_t num_slabs = desc->num_slabs;
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

        md_coord_stream_t pts = md_coord_stream_from_soa(qx, qy, qz, NULL, (size_t)n);
        md_spatial_acc_query_nearest(desc->acc, &pts, desc->max_dist, NULL, dist);

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

            const uint32_t slab = void_profile_slab_of(vk[p], grid.dim[2], num_slabs);

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
    if (!accum || !desc || !desc->grid || desc->num_slabs == 0 || desc->num_bins == 0) return p;
    const md_grid_t& grid = *desc->grid;
    p.hist         = accum->hist;
    p.solid        = accum->solid;
    p.total        = accum->total;
    p.num_slabs    = desc->num_slabs;
    p.num_bins     = desc->num_bins;
    p.bin_width    = desc->max_dist / (double)desc->num_bins;
    p.z_min        = (double)grid.origin.z;
    p.slab_height  = (double)grid.spacing.z * (double)grid.dim[2] / (double)desc->num_slabs;
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
        md_coord_stream_t pts = md_coord_stream_from_soa(qx, qy, qz, NULL, (size_t)num_act);
        md_spatial_acc_query_nearest(desc->acc, &pts, desc->max_dist, NULL, dist);

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
    if (!desc || !desc->acc || !desc->grid) return;
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

// ---------------------------------------------------------------------------------------------
// Tracing a route
// ---------------------------------------------------------------------------------------------

namespace {

const int PATH_DIR[6][3] = { {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1} };
const uint8_t PATH_SEED  = 7;

// Voxel index of an integer lattice point, wrapping x and y where the field is periodic. Returns
// false for a point off a non periodic face.
inline bool path_voxel(const channel_field_t* f, long ix, long iy, long ik, size_t* out) {
    const long nx = f->dim[0], ny = f->dim[1], nz = f->dim[2];
    if (f->pbc[0]) { ix = ((ix % nx) + nx) % nx; } else if (ix < 0 || ix >= nx) return false;
    if (f->pbc[1]) { iy = ((iy % ny) + ny) % ny; } else if (iy < 0 || iy >= ny) return false;
    if (ik < 0 || ik >= nz) return false;
    *out = ((size_t)ik * (size_t)ny + (size_t)iy) * (size_t)nx + (size_t)ix;
    return true;
}

// Does the straight segment between two lattice points stay inside {d >= r}? Sampled at half a
// voxel, which is the finest statement the field itself supports.
bool path_segment_clear(const channel_field_t* f, float r, const double a[3], const double b[3]) {
    const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
    const double len = sqrt(dx * dx + dy * dy + dz * dz);
    const int steps = (int)(2.0 * len) + 1;
    for (int i = 0; i <= steps; ++i) {
        const double t = (double)i / (double)steps;
        size_t v;
        if (!path_voxel(f, lround(a[0] + t * dx), lround(a[1] + t * dy), lround(a[2] + t * dz), &v)) return false;
        if (f->data[v] < r) return false;
    }
    return true;
}

// One sample of the route: the wrapped world position, with the clearance the field actually
// reports there rather than an interpolation between two distant anchors.
void path_emit(md_array(vec4_t)* out, const channel_field_t* f, double sx, double sy, double sz, struct md_allocator_i* alloc) {
    size_t v = 0;
    if (!path_voxel(f, lround(sx), lround(sy), lround(sz), &v)) return;
    double wx = sx, wy = sy;
    if (f->pbc[0]) { wx = fmod(sx, (double)f->dim[0]); if (wx < 0.0) wx += (double)f->dim[0]; }
    if (f->pbc[1]) { wy = fmod(sy, (double)f->dim[1]); if (wy < 0.0) wy += (double)f->dim[1]; }
    vec4_t pt;
    pt.x = f->origin[0] + (float)(wx + 0.5) * f->spacing[0];
    pt.y = f->origin[1] + (float)(wy + 0.5) * f->spacing[1];
    pt.z = f->origin[2] + (float)(sz + 0.5) * f->spacing[2];
    pt.w = f->data[v];
    md_array_push(*out, pt, alloc);
}

}  // namespace

bool channel_trace_path(md_array(vec4_t)* out_path, double* out_length, const channel_field_t* field, double r, struct md_allocator_i* alloc) {
    ASSERT(out_path);
    ASSERT(field);
    ASSERT(alloc);

    if (!field->data) return false;
    if (field->pbc[2]) return false;
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    const size_t plane = (size_t)nx * (size_t)ny;
    const size_t N     = plane * (size_t)nz;
    const float  rr    = (float)r;

    // One byte per voxel: which way we arrived. A parent index would be four, and the direction is
    // all a backtrace needs.
    uint8_t* from = (uint8_t*)md_alloc(alloc, N);
    if (!from) return false;
    MEMSET(from, 0, N);

    md_array(uint32_t) queue = 0;
    size_t head = 0;
    size_t found = N;

    for (size_t i = 0; i < plane; ++i) {
        const size_t v = (size_t)(nz - 1) * plane + i;
        if (field->data[v] < rr) continue;
        from[v] = PATH_SEED;
        if (nz == 1) { found = v; break; }
        md_array_push(queue, (uint32_t)v, alloc);
    }

    while (found == N && head < md_array_size(queue)) {
        const size_t v = (size_t)queue[head++];
        const int k = (int)(v / plane);
        const int y = (int)((v % plane) / (size_t)nx);
        const int x = (int)(v % (size_t)nx);

        for (int d = 0; d < 6; ++d) {
            size_t w;
            if (!path_voxel(field, (long)x + PATH_DIR[d][0], (long)y + PATH_DIR[d][1], (long)k + PATH_DIR[d][2], &w)) continue;
            if (from[w] || field->data[w] < rr) continue;
            from[w] = (uint8_t)(d + 1);
            if (w < plane) { found = w; break; }        // Reached the bottom face
            md_array_push(queue, (uint32_t)w, alloc);
        }
    }

    md_array_free(queue, alloc);

    if (found == N) {
        md_free(alloc, from, N);
        return false;
    }

    // Backtrace, then reverse, carrying unwrapped lattice coordinates so a route that leaves through
    // a periodic face keeps going in a straight line instead of jumping the box.
    md_array(int32_t) chain = 0;                        // Directions, bottom to top
    {
        size_t v = found;
        while (from[v] != PATH_SEED) {
            const int d = (int)from[v] - 1;
            md_array_push(chain, (int32_t)d, alloc);
            size_t w;
            const int k = (int)(v / plane);
            const int y = (int)((v % plane) / (size_t)nx);
            const int x = (int)(v % (size_t)nx);
            if (!path_voxel(field, (long)x - PATH_DIR[d][0], (long)y - PATH_DIR[d][1], (long)k - PATH_DIR[d][2], &w)) break;
            v = w;
        }
        // v is now the seed at the top face; rebuild the route downward from it
        md_array(double) pts = 0;
        long ix = (long)(v % (size_t)nx);
        long iy = (long)((v % plane) / (size_t)nx);
        long ik = (long)(v / plane);
        md_array_push(pts, (double)ix, alloc);
        md_array_push(pts, (double)iy, alloc);
        md_array_push(pts, (double)ik, alloc);
        for (size_t i = md_array_size(chain); i > 0; --i) {
            const int d = (int)chain[i - 1];
            ix += PATH_DIR[d][0];
            iy += PATH_DIR[d][1];
            ik += PATH_DIR[d][2];
            md_array_push(pts, (double)ix, alloc);
            md_array_push(pts, (double)iy, alloc);
            md_array_push(pts, (double)ik, alloc);
        }
        md_array_free(chain, alloc);
        md_free(alloc, from, N);

        const size_t n = md_array_size(pts) / 3;

        // Straighten it. Breadth first through a six connected grid returns a staircase: it is the
        // shortest route in voxel steps, which overstates the length of anything not axis aligned by
        // up to sqrt(3). Replacing a run by the straight segment between its ends, where that segment
        // stays inside the set, recovers a length worth quoting and a route worth looking at.
        const size_t LOOKAHEAD = 64;
        md_array(double) anchor = 0;
        for (size_t i = 0; ; ) {
            md_array_push(anchor, pts[i*3+0], alloc);
            md_array_push(anchor, pts[i*3+1], alloc);
            md_array_push(anchor, pts[i*3+2], alloc);
            if (i + 1 >= n) break;

            size_t best = i + 1;
            const size_t far = MIN(n - 1, i + LOOKAHEAD);
            for (size_t j = far; j > i + 1; --j) {
                const double a[3] = { pts[i*3+0], pts[i*3+1], pts[i*3+2] };
                const double b[3] = { pts[j*3+0], pts[j*3+1], pts[j*3+2] };
                if (path_segment_clear(field, rr, a, b)) { best = j; break; }
            }
            i = best;
        }
        md_array_free(pts, alloc);

        const size_t na = md_array_size(anchor) / 3;

        double length = 0.0;
        for (size_t i = 1; i < na; ++i) {
            const double dx = (anchor[i*3+0] - anchor[(i-1)*3+0]) * (double)field->spacing[0];
            const double dy = (anchor[i*3+1] - anchor[(i-1)*3+1]) * (double)field->spacing[1];
            const double dz = (anchor[i*3+2] - anchor[(i-1)*3+2]) * (double)field->spacing[2];
            length += sqrt(dx * dx + dy * dy + dz * dz);
        }
        if (out_length) *out_length = length;

        // Resample the straightened route at half a voxel, reading the clearance at each sample.
        //
        // The anchors alone are not something to draw with. Straightening is free to leave them
        // sixty voxels apart, and the clearance between two of them is whatever the field says, not
        // the interpolation of its endpoints - so anything drawn from the anchors is at its widest
        // exactly where its width was never measured. Sampling also keeps every step short, which is
        // what lets a caller drop the one segment that crosses a periodic face without leaving a
        // visible gap in the rest.
        md_array_shrink(*out_path, 0);
        for (size_t i = 0; i + 1 < na; ++i) {
            const double ax = anchor[i*3+0],     ay = anchor[i*3+1],     az = anchor[i*3+2];
            const double dx = anchor[(i+1)*3+0] - ax;
            const double dy = anchor[(i+1)*3+1] - ay;
            const double dz = anchor[(i+1)*3+2] - az;
            const int steps = MAX(1, (int)(2.0 * sqrt(dx * dx + dy * dy + dz * dz)));
            for (int t = 0; t < steps; ++t) {
                const double u = (double)t / (double)steps;
                path_emit(out_path, field, ax + u * dx, ay + u * dy, az + u * dz, alloc);
            }
        }
        path_emit(out_path, field, anchor[(na-1)*3+0], anchor[(na-1)*3+1], anchor[(na-1)*3+2], alloc);
        md_array_free(anchor, alloc);
    }

    return md_array_size(*out_path) > 1;
}

float channel_tree_layout(float* out_slot, const channel_tree_t* tree, struct md_allocator_i* alloc) {
    ASSERT(out_slot);
    ASSERT(tree);

    const size_t num_nodes = md_array_size(tree->nodes);
    for (size_t i = 0; i < num_nodes; ++i) out_slot[i] = -1.0f;

    // An explicit traversal stack. Post order: a node's column is only known once every child has
    // one, so a frame stays on the stack accumulating its children's columns until they are done.
    typedef struct {
        uint32_t node;
        uint32_t next_child;
        float    sum;
        uint32_t count;
    } frame_t;

    md_array(frame_t) stack = 0;
    float slot = 0.0f;

    for (size_t i = 0; i < md_array_size(tree->roots); ++i) {
        const uint32_t root = tree->roots[i];
        if (!(tree->nodes[root].reaches_top && tree->nodes[root].reaches_bottom)) continue;

        md_array_shrink(stack, 0);
        frame_t root_frame = { root, tree->nodes[root].first_child, 0.0f, 0 };
        md_array_push(stack, root_frame, alloc);

        while (md_array_size(stack) > 0) {
            frame_t* top = md_array_last(stack);

            if (top->next_child != CHANNEL_INVALID_INDEX) {
                const uint32_t child = top->next_child;
                top->next_child = tree->nodes[child].next_sibling;
                frame_t child_frame = { child, tree->nodes[child].first_child, 0.0f, 0 };
                // The push may move the array, so nothing may be read through top after this
                md_array_push(stack, child_frame, alloc);
                continue;
            }

            const float column = top->count ? (top->sum / (float)top->count) : slot++;
            out_slot[top->node] = column;
            md_array_pop(stack);

            if (md_array_size(stack) > 0) {
                frame_t* parent = md_array_last(stack);
                parent->sum   += column;
                parent->count += 1;
            }
        }
    }

    md_array_free(stack, alloc);
    return slot;
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
