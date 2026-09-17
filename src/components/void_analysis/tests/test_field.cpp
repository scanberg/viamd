#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <core/md_allocator.h>
#include <core/md_coord_stream.h>
#include <core/md_grid.h>
#include <core/md_spatial_acc.h>
#include <md_types.h>

#include <float.h>
#include <math.h>
#include <vector>

// The distance field pass, which is the one stage the porosity and channel tests do not reach: they
// supply the field analytically. Here the field is computed by the pass the component runs and held
// against a brute force answer, and the profile it fills against geometry with a closed form.

namespace {

const double PI_D = 3.14159265358979323846;

// Everything one run of the pass needs, owned in one place.
struct Run {
    std::vector<uint64_t> hist, solid, total;
    void_field_accum_t accum = {};

    Run(uint32_t num_slabs, uint32_t num_bins) {
        hist.assign((size_t)num_slabs * num_bins, 0);
        solid.assign(num_slabs, 0);
        total.assign(num_slabs, 0);
        accum.hist  = hist.data();
        accum.solid = solid.data();
        accum.total = total.data();
        void_field_accum_reset(&accum, num_slabs, num_bins);
    }
};

struct Beads {
    std::vector<float> x, y, z, r;
    md_coord_stream_t stream = {};
    md_spatial_acc_t  acc = {};

    void build(const md_unitcell_t* cell) {
        stream = md_coord_stream_from_soa(x.data(), y.data(), z.data(), NULL, x.size());
        acc = {};
        acc.alloc = md_get_heap_allocator();
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &stream;
        desc.radii    = r.data();
        desc.cell_ext = 5.0;
        desc.unitcell = cell;
        md_spatial_acc_init_desc(&acc, &desc);
    }
    ~Beads() { md_spatial_acc_free(&acc); }
};

// Weighted distance by brute force, minimum image on the periodic axes of an orthorhombic box.
double brute_distance(const Beads& b, double px, double py, double pz, const double ext[3], const bool pbc[3]) {
    double best = DBL_MAX;
    for (size_t i = 0; i < b.x.size(); ++i) {
        double d[3] = { px - b.x[i], py - b.y[i], pz - b.z[i] };
        for (int a = 0; a < 3; ++a) {
            if (pbc[a]) d[a] -= ext[a] * round(d[a] / ext[a]);
        }
        const double s = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) - (double)b.r[i];
        if (s < best) best = s;
    }
    return best;
}

}  // namespace

UTEST(viamd_void_field, grid_tiles_the_box_exactly) {
    // A spacing that divides nothing. The grid must still end exactly on the cell faces, or a
    // periodic field seams where it wraps.
    const md_unitcell_t cell = md_unitcell_from_extent(103.0, 57.0, 41.0);
    md_grid_t grid = {};
    ASSERT_TRUE(void_field_grid(&grid, &cell, NULL, NULL, NULL, 0, 5.0f));

    const double ext[3] = { 103.0, 57.0, 41.0 };
    const int    dim[3] = { 21, 11, 8 };
    for (int a = 0; a < 3; ++a) {
        EXPECT_EQ(dim[a], grid.dim[a]);
        EXPECT_NEAR(0.0, (double)grid.origin.elem[a], 1e-6);
        EXPECT_NEAR(ext[a], (double)grid.spacing.elem[a] * grid.dim[a], 1e-4);
    }

    // No cell: the points plus two voxels of margin on every side
    const float x[2] = { 0.0f, 10.0f };
    const float y[2] = { 0.0f, 20.0f };
    const float z[2] = { 0.0f, 30.0f };
    ASSERT_TRUE(void_field_grid(&grid, NULL, x, y, z, 2, 1.0f));
    EXPECT_EQ(14, grid.dim[0]);
    EXPECT_EQ(24, grid.dim[1]);
    EXPECT_EQ(34, grid.dim[2]);
    EXPECT_NEAR(-2.0, (double)grid.origin.x, 1e-6);

    // Nothing to cover
    EXPECT_FALSE(void_field_grid(&grid, NULL, x, y, z, 0, 1.0f));
}

UTEST(viamd_void_field, field_matches_brute_force_and_does_not_depend_on_the_tile_split) {
    // A film-like box, periodic in x and y and open in z, whose dimensions are not multiples of the
    // tile size, so partial tiles on every far face are exercised.
    md_unitcell_t cell = md_unitcell_from_extent(37.0, 29.0, 23.0);
    cell.flags = (md_unitcell_flags_t)(MD_UNITCELL_ORTHO | MD_UNITCELL_PBC_X | MD_UNITCELL_PBC_Y);
    const double ext[3] = { 37.0, 29.0, 23.0 };
    const bool   pbc[3] = { true, true, false };

    Beads beads;
    // Deterministic scatter, including beads near the periodic faces so that images matter
    uint32_t s = 12345;
    auto next = [&s]() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / (double)(1u << 24); };
    for (int i = 0; i < 12; ++i) {
        beads.x.push_back((float)(next() * ext[0]));
        beads.y.push_back((float)(next() * ext[1]));
        beads.z.push_back((float)(2.0 + next() * (ext[2] - 4.0)));
        beads.r.push_back((float)(2.0 + 2.0 * next()));
    }
    beads.x[0] = 0.5f;  beads.y[0] = 28.5f;
    beads.build(&cell);

    md_grid_t grid = {};
    ASSERT_TRUE(void_field_grid(&grid, &cell, NULL, NULL, NULL, 0, 1.3f));
    ASSERT_NE(0, grid.dim[0] % VOID_FIELD_TILE_DIM);
    ASSERT_NE(0, grid.dim[1] % VOID_FIELD_TILE_DIM);
    ASSERT_NE(0, grid.dim[2] % VOID_FIELD_TILE_DIM);

    const size_t   num_voxels = md_grid_num_points(&grid);
    const uint32_t num_slabs  = 7;
    const uint32_t num_bins   = 64;
    std::vector<float> field(num_voxels, -1.0f);

    void_field_desc_t desc = {};
    desc.acc       = &beads.acc;
    desc.cell      = &cell;
    desc.grid      = &grid;
    desc.max_dist  = 60.0;      // Beyond anything in the box, so nothing clamps
    desc.num_slabs = num_slabs;
    desc.num_bins  = num_bins;
    desc.field     = field.data();

    const uint32_t num_tiles = void_field_num_tiles(&grid);
    ASSERT_EQ((uint32_t)(4 * 3 * 3), num_tiles);

    Run whole(num_slabs, num_bins);
    void_field_eval_tiles(&whole.accum, &desc, 0, num_tiles);

    // Every voxel written, and written with the brute force answer
    Run rebuilt(num_slabs, num_bins);
    const double bin_width = desc.max_dist / num_bins;
    double worst = 0.0;
    for (int k = 0; k < grid.dim[2]; ++k) {
        for (int j = 0; j < grid.dim[1]; ++j) {
            for (int i = 0; i < grid.dim[0]; ++i) {
                const size_t idx = ((size_t)k * grid.dim[1] + j) * grid.dim[0] + i;
                const double px = grid.origin.x + (i + 0.5) * grid.spacing.x;
                const double py = grid.origin.y + (j + 0.5) * grid.spacing.y;
                const double pz = grid.origin.z + (k + 0.5) * grid.spacing.z;
                const double ref = fmax(0.0, brute_distance(beads, px, py, pz, ext, pbc));
                worst = fmax(worst, fabs(ref - (double)field[idx]));

                // The profile rebuilt from the field, binned by the same helpers
                const uint32_t slab = void_profile_slab_of(k, grid.dim[2], num_slabs);
                rebuilt.total[slab] += 1;
                if (field[idx] <= 0.0f) {
                    rebuilt.solid[slab] += 1;
                } else {
                    rebuilt.hist[(size_t)slab * num_bins + void_profile_bin_of(field[idx], bin_width, num_bins)] += 1;
                }
            }
        }
    }
    EXPECT_LT(worst, 1e-3);
    EXPECT_TRUE(whole.hist  == rebuilt.hist);
    EXPECT_TRUE(whole.solid == rebuilt.solid);
    EXPECT_TRUE(whole.total == rebuilt.total);
    EXPECT_EQ(0u, whole.accum.num_clamped);

    // The same pass split across uneven ranges into separate accumulators, as the task system splits
    // it, merged back. Must be identical, not close.
    desc.field = NULL;
    const uint32_t cuts[] = { 0, 5, 6, 23, num_tiles };
    Run merged(num_slabs, num_bins);
    for (size_t c = 0; c + 1 < sizeof(cuts) / sizeof(cuts[0]); ++c) {
        Run part(num_slabs, num_bins);
        void_field_eval_tiles(&part.accum, &desc, cuts[c], cuts[c + 1]);
        void_field_accum_merge(&merged.accum, &part.accum, num_slabs, num_bins);
    }
    EXPECT_TRUE(merged.hist  == whole.hist);
    EXPECT_TRUE(merged.solid == whole.solid);
    EXPECT_TRUE(merged.total == whole.total);
    EXPECT_EQ(whole.accum.d_min, merged.accum.d_min);
    EXPECT_EQ(whole.accum.d_max, merged.accum.d_max);

    // A range past the end is clamped rather than read out of bounds
    Run past(num_slabs, num_bins);
    void_field_eval_tiles(&past.accum, &desc, num_tiles, num_tiles + 100);
    const void_profile_t past_prof = void_field_profile(&past.accum, &desc);
    EXPECT_EQ(0u, void_profile_num_total(&past_prof, 0, num_slabs));
}

UTEST(viamd_void_field, voxels_out_of_range_report_max_dist) {
    // One bead and no cell. Everything further than max_dist from its surface comes back at exactly
    // max_dist and is counted as clamped, in the last bin.
    Beads beads;
    beads.x = { 0.0f }; beads.y = { 0.0f }; beads.z = { 0.0f }; beads.r = { 1.0f };
    beads.build(NULL);

    md_grid_t grid = {};
    grid.orientation = mat3_ident();
    grid.origin  = { -10.0f, -10.0f, -10.0f };
    grid.spacing = { 1.0f, 1.0f, 1.0f };
    grid.dim[0] = grid.dim[1] = grid.dim[2] = 20;

    const uint32_t num_slabs = 4, num_bins = 10;
    std::vector<float> field(md_grid_num_points(&grid), -1.0f);

    void_field_desc_t desc = {};
    desc.acc = &beads.acc;
    desc.grid = &grid;
    desc.max_dist = 5.0;
    desc.num_slabs = num_slabs;
    desc.num_bins = num_bins;
    desc.field = field.data();

    Run run(num_slabs, num_bins);
    void_field_eval_tiles(&run.accum, &desc, 0, void_field_num_tiles(&grid));

    // Half integer coordinates never put a voxel centre exactly on |p| = 6
    uint64_t expect_clamped = 0, expect_solid = 0, expect_last_bin = 0;
    for (int k = 0; k < 20; ++k) for (int j = 0; j < 20; ++j) for (int i = 0; i < 20; ++i) {
        const double px = -9.5 + i, py = -9.5 + j, pz = -9.5 + k;
        const double d = sqrt(px * px + py * py + pz * pz) - 1.0;
        if (d >= 5.0) expect_clamped += 1;
        if (d <= 0.0) expect_solid += 1;
        if (d >= 4.5) expect_last_bin += 1;
        const float f = field[((size_t)k * 20 + j) * 20 + i];
        if (d >= 5.0) { EXPECT_NEAR(5.0, (double)f, 1e-4); }
    }

    EXPECT_EQ(expect_clamped, run.accum.num_clamped);
    EXPECT_NEAR(5.0, (double)run.accum.d_max, 1e-4);

    uint64_t solid = 0, last_bin = 0;
    for (uint32_t sl = 0; sl < num_slabs; ++sl) {
        solid    += run.solid[sl];
        last_bin += run.hist[(size_t)sl * num_bins + num_bins - 1];
    }
    EXPECT_EQ(expect_solid, solid);
    EXPECT_EQ(expect_last_bin, last_bin);
}

UTEST(viamd_void_field, triclinic_cell_is_counted_once) {
    // The grid covers the bounding box of a triclinic cell, so its corners are periodic images of
    // voxels already inside. The statistics must add up to the volume of the cell, not of the box.
    const md_unitcell_t cell = md_unitcell_from_basis_parameters(40.0, 36.0, 30.0, 12.0, -8.0, 6.0);
    ASSERT_TRUE(md_unitcell_is_triclinic(&cell));
    const double cell_volume = 40.0 * 36.0 * 30.0;

    Beads beads;
    beads.x = { 20.0f }; beads.y = { 18.0f }; beads.z = { 15.0f }; beads.r = { 3.0f };
    beads.build(&cell);

    md_grid_t grid = {};
    ASSERT_TRUE(void_field_grid(&grid, &cell, NULL, NULL, NULL, 0, 0.5f));

    const uint32_t num_slabs = 16, num_bins = 64;
    void_field_desc_t desc = {};
    desc.acc = &beads.acc;
    desc.cell = &cell;
    desc.grid = &grid;
    desc.max_dist = 40.0;
    desc.num_slabs = num_slabs;
    desc.num_bins = num_bins;

    Run run(num_slabs, num_bins);
    void_field_eval_tiles(&run.accum, &desc, 0, void_field_num_tiles(&grid));

    const void_profile_t prof = void_field_profile(&run.accum, &desc);
    const double counted = (double)void_profile_num_total(&prof, 0, num_slabs) * prof.voxel_volume;
    const double boxed   = (double)md_grid_num_points(&grid) * prof.voxel_volume;

    EXPECT_GT(boxed, 1.2 * cell_volume);                    // The mask has something to do
    EXPECT_NEAR(1.0, counted / cell_volume, 5e-3);          // and does it
}

UTEST(viamd_void_field, porosity_of_a_sphere_in_a_periodic_cube) {
    // End to end: bead, pass, profile, reduction. One sphere of radius a in a periodic cube of side L
    // has porosity 1 - (4/3) pi a^3 / L^3, and the volume a probe of radius R can occupy is the same
    // with a + R, as long as the grown sphere still fits.
    const double L = 40.0, a = 12.0, R = 4.0;
    const md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    Beads beads;
    beads.x = { 20.0f }; beads.y = { 20.0f }; beads.z = { 20.0f }; beads.r = { (float)a };
    beads.build(&cell);

    md_grid_t grid = {};
    ASSERT_TRUE(void_field_grid(&grid, &cell, NULL, NULL, NULL, 0, 0.5f));

    const uint32_t num_slabs = 32, num_bins = 512;
    void_field_desc_t desc = {};
    desc.acc = &beads.acc;
    desc.cell = &cell;
    desc.grid = &grid;
    desc.max_dist = L;
    desc.num_slabs = num_slabs;
    desc.num_bins = num_bins;

    Run run(num_slabs, num_bins);
    void_field_eval_tiles(&run.accum, &desc, 0, void_field_num_tiles(&grid));
    const void_profile_t prof = void_field_profile(&run.accum, &desc);
    ASSERT_TRUE(void_profile_valid(&prof));

    const double V = L * L * L;
    const double phi_ref = 1.0 - (4.0 / 3.0) * PI_D * a * a * a / V;
    const double acc_ref = 1.0 - (4.0 / 3.0) * PI_D * (a + R) * (a + R) * (a + R) / V;

    EXPECT_NEAR(phi_ref, void_profile_porosity(&prof, 0, num_slabs), 2e-3);
    EXPECT_NEAR(acc_ref, void_profile_accessible_fraction(&prof, 0, num_slabs, R), 2e-3);

    // The deepest point of the void is the cube corner, at half the body diagonal from the centre
    EXPECT_NEAR(sqrt(3.0) * 0.5 * L - a, (double)run.accum.d_max, 0.5);
}
