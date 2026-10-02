#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <core/md_allocator.h>
#include <core/md_coord_stream.h>
#include <core/md_grid.h>
#include <core/md_spatial_acc.h>
#include <md_types.h>

#include <math.h>
#include <vector>

// Surface topography and the mass profile, against geometry with a closed form answer.

namespace {

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
        desc.cutoff   = 5.0;
        desc.unitcell = cell;
        md_spatial_acc_init(&acc, &desc);
    }
    ~Beads() { md_spatial_acc_free(&acc); }

    void_beads_t view() const { return void_beads(&acc, r.data(), r.size()); }
};

md_grid_t make_grid(int nx, int ny, int nz, float h) {
    md_grid_t g = {};
    g.orientation = mat3_ident();
    g.origin  = vec3_t{0.0f, 0.0f, 0.0f};
    g.spacing = vec3_t{h, h, h};
    g.dim[0] = nx;
    g.dim[1] = ny;
    g.dim[2] = nz;
    return g;
}

struct Maps {
    std::vector<float> top, bot;
    Maps(const md_grid_t& g, const void_beads_t& beads, double probe, double max_dist, uint32_t split) {
        const size_t n = (size_t)g.dim[0] * g.dim[1];
        top.assign(n, 12345.0f);
        bot.assign(n, 12345.0f);
        void_heightmap_desc_t d = {};
        d.beads = beads;
        d.grid = &g;
        d.probe_radius = probe;
        d.max_dist = max_dist;
        const uint32_t np = void_heightmap_num_patches(&g);
        // Evaluate in uneven chunks: the result must not depend on how the patches were handed out
        for (uint32_t b = 0; b < np; b += split) {
            void_heightmap_eval_patches(top.data(), bot.data(), &d, b, b + split);
        }
    }
};

}  // namespace

UTEST(viamd_void_heightmap, sphere_is_seen_dilated_by_the_probe) {
    // One bead of radius Rs in an open box. A probe of radius R lowered onto it stops with its centre
    // at |p - c| = Rs + R, so the apex sits at cz + sqrt((Rs + R)^2 - rho^2) - R, and every column
    // further out than Rs + R falls straight through.
    const double Rs = 10.0;
    const double c[3] = { 30.0, 30.0, 30.0 };

    Beads beads;
    beads.x = { (float)c[0] };
    beads.y = { (float)c[1] };
    beads.z = { (float)c[2] };
    beads.r = { (float)Rs };
    beads.build(NULL);

    const md_grid_t g = make_grid(60, 60, 60, 1.0f);

    for (double R : { 0.0, 3.0 }) {
        Maps m(g, beads.view(), R, 100.0, 5);
        const double Rt = Rs + R;
        int num_hit = 0, num_miss = 0;
        for (int j = 0; j < g.dim[1]; ++j) {
            for (int i = 0; i < g.dim[0]; ++i) {
                const size_t k = (size_t)j * g.dim[0] + i;
                const double dx = (i + 0.5) - c[0];
                const double dy = (j + 0.5) - c[1];
                const double rho2 = dx * dx + dy * dy;
                if (rho2 < Rt * Rt - 1.0e-3) {
                    const double s = sqrt(Rt * Rt - rho2);
                    EXPECT_NEAR(c[2] + s - R, (double)m.top[k], 2.0e-2);
                    EXPECT_NEAR(c[2] - s + R, (double)m.bot[k], 2.0e-2);
                    num_hit += 1;
                } else if (rho2 > Rt * Rt + 1.0e-3) {
                    EXPECT_TRUE(isnan(m.top[k]));
                    EXPECT_TRUE(isnan(m.bot[k]));
                    num_miss += 1;
                }
            }
        }
        EXPECT_GT(num_hit, 100);
        EXPECT_GT(num_miss, 100);
    }
}

UTEST(viamd_void_heightmap, probe_bridges_a_gap_narrower_than_its_diameter) {
    // A square layer of beads at one height, with a gap between neighbours. A point probe drops
    // through the middle of four beads; one wider than the clearance there rests on them.
    const double Rb = 4.0, pitch = 10.0, zc = 20.0;     // Gap between surfaces: 2 A
    Beads beads;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            beads.x.push_back((float)(5.0 + n * pitch));
            beads.y.push_back((float)(5.0 + m * pitch));
            beads.z.push_back((float)zc);
            beads.r.push_back((float)Rb);
        }
    }
    beads.build(NULL);

    // Columns through the middle of a gap between four beads, where the clearance is largest
    const md_grid_t g = make_grid(80, 80, 80, 0.5f);
    const int gi = (int)(10.0 / 0.5);    // x = 10.25, between the beads at 5 and 15
    const size_t k = (size_t)gi * g.dim[0] + gi;

    {
        Maps m(g, beads.view(), 0.0, 100.0, 7);
        EXPECT_TRUE(isnan(m.top[k]));
    }
    {
        // Clearance in the middle of four beads is sqrt(2) * 5 - 4 = 3.07 A, so a 3.5 A probe rests
        // on them: centre at |p - c| = 7.5 with rho = sqrt(2) * 5.25 from the nearest bead.
        const double R = 3.5;
        Maps m(g, beads.view(), R, 100.0, 7);
        ASSERT_FALSE(isnan(m.top[k]));
        const double x = (gi + 0.5) * 0.5;
        const double dx = x - 5.0, dy = x - 5.0;
        // The nearest bead centre to (x, x) is at (5, 5) or (15, 15); x = 10.25 is nearer (15, 15)
        const double ex = 15.0 - x, ey = 15.0 - x;
        const double rho2 = MIN(dx * dx + dy * dy, ex * ex + ey * ey);
        const double expect = zc + sqrt((Rb + R) * (Rb + R) - rho2) - R;
        EXPECT_NEAR(expect, (double)m.top[k], 2.0e-2);
    }
}

UTEST(viamd_void_heightmap, stats_skip_open_columns) {
    const float h[6] = { 1.0f, 3.0f, NAN, 1.0f, 3.0f, NAN };
    void_heightmap_stats_t s;
    void_heightmap_stats(&s, h, 6);
    EXPECT_EQ(4u, (unsigned)s.num_valid);
    EXPECT_EQ(2u, (unsigned)s.num_open);
    EXPECT_NEAR(2.0, s.mean, 1e-12);
    EXPECT_NEAR(1.0, s.min,  1e-12);
    EXPECT_NEAR(3.0, s.max,  1e-12);
    EXPECT_NEAR(1.0, s.rq,   1e-12);
    EXPECT_NEAR(1.0, s.ra,   1e-12);
}

UTEST(viamd_void_heightmap, mass_profile_bins_by_slab_and_wraps_only_when_periodic) {
    uint64_t hist[4 * 2] = {};
    uint64_t solid[4] = {};
    uint64_t total[4] = { 100, 100, 100, 100 };
    void_profile_t p = {};
    p.hist = hist; p.solid = solid; p.total = total;
    p.num_slabs = 4; p.num_bins = 2; p.bin_width = 1.0;
    p.z_min = 10.0; p.z_max = 30.0; p.slab_height = 5.0; p.voxel_volume = 0.5;

    const vec3_t xyz[5] = { {0, 0, 11.0f}, {0, 0, 16.0f}, {0, 0, 29.9f}, {0, 0, 31.0f}, {0, 0, 8.0f} };   // last two outside [10, 30)
    const float mass[5] = { 1.0f, 2.0f, 4.0f, 8.0f, 16.0f };
    double m[4];

    EXPECT_NEAR(7.0, void_profile_bin_mass(m, &p, xyz, mass, 5, false), 1e-12);
    EXPECT_NEAR(1.0, m[0], 1e-12);
    EXPECT_NEAR(2.0, m[1], 1e-12);
    EXPECT_NEAR(0.0, m[2], 1e-12);
    EXPECT_NEAR(4.0, m[3], 1e-12);

    // Periodic: 31 wraps to 11 (slab 0), 8 wraps to 28 (slab 3)
    EXPECT_NEAR(31.0, void_profile_bin_mass(m, &p, xyz, mass, 5, true), 1e-12);
    EXPECT_NEAR(9.0,  m[0], 1e-12);
    EXPECT_NEAR(20.0, m[3], 1e-12);

    // Density over slabs [0, 2): 11 / (200 * 0.5)
    EXPECT_NEAR(11.0 / 100.0, void_profile_mass_density(&p, m, 0, 2), 1e-12);

    // Number density with no masses
    EXPECT_NEAR(5.0, void_profile_bin_mass(m, &p, xyz, NULL, 5, true), 1e-12);
}
