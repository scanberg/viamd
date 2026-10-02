#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <core/md_allocator.h>

#include <math.h>
#include <vector>

// Channels through a structure along z, checked against geometry with a closed form answer. A
// disordered network cannot tell a bug from physics, so every case here is one that can be worked
// out on paper.

namespace {

const double PI_D = 3.14159265358979323846;

}  // namespace

UTEST(viamd_channels, drilled_holes_count_and_critical_radius) {
    // A solid block with three cylindrical holes along z, radii 3, 5 and 7. The channel count is
    // then a step function of the probe radius and the critical radius is the widest hole.
    md_allocator_i* alloc = md_get_heap_allocator();

    const float h = 0.5f;
    const int dim[3] = { 180, 180, 120 };
    const double ax[3] = { 20.25, 45.25, 70.25 };   // on voxel centres, so the maxima are exact
    const double ay     = 45.25;
    const double rad[3] = { 3.0, 5.0, 7.0 };

    std::vector<float> f((size_t)dim[0]*dim[1]*dim[2], 0.0f);
    for (int k = 0; k < dim[2]; ++k)
    for (int j = 0; j < dim[1]; ++j)
    for (int i = 0; i < dim[0]; ++i) {
        const double x = (i + 0.5) * h, y = (j + 0.5) * h;
        double best = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double rho = sqrt((x - ax[c])*(x - ax[c]) + (y - ay)*(y - ay));
            best = fmax(best, rad[c] - rho);
        }
        f[(size_t)k*dim[0]*dim[1] + (size_t)j*dim[0] + i] = (float)fmax(0.0, best);
    }

    channel_field_t field = {};
    field.data = f.data();
    for (int a = 0; a < 3; ++a) {
        field.dim[a] = dim[a];
        field.spacing[a] = h;
        field.origin[a] = 0.0f;
        field.pbc[a] = false;
    }

    const double probes[4] = { 2.0, 4.0, 6.0, 8.0 };
    const uint32_t want[4] = { 3, 2, 1, 0 };
    uint32_t counts[4];
    channel_spanning_counts(counts, probes, 4, &field, alloc);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(want[i], counts[i]);
    }

    EXPECT_NEAR(7.0, channel_critical_radius(&field, 0.5, 12.0, 0.01, alloc), 0.05);

    channel_tree_t tree;
    ASSERT_TRUE(channel_sweep(&tree, &field, 2.0, true, alloc));
    EXPECT_EQ(3u, tree.num_spanning);

    int spanning_roots = 0;
    for (size_t i = 0; i < md_array_size(tree.roots); ++i) {
        const channel_node_t& n = tree.nodes[tree.roots[i]];
        if (!(n.reaches_top && n.reaches_bottom)) continue;
        spanning_roots += 1;
        // A drilled hole is one unbranched channel running the full depth
        EXPECT_EQ(CHANNEL_INVALID_INDEX, n.first_child);
        EXPECT_EQ((size_t)dim[2], md_array_size(n.path));
    }
    EXPECT_EQ(3, spanning_roots);
    channel_tree_free(&tree);
}

UTEST(viamd_channels, hex_array_gives_the_throat_not_the_cavity) {
    // The case from docs/cnf-void-analysis-design.md, turned so the sweep crosses the cylinder axes:
    // a hexagonal array of parallel cylinders, radius R, lattice constant a.
    //   throat between neighbouring cylinders : a/2 - R           = 5.0
    //   interstitial cavity radius            : a/sqrt(3) - R     = 6.55
    // The cavity is strictly the larger. A pipeline that reports the cavity as the percolation
    // radius has the pore size / percolation confusion baked in, and this is where it shows.
    md_allocator_i* alloc = md_get_heap_allocator();

    const double R = 5.0, a = 20.0;
    const double row_dz = a * sqrt(3.0) / 2.0;
    const double throat = a / 2.0 - R;
    const double cavity = a / sqrt(3.0) - R;

    const float h = 0.25f;
    const double Ly = 4 * a;          // periodic
    const double Lz = 4 * row_dz;     // open, and the sweep axis
    const int dim[3] = { 16, (int)(Ly / h), (int)(Lz / h) };

    std::vector<float> f((size_t)dim[0]*dim[1]*dim[2], 0.0f);
    double max_clear = 0.0;
    for (int k = 0; k < dim[2]; ++k) {
        const double z = (k + 0.5) * h;
        for (int j = 0; j < dim[1]; ++j) {
            const double y = (j + 0.5) * h;
            double dmin = 1.0e30;
            for (int m = 0; m < 4; ++m) {
                const double zc = row_dz * 0.5 + m * row_dz;
                for (int n = 0; n < 4; ++n) {
                    const double yc = n * a + ((m & 1) ? 0.5 * a : 0.0);
                    double dy = y - yc;
                    dy -= Ly * round(dy / Ly);
                    const double dz = z - zc;
                    dmin = fmin(dmin, sqrt(dy*dy + dz*dz));
                }
            }
            const double clear = fmax(0.0, dmin - R);
            // Only inside the lattice: beyond the outermost rows there is open space, whose
            // clearance says nothing about the interstitial cavity
            if (z > row_dz * 0.5 && z < row_dz * 3.5) max_clear = fmax(max_clear, clear);
            for (int i = 0; i < dim[0]; ++i) {
                f[(size_t)k*dim[0]*dim[1] + (size_t)j*dim[0] + i] = (float)clear;
            }
        }
    }

    channel_field_t field = {};
    field.data = f.data();
    for (int a2 = 0; a2 < 3; ++a2) {
        field.dim[a2] = dim[a2];
        field.spacing[a2] = h;
        field.origin[a2] = 0.0f;
    }
    field.pbc[0] = true; field.pbc[1] = true; field.pbc[2] = false;

    // Half a voxel diagonal, since the cavity centre does not land on a voxel centre
    EXPECT_NEAR(cavity, max_clear, 0.2);

    const double rc = channel_critical_radius(&field, 0.1, 10.0, 0.005, alloc);
    EXPECT_NEAR(throat, rc, 0.2);
    EXPECT_LT(rc, cavity - 1.0);
}

UTEST(viamd_channels, converging_tubes_make_one_channel_with_two_branches) {
    md_allocator_i* alloc = md_get_heap_allocator();

    const float h = 0.5f;
    const int dim[3] = { 180, 120, 120 };
    const double R = 4.0, spread = 10.0, probe = 2.0;
    const double zmax = dim[2] * h;

    std::vector<float> f((size_t)dim[0]*dim[1]*dim[2], 0.0f);
    for (int k = 0; k < dim[2]; ++k) {
        const double z = (k + 0.5) * h;
        const double off = spread * (z / zmax);
        for (int j = 0; j < dim[1]; ++j) {
            const double y = (j + 0.5) * h;
            for (int i = 0; i < dim[0]; ++i) {
                const double x = (i + 0.5) * h;
                const double da = R - sqrt((x - (45.0 - off))*(x - (45.0 - off)) + (y - 30.0)*(y - 30.0));
                const double db = R - sqrt((x - (45.0 + off))*(x - (45.0 + off)) + (y - 30.0)*(y - 30.0));
                f[(size_t)k*dim[0]*dim[1] + (size_t)j*dim[0] + i] = (float)fmax(0.0, fmax(da, db));
            }
        }
    }

    channel_field_t field = {};
    field.data = f.data();
    for (int a = 0; a < 3; ++a) {
        field.dim[a] = dim[a];
        field.spacing[a] = h;
        field.origin[a] = 0.0f;
        field.pbc[a] = false;
    }

    channel_tree_t tree;
    ASSERT_TRUE(channel_sweep(&tree, &field, probe, true, alloc));
    EXPECT_EQ(1u, tree.num_spanning);

    int spanning_roots = 0, branches_at_top = 0;
    double junction = 0.0;
    for (size_t i = 0; i < md_array_size(tree.roots); ++i) {
        const channel_node_t& n = tree.nodes[tree.roots[i]];
        if (!(n.reaches_top && n.reaches_bottom)) continue;
        spanning_roots += 1;
        junction = n.z_top;
        for (uint32_t c = n.first_child; c != CHANNEL_INVALID_INDEX; c = tree.nodes[c].next_sibling) {
            if (tree.nodes[c].reaches_top) branches_at_top += 1;
        }
    }
    EXPECT_EQ(1, spanning_roots);
    EXPECT_EQ(2, branches_at_top);

    // At probe radius r the accessible tubes have radius R - r, so they touch where the axis
    // separation reaches 2(R - r), that is at z = (R - r) * zmax / spread.
    //
    // The sweep is allowed to report the junction slightly above that, and does: once the gap
    // between the two cross sections is narrower than a voxel, no voxel centre lands inside it and
    // the grid cannot tell them apart, which closes the gap early by up to h/2 of separation. That
    // is the resolution limit of the whole approach, so it is pinned rather than tolerated - a
    // junction outside this window means something other than discretization.
    const double z_exact = (R - probe) * zmax / spread;
    const double z_grid  = (R - probe + 0.5 * h) * zmax / spread;
    EXPECT_GE(junction, z_exact - 0.5);
    EXPECT_LE(junction, z_grid + 0.5);

    channel_tree_free(&tree);
}

// --- Percolation -------------------------------------------------------------------------------
//
// The same questions the sweeps above answer by repetition, answered once in order of decreasing
// clearance. The cross checks below are the point: two algorithms with almost nothing in common
// have to land on the same critical radius, and if they do not, one of them is wrong.

namespace {

// A solid block with three cylindrical holes along z, radii 3, 5 and 7 - the same geometry as the
// first case above, shared so the two families are checked against the identical field.
struct DrilledBlock {
    static const int DIM[3];
    static constexpr float H = 0.5f;
    static constexpr double AX[3] = { 20.25, 45.25, 70.25 };
    static constexpr double AY = 45.25;
    static constexpr double RAD[3] = { 3.0, 5.0, 7.0 };

    std::vector<float> f;
    channel_field_t field = {};

    DrilledBlock() : f((size_t)DIM[0]*DIM[1]*DIM[2], 0.0f) {
        for (int k = 0; k < DIM[2]; ++k)
        for (int j = 0; j < DIM[1]; ++j)
        for (int i = 0; i < DIM[0]; ++i) {
            const double x = (i + 0.5) * H, y = (j + 0.5) * H;
            double best = 0.0;
            for (int c = 0; c < 3; ++c) {
                const double rho = sqrt((x - AX[c])*(x - AX[c]) + (y - AY)*(y - AY));
                best = fmax(best, RAD[c] - rho);
            }
            f[(size_t)k*DIM[0]*DIM[1] + (size_t)j*DIM[0] + i] = (float)fmax(0.0, best);
        }
        field.data = f.data();
        for (int a = 0; a < 3; ++a) {
            field.dim[a] = DIM[a];
            field.spacing[a] = H;
            field.origin[a] = 0.0f;
            field.pbc[a] = false;
        }
    }
};
const int DrilledBlock::DIM[3] = { 180, 180, 120 };
constexpr double DrilledBlock::AX[3];
constexpr double DrilledBlock::RAD[3];

}  // namespace

UTEST(viamd_channels, percolation_finds_the_widest_hole_and_where_it_is) {
    // The critical radius is the widest hole, and - unlike a bisection, which returns a number and
    // nothing else - the pass also says which voxel it belongs to. That voxel has to sit on the axis
    // of the widest hole, because that is the only place with clearance 7.
    md_allocator_i* alloc = md_get_heap_allocator();
    DrilledBlock b;

    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &b.field, 0.5, 64, alloc));

    ASSERT_TRUE(perc.has_r_c);
    EXPECT_NEAR(7.0, perc.r_c, 0.02);
    EXPECT_NEAR(DrilledBlock::AX[2], perc.throat[0], 0.6);
    EXPECT_NEAR(DrilledBlock::AY,    perc.throat[1], 0.6);

    // Agreement with the bisection, which asks a threshold question about the same field and shares
    // no code with this.
    const double rc_bisect = channel_critical_radius(&b.field, 0.5, 12.0, 0.005, alloc);
    EXPECT_NEAR(rc_bisect, perc.r_c, 0.02);

    // The count curve is the step function the sweep reports, recovered from the same pass.
    for (size_t i = 0; i < md_array_size(perc.radius); ++i) {
        const double r = perc.radius[i];
        uint32_t want = 0;
        if (r <= 3.0) want = 3; else if (r <= 5.0) want = 2; else if (r <= 7.0) want = 1;
        // Only away from the radii themselves, where a voxel either side is a coin toss
        if (fabs(r - 3.0) < 0.2 || fabs(r - 5.0) < 0.2 || fabs(r - 7.0) < 0.2) continue;
        EXPECT_EQ(want, perc.num_spanning[i]);
    }

    channel_percolation_free(&perc);
}

UTEST(viamd_channels, percolation_curves_are_ordered_and_monotone) {
    // Properties that hold whatever the geometry. A spanning component is open and an open one is
    // void, so the three fractions nest; all three can only grow as the probe shrinks; and the two
    // penetration depths are what the spanning flag is made of, so they must meet exactly where it
    // turns on.
    md_allocator_i* alloc = md_get_heap_allocator();
    DrilledBlock b;

    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &b.field, 0.5, 96, alloc));

    const size_t n = md_array_size(perc.radius);
    ASSERT_TRUE(n > 8);

    for (size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(perc.frac_void[i] >= perc.frac_open[i] - 1.0e-12);
        EXPECT_TRUE(perc.frac_open[i] >= perc.frac_spanning[i] - 1.0e-12);
        EXPECT_TRUE(perc.num_spanning[i] <= perc.num_components[i]);
        if (i > 0) {
            EXPECT_TRUE(perc.radius[i] > perc.radius[i - 1]);
            // Ascending in radius, so the sets shrink
            EXPECT_TRUE(perc.frac_void[i] <= perc.frac_void[i - 1] + 1.0e-12);
            EXPECT_TRUE(perc.frac_open[i] <= perc.frac_open[i - 1] + 1.0e-12);
            EXPECT_TRUE(perc.z_from_top[i]    >= perc.z_from_top[i - 1]    - 1.0e-9);
            EXPECT_TRUE(perc.z_from_bottom[i] <= perc.z_from_bottom[i - 1] + 1.0e-9);
        }
        // Spanning is exactly "the two penetration fronts have met"
        const bool fronts_met = perc.z_from_top[i] <= perc.z_from_bottom[i] + 1.0e-9;
        EXPECT_EQ(fronts_met, perc.num_spanning[i] > 0);
        EXPECT_EQ(perc.radius[i] <= perc.r_c, perc.num_spanning[i] > 0);
    }

    channel_percolation_free(&perc);
}

UTEST(viamd_channels, percolation_separates_a_sealed_cavity_from_the_open_pore) {
    // One tube through the block plus a sealed spherical cavity that touches nothing. The cavity is
    // wider than the tube, so an accessible-volume measure counts it first and a probe can never
    // reach it. Open porosity is what tells them apart, and this is the flood fill the design
    // document lists as owed - it falls out of the ordering rather than needing a second pass.
    md_allocator_i* alloc = md_get_heap_allocator();

    const float h = 0.5f;
    const int dim[3] = { 120, 120, 120 };
    const double tube_x = 30.25, tube_y = 30.25, tube_r = 4.0;
    const double cav_x = 45.25, cav_y = 45.25, cav_z = 30.25, cav_r = 6.0;

    std::vector<float> f((size_t)dim[0]*dim[1]*dim[2], 0.0f);
    for (int k = 0; k < dim[2]; ++k) {
        const double z = (k + 0.5) * h;
        for (int j = 0; j < dim[1]; ++j) {
            const double y = (j + 0.5) * h;
            for (int i = 0; i < dim[0]; ++i) {
                const double x = (i + 0.5) * h;
                const double tube = tube_r - sqrt((x-tube_x)*(x-tube_x) + (y-tube_y)*(y-tube_y));
                const double cav  = cav_r  - sqrt((x-cav_x)*(x-cav_x) + (y-cav_y)*(y-cav_y) + (z-cav_z)*(z-cav_z));
                f[(size_t)k*dim[0]*dim[1] + (size_t)j*dim[0] + i] = (float)fmax(0.0, fmax(tube, cav));
            }
        }
    }

    channel_field_t field = {};
    field.data = f.data();
    for (int a = 0; a < 3; ++a) {
        field.dim[a] = dim[a]; field.spacing[a] = h; field.origin[a] = 0.0f; field.pbc[a] = false;
    }

    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &field, 0.25, 128, alloc));

    // The cavity is the widest thing in the box, but nothing gets through above the tube radius.
    ASSERT_TRUE(perc.has_r_c);
    EXPECT_NEAR(tube_r, perc.r_c, 0.2);
    EXPECT_TRUE(perc.r_c < cav_r - 1.0);

    // Between the tube radius and the cavity radius the only void left is the sealed cavity: void
    // volume with nothing open and nothing spanning.
    bool checked = false;
    for (size_t i = 0; i < md_array_size(perc.radius); ++i) {
        const double r = perc.radius[i];
        if (r < tube_r + 0.5 || r > cav_r - 0.5) continue;
        EXPECT_TRUE(perc.frac_void[i] > 0.0);
        EXPECT_NEAR(0.0, perc.frac_open[i], 1.0e-12);
        EXPECT_EQ(0u, perc.num_spanning[i]);
        checked = true;
    }
    EXPECT_TRUE(checked);

    // And below the tube radius the cavity is still closed: the gap between void and open is it.
    const double probe = 2.0;
    size_t at = 0;
    for (size_t i = 0; i < md_array_size(perc.radius); ++i) if (perc.radius[i] <= probe) at = i;
    EXPECT_TRUE(perc.frac_void[at] - perc.frac_open[at] > 0.0);

    channel_percolation_free(&perc);
}

UTEST(viamd_channels, percolation_crosses_the_periodic_seam) {
    // The box is periodic in x and the only tube through it drifts far enough to leave one x face
    // and come back through the other. If the wrap is not honoured the tube is two blind pores and
    // nothing gets through - so this pins the periodic neighbour handling in the percolation pass. It is the same wrap md_spatial_acc applies when it builds the field
    // in the first place; here it is asserted rather than assumed.
    md_allocator_i* alloc = md_get_heap_allocator();

    const float h = 0.5f;
    const int dim[3] = { 120, 40, 120 };
    const double Lx = dim[0] * h, Ly = dim[1] * h;
    const double R = 4.0, x0 = 50.0, drift = 40.0;
    const double zmax = dim[2] * h;

    std::vector<float> f((size_t)dim[0]*dim[1]*dim[2], 0.0f);
    for (int k = 0; k < dim[2]; ++k) {
        const double z = (k + 0.5) * h;
        const double cx = x0 + drift * (z / zmax);       // 50 -> 90, i.e. straight across the seam at 60
        for (int j = 0; j < dim[1]; ++j) {
            const double y = (j + 0.5) * h;
            for (int i = 0; i < dim[0]; ++i) {
                const double x = (i + 0.5) * h;
                double dx = x - cx;    dx -= Lx * round(dx / Lx);
                double dy = y - 10.25; dy -= Ly * round(dy / Ly);
                f[(size_t)k*dim[0]*dim[1] + (size_t)j*dim[0] + i] = (float)fmax(0.0, R - sqrt(dx*dx + dy*dy));
            }
        }
    }

    channel_field_t field = {};
    field.data = f.data();
    for (int a = 0; a < 3; ++a) { field.dim[a] = dim[a]; field.spacing[a] = h; field.origin[a] = 0.0f; }
    field.pbc[0] = true; field.pbc[1] = true; field.pbc[2] = false;

    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &field, 0.25, 64, alloc));
    ASSERT_TRUE(perc.has_r_c);

    // Comfortably inside the tube, so the wrap is genuinely being followed, and never wider than the
    // tube itself. It measures ~3.58 rather than 4, and that is the grid rather than an error: the
    // axis drifts a third of a voxel in x per step in z, so a six connected route cannot follow it
    // without stepping sideways onto a voxel further from the axis than the one it left. The same
    // discretization that closes narrow gaps early - biasing r_c up - biases it down whenever the
    // route is not axis aligned, and a tube through a real network is never axis aligned.
    EXPECT_TRUE(perc.r_c > 3.0);
    EXPECT_TRUE(perc.r_c <= R + 0.05);

    channel_percolation_free(&perc);

    // Turn the wrap off and the same field is cut at the seam: nothing through.
    field.pbc[0] = false;
    channel_percolation_t flat;
    ASSERT_TRUE(channel_percolate(&flat, &field, 0.25, 64, alloc));
    EXPECT_FALSE(flat.has_r_c);
    channel_percolation_free(&flat);
}
