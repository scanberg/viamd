#include "utest.h"

#include <channels.h>

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

UTEST(viamd_channels, layout_columns_a_small_tree) {
    // One merge node with two leaves below it: the leaves take columns 0 and 1, and the merge sits
    // between them at 0.5.
    md_allocator_i* alloc = md_get_heap_allocator();

    channel_tree_t tree = {};
    tree.alloc = alloc;

    channel_node_t n = {};
    n.parent = n.first_child = n.next_sibling = CHANNEL_INVALID_INDEX;

    channel_node_t root = n;   // 0
    channel_node_t a    = n;   // 1
    channel_node_t b    = n;   // 2
    root.reaches_top = root.reaches_bottom = true;
    a.reaches_top = b.reaches_top = true;
    root.first_child = 1;
    a.next_sibling   = 2;
    a.parent = b.parent = 0;

    md_array_push(tree.nodes, root, alloc);
    md_array_push(tree.nodes, a, alloc);
    md_array_push(tree.nodes, b, alloc);
    md_array_push(tree.roots, 0u, alloc);

    float slot[3];
    const float used = channel_tree_layout(slot, &tree, alloc);

    EXPECT_NEAR(2.0, (double)used, 1.0e-6);
    EXPECT_NEAR(0.0, (double)slot[1], 1.0e-6);
    EXPECT_NEAR(1.0, (double)slot[2], 1.0e-6);
    EXPECT_NEAR(0.5, (double)slot[0], 1.0e-6);

    channel_tree_free(&tree);
}

UTEST(viamd_channels, layout_survives_a_deep_tree) {
    // A merge tree carries one node per pair of branches that join, so a dense network produces a
    // chain many thousands of levels deep. Laying that out by recursing over the children overflows
    // the stack on a real system, so the traversal carries its own.
    //
    // The tree is built directly rather than swept out of a field: reaching this depth through a
    // field would need a grid with as many slabs, and the depth is the whole point of the case.
    md_allocator_i* alloc = md_get_heap_allocator();

    // Deep enough to overrun a default stack on every platform the project builds for: recursing
    // here costs roughly 50 bytes a level, so this is tens of megabytes of frames against the 8 MB
    // a linux thread gets and the 1 MB a windows one does. The window between those two is exactly
    // where a version that only crashes for some users lives.
    const uint32_t DEPTH = 250000;

    channel_tree_t tree = {};
    tree.alloc = alloc;

    // Level i is a merge node at index 2i whose children are level i+1 and a leaf at index 2i+1.
    // The deepest level has only its leaf.
    for (uint32_t i = 0; i < DEPTH; ++i) {
        channel_node_t merge = {};
        merge.parent = (i == 0) ? CHANNEL_INVALID_INDEX : 2 * (i - 1);
        merge.next_sibling = CHANNEL_INVALID_INDEX;
        merge.first_child = (i + 1 < DEPTH) ? 2 * (i + 1) : (2 * i + 1);
        merge.reaches_top = true;
        merge.reaches_bottom = (i == 0);
        md_array_push(tree.nodes, merge, alloc);

        channel_node_t leaf = {};
        leaf.parent = 2 * i;
        leaf.first_child = CHANNEL_INVALID_INDEX;
        leaf.next_sibling = CHANNEL_INVALID_INDEX;
        leaf.reaches_top = true;
        md_array_push(tree.nodes, leaf, alloc);
    }
    // Wire the leaf of level i as the sibling of level i+1
    for (uint32_t i = 0; i + 1 < DEPTH; ++i) {
        tree.nodes[2 * (i + 1)].next_sibling = 2 * i + 1;
    }
    md_array_push(tree.roots, 0u, alloc);

    md_array(float) slot = 0;
    md_array_resize(slot, md_array_size(tree.nodes), alloc);
    const float used = channel_tree_layout(slot, &tree, alloc);

    // One column per leaf, and every leaf got a distinct one
    EXPECT_NEAR((double)DEPTH, (double)used, 1.0e-6);
    md_array(uint8_t) taken = 0;
    md_array_resize(taken, DEPTH, alloc);
    MEMSET(taken, 0, DEPTH);
    uint32_t bad = 0;
    for (uint32_t i = 0; i < DEPTH; ++i) {
        const float s = slot[2 * i + 1];
        if (!(s >= 0.0f && s < (float)DEPTH)) { bad += 1; continue; }
        const uint32_t c = (uint32_t)s;
        if (taken[c]) bad += 1;
        taken[c] = 1;
    }
    EXPECT_EQ(0u, bad);

    // A merge node sits at the mean of its children
    for (uint32_t i = 0; i + 1 < DEPTH; ++i) {
        const double want = 0.5 * ((double)slot[2 * (i + 1)] + (double)slot[2 * i + 1]);
        if (fabs(want - (double)slot[2 * i]) > 1.0e-3) {
            EXPECT_NEAR(want, (double)slot[2 * i], 1.0e-3);
            break;
        }
    }

    md_array_free(slot, alloc);
    md_array_free(taken, alloc);
    channel_tree_free(&tree);
}
