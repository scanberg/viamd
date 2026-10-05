#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <core/md_allocator.h>

#include <math.h>
#include <algorithm>
#include <vector>

// The pore network, checked on fields where the pores and throats can be named on paper, and on a
// disordered one against channel_percolate, which reaches r_c by a different route through the same
// voxels.

namespace {

struct Field {
    std::vector<float> f;
    channel_field_t field = {};

    template <typename F>
    Field(int nx, int ny, int nz, float h, bool pbc_xy, F dist) : f((size_t)nx * ny * nz, 0.0f) {
        for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double v = dist((i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h);
            f[((size_t)k * ny + j) * nx + i] = (float)fmax(0.0, v);
        }
        field.data = f.data();
        const int dim[3] = { nx, ny, nz };
        for (int a = 0; a < 3; ++a) {
            field.dim[a]     = dim[a];
            field.spacing[a] = h;
            field.origin[a]  = 0.0f;
        }
        field.pbc[0] = field.pbc[1] = pbc_xy;
        field.pbc[2] = false;
    }
};

double dist3(double x, double y, double z, double cx, double cy, double cz) {
    return sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));
}

}  // namespace

UTEST(viamd_pores, two_cavities_and_the_throat_between_them) {
    // Two spherical cavities on a tube along z that runs from face to face. Two pores, the size of
    // the cavities, and one throat between them, the size of the tube - nothing else, since the tube
    // on either side is a ridge with no maximum of its own and has to fold into the cavity it leads to.
    md_allocator_i* alloc = md_get_heap_allocator();
    const double ax = 15.25, ay = 15.25;            // On voxel centres, so the maxima are exact
    const double z1 = 25.25, z2 = 55.25;
    const double R1 = 10.0, R2 = 7.0, rt = 2.5;

    Field fd(60, 60, 160, 0.5f, false, [&](double x, double y, double z) {
        const double rho = sqrt((x - ax) * (x - ax) + (y - ay) * (y - ay));
        return fmax(fmax(R1 - dist3(x, y, z, ax, ay, z1), R2 - dist3(x, y, z, ax, ay, z2)), rt - rho);
    });

    pore_network_t net;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 0.5, alloc, NULL, NULL));
    ASSERT_EQ((size_t)2, md_array_size(net.vertices));
    ASSERT_EQ((size_t)1, md_array_size(net.edges));

    const uint32_t big   = (net.vertices[0].radius > net.vertices[1].radius) ? 0 : 1;
    const uint32_t small = 1 - big;
    const pore_vertex_t& A = net.vertices[big];
    const pore_vertex_t& B = net.vertices[small];
    EXPECT_NEAR(R1, A.radius, 1e-5);
    EXPECT_NEAR(R2, B.radius, 1e-5);
    EXPECT_NEAR(ax, A.pos[0], 1e-4); EXPECT_NEAR(ay, A.pos[1], 1e-4); EXPECT_NEAR(z1, A.pos[2], 1e-4);
    EXPECT_NEAR(z2, B.pos[2], 1e-4);

    // The lower cavity owns the bottom face, the upper one the top, each through the tube
    EXPECT_NEAR(rt, A.face_bottom, 1e-5);
    EXPECT_LT(A.face_top, 0.0f);
    EXPECT_NEAR(rt, B.face_top, 1e-5);
    EXPECT_LT(B.face_bottom, 0.0f);

    const pore_edge_t& e = net.edges[0];
    EXPECT_NEAR(rt, e.radius, 1e-5);
    EXPECT_NEAR(ax, e.pos[0], 1e-4);
    EXPECT_TRUE(e.pos[2] > z1 && e.pos[2] < z2);
    EXPECT_EQ(0u, pore_network_find_edge(&net, big, small));
    EXPECT_EQ(0u, pore_network_find_edge(&net, small, big));

    // Every voxel the pass took in belongs to exactly one pore
    EXPECT_EQ(net.num_active, (size_t)A.num_voxels + (size_t)B.num_voxels);

    // r_c is the tube, and agrees with the percolation pass on the same voxels
    ASSERT_TRUE(net.has_r_c);
    EXPECT_NEAR(rt, net.r_c, 1e-5);
    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &fd.field, 0.25, 64, alloc));
    ASSERT_TRUE(perc.has_r_c);
    EXPECT_NEAR(perc.r_c, net.r_c, 1e-5);
    channel_percolation_free(&perc);

    // What each pore is to a probe of a given size
    uint8_t cls[2];
    pore_network_classify(cls, &net, 2.0, alloc);
    EXPECT_EQ(PORE_CLASS_SPANNING, cls[big]);
    EXPECT_EQ(PORE_CLASS_SPANNING, cls[small]);
    pore_network_classify(cls, &net, 3.0, alloc);
    EXPECT_EQ(PORE_CLASS_CLOSED, cls[big]);
    EXPECT_EQ(PORE_CLASS_CLOSED, cls[small]);
    pore_network_classify(cls, &net, 8.0, alloc);
    EXPECT_EQ(PORE_CLASS_CLOSED, cls[big]);
    EXPECT_EQ(PORE_CLASS_SMALL, cls[small]);

    // The route enters at the top, so through the upper cavity first
    md_array(uint32_t) route = 0;
    double bottleneck = 0.0;
    ASSERT_TRUE(pore_network_widest_route(&route, &bottleneck, &net, alloc));
    ASSERT_EQ((size_t)2, md_array_size(route));
    EXPECT_EQ(small, route[0]);
    EXPECT_EQ(big, route[1]);
    EXPECT_NEAR(rt, bottleneck, 1e-5);
    md_array_free(route, alloc);

    pore_network_free(&net);
}

UTEST(viamd_pores, a_shallow_dip_is_one_pore) {
    // Two overlapping cavities. Along the line between the centres the field dips to 6 between
    // maxima of 8 and 7.5, so the smaller one stands 1.5 above the saddle: two pores below that
    // persistence, one above it.
    md_allocator_i* alloc = md_get_heap_allocator();
    Field fd(60, 60, 60, 0.5f, false, [&](double x, double y, double z) {
        return fmax(8.0 - dist3(x, y, z, 12.25, 15.25, 15.25), 7.5 - dist3(x, y, z, 16.25, 15.25, 15.25));
    });

    pore_network_t net;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 1.0, alloc, NULL, NULL));
    ASSERT_EQ((size_t)2, md_array_size(net.vertices));
    ASSERT_EQ((size_t)1, md_array_size(net.edges));
    EXPECT_NEAR(6.0, net.edges[0].radius, 1e-5);
    EXPECT_FALSE(net.has_r_c);                      // Neither touches a face
    pore_network_free(&net);

    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 2.0, alloc, NULL, NULL));
    ASSERT_EQ((size_t)1, md_array_size(net.vertices));
    EXPECT_EQ((size_t)0, md_array_size(net.edges));
    EXPECT_NEAR(8.0, net.vertices[0].radius, 1e-5);
    EXPECT_NEAR(12.25, net.vertices[0].pos[0], 1e-4);
    EXPECT_EQ(net.num_active, (size_t)net.vertices[0].num_voxels);

    md_array(uint32_t) route = 0;
    EXPECT_FALSE(pore_network_widest_route(&route, nullptr, &net, alloc));
    md_array_free(route, alloc);
    pore_network_free(&net);
}

UTEST(viamd_pores, disordered_field_agrees_with_percolation) {
    // Solid spheres scattered through a box periodic in x and y. Nothing here has a closed form, so
    // the checks are the ones that have to hold whatever the geometry: the network's invariants, and
    // r_c - from the network's own union find, and from the widest route through the graph when
    // nothing is merged - against channel_percolate.
    md_allocator_i* alloc = md_get_heap_allocator();
    const int n = 48;
    const float h = 0.5f;
    const double L = n * h;

    uint32_t s = 2024;
    auto next = [&s]() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / (double)(1u << 24); };
    std::vector<double> cx, cy, cz, cr;
    for (int i = 0; i < 30; ++i) {
        cx.push_back(next() * L); cy.push_back(next() * L); cz.push_back(next() * L);
        cr.push_back(1.5 + 1.5 * next());
    }

    Field fd(n, n, n, h, true, [&](double x, double y, double z) {
        double best = 1.0e30;
        for (size_t i = 0; i < cx.size(); ++i) {
            double dx = x - cx[i]; dx -= L * round(dx / L);
            double dy = y - cy[i]; dy -= L * round(dy / L);
            const double dz = z - cz[i];
            best = fmin(best, sqrt(dx * dx + dy * dy + dz * dz) - cr[i]);
        }
        return best;
    });

    const double r_min = 0.25;
    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &fd.field, r_min, 64, alloc));
    ASSERT_TRUE(perc.has_r_c);

    pore_network_t net;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, r_min, 0.0, alloc, NULL, NULL));
    ASSERT_TRUE(net.has_r_c);
    EXPECT_NEAR(perc.r_c, net.r_c, 1e-2);
    EXPECT_TRUE(md_array_size(net.vertices) > 2);

    const size_t V = md_array_size(net.vertices);
    const size_t E = md_array_size(net.edges);

    size_t total = 0;
    for (size_t i = 0; i < V; ++i) {
        const pore_vertex_t& v = net.vertices[i];
        total += v.num_voxels;
        // The pore centre is a voxel of the field, and its clearance is the pore radius
        const int ix = (int)(v.pos[0] / h), iy = (int)(v.pos[1] / h), iz = (int)(v.pos[2] / h);
        EXPECT_NEAR((double)fd.f[((size_t)iz * n + iy) * n + ix], (double)v.radius, 1e-6);
        EXPECT_TRUE(v.face_top <= v.radius && v.face_bottom <= v.radius);
        if (net.adj_offset[i + 1] < net.adj_offset[i]) { EXPECT_TRUE(false); break; }
    }
    EXPECT_EQ(net.num_active, total);
    EXPECT_EQ(2 * E, (size_t)net.adj_offset[V]);

    for (size_t e = 0; e < E; ++e) {
        const pore_edge_t& edge = net.edges[e];
        EXPECT_TRUE(edge.a < edge.b && edge.b < V);
        // A throat is never wider than either pore it joins, and the list runs widest first
        EXPECT_TRUE(edge.radius <= net.vertices[edge.a].radius && edge.radius <= net.vertices[edge.b].radius);
        if (e > 0) EXPECT_TRUE(edge.radius <= net.edges[e - 1].radius);
        EXPECT_EQ((uint32_t)e, pore_network_find_edge(&net, edge.a, edge.b));

        // Drawn as a -> throat -> b, nearest images: no leg longer than half the box
        vec3_t pa, ps, pb;
        pore_network_edge_points(&pa, &ps, &pb, &net, (uint32_t)e);
        EXPECT_TRUE(fabsf(pa.x - ps.x) <= 0.5f * (float)L && fabsf(pb.y - ps.y) <= 0.5f * (float)L);
    }

    // Unmerged, the graph keeps every throat, so the widest route is r_c
    md_array(uint32_t) route = 0;
    double bottleneck = 0.0;
    ASSERT_TRUE(pore_network_widest_route(&route, &bottleneck, &net, alloc));
    EXPECT_NEAR(net.r_c, bottleneck, 1e-2);
    for (size_t i = 0; i + 1 < md_array_size(route); ++i) {
        const uint32_t ei = pore_network_find_edge(&net, route[i], route[i + 1]);
        ASSERT_NE(PORE_INVALID, ei);
        EXPECT_TRUE((double)net.edges[ei].radius >= bottleneck - 1e-6);
    }
    EXPECT_TRUE(net.vertices[route[0]].face_top >= (float)bottleneck - 1e-6f);
    EXPECT_TRUE(net.vertices[*md_array_last(route)].face_bottom >= (float)bottleneck - 1e-6f);
    md_array_free(route, alloc);

    // Just below r_c something spans, just above nothing does
    std::vector<uint8_t> cls(V);
    pore_network_classify(cls.data(), &net, net.r_c - 0.05, alloc);
    size_t spanning = 0;
    for (uint8_t c : cls) spanning += (c == PORE_CLASS_SPANNING);
    EXPECT_TRUE(spanning > 0);
    pore_network_classify(cls.data(), &net, net.r_c + 0.05, alloc);
    spanning = 0;
    for (uint8_t c : cls) spanning += (c == PORE_CLASS_SPANNING);
    EXPECT_EQ((size_t)0, spanning);

    // Merging only ever removes pores, and leaves r_c alone since it never came from the graph
    pore_network_t merged;
    ASSERT_TRUE(pore_network_build(&merged, &fd.field, r_min, 1.0, alloc, NULL, NULL));
    EXPECT_TRUE(md_array_size(merged.vertices) < V);
    EXPECT_EQ(net.r_c, merged.r_c);
    EXPECT_EQ(net.num_active, merged.num_active);
    pore_network_free(&merged);

    pore_network_free(&net);
    channel_percolation_free(&perc);
}

UTEST(viamd_pores, class_sweep_agrees_with_classify_at_every_radius) {
    // The sweep is one union-find over decreasing width standing in for a classification per radius,
    // so it has to give classify's counts exactly - at radii between events and at radii equal to
    // one, where the >= and < of the two have to agree. On a disordered field, merged and not.
    md_allocator_i* alloc = md_get_heap_allocator();
    const int n = 48;
    const float h = 0.5f;
    const double L = n * h;

    uint32_t s = 2024;
    auto next = [&s]() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / (double)(1u << 24); };
    std::vector<double> cx, cy, cz, cr;
    for (int i = 0; i < 30; ++i) {
        cx.push_back(next() * L); cy.push_back(next() * L); cz.push_back(next() * L);
        cr.push_back(1.5 + 1.5 * next());
    }
    Field fd(n, n, n, h, true, [&](double x, double y, double z) {
        double best = 1.0e30;
        for (size_t i = 0; i < cx.size(); ++i) {
            double dx = x - cx[i]; dx -= L * round(dx / L);
            double dy = y - cy[i]; dy -= L * round(dy / L);
            const double dz = z - cz[i];
            best = fmin(best, sqrt(dx * dx + dy * dy + dz * dz) - cr[i]);
        }
        return best;
    });

    for (double merge : { 0.0, 1.0 }) {
        pore_network_t net;
        ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, merge, alloc, NULL, NULL));
        const size_t V = md_array_size(net.vertices);
        ASSERT_TRUE(V > 2);

        // A uniform sweep past the widest pore, and every distinct event width on top of it
        std::vector<double> r;
        for (int k = 0; k <= 64; ++k) r.push_back(0.25 + 8.0 * k / 64.0);
        for (size_t i = 0; i < V; ++i) {
            r.push_back(net.vertices[i].radius);
            if (net.vertices[i].face_top    >= 0.0f) r.push_back(net.vertices[i].face_top);
            if (net.vertices[i].face_bottom >= 0.0f) r.push_back(net.vertices[i].face_bottom);
        }
        for (size_t e = 0; e < md_array_size(net.edges); ++e) r.push_back(net.edges[e].radius);
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        r.push_back(r.back() + 1.0);

        std::vector<uint32_t> counts(r.size() * PORE_CLASS_COUNT);
        pore_network_class_sweep(counts.data(), &net, r.data(), r.size(), alloc);

        std::vector<uint8_t> cls(V);
        size_t mismatches = 0, spanning_seen = 0;
        for (size_t k = 0; k < r.size(); ++k) {
            pore_network_classify(cls.data(), &net, r[k], alloc);
            uint32_t ref[PORE_CLASS_COUNT] = {};
            for (uint8_t c : cls) ref[c] += 1;
            for (int c = 0; c < PORE_CLASS_COUNT; ++c) mismatches += (ref[c] != counts[k * PORE_CLASS_COUNT + c]);
            spanning_seen += ref[PORE_CLASS_SPANNING] > 0;
        }
        EXPECT_EQ((size_t)0, mismatches);
        EXPECT_TRUE(spanning_seen > 0);     // The field is open enough that the sweep is not all zeros

        // Past the widest pore everything is too narrow
        EXPECT_EQ((uint32_t)V, counts[(r.size() - 1) * PORE_CLASS_COUNT + PORE_CLASS_SMALL]);
        pore_network_free(&net);
    }
}

UTEST(viamd_pores, a_tube_through_the_periodic_seam) {
    // The tube drifts out of one x face and back in through the other. Without the wrap it is two
    // blind pores and nothing gets through.
    md_allocator_i* alloc = md_get_heap_allocator();
    const int nx = 120, ny = 40, nz = 120;
    const float h = 0.5f;
    const double Lx = nx * h, Ly = ny * h, zmax = nz * h;

    Field fd(nx, ny, nz, h, true, [&](double x, double y, double z) {
        const double c = 50.0 + 40.0 * (z / zmax);
        double dx = x - c;     dx -= Lx * round(dx / Lx);
        double dy = y - 10.25; dy -= Ly * round(dy / Ly);
        return 4.0 - sqrt(dx * dx + dy * dy);
    });

    channel_percolation_t perc;
    ASSERT_TRUE(channel_percolate(&perc, &fd.field, 0.25, 64, alloc));
    ASSERT_TRUE(perc.has_r_c);

    pore_network_t net;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 0.5, alloc, NULL, NULL));
    ASSERT_TRUE(net.has_r_c);
    EXPECT_NEAR(perc.r_c, net.r_c, 1e-2);

    md_array(uint32_t) route = 0;
    double bottleneck = 0.0;
    EXPECT_TRUE(pore_network_widest_route(&route, &bottleneck, &net, alloc));
    md_array_free(route, alloc);
    pore_network_free(&net);
    channel_percolation_free(&perc);

    fd.field.pbc[0] = false;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 0.5, alloc, NULL, NULL));
    EXPECT_FALSE(net.has_r_c);
    pore_network_free(&net);
}

UTEST(viamd_pores, progress_is_reported_and_a_stop_is_honoured) {
    // The build runs in the background in the component, which watches it through the callback and
    // stops it through the same. Reported fractions only rise, and a stop leaves nothing behind.
    md_allocator_i* alloc = md_get_heap_allocator();
    Field fd(64, 64, 64, 0.5f, true, [&](double x, double y, double z) {
        return 3.0 - sqrt((x - 16.25) * (x - 16.25) + (y - 16.25) * (y - 16.25)) + 0.0 * z;
    });

    struct Log { std::vector<float> f; int stop_after; };
    auto cb = [](float fraction, void* user) {
        Log* l = (Log*)user;
        l->f.push_back(fraction);
        return l->stop_after < 0 || (int)l->f.size() < l->stop_after;
    };

    Log all = { {}, -1 };
    pore_network_t net;
    ASSERT_TRUE(pore_network_build(&net, &fd.field, 0.25, 0.5, alloc, cb, &all));
    ASSERT_TRUE(all.f.size() >= 2);
    for (size_t i = 1; i < all.f.size(); ++i) EXPECT_TRUE(all.f[i] >= all.f[i - 1]);
    EXPECT_TRUE(md_array_size(net.vertices) > 0);
    pore_network_free(&net);

    Log stop = { {}, 2 };
    EXPECT_FALSE(pore_network_build(&net, &fd.field, 0.25, 0.5, alloc, cb, &stop));
    EXPECT_EQ((size_t)2, stop.f.size());
    EXPECT_EQ((size_t)0, md_array_size(net.vertices));
    EXPECT_EQ((size_t)0, md_array_size(net.edges));
    EXPECT_FALSE(net.has_r_c);
}
