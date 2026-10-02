#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <core/md_allocator.h>
#include <core/md_coord_stream.h>
#include <core/md_spatial_acc.h>
#include <md_types.h>
#include <md_unitcell.h>

#include <float.h>
#include <math.h>
#include <vector>

// The weighted nearest query everything in this component stands on, held against brute force over
// the cases which have broken nearest queries before: images across periodic faces, a skewed cell,
// a search which outgrows the period, points outside the cell, an open axis, and coordinates far
// from the origin.

namespace {

struct Rng {
    uint32_t s;
    double next() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / (double)(1u << 24); }
};

struct Case {
    md_unitcell_t cell = {};
    bool has_cell = false;
    double A[3][3] = {};
    bool pbc[3] = {};
    std::vector<float> x, y, z, r;
    md_coord_stream_t stream = {};
    md_spatial_acc_t acc = {};

    void build() {
        stream = md_coord_stream_from_soa(x.data(), y.data(), z.data(), NULL, x.size());
        acc = {};
        acc.alloc = md_get_heap_allocator();
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &stream;
        desc.cutoff   = 5.0;
        desc.unitcell = has_cell ? &cell : NULL;
        md_spatial_acc_init(&acc, &desc);
        if (has_cell) {
            md_unitcell_A_extract_double(A, &cell);
            const uint32_t f = md_unitcell_flags(&cell);
            pbc[0] = (f & MD_UNITCELL_PBC_X) != 0;
            pbc[1] = (f & MD_UNITCELL_PBC_Y) != 0;
            pbc[2] = (f & MD_UNITCELL_PBC_Z) != 0;
        }
    }
    ~Case() { md_spatial_acc_free(&acc); }

    // Every image within two periods, which is far more than any of these cells need
    void brute(double* out_d, uint32_t* out_i, double* out_second, double px, double py, double pz, double max_dist) const {
        double best = max_dist, second = max_dist;
        uint32_t bi = VOID_BEAD_NONE;
        const int n[3] = { pbc[0] ? 2 : 0, pbc[1] ? 2 : 0, pbc[2] ? 2 : 0 };
        for (size_t i = 0; i < x.size(); ++i) {
            double bead_best = DBL_MAX;
            for (int k = -n[2]; k <= n[2]; ++k) for (int j = -n[1]; j <= n[1]; ++j) for (int l = -n[0]; l <= n[0]; ++l) {
                const double cx = x[i] + A[0][0] * l + A[1][0] * j + A[2][0] * k;
                const double cy = y[i] + A[0][1] * l + A[1][1] * j + A[2][1] * k;
                const double cz = z[i] + A[0][2] * l + A[1][2] * j + A[2][2] * k;
                const double d = sqrt((px - cx) * (px - cx) + (py - cy) * (py - cy) + (pz - cz) * (pz - cz)) - r[i];
                bead_best = fmin(bead_best, d);
            }
            if (bead_best < best) { second = best; best = bead_best; bi = (uint32_t)i; }
            else if (bead_best < second) { second = bead_best; }
        }
        *out_d = best;
        *out_i = bi;
        *out_second = second;
    }
};

// Query points: compact blocks the way the field and the topography supply them, and a scatter.
// Counts which are not multiples of a vector, a group or a batch.
void make_points(std::vector<float>& qx, std::vector<float>& qy, std::vector<float>& qz, Rng& rng, const double lo[3], const double hi[3]) {
    qx.clear(); qy.clear(); qz.clear();
    for (int blk = 0; blk < 6; ++blk) {
        const double h = 0.3 + 1.2 * rng.next();
        const int nb[3] = { 3 + (int)(rng.next() * 9), 3 + (int)(rng.next() * 9), 1 + (int)(rng.next() * 9) };
        double o[3];
        for (int a = 0; a < 3; ++a) o[a] = lo[a] + rng.next() * (hi[a] - lo[a]);
        for (int k = 0; k < nb[2]; ++k) for (int j = 0; j < nb[1]; ++j) for (int i = 0; i < nb[0]; ++i) {
            qx.push_back((float)(o[0] + i * h));
            qy.push_back((float)(o[1] + j * h));
            qz.push_back((float)(o[2] + k * h));
        }
    }
    for (int i = 0; i < 333; ++i) {
        qx.push_back((float)(lo[0] + rng.next() * (hi[0] - lo[0])));
        qy.push_back((float)(lo[1] + rng.next() * (hi[1] - lo[1])));
        qz.push_back((float)(lo[2] + rng.next() * (hi[2] - lo[2])));
    }
}

// Runs the query over the points, one call each for a few batch sizes, and compares with brute force.

void check(int* utest_result, const Case& c, const std::vector<float>& qx, const std::vector<float>& qy, const std::vector<float>& qz, double max_dist, double tol) {
    const void_beads_t beads = void_beads(&c.acc, c.r.data(), c.r.size());
    const size_t n = qx.size();
    std::vector<float> dist(n);
    std::vector<uint32_t> idx(n);

    std::vector<double>   ref_d(n), ref_second(n);
    std::vector<uint32_t> ref_i(n);
    for (size_t i = 0; i < n; ++i) {
        c.brute(&ref_d[i], &ref_i[i], &ref_second[i], qx[i], qy[i], qz[i], max_dist);
    }

    size_t num_bad = 0, num_idx_bad = 0;
    for (size_t chunk : { (size_t)7, (size_t)64, (size_t)513, n }) {
        for (size_t b = 0; b < n; b += chunk) {
            const size_t m = MIN(chunk, n - b);
            void_beads_nearest(&beads, qx.data() + b, qy.data() + b, qz.data() + b, m, max_dist, idx.data() + b, dist.data() + b);
        }
        for (size_t i = 0; i < n; ++i) {
            if (fabs(ref_d[i] - (double)dist[i]) > tol) {
                if (num_bad < 5) printf("  point %zu (%g %g %g): got %.6f (bead %u), expected %.6f (bead %u)\n", i, qx[i], qy[i], qz[i], dist[i], idx[i], ref_d[i], ref_i[i]);
                num_bad += 1;
            }
            // The bead only where the answer is not a near tie
            if (ref_second[i] - ref_d[i] > 2.0 * tol && idx[i] != ref_i[i]) num_idx_bad += 1;
        }
    }
    EXPECT_EQ((size_t)0, num_bad);
    EXPECT_EQ((size_t)0, num_idx_bad);
}

void scatter_beads(Case& c, Rng& rng, size_t num, const double lo[3], const double hi[3], double r_lo, double r_hi) {
    for (size_t i = 0; i < num; ++i) {
        c.x.push_back((float)(lo[0] + rng.next() * (hi[0] - lo[0])));
        c.y.push_back((float)(lo[1] + rng.next() * (hi[1] - lo[1])));
        c.z.push_back((float)(lo[2] + rng.next() * (hi[2] - lo[2])));
        c.r.push_back((float)(r_lo + rng.next() * (r_hi - r_lo)));
    }
}

}  // namespace

UTEST(viamd_void_nearest, open_and_periodic_boxes_against_brute_force) {
    // Small max_dist (searches which stop at once), medium, and one larger than the box, which takes the
    // search past the period. Beads include zero radius ones and a cluster near a corner.
    const double max_dists[] = { 3.0, 12.0, 80.0 };
    const double ext[3] = { 31.0, 26.0, 22.0 };

    for (int mode = 0; mode < 3; ++mode) {
        for (double max_dist : max_dists) {
            Rng rng = { 777u + (uint32_t)mode * 31u + (uint32_t)max_dist };
            Case c;
            if (mode > 0) {
                c.cell = md_unitcell_from_extent(ext[0], ext[1], ext[2]);
                if (mode == 2) c.cell.flags = (md_unitcell_flags_t)(MD_UNITCELL_ORTHO | MD_UNITCELL_PBC_X | MD_UNITCELL_PBC_Y);
                c.has_cell = true;
            }
            const double lo[3] = { 0.0, 0.0, 0.0 };
            scatter_beads(c, rng, 150, lo, ext, 0.5, 3.0);
            const double clo[3] = { 0.0, 0.0, 0.0 }, chi[3] = { 3.0, 3.0, 3.0 };
            scatter_beads(c, rng, 30, clo, chi, 0.0, 1.5);
            c.build();

            // Points reaching outside the box on every side
            const double plo[3] = { -6.0, -6.0, -6.0 };
            const double phi[3] = { ext[0] + 6.0, ext[1] + 6.0, ext[2] + 6.0 };
            std::vector<float> qx, qy, qz;
            make_points(qx, qy, qz, rng, plo, phi);
            check(utest_result, c, qx, qy, qz, max_dist, 1e-3);
        }
    }
}

UTEST(viamd_void_nearest, triclinic_cell_against_brute_force) {
    const double max_dists[] = { 4.0, 15.0, 70.0 };
    for (double max_dist : max_dists) {
        Rng rng = { 4242u + (uint32_t)max_dist };
        Case c;
        c.cell = md_unitcell_from_basis_parameters(34.0, 29.0, 25.0, 9.0, -7.0, 5.0);
        c.has_cell = true;
        ASSERT_TRUE(md_unitcell_is_triclinic(&c.cell));
        // Beads anywhere in the cell, by fractional coordinate
        double A[3][3];
        md_unitcell_A_extract_double(A, &c.cell);
        for (int i = 0; i < 200; ++i) {
            const double f[3] = { rng.next(), rng.next(), rng.next() };
            c.x.push_back((float)(A[0][0] * f[0] + A[1][0] * f[1] + A[2][0] * f[2]));
            c.y.push_back((float)(A[0][1] * f[0] + A[1][1] * f[1] + A[2][1] * f[2]));
            c.z.push_back((float)(A[0][2] * f[0] + A[1][2] * f[1] + A[2][2] * f[2]));
            c.r.push_back((float)(0.5 + 2.5 * rng.next()));
        }
        c.build();

        const double plo[3] = { -8.0, -8.0, -8.0 };
        const double phi[3] = { 50.0, 40.0, 33.0 };
        std::vector<float> qx, qy, qz;
        make_points(qx, qy, qz, rng, plo, phi);
        check(utest_result, c, qx, qy, qz, max_dist, 1e-3);
    }
}

UTEST(viamd_void_nearest, sparse_beads_in_a_large_open_space) {
    // Few beads, far apart, with points in the open between them: the search has to grow over several
    // rounds before it reaches anything, and stops at max_dist where there is nothing.
    for (int mode = 0; mode < 2; ++mode) {
        Rng rng = { 99u + (uint32_t)mode };
        Case c;
        const double ext[3] = { 400.0, 300.0, 250.0 };
        if (mode == 1) {
            c.cell = md_unitcell_from_extent(ext[0], ext[1], ext[2]);
            c.has_cell = true;
        }
        const double lo[3] = { 0.0, 0.0, 0.0 };
        scatter_beads(c, rng, 40, lo, ext, 1.0, 6.0);
        c.build();
        std::vector<float> qx, qy, qz;
        make_points(qx, qy, qz, rng, lo, ext);
        check(utest_result, c, qx, qy, qz, 90.0, 2e-3);
        check(utest_result, c, qx, qy, qz, 1000.0, 2e-3);
    }
}

UTEST(viamd_void_nearest, far_from_the_origin) {
    // A 20 um periodic box with the beads in one corner region and points around them, near the far face
    // and across it. The distances are relative to the batch, so they keep their precision out here.
    Rng rng = { 5150u };
    Case c;
    const double L[3] = { 20000.0, 20000.0, 400.0 };
    c.cell = md_unitcell_from_extent(L[0], L[1], L[2]);
    c.has_cell = true;
    const double lo[3] = { L[0] - 30.0, L[1] - 30.0, 0.0 };
    const double hi[3] = { L[0],        L[1],        40.0 };
    scatter_beads(c, rng, 300, lo, hi, 1.0, 2.5);
    c.build();
    const double plo[3] = { L[0] - 35.0, L[1] - 35.0, -5.0 };
    const double phi[3] = { L[0] + 10.0, L[1] + 10.0, 45.0 };
    std::vector<float> qx, qy, qz;
    make_points(qx, qy, qz, rng, plo, phi);
    // Floats carry ~2e-3 A at 2e4 A, which is what the query points themselves are good to
    check(utest_result, c, qx, qy, qz, 25.0, 5e-3);
}

UTEST(viamd_void_nearest, empty_and_out_of_range) {
    Case c;
    c.x = { 0.0f }; c.y = { 0.0f }; c.z = { 0.0f }; c.r = { 1.0f };
    c.build();
    const void_beads_t beads = void_beads(&c.acc, c.r.data(), c.r.size());

    const float qx[3] = { 0.0f, 3.0f, 100.0f }, qy[3] = {}, qz[3] = {};
    float d[3];
    uint32_t i[3];
    void_beads_nearest(&beads, qx, qy, qz, 3, 10.0, i, d);
    EXPECT_NEAR(-1.0, (double)d[0], 1e-6);
    EXPECT_EQ(0u, i[0]);
    EXPECT_NEAR(2.0, (double)d[1], 1e-6);
    EXPECT_EQ(0u, i[1]);
    EXPECT_EQ(10.0f, d[2]);
    EXPECT_EQ(VOID_BEAD_NONE, i[2]);

    // No beads at all
    void_beads_t none = {};
    void_beads_nearest(&none, qx, qy, qz, 3, 7.0, i, d);
    EXPECT_EQ(7.0f, d[0]);
    EXPECT_EQ(VOID_BEAD_NONE, i[0]);
}
