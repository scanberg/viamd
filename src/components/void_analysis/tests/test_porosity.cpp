#include "utest.h"

#include <void_analysis/void_analysis_core.h>

#include <math.h>
#include <vector>

// Porosity and accessible volume, checked against geometry with a closed form answer. A disordered
// network cannot tell a bug from physics, so every case here is one that can be worked out on paper.
//
// The histogram is filled here with the same void_profile_bin_of and void_profile_slab_of the field
// pass uses, so the binning convention is under test rather than restated. What is not under test is
// the distance field itself: these cases supply d analytically, which is the point - a reduction
// checked against a field checked against nothing proves only that they agree.

namespace {

const double PI_D = 3.14159265358979323846;

// A profile built by sampling an analytic distance function on a grid, binned exactly the way the
// field pass bins it.
struct Sampled {
    std::vector<uint64_t> hist, solid, total;
    void_profile_t prof = {};

    template <typename F>
    Sampled(int nx, int ny, int nz, double h, uint32_t num_slabs, uint32_t num_bins, double max_dist, F dist) {
        const double bin_width = max_dist / (double)num_bins;

        hist.assign((size_t)num_slabs * num_bins, 0);
        solid.assign(num_slabs, 0);
        total.assign(num_slabs, 0);

        for (int k = 0; k < nz; ++k) {
            const uint32_t slab = void_profile_slab_of(k, nz, num_slabs);
            const double z = ((double)k + 0.5) * h;
            for (int j = 0; j < ny; ++j) {
                const double y = ((double)j + 0.5) * h;
                for (int i = 0; i < nx; ++i) {
                    const double x = ((double)i + 0.5) * h;
                    double d = dist(x, y, z);
                    if (d > max_dist) d = max_dist;
                    total[slab] += 1;
                    if (d <= 0.0) {
                        solid[slab] += 1;
                    } else {
                        hist[(size_t)slab * num_bins + void_profile_bin_of(d, bin_width, num_bins)] += 1;
                    }
                }
            }
        }

        prof.hist         = hist.data();
        prof.solid        = solid.data();
        prof.total        = total.data();
        prof.num_slabs    = num_slabs;
        prof.num_bins     = num_bins;
        prof.bin_width    = bin_width;
        prof.z_min        = 0.0;
        prof.slab_height  = h * (double)nz / (double)num_slabs;
        prof.voxel_volume = h * h * h;
    }
};

}  // namespace

UTEST(viamd_void_profile, binning_conventions) {
    // The two mappings the field pass and every reader share. They are three lines each, which is
    // exactly why they are worth pinning: a one plane shift in the slab mapping or a half bin shift
    // in the distance mapping changes every reported profile a little and nothing enough to notice.

    // One slab per plane of voxels is the identity. Anything else is an off by one.
    for (int k = 0; k < 37; ++k) {
        EXPECT_EQ((uint32_t)k, void_profile_slab_of(k, 37, 37));
    }

    // Otherwise a plane belongs to the slab its centre falls in, and the slabs partition the grid:
    // non decreasing, starting at the first and ending at the last.
    const int dim_z = 100;
    const uint32_t ns = 7;
    uint32_t prev = 0;
    for (int k = 0; k < dim_z; ++k) {
        const uint32_t got  = void_profile_slab_of(k, dim_z, ns);
        const uint32_t want = (uint32_t)(((double)k + 0.5) / (double)dim_z * (double)ns);
        EXPECT_EQ(want, got);
        EXPECT_TRUE(got >= prev);
        prev = got;
    }
    EXPECT_EQ(0u, void_profile_slab_of(0, dim_z, ns));
    EXPECT_EQ(ns - 1, void_profile_slab_of(dim_z - 1, dim_z, ns));

    // Bin b is [b*w, (b+1)*w): closed below, open above.
    const double w = 0.25;
    EXPECT_EQ(0u, void_profile_bin_of(0.01, w, 16));
    EXPECT_EQ(0u, void_profile_bin_of(0.2499, w, 16));
    EXPECT_EQ(1u, void_profile_bin_of(0.25, w, 16));
    EXPECT_EQ(3u, void_profile_bin_of(0.99, w, 16));

    // The query clamps an unreached voxel to max_dist, which is the top edge of the last bin. It has
    // to land in that bin rather than one past the end.
    EXPECT_EQ(15u, void_profile_bin_of(16.0 * w, w, 16));
    EXPECT_EQ(15u, void_profile_bin_of(1000.0, w, 16));
}

UTEST(viamd_void_profile, slab_porosity_and_linear_accessible_volume) {
    // A solid slab spanning the middle of the box in z and everything in x and y. The distance to it
    // is |z - z0| - t/2, so the accessible volume falls exactly linearly with the probe radius:
    //
    //     V(R)/V = (H - t - 2R) / H,     porosity = (H - t) / H.
    //
    // Linear is the case which catches an off by half a bin: a reader that takes the bin containing
    // R whole, or reads at bin edges instead of splitting them, bends a straight line.
    const double h  = 0.25;
    const int    nz = 320;                   // H = 80
    const double H  = (double)nz * h;
    const double t  = 20.0;                  // Slab thickness
    const double z0 = 0.5 * H;

    Sampled s(4, 4, nz, h, 64, 256, 64.0, [&](double, double, double z) {
        return fabs(z - z0) - 0.5 * t;
    });

    const uint32_t ns = s.prof.num_slabs;

    EXPECT_NEAR((H - t) / H, void_profile_porosity(&s.prof, 0, ns), 1.0e-9);

    // Deliberately off the bin edges. A reader that counts the bin containing R whole agrees with a
    // correct one at every multiple of the bin width and nowhere else, so radii chosen on the grid
    // test nothing - which is exactly the blind spot this line exists to close.
    for (double R = 0.137; R <= 0.5 * (H - t) - 1.0; R += 2.313) {
        const double want = (H - t - 2.0 * R) / H;
        EXPECT_NEAR(want, void_profile_accessible_fraction(&s.prof, 0, ns, R), 2.0e-3);
    }

    // Beyond the widest gap nothing fits, and it must reach zero rather than stopping at a bin.
    EXPECT_NEAR(0.0, void_profile_accessible_fraction(&s.prof, 0, ns, 0.5 * (H - t) + 0.5), 1.0e-9);

    // -dV/dR is constant for a linear V, and equals 2/H of the box per unit R.
    const double n_box = (double)void_profile_num_total(&s.prof, 0, ns);
    EXPECT_NEAR(2.0 / H, void_profile_density(&s.prof, 0, ns, 10.0, h) / n_box, 5.0e-3);

    // Per slab, the answer is 0 or 1 and the interior is solid throughout.
    EXPECT_NEAR(1.0, void_profile_porosity(&s.prof, 0, 1), 1.0e-9);
    EXPECT_NEAR(0.0, void_profile_porosity(&s.prof, ns / 2, ns / 2 + 1), 1.0e-9);
}

UTEST(viamd_void_profile, hex_cylinder_array_accessible_volume) {
    // A hexagonal array of parallel cylinders along z, the idealization the CNF work is checked
    // against. The set where a probe of radius R fits is the complement of the cylinders dilated by
    // R, and while the dilated cylinders stay apart - Rc + R < a/2 - that has a closed form:
    //
    //     V(R)/V = 1 - (2 pi / sqrt 3) ((Rc + R) / a)^2.
    //
    // This is the case that separates an accessible volume from a pore size: the answer depends on
    // the dilated area alone and not at all on how the space between the cylinders is shaped.
    const double a  = 20.0;                  // Lattice constant
    const double Rc = 6.0;                   // Cylinder radius
    const double sq3 = sqrt(3.0);
    const double h  = 0.05;

    // Rectangular supercell of the hexagonal lattice: a x a*sqrt(3), holding exactly two lattice
    // points, at (0,0) and (a/2, a*sqrt(3)/2).
    const int nx = (int)llround(a / h);
    const int ny = (int)llround(a * sq3 / h);

    Sampled s(nx, ny, 4, h, 2, 512, 20.0, [&](double x, double y, double) {
        double best = 1.0e30;
        for (int i = -2; i <= 2; ++i) {
            for (int j = -2; j <= 4; ++j) {
                const double cx = (double)i * a + (double)j * 0.5 * a;
                const double cy = (double)j * 0.5 * sq3 * a;
                const double rho = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                if (rho - Rc < best) best = rho - Rc;
            }
        }
        return best;
    });

    const uint32_t ns = s.prof.num_slabs;

    for (double R = 0.0; R <= 3.5; R += 0.5) {
        const double r = Rc + R;
        const double want = 1.0 - (2.0 * PI_D / sq3) * (r / a) * (r / a);
        const double got  = void_profile_accessible_fraction(&s.prof, 0, ns, R);
        EXPECT_NEAR(want, got, 5.0e-3);
    }

    // The R = 0 endpoint of that curve is the porosity, by construction and not by coincidence.
    EXPECT_NEAR(void_profile_accessible_fraction(&s.prof, 0, ns, 0.0),
                void_profile_porosity(&s.prof, 0, ns), 1.0e-12);

    // Cylinders run along z, so every slab reports the same thing.
    EXPECT_NEAR(void_profile_porosity(&s.prof, 0, 1), void_profile_porosity(&s.prof, 1, 2), 1.0e-12);

    // The coarea derivative of the closed form is the dilated perimeter per unit area,
    // (4 pi / sqrt 3) (Rc + R) / a^2, which the histogram density has to reproduce.
    const double R = 2.0;
    const double n_box = (double)void_profile_num_total(&s.prof, 0, ns);
    const double want_dens = (4.0 * PI_D / sq3) * (Rc + R) / (a * a);
    EXPECT_NEAR(want_dens, void_profile_density(&s.prof, 0, ns, R, h) / n_box, 0.05 * want_dens);
}

UTEST(viamd_void_profile, film_extent_finds_the_half_density_surface) {
    // A film with smooth surfaces: the solid fraction follows a symmetric pair of tanh steps, which
    // is what a rough free surface averages to. The half density convention puts the boundary where
    // the profile has fallen to half its interior value, and for a symmetric step that is the
    // inflection of the step - so the answer is the step centre, to within a slab.
    const uint32_t ns = 200;
    const double slab_h = 1.0;               // 200 units of z
    const double z_lo = 50.0, z_hi = 150.0, w = 4.0, phi0 = 0.42;

    std::vector<uint64_t> hist((size_t)ns * 4, 0), solid(ns, 0), total(ns, 0);
    for (uint32_t s = 0; s < ns; ++s) {
        const double z = ((double)s + 0.5) * slab_h;
        const double f = phi0 * 0.5 * (tanh((z - z_lo) / w) - tanh((z - z_hi) / w));
        total[s] = 10000;
        solid[s] = (uint64_t)llround(f * 10000.0);
        hist[(size_t)s * 4 + 1] = total[s] - solid[s];
    }

    void_profile_t p = {};
    p.hist = hist.data(); p.solid = solid.data(); p.total = total.data();
    p.num_slabs = ns; p.num_bins = 4; p.bin_width = 1.0;
    p.z_min = 0.0; p.slab_height = slab_h; p.voxel_volume = 1.0;

    uint32_t beg = 0, end = 0;
    double interior = 0.0;
    ASSERT_TRUE(void_profile_film_extent(&p, 0.5, &beg, &end, &interior));

    EXPECT_NEAR(phi0, interior, 1.0e-3);
    EXPECT_NEAR(z_lo, void_profile_z_lo(&p, beg),     1.5 * slab_h);
    EXPECT_NEAR(z_hi, void_profile_z_hi(&p, end - 1), 1.5 * slab_h);

    // Raising the threshold can only shrink the film, never grow it.
    uint32_t beg2 = 0, end2 = 0;
    ASSERT_TRUE(void_profile_film_extent(&p, 0.9, &beg2, &end2, nullptr));
    EXPECT_TRUE(beg2 >= beg);
    EXPECT_TRUE(end2 <= end);

    // A profile with no solid has no film, rather than a film covering everything.
    std::vector<uint64_t> no_solid(ns, 0);
    p.solid = no_solid.data();
    EXPECT_FALSE(void_profile_film_extent(&p, 0.5, &beg, &end, nullptr));
}

UTEST(viamd_void_profile, film_extent_spans_an_internal_void) {
    // Two dense layers with a gap between them wide enough to drop the solid fraction to zero. That
    // is one film with a void inside it, not two films, so the extent has to reach from the first
    // crossing to the last rather than stopping at the end of the first run.
    const uint32_t ns = 120;
    std::vector<uint64_t> hist((size_t)ns * 2, 0), solid(ns, 0), total(ns, 0);
    for (uint32_t s = 0; s < ns; ++s) {
        const bool dense = (s >= 20 && s < 45) || (s >= 75 && s < 100);
        total[s] = 1000;
        solid[s] = dense ? 500 : 0;
        hist[(size_t)s * 2 + 1] = total[s] - solid[s];
    }

    void_profile_t p = {};
    p.hist = hist.data(); p.solid = solid.data(); p.total = total.data();
    p.num_slabs = ns; p.num_bins = 2; p.bin_width = 1.0;
    p.z_min = 0.0; p.slab_height = 1.0; p.voxel_volume = 1.0;

    uint32_t beg = 0, end = 0;
    ASSERT_TRUE(void_profile_film_extent(&p, 0.5, &beg, &end, nullptr));
    EXPECT_TRUE(beg >= 19 && beg <= 21);
    EXPECT_TRUE(end >= 99 && end <= 101);

    // The gap is pure void, and the film porosity has to sit between the dense value and 1.
    const double phi_film  = void_profile_porosity(&p, beg, end);
    const double phi_dense = void_profile_porosity(&p, 25, 40);
    EXPECT_TRUE(phi_film > phi_dense);
    EXPECT_TRUE(phi_film < 1.0);
}

UTEST(viamd_void_profile, reduction_identities) {
    // Properties the reductions have to satisfy whatever the geometry, checked on a field with no
    // closed form at all: a lattice of spheres of three different radii. These are the invariants
    // that a fast path or a cached cumulative would be at risk of breaking later.
    const double h = 0.4;
    const int n = 60;
    const double rad[3] = { 3.0, 5.0, 7.0 };
    const double cx[3]  = { 6.0, 14.0, 18.0 };
    const double cy[3]  = { 5.0, 17.0, 8.0 };
    const double cz[3]  = { 7.0, 9.0, 19.0 };

    Sampled s(n, n, n, h, 12, 256, 16.0, [&](double x, double y, double z) {
        double best = 1.0e30;
        for (int c = 0; c < 3; ++c) {
            const double d = sqrt((x - cx[c]) * (x - cx[c]) + (y - cy[c]) * (y - cy[c]) + (z - cz[c]) * (z - cz[c])) - rad[c];
            if (d < best) best = d;
        }
        return best;
    });

    const uint32_t ns = s.prof.num_slabs;
    const uint64_t n_tot   = void_profile_num_total(&s.prof, 0, ns);
    const uint64_t n_solid = void_profile_num_solid(&s.prof, 0, ns);

    // Every voxel is solid or binned, so R = 0 recovers the complement of the solid exactly. This is
    // what makes porosity and accessible volume the same quantity rather than two that agree.
    EXPECT_NEAR((double)(n_tot - n_solid), void_profile_count_above(&s.prof, 0, ns, 0.0), 1.0e-6);
    EXPECT_NEAR(1.0 - (double)n_solid / (double)n_tot, void_profile_porosity(&s.prof, 0, ns), 1.0e-12);

    // Slabs partition the box: the whole is the sum of its parts, for counts and for volumes.
    double sum_count = 0.0, sum_vol = 0.0;
    uint64_t sum_tot = 0;
    for (uint32_t sl = 0; sl < ns; ++sl) {
        sum_count += void_profile_count_above(&s.prof, sl, sl + 1, 1.5);
        sum_vol   += void_profile_accessible_volume(&s.prof, sl, sl + 1, 1.5);
        sum_tot   += void_profile_num_total(&s.prof, sl, sl + 1);
    }
    EXPECT_NEAR(void_profile_count_above(&s.prof, 0, ns, 1.5), sum_count, 1.0e-6);
    EXPECT_NEAR(void_profile_accessible_volume(&s.prof, 0, ns, 1.5), sum_vol, 1.0e-6);
    EXPECT_EQ(n_tot, sum_tot);

    // Monotone in R, and never increasing.
    double prev = void_profile_count_above(&s.prof, 0, ns, 0.0);
    for (double R = 0.01; R < 16.0; R += 0.01) {
        const double cur = void_profile_count_above(&s.prof, 0, ns, R);
        EXPECT_TRUE(cur <= prev + 1.0e-9);
        prev = cur;
    }
    EXPECT_NEAR(0.0, void_profile_count_above(&s.prof, 0, ns, 16.0), 1.0e-9);

    // Continuous across a bin, not a staircase: the value at a bin's midpoint has to sit between
    // its two edges, near the mean of them. A reader that counts the bin containing R whole returns
    // the lower edge's value right across the bin and fails this while still passing every test
    // whose radii happen to land on the grid.
    const double bw = s.prof.bin_width;
    for (int b = 2; b < 24; ++b) {
        const double lo  = void_profile_count_above(&s.prof, 0, ns, (double)b * bw);
        const double hi  = void_profile_count_above(&s.prof, 0, ns, (double)(b + 1) * bw);
        const double mid = void_profile_count_above(&s.prof, 0, ns, ((double)b + 0.5) * bw);
        EXPECT_NEAR(0.5 * (lo + hi), mid, 0.02 * (lo - hi) + 1.0e-9);
    }

    // The density is the derivative of that curve, so integrating it back over the whole range has
    // to return every void voxel.
    const double w = s.prof.bin_width;
    double integral = 0.0;
    for (int b = 0; b < (int)s.prof.num_bins; ++b) {
        integral += void_profile_density(&s.prof, 0, ns, ((double)b + 0.5) * w, 0.0) * w;
    }
    EXPECT_NEAR((double)(n_tot - n_solid), integral, 1.0e-3 * (double)(n_tot - n_solid));

    // An empty or inverted slab range is a question with no answer, not a crash or a division.
    EXPECT_NEAR(0.0, void_profile_count_above(&s.prof, 5, 5, 1.0), 1.0e-12);
    EXPECT_NEAR(0.0, void_profile_porosity(&s.prof, 7, 3), 1.0e-12);
    EXPECT_NEAR(0.0, void_profile_accessible_fraction(&s.prof, ns + 10, ns + 20, 1.0), 1.0e-12);
}
