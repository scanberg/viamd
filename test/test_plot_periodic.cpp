#include "utest.h"

#include <plot_periodic_core.h>

#include <math.h>
#include <vector>

static bool is_nan(double v) { return v != v; }

UTEST(plot_periodic, wrap_and_near) {
    EXPECT_NEAR(-170.0, periodic_wrap(190, -180, 360), 1e-9);
    EXPECT_NEAR(-180.0, periodic_wrap(-180, -180, 360), 1e-9);
    // The end of the period is its beginning
    EXPECT_NEAR(-180.0, periodic_wrap(180, -180, 360), 1e-9);
    EXPECT_NEAR(10.0, periodic_wrap(370 + 720, -180, 360), 1e-9);
    // As near as it comes: -190 is 20 from -170, where 170 is 340
    EXPECT_NEAR(-190.0, periodic_near(170, -170, 360), 1e-9);
}

// 170 to -170 is 20 the short way, across 180: out at the top of [-180, 180[ and in at the bottom, where it crosses
UTEST(plot_periodic, polyline_cut_where_it_crosses) {
    const float x[] = {0, 1, 2};
    const float y[] = {170, -170, -160};
    double px[16], py[16];
    const int m = periodic_polyline(px, py, 16, x, y, 3, 1, 0, 1.0, -180, 360, 360 * PERIODIC_MAX_STEP);
    ASSERT_EQ(6, m);
    EXPECT_EQ(0.0, px[0]); EXPECT_EQ(170.0, py[0]);
    EXPECT_NEAR(0.5, px[1], 1e-9); EXPECT_NEAR(180.0, py[1], 1e-9);
    EXPECT_TRUE(is_nan(py[2]));
    EXPECT_NEAR(0.5, px[3], 1e-9); EXPECT_NEAR(-180.0, py[3], 1e-9);
    EXPECT_EQ(-170.0, py[4]);
    EXPECT_EQ(-160.0, py[5]);
}

// The same values cut at 0 instead: [0, 360[, through which they run unbroken
UTEST(plot_periodic, polyline_unbroken_away_from_the_cut) {
    const float x[] = {0, 1, 2};
    const float y[] = {170, -170, -160};
    double px[16], py[16];
    const int m = periodic_polyline(px, py, 16, x, y, 3, 1, 0, 1.0, 0, 360, 360 * PERIODIC_MAX_STEP);
    ASSERT_EQ(3, m);
    EXPECT_EQ(170.0, py[0]);
    EXPECT_EQ(190.0, py[1]);
    EXPECT_EQ(200.0, py[2]);
}

// Downwards across the cut; a step too long to tell which way it went, and a frame not evaluated, are gaps
UTEST(plot_periodic, polyline_gaps) {
    const float x[] = {0, 1, 2, 3, 4};
    const float y[] = {-175, 175, 0, NAN, 10};
    double px[32], py[32];
    const int m = periodic_polyline(px, py, 32, x, y, 5, 1, 0, 1.0, -180, 360, 360 * PERIODIC_MAX_STEP);
    // (0,-175) (0.5,-180) gap (0.5,180) (1,175) gap (2,0) gap (4,10)
    ASSERT_EQ(9, m);
    EXPECT_NEAR(-180.0, py[1], 1e-9);
    EXPECT_TRUE(is_nan(py[2]));
    EXPECT_NEAR(180.0, py[3], 1e-9);
    EXPECT_TRUE(is_nan(py[5]));
    EXPECT_EQ(0.0, py[6]);
    EXPECT_TRUE(is_nan(py[7]));
    EXPECT_EQ(10.0, py[8]);
}

// A member of a population (stride, k), in radians shown in degrees
UTEST(plot_periodic, polyline_member_and_scale) {
    const float x[] = {0, 1};
    const float y[] = {0, 3.0f, 0, -3.0f};
    double px[16], py[16];
    const int m = periodic_polyline(px, py, 16, x, y, 2, 2, 1, 180.0 / M_PI, -180, 360, 360 * PERIODIC_MAX_STEP);
    ASSERT_EQ(5, m);
    EXPECT_NEAR(3.0 * 180.0 / M_PI, py[0], 1e-3);
    EXPECT_NEAR(180.0, py[1], 1e-9);
}

// A band whose centre crosses the end of the window stays one run, unwrapped past it
UTEST(plot_periodic, band_runs_unwrapped_across_the_cut) {
    const double x[]  = {0, 1, 2, 3};
    const double lo[] = {160, 170, -180, -170};
    const double hi[] = {180, 190, -160, -150};
    double ox[8], olo[8], ohi[8];
    int rb[4], rc[4];
    const int runs = periodic_band_runs(ox, olo, ohi, rb, rc, x, lo, hi, 4, -180, 360, 360 * PERIODIC_MAX_STEP);
    ASSERT_EQ(1, runs);
    ASSERT_EQ(4, rc[0]);
    EXPECT_NEAR(200.0, 0.5 * (olo[3] + ohi[3]), 1e-9);
    EXPECT_NEAR(20.0, ohi[3] - olo[3], 1e-9);
}

// A centre which goes around and around (a rotor) is taken back towards the window, without a break in the band
UTEST(plot_periodic, band_runs_of_a_rotor) {
    const int n = 2000;
    std::vector<double> x(n), lo(n), hi(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i;
        double c = fmod(i * 10.0, 360.0);
        if (c >= 180) c -= 360;
        lo[i] = c - 5;
        hi[i] = c + 5;
    }
    std::vector<double> ox(2 * n), olo(2 * n), ohi(2 * n);
    std::vector<int> rb(n), rc(n);
    const int runs = periodic_band_runs(ox.data(), olo.data(), ohi.data(), rb.data(), rc.data(), x.data(), lo.data(), hi.data(), n, -180, 360, 360 * PERIODIC_MAX_STEP);
    ASSERT_GT(runs, 1);
    int total = 0;
    for (int r = 0; r < runs; ++r) {
        for (int i = 0; i < rc[r]; ++i) {
            const double c = 0.5 * (olo[rb[r] + i] + ohi[rb[r] + i]);
            EXPECT_GE(c, -180.0 - 180.0 - 1e-9);
            EXPECT_LT(c, -180.0 + 1.5 * 360.0 + 10.0);
            total += 1;
        }
        // Each run takes over from the sample before, a period away
        if (r > 0) EXPECT_EQ(ox[rb[r - 1] + rc[r - 1] - 1], ox[rb[r]]);
    }
    EXPECT_GE(total, n);
    EXPECT_LE(total, 2 * n);
}

// The least populated place is found around the circle, and not fooled by one empty bin among full ones
UTEST(plot_periodic, least_populated) {
    float d[36];
    for (int i = 0; i < 36; ++i) d[i] = 1.0f;
    // A dip across the end of the period, and a single empty bin elsewhere
    d[35] = 0.0f; d[0] = 0.0f; d[1] = 0.0f;
    d[18] = 0.0f;
    const int b = periodic_least_populated(d, 36);
    EXPECT_EQ(0, b);
    EXPECT_EQ(-1, periodic_least_populated(d, 0));
}
