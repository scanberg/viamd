#include "utest.h"

#include <correlation/correlation_core.h>

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdint.h>
#include <vector>

using namespace correlation;

UTEST(correlation, rank_with_ties) {
    const double v[] = {10, 30, 20, 30, 5};
    double r[5];
    corr_rank(r, v, 5);
    EXPECT_EQ(2.0, r[0]);
    EXPECT_EQ(4.5, r[1]);
    EXPECT_EQ(3.0, r[2]);
    EXPECT_EQ(4.5, r[3]);
    EXPECT_EQ(1.0, r[4]);
}

UTEST(correlation, pearson) {
    const double x[] = {1, 2, 3, 4, 5};
    const double up[] = {3, 5, 7, 9, 11};
    const double down[] = {5, 4, 3, 2, 1};
    const double y[] = {2, 1, 4, 3, 5};
    double r = 0;
    ASSERT_TRUE(corr_pearson(&r, x, up, 5));
    EXPECT_NEAR(1.0, r, 1e-12);
    ASSERT_TRUE(corr_pearson(&r, x, down, 5));
    EXPECT_NEAR(-1.0, r, 1e-12);
    // Sxy = 8, Sxx = Syy = 10
    ASSERT_TRUE(corr_pearson(&r, x, y, 5));
    EXPECT_NEAR(0.8, r, 1e-12);
}

UTEST(correlation, pearson_undefined) {
    const double x[] = {1, 2, 3, 4};
    const double flat[] = {2, 2, 2, 2};
    double r = 42;
    EXPECT_FALSE(corr_pearson(&r, x, flat, 4));
    EXPECT_FALSE(corr_pearson(&r, x, x, 2));
    EXPECT_EQ(42.0, r);
}

// Monotonic but not linear: Spearman's sees it whole, Pearson's does not
UTEST(correlation, spearman_monotonic) {
    std::vector<double> x, y;
    for (int i = 0; i < 20; ++i) {
        x.push_back(i);
        y.push_back(exp(0.5 * i));
    }
    double rho = 0, r = 0;
    ASSERT_TRUE(corr_spearman(&rho, x.data(), y.data(), x.size()));
    EXPECT_NEAR(1.0, rho, 1e-12);
    ASSERT_TRUE(corr_pearson(&r, x.data(), y.data(), x.size()));
    EXPECT_LT(r, 0.9);
}

// Of ranks without ties: 1 - 6 sum d^2 / (n (n^2 - 1))
UTEST(correlation, spearman_closed_form) {
    const double x[] = {1, 2, 3, 4, 5};
    const double y[] = {2, 1, 4, 3, 5};
    double rho = 0;
    ASSERT_TRUE(corr_spearman(&rho, x, y, 5));
    EXPECT_NEAR(1.0 - 6.0 * 4.0 / (5.0 * 24.0), rho, 1e-12);
}

// Angles that turn together: 1, however they wrap around the circle; against each other: -1
UTEST(correlation, circular_together_across_the_cut) {
    std::vector<double> a, b, c;
    for (int i = 0; i < 12; ++i) {
        const double t = -0.6 + 0.1 * i;
        a.push_back(remainder(M_PI + t, 2 * M_PI));          // About pi: across the cut of [-pi, pi[
        b.push_back(remainder(M_PI + t + 0.3, 2 * M_PI));
        c.push_back(remainder(-t, 2 * M_PI));
    }
    double r = 0;
    ASSERT_TRUE(corr_circular(&r, a.data(), b.data(), a.size()));
    EXPECT_NEAR(1.0, r, 1e-9);
    ASSERT_TRUE(corr_circular(&r, a.data(), c.data(), a.size()));
    EXPECT_NEAR(-1.0, r, 1e-9);
}

// A rotation of either set of angles does not change it
UTEST(correlation, circular_rotation_invariant) {
    const double a[] = {0.1, 0.9, 2.0, 2.5, 4.0, 5.5, 6.0};
    const double b[] = {0.4, 0.2, 1.9, 3.0, 3.1, 5.0, 0.3};
    double b2[7];
    for (int i = 0; i < 7; ++i) b2[i] = b[i] + 1.234;
    double r1 = 0, r2 = 0;
    ASSERT_TRUE(corr_circular(&r1, a, b, 7));
    ASSERT_TRUE(corr_circular(&r2, a, b2, 7));
    EXPECT_NEAR(r1, r2, 1e-12);
    EXPECT_LE(fabs(r1), 1.0);
}

// x a function of the cosine and sine of the angle: 1, whatever the phase; unsigned, so -x as well
UTEST(correlation, circular_linear_exact) {
    std::vector<double> a, x, nx;
    for (int i = 0; i < 36; ++i) {
        const double t = i * (2 * M_PI / 36);
        a.push_back(t);
        x.push_back(3.0 * cos(t - 1.0) + 2.0);
        nx.push_back(-x.back());
    }
    double r = 0;
    ASSERT_TRUE(corr_circular_linear(&r, a.data(), x.data(), a.size()));
    EXPECT_NEAR(1.0, r, 1e-9);
    ASSERT_TRUE(corr_circular_linear(&r, a.data(), nx.data(), a.size()));
    EXPECT_NEAR(1.0, r, 1e-9);
}

// x of the angle, but not through its cosine or sine: none
UTEST(correlation, circular_linear_unrelated) {
    std::vector<double> a, x;
    for (int i = 0; i < 64; ++i) {
        a.push_back(i * (2 * M_PI / 64));
        x.push_back(cos(2.0 * a.back()));   // Twice around for once around the angle
    }
    double r = 1;
    ASSERT_TRUE(corr_circular_linear(&r, a.data(), x.data(), a.size()));
    EXPECT_NEAR(0.0, r, 1e-9);
}

UTEST(correlation, corr_of_dispatches) {
    const double x[] = {0.1, 0.5, 1.0, 1.4, 2.2, 3.0};
    const double y[] = {1.0, 1.3, 2.1, 2.0, 3.5, 3.9};
    double a = 0, b = 0;
    ASSERT_TRUE(corr_of(&a, CorrKind_Linear, x, y, 6));
    ASSERT_TRUE(corr_pearson(&b, x, y, 6));
    EXPECT_EQ(b, a);
    ASSERT_TRUE(corr_of(&a, CorrKind_CircularX, x, y, 6));
    ASSERT_TRUE(corr_circular_linear(&b, x, y, 6));
    EXPECT_EQ(b, a);
    ASSERT_TRUE(corr_of(&a, CorrKind_CircularY, y, x, 6));
    EXPECT_EQ(b, a);
    ASSERT_TRUE(corr_of(&a, CorrKind_Circular, x, y, 6));
    ASSERT_TRUE(corr_circular(&b, x, y, 6));
    EXPECT_EQ(b, a);
    EXPECT_EQ(CorrKind_CircularY, corr_kind(false, true));
}

UTEST(correlation, lag_match) {
    // Frames with gaps: what is lag away is found, what falls in a gap is not
    const uint32_t f[] = {0, 1, 2, 4, 5, 8};
    uint32_t ia[6], ib[6];
    size_t n = corr_lag_match(ia, ib, f, 6, 0);
    EXPECT_EQ((size_t)6, n);
    n = corr_lag_match(ia, ib, f, 6, 1);
    ASSERT_EQ((size_t)3, n);    // 0-1, 1-2, 4-5
    EXPECT_EQ(0u, ia[0]); EXPECT_EQ(1u, ib[0]);
    EXPECT_EQ(1u, ia[1]); EXPECT_EQ(2u, ib[1]);
    EXPECT_EQ(3u, ia[2]); EXPECT_EQ(4u, ib[2]);
    n = corr_lag_match(ia, ib, f, 6, -3);
    ASSERT_EQ((size_t)3, n);    // 4-1, 5-2, 8-5
    EXPECT_EQ(3u, ia[0]); EXPECT_EQ(1u, ib[0]);
    EXPECT_EQ(4u, ia[1]); EXPECT_EQ(2u, ib[1]);
    EXPECT_EQ(5u, ia[2]); EXPECT_EQ(4u, ib[2]);
    EXPECT_EQ((size_t)0, corr_lag_match(ia, ib, f, 6, 100));
}

// y following x by 3 frames: at lag 3 they are the same
UTEST(correlation, lagged_series) {
    const int n = 200;
    std::vector<double> x(n), y(n);
    std::vector<uint32_t> f(n);
    uint32_t seed = 12345;
    for (int i = 0; i < n; ++i) {
        seed = seed * 1664525u + 1013904223u;
        x[i] = (double)(seed >> 8) / (double)(1u << 24);
        f[i] = (uint32_t)i;
    }
    for (int i = 0; i < n; ++i) y[i] = i >= 3 ? x[i - 3] : 0.5;
    std::vector<uint32_t> ia(n), ib(n);
    std::vector<double> a(n), b(n);
    auto at = [&](int lag) {
        const size_t m = corr_lag_match(ia.data(), ib.data(), f.data(), n, lag);
        for (size_t k = 0; k < m; ++k) { a[k] = x[ia[k]]; b[k] = y[ib[k]]; }
        double r = 0;
        corr_pearson(&r, a.data(), b.data(), m);
        return r;
    };
    EXPECT_NEAR(1.0, at(3), 1e-12);
    EXPECT_LT(fabs(at(0)), 0.3);
    EXPECT_LT(fabs(at(-3)), 0.3);
}

UTEST(correlation, pairing_shapes) {
    EXPECT_EQ(CorrShape_Empty, corr_shape(0, 3));
    EXPECT_EQ(CorrShape_Empty, corr_shape(2, 0));
    EXPECT_EQ(CorrShape_Mismatch, corr_shape(2, 3));
    EXPECT_EQ(CorrShape_Valid, corr_shape(1, 4));
    EXPECT_EQ(CorrShape_Valid, corr_shape(3, 3));
    EXPECT_EQ(0, corr_num_pairs(2, 3));
    EXPECT_EQ(4, corr_num_pairs(1, 4));
    EXPECT_EQ(4, corr_num_pairs(4, 1));
    EXPECT_EQ(3, corr_num_pairs(3, 3));
    // 1 x n: the one X in every pair, the p'th Y
    EXPECT_EQ(0, corr_pair_slot(1, 2));
    EXPECT_EQ(2, corr_pair_slot(4, 2));
}

UTEST(correlation, pairing_members) {
    bool truncated = true;
    EXPECT_EQ(1, corr_members(1, 1, &truncated));
    EXPECT_FALSE(truncated);
    // A single value against a population: against each member
    EXPECT_EQ(5, corr_members(1, 5, &truncated));
    EXPECT_FALSE(truncated);
    EXPECT_EQ(5, corr_members(5, 1, &truncated));
    // Two populations: member by member, up to the fewer
    EXPECT_EQ(3, corr_members(3, 7, &truncated));
    EXPECT_TRUE(truncated);
    EXPECT_EQ(4, corr_members(4, 4, &truncated));
    EXPECT_FALSE(truncated);
    EXPECT_EQ(0, corr_member_slot(1, 3));
    EXPECT_EQ(3, corr_member_slot(5, 3));
}

UTEST(correlation, interp_linear) {
    // Two members per sample
    const float v[] = {0, 10, 2, 20, 4, 30};
    EXPECT_NEAR(1.0, corr_interp(v, 2, 3, 1.0, 0.5, 0, 0.0), 1e-9);
    EXPECT_NEAR(25.0, corr_interp(v, 2, 3, 1.0, 1.5, 1, 0.0), 1e-9);
    EXPECT_NEAR(4.0, corr_interp(v, 2, 3, 2.0, 1.0, 0, 0.0), 1e-9);
    // Clamped to the ends
    EXPECT_NEAR(0.0, corr_interp(v, 2, 3, 1.0, -1.0, 0, 0.0), 1e-9);
    EXPECT_NEAR(30.0, corr_interp(v, 2, 3, 1.0, 9.0, 1, 0.0), 1e-9);
}

UTEST(correlation, interp_periodic) {
    // 170 to -170 degrees: the short way, through 180
    const float v[] = {170, -170};
    EXPECT_NEAR(180.0, corr_interp(v, 1, 2, 1.0, 0.5, 0, 360.0), 1e-9);
    EXPECT_NEAR(175.0, corr_interp(v, 1, 2, 1.0, 0.25, 0, 360.0), 1e-9);
    // Not periodic: through 0
    EXPECT_NEAR(0.0, corr_interp(v, 1, 2, 1.0, 0.5, 0, 0.0), 1e-9);
}

UTEST(correlation, lag1_autocorrelation) {
    // Alternating: -1-ish; a slow ramp: close to 1
    std::vector<double> alt(100), ramp(100);
    std::vector<uint32_t> g(100, 0);
    for (int i = 0; i < 100; ++i) {
        alt[i] = (i & 1) ? 1.0 : -1.0;
        ramp[i] = (double)i;
    }
    EXPECT_LT(corr_lag1(alt.data(), g.data(), 100, 1), -0.95);
    EXPECT_GT(corr_lag1(ramp.data(), g.data(), 100, 1), 0.95);
    // Two ramps interleaved, each its own sequence: still alike from one to the next
    std::vector<double> two(100);
    for (int i = 0; i < 100; ++i) {
        g[i] = (uint32_t)(i & 1);
        two[i] = (i & 1) ? 1000.0 - i : (double)i;
    }
    EXPECT_GT(corr_lag1(two.data(), g.data(), 100, 2), 0.9);
    // Taken as one sequence, they alternate
    std::vector<uint32_t> g0(100, 0);
    EXPECT_LT(corr_lag1(two.data(), g0.data(), 100, 1), 0.0);
}

UTEST(correlation, effective_n) {
    EXPECT_NEAR(1000.0, corr_effective_n(1000.0, 0.0, 0.9), 1e-9);
    // (1 - 0.81) / (1 + 0.81)
    EXPECT_NEAR(1000.0 * 0.19 / 1.81, corr_effective_n(1000.0, 0.9, 0.9), 1e-6);
    // Never more than there are, nor fewer than 2
    EXPECT_NEAR(100.0, corr_effective_n(100.0, -0.5, 0.5), 1e-9);
    EXPECT_NEAR(2.0, corr_effective_n(100.0, 1.0, 1.0), 1e-9);
}

UTEST(correlation, fisher_interval) {
    double lo = 0, hi = 0;
    ASSERT_TRUE(corr_fisher_interval(&lo, &hi, 0.5, 103.0, 1.959964));
    // atanh(0.5) = 0.549306, se = 0.1
    EXPECT_NEAR(tanh(0.549306 - 0.1959964), lo, 1e-5);
    EXPECT_NEAR(tanh(0.549306 + 0.1959964), hi, 1e-5);
    EXPECT_LT(lo, 0.5);
    EXPECT_GT(hi, 0.5);
    // Few samples: wide, around 0
    ASSERT_TRUE(corr_fisher_interval(&lo, &hi, 0.3, 8.0, 1.959964));
    EXPECT_LT(lo, 0.0);
    EXPECT_FALSE(corr_fisher_interval(&lo, &hi, 0.3, 3.0, 1.959964));
    // Perfect: finite
    ASSERT_TRUE(corr_fisher_interval(&lo, &hi, 1.0, 50.0, 1.959964));
    EXPECT_TRUE(isfinite(lo) && isfinite(hi));
}

UTEST(correlation, circular_linear_parts_signed) {
    // x = -cos(a): against the cosine -1, Mardia's 1 without a sign
    const int n = 72;
    std::vector<double> a(n), x(n);
    for (int i = 0; i < n; ++i) {
        a[i] = 2.0 * M_PI * i / n;
        x[i] = -cos(a[i]);
    }
    double rc = 0, rs = 0, R = 0;
    ASSERT_TRUE(corr_circular_linear_parts(&rc, &rs, a.data(), x.data(), n));
    EXPECT_NEAR(-1.0, rc, 1e-9);
    EXPECT_NEAR(0.0, rs, 1e-9);
    ASSERT_TRUE(corr_circular_linear(&R, a.data(), x.data(), n));
    EXPECT_NEAR(1.0, R, 1e-9);
}
