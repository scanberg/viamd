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
