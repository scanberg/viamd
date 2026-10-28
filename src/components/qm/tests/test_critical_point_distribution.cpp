#include "utest.h"
#include <qm/critical_point_distribution_core.h>

#include <algorithm>
#include <cmath>
#include <limits>

using critical_point_distribution::estimate;

UTEST(critical_point_distribution, focus_radius_is_three_angstrom_in_3d) {
    using namespace critical_point_distribution;
    const Position centers[] = {{0.0, 0.0, 0.0}};
    EXPECT_TRUE(within_focus({0.0, 0.0, 0.0}, centers));
    EXPECT_TRUE(within_focus({3.0, 0.0, 0.0}, centers));
    EXPECT_TRUE(within_focus({1.0, 2.0, 2.0}, centers));
    EXPECT_FALSE(within_focus({3.0001, 0.0, 0.0}, centers));
    EXPECT_FALSE(within_focus({0.0, 0.0, 3.0001}, centers));
    EXPECT_FALSE(within_focus({2.0, 2.0, 2.0}, centers));
}

UTEST(critical_point_distribution, grouped_focus_is_union_of_neighborhoods) {
    using namespace critical_point_distribution;
    const Position centers[] = {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}};
    EXPECT_TRUE(within_focus({-3.0, 0.0, 0.0}, centers));
    EXPECT_TRUE(within_focus({13.0, 0.0, 0.0}, centers));
    EXPECT_FALSE(within_focus({5.0, 0.0, 0.0}, centers));
    EXPECT_FALSE(within_focus({0.0, 0.0, 0.0}, {}));
}

UTEST(critical_point_distribution, density_filter_is_inclusive_and_optional) {
    critical_point_distribution::Filter filter;
    EXPECT_TRUE(filter.accepts(1.0e-9));
    EXPECT_TRUE(filter.accepts(100.0));
    filter.enabled = true;
    filter.ranges = {{0.001, 0.1}};
    EXPECT_FALSE(filter.accepts(0.0009));
    EXPECT_TRUE(filter.accepts(0.001));
    EXPECT_TRUE(filter.accepts(0.05));
    EXPECT_TRUE(filter.accepts(0.1));
    EXPECT_FALSE(filter.accepts(0.1001));
    filter.ranges = {{0.1, 0.1}};
    EXPECT_TRUE(filter.accepts(0.1));
    EXPECT_FALSE(filter.accepts(0.05));
    filter.enabled = false;
    EXPECT_TRUE(filter.accepts(0.05));
    EXPECT_FALSE(filter.accepts(std::numeric_limits<double>::quiet_NaN()));
}

UTEST(critical_point_distribution, category_filters_are_independent) {
    critical_point_distribution::Filter filters[4];
    filters[1] = {true, {{0.001, 0.01}}};
    EXPECT_FALSE(filters[1].accepts(0.1));
    EXPECT_TRUE(filters[0].accepts(0.1));
    EXPECT_TRUE(filters[2].accepts(0.1));
    EXPECT_TRUE(filters[3].accepts(0.1));
    filters[2] = {true, {{0.05, 0.2}}};
    EXPECT_TRUE(filters[2].accepts(0.1));
    EXPECT_FALSE(filters[1].accepts(0.1));
    filters[1].enabled = false;
    EXPECT_TRUE(filters[1].accepts(0.1));
    EXPECT_FALSE(filters[2].accepts(0.01));
}

UTEST(critical_point_distribution, density_filter_accepts_union_of_ranges) {
    critical_point_distribution::Filter filter;
    filter.enabled = true;
    EXPECT_FALSE(filter.accepts(0.1));
    filter.ranges = {{0.001, 0.01}, {0.1, 0.2}};
    EXPECT_TRUE(filter.accepts(0.005));
    EXPECT_TRUE(filter.accepts(0.15));
    EXPECT_FALSE(filter.accepts(0.05));
    EXPECT_FALSE(filter.accepts(0.3));
    filter.ranges.push_back({0.04, 0.06});
    EXPECT_TRUE(filter.accepts(0.05));
}

UTEST(critical_point_distribution, components_sum_to_total) {
    using critical_point_distribution::estimate;
    const double values[] = {0.3, 0.001, 0.25, 0.002, 0.28, 0.0015};
    const uint32_t component[] = {0, 1, 0, 1, 0, 1};
    const auto split = estimate(values, component, 2);
    const auto whole = estimate(values);
    ASSERT_TRUE(split.valid);
    ASSERT_EQ(split.components.size(), size_t(2));
    EXPECT_TRUE(whole.components.empty());
    double low = 0.0, high = 0.0;
    for (int j = 0; j < split.sample_count; ++j) {
        EXPECT_NEAR(split.probability[j], whole.probability[j], 1.0e-12);
        EXPECT_NEAR(split.components[0][j] + split.components[1][j], split.probability[j], 1.0e-12);
        (split.density[j] < 0.02 ? low : high) += split.components[0][j];
    }
    EXPECT_LT(low, high * 1.0e-3);
    const uint32_t bad[] = {0, 1, 0, 1, 0, 2};
    EXPECT_FALSE(estimate(values, bad, 2).valid);
}

UTEST(critical_point_distribution, empty) {
    EXPECT_FALSE(estimate({}).valid);
}

UTEST(critical_point_distribution, invalid_density_is_reported) {
    const double values[] = {std::numeric_limits<double>::quiet_NaN()};
    EXPECT_FALSE(estimate(values).valid);
}

static void check_curve(int* utest_result, std::span<const double> values, bool logarithmic) {
    const auto curve = estimate(values);
    ASSERT_TRUE(curve.valid);
    EXPECT_EQ(curve.logarithmic, logarithmic);
    double area = 0.0;
    for (int i = 0; i < curve.sample_count; ++i) {
        EXPECT_TRUE(std::isfinite(curve.density[i]));
        EXPECT_TRUE(std::isfinite(curve.probability[i]));
        EXPECT_GE(curve.probability[i], 0.0);
        if (i > 0) {
            EXPECT_GT(curve.density[i], curve.density[i - 1]);
            const double dx = logarithmic ? std::log10(curve.density[i] / curve.density[i - 1])
                                          : curve.density[i] - curve.density[i - 1];
            area += dx * 0.5 * (curve.probability[i] + curve.probability[i - 1]);
        }
    }
    EXPECT_NEAR(area, 1.0, 0.002);
}

UTEST(critical_point_distribution, single_and_coincident) {
    const double single[] = {0.003};
    const double coincident[] = {0.003, 0.003, 0.003, 0.003};
    check_curve(utest_result, single, true);
    check_curve(utest_result, coincident, true);
    const auto a = estimate(single), b = estimate(coincident);
    for (int i = 0; i < a.sample_count; ++i) {
        EXPECT_NEAR(a.probability[i], b.probability[i], 1.0e-12);
    }
}

UTEST(critical_point_distribution, wide_dynamic_range) {
    const double values[] = {100.0, 0.0001, 0.00011, 0.0001, 0.02, 0.021, 0.023, 5.0};
    check_curve(utest_result, values, true);
}

UTEST(critical_point_distribution, cluster_is_narrower_than_isolated_point) {
    const double values[] = {0.001, 0.00101, 0.00102, 0.00103, 0.00104, 0.00105, 0.00106, 0.00107, 1.0};
    const auto curve = estimate(values);
    double cluster_peak = 0.0, isolated_peak = 0.0;
    for (int i = 0; i < curve.sample_count; ++i) {
        if (curve.density[i] > 0.0005 && curve.density[i] < 0.002) {
            cluster_peak = std::max(cluster_peak, curve.probability[i]);
        }
        if (curve.density[i] > 0.5 && curve.density[i] < 2.0) {
            isolated_peak = std::max(isolated_peak, curve.probability[i]);
        }
    }
    EXPECT_GT(cluster_peak, 10.0 * isolated_peak);
    check_curve(utest_result, values, true);
}

UTEST(critical_point_distribution, nonpositive_values_use_linear_axis) {
    const double values[] = {-0.001, 0.0, 0.002};
    const double zero[] = {0.0};
    check_curve(utest_result, values, false);
    check_curve(utest_result, zero, false);
}

UTEST(critical_point_distribution, log_distribution_is_scale_invariant) {
    const double values[] = {0.001, 0.002, 0.003, 0.1, 2.0};
    const double scaled[] = {1.0, 2.0, 3.0, 100.0, 2000.0};
    const auto a = estimate(values), b = estimate(scaled);
    for (int i = 0; i < a.sample_count; ++i) {
        EXPECT_NEAR(b.density[i] / a.density[i], 1000.0, 1.0e-8);
        EXPECT_NEAR(a.probability[i], b.probability[i], 1.0e-10);
    }
}
