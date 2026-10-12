#include "utest.h"
#include <components/analysis_density.h>
#include <limits>

using analysis_density::Triangle;

UTEST(analysis_density, normalized_independent_of_frame_and_group_count) {
    const vec2_t points[] = {{0,0}, {1,0}, {0.5f,Triangle::height},
                            {0,0}, {1,0}, {0.5f,Triangle::height}};
    Triangle density;
    density.update(points, 3, 1, 0, 3, 2.0f);
    const auto original = density.values;
    double sum = 0;
    for (float v : original) {
        EXPECT_TRUE(std::isfinite(v));
        EXPECT_GE(v, 0.0f);
        sum += v;
    }
    EXPECT_NEAR(1.0, sum / original.size(), 1e-5);
    density.update(points, 3, 2, 0, 3, 2.0f);
    for (size_t i = 0; i < original.size(); ++i) {
        EXPECT_NEAR(original[i], density.values[i], 1e-4f);
    }
    density.update(points, 6, 1, 0, 6, 2.0f);
    for (size_t i = 0; i < original.size(); ++i) {
        EXPECT_NEAR(original[i], density.values[i], 1e-4f);
    }
    density.update(points, 6, 1, 3, 6, 2.0f);
    for (size_t i = 0; i < original.size(); ++i) {
        EXPECT_NEAR(original[i], density.values[i], 1e-4f);
    }
}

UTEST(analysis_density, triangle_edges_do_not_wrap_and_smoothing_preserves_mass) {
    const vec2_t point[] = {{0,0}};
    Triangle density;
    density.update(point, 1, 1, 0, 1, 1.0f);
    const float peak = density.values[0];
    EXPECT_GT(peak, 0.0f);
    EXPECT_EQ(0.0f, density.values[Triangle::index(Triangle::resolution,0)]);
    EXPECT_EQ(0.0f, density.values[Triangle::index(0,Triangle::resolution)]);
    density.update(point, 1, 1, 0, 1, 4.0f);
    EXPECT_LT(density.values[0], peak);
    double sum = 0;
    for (float v : density.values) sum += v;
    EXPECT_NEAR(1.0, sum / density.values.size(), 1e-5);
    const auto texels = density.texture_data();
    EXPECT_EQ((size_t)65 * 65, texels.size());
    EXPECT_GT(texels[0], 0.0f);
    EXPECT_EQ(0.0f, texels[64]);
    for (float v : texels) {
        EXPECT_TRUE(std::isfinite(v));
        EXPECT_GE(v, 0.0f);
        EXPECT_LE(v, 1.0f);
    }
}

UTEST(analysis_density, empty_subset_and_invalid_samples_are_empty) {
    const vec2_t points[] = {{std::numeric_limits<float>::quiet_NaN(),0}, {5,5}};
    Triangle density;
    density.update(points, 2, 1, 0, 2, 2.0f);
    for (float v : density.values) EXPECT_EQ(0.0f, v);
    const vec2_t valid[] = {{0.5f,0.25f}};
    density.update(valid, 1, 1, 0, 0, 2.0f);
    for (float v : density.values) EXPECT_EQ(0.0f, v);
    density.clear();
    for (float v : density.texture_data()) EXPECT_EQ(0.0f, v);
}
