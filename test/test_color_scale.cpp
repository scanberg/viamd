#include "utest.h"

#include <color_scale.h>

// The colour scale's arithmetic: what an automatic range takes from a span, and where a value
// lands on the scale. The controls and the legend are ImGui and are not tested here.

UTEST(color_scale, auto_range_takes_the_span) {
    ColorScale scale;
    scale.auto_range = true;
    scale.symmetric  = false;
    ColorScaleSpan span = color_scale_span(-0.2f, 0.6f);
    color_scale_update_range(&scale, span);
    EXPECT_EQ(-0.2f, scale.range_beg);
    EXPECT_EQ( 0.6f, scale.range_end);

    // Symmetric: centred on zero, reaching the larger magnitude
    scale.symmetric = true;
    color_scale_update_range(&scale, span);
    EXPECT_EQ(-0.6f, scale.range_beg);
    EXPECT_EQ( 0.6f, scale.range_end);
}

UTEST(color_scale, auto_range_takes_the_robust_part) {
    // lo / hi, not the extremes: a spike next to a surface must not become the range
    ColorScaleSpan span;
    span.valid = true;
    span.min = -5.0f; span.lo = 0.15f;
    span.max =  0.4f; span.hi = 0.31f;
    ColorScale scale;
    scale.symmetric = false;
    color_scale_update_range(&scale, span);
    EXPECT_EQ(0.15f, scale.range_beg);
    EXPECT_EQ(0.31f, scale.range_end);
}

UTEST(color_scale, manual_range_and_empty_span_are_left_alone) {
    ColorScale scale;
    scale.range_beg = 1.0f;
    scale.range_end = 2.0f;
    scale.auto_range = false;
    color_scale_update_range(&scale, color_scale_span(-3.0f, 3.0f));
    EXPECT_EQ(1.0f, scale.range_beg);
    EXPECT_EQ(2.0f, scale.range_end);

    scale.auto_range = true;
    color_scale_update_range(&scale, ColorScaleSpan{});
    EXPECT_EQ(1.0f, scale.range_beg);
    EXPECT_EQ(2.0f, scale.range_end);

    EXPECT_FALSE(color_scale_span(1.0f, -1.0f).valid);
}

UTEST(color_scale, values_land_on_the_range) {
    ColorScale scale;
    scale.range_beg = -1.0f;
    scale.range_end =  1.0f;
    EXPECT_NEAR(0.0f,  color_scale_param(scale, -1.0f), 1e-6f);
    EXPECT_NEAR(0.5f,  color_scale_param(scale,  0.0f), 1e-6f);
    EXPECT_NEAR(0.75f, color_scale_param(scale,  0.5f), 1e-6f);
    EXPECT_EQ(0.0f, color_scale_param(scale, -7.0f));   // clamped
    EXPECT_EQ(1.0f, color_scale_param(scale,  7.0f));

    // Small values are values: no floor on the extent (the old atom colouring had one of 0.001)
    scale.range_beg = 0.0f;
    scale.range_end = 1e-5f;
    EXPECT_NEAR(0.5f, color_scale_param(scale, 5e-6f), 1e-4f);

    // Backwards reverses the map; no extent is the middle
    scale.range_beg = 1.0f;
    scale.range_end = 0.0f;
    EXPECT_NEAR(0.25f, color_scale_param(scale, 0.75f), 1e-6f);
    scale.range_end = 1.0f;
    EXPECT_EQ(0.5f, color_scale_param(scale, 3.0f));
}

UTEST(color_scale, signed_spans) {
    EXPECT_TRUE (color_scale_span_signed(color_scale_span(-0.8f, 0.4f)));
    EXPECT_FALSE(color_scale_span_signed(color_scale_span( 0.0f, 0.4f)));
    EXPECT_FALSE(color_scale_span_signed(color_scale_span( 1.0f, 16.0f)));
    EXPECT_FALSE(color_scale_span_signed(ColorScaleSpan{}));
}
