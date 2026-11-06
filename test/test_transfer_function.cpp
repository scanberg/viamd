#include "utest.h"

#include <transfer_function.h>

// The transfer function's arithmetic (transfer_function.h): the curve through its points, the alpha baked from it,
// editing the points, the range. The editor is ImGui and ImPlot and is not tested here.

UTEST(transfer_function, curve_is_piecewise_linear_and_flat_beyond_its_ends) {
    TransferFunction tf;
    tf.num_points = 0;
    transfer_function_insert(&tf, 0.2f, 0.0f);
    transfer_function_insert(&tf, 0.6f, 1.0f);
    EXPECT_EQ(0.0f, transfer_function_curve(tf, 0.0f));     // before the first: its value
    EXPECT_EQ(0.0f, transfer_function_curve(tf, 0.2f));
    EXPECT_NEAR(0.5f, transfer_function_curve(tf, 0.4f), 1e-6f);
    EXPECT_EQ(1.0f, transfer_function_curve(tf, 0.6f));
    EXPECT_EQ(1.0f, transfer_function_curve(tf, 1.0f));     // after the last: its value
}

UTEST(transfer_function, alpha_is_the_curve_times_the_opacity) {
    TransferFunction tf;    // 0 at t = 0 rising to 1 at t = 1
    tf.opacity = 0.2f;
    float a[5];
    transfer_function_bake_alpha(tf, a, 5, 1);
    EXPECT_EQ(0.0f, a[0]);
    EXPECT_NEAR(0.1f, a[2], 1e-6f);
    EXPECT_NEAR(0.2f, a[4], 1e-6f);

    // into RGBA, every fourth
    float rgba[3 * 4] = {};
    transfer_function_bake_alpha(tf, rgba + 3, 3, 4);
    EXPECT_NEAR(0.1f, rgba[7], 1e-6f);
    EXPECT_EQ(0.0f, rgba[4]);   // the colours left alone
}

UTEST(transfer_function, points_stay_in_order_and_in_range) {
    TransferFunction tf;
    const int i = transfer_function_insert(&tf, 0.5f, 2.0f);   // between the two, alpha clamped
    EXPECT_EQ(1, i);
    EXPECT_EQ(3, tf.num_points);
    EXPECT_EQ(1.0f, tf.points[1].alpha);

    // moved past a neighbour: it stops at it
    transfer_function_move(&tf, 1, 1.5f, 0.3f);
    EXPECT_EQ(tf.points[2].t, tf.points[1].t);
    EXPECT_EQ(0.3f, tf.points[1].alpha);

    // the last one is not removed
    EXPECT_TRUE(transfer_function_remove(&tf, 0));
    EXPECT_TRUE(transfer_function_remove(&tf, 0));
    EXPECT_FALSE(transfer_function_remove(&tf, 0));
    EXPECT_EQ(1, tf.num_points);

    // full
    for (int k = 0; k < TRANSFER_FUNCTION_MAX_POINTS + 2; ++k) transfer_function_insert(&tf, 0.01f * (float)k, 0.5f);
    EXPECT_EQ(TRANSFER_FUNCTION_MAX_POINTS, tf.num_points);
    EXPECT_EQ(-1, transfer_function_insert(&tf, 0.5f, 0.5f));
}

UTEST(transfer_function, range_and_param) {
    TransferFunction tf;
    transfer_function_set_range(&tf, -0.1f, 0.3f);
    EXPECT_NEAR(0.25f, transfer_function_param(tf, 0.0f), 1e-6f);
    EXPECT_EQ(0.0f, transfer_function_param(tf, -1.0f));  // clamped
    EXPECT_NEAR(0.3f, transfer_function_value(tf, 1.0f), 1e-6f);

    // symmetric: about 0, out to the larger magnitude
    tf.symmetric = true;
    transfer_function_set_range(&tf, -0.1f, 0.3f);
    EXPECT_NEAR(-0.3f, tf.range_beg, 1e-6f);
    EXPECT_NEAR(0.3f, tf.range_end, 1e-6f);
    EXPECT_NEAR(0.5f, transfer_function_param(tf, 0.0f), 1e-6f);

    // of no extent: every value in the middle, not a division by 0
    transfer_function_set_range(&tf, 0.0f, 0.0f);
    EXPECT_EQ(0.5f, transfer_function_param(tf, 3.0f));
}

UTEST(transfer_function, ramps_and_hash) {
    TransferFunction tf;
    transfer_function_ramp_magnitude(&tf, 0.2f);   // 1 at both ends, 0 within 0.1 of the middle
    EXPECT_EQ(1.0f, transfer_function_curve(tf, 0.0f));
    EXPECT_EQ(0.0f, transfer_function_curve(tf, 0.5f));
    EXPECT_EQ(0.0f, transfer_function_curve(tf, 0.6f));
    EXPECT_NEAR(0.5f, transfer_function_curve(tf, 0.8f), 1e-6f);

    const uint64_t h = transfer_function_hash(tf);
    transfer_function_ramp(&tf, 0.5f);
    EXPECT_NE(h, transfer_function_hash(tf));
    EXPECT_EQ(0.0f, transfer_function_curve(tf, 0.5f));
    EXPECT_NEAR(0.5f, transfer_function_curve(tf, 0.75f), 1e-6f);
}
