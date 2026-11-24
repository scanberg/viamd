#include "utest.h"

#include <clock.h>

// Two runs of one simulation: every frame (200 ps apart) and every fourth (800 ps), placed on a clock in ps. At a time,
// each shows its frame of that time, not the frame of the same number.
UTEST(clock, runs_sampled_differently_meet_in_time) {
    double every[501], fourth[126];
    for (int i = 0; i < 501; ++i) every[i]  = 200.0 * i;
    for (int i = 0; i < 126; ++i) fourth[i] = 800.0 * i;
    const ClockAxis a = {every, 501, 1.0};
    const ClockAxis b = {fourth, 126, 1.0};

    const double t = clock_axis_time(a, 40.0);  // 8000 ps
    EXPECT_NEAR(8000.0, t, 1e-9);
    EXPECT_NEAR(10.0, clock_axis_frame(b, t), 1e-9);
    // Between two of b's frames: the fraction of the way
    EXPECT_NEAR(10.25, clock_axis_frame(b, clock_axis_time(a, 41.0)), 1e-9);
    EXPECT_NEAR(200.0, clock_axis_frame_duration(a), 1e-9);
    EXPECT_NEAR(800.0, clock_axis_frame_duration(b), 1e-9);
}

// A run in ns on a clock in ps: its times taken into the clock's unit
UTEST(clock, scaled_into_the_clock_unit) {
    const double ns[3] = {0.0, 0.5, 1.0};
    const ClockAxis a = {ns, 3, 1000.0};
    EXPECT_NEAR(500.0, clock_axis_time(a, 1.0), 1e-9);
    EXPECT_NEAR(1.5, clock_axis_frame(a, 750.0), 1e-9);
    // Clamped to its frames: before its first, after its last
    EXPECT_NEAR(0.0, clock_axis_frame(a, -10.0), 1e-9);
    EXPECT_NEAR(2.0, clock_axis_frame(a, 5000.0), 1e-9);
}

UTEST(clock, untimed) {
    const double t[2] = {0.0, 1.0};
    EXPECT_FALSE(clock_axis_timed(ClockAxis{t, 2, 0.0}));
    EXPECT_FALSE(clock_axis_timed(ClockAxis{nullptr, 0, 1.0}));
    EXPECT_TRUE(clock_axis_timed(ClockAxis{t, 2, 1.0}));
    // Fewer than two frames, or none between them: a frame is one unit of time
    EXPECT_NEAR(1.0, clock_axis_frame_duration(ClockAxis{t, 1, 1.0}), 1e-12);
    const double same[2] = {5.0, 5.0};
    EXPECT_NEAR(1.0, clock_axis_frame_duration(ClockAxis{same, 2, 1.0}), 1e-12);
}

// Stepping by the reference's frames: to its own frames within them, by its mean time per frame past them, within the
// span of every run shown
UTEST(clock, step) {
    double ref[100];
    for (int i = 0; i < 100; ++i) ref[i] = 10.0 * i;     // 0 .. 990
    const ClockAxis r = {ref, 100, 1.0};

    EXPECT_NEAR(20.0, clock_step_time(r, 10.0, 1.0, 0.0, 2000.0), 1e-9);
    EXPECT_NEAR(110.0, clock_step_time(r, 10.0, 10.0, 0.0, 2000.0), 1e-9);
    EXPECT_NEAR(0.0, clock_step_time(r, 10.0, -1.0, 0.0, 2000.0), 1e-9);
    // Fractional frames stay fractional
    EXPECT_NEAR(25.0, clock_step_time(r, 15.0, 1.0, 0.0, 2000.0), 1e-9);
    // Past its last frame, a longer run goes on: by the mean time per frame, and back again
    EXPECT_NEAR(1000.0, clock_step_time(r, 990.0, 1.0, 0.0, 2000.0), 1e-9);
    EXPECT_NEAR(1510.0, clock_step_time(r, 1500.0, 1.0, 0.0, 2000.0), 1e-9);
    EXPECT_NEAR(1490.0, clock_step_time(r, 1500.0, -1.0, 0.0, 2000.0), 1e-9);
    // Within the span
    EXPECT_NEAR(2000.0, clock_step_time(r, 1995.0, 1.0, 0.0, 2000.0), 1e-9);
    EXPECT_NEAR(0.0, clock_step_time(r, 0.0, -1.0, 0.0, 2000.0), 1e-9);
}
