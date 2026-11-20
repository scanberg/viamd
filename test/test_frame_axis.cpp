#include "utest.h"

#include <frame_axis.h>

#include <math.h>

// A frame between two others lies on the line between their times, and back
UTEST(frame_axis, between_frames) {
    const double t[] = {0.0, 10.0, 30.0, 60.0};
    const size_t n = 4;

    EXPECT_EQ(0.0,  time_at_frame(t, n, 0.0));
    EXPECT_EQ(30.0, time_at_frame(t, n, 2.0));
    EXPECT_EQ(5.0,  time_at_frame(t, n, 0.5));
    EXPECT_EQ(45.0, time_at_frame(t, n, 2.5));

    EXPECT_EQ(2.0,  frame_at_time(t, n, 30.0));
    EXPECT_EQ(0.5,  frame_at_time(t, n, 5.0));
    EXPECT_EQ(2.5,  frame_at_time(t, n, 45.0));
    EXPECT_EQ(3.0,  frame_at_time(t, n, 60.0));

    // Each other's inverse where the time increases
    for (double f = 0.0; f <= 3.0; f += 0.125) {
        EXPECT_NEAR(f, frame_at_time(t, n, time_at_frame(t, n, f)), 1e-12);
    }
}

// Outside the axis, its ends; nothing on an empty axis
UTEST(frame_axis, clamped) {
    const double t[] = {100.0, 110.0, 120.0};
    EXPECT_EQ(100.0, time_at_frame(t, 3, -1.0));
    EXPECT_EQ(120.0, time_at_frame(t, 3, 7.0));
    EXPECT_EQ(0.0,   frame_at_time(t, 3, 50.0));
    EXPECT_EQ(2.0,   frame_at_time(t, 3, 500.0));
    EXPECT_EQ(0.0,   frame_at_time(t, 3, NAN));
    EXPECT_EQ(100.0, time_at_frame(t, 3, NAN));

    const double one[] = {42.0};
    EXPECT_EQ(42.0, time_at_frame(one, 1, 0.3));
    EXPECT_EQ(0.0,  frame_at_time(one, 1, 42.0));
    EXPECT_EQ(0.0,  frame_at_time(one, 1, 43.0));

    EXPECT_EQ(0.0, time_at_frame(nullptr, 0, 1.0));
    EXPECT_EQ(0.0, frame_at_time(nullptr, 0, 1.0));
}

// Frames which share a time: the first of them, and never a division by the zero between them
UTEST(frame_axis, shared_times) {
    const double t[] = {0.0, 10.0, 10.0, 10.0, 20.0, 20.0};
    const size_t n = 6;
    EXPECT_EQ(1.0, frame_at_time(t, n, 10.0));
    EXPECT_EQ(0.5, frame_at_time(t, n, 5.0));
    EXPECT_EQ(3.5, frame_at_time(t, n, 15.0));
    EXPECT_EQ(4.0, frame_at_time(t, n, 20.0));
    EXPECT_EQ(10.0, time_at_frame(t, n, 2.5));
    for (double x = 0.0; x <= 20.0; x += 0.5) {
        const double f = frame_at_time(t, n, x);
        EXPECT_TRUE(f == f);
        EXPECT_NEAR(x, time_at_frame(t, n, f), 1e-12);
    }
}

// The time of a long run in ps, where a float would have lost the fraction of a frame
UTEST(frame_axis, long_runs_keep_the_fraction) {
    enum { N = 1000 };
    double t[N];
    for (int i = 0; i < N; ++i) t[i] = 1.0e7 + 2.0 * i;    // 10 us, a frame every 2 ps
    EXPECT_NEAR(500.25, frame_at_time(t, N, time_at_frame(t, N, 500.25)), 1e-9);
    EXPECT_EQ(1.0e7 + 1000.5, time_at_frame(t, N, 500.25));
}
