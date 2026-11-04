#include "utest.h"

#include <plane.h>

#include <math.h>

// The plane's geometry (plane.h): the fit through points, the polygon a plane cuts out of a box, and the axes laid
// out on it. The gizmo and the drawing are ImGui and OpenGL and are not tested here.

static float max_distance(const Plane& p, const vec3_t* pts, size_t n) {
    float m = 0.0f;
    for (size_t i = 0; i < n; ++i) m = fmaxf(m, fabsf(plane_distance(p, pts[i])));
    return m;
}

UTEST(plane, fit_three_points_is_exact) {
    const vec3_t pts[3] = {{1.3f, -0.2f, 0.7f}, {-2.1f, 0.9f, 1.4f}, {0.4f, 2.6f, -1.1f}};
    const Plane p = plane_fit(Plane{}, pts, 3);
    EXPECT_LT(max_distance(p, pts, 3), 1.0e-5f);
    EXPECT_NEAR(1.0f, vec3_length(p.normal), 1.0e-6f);
}

UTEST(plane, fit_keeps_the_side_of_the_previous_plane) {
    // z = 2, as the plane facing -z: the fit faces -z too, so what it cut away stays cut away
    const vec3_t pts[4] = {{0, 0, 2}, {1, 0, 2}, {0, 1, 2}, {3, -2, 2}};
    const Plane prev = {{0.1f, 0.0f, -1.0f}, 5.0f};
    const Plane p = plane_fit(prev, pts, 4);
    EXPECT_NEAR(-1.0f, p.normal.z, 1.0e-6f);
    EXPECT_NEAR(2.0f, p.d, 1.0e-5f);   // -n . x = -(-1)(2)

    const Plane q = plane_fit(Plane{{0, 0, 1}, 0}, pts, 4);
    EXPECT_NEAR(1.0f, q.normal.z, 1.0e-6f);
    EXPECT_NEAR(-2.0f, q.d, 1.0e-5f);
}

UTEST(plane, fit_many_points_is_the_least_squares_plane) {
    // points on a tilted plane, pushed off it symmetrically: the fit is the plane itself
    const vec3_t n = vec3_normalize(vec3_set(0.3f, -0.5f, 0.8f));
    const vec3_t o = vec3_set(1.0f, 2.0f, -3.0f);
    vec3_t u, v;
    plane_frame(&u, &v, Plane{n, 0.0f}, vec3_set(1, 0, 0));
    vec3_t pts[40];
    size_t k = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) {
            const vec3_t in_plane = vec3_add(o, vec3_add(vec3_mul1(u, (float)i - 2.0f), vec3_mul1(v, 1.5f * (float)j - 2.25f)));
            pts[k++] = vec3_add(in_plane, vec3_mul1(n, 0.05f));
            pts[k++] = vec3_sub(in_plane, vec3_mul1(n, 0.05f));
        }
    }
    const Plane p = plane_fit(Plane{n, 0.0f}, pts, k);
    EXPECT_GT(vec3_dot(p.normal, n), 0.99999f);
    EXPECT_LT(fabsf(plane_distance(p, o)), 1.0e-4f);
}

UTEST(plane, fit_a_line_turns_the_previous_plane_onto_it) {
    // along x: the plane through the line nearest to the previous one, its normal without the part along x
    const vec3_t pts[3] = {{-1, 1, 1}, {0, 1, 1}, {2, 1, 1}};
    const Plane prev = {vec3_normalize(vec3_set(0.5f, 1.0f, 1.0f)), 0.0f};
    const Plane p = plane_fit(prev, pts, 3);
    const vec3_t want = vec3_normalize(vec3_set(0.0f, 1.0f, 1.0f));
    EXPECT_GT(vec3_dot(p.normal, want), 0.99999f);
    EXPECT_LT(max_distance(p, pts, 3), 1.0e-5f);

    // two points are a line too
    const Plane q = plane_fit(prev, pts, 2);
    EXPECT_GT(vec3_dot(q.normal, want), 0.99999f);
    EXPECT_LT(max_distance(q, pts, 2), 1.0e-5f);

    // a previous plane across the line: still some plane through it
    const Plane r = plane_fit(Plane{{1, 0, 0}, 0}, pts, 3);
    EXPECT_LT(fabsf(r.normal.x), 1.0e-4f);
    EXPECT_LT(max_distance(r, pts, 3), 1.0e-5f);
}

UTEST(plane, fit_one_point_moves_the_plane_through_it) {
    const vec3_t pt = {1.0f, -2.0f, 4.0f};
    const Plane prev = {vec3_normalize(vec3_set(1, 1, 0)), 3.0f};
    const Plane p = plane_fit(prev, &pt, 1);
    EXPECT_GT(vec3_dot(p.normal, prev.normal), 0.999999f);
    EXPECT_LT(fabsf(plane_distance(p, pt)), 1.0e-5f);

    // and no points leave it as it was
    const Plane q = plane_fit(prev, &pt, 0);
    EXPECT_EQ(prev.d, q.d);
}

UTEST(plane, box_polygon) {
    const mat4_t box = mat4_ident();
    vec3_t poly[6];
    EXPECT_EQ(4, plane_box_polygon(poly, plane_from_point_normal(vec3_set(0.5f, 0.5f, 0.5f), vec3_set(0, 0, 1)), box));
    // through the middle across the main diagonal: a regular hexagon
    const Plane diag = plane_from_point_normal(vec3_set(0.5f, 0.5f, 0.5f), vec3_set(1, 1, 1));
    const int n = plane_box_polygon(poly, diag, box);
    EXPECT_EQ(6, n);
    for (int i = 0; i < n; ++i) {
        EXPECT_LT(fabsf(plane_distance(diag, poly[i])), 1.0e-5f);
        // in order around it: neighbours are a side apart, sqrt(2)/2
        EXPECT_NEAR(0.70710678f, vec3_distance(poly[i], poly[(i + 1) % n]), 1.0e-5f);
    }
    // missing it
    EXPECT_EQ(0, plane_box_polygon(poly, plane_from_point_normal(vec3_set(0, 0, 2), vec3_set(0, 0, 1)), box));
}

UTEST(plane, frame_follows_right_and_is_right_handed) {
    const Plane p = {vec3_normalize(vec3_set(0.2f, -0.3f, -0.93f)), 1.0f};
    vec3_t u, v;
    plane_frame(&u, &v, p, vec3_set(1, 0, 0));
    EXPECT_NEAR(1.0f, vec3_length(u), 1.0e-6f);
    EXPECT_NEAR(0.0f, vec3_dot(u, p.normal), 1.0e-6f);
    EXPECT_NEAR(0.0f, vec3_dot(u, v), 1.0e-6f);
    EXPECT_GT(vec3_dot(vec3_cross(u, v), p.normal), 0.99999f);
    EXPECT_GT(u.x, 0.9f);   // as near to 'right' as the plane allows

    // 'right' along the normal: another axis, still a frame
    plane_frame(&u, &v, Plane{{1, 0, 0}, 0}, vec3_set(2, 0, 0));
    EXPECT_NEAR(1.0f, vec3_length(u), 1.0e-6f);
    EXPECT_NEAR(0.0f, u.x, 1.0e-6f);
}
