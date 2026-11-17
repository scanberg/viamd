#include "utest.h"

#include <plane.h>

#include <math.h>
#include <core/md_vec_math.h>

// The plane's geometry (plane.h): the fit through points, the polygon a plane cuts out of a box, and the axes laid
// out on it. The gizmo and the drawing are ImGui and OpenGL and are not tested here.

static float max_distance(const Plane& p, const md_float3_t* pts, size_t n) {
    float m = 0.0f;
    for (size_t i = 0; i < n; ++i) m = fmaxf(m, fabsf(plane_distance(p, pts[i])));
    return m;
}

UTEST(plane, fit_three_points_is_exact) {
    const md_float3_t pts[3] = {{1.3f, -0.2f, 0.7f}, {-2.1f, 0.9f, 1.4f}, {0.4f, 2.6f, -1.1f}};
    const Plane p = plane_fit(Plane{}, pts, 3);
    EXPECT_LT(max_distance(p, pts, 3), 1.0e-5f);
    EXPECT_NEAR(1.0f, md_float3_length(p.normal), 1.0e-6f);
}

UTEST(plane, fit_keeps_the_side_of_the_previous_plane) {
    // z = 2, as the plane facing -z: the fit faces -z too, so what it cut away stays cut away
    const md_float3_t pts[4] = {{0, 0, 2}, {1, 0, 2}, {0, 1, 2}, {3, -2, 2}};
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
    const md_float3_t n = md_float3_normalize(md_float3_set(0.3f, -0.5f, 0.8f));
    const md_float3_t o = md_float3_set(1.0f, 2.0f, -3.0f);
    md_float3_t u, v;
    plane_frame(&u, &v, Plane{n, 0.0f}, md_float3_set(1, 0, 0));
    md_float3_t pts[40];
    size_t k = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) {
            const md_float3_t in_plane = md_float3_add(o, md_float3_add(md_float3_mul1(u, (float)i - 2.0f), md_float3_mul1(v, 1.5f * (float)j - 2.25f)));
            pts[k++] = md_float3_add(in_plane, md_float3_mul1(n, 0.05f));
            pts[k++] = md_float3_sub(in_plane, md_float3_mul1(n, 0.05f));
        }
    }
    const Plane p = plane_fit(Plane{n, 0.0f}, pts, k);
    EXPECT_GT(md_float3_dot(p.normal, n), 0.99999f);
    EXPECT_LT(fabsf(plane_distance(p, o)), 1.0e-4f);
}

UTEST(plane, fit_a_line_turns_the_previous_plane_onto_it) {
    // along x: the plane through the line nearest to the previous one, its normal without the part along x
    const md_float3_t pts[3] = {{-1, 1, 1}, {0, 1, 1}, {2, 1, 1}};
    const Plane prev = {md_float3_normalize(md_float3_set(0.5f, 1.0f, 1.0f)), 0.0f};
    const Plane p = plane_fit(prev, pts, 3);
    const md_float3_t want = md_float3_normalize(md_float3_set(0.0f, 1.0f, 1.0f));
    EXPECT_GT(md_float3_dot(p.normal, want), 0.99999f);
    EXPECT_LT(max_distance(p, pts, 3), 1.0e-5f);

    // two points are a line too
    const Plane q = plane_fit(prev, pts, 2);
    EXPECT_GT(md_float3_dot(q.normal, want), 0.99999f);
    EXPECT_LT(max_distance(q, pts, 2), 1.0e-5f);

    // a previous plane across the line: still some plane through it
    const Plane r = plane_fit(Plane{{1, 0, 0}, 0}, pts, 3);
    EXPECT_LT(fabsf(r.normal.x), 1.0e-4f);
    EXPECT_LT(max_distance(r, pts, 3), 1.0e-5f);
}

UTEST(plane, fit_one_point_moves_the_plane_through_it) {
    const md_float3_t pt = {1.0f, -2.0f, 4.0f};
    const Plane prev = {md_float3_normalize(md_float3_set(1, 1, 0)), 3.0f};
    const Plane p = plane_fit(prev, &pt, 1);
    EXPECT_GT(md_float3_dot(p.normal, prev.normal), 0.999999f);
    EXPECT_LT(fabsf(plane_distance(p, pt)), 1.0e-5f);

    // and no points leave it as it was
    const Plane q = plane_fit(prev, &pt, 0);
    EXPECT_EQ(prev.d, q.d);
}

UTEST(plane, box_polygon) {
    const md_float4x4_t box = md_float4x4_ident();
    md_float3_t poly[6];
    EXPECT_EQ(4, plane_box_polygon(poly, plane_from_point_normal(md_float3_set(0.5f, 0.5f, 0.5f), md_float3_set(0, 0, 1)), box));
    // through the middle across the main diagonal: a regular hexagon
    const Plane diag = plane_from_point_normal(md_float3_set(0.5f, 0.5f, 0.5f), md_float3_set(1, 1, 1));
    const int n = plane_box_polygon(poly, diag, box);
    EXPECT_EQ(6, n);
    for (int i = 0; i < n; ++i) {
        EXPECT_LT(fabsf(plane_distance(diag, poly[i])), 1.0e-5f);
        // in order around it: neighbours are a side apart, sqrt(2)/2
        EXPECT_NEAR(0.70710678f, md_float3_distance(poly[i], poly[(i + 1) % n]), 1.0e-5f);
    }
    // missing it
    EXPECT_EQ(0, plane_box_polygon(poly, plane_from_point_normal(md_float3_set(0, 0, 2), md_float3_set(0, 0, 1)), box));
}

UTEST(plane, frame_follows_right_and_is_right_handed) {
    const Plane p = {md_float3_normalize(md_float3_set(0.2f, -0.3f, -0.93f)), 1.0f};
    md_float3_t u, v;
    plane_frame(&u, &v, p, md_float3_set(1, 0, 0));
    EXPECT_NEAR(1.0f, md_float3_length(u), 1.0e-6f);
    EXPECT_NEAR(0.0f, md_float3_dot(u, p.normal), 1.0e-6f);
    EXPECT_NEAR(0.0f, md_float3_dot(u, v), 1.0e-6f);
    EXPECT_GT(md_float3_dot(md_float3_cross(u, v), p.normal), 0.99999f);
    EXPECT_GT(u.x, 0.9f);   // as near to 'right' as the plane allows

    // 'right' along the normal: another axis, still a frame
    plane_frame(&u, &v, Plane{{1, 0, 0}, 0}, md_float3_set(2, 0, 0));
    EXPECT_NEAR(1.0f, md_float3_length(u), 1.0e-6f);
    EXPECT_NEAR(0.0f, u.x, 1.0e-6f);
}

UTEST(plane, fit_reports_how_well_the_points_determine_it) {
    PlaneFitInfo info;
    const md_float3_t ring[6] = {{1, 0, 0}, {0.5f, 0.866f, 0}, {-0.5f, 0.866f, 0}, {-1, 0, 0}, {-0.5f, -0.866f, 0}, {0.5f, -0.866f, 0}};
    plane_fit(Plane{}, ring, 6, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::Fit);
    EXPECT_LT(info.rms, 1.0e-5f);

    plane_fit(Plane{}, ring, 1, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::Point);
    plane_fit(Plane{}, ring, 0, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::None);
}

UTEST(plane, fit_a_tetrahedron_keeps_the_normal) {
    // spread alike every way (methane's hydrogens): no plane fits better than another, the normal is kept
    const md_float3_t tet[4] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    const Plane prev = {md_float3_normalize(md_float3_set(0.3f, 0.4f, 0.87f)), 7.0f};
    PlaneFitInfo info;
    const Plane p = plane_fit(prev, tet, 4, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::Kept);
    EXPECT_GT(md_float3_dot(p.normal, prev.normal), 0.999999f);
    EXPECT_LT(fabsf(p.d), 1.0e-5f);   // through the centroid, the origin
}

UTEST(plane, fit_a_loose_set_turns_as_little_as_it_allows) {
    // long along x, alike along y and z: every normal in the yz plane fits as well, and the one nearest the previous
    // normal is taken - not whichever of them the eigen solver gives, which turns with the smallest change
    const md_float3_t pts[6] = {{-3, 0, 0}, {3, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    const Plane prev = {md_float3_normalize(md_float3_set(0.5f, 0.6f, -0.8f)), 0.0f};
    PlaneFitInfo info;
    const Plane p = plane_fit(prev, pts, 6, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::Loose);
    const md_float3_t want = md_float3_normalize(md_float3_set(0.0f, 0.6f, -0.8f));
    EXPECT_GT(md_float3_dot(p.normal, want), 0.99999f);

    // squashed along z, it is a plane again: the xy one, whatever the previous normal
    md_float3_t flat[6];
    for (int i = 0; i < 6; ++i) flat[i] = md_float3_set(pts[i].x, pts[i].y, pts[i].z * 0.5f);
    const Plane q = plane_fit(prev, flat, 6, &info);
    EXPECT_TRUE(info.kind == PlaneFitKind::Fit);
    EXPECT_GT(fabsf(q.normal.z), 0.99999f);
    EXPECT_LT(q.normal.z, 0.0f);   // on prev's side
}
