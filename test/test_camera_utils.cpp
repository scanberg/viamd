#include "utest.h"

#include <gfx/camera_utils.h>
#include <gfx/camera.h>

#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>

/* camera_utils is where a sign or a transposition goes unnoticed for a long time: a projection that
 * is subtly wrong still renders something, and the picture only looks off at the extremes. So these
 * tests check the properties rather than the entries - a matrix and its stated inverse must compose
 * to the identity, and a point sent through a projection must come back.
 *
 * Every camera here is off-axis and off-origin. With an identity orientation at the origin most of
 * these compositions hold for the wrong reasons. */

static float max_abs_diff(md_float4x4_t A, md_float4x4_t B) {
    float m = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            const float d = fabsf(A.elem[i][j] - B.elem[i][j]);
            if (d > m) m = d;
        }
    return m;
}

static Camera test_camera(void) {
    Camera c;
    c.orientation = md_quat_normalize(md_quat_axis_angle(md_float3_normalize(md_float3_set(0.3f, 1.0f, -0.2f)), 0.7f));
    c.position    = md_float3_set(2.5f, -1.25f, 8.0f);
    c.distance    = 12.0f;
    c.near_plane  = 0.5f;
    c.far_plane   = 250.0f;
    c.fov_y       = 1.0471975512f;   /* 60 degrees */
    return c;
}

UTEST(viamd_camera, view_and_world_transforms_are_inverses) {
    const Camera c = test_camera();

    const md_float4x4_t w2v = camera_world_to_view_matrix(c);
    const md_float4x4_t v2w = camera_view_to_world_matrix(c);

    EXPECT_LT(max_abs_diff(md_float4x4_mul(w2v, v2w), md_float4x4_ident()), 1.0e-4f);
    EXPECT_LT(max_abs_diff(md_float4x4_mul(v2w, w2v), md_float4x4_ident()), 1.0e-4f);

    /* The camera sits at the origin of its own view space. */
    const md_float4_t eye = md_float4x4_mul_float4(w2v, md_float4_set(c.position.x, c.position.y, c.position.z, 1.0f));
    EXPECT_NEAR(0.0f, eye.x, 1.0e-4f);
    EXPECT_NEAR(0.0f, eye.y, 1.0e-4f);
    EXPECT_NEAR(0.0f, eye.z, 1.0e-4f);
}

UTEST(viamd_camera, perspective_clip_and_view_are_inverses) {
    const Camera c = test_camera();
    const float aspect = 16.0f / 9.0f;

    const md_float4x4_t v2c = camera_view_to_clip_matrix_persp(c, aspect);
    const md_float4x4_t c2v = camera_clip_to_view_matrix_persp(c, aspect);

    EXPECT_LT(max_abs_diff(md_float4x4_mul(v2c, c2v), md_float4x4_ident()), 1.0e-4f);
    EXPECT_LT(max_abs_diff(md_float4x4_mul(c2v, v2c), md_float4x4_ident()), 1.0e-4f);
}

UTEST(viamd_camera, perspective_maps_the_clip_planes_to_the_depth_range) {
    const Camera c = test_camera();
    const md_float4x4_t v2c = camera_view_to_clip_matrix_persp(c, 1.0f);

    /* Looking down -z in view space. The near plane lands at one end of the depth range and the far
     * plane at the other. The span is what identifies the convention: 2 is OpenGL's [-1,1] clip
     * space, 1 would be the [0,1] that D3D and Vulkan use. mdlib's md_float4x4_persp builds the GL form and
     * viamd's shaders read it back that way, so the two have to keep agreeing. */
    const md_float4_t n = md_float4x4_mul_float4(v2c, md_float4_set(0, 0, -c.near_plane, 1.0f));
    const md_float4_t f = md_float4x4_mul_float4(v2c, md_float4_set(0, 0, -c.far_plane,  1.0f));

    ASSERT_GT(fabsf(n.w), 1.0e-6f);
    ASSERT_GT(fabsf(f.w), 1.0e-6f);
    const float zn = n.z / n.w;
    const float zf = f.z / f.w;

    EXPECT_TRUE(isfinite(zn));
    EXPECT_TRUE(isfinite(zf));
    EXPECT_NEAR(2.0f, fabsf(zf - zn), 1.0e-3f);   /* [-1,1], not [0,1] */
    EXPECT_NEAR(1.0f, fabsf(zn), 1.0e-3f);
    EXPECT_NEAR(1.0f, fabsf(zf), 1.0e-3f);

    /* A point on the axis stays on the axis. */
    EXPECT_NEAR(0.0f, n.x / n.w, 1.0e-5f);
    EXPECT_NEAR(0.0f, n.y / n.w, 1.0e-5f);
}

UTEST(viamd_camera, perspective_by_aspect_and_by_resolution_agree) {
    const Camera c = test_camera();

    /* Two spellings of the same projection: one takes an aspect ratio, the other a pixel size with
     * a jitter offset. With no jitter they have to produce the same matrix, or temporal
     * antialiasing silently renders a different frustum from everything else. */
    const md_float4x4_t by_aspect = camera_view_to_clip_matrix_persp(c, 1920.0f / 1080.0f);
    const md_float4x4_t by_size   = camera_view_to_clip_matrix_persp(c, 1920, 1080, 0.0f, 0.0f);

    EXPECT_LT(max_abs_diff(by_aspect, by_size), 1.0e-5f);
}

UTEST(viamd_camera, orthographic_clip_and_view_are_inverses) {
    /* Deliberately asymmetric bounds - a symmetric box hides a swapped pair. */
    const float l = -3.0f, r = 5.0f, b = -1.5f, t = 4.5f;

    {
        const md_float4x4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t);
        const md_float4x4_t c2v = camera_clip_to_view_matrix_ortho(l, r, b, t);
        EXPECT_LT(max_abs_diff(md_float4x4_mul(v2c, c2v), md_float4x4_ident()), 1.0e-5f);
        EXPECT_LT(max_abs_diff(md_float4x4_mul(c2v, v2c), md_float4x4_ident()), 1.0e-5f);
    }
    {
        const md_float4x4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t, 0.25f, 100.0f);
        const md_float4x4_t c2v = camera_clip_to_view_matrix_ortho(l, r, b, t, 0.25f, 100.0f);
        EXPECT_LT(max_abs_diff(md_float4x4_mul(v2c, c2v), md_float4x4_ident()), 1.0e-5f);
        EXPECT_LT(max_abs_diff(md_float4x4_mul(c2v, v2c), md_float4x4_ident()), 1.0e-5f);
    }
}

UTEST(viamd_camera, orthographic_corners_map_to_the_clip_cube) {
    const float l = -3.0f, r = 5.0f, b = -1.5f, t = 4.5f;
    const md_float4x4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t);

    const md_float4_t lo = md_float4x4_mul_float4(v2c, md_float4_set(l, b, 0, 1));
    const md_float4_t hi = md_float4x4_mul_float4(v2c, md_float4_set(r, t, 0, 1));

    EXPECT_NEAR(-1.0f, lo.x, 1.0e-5f);
    EXPECT_NEAR(-1.0f, lo.y, 1.0e-5f);
    EXPECT_NEAR( 1.0f, hi.x, 1.0e-5f);
    EXPECT_NEAR( 1.0f, hi.y, 1.0e-5f);
}

UTEST(viamd_camera, look_at_and_position_are_inverses) {
    const Camera c = test_camera();

    /* The camera is described by a position plus a distance to what it orbits; the look-at point is
     * derived from those, and the position is derivable back from it. Trackball navigation leans on
     * that round trip every frame. */
    const md_float3_t look_at = camera_get_look_at(c);
    const md_float3_t pos     = camera_position_from_look_at(look_at, c.orientation, c.distance);

    EXPECT_NEAR(c.position.x, pos.x, 1.0e-3f);
    EXPECT_NEAR(c.position.y, pos.y, 1.0e-3f);
    EXPECT_NEAR(c.position.z, pos.z, 1.0e-3f);

    /* And the look-at point really is 'distance' away. */
    EXPECT_NEAR(c.distance, md_float3_length(md_float3_sub(look_at, c.position)), 1.0e-3f);
}

UTEST(viamd_camera, interpolation_reproduces_its_endpoints) {
    const Camera a = test_camera();
    Camera b = test_camera();
    b.position    = md_float3_set(-4.0f, 6.0f, 1.0f);
    b.distance    = 30.0f;
    b.orientation = md_quat_normalize(md_quat_axis_angle(md_float3_set(0, 1, 0), 2.1f));

    md_float3_t in_pos[2] = { a.position,    b.position };
    md_quat_t in_ori[2] = { a.orientation, b.orientation };
    float  in_dst[2] = { a.distance,    b.distance };

    md_float3_t pos; md_quat_t ori; float dst;

    camera_interpolate_look_at(&pos, &ori, &dst, in_pos, in_ori, in_dst, 0.0);
    EXPECT_NEAR(a.position.x, pos.x, 1.0e-4f);
    EXPECT_NEAR(a.distance,   dst,   1.0e-4f);

    camera_interpolate_look_at(&pos, &ori, &dst, in_pos, in_ori, in_dst, 1.0);
    EXPECT_NEAR(b.position.x, pos.x, 1.0e-4f);
    EXPECT_NEAR(b.distance,   dst,   1.0e-4f);

    /* Halfway stays between the two, and the orientation stays a unit quaternion - a normalisation
     * dropped from the blend shows up as a slowly growing scale in the view matrix. */
    camera_interpolate_look_at(&pos, &ori, &dst, in_pos, in_ori, in_dst, 0.5);
    EXPECT_GT(dst, a.distance);
    EXPECT_LT(dst, b.distance);
    const float len = sqrtf(ori.x*ori.x + ori.y*ori.y + ori.z*ori.z + ori.w*ori.w);
    EXPECT_NEAR(1.0f, len, 1.0e-4f);
}

/* Default view. These pin the behaviour rather than the numbers: what is up, what faces the viewer, and
 * that everything fits. Synthetic systems, deterministic (a fixed LCG, not rand()). */

struct DvSys {
    md_float3_t xyz[4096];
    size_t n = 0;
    void add(md_float3_t p) { xyz[n++] = p; }
};

static float dv_rand(uint32_t* s) { *s = *s * 1664525u + 1013904223u; return (float)(*s >> 8) / 16777216.0f; }

static md_float3_t dv_view_dir(const ViewTransform& v) { return md_quat_mul_float3(v.orientation, md_float3_set(0, 0, 1)); }
static md_float3_t dv_view_up (const ViewTransform& v) { return md_quat_mul_float3(v.orientation, md_float3_set(0, 1, 0)); }

static bool dv_all_in_view(const DvSys& s, const ViewTransform& v, float fov_y) {
    const md_float4x4_t V = camera_world_to_view_matrix(v);
    const float t = tanf(fov_y * 0.5f);
    for (size_t i = 0; i < s.n; ++i) {
        const md_float4_t p = md_float4x4_mul_float4(V, md_float4_from_float3(s.xyz[i], 1.0f));
        if (p.z >= 0.0f) return false;
        if (fabsf(p.x) > t * -p.z || fabsf(p.y) > t * -p.z) return false;
    }
    return true;
}

UTEST(viamd_camera, default_view_of_a_filled_box_is_the_world_view) {
    static DvSys s; s.n = 0;
    uint32_t seed = 1;
    for (int i = 0; i < 4000; ++i) s.add(md_float3_set(60 * dv_rand(&seed), 60 * dv_rand(&seed), 60 * dv_rand(&seed)));
    md_float3x3_t A = {};
    A.elem[0][0] = A.elem[1][1] = A.elem[2][2] = 60.0f;
    const float fov = 0.785f;
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, &A, fov);
    const md_float3_t b = dv_view_dir(v);   /* towards the camera */
    const md_float3_t u = dv_view_up(v);
    EXPECT_GT(u.z, 0.9f);              /* Z up */
    EXPECT_GT(b.x, 0.8f);              /* X towards the viewer... */
    EXPECT_GT(b.y, 0.0f);              /* ...and to the left: the camera is swung towards +Y */
    EXPECT_GT(b.z, 0.0f);              /* seen a little from above */
    EXPECT_TRUE(dv_all_in_view(s, v, fov));
}

UTEST(viamd_camera, default_view_of_a_membrane_without_a_cell_keeps_its_normal_up) {
    /* A bilayer-sized slab whose normal is 10 degrees off Z: snapped to Z, seen from the side */
    static DvSys s; s.n = 0;
    uint32_t seed = 2;
    const float c = cosf(0.1745f), sn = sinf(0.1745f);
    for (int i = 0; i < 4000; ++i) {
        const float px = 80 * dv_rand(&seed), py = 80 * dv_rand(&seed) - 40, pz = 40 * dv_rand(&seed) - 20;
        s.add(md_float3_set(px, c * py - sn * pz, sn * py + c * pz));
    }
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_GT(dv_view_up(v).z, 0.9f);
    EXPECT_LT(fabsf(dv_view_dir(v).z), 0.5f);
}

UTEST(viamd_camera, default_view_of_a_planar_molecule_is_face_on) {
    /* Benzene, in an arbitrary orientation */
    static DvSys s; s.n = 0;
    const md_float3_t n  = md_float3_normalize(md_float3_set(0.3f, -0.5f, 0.8f));
    const md_float3_t e1 = md_float3_normalize(md_float3_cross(n, md_float3_set(1, 0, 0)));
    const md_float3_t e2 = md_float3_cross(n, e1);
    for (int k = 0; k < 6; ++k) {
        const float a = k * 1.0471976f;
        s.add(md_float3_add(md_float3_mul1(e1, 1.39f * cosf(a)), md_float3_mul1(e2, 1.39f * sinf(a))));
        s.add(md_float3_add(md_float3_mul1(e1, 2.47f * cosf(a)), md_float3_mul1(e2, 2.47f * sinf(a))));
    }
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_GT(fabsf(md_float3_dot(dv_view_dir(v), n)), 0.99f);
}

UTEST(viamd_camera, default_view_does_not_look_down_a_bond) {
    /* Staggered ethane: looking down C-C hides one carbon behind the other */
    static DvSys s; s.n = 0;
    const md_float3_t ax = md_float3_normalize(md_float3_set(0.6f, 0.2f, -0.7f));
    const md_float3_t e1 = md_float3_normalize(md_float3_cross(ax, md_float3_set(0, 1, 0)));
    const md_float3_t e2 = md_float3_cross(ax, e1);
    for (int c = 0; c < 2; ++c) {
        const float h = c ? 0.765f : -0.765f;
        s.add(md_float3_mul1(ax, h));
        for (int k = 0; k < 3; ++k) {
            const float a = k * 2.0943951f + (c ? 1.0471976f : 0.0f);
            s.add(md_float3_add(md_float3_mul1(ax, h + (c ? 0.36f : -0.36f)), md_float3_add(md_float3_mul1(e1, 1.03f * cosf(a)), md_float3_mul1(e2, 1.03f * sinf(a)))));
        }
    }
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_LT(fabsf(md_float3_dot(dv_view_dir(v), ax)), 0.8f);
}

UTEST(viamd_camera, default_view_is_independent_of_the_structure_sign) {
    /* The same elongated structure mirrored through its center must not flip the view upside down */
    static DvSys s; s.n = 0;
    static DvSys m; m.n = 0;
    uint32_t seed = 3;
    for (int i = 0; i < 2000; ++i) {
        const md_float3_t p = md_float3_set(40 * (dv_rand(&seed) - 0.5f), 20 * (dv_rand(&seed) - 0.5f), 10 * (dv_rand(&seed) - 0.5f));
        s.add(p);
        m.add(md_float3_mul1(p, -1.0f));
    }
    const ViewTransform a = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    const ViewTransform b = camera_compute_default_view(m.xyz, m.n, NULL, m.n, NULL, 0.785f);
    EXPECT_GT(md_float3_dot(dv_view_up(a), dv_view_up(b)), 0.99f);
    EXPECT_GT(md_float3_dot(dv_view_dir(a), dv_view_dir(b)), 0.99f);
}

/* Framing in a region of the view, the part that windows leave uncovered. The camera looks at what it frames
 * and the projection is shifted so that it is seen at the region's center, so these project through that. */

static md_float2_t dv_project(const ViewTransform& v, const md_float4x4_t& P, md_float3_t x) {
    const md_float4_t p = md_float4x4_mul_float4(P * camera_world_to_view_matrix(v), md_float4_from_float3(x, 1.0f));
    return md_float2_set(p.x / p.w, p.y / p.w);
}

UTEST(viamd_camera, shifted_projection_shows_the_look_at_at_the_region_center) {
    const Camera c = test_camera();
    const float aspect = 16.0f / 9.0f;
    const ViewRegion region = {.min = {0.1f, -0.9f}, .max = {0.9f, 0.3f}, .aspect_ratio = aspect};
    const md_float2_t center = view_region_center(region);
    const md_float4x4_t shift = md_float4x4_translate(center.x, center.y, 0.0f);
    const md_float3_t look_at = camera_get_look_at(c);

    const md_float2_t persp = dv_project(c, shift * camera_view_to_clip_matrix_persp(c, aspect), look_at);
    EXPECT_NEAR( 0.5f, persp.x, 1.0e-4f);
    EXPECT_NEAR(-0.3f, persp.y, 1.0e-4f);

    const float h = c.distance * tanf(c.fov_y * 0.5f);
    const float w = aspect * h;
    const md_float2_t ortho = dv_project(c, shift * camera_view_to_clip_matrix_ortho(-w, w, -h, h, c.near_plane, c.far_plane), look_at);
    EXPECT_NEAR( 0.5f, ortho.x, 1.0e-4f);
    EXPECT_NEAR(-0.3f, ortho.y, 1.0e-4f);

    /* And back: the shifted inverse sends the region center to the view axis */
    const md_float4x4_t inv = camera_clip_to_view_matrix_persp(c, aspect) * md_float4x4_translate(-center.x, -center.y, 0.0f);
    const md_float4_t v = md_float4x4_mul_float4(inv, md_float4_set(center.x, center.y, 0.5f, 1.0f));
    EXPECT_NEAR(0.0f, v.x / v.w, 1.0e-4f);
    EXPECT_NEAR(0.0f, v.y / v.w, 1.0e-4f);
}

UTEST(viamd_camera, default_view_in_a_region_fits_in_it) {
    static DvSys s; s.n = 0;
    uint32_t seed = 4;
    for (int i = 0; i < 2000; ++i) {
        s.add(md_float3_set(60 * (dv_rand(&seed) - 0.5f), 25 * (dv_rand(&seed) - 0.5f), 25 * (dv_rand(&seed) - 0.5f)));
    }
    const float fov = 0.785f;
    const float aspect = 16.0f / 9.0f;
    const md_float4x4_t P = md_float4x4_persp(fov, aspect, 1.0f, 10000.0f);

    /* Windows to the left and at the bottom; then windows over most of the view, which leave a tall column */
    const ViewRegion regions[2] = {
        {.min = {-0.3f, -0.5f}, .max = {1.0f, 1.0f}, .aspect_ratio = aspect},
        {.min = { 0.5f, -1.0f}, .max = {1.0f, 1.0f}, .aspect_ratio = aspect},
    };
    for (const ViewRegion& region : regions) {
        const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, fov, region);
        const md_float2_t center = view_region_center(region);
        const md_float4x4_t shifted = md_float4x4_translate(center.x, center.y, 0.0f) * P;
        md_float2_t lo = md_float2_set( FLT_MAX,  FLT_MAX);
        md_float2_t hi = md_float2_set(-FLT_MAX, -FLT_MAX);
        for (size_t i = 0; i < s.n; ++i) {
            const md_float2_t p = dv_project(v, shifted, s.xyz[i]);
            lo = md_float2_set(fminf(lo.x, p.x), fminf(lo.y, p.y));
            hi = md_float2_set(fmaxf(hi.x, p.x), fmaxf(hi.y, p.y));
        }
        EXPECT_GE(lo.x, region.min.x);
        EXPECT_GE(lo.y, region.min.y);
        EXPECT_LE(hi.x, region.max.x);
        EXPECT_LE(hi.y, region.max.y);
        /* In the middle of it, and not lost in it */
        const md_float2_t mid = (lo + hi) * 0.5f;
        EXPECT_LT(fabsf(mid.x - center.x), 0.1f);
        EXPECT_LT(fabsf(mid.y - center.y), 0.1f);
        EXPECT_GT(fmaxf((hi.x - lo.x) / (region.max.x - region.min.x), (hi.y - lo.y) / (region.max.y - region.min.y)), 0.5f);
    }
}

UTEST(viamd_camera, trackball_rotates_about_the_shifted_center) {
    /* A vertical drag straight through where the look-at is seen tilts the view about its right axis, and does
     * not roll it - which a drag off the trackball's center would */
    const Camera c = test_camera();
    TrackballControllerInput input = {};
    input.rotate_button    = true;
    input.screen_size      = md_float2_set(200, 100);
    input.ndc_center       = md_float2_set(0.5f, 0.0f);   /* pixel x 150 */
    input.mouse_coord_prev = md_float2_set(150, 40);
    input.mouse_coord_curr = md_float2_set(150, 60);

    ViewTransform v = c;
    EXPECT_TRUE(camera_controller_trackball(&v, input, {}, TrackballFlags_RotateEnabled | TrackballFlags_RotateReturnsTrue));
    const md_float3_t right_before = md_quat_mul_float3(c.orientation, md_float3_set(1, 0, 0));
    const md_float3_t right_after  = md_quat_mul_float3(v.orientation, md_float3_set(1, 0, 0));
    EXPECT_NEAR(1.0f, md_float3_dot(right_before, right_after), 1.0e-5f);
    EXPECT_LT(md_float3_dot(dv_view_dir(c), dv_view_dir(v)), 0.9999f);   /* it did turn */
    /* About the look-at, which stays put */
    const md_float3_t a = camera_get_look_at(c), b = camera_get_look_at(v);
    EXPECT_NEAR(a.x, b.x, 1.0e-3f);
    EXPECT_NEAR(a.y, b.y, 1.0e-3f);
    EXPECT_NEAR(a.z, b.z, 1.0e-3f);
}
