#include "utest.h"

#include <gfx/camera_utils.h>
#include <gfx/camera.h>

#include <core/md_vec_math.h>

#include <math.h>

/* camera_utils is where a sign or a transposition goes unnoticed for a long time: a projection that
 * is subtly wrong still renders something, and the picture only looks off at the extremes. So these
 * tests check the properties rather than the entries - a matrix and its stated inverse must compose
 * to the identity, and a point sent through a projection must come back.
 *
 * Every camera here is off-axis and off-origin. With an identity orientation at the origin most of
 * these compositions hold for the wrong reasons. */

static float max_abs_diff(mat4_t A, mat4_t B) {
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
    c.orientation = quat_normalize(quat_axis_angle(vec3_normalize(vec3_set(0.3f, 1.0f, -0.2f)), 0.7f));
    c.position    = vec3_set(2.5f, -1.25f, 8.0f);
    c.distance    = 12.0f;
    c.near_plane  = 0.5f;
    c.far_plane   = 250.0f;
    c.fov_y       = 1.0471975512f;   /* 60 degrees */
    return c;
}

UTEST(viamd_camera, view_and_world_transforms_are_inverses) {
    const Camera c = test_camera();

    const mat4_t w2v = camera_world_to_view_matrix(c);
    const mat4_t v2w = camera_view_to_world_matrix(c);

    EXPECT_LT(max_abs_diff(mat4_mul(w2v, v2w), mat4_ident()), 1.0e-4f);
    EXPECT_LT(max_abs_diff(mat4_mul(v2w, w2v), mat4_ident()), 1.0e-4f);

    /* The camera sits at the origin of its own view space. */
    const vec4_t eye = mat4_mul_vec4(w2v, vec4_set(c.position.x, c.position.y, c.position.z, 1.0f));
    EXPECT_NEAR(0.0f, eye.x, 1.0e-4f);
    EXPECT_NEAR(0.0f, eye.y, 1.0e-4f);
    EXPECT_NEAR(0.0f, eye.z, 1.0e-4f);
}

UTEST(viamd_camera, perspective_clip_and_view_are_inverses) {
    const Camera c = test_camera();
    const float aspect = 16.0f / 9.0f;

    const mat4_t v2c = camera_view_to_clip_matrix_persp(c, aspect);
    const mat4_t c2v = camera_clip_to_view_matrix_persp(c, aspect);

    EXPECT_LT(max_abs_diff(mat4_mul(v2c, c2v), mat4_ident()), 1.0e-4f);
    EXPECT_LT(max_abs_diff(mat4_mul(c2v, v2c), mat4_ident()), 1.0e-4f);
}

UTEST(viamd_camera, perspective_maps_the_clip_planes_to_the_depth_range) {
    const Camera c = test_camera();
    const mat4_t v2c = camera_view_to_clip_matrix_persp(c, 1.0f);

    /* Looking down -z in view space. The near plane lands at one end of the depth range and the far
     * plane at the other. The span is what identifies the convention: 2 is OpenGL's [-1,1] clip
     * space, 1 would be the [0,1] that D3D and Vulkan use. mdlib's mat4_persp builds the GL form and
     * viamd's shaders read it back that way, so the two have to keep agreeing. */
    const vec4_t n = mat4_mul_vec4(v2c, vec4_set(0, 0, -c.near_plane, 1.0f));
    const vec4_t f = mat4_mul_vec4(v2c, vec4_set(0, 0, -c.far_plane,  1.0f));

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
    const mat4_t by_aspect = camera_view_to_clip_matrix_persp(c, 1920.0f / 1080.0f);
    const mat4_t by_size   = camera_view_to_clip_matrix_persp(c, 1920, 1080, 0.0f, 0.0f);

    EXPECT_LT(max_abs_diff(by_aspect, by_size), 1.0e-5f);
}

UTEST(viamd_camera, orthographic_clip_and_view_are_inverses) {
    /* Deliberately asymmetric bounds - a symmetric box hides a swapped pair. */
    const float l = -3.0f, r = 5.0f, b = -1.5f, t = 4.5f;

    {
        const mat4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t);
        const mat4_t c2v = camera_clip_to_view_matrix_ortho(l, r, b, t);
        EXPECT_LT(max_abs_diff(mat4_mul(v2c, c2v), mat4_ident()), 1.0e-5f);
        EXPECT_LT(max_abs_diff(mat4_mul(c2v, v2c), mat4_ident()), 1.0e-5f);
    }
    {
        const mat4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t, 0.25f, 100.0f);
        const mat4_t c2v = camera_clip_to_view_matrix_ortho(l, r, b, t, 0.25f, 100.0f);
        EXPECT_LT(max_abs_diff(mat4_mul(v2c, c2v), mat4_ident()), 1.0e-5f);
        EXPECT_LT(max_abs_diff(mat4_mul(c2v, v2c), mat4_ident()), 1.0e-5f);
    }
}

UTEST(viamd_camera, orthographic_corners_map_to_the_clip_cube) {
    const float l = -3.0f, r = 5.0f, b = -1.5f, t = 4.5f;
    const mat4_t v2c = camera_view_to_clip_matrix_ortho(l, r, b, t);

    const vec4_t lo = mat4_mul_vec4(v2c, vec4_set(l, b, 0, 1));
    const vec4_t hi = mat4_mul_vec4(v2c, vec4_set(r, t, 0, 1));

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
    const vec3_t look_at = camera_get_look_at(c);
    const vec3_t pos     = camera_position_from_look_at(look_at, c.orientation, c.distance);

    EXPECT_NEAR(c.position.x, pos.x, 1.0e-3f);
    EXPECT_NEAR(c.position.y, pos.y, 1.0e-3f);
    EXPECT_NEAR(c.position.z, pos.z, 1.0e-3f);

    /* And the look-at point really is 'distance' away. */
    EXPECT_NEAR(c.distance, vec3_length(vec3_sub(look_at, c.position)), 1.0e-3f);
}

UTEST(viamd_camera, interpolation_reproduces_its_endpoints) {
    const Camera a = test_camera();
    Camera b = test_camera();
    b.position    = vec3_set(-4.0f, 6.0f, 1.0f);
    b.distance    = 30.0f;
    b.orientation = quat_normalize(quat_axis_angle(vec3_set(0, 1, 0), 2.1f));

    vec3_t in_pos[2] = { a.position,    b.position };
    quat_t in_ori[2] = { a.orientation, b.orientation };
    float  in_dst[2] = { a.distance,    b.distance };

    vec3_t pos; quat_t ori; float dst;

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
    float x[4096], y[4096], z[4096];
    size_t n = 0;
    void add(vec3_t p) { x[n] = p.x; y[n] = p.y; z[n] = p.z; ++n; }
};

static float dv_rand(uint32_t* s) { *s = *s * 1664525u + 1013904223u; return (float)(*s >> 8) / 16777216.0f; }

static vec3_t dv_view_dir(const ViewTransform& v) { return quat_mul_vec3(v.orientation, vec3_set(0, 0, 1)); }
static vec3_t dv_view_up (const ViewTransform& v) { return quat_mul_vec3(v.orientation, vec3_set(0, 1, 0)); }

static bool dv_all_in_view(const DvSys& s, const ViewTransform& v, float fov_y) {
    const mat4_t V = camera_world_to_view_matrix(v);
    const float t = tanf(fov_y * 0.5f);
    for (size_t i = 0; i < s.n; ++i) {
        const vec4_t p = mat4_mul_vec4(V, vec4_set(s.x[i], s.y[i], s.z[i], 1.0f));
        if (p.z >= 0.0f) return false;
        if (fabsf(p.x) > t * -p.z || fabsf(p.y) > t * -p.z) return false;
    }
    return true;
}

UTEST(viamd_camera, default_view_of_a_filled_box_is_the_world_view) {
    static DvSys s; s.n = 0;
    uint32_t seed = 1;
    for (int i = 0; i < 4000; ++i) s.add(vec3_set(60 * dv_rand(&seed), 60 * dv_rand(&seed), 60 * dv_rand(&seed)));
    mat3_t A = {};
    A.elem[0][0] = A.elem[1][1] = A.elem[2][2] = 60.0f;
    const float fov = 0.785f;
    const ViewTransform v = camera_compute_default_view(s.x, s.y, s.z, s.n, NULL, s.n, &A, fov);
    const vec3_t b = dv_view_dir(v);   /* towards the camera */
    const vec3_t u = dv_view_up(v);
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
        s.add(vec3_set(px, c * py - sn * pz, sn * py + c * pz));
    }
    const ViewTransform v = camera_compute_default_view(s.x, s.y, s.z, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_GT(dv_view_up(v).z, 0.9f);
    EXPECT_LT(fabsf(dv_view_dir(v).z), 0.5f);
}

UTEST(viamd_camera, default_view_of_a_planar_molecule_is_face_on) {
    /* Benzene, in an arbitrary orientation */
    static DvSys s; s.n = 0;
    const vec3_t n  = vec3_normalize(vec3_set(0.3f, -0.5f, 0.8f));
    const vec3_t e1 = vec3_normalize(vec3_cross(n, vec3_set(1, 0, 0)));
    const vec3_t e2 = vec3_cross(n, e1);
    for (int k = 0; k < 6; ++k) {
        const float a = k * 1.0471976f;
        s.add(vec3_add(vec3_mul1(e1, 1.39f * cosf(a)), vec3_mul1(e2, 1.39f * sinf(a))));
        s.add(vec3_add(vec3_mul1(e1, 2.47f * cosf(a)), vec3_mul1(e2, 2.47f * sinf(a))));
    }
    const ViewTransform v = camera_compute_default_view(s.x, s.y, s.z, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_GT(fabsf(vec3_dot(dv_view_dir(v), n)), 0.99f);
}

UTEST(viamd_camera, default_view_does_not_look_down_a_bond) {
    /* Staggered ethane: looking down C-C hides one carbon behind the other */
    static DvSys s; s.n = 0;
    const vec3_t ax = vec3_normalize(vec3_set(0.6f, 0.2f, -0.7f));
    const vec3_t e1 = vec3_normalize(vec3_cross(ax, vec3_set(0, 1, 0)));
    const vec3_t e2 = vec3_cross(ax, e1);
    for (int c = 0; c < 2; ++c) {
        const float h = c ? 0.765f : -0.765f;
        s.add(vec3_mul1(ax, h));
        for (int k = 0; k < 3; ++k) {
            const float a = k * 2.0943951f + (c ? 1.0471976f : 0.0f);
            s.add(vec3_add(vec3_mul1(ax, h + (c ? 0.36f : -0.36f)), vec3_add(vec3_mul1(e1, 1.03f * cosf(a)), vec3_mul1(e2, 1.03f * sinf(a)))));
        }
    }
    const ViewTransform v = camera_compute_default_view(s.x, s.y, s.z, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_LT(fabsf(vec3_dot(dv_view_dir(v), ax)), 0.8f);
}

UTEST(viamd_camera, default_view_is_independent_of_the_structure_sign) {
    /* The same elongated structure mirrored through its center must not flip the view upside down */
    static DvSys s; s.n = 0;
    static DvSys m; m.n = 0;
    uint32_t seed = 3;
    for (int i = 0; i < 2000; ++i) {
        const vec3_t p = vec3_set(40 * (dv_rand(&seed) - 0.5f), 20 * (dv_rand(&seed) - 0.5f), 10 * (dv_rand(&seed) - 0.5f));
        s.add(p);
        m.add(vec3_mul1(p, -1.0f));
    }
    const ViewTransform a = camera_compute_default_view(s.x, s.y, s.z, s.n, NULL, s.n, NULL, 0.785f);
    const ViewTransform b = camera_compute_default_view(m.x, m.y, m.z, m.n, NULL, m.n, NULL, 0.785f);
    EXPECT_GT(vec3_dot(dv_view_up(a), dv_view_up(b)), 0.99f);
    EXPECT_GT(vec3_dot(dv_view_dir(a), dv_view_dir(b)), 0.99f);
}
