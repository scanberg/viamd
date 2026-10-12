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

/* Keyframe path. Uneven spacing and an off-axis, off-origin pose per key, so a time that is not used
 * as a time (or a rotation that is not taken the short way) shows up. */

static void kf_set(CameraKeyframe& k, double time, vec3_t axis, float angle, vec3_t look_at, float dist, float fov) {
    k.time = time;
    k.transform.orientation = quat_normalize(quat_axis_angle(vec3_normalize(axis), angle));
    k.transform.distance = dist;
    k.transform.position = camera_position_from_look_at(look_at, k.transform.orientation, dist);
    k.fov_y = fov;
}

static void kf_make3(CameraKeyframe* k) {
    kf_set(k[0], 0.0, vec3_set(0.3f, 1.0f, 0.1f), 0.2f, vec3_set(0, 0, 0),     10.0f, 0.8f);
    kf_set(k[1], 1.0, vec3_set(0.3f, 1.0f, 0.1f), 1.0f, vec3_set(10, 5, -3),   20.0f, 0.6f);
    kf_set(k[2], 3.0, vec3_set(0.3f, 1.0f, 0.1f), 1.8f, vec3_set(30, -5, 12),  15.0f, 1.0f);
}

UTEST(viamd_camera, depth_of_point_is_measured_along_the_view_direction) {
    const Camera c = test_camera();
    /* The look-at point is `distance` in front of the eye, however the camera is turned */
    EXPECT_NEAR(c.distance, camera_depth_of_point(c, camera_get_look_at(c)), 1.0e-3f);

    /* A sideways offset does not change the depth, a move along the view direction does */
    const vec3_t fwd   = quat_mul_vec3(c.orientation, vec3_set(0, 0, -1));
    const vec3_t right = quat_mul_vec3(c.orientation, vec3_set(1, 0, 0));
    const vec3_t p = c.position + fwd * 7.0f + right * 3.0f;
    EXPECT_NEAR(7.0f, camera_depth_of_point(c, p), 1.0e-3f);
    EXPECT_LT(camera_depth_of_point(c, c.position - fwd * 2.0f), 0.0f);
}

UTEST(viamd_camera, keyframes_hit_their_poses_and_hold_outside) {
    CameraKeyframe k[3];
    kf_make3(k);
    for (int i = 0; i < 3; ++i) {
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 3, k[i].time);
        EXPECT_NEAR(k[i].transform.position.x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.z, v.position.z, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.distance, v.distance, 1.0e-4f);
        EXPECT_NEAR(k[i].fov_y, fov, 1.0e-5f);
        EXPECT_NEAR(1.0f, fabsf(quat_dot(k[i].transform.orientation, v.orientation)), 1.0e-5f);
    }
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 3, -5.0);
    EXPECT_NEAR(k[0].transform.distance, v.distance, 1.0e-4f);
    camera_keyframes_evaluate(&v, &fov, k, 3, 99.0);
    EXPECT_NEAR(k[2].transform.distance, v.distance, 1.0e-4f);
    EXPECT_NEAR(k[2].fov_y, fov, 1.0e-5f);
    camera_keyframes_evaluate(&v, &fov, k, 1, 0.5);
    EXPECT_NEAR(k[0].transform.distance, v.distance, 1.0e-4f);
}

UTEST(viamd_camera, a_following_key_keeps_its_look_at_relative_to_the_target) {
    CameraKeyframe k[3];
    kf_make3(k);
    const vec3_t c0 = vec3_set(1, 2, 3);
    const vec3_t c1 = vec3_set(-4, 0, 5);
    k[0].follow = true; k[0].follow_center = c0;
    k[1].follow = true; k[1].follow_center = c1;
    k[2].follow = true; k[2].follow_center = vec3_set(0, 0, 0);

    /* With the target where a key saw it, the key is hit exactly. */
    for (int i = 0; i < 3; ++i) {
        const vec3_t now = k[i].follow_center;
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 3, k[i].time, false, &now);
        EXPECT_NEAR(k[i].transform.position.x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.z, v.position.z, 1.0e-3f);
    }

    /* A target that has moved takes the look-at point with it, by the same amount at every key. */
    const vec3_t shift = vec3_set(7, -3, 2);
    for (int i = 0; i < 3; ++i) {
        const vec3_t now = k[i].follow_center + shift;
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 3, k[i].time, false, &now);
        EXPECT_NEAR(k[i].transform.position.x + shift.x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.y + shift.y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.z + shift.z, v.position.z, 1.0e-3f);
    }

    /* Without a target the keys are fixed points as before. */
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 3, k[1].time, false, nullptr);
    EXPECT_NEAR(k[1].transform.position.x, v.position.x, 1.0e-3f);
}

UTEST(viamd_camera, a_fixed_key_ignores_the_target) {
    CameraKeyframe k[2];
    kf_set(k[0], 0.0, vec3_set(0.3f, 1.0f, 0.1f), 0.2f, vec3_set(0, 0, 0),  10.0f, 0.8f);
    kf_set(k[1], 2.0, vec3_set(0.3f, 1.0f, 0.1f), 0.2f, vec3_set(10, 0, 0), 10.0f, 0.8f);
    k[1].follow = true;
    k[1].follow_center = vec3_set(1, 1, 1);

    const vec3_t now = vec3_set(5, 5, 5);
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 2, 0.0, false, &now);
    EXPECT_NEAR(k[0].transform.position.x, v.position.x, 1.0e-3f);
    EXPECT_NEAR(k[0].transform.position.y, v.position.y, 1.0e-3f);

    /* The follow weight grows from 0 to 1 over the segment, so halfway the target counts for half. */
    camera_keyframes_evaluate(&v, &fov, k, 2, 1.0, false, &now);
    const vec3_t plain = camera_get_look_at(k[0].transform) * 0.5f + camera_get_look_at(k[1].transform) * 0.5f;
    const vec3_t look = camera_get_look_at(v);
    EXPECT_NEAR(plain.x + 0.5f * (now.x - 0.0f) - 0.5f * (k[1].follow_center.x - 0.0f), look.x, 5.0e-2f);
}

UTEST(viamd_camera, keyframes_pass_through_a_key_without_stopping) {
    CameraKeyframe k[3];
    kf_make3(k);
    const double eps = 1.0e-3;
    ViewTransform a, b, c; float fov;
    camera_keyframes_evaluate(&a, &fov, k, 3, 1.0 - eps);
    camera_keyframes_evaluate(&b, &fov, k, 3, 1.0);
    camera_keyframes_evaluate(&c, &fov, k, 3, 1.0 + eps);

    /* Same velocity on both sides (no corner) ... */
    const vec3_t v0 = (b.position - a.position) * (float)(1.0 / eps);
    const vec3_t v1 = (c.position - b.position) * (float)(1.0 / eps);
    EXPECT_LT(vec3_length(v1 - v0), 0.05f * vec3_length(v0));
    /* ... and it is moving: the average speed over the path is ~14 per second; a per-segment ease would be 0 here. */
    EXPECT_GT(vec3_length(v0), 5.0f);

    /* Rotation speed through the key is likewise non-zero and continuous */
    const double re = 0.05;   /* acos of a float dot is too coarse for a smaller step */
    camera_keyframes_evaluate(&a, &fov, k, 3, 1.0 - re);
    camera_keyframes_evaluate(&c, &fov, k, 3, 1.0 + re);
    const float ang0 = 2.0f * acosf(fminf(1.0f, fabsf(quat_dot(a.orientation, b.orientation)))) / (float)re;
    const float ang1 = 2.0f * acosf(fminf(1.0f, fabsf(quat_dot(b.orientation, c.orientation)))) / (float)re;
    EXPECT_GT(ang0, 0.1f);
    EXPECT_NEAR(ang0, ang1, 0.15f * ang0);
}

UTEST(viamd_camera, keyframes_stay_valid_and_do_not_overshoot) {
    CameraKeyframe k[3];
    kf_make3(k);
    /* The middle key is the extreme of both distance and fov: the path must not exceed it */
    k[1].transform.distance = 30.0f;
    k[1].fov_y = 1.2f;
    /* A quaternion of opposite sign is the same rotation; the path must not spin the long way around */
    k[1].transform.orientation = quat_t{-k[1].transform.orientation.x, -k[1].transform.orientation.y, -k[1].transform.orientation.z, -k[1].transform.orientation.w};
    for (double t = 0.0; t <= 3.0; t += 0.01) {
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 3, t);
        EXPECT_LE(v.distance, 30.0f + 1.0e-3f);
        EXPECT_GE(v.distance, 10.0f - 1.0e-3f);
        EXPECT_LE(fov, 1.2f + 1.0e-4f);
        const quat_t q = v.orientation;
        EXPECT_NEAR(1.0f, sqrtf(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w), 1.0e-4f);
        /* Rotation about one axis, 0.2 -> 1.0 -> 1.8 rad: the angle stays within that range */
        const float angle = 2.0f * acosf(fminf(1.0f, fabsf(q.w)));
        EXPECT_GT(angle, 0.2f - 0.05f);
        EXPECT_LT(angle, 1.8f + 0.05f);
    }
}

/* Extra turns. Two keys with the same pose, so whatever the camera does in between is the spin. */

static void kf_orbit(CameraKeyframe* k, int turns, SpinAxis axis, bool constant_speed) {
    kf_set(k[0], 0.0, vec3_set(0.3f, 1.0f, 0.1f), 0.4f, vec3_set(5, -2, 3), 12.0f, 0.8f);
    k[1] = k[0];
    k[1].time = 4.0;
    k[1].spin_turns = turns;
    k[1].spin_axis = axis;
    k[1].spin_constant_speed = constant_speed;
}

static vec3_t kf_relative_eye(const CameraKeyframe& k, const ViewTransform& v) {
    return v.position - camera_get_look_at(k.transform);
}

UTEST(viamd_camera, a_spin_orbits_the_look_at_point_and_arrives_where_it_started) {
    CameraKeyframe k[2];
    kf_orbit(k, 1, SpinAxis::WorldY, true);
    const vec3_t r0 = kf_relative_eye(k[0], k[0].transform);

    ViewTransform v; float fov;
    for (double t = 0.0; t <= 4.0; t += 0.25) {
        camera_keyframes_evaluate(&v, &fov, k, 2, t);
        EXPECT_NEAR(vec3_length(r0), vec3_length(kf_relative_eye(k[0], v)), 1.0e-3f);
        EXPECT_NEAR(1.0f, sqrtf(v.orientation.x*v.orientation.x + v.orientation.y*v.orientation.y + v.orientation.z*v.orientation.z + v.orientation.w*v.orientation.w), 1.0e-4f);
    }

    /* Both ends are the key's pose, whatever the sign of the quaternion */
    for (int e = 0; e < 2; ++e) {
        camera_keyframes_evaluate(&v, &fov, k, 2, k[e].time);
        EXPECT_NEAR(k[0].transform.position.x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[0].transform.position.y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[0].transform.position.z, v.position.z, 1.0e-3f);
        EXPECT_NEAR(1.0f, fabsf(quat_dot(k[0].transform.orientation, v.orientation)), 1.0e-5f);
    }

    /* A quarter of the way at constant speed is a quarter turn: counter-clockwise seen from +Y is
     * (x, y, z) -> (z, y, -x). Halfway it is on the other side. */
    camera_keyframes_evaluate(&v, &fov, k, 2, 1.0);
    vec3_t r = kf_relative_eye(k[0], v);
    EXPECT_NEAR( r0.z, r.x, 2.0e-3f);
    EXPECT_NEAR( r0.y, r.y, 2.0e-3f);
    EXPECT_NEAR(-r0.x, r.z, 2.0e-3f);

    camera_keyframes_evaluate(&v, &fov, k, 2, 2.0);
    r = kf_relative_eye(k[0], v);
    EXPECT_NEAR(-r0.x, r.x, 2.0e-3f);
    EXPECT_NEAR( r0.y, r.y, 2.0e-3f);
    EXPECT_NEAR(-r0.z, r.z, 2.0e-3f);
}

UTEST(viamd_camera, a_negative_spin_goes_the_other_way) {
    CameraKeyframe k[2];
    kf_orbit(k, -1, SpinAxis::WorldY, true);
    const vec3_t r0 = kf_relative_eye(k[0], k[0].transform);

    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 2, 1.0);
    const vec3_t r = kf_relative_eye(k[0], v);
    EXPECT_NEAR(-r0.z, r.x, 2.0e-3f);
    EXPECT_NEAR( r0.x, r.z, 2.0e-3f);
}

UTEST(viamd_camera, a_spin_eases_in_and_out_unless_it_is_constant) {
    CameraKeyframe k[2];
    kf_orbit(k, 1, SpinAxis::WorldY, false);
    const vec3_t r0 = kf_relative_eye(k[0], k[0].transform);

    /* A quarter of the way the eased turn has covered smoothstep(0.25) of it */
    const double s = 0.25 * 0.25 * (3.0 - 2.0 * 0.25);
    const float a = (float)(6.283185307179586 * s);
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 2, 1.0);
    const vec3_t r = kf_relative_eye(k[0], v);
    EXPECT_NEAR(r0.x * cosf(a) + r0.z * sinf(a), r.x, 2.0e-3f);
    EXPECT_NEAR(-r0.x * sinf(a) + r0.z * cosf(a), r.z, 2.0e-3f);

    /* and it starts slowly */
    ViewTransform a0, a1; 
    camera_keyframes_evaluate(&a0, &fov, k, 2, 0.0);
    camera_keyframes_evaluate(&a1, &fov, k, 2, 0.01);
    EXPECT_LT(vec3_length(a1.position - a0.position), 0.05f);
}

UTEST(viamd_camera, a_spin_about_the_view_up_keeps_the_camera_level) {
    CameraKeyframe k[2];
    kf_orbit(k, 2, SpinAxis::ViewUp, true);
    const vec3_t up0 = k[0].transform.orientation * vec3_set(0, 1, 0);
    const vec3_t r0 = kf_relative_eye(k[0], k[0].transform);

    ViewTransform v; float fov;
    for (double t = 0.0; t <= 4.0; t += 0.2) {
        camera_keyframes_evaluate(&v, &fov, k, 2, t);
        const vec3_t up = v.orientation * vec3_set(0, 1, 0);
        EXPECT_NEAR(up0.x, up.x, 1.0e-3f);
        EXPECT_NEAR(up0.y, up.y, 1.0e-3f);
        EXPECT_NEAR(up0.z, up.z, 1.0e-3f);
        /* and it stays in the plane it started in */
        EXPECT_NEAR(vec3_dot(r0, up0), vec3_dot(kf_relative_eye(k[0], v), up0), 2.0e-3f);
    }
}

UTEST(viamd_camera, a_spin_is_not_applied_to_the_first_key) {
    CameraKeyframe k[2];
    kf_orbit(k, 1, SpinAxis::WorldY, true);
    k[0].spin_turns = 3;
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 2, 0.0);
    EXPECT_NEAR(k[0].transform.position.x, v.position.x, 1.0e-3f);
    EXPECT_NEAR(k[0].transform.position.z, v.position.z, 1.0e-3f);
}

/* Trajectory frame from keyframes. */

static void kf_frame(CameraKeyframe& k, double time, double frame) {
    k = CameraKeyframe{};
    k.time = time;
    k.use_frame = true;
    k.frame = frame;
}

UTEST(viamd_camera, no_frame_keys_gives_no_frame_curve) {
    CameraKeyframe k[2];
    k[0] = CameraKeyframe{};
    k[1] = CameraKeyframe{};
    k[1].time = 1.0;
    double f = -1.0;
    EXPECT_FALSE(camera_keyframes_evaluate_frame(&f, k, 2, 0.5));
    EXPECT_FALSE(camera_keyframes_evaluate_frame(&f, k, 0, 0.5));
}

UTEST(viamd_camera, two_frame_keys_play_at_constant_speed_and_hold_outside) {
    CameraKeyframe k[2];
    kf_frame(k[0], 1.0, 0.0);
    kf_frame(k[1], 11.0, 100.0);
    double f = 0.0;
    ASSERT_TRUE(camera_keyframes_evaluate_frame(&f, k, 2, 6.0));
    EXPECT_NEAR(50.0, f, 1.0e-9);
    camera_keyframes_evaluate_frame(&f, k, 2, 3.0);
    EXPECT_NEAR(20.0, f, 1.0e-9);
    camera_keyframes_evaluate_frame(&f, k, 2, -5.0);
    EXPECT_NEAR(0.0, f, 1.0e-12);
    camera_keyframes_evaluate_frame(&f, k, 2, 99.0);
    EXPECT_NEAR(100.0, f, 1.0e-12);
}

UTEST(viamd_camera, a_single_frame_key_holds_that_frame) {
    CameraKeyframe k[1];
    kf_frame(k[0], 3.0, 42.0);
    double f = 0.0;
    ASSERT_TRUE(camera_keyframes_evaluate_frame(&f, k, 1, 0.0));
    EXPECT_NEAR(42.0, f, 1.0e-12);
    camera_keyframes_evaluate_frame(&f, k, 1, 10.0);
    EXPECT_NEAR(42.0, f, 1.0e-12);
}

UTEST(viamd_camera, a_slow_then_fast_frame_curve_is_forward_only_and_hits_its_keys) {
    CameraKeyframe k[3];
    kf_frame(k[0], 0.0, 0.0);
    kf_frame(k[1], 5.0, 10.0);
    kf_frame(k[2], 10.0, 110.0);
    double prev = -1.0, f = 0.0;
    for (double t = 0.0; t <= 10.0; t += 0.01) {
        ASSERT_TRUE(camera_keyframes_evaluate_frame(&f, k, 3, t));
        EXPECT_GE(f, prev - 1.0e-9);
        EXPECT_GE(f, -1.0e-9);
        EXPECT_LE(f, 110.0 + 1.0e-9);
        prev = f;
    }
    for (int i = 0; i < 3; ++i) {
        camera_keyframes_evaluate_frame(&f, k, 3, k[i].time);
        EXPECT_NEAR(k[i].frame, f, 1.0e-9);
    }

    /* Slow in the first half, fast in the second */
    double a0, a1, b0, b1;
    camera_keyframes_evaluate_frame(&a0, k, 3, 2.4);
    camera_keyframes_evaluate_frame(&a1, k, 3, 2.5);
    camera_keyframes_evaluate_frame(&b0, k, 3, 7.4);
    camera_keyframes_evaluate_frame(&b1, k, 3, 7.5);
    EXPECT_LT(a1 - a0, 0.5 * (b1 - b0));
}

UTEST(viamd_camera, a_frame_behind_the_previous_plays_backward_and_hits_its_keys) {
    CameraKeyframe k[4];
    kf_frame(k[0], 0.0, 0.0);
    kf_frame(k[1], 2.0, 100.0);
    kf_frame(k[2], 4.0, 40.0);     /* backward from 100 */
    kf_frame(k[3], 6.0, 200.0);
    double f = 0.0;
    for (int i = 0; i < 4; ++i) {
        camera_keyframes_evaluate_frame(&f, k, 4, k[i].time);
        EXPECT_NEAR(k[i].frame, f, 1.0e-9);
    }

    /* Going down between 2 and 4, never leaving the range of its two keys */
    double prev = 100.0;
    for (double t = 2.0; t <= 4.0; t += 0.05) {
        camera_keyframes_evaluate_frame(&f, k, 4, t);
        EXPECT_LE(f, prev + 1.0e-9);
        EXPECT_GE(f, 40.0 - 1.0e-9);
        prev = f;
    }
    camera_keyframes_evaluate_frame(&f, k, 4, 3.0);
    EXPECT_LT(f, 100.0);
    EXPECT_GT(f, 40.0);
}

UTEST(viamd_camera, equal_frames_hold) {
    CameraKeyframe k[3];
    kf_frame(k[0], 0.0, 0.0);
    kf_frame(k[1], 2.0, 100.0);
    kf_frame(k[2], 4.0, 100.0);
    double f = 0.0;
    for (double t = 2.0; t <= 4.0; t += 0.05) {
        camera_keyframes_evaluate_frame(&f, k, 3, t);
        EXPECT_NEAR(100.0, f, 1.0e-9);
    }
}

UTEST(viamd_camera, frame_keys_need_not_be_sorted_or_alone) {
    CameraKeyframe sorted[3], mixed[5];
    kf_frame(sorted[0], 0.0, 0.0);
    kf_frame(sorted[1], 5.0, 10.0);
    kf_frame(sorted[2], 10.0, 110.0);

    mixed[0] = sorted[2];
    mixed[1] = CameraKeyframe{};     /* camera only, no frame */
    mixed[1].time = 3.0;
    mixed[2] = sorted[0];
    mixed[3] = CameraKeyframe{};
    mixed[3].time = 8.0;
    mixed[4] = sorted[1];

    for (double t = 0.0; t <= 10.0; t += 0.37) {
        double a = 0.0, b = 0.0;
        camera_keyframes_evaluate_frame(&a, sorted, 3, t);
        camera_keyframes_evaluate_frame(&b, mixed, 5, t);
        EXPECT_NEAR(a, b, 1.0e-12);
    }
}

/* Easing of a segment, and loops */

static void kf_two_poses(CameraKeyframe* k, KeyEase ease) {
    kf_set(k[0], 0.0, vec3_set(0.3f, 1.0f, 0.1f), 0.2f, vec3_set(0, 0, 0),   10.0f, 0.8f);
    kf_set(k[1], 4.0, vec3_set(0.3f, 1.0f, 0.1f), 1.2f, vec3_set(40, 0, 0),  20.0f, 0.6f);
    k[1].ease = ease;
}

UTEST(viamd_camera, a_hold_stays_in_the_pose_until_the_key_and_then_jumps) {
    CameraKeyframe k[2];
    kf_two_poses(k, KeyEase::Hold);
    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, k, 2, 3.99);
    EXPECT_NEAR(k[0].transform.distance, v.distance, 1.0e-4f);
    EXPECT_NEAR(k[0].fov_y, fov, 1.0e-5f);
    EXPECT_NEAR(k[0].transform.position.x, v.position.x, 1.0e-3f);
    EXPECT_NEAR(1.0f, fabsf(quat_dot(k[0].transform.orientation, v.orientation)), 1.0e-5f);
    camera_keyframes_evaluate(&v, &fov, k, 2, 4.0);
    EXPECT_NEAR(k[1].transform.distance, v.distance, 1.0e-4f);
    EXPECT_NEAR(k[1].fov_y, fov, 1.0e-5f);
}

UTEST(viamd_camera, a_linear_segment_moves_at_constant_speed_where_a_smooth_one_eases) {
    CameraKeyframe lin[2], smooth[2];
    kf_two_poses(lin, KeyEase::Linear);
    kf_two_poses(smooth, KeyEase::Smooth);

    ViewTransform v; float fov;
    camera_keyframes_evaluate(&v, &fov, lin, 2, 1.0);
    EXPECT_NEAR(12.5f, v.distance, 1.0e-3f);
    EXPECT_NEAR(0.75f, fov, 1.0e-4f);
    camera_keyframes_evaluate(&v, &fov, smooth, 2, 1.0);
    EXPECT_NEAR(10.0f + 10.0f * 0.15625f, v.distance, 1.0e-3f);   /* two keys alone ease at both ends */

    /* The orientation is turned at a constant rate too: halfway is halfway */
    camera_keyframes_evaluate(&v, &fov, lin, 2, 2.0);
    const float full = 2.0f * acosf(fminf(1.0f, fabsf(quat_dot(lin[0].transform.orientation, lin[1].transform.orientation))));
    const float half = 2.0f * acosf(fminf(1.0f, fabsf(quat_dot(lin[0].transform.orientation, v.orientation))));
    EXPECT_NEAR(0.5f * full, half, 2.0e-3f);
}

UTEST(viamd_camera, ease_in_and_out_starts_and_ends_slowly) {
    CameraKeyframe k[2];
    kf_two_poses(k, KeyEase::EaseInOut);
    ViewTransform a, b; float fov;
    camera_keyframes_evaluate(&a, &fov, k, 2, 0.0);
    camera_keyframes_evaluate(&b, &fov, k, 2, 0.02);
    const float start = vec3_length(b.position - a.position);
    camera_keyframes_evaluate(&a, &fov, k, 2, 1.99);
    camera_keyframes_evaluate(&b, &fov, k, 2, 2.01);
    const float middle = vec3_length(b.position - a.position);
    EXPECT_LT(start * 20.0f, middle);
}

static void kf_loop_keys(CameraKeyframe* k) {
    kf_set(k[0], 0.0, vec3_set(0.3f, 1.0f, 0.1f), 0.2f, vec3_set(0, 0, 0),    10.0f, 0.8f);
    kf_set(k[1], 2.0, vec3_set(0.3f, 1.0f, 0.1f), 0.8f, vec3_set(20, 10, 0),  14.0f, 0.8f);
    kf_set(k[2], 4.0, vec3_set(0.3f, 1.0f, 0.1f), 1.4f, vec3_set(40, -5, 8),  12.0f, 0.8f);
    k[3] = k[0];
    k[3].time = 6.0;
}

UTEST(viamd_camera, a_loop_has_no_corner_at_the_seam) {
    CameraKeyframe k[4];
    kf_loop_keys(k);
    const double e = 1.0e-2;   /* a smaller step drowns in the float noise of slerp over tiny angles */

    /* Not a loop: the camera eases away from the first key and into the last, it stands still at the seam */
    ViewTransform a, b; float fov;
    camera_keyframes_evaluate(&a, &fov, k, 4, 0.0, false);
    camera_keyframes_evaluate(&b, &fov, k, 4, e, false);
    const vec3_t open_start = (b.position - a.position) * (float)(1.0 / e);

    camera_keyframes_evaluate(&a, &fov, k, 4, 0.0, true);
    camera_keyframes_evaluate(&b, &fov, k, 4, e, true);
    const vec3_t v_start = (b.position - a.position) * (float)(1.0 / e);

    camera_keyframes_evaluate(&a, &fov, k, 4, 6.0 - e, true);
    camera_keyframes_evaluate(&b, &fov, k, 4, 6.0, true);
    const vec3_t v_end = (b.position - a.position) * (float)(1.0 / e);

    EXPECT_GT(vec3_length(v_start), 5.0f);
    EXPECT_GT(vec3_length(v_start), 5.0f * vec3_length(open_start));
    EXPECT_LT(vec3_length(v_end - v_start), 0.15f * vec3_length(v_start));
}

UTEST(viamd_camera, a_loop_still_hits_every_key) {
    CameraKeyframe k[4];
    kf_loop_keys(k);
    for (int i = 0; i < 4; ++i) {
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 4, k[i].time, true);
        EXPECT_NEAR(k[i].transform.position.x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.z, v.position.z, 1.0e-3f);
    }
}

/* Default view. These pin the behaviour rather than the numbers: what is up, what faces the viewer, and
 * that everything fits. Synthetic systems, deterministic (a fixed LCG, not rand()). */

struct DvSys {
    vec3_t xyz[4096];
    size_t n = 0;
    void add(vec3_t p) { xyz[n++] = p; }
};

static float dv_rand(uint32_t* s) { *s = *s * 1664525u + 1013904223u; return (float)(*s >> 8) / 16777216.0f; }

static vec3_t dv_view_dir(const ViewTransform& v) { return quat_mul_vec3(v.orientation, vec3_set(0, 0, 1)); }
static vec3_t dv_view_up (const ViewTransform& v) { return quat_mul_vec3(v.orientation, vec3_set(0, 1, 0)); }

static bool dv_all_in_view(const DvSys& s, const ViewTransform& v, float fov_y) {
    const mat4_t V = camera_world_to_view_matrix(v);
    const float t = tanf(fov_y * 0.5f);
    for (size_t i = 0; i < s.n; ++i) {
        const vec4_t p = mat4_mul_vec4(V, vec4_from_vec3(s.xyz[i], 1.0f));
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
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, &A, fov);
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
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    EXPECT_GT(dv_view_up(v).z, 0.9f);
    EXPECT_LT(fabsf(dv_view_dir(v).z), 0.5f);
}

UTEST(viamd_camera, default_view_of_a_planar_molecule_is_face_on) {
    /* Benzene, in an arbitrary orientation. Real coordinates are never exactly symmetric: with the ideal
     * geometry the visibility scores of nearby directions tie and the view can settle ~10 degrees off the
     * normal, so a deterministic jitter of 0.001 A (what a coordinate file's rounding gives) is added. */
    static DvSys s; s.n = 0;
    uint32_t seed = 5;
    const vec3_t n  = vec3_normalize(vec3_set(0.3f, -0.5f, 0.8f));
    const vec3_t e1 = vec3_normalize(vec3_cross(n, vec3_set(1, 0, 0)));
    const vec3_t e2 = vec3_cross(n, e1);
    auto jitter = [&]() { return vec3_set(0.001f * (dv_rand(&seed) - 0.5f), 0.001f * (dv_rand(&seed) - 0.5f), 0.001f * (dv_rand(&seed) - 0.5f)); };
    for (int k = 0; k < 6; ++k) {
        const float a = k * 1.0471976f;
        s.add(vec3_add(vec3_add(vec3_mul1(e1, 1.39f * cosf(a)), vec3_mul1(e2, 1.39f * sinf(a))), jitter()));
        s.add(vec3_add(vec3_add(vec3_mul1(e1, 2.47f * cosf(a)), vec3_mul1(e2, 2.47f * sinf(a))), jitter()));
    }
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
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
    const ViewTransform v = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
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
    const ViewTransform a = camera_compute_default_view(s.xyz, s.n, NULL, s.n, NULL, 0.785f);
    const ViewTransform b = camera_compute_default_view(m.xyz, m.n, NULL, m.n, NULL, 0.785f);
    EXPECT_GT(vec3_dot(dv_view_up(a), dv_view_up(b)), 0.99f);
    EXPECT_GT(vec3_dot(dv_view_dir(a), dv_view_dir(b)), 0.99f);
}

UTEST(viamd_camera, aim_at_keeps_the_eye_and_looks_at_the_point) {
    Camera c = test_camera();
    const vec3_t eye = c.position;
    const vec3_t p = eye + vec3_set(4.0f, -2.0f, -9.0f);
    ASSERT_TRUE(camera_aim_at(&c, p));
    EXPECT_NEAR(eye.x, c.position.x, 1.0e-6f);
    EXPECT_NEAR(eye.y, c.position.y, 1.0e-6f);
    EXPECT_NEAR(eye.z, c.position.z, 1.0e-6f);
    const vec3_t l = camera_get_look_at(c);
    EXPECT_NEAR(p.x, l.x, 1.0e-3f);
    EXPECT_NEAR(p.y, l.y, 1.0e-3f);
    EXPECT_NEAR(p.z, l.z, 1.0e-3f);
    EXPECT_NEAR(vec3_length(p - eye), c.distance, 1.0e-3f);
    EXPECT_FALSE(camera_aim_at(&c, eye));
}

UTEST(viamd_camera, keys_that_track_their_own_atoms_each_follow_their_atom) {
    CameraKeyframe k[3];
    kf_make3(k);
    for (int i = 0; i < 3; ++i) {
        k[i].follow = true;
        k[i].follow_atom = i;
        k[i].follow_center = vec3_set((float)i, 0, 0);
    }
    /* Each atom has moved by its own amount, and every key is hit at its own time */
    const vec3_t shifts[3] = {vec3_set(1, 0, 0), vec3_set(0, 2, 0), vec3_set(0, 0, 3)};
    vec3_t now[3];
    for (int i = 0; i < 3; ++i) now[i] = k[i].follow_center + shifts[i];
    for (int i = 0; i < 3; ++i) {
        ViewTransform v; float fov;
        camera_keyframes_evaluate(&v, &fov, k, 3, k[i].time, false, nullptr, now);
        EXPECT_NEAR(k[i].transform.position.x + shifts[i].x, v.position.x, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.y + shifts[i].y, v.position.y, 1.0e-3f);
        EXPECT_NEAR(k[i].transform.position.z + shifts[i].z, v.position.z, 1.0e-3f);
    }
}

UTEST(viamd_camera, trajectory_frame_goes_through_the_anchors_and_the_keys) {
    /* No keys: linear between the anchors, held outside */
    EXPECT_NEAR(0.0,   camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 12.0, 100.0, 0.0), 1.0e-9);
    EXPECT_NEAR(0.0,   camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 12.0, 100.0, 2.0), 1.0e-9);
    EXPECT_NEAR(50.0,  camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 12.0, 100.0, 7.0), 1.0e-9);
    EXPECT_NEAR(100.0, camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 12.0, 100.0, 20.0), 1.0e-9);

    /* A key in between is hit and the trajectory still reaches the end anchor */
    CameraKeyframe k = {};
    k.time = 4.0; k.use_frame = true; k.frame = 80.0;
    EXPECT_NEAR(80.0,  camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 4.0), 1.0e-9);
    EXPECT_NEAR(100.0, camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 12.0), 1.0e-9);
    const double mid = camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 8.0);
    EXPECT_GT(mid, 80.0);
    EXPECT_LT(mid, 100.0);

    /* A key before the begin anchor takes its place */
    k.time = 1.0; k.frame = 30.0;
    EXPECT_NEAR(30.0, camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 0.5), 1.0e-9);
    EXPECT_NEAR(30.0, camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 1.0), 1.0e-9);
    k.time = 12.0; k.frame = 60.0;
    EXPECT_NEAR(60.0, camera_keyframes_frame_with_anchors(&k, 1, 2.0, 0.0, 12.0, 100.0, 14.0), 1.0e-9);
    EXPECT_NEAR(75.0, camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 100.0, 12.0, 0.0, 4.5), 1.0e-9);
    EXPECT_NEAR(0.0, camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 2.0, 100.0, 1.0), 1.0e-9);
    EXPECT_NEAR(100.0, camera_keyframes_frame_with_anchors(nullptr, 0, 2.0, 0.0, 2.0, 100.0, 2.0), 1.0e-9);
}

static float v3_dist(vec3_t a, vec3_t b) { return vec3_length(a - b); }

UTEST(viamd_camera, levelling_keeps_where_the_camera_looks_and_makes_it_upright) {
    Camera c = test_camera();
    // Upside down: rolled half a turn about where it looks
    c.orientation = quat_normalize(c.orientation * quat_axis_angle(vec3_set(0, 0, 1), 3.0f));
    ViewTransform t = c;
    const vec3_t up = {0, 1, 0};
    const vec3_t look = camera_get_look_at(t);
    const vec3_t fwd = t.orientation * vec3_t{0, 0, -1};

    camera_level(&t, up);
    EXPECT_LT(v3_dist(t.position, c.position), 1.0e-4f);
    EXPECT_LT(v3_dist(camera_get_look_at(t), look), 1.0e-4f);
    EXPECT_LT(v3_dist(t.orientation * vec3_t{0, 0, -1}, fwd), 1.0e-5f);
    EXPECT_LT(fabsf(vec3_dot(t.orientation * vec3_t{1, 0, 0}, up)), 1.0e-5f);  // The horizon is level
    EXPECT_GT(vec3_dot(t.orientation * vec3_t{0, 1, 0}, up), 0.0f);           // and the right way up
    EXPECT_NEAR(camera_roll(t, up), 0.0f, 1.0e-4f);
}

UTEST(viamd_camera, the_roll_that_is_set_is_the_roll_that_is_read) {
    const vec3_t up = {0, 0, 1};
    for (float roll : {-2.5f, -0.4f, 0.0f, 0.7f, 3.0f}) {
        ViewTransform t = test_camera();
        camera_level(&t, up, roll);
        EXPECT_NEAR(camera_roll(t, up), roll, 1.0e-4f);
    }
    // Positive leans the camera's up to its left
    ViewTransform t = test_camera();
    camera_level(&t, up, 0.3f);
    ViewTransform level = test_camera();
    camera_level(&level, up);
    EXPECT_LT(vec3_dot(t.orientation * vec3_t{0, 1, 0}, level.orientation * vec3_t{1, 0, 0}), 0.0f);
}

UTEST(viamd_camera, an_upright_path_stays_level_through_a_spin_and_goes_through_the_rolls) {
    const vec3_t up = {0, 1, 0};
    CameraKeyframe keys[2];
    for (int i = 0; i < 2; ++i) {
        ViewTransform t = test_camera();
        t.position = vec3_set(10.0f * i, 2, 10);
        camera_aim_at(&t, vec3_set(0, 0, 0));
        keys[i].transform = t;
        keys[i].time = 4.0 * i;
    }
    // The second key is upside down and tilted, but only its roll counts
    keys[1].transform.orientation = quat_normalize(keys[1].transform.orientation * quat_axis_angle(vec3_set(0, 0, 1), 2.9f));
    keys[1].spin_turns = 1;
    keys[0].roll = 0.0f;
    keys[1].roll = 0.4f;
    for (int s = 0; s <= 8; ++s) {
        ViewTransform vt;
        float fov;
        camera_keyframes_evaluate(&vt, &fov, keys, 2, 0.5 * s, false, nullptr, nullptr, &up);
        const float r = camera_roll(vt, up);
        EXPECT_GE(r, -1.0e-4f);
        EXPECT_LE(r, 0.4f + 1.0e-4f);
    }
    ViewTransform end;
    float fov;
    camera_keyframes_evaluate(&end, &fov, keys, 2, 4.0, false, nullptr, nullptr, &up);
    EXPECT_NEAR(camera_roll(end, up), 0.4f, 1.0e-4f);
    EXPECT_LT(v3_dist(camera_get_look_at(end), camera_get_look_at(keys[1].transform)), 1.0e-3f);
}
