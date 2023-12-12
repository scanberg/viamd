#include "camera_utils.h"
#include <core/md_common.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>

static vec4_t projection_extents(float fov_y, int width, int height, float texel_offset_x, float texel_offset_y) {
    const float aspect_ratio = (float)width / (float)height;
    const float half_h = tanf(fov_y * 0.5f);
    const float half_w = aspect_ratio * half_h;
    const float texel_size_x = half_w / (float)(0.5f * width);
    const float texel_size_y = half_h / (float)(0.5f * height);
    const float jitter_x = texel_size_x * texel_offset_x;
    const float jitter_y = texel_size_y * texel_offset_y;

    // xy = frustum extents at distance 1, zw = jitter at distance 1
    return {half_w, half_h, jitter_x, jitter_y};
}

mat4_t camera_view_to_world_matrix(const Camera& camera) {
    mat4_t M = mat4_from_quat(camera.orientation);
    M.col[3] = vec4_from_vec3(camera.position, 1);
    return M;
}

mat4_t camera_world_to_view_matrix(const Camera& camera) {
    const mat4_t R = mat4_from_quat(quat_conj(camera.orientation));
    const mat4_t T = mat4_translate(-camera.position.x, -camera.position.y, -camera.position.z);
    return R * T;
}

mat4_t camera_perspective_projection_matrix(const Camera& camera, float aspect_ratio) {
    return mat4_persp(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

mat4_t camera_inverse_perspective_projection_matrix(const Camera& camera, float aspect_ratio) {
    return mat4_persp_inv(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

mat4_t camera_perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const vec4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return mat4_frustum(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

mat4_t camera_inverse_perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const vec4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return mat4_frustum_inv(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

mat4_t camera_orthographic_projection_matrix(float l, float r, float b, float t) {
    return mat4_ortho_2d(l, r, b, t);
}

mat4_t camera_inverse_orthographic_projection_matrix(float l, float r, float b, float t) {
    return mat4_ortho_2d_inv(l, r, b, t);
}

mat4_t camera_orthographic_projection_matrix(float l, float r, float b, float t, float n, float f) {
    return mat4_ortho(l, r, b, t, n, f);
}

mat4_t camera_inverse_orthographic_projection_matrix(float l, float r, float b, float t, float n, float f) {
    return mat4_ortho_inv(l, r, b, t, n, f);
}

mat4_t look_at(vec3_t look_from, vec3_t look_at, vec3_t look_up) {
    const vec3_t f = vec3_normalize(look_at - look_from);
    const vec3_t s = vec3_normalize(vec3_cross(f, look_up));
    const vec3_t u = vec3_cross(s, f);
	const mat4_t M = {
            s.x, u.x, -f.x, 0.0f,
            s.y, u.y, -f.y, 0.0f,
            s.z, u.z, -f.z, 0.0f,
            -vec3_dot(s, look_from), -vec3_dot(u, look_from), vec3_dot(f, look_from), 1.0f,
    };
    return M;
}

static inline float project_to_sphere(float r, vec2_t v) {
    const float d = vec2_length(v);
    if (d < r * 0.70710678118654752440f) {
        // On sphere
        return sqrtf(r * r - d * d);
    } else {
        // On hyperbola
        float t = r / 1.41421356237309504880f;
        return t * t / d;
    }
}

static inline quat_t trackball(vec2_t prev_ndc, vec2_t curr_ndc) {
    static const float TRACKBALLSIZE = 0.8f;
    if (vec2_dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return quat_t{1, 0, 0, 0};

    const vec3_t p1 = vec3_from_vec2(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
    const vec3_t p2 = vec3_from_vec2(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));

    const vec3_t axis = vec3_normalize(vec3_cross(p2, p1));
    float t = CLAMP(vec3_length(p1 - p2) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
    float angle = 2.f * asinf(t);

    return quat_angle_axis(angle, axis);
}

void camera_trackball(Camera* camera, vec2_t prev_ndc, vec2_t curr_ndc) {
    ASSERT(camera);
    const quat_t q = trackball(prev_ndc, curr_ndc);
    camera->orientation = camera->orientation * q;
}

void camera_move(Camera* camera, vec3_t t) {
    ASSERT(camera);
    camera->position = camera->position + camera->orientation * t;
}

bool camera_controller_trackball(vec3_t* position, quat_t* orientation, float* distance, TrackballControllerInput input, TrackballControllerParam param, TrackballFlags flags) {
    ASSERT(position);
    ASSERT(orientation);
    ASSERT(distance);

    const vec2_t half_res = input.screen_size * 0.5f;
    const vec2_t ndc_prev = (vec2_t{input.mouse_coord_prev.x, input.screen_size.y - input.mouse_coord_prev.y} - half_res) / half_res;
    const vec2_t ndc_curr = (vec2_t{input.mouse_coord_curr.x, input.screen_size.y - input.mouse_coord_curr.y} - half_res) / half_res;
    const vec2_t mouse_coord_delta = input.mouse_coord_curr - input.mouse_coord_prev;
    const bool mouse_move = mouse_coord_delta != vec2_t{0, 0};

    if ((flags & TrackballFlags_RotateEnabled) && input.rotate_button && mouse_move) {
        const quat_t q = trackball(ndc_prev, ndc_curr);
        const vec3_t look_at = *position - *orientation * vec3_t{0, 0, *distance};
        *orientation = quat_normalize(*orientation * q);
        *position = look_at + *orientation * vec3_t{0, 0, *distance};
        if (flags & TrackballFlags_RotateReturnsTrue) return true;
    } else if ((flags & TrackballFlags_PanEnabled) && input.pan_button && mouse_move) {
        const float aspect_ratio = input.screen_size.x / input.screen_size.y;
        const float scl = tanf(input.fov_y * 0.5f);
        const vec2_t delta = (ndc_curr - ndc_prev) * vec2_t{1, -1} * vec2_t{aspect_ratio * scl, scl};
        const vec3_t move = *orientation * vec3_t{-delta.x, delta.y, 0} * powf(*distance * param.pan_scale, param.pan_exponent);
        *position = *position + move;
        if (flags & TrackballFlags_PanReturnsTrue) return true;
    } else if ((flags & TrackballFlags_DollyEnabled) && ((input.dolly_button && mouse_move) || input.dolly_delta != 0.f)) {
        float delta = -(input.mouse_coord_curr.y - input.mouse_coord_prev.y) * powf(*distance * param.dolly_drag_scale, param.dolly_drag_exponent);
        delta -= input.dolly_delta * powf(*distance * param.dolly_delta_scale, param.dolly_delta_exponent);
        const vec3_t look_at = *position - *orientation * vec3_t{0, 0, *distance};
        *distance = CLAMP(*distance + delta, param.min_distance, param.max_distance);
        *position = look_at + *orientation * vec3_t{0, 0, *distance};
        if (flags & TrackballFlags_DollyReturnsTrue) return true;
    }
    return false;
}
