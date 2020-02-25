#include <core/common.h>
#include <core/math_utils.h>
#include "camera_utils.h"

static mat4 frustum(float left, float right, float bottom, float top, float near, float far) {
    mat4 M{0};
    M[0][0] = (2.0f * near) / (right - left);
    M[1][1] = (2.0f * near) / (top - bottom);
    M[2][0] = (right + left) / (right - left);
    M[2][1] = (top + bottom) / (top - bottom);
    M[2][2] = -(far + near) / (far - near);
    M[2][3] = -1.0f;
    M[3][2] = -(2.0f * far * near) / (far - near);
    return M;
}

mat4 perspective(float fovy, float aspect, float near, float far) {
    ASSERT(math::abs(aspect - FLT_EPSILON) > 0.0f);
    const float tan_half_fovy = math::tan(fovy * 0.5f);

    mat4 M{0}; 
    M[0][0] = 1.0f / (aspect * tan_half_fovy);
    M[1][1] = 1.0f / (tan_half_fovy);
    M[2][2] = -(far + near) / (far - near);
    M[2][3] = -1.0f;
    M[3][2] = -(2.0f * far * near) / (far - near);
    return M;
}

static vec4 projection_extents(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const float aspect_ratio = (float)width / (float)height;
    const float half_h = math::tan(camera.fov_y * 0.5f);
    const float half_w = aspect_ratio * half_h;
    const float texel_size_x = half_w / (float)(0.5f * width);
    const float texel_size_y = half_h / (float)(0.5f * height);
    const float jitter_x = texel_size_x * texel_offset_x;
    const float jitter_y = texel_size_y * texel_offset_y;

    // xy = frustum extents at distance 1, zw = jitter at distance 1
    return vec4(half_w, half_h, jitter_x, jitter_y);
}

mat4 compute_view_to_world_matrix(const Camera& camera) {
    auto m = math::mat4_cast(camera.orientation);
    m[3] = vec4(camera.orientation * camera.position, 1);
    return m;
    // auto t = math::translate(mat4(1), camera.position);
    // return t * r;
}

mat4 compute_world_to_view_matrix(const Camera& camera) {
    const mat4 R = glm::mat4_cast(glm::conjugate(camera.orientation));
    const mat4 T = glm::translate(mat4(1), -camera.position);
    return R * T;
}

mat4 perspective_projection_matrix(const Camera& camera, int width, int height) {
    const float aspect_ratio = (float)width / (float)height;
    return perspective(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

mat4 perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const vec4 ext = projection_extents(camera, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return frustum(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

mat4 orthographic_projection_matrix(float l, float r, float b, float t) {
    mat4 M{};
    M[0][0] = 2.0f / (r - l);
    M[1][1] = 2.0f / (t - b);
    M[2][2] = -1.0f;
    M[3][0] = -(r + l) / (r - l);
    M[3][1] = -(t + b) / (t - b);
    M[3][3] = 1.0f;
    return M;
}

mat4 orthographic_projection_matrix(float l, float r, float b, float t, float n, float f) {
    mat4 M{};
    M[0][0] = 2.0f / (r - l);
    M[1][1] = 2.0f / (t - b);
    M[2][2] = -2.0f / (f - n);
    M[3][0] = -(r + l) / (r - l);
    M[3][1] = -(t + b) / (t - b);
    M[3][2] = -(f + n) / (f - n);
    M[3][3] = 1.0f;
    return M;
}

mat3 look_at(const vec3& look_from, const vec3& look_at, const vec3& look_up) {
    const vec3 f(math::normalize(look_at - look_from));
    const vec3 s(math::normalize(math::cross(f, look_up)));
    const vec3 u(math::cross(s, f));
	const mat4 M = {{s.x, u.x, -f.x, 0.0f}, {s.y, u.y, -f.y, 0.0f}, {s.z, u.z, -f.z, 0.0f}, {-math::dot(s, look_from), -math::dot(u, look_from), math::dot(f, look_from), 1.0f}};
    return M;
}

static inline float project_to_sphere(float r, vec2 v) {
    const float d = math::length(v);
    if (d < r * 0.70710678118654752440f) {
        // On sphere
        return math::sqrt(r * r - d * d);
    } else {
        // On hyperbola
        float t = r / 1.41421356237309504880f;
        return t * t / d;
    }
}

static inline quat trackball(vec2 prev_ndc, vec2 curr_ndc) {
    constexpr float TRACKBALLSIZE = 0.8f;
    if (dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return quat(1, 0, 0, 0);

    const vec3 p1 = vec3(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
    const vec3 p2 = vec3(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));

    const vec3 axis = math::normalize(cross(p2, p1));
    float t = math::clamp(math::length(p1 - p2) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
    float angle = 2.f * math::asin(t);

    return math::angle_axis(angle, axis);
}

void camera_trackball(Camera* camera, vec2 prev_ndc, vec2 curr_ndc) {
    ASSERT(camera);
    const quat q = trackball(prev_ndc, curr_ndc);
    camera->orientation = camera->orientation * q;
}

void camera_move(Camera* camera, vec3 vec) {
    ASSERT(camera);
    const mat3 m = compute_view_to_world_matrix(*camera);
    camera->position += m * vec;
}

bool camera_controller_trackball(vec3* position, quat* orientation, TrackballControllerState* state, TrackballFlags flags) {
    ASSERT(position);
    ASSERT(orientation);
    ASSERT(state);

    const vec2 half_res = state->input.screen_size * 0.5f;
    const vec2 ndc_prev = (vec2(state->input.mouse_coord_prev.x, state->input.screen_size.y - state->input.mouse_coord_prev.y) - half_res) / half_res;
    const vec2 ndc_curr = (vec2(state->input.mouse_coord_curr.x, state->input.screen_size.y - state->input.mouse_coord_curr.y) - half_res) / half_res;
    const vec2 mouse_coord_delta = state->input.mouse_coord_curr - state->input.mouse_coord_prev;
    const bool mouse_move = mouse_coord_delta != vec2(0, 0);

    if (state->input.rotate_button && mouse_move) {
        const quat q = trackball(ndc_prev, ndc_curr);
        const vec3 look_at = *position - *orientation * vec3(0, 0, state->distance);
        *orientation = glm::normalize(*orientation * q);
        *position = look_at + *orientation * vec3(0, 0, state->distance);
        if (flags & TrackballFlags_RotateReturnsTrue) return true;
    } else if (state->input.pan_button && mouse_move) {
        const float aspect_ratio = state->input.screen_size.x / state->input.screen_size.y;
        const float scl = math::tan(state->input.fov_y * 0.5f);
        const vec2 delta = (ndc_curr - ndc_prev) * vec2(1, -1) * vec2(aspect_ratio * scl, scl);
        const vec3 move = *orientation * vec3(-delta.x, delta.y, 0) * math::pow(state->distance * state->params.pan_scale, state->params.pan_exponent);
        *position += move;
        if (flags & TrackballFlags_PanReturnsTrue) return true;
    } else if ((state->input.dolly_button && mouse_move) || state->input.dolly_delta != 0.f) {
        float delta = -(state->input.mouse_coord_curr.y - state->input.mouse_coord_prev.y) * math::pow(state->distance * state->params.dolly_drag_scale, state->params.dolly_drag_exponent);
        delta -= state->input.dolly_delta * math::pow(state->distance * state->params.dolly_delta_scale, state->params.dolly_delta_exponent);
        const vec3 look_at = *position - *orientation * vec3(0, 0, state->distance);
        state->distance = math::clamp(state->distance + delta, state->params.min_distance, state->params.max_distance);
        *position = look_at + *orientation * vec3(0, 0, state->distance);
        if (flags & TrackballFlags_DollyReturnsTrue) return true;
    }
    return false;
}
