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

    return quat_axis_angle(axis, angle);
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

// high precision version
static inline vec3_t highp_quat_vec3_mul(const quat_t& q, const vec3_t& v) {
    double u[4] = {q.x, q.y, q.z, q.w};
    double t[3] = {
        2.0 * (u[1] * (double)v.z - (double)v.y * u[2]),
        2.0 * (u[2] * (double)v.x - (double)v.z * u[0]),
        2.0 * (u[0] * (double)v.y - (double)v.x * u[1]),
    };
    double res[3] = {
        v.x + t[0] * u[3] + (u[1] * t[2] - t[1] * u[2]),
        v.y + t[1] * u[3] + (u[2] * t[0] - t[2] * u[0]),
        v.z + t[2] * u[3] + (u[0] * t[1] - t[0] * u[1]),
    };
    return vec3_set((float)res[0], (float)res[1], (float)res[2]);
}

// We want to interpolate along an arc which is formed by maintaining a distance to the look_at position and smoothly interpolating the orientation,
// We linearly interpolate a look_at position which is implicitly defined by position, orientation and distance
// There is some precision errors creeping into the posision because we transform back and forth to look at using the orientation
void camera_interpolate_look_at(vec3_t* out_pos, quat_t* out_ori, float* out_dist, vec3_t in_pos[2], quat_t in_ori[2], float in_dist[2], double t) {

    // This is to combat floating point inaccuracies when interpolating in this arc like fashion
    // If we are close enough, we just set it to that
    const double quat_epsilon = 1.0e-7;
    const double pos_epsilon = 1.0e-7;
    const double dist_epsilon = 1.0e-7;

    quat_t ori;
    vec3_t pos;
    float dist;

    const float qd = quat_dot(in_ori[0], in_ori[1]);
    if (fabsf(qd) > 1.0 - quat_epsilon) {
        ori = in_ori[1];
    } else {
        // Due to the inherent rotational duality of quaternions we want to make sure we rotate along the shortest 'path'
        quat_t qa = qd < 0.0f ? quat_conj(in_ori[0]) : in_ori[0];
        quat_t qb = in_ori[1];
        ori = quat_normalize(quat_slerp(qa, qb, (float)t));
    }

    if (fabs((double)in_dist[0] - (double)in_dist[1]) < dist_epsilon) {
        dist = in_dist[1];
    } else {
        dist = (float)lerp(in_dist[0], in_dist[1], t);
    }

    double dx = in_pos[0].x - in_pos[1].x;
    double dy = in_pos[0].y - in_pos[1].y;
    double dz = in_pos[0].z - in_pos[1].z;
    double d2 = sqrt(dx*dx + dy*dy + dz*dz);
    if (d2 < pos_epsilon) {
        pos = in_pos[1];
    } else {
        vec3_t l0 = in_pos[0] - highp_quat_vec3_mul(in_ori[0], vec3_set(0, 0, in_dist[0]));
        vec3_t l1 = in_pos[1] - highp_quat_vec3_mul(in_ori[1], vec3_set(0, 0, in_dist[1]));
        vec3_t look_at = {
            lerp(l0.x, l1.x, t),
            lerp(l0.y, l1.y, t),
            lerp(l0.z, l1.z, t),
        };
        pos = look_at + ori * vec3_t{0, 0, dist};
    }

    *out_pos = pos;
    *out_ori = ori;
    *out_dist = dist;
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

void camera_compute_optimal_view(vec3_t* out_pos, quat_t* out_ori, float* out_dist, mat3_t in_basis, vec3_t in_min_ext, vec3_t in_max_ext, float distance_scale) {
    const vec3_t ext = in_max_ext - in_min_ext;
    const float len = MAX(vec3_length(ext * 0.5f), 5.0f);

    const float max_ext = MAX(MAX(ext.x, ext.y), ext.z);
    const float min_ext = MIN(MIN(ext.x, ext.y), ext.z);
    const float aniso_ext = max_ext / min_ext;

    // We want to align the view such that we the longest axis of the aabb align with the X-axis, the mid axis with the Y-axis

    int l[3] = { 0, 1, 2 };

    auto swap = [](int& a, int& b) { int t = a; a = b; b = t; };

    if (aniso_ext > 1.1f) {
        // The aabb is not uniform, so we sort the axes by length
        if (ext[l[0]] < ext[l[1]]) swap(l[0], l[1]);
        if (ext[l[1]] < ext[l[2]]) swap(l[1], l[2]);
        if (ext[l[0]] < ext[l[1]]) swap(l[0], l[1]);
        // Now the axes are sorted with respect to the length l[0] > l[1] > l[2]
    }

    const vec3_t right = in_basis[l[0]];
    const vec3_t up    = in_basis[l[1]];
    const vec3_t out   = in_basis[l[2]];

    const vec3_t dir = vec3_normalize(right * 0.6f + up * 0.5f + out * 1.0f);

    const vec3_t cen = in_basis * ((in_min_ext + in_max_ext) * 0.5f);
    const vec3_t pos = cen + dir * len * distance_scale;

    *out_ori = quat_from_mat4(mat4_look_at(pos, cen, up));
    *out_pos = pos;
    *out_dist = vec3_length(pos - cen);
}

void camera_animate(Camera* camera, quat_t target_ori, vec3_t target_pos, float target_dist, double dt, double target_factor) {
    ASSERT(camera);

    dt = CLAMP(dt, 1.0 / 1000.0, 1.0 / 20.0);

    // We use an exponential interpolation of the deltas with a common factor
    const double INV_TARGET_DT = 100.0;
    double interpolation_factor = target_factor * dt * INV_TARGET_DT;

    vec3_t pos[2] = {camera->position, target_pos};
    quat_t ori[2] = {camera->orientation, target_ori};
    float dist[2] = {camera->focus_distance, target_dist};
    camera_interpolate_look_at(&camera->position, &camera->orientation, &camera->focus_distance, pos, ori, dist, interpolation_factor);
}