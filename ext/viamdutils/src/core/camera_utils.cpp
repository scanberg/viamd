#include <core/camera_utils.h>
#include <core/common.h>
#include <core/math_utils.h>

mat4 compute_view_to_world_matrix(const Camera& camera) {
    auto m = math::mat4_cast(camera.orientation);
    m[3] = vec4(camera.orientation * camera.position, 1);
    return m;
    // auto t = math::translate(mat4(1), camera.position);
    // return t * r;
}

mat4 compute_world_to_view_matrix(const Camera& camera) {
    auto r = glm::mat4_cast(glm::conjugate(camera.orientation));
    auto t = glm::translate(mat4(1), -camera.position);
    return r * t;
}

mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height) {
    float aspect = (float)width / (float)height;
    return glm::perspective(camera.fov_y, aspect, camera.near_plane, camera.far_plane);
}

// @TODO: This is messed up... what values should one use to control the zoomlevel?
mat4 compute_orthographic_projection_matrix(const Camera& camera, int width, int height) {
    float h_w = width * 0.05f;
    float h_h = height * 0.05f;
    return glm::ortho(-h_w, h_w, -h_h, h_h, camera.near_plane, camera.far_plane);
}

void look_at(vec3* position, quat* orientation, vec3 look_at, vec3 look_up) {
    ASSERT(position);
    ASSERT(orientation);

    vec3 look = *position - look_at;
    float len2 = dot(look, look);
    if (len2 > 0.f) {
        look /= sqrtf(len2);
        vec3 right = normalize(cross(look_up, look));
        // @TODO: Make sure look and look_up dont coincide
        look_up = cross(look, right);

        *orientation = glm::quat_cast(mat3(right, look_up, look));
    }

    // @TODO: make sure look_up is kept.
}

void camera_controller_fps(Camera* camera, const FpsControllerState& state) {
    (void)camera;
    (void)state;
    /*
ASSERT(camera);

mouse_vel *= 0.01f;
quat pitch = quat(vec3(-mouse_vel.y, 0, 0));
quat yaw = quat(vec3(0, mouse_vel.x, 0));
// glm::clamp(euler.x, -math::PI * 0.5f, math::PI * 0.5f);
camera->orientation = yaw * camera->orientation * pitch;
camera->orientation = math::normalize(camera->orientation);
// camera->orientation = quat(mouse_vel.x, vec3(0, 1, 0)) * camera->orientation * quat(mouse_vel.y, vec3(1, 0, 0));

vec3 move_vec(0);
move_speed *= 10.f * delta_time;

if (key_fwd) move_vec.z -= move_speed;
if (key_bwd) move_vec.z += move_speed;
if (key_lft) move_vec.x -= move_speed;
if (key_rht) move_vec.x += move_speed;

mat3 mat = math::mat3_cast(camera->orientation);
camera->position += mat * move_vec;
    */
}

static float project_to_sphere(float r, vec2 v) {
    float d = math::length(v);
    if (d < r * 0.70710678118654752440f) {
        // On sphere
        return math::sqrt(r * r - d * d);
    } else {
        // On hyperbola
        float t = r / 1.41421356237309504880f;
        return t * t / d;
    }
}

static quat trackball(vec2 prev_ndc, vec2 curr_ndc) {
    constexpr float TRACKBALLSIZE = 0.8f;
    if (dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return quat(1, 0, 0, 0);

    vec3 p1 = vec3(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
    vec3 p2 = vec3(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));
    vec3 d = p1 - p2;

    vec3 axis = math::normalize(cross(p2, p1));
    float t = math::clamp(math::length(d) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
    float angle = 2.f * math::asin(t);

    return math::angle_axis(angle, axis);
}

void camera_trackball(Camera* camera, vec2 prev_ndc, vec2 curr_ndc) {
    ASSERT(camera);
    quat q = trackball(prev_ndc, curr_ndc);
    camera->orientation = camera->orientation * q;
}

void camera_move(Camera* camera, vec3 vec) {
    ASSERT(camera);
    mat3 m = compute_view_to_world_matrix(*camera);
    camera->position += m * vec;
}

static inline vec3 transform_vec(quat q, vec3 v) { return q * v; }

bool camera_controller_trackball(vec3* position, quat* orientation, TrackballControllerState* state) {
    ASSERT(camera);
    ASSERT(state);

    const vec2 half_res = state->input.screen_size * 0.5f;
    vec2 ndc_prev = (vec2(state->input.mouse_coord_prev.x, state->input.screen_size.y - state->input.mouse_coord_prev.y) - half_res) / half_res;
    vec2 ndc_curr = (vec2(state->input.mouse_coord_curr.x, state->input.screen_size.y - state->input.mouse_coord_curr.y) - half_res) / half_res;

    if (state->input.rotate_button) {
        const quat q = trackball(ndc_prev, ndc_curr);
        const vec3 look_at = *position - *orientation * vec3(0, 0, state->distance);
        *orientation = glm::normalize(*orientation * q);
        *position = look_at + *orientation * vec3(0, 0, state->distance);
        return true;
    } else if (state->input.pan_button) {
        const vec2 delta = (ndc_curr - ndc_prev) * vec2(1, -1);
        const vec3 move =
            *orientation * vec3(-delta.x, delta.y, 0) * math::pow(state->distance * state->params.pan_scale, state->params.pan_exponent);
        *position += move;
        return true;
    } else if (state->input.dolly_button || state->input.dolly_delta != 0.f) {
        float delta = -(state->input.mouse_coord_curr.y - state->input.mouse_coord_prev.y) *
                      math::pow(state->distance * state->params.dolly_drag_scale, state->params.dolly_drag_exponent);
        delta -= state->input.dolly_delta * math::pow(state->distance * state->params.dolly_delta_scale, state->params.dolly_delta_exponent);
        const vec3 look_at = *position - *orientation * vec3(0, 0, state->distance);
        state->distance = math::clamp(state->distance + delta, state->params.min_distance, state->params.max_distance);
        *position = look_at + transform_vec(*orientation, vec3(0, 0, state->distance));
        return true;
    }
    return false;
}
