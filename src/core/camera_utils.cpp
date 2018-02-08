#include <core/camera_utils.h>
#include <core/common.h>
#include <core/math_utils.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/euler_angles.hpp>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel, float delta_time, float move_speed,
                           float rot_speed) {
    ASSERT(camera);

	mouse_vel *= 0.01f;
	quat pitch = glm::quat(glm::vec3(-mouse_vel.y, 0, 0));
	quat yaw = glm::quat(glm::vec3(0, mouse_vel.x, 0));
	//glm::clamp(euler.x, -math::PI * 0.5f, math::PI * 0.5f);
	camera->orientation = yaw * camera->orientation * pitch;
	camera->orientation = glm::normalize(camera->orientation);
    //camera->orientation = quat(mouse_vel.x, vec3(0, 1, 0)) * camera->orientation * quat(mouse_vel.y, vec3(1, 0, 0));

    vec3 move_vec(0);
	move_speed *= 10.f * delta_time;

    if (key_fwd) move_vec.z -= move_speed;
    if (key_bwd) move_vec.z += move_speed;
    if (key_lft) move_vec.x -= move_speed;
    if (key_rht) move_vec.x += move_speed;

	mat3 mat = glm::mat3_cast(camera->orientation);
	camera->position += mat * move_vec;
}

static float project_to_sphere(float r, vec2 v) {
	float d = length(v);
	if (d < r * 0.70710678118654752440f) { 
		// On sphere
		return sqrt(r*r - d*d);
	}
	else {
		// On hyperbola
		float t = r / 1.41421356237309504880f;
		return t*t / d;
	}
}

static quat trackball(vec2 prev_ndc, vec2 curr_ndc) {
	constexpr float TRACKBALLSIZE = 0.8f;
	if (dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return quat(1,0,0,0);

	vec3 p1 = vec3(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
	vec3 p2 = vec3(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));
	vec3 d = p1 - p2;

	vec3 axis = glm::normalize(cross(p2, p1));
	float t = math::clamp(length(d) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
	float angle = 2.f * math::asin(t);

	return glm::angleAxis(angle, axis);
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

void camera::TrackballController::update() {
	constexpr float PAN_SCL = 10.f;
	constexpr float DOLLY_DRAG_SCL = 50.f;
	constexpr float DOLLY_DELTA_SCL = 0.8f;

	if (input.rotate_button) {
		quat q = trackball(input.prev_mouse_ndc, input.curr_mouse_ndc);
		vec3 look_at = compute_look_at();
		orientation = glm::normalize(orientation * q);
		position = look_at + glm::mat3_cast(orientation) * vec3(0,0,distance);
	}
	else if (input.pan_button) {
		vec3 move = vec3(-(input.curr_mouse_ndc - input.prev_mouse_ndc), 0);
		move = glm::mat3_cast(orientation) * move * sqrtf(distance * PAN_SCL);
		position += move;
	}
	else if (input.dolly_button || input.dolly_delta != 0.f) {
		float delta = -(input.curr_mouse_ndc.y - input.prev_mouse_ndc.y) * sqrtf(distance * DOLLY_DRAG_SCL);
		delta -= input.dolly_delta * sqrtf(distance * DOLLY_DELTA_SCL);
		vec3 look_at = compute_look_at();
		distance = math::clamp(distance + delta, min_distance, max_distance);
		position = look_at + glm::mat3_cast(orientation) * vec3(0, 0, distance);
	}
}
