#include <core/camera_utils.h>
#include <core/common.h>
#include <core/math_utils.h>

mat4 compute_view_to_world_matrix(const Camera& camera) {
	auto r = glm::mat4_cast(camera.orientation);
	auto t = glm::translate(mat4(1), camera.position);
	return  t * r;
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
mat4 compute_orthographic_projection_matrix(const Camera & camera, int width, int height) {
	float h_w = width * 0.05f;
	float h_h = height * 0.05f;
	return glm::ortho(-h_w, h_w, -h_h, h_h, camera.near_plane, camera.far_plane);
}

void camera_look_at(Camera* camera, vec3 look_at, vec3 look_up) {
	ASSERT(camera);
	vec3 look = camera->position - look_at;
	float len2 = dot(look, look);
	if (len2 > 0.f) {
		look /= sqrtf(len2);
		vec3 right = normalize(cross(look_up, look));
		// @TODO: Make sure look and look_up dont coincide
		look_up = cross(look, right);

		camera->orientation = glm::quat_cast(mat3(right, look_up, look));
	}

	// @TODO: make sure look_up is kept.
}

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel, float delta_time, float move_speed,
                           float rot_speed) {
	(void)rot_speed;
    ASSERT(camera);

	mouse_vel *= 0.01f;
	quat pitch = quat(vec3(-mouse_vel.y, 0, 0));
	quat yaw = quat(vec3(0, mouse_vel.x, 0));
	//glm::clamp(euler.x, -math::PI * 0.5f, math::PI * 0.5f);
	camera->orientation = yaw * camera->orientation * pitch;
	camera->orientation = math::normalize(camera->orientation);
    //camera->orientation = quat(mouse_vel.x, vec3(0, 1, 0)) * camera->orientation * quat(mouse_vel.y, vec3(1, 0, 0));

    vec3 move_vec(0);
	move_speed *= 10.f * delta_time;

    if (key_fwd) move_vec.z -= move_speed;
    if (key_bwd) move_vec.z += move_speed;
    if (key_lft) move_vec.x -= move_speed;
    if (key_rht) move_vec.x += move_speed;

	mat3 mat = math::mat3_cast(camera->orientation);
	camera->position += mat * move_vec;
}

static float project_to_sphere(float r, vec2 v) {
	float d = math::length(v);
	if (d < r * 0.70710678118654752440f) { 
		// On sphere
		return math::sqrt(r*r - d*d);
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

static inline vec3 transform_vec(quat q, vec3 v) {
	return q * v;
}

void TrackballController::update() {
	constexpr float PAN_SCL = 0.0008f;
	constexpr float PAN_EXP = 1.f;
	constexpr float DOLLY_DRAG_SCL = 0.01f;
	constexpr float DOLLY_DRAG_EXP = 1.1f;
	constexpr float DOLLY_DELTA_SCL = 0.1f;
	constexpr float DOLLY_DELTA_EXP = 1.1f;

	if (input.rotate_button) {
		const vec2 half_res = input.screen_size * 0.5f;
		vec2 ndc_prev = (vec2(input.mouse_coord_prev.x, input.screen_size.y - input.mouse_coord_prev.y) - half_res) / half_res;
		vec2 ndc_curr = (vec2(input.mouse_coord_curr.x, input.screen_size.y - input.mouse_coord_curr.y) - half_res) / half_res;
		quat q = trackball(ndc_prev, ndc_curr);
		vec3 look_at = compute_look_at();
		orientation = glm::normalize(orientation * q);
		position = look_at + transform_vec(orientation, vec3(0, 0, distance));
	}
	else if (input.pan_button) {
		vec2 delta = input.mouse_coord_curr - input.mouse_coord_prev;
		vec3 move = vec3(-delta.x, delta.y, 0);
		move = transform_vec(orientation, move) * math::pow(distance * PAN_SCL, PAN_EXP);
		position += move;
	}
	else if (input.dolly_button || input.dolly_delta != 0.f) {
		float delta = -(input.mouse_coord_curr.y - input.mouse_coord_prev.y) * math::pow(distance * DOLLY_DRAG_SCL, DOLLY_DRAG_EXP);
		delta -= input.dolly_delta * math::pow(distance * DOLLY_DELTA_SCL, DOLLY_DELTA_EXP);
		vec3 look_at = compute_look_at();
		distance = math::clamp(distance + delta, min_distance, max_distance);
		position = look_at + transform_vec(orientation, vec3(0, 0, distance));
	}
}
