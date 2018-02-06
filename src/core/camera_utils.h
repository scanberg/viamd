#pragma once

#include <core/camera.h>
#include <core/math_utils.h>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel,
                           float delta_time, float move_speed = 1.f, float rot_speed = 1.f);

void camera_controller_trackball(Camera* camera, bool lft_mouse, bool mdl_mouse, bool rht_mouse, vec2 mouse_vel, int screen_w, int screen_h);

struct Trackball {
	float angle;
	float prev_angle;

	vec2 curr_move;
	vec2 prev_move;

	vec3 last_axis;

	vec3 target = vec3(0, 0, 0);

	void update(Camera* camera, bool rot_mouse_btn, bool zoom_mouse_btn, bool pan_mouse_btn, vec2 mouse_vel, float mouse_zoom_delta, int screen_w, int screen_h) {
		vec3 move_dir{ mouse_vel.x / (float)screen_w, -mouse_vel.y / (float)screen_h, 0 };
		angle = glm::length(move_dir);
		vec3 up = glm::normalize(glm::mat3_cast(glm::conjugate(camera->orientation))[1]);
		vec3 eye = camera->position - target;
		vec3 eye_dir = glm::normalize(eye);

		constexpr float ROT_SPEED = 1.f;
		constexpr float PAN_SPEED = 0.01f;
		constexpr float ZOOM_SPEED = 0.05f;

		if (angle && rot_mouse_btn) {

			vec3 object_up_dir = up;
			vec3 object_sideways_dir = glm::cross(object_up_dir, eye_dir);

			object_up_dir *= move_dir.y;
			object_sideways_dir *= move_dir.x;

			move_dir = object_up_dir + object_sideways_dir;
			vec3 axis = glm::normalize(glm::cross(move_dir, eye));

			angle *= ROT_SPEED;
			quat quaternion = glm::angleAxis(angle, axis);

			eye = quaternion * eye * glm::conjugate(quaternion);
			up = quaternion * up * glm::conjugate(quaternion);
		}
		else if (zoom_mouse_btn) {

		}
		else if (pan_mouse_btn) {
			if (glm::dot(mouse_vel, mouse_vel) > 0.f) {
				vec2 mouse_change = mouse_vel * glm::length(eye) * PAN_SPEED;

				vec3 pan = cross(eye_dir, up) * mouse_change.x + up * mouse_change.y;
				target += pan;
			}
		}

		if (mouse_zoom_delta != 0) {
			float scl = glm::distance(target, eye);
			eye -= eye_dir * mouse_zoom_delta * scl * ZOOM_SPEED;
		}

		camera->position = target + eye;
		auto m = glm::inverse(glm::lookAt(eye, target, up));
		camera->orientation = glm::quat_cast(m);
		camera->orientation = glm::normalize(camera->orientation);


		//camera_look_at(camera, target, object_up_dir);
	}
};