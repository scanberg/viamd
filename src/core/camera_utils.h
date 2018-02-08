#pragma once

#include <core/camera.h>
#include <core/math_utils.h>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel,
                           float delta_time, float move_speed = 1.f, float rot_speed = 1.f);

void camera_trackball(Camera* camera, vec2 prev_ndc, vec2 curr_ndc);
void camera_move(Camera* camera, vec3 vec);

namespace camera {
	struct InputState {
		bool rotate_button = false;
		bool pan_button = false;
		bool dolly_button = false;
		vec2 prev_mouse_ndc = vec2(0,0);
		vec2 curr_mouse_ndc = vec2(0,0);
		float dolly_delta = 0;
	};

	struct TrackballController {
		InputState input;
		
		vec3 position = vec3(0,10,10);
		float distance = glm::length(position);
		quat orientation = quat(1,0,0,0);

		float min_distance = 1.f;
		float max_distance = 1000.f;

		vec3 compute_look_at() const {
			return position - glm::mat3_cast(orientation) * vec3(0, 0, distance);
		}

		void update();
	};
}
