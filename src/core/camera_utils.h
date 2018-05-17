#pragma once

#include <core/camera.h>
#include <core/math_utils.h>

/*
void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel,
                           float delta_time, float move_speed = 1.f, float rot_speed = 1.f);
						   */

void camera_trackball(Camera* camera, vec2 prev_ndc, vec2 curr_ndc);
void camera_move(Camera* camera, vec3 vec);
void camera_look_at(Camera* camera, vec3 look_at, vec3 look_up);

mat4 compute_world_to_view_matrix(const Camera& camera);
mat4 compute_view_to_world_matrix(const Camera& camera);
mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height);
mat4 compute_orthographic_projection_matrix(const Camera& camera, int width, int height);

struct TrackballController {
    struct InputState {
        bool rotate_button = false;
        bool pan_button = false;
        bool dolly_button = false;
		vec2 mouse_coord_prev = vec2(0, 0);
		vec2 mouse_coord_curr = vec2(0, 0);
		vec2 screen_size = vec2(0, 0);
        float dolly_delta = 0;
    };
	InputState input;
	
	vec3 position = vec3(0,10,10);
	float distance = glm::length(position);
	quat orientation = quat(1,0,0,0);

	float min_distance = 1.f;
	float max_distance = 1000.f;

	vec3 compute_look_at() const {
		return position - glm::mat3_cast(orientation) * vec3(0, 0, distance);
	}

    void look_at(const vec3& dst, const vec3& src) {
        const vec3 FORWARD = vec3(0,0,-1);
        const vec3 UP = vec3(0,1,0);

        position = src;
        distance = math::distance(dst, src);

        vec3 fwd = math::normalize(dst - src);
        float dot = math::dot(FORWARD, fwd);

        if (math::abs(dot - (-1.0f)) < 0.000001f) {
            orientation = quat(math::PI, UP.x, UP.y, UP.z);
            return;
        }
        if (math::abs(dot - (1.0f)) < 0.000001f) {
            orientation = quat(1,0,0,0);
            return;
        }

        float angle = math::acos(dot);
        vec3 axis = math::normalize(math::cross(FORWARD, fwd));
        orientation = math::angle_axis(angle, axis);
    }

	void update();
};
