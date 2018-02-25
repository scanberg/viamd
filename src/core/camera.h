#pragma once

#include <core/types.h>

struct Camera {
	vec3 position {};
	quat orientation = quat(0,0,0,1);

    float near_plane = 1.f;
    float far_plane = 1000.f;
    float fov_y = 3.1415926534f / 4.f;
};

mat4 compute_world_to_view_matrix(const Camera& camera);
mat4 compute_view_to_world_matrix(const Camera& camera);
mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height);
mat4 compute_orthographic_projection_matrix(const Camera& camera, int width, int height);
void camera_look_at(Camera* camera, vec3 look_at, vec3 look_up);