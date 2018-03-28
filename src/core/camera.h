#pragma once

#include <core/types.h>

struct Camera {
	vec3 position {};
	quat orientation = quat(0,0,0,1);

    float near_plane = 1.f;
    float far_plane = 1000.f;
    float fov_y = 3.1415926534f / 4.f;
};