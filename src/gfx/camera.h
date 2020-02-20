#pragma once

#include <core/vector_types.h>

struct Camera {
    vec3 position = {0, 0, 0};
    quat orientation = quat(0, 0, 0, 1);

    float near_plane = 1.0f;
    float far_plane = 1000.0f;
    float fov_y = (3.1415926534f / 4.0f);
};
