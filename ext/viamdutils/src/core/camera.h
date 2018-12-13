#pragma once

#include <core/types.h>
#include <core/vector_types.h>

struct Camera {
    vec3 position = {0, 0, 0};
    quat orientation = quat(0, 0, 0, 1);

    float32 near_plane = 1.f;
    float32 far_plane = 1000.f;
    float32 fov_y = 3.1415926534f / 4.f;
};
