#pragma once

#include <core/md_vec_math.h>

struct Camera {
    quat_t orientation = {0, 0, 0, 1};
    vec3_t position = {0, 0, 0};

    float focus_distance = 10.0f;
    float near_plane = 1.0f;
    float far_plane = 10000.0f;
    float fov_y = (3.1415926534f / 4.0f);
};
