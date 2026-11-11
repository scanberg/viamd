#pragma once

#include <core/md_vec.h>

struct ViewTransform {
    md_quat_t orientation = {0, 0, 0, 1};
    md_f32x3_t position = {0, 0, 0};
	float  distance = 10.0f;
};

struct Camera : ViewTransform {
    float near_plane = 1.0f;
    float far_plane = 10000.0f;
    float fov_y = (3.1415926534f / 4.0f);

    Camera& operator = (const ViewTransform& t) {
        this->orientation = t.orientation;
        this->position = t.position;
        this->distance = t.distance;
        return *this;
	}
};
