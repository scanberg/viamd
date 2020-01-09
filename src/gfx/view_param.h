#pragma once

#include <core/vector_types.h>

struct ViewParam {
    struct Block {
        mat4 view;
        mat4 proj;
        mat4 proj_jittered;
        mat4 view_proj;
        mat4 view_proj_jittered;
        mat4 norm;
    };

    struct {
        Block current;
        Block inverse;
        Block previous;
    } matrix;

    struct {
        vec2 current;
        vec2 previous;
    } jitter;

    struct {
        float near;
        float far;
    } clip_planes;

    float fov_y;
    vec2 resolution;
};
