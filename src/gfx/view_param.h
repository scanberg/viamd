#pragma once

#include <core/md_vec_math.h>

struct ViewParam {
    struct Block {
        mat4_t view;
        mat4_t proj;
        mat4_t norm;
    };

    struct {
        Block curr;
        Block inv;
        Block prev;
    } matrix;

    struct {
        vec2_t curr;
        vec2_t prev;
    } jitter;

    struct {
        float near;
        float far;
    } clip_planes;

    vec2_t resolution;
    float fov_y;
};
