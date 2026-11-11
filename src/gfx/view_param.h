#pragma once

#include <core/md_vec.h>

struct ViewParam {
    struct Block {
        md_f32x4x4_t view;
        md_f32x4x4_t proj;
        md_f32x4x4_t norm;
        md_f32x4x4_t proj_no_jitter;
    };

    struct {
        Block curr;
        Block inv;
        Block prev;
    } matrix;

    struct {
        md_f32x2_t curr;
        md_f32x2_t prev;
    } jitter;

    struct {
        float near;
        float far;
    } clip_planes;

    md_f32x2_t resolution;
    float fov_y;
};
