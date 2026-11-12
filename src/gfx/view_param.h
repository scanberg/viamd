#pragma once

#include <core/md_vec.h>

struct ViewParam {
    struct Block {
        md_float4x4_t view;
        md_float4x4_t proj;
        md_float4x4_t norm;
        md_float4x4_t proj_no_jitter;
    };

    struct {
        Block curr;
        Block inv;
        Block prev;
    } matrix;

    struct {
        md_float2_t curr;
        md_float2_t prev;
    } jitter;

    struct {
        float near;
        float far;
    } clip_planes;

    md_float2_t resolution;
    float fov_y;
};
