#pragma once

#include <core/vector_types.h>

struct ViewParam {
    struct {
        mat4 view;
        mat4 proj;
        mat4 view_proj;
        mat4 norm;

        struct {
            mat4 view;
            mat4 proj;
            mat4 view_proj;
            mat4 norm;
        } inverse;
    } matrix;

    vec2 resolution;
    vec2 jitter;

    struct {
        struct {
            mat4 view;
            mat4 proj;
            mat4 view_proj;
        } matrix;

        vec2 jitter;
    } previous;
};
