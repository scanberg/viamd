#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/md_vec_math.h>

namespace postprocessing {
// Legacy utility helpers used by callsites that manipulate the current bound target.

void blit_static_velocity(GLuint tex_depth, const ViewParam& view_param);

void scale_hsv(GLuint color_tex, vec3_t hsv_scale);

void blit_texture(GLuint tex);
void blit_color(vec4_t color);

}  // namespace postprocessing

namespace postprocess_pipeline {

enum Tonemapper {
    Tonemapper_Passthrough,
    Tonemapper_ExposureGamma,
    Tonemapper_Filmic,
    Tonemapper_ACES,
};

struct Inputs {
    GLuint depth = 0;
    GLuint color = 0;
    GLuint normal = 0;
    GLuint velocity = 0;
    GLuint transparency = 0;
    GLuint history = 0;
};

struct Settings {
    vec3_t background_color = {20.f, 20.f, 20.f};

    struct {
        bool enabled = true;
        Tonemapper mode = Tonemapper_ACES;
        float exposure = 1.0f;
        float gamma = 2.4f;
    } tonemap;

    struct {
        bool enabled = true;
        float radius = 6.0f;
        float intensity = 3.0f;
        float bias = 0.1f;
    } ssao;

    struct {
        bool enabled = true;
        float focus_depth = 0.5f;
        float focus_scale = 10.f;
    } dof;

    struct {
        bool enabled = true;
    } fxaa;

    struct {
        bool enabled = true;
        float feedback_min = 0.80f;
        float feedback_max = 0.95f;
        struct {
            bool enabled = true;
            float motion_scale = 0.5f;
        } motion_blur;
    } taa;

    struct {
        bool enabled = true;
        float weight = 1.0f;
    } sharpen;
};

void initialize(int width, int height);
// Depth is recorded separately to allow for depth buffer reuse.
// This will internally store a linear depth representation that is required later in execute.
void record_depth(GLuint in_depth, const ViewParam& view);
void execute(const Inputs& in, const Settings& settings, const ViewParam& view);
void shutdown();

}  // namespace postprocess_pipeline