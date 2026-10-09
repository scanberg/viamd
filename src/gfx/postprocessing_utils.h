#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/md_vec_math.h>

namespace postprocessing {
// Legacy utility helpers used by callsites that manipulate the current bound target.

void blit_static_velocity(GLuint tex_depth, const ViewParam& view_param);

void scale_hsv(GLuint color_tex, vec3_t hsv_scale);

void blit_texture(GLuint tex);
// Writes the colour premultiplied by its alpha, (rgb * a, a): meant for the transparency buffer
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
    GLuint transparency = 0;        // LDR premultiplied colour, blended over the tone mapped image (overlays, selection, DVR)
    GLuint transparency_hdr = 0;    // optional: premultiplied linear radiance, blended over the HDR image before tone mapping (isosurfaces)
    bool   transparency_under_hdr = false;  // transparency lies under transparency_hdr (a selection tint under isosurfaces): it is covered by its alpha
    GLuint history = 0;         // TAA history target, written this frame
    GLuint history_prev = 0;    // optional: last frame's history (ping-pong). If 0, history is copied internally
};

struct Settings {
    vec3_t background_color = {20.f, 20.f, 20.f};

    struct {
        bool enabled = true;
        Tonemapper mode = Tonemapper_ACES;
        float exposure = 1.0f;
        float gamma = 2.4f;
    } tonemap;

    // Scale-free: there is no world-space radius. The sampling range is a fixed fraction of the viewport.
    struct {
        bool enabled = true;
        float intensity = 5.0f;
    } ssao;

    struct {
        bool enabled = true;
        float focus_depth = 0.5f;   // view-space distance to the focus plane
        float aperture = 0.01f;     // CoC of an object at infinity, as a fraction of the viewport height
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
void execute(const Inputs& in, const Settings& settings, const ViewParam& view);
void shutdown();

}  // namespace postprocess_pipeline