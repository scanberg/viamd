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

// Stored in workspaces by value: new ones go at the end
enum Tonemapper {
    Tonemapper_Passthrough,
    Tonemapper_ExposureGamma,
    Tonemapper_Filmic,      // Uncharted 2
    Tonemapper_ACES,        // fitted ACES (RRT + ODT): contrast, but it shifts and desaturates colours
    Tonemapper_Neutral,     // Khronos PBR Neutral: the colours as they are, only the highlights compressed
    Tonemapper_Count
};

// How the compose pass, the first step of the pipeline, lights the colours of the geometry
enum ShadingMode {
    // A simplified Cook-Torrance model of the colours (sRGB, decoded to linear albedo) under a key light and an
    // environment of sky and ground (lighting()), then tone mapped (Settings::tonemap)
    ShadingMode_PhysicallyBased,
    // The colours as given where a surface faces the light, darkened away from it, with an optional white
    // highlight (Settings::shading). Shown through the exact sRGB encoding rather than a tone mapper
    // (Settings::tonemap is not used), so a colour in full light appears as it was picked (to within the
    // precision of the HDR targets: one or two steps of 8 bits).
    ShadingMode_Basic,
    ShadingMode_Count
};

// Estimator of the ambient occlusion. Both feed the same temporal filter, blur and upsample.
enum SsaoMode {
    SsaoMode_Performance,   // point samples (shaders/ssao/ssao.frag)
    SsaoMode_Quality,       // screen-space slices with a visibility bitmask (shaders/ssao/ssao_bitmask.frag)
    SsaoMode_Count
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
    // The background as picked (sRGB) and its intensity. Physically based: the background is radiance, the
    // colour decoded to linear times the intensity, which is then tone mapped (24 is the white point of the
    // tone mappers). Basic: it shows as the colour itself and the intensity is not used.
    vec3_t background_color = {1.f, 1.f, 1.f};
    float  background_intensity = 24.f;

    struct {
        ShadingMode mode = ShadingMode_PhysicallyBased;
        // Basic only
        float ambient  = 0.3f;      // the fraction of the colour left on the side facing away from the light
        float specular = 0.25f;     // strength of the white highlight, 0 for none
    } shading;

    struct {
        bool enabled = true;
        Tonemapper mode = Tonemapper_Neutral;
        float exposure = 1.0f;
        float gamma = 2.2f;     // 2.2: the sRGB curve for Neutral
    } tonemap;

    // Scale-free: there is no world-space radius. The sampling range is a fixed fraction of the viewport.
    struct {
        bool enabled = true;
        float intensity = 5.0f;
        bool temporal = true;       // accumulate over frames (needs Inputs::velocity); off gives a fixed sample pattern
        SsaoMode mode = SsaoMode_Performance;
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

// GPU time of the ambient occlusion, averaged over a number of frames (0 until the first results are in)
struct SsaoTimings {
    float estimator_ms = 0.0f;  // the AO pass alone
    float total_ms = 0.0f;      // estimator, temporal filter, blur and upsample
};

void initialize(int width, int height);
void execute(const Inputs& in, const Settings& settings, const ViewParam& view);
void shutdown();
SsaoTimings ssao_timings();

// The lights of the compose pass for a shading mode, in view space: they move with the camera. The isosurfaces
// (volume::render_isosurfaces) are lit by the same. The radiances are those of the physically based mode, in the
// units the Neutral and ACES tone mappers expect (their exposure bias takes 4 to 1); they do not depend on the
// background.
struct Lighting {
    vec3_t light_dir;           // toward the key light, normalized
    vec3_t light_irradiance;    // the key light: irradiance at normal incidence
    vec3_t sky_radiance;        // the environment above the horizon of the view (view space y > 0)
    vec3_t ground_radiance;     // the environment below it
};
Lighting lighting(ShadingMode mode);

}  // namespace postprocess_pipeline