#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/md_vec_math.h>

namespace postprocessing {
// Legacy utility helpers used by callsites that manipulate the current bound target.

void blit_static_velocity(GLuint tex_depth, const ViewParam& view_param);

void blit_texture(GLuint tex);

}  // namespace postprocessing

namespace postprocess_pipeline {

// The curve of the display transform (shaders/tonemap/display.frag). Stored in workspaces by value: new ones go at
// the end.
enum Tonemapper {
    Tonemapper_Passthrough,     // no curve: clipped to the display
    Tonemapper_ExposureGamma,   // the same as Passthrough (exposure and gamma apply to every curve), kept for workspaces
    Tonemapper_Filmic,          // John Hable's Uncharted 2 curve
    Tonemapper_ACES,            // fitted ACES (RRT + ODT): contrast, but it shifts and desaturates colours
    Tonemapper_Neutral,         // Khronos PBR Neutral: the colours as they are, only the highlights compressed
    Tonemapper_Count
};

// How the compose pass, the first step of the pipeline, lights the colours of the geometry
enum ShadingMode {
    // A simplified Cook-Torrance model of the colours (sRGB, decoded to linear albedo) under a key light and an
    // environment of sky and ground (lighting()), then tone mapped (Settings::tonemap)
    ShadingMode_PhysicallyBased,
    // The colours as given where a surface faces the light, darkened away from it, with an optional white
    // highlight (Settings::shading). Shown through the display transform without a curve or exposure, the exact
    // sRGB encoding (Settings::tonemap is not used), so a colour in full light appears as it was picked (to within
    // the precision of the HDR targets: one or two steps of 8 bits).
    ShadingMode_Basic,
    ShadingMode_Count
};

// How the colours of the geometry (and of the isosurfaces) are read. Stored in workspaces by value. Linear by
// default, the way the colours were read before.
enum ColorSpace {
    // As picked and seen: the sRGB colours the pickers show, decoded to linear for the lighting. A colour in full
    // light shows as it was picked.
    ColorSpace_sRGB,
    // As the fraction of the light reflected per channel, used directly. Lighter and paler on screen than the swatch
    // the picker shows (0.5 reflects half the light, which the display shows as 0.735).
    ColorSpace_Linear,
    ColorSpace_Count
};

// Estimator of the ambient occlusion. Both feed the same temporal filter, blur and upsample.
enum SsaoMode {
    SsaoMode_Performance,   // point samples (shaders/ssao/ssao.frag)
    SsaoMode_Quality,       // screen-space slices with a visibility bitmask (shaders/ssao/ssao_bitmask.frag)
    SsaoMode_Count
};

struct Inputs {
    // The G-buffer's size, when it is smaller than the image made (the viewport's size): upscaled to it, by the TAA
    // where it is on, else a frame on its own. 0: the viewport's size.
    int    width = 0;
    int    height = 0;
    GLuint depth = 0;
    GLuint color = 0;
    GLuint normal = 0;
    GLuint velocity = 0;
    GLuint transparency = 0;        // LDR premultiplied colour, blended over the tone mapped image (overlays, DVR)
    GLuint transparency_hdr = 0;    // optional: premultiplied linear radiance, blended over the HDR image before tone mapping (isosurfaces); its alpha is part of the coverage that puts the scene over the background
    GLuint history = 0;         // TAA history target, written this frame, at the size of the image made
    GLuint history_prev = 0;    // optional: last frame's history (ping-pong). If 0, history is copied internally
    // optional: the marks of the selection and the highlight at the size of the image made, put over it as it is presented
    // (Settings::marks), filtered linearly: the selection in view (r), there at all (g), the highlight the same (b, a)
    GLuint marks = 0;
    // with marks: the group of the highlight (Settings::marks.groups) at the same size, R8 in units of 1/255: 1 + the
    // group, 8 more where it is in view, 0 none
    GLuint marks_group = 0;
};

struct Settings {
    // The background (sRGB). Not part of the scene: it lights nothing and is not tone mapped. The scene is put over
    // it after the display transform, by its coverage (the geometry, through the depth of field, and the
    // isosurfaces), so it shows as picked in every mode and with every curve, and a translucent surface over it is
    // blended with it.
    vec3_t background_color = {1.f, 1.f, 1.f};

    struct {
        ShadingMode mode = ShadingMode_PhysicallyBased;
        ColorSpace colors = ColorSpace_Linear;  // how the colours of the geometry are read, in either mode
        float light = 1.0f;         // strength of the directional light, relative to lighting(): in either mode
        float environment = 1.0f;   // strength of the light from the surroundings (sky and ground), relative to lighting(): physically based only
        // Basic only
        float ambient  = 0.3f;      // the fraction of the colour left on the side facing away from the light
        float specular = 0.25f;     // strength of the white highlight, 0 for none
    } shading;

    // The display transform: x = exposure * radiance (in units of reference white, see lighting()), shown as
    // curve(x) / curve(W), W the white point: 8, three stops above reference white (Hable's curve: its own, 5.6);
    // then the sRGB encoding raised to 2.2 / gamma (at 2.2 the sRGB curve itself). Off: no curve, exposure 1 and
    // gamma 2.2.
    struct {
        bool enabled = true;
        Tonemapper mode = Tonemapper_ACES;
        float exposure = 1.0f;
        float gamma = 2.2f;
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

    // The marks (Inputs::marks), put over the finished image as it is presented, in display colours: [0] of the
    // selection, [1] of the highlight. A tint where they are in view, a fill where they are behind something, and the
    // outline of where they are at all. A colour's alpha is how much of it. The highlight's group 0 has its colours,
    // group k > 0 the colour groups[k - 1] with their alphas.
    struct {
        vec4_t visible[2] = {};
        vec4_t hidden[2]  = {};
        vec4_t outline[2] = {};
        vec3_t groups[7]  = {};
        float highlight   = 1.0f;   // how much of the highlight shows: it fades in
        int   outline_step = 1;     // pixels between the taps that find the outline, which is about twice as wide
    } marks;
};


// width / height: the G-buffer's size; out_width / out_height: the image made, when it is upscaled to it (0: the same)
void initialize(int width, int height, int out_width = 0, int out_height = 0);
void execute(const Inputs& in, const Settings& settings, const ViewParam& view);
void shutdown();

// The lights of the compose pass for a shading mode, in view space: they move with the camera. The isosurfaces
// (volume::render_isosurfaces) are lit by the same. They do not depend on the background.
// Units: radiance in units of reference white, the radiance of a white Lambertian surface under irradiance pi
// (BRDF 1/pi, so its radiance is the irradiance over pi). The physically based lights are scaled so that the diffuse
// reflection of a white surface facing the viewer is 1: the display transform at exposure 1 then shows a colour on a
// surface facing the viewer close to as it was picked (with the Neutral curve, within about 0.02 of reference white:
// its offset takes off 0.04, the specular reflection of these surroundings adds about 0.02).
struct Lighting {
    vec3_t light_dir;           // toward the key light, normalized
    vec3_t light_irradiance;    // the key light: irradiance at normal incidence (pi: a white surface facing it is 1)
    vec3_t sky_radiance;        // the environment above the horizon of the view (view space y > 0)
    vec3_t ground_radiance;     // the environment below it
};
Lighting lighting(ShadingMode mode);

}  // namespace postprocess_pipeline