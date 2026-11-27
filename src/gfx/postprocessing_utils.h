#pragma once

#include <core/md_vec_math.h>

// THE LOOK of a view: how md_render_bridge.cpp shades what md_render drew (md_compose), and how md_post shows it (the
// display transform, the temporal antialiasing, the marks of the selection and the highlight). The main view takes it
// from the visuals (main.cpp's frame_settings), the components' views from the same (panel_view_settings).

namespace postprocess_pipeline {

// The curve of the display transform (md_post_display). Stored in workspaces by value: new ones go at the end.
enum Tonemapper {
    Tonemapper_Passthrough,     // no curve: clipped to the display
    Tonemapper_ExposureGamma,   // the same as Passthrough (exposure and gamma apply to every curve), kept for workspaces
    Tonemapper_Filmic,          // John Hable's Uncharted 2 curve
    Tonemapper_ACES,            // fitted ACES (RRT + ODT): contrast, but it shifts and desaturates colours
    Tonemapper_Neutral,         // Khronos PBR Neutral: the colours as they are, only the highlights compressed
    Tonemapper_Count
};

// How the colours of the geometry are lit
enum ShadingMode {
    // A simplified Cook-Torrance model of the colours (sRGB, decoded to linear albedo) under a key light and an
    // environment of sky and ground (lighting()), then tone mapped (Settings::tonemap)
    ShadingMode_PhysicallyBased,
    // The colours as given where a surface faces the light, darkened away from it, with an optional white
    // highlight (Settings::shading). Shown through the display transform without a curve or exposure, the exact
    // sRGB encoding (Settings::tonemap is not used), so a colour in full light appears as it was picked.
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

// Estimator of the ambient occlusion (md_compose's ao_mode)
enum SsaoMode {
    SsaoMode_Performance,
    SsaoMode_Quality,
    SsaoMode_Count
};

struct Settings {
    // The background (sRGB). Not part of the scene: it lights nothing and is not tone mapped. The scene is put over
    // it after the display transform, by its coverage, so it shows as picked in every mode and with every curve, and
    // a translucent surface over it is blended with it.
    md_float3_t background_color = {1.f, 1.f, 1.f};

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

    // Shadows of the directional light (both modes): the visibility of the light, as the ambient occlusion is the
    // visibility of the surroundings. The light is a disk: softness is its angular radius (degrees), which sets how
    // the shadows soften with distance.
    struct {
        bool enabled = true;
        float softness = 10.0f;
    } shadows;

    // The visibility of the surroundings, which darkens their light (physically based) or all of it (basic)
    struct {
        bool enabled = true;
        float intensity = 5.0f;
        SsaoMode mode = SsaoMode_Performance;
    } ssao;

    // The shadows and the ambient occlusion accumulated over frames, their sample patterns turning with the frame
    // (md_compose's history); else the patterns held still
    bool occlusion_temporal = true;

    // The temporal antialiasing (md_post_resolve): how much of the history is kept, from where the image moves to
    // where it is still
    struct {
        bool enabled = true;
        float feedback_min = 0.80f;
        float feedback_max = 0.95f;
        // The image shown smeared along the motion (md_post_resolve's motion_blur): of the frame's motion, how much
        struct {
            bool enabled = true;
            float motion_scale = 1.0f;
        } motion_blur;
    } taa;

    // The image shown sharpened after the TAA (md_post_present's): weight is its unsharp mask's, (1 + w) of a pixel
    // less w / 4 of its four neighbours. With the TAA only: the softness it is there for is the TAA's.
    struct {
        bool enabled = true;
        float weight = 0.5f;
    } sharpen;

    // The marks, put over the finished image as it is presented, in display colours: [0] of the selection, [1] of the
    // highlight. A tint where they are in view, a fill where they are behind something, and the outline of where they
    // are at all. A colour's alpha is how much of it. The highlight's group 0 has its colours, group k > 0 the colour
    // groups[k - 1] with their alphas.
    struct {
        md_float4_t visible[2] = {};
        md_float4_t hidden[2]  = {};
        md_float4_t outline[2] = {};
        md_float3_t groups[7]  = {};
        float highlight   = 1.0f;   // how much of the highlight shows: it fades in
        int   outline_step = 1;     // pixels between the taps that find the outline, which is about twice as wide
    } marks;
};

// The lights for a shading mode, in view space: they move with the camera. They do not depend on the background.
// Units: radiance in units of reference white, the radiance of a white Lambertian surface under irradiance pi
// (BRDF 1/pi, so its radiance is the irradiance over pi). The physically based lights are scaled so that the diffuse
// reflection of a white surface facing the viewer is 1: the display transform at exposure 1 then shows a colour on a
// surface facing the viewer close to as it was picked (with the Neutral curve, within about 0.02 of reference white:
// its offset takes off 0.04, the specular reflection of these surroundings adds about 0.02).
struct Lighting {
    md_float3_t light_dir;           // toward the key light, normalized
    md_float3_t light_irradiance;    // the key light: irradiance at normal incidence (pi: a white surface facing it is 1)
    md_float3_t sky_radiance;        // the environment above the horizon of the view (view space y > 0)
    md_float3_t ground_radiance;     // the environment below it
};

// Physically based: the key light from the upper right, 55 degrees off the view axis, which models the shapes, and
// light from all around, nearly uniform: a sky a little brighter than the ground, which keeps a hint of up and down.
// The side away from the key light keeps about 0.37 of the lit side. A white surface facing the viewer:
//   key: 0.87 * cos(55 deg) = 0.50, environment: (0.59 + 0.49) / 2 = 0.54
// which, less what the specular reflection takes off the diffuse (Fresnel, about 3.5%), is a diffuse reflection of 1:
// reference white. Facing the key light it is about 1.4. Settings::shading.light and .environment scale the two.
// Basic: from the upper right too, but 27 degrees off the view axis, so the surfaces facing the viewer keep close to
// their full colour. The radiances are not used.
static inline Lighting lighting(ShadingMode mode) {
    Lighting l = {
        .light_dir        = md_float3_normalize(md_float3_set(1.0f, 1.0f, 1.0f)),
        .light_irradiance = md_float3_set1(0.87f * (float)MD_PI),
        .sky_radiance     = md_float3_set1(0.59f),
        .ground_radiance  = md_float3_set1(0.49f),
    };
    if (mode == ShadingMode_Basic) {
        l.light_dir = md_float3_normalize(md_float3_set(0.3f, 0.4f, 1.0f));
    }
    return l;
}

}  // namespace postprocess_pipeline
