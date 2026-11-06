#version 150 core

// The display transform: scene referred linear radiance -> display, and the composite over the background.
//
// Units: 1 is reference white, the radiance of a white Lambertian surface under irradiance pi (see
// postprocess_pipeline::lighting). The transform is that of a camera and a display, with nothing in it to
// compensate for the scale of the lights:
//   x       = exposure * L                              (exposure 1: reference white is 1)
//   display = curve(x) / curve(W)                       (W, the white point: the radiance shown as white)
//   out     = srgb_encode(display) ^ (2.2 / gamma)      (gamma 2.2 is the sRGB curve itself)
// The background is not part of the scene: it lights nothing and the curve does not apply to it. The scene arrives
// premultiplied by its coverage, a = 1 - (1 - geometry) (1 - isosurfaces), with nothing behind it, and is put over the
// background after the curve, in display referred linear values (as a render is composited over a plate):
//   display = a * curve(x / a) / curve(W) + (1 - a) * background
// so the background shows as picked whatever the curve, and a translucent surface over it is blended with it, not
// with a luminous wall bright enough to come out white.
// The curves:
//   0 none:    clip to [0, 1], no W (the basic shading mode, and tone mapping off)
//   1 neutral: Khronos PBR Neutral. Below its highlight compression the colours pass through as they are, but
//              for a small offset in the darkest ones (the 4% a dielectric reflects of a white environment),
//              so the base colours are reproduced; only the highlights are compressed, toward white.
//   2 ACES:    the fit of the ACES RRT + ODT by Stephen Hill: more contrast, shifts and desaturates colours
//   3 filmic:  John Hable's Uncharted 2 curve, with its own published input scale (2) and white point (11.2)

uniform sampler2D u_texture;            // premultiplied radiance of the scene
uniform sampler2D u_texture_coverage;   // r: coverage of the geometry (blurred by the depth of field when it is on)
uniform sampler2D u_texture_iso;        // a: coverage of the isosurfaces (premultiplied radiance, added after DOF)
uniform int   u_use_iso = 0;
uniform vec3  u_background = vec3(1.0); // display referred linear: what the encoding shows as the colour picked
uniform int   u_curve = 1;
uniform float u_exposure = 1.0;
uniform float u_gamma = 2.2;

const float WHITE_POINT = 8.0;          // three stops above reference white
const float HABLE_SCALE = 2.0;          // Hable's published input scale and white point
const float HABLE_WHITE_POINT = 11.2;

out vec4 out_frag;

vec3 neutral(vec3 color) {
    const float start_compression = 0.8 - 0.04;
    const float desaturation = 0.15;

    float x = min(color.r, min(color.g, color.b));
    float offset = x < 0.08 ? x - 6.25 * x * x : 0.04;
    color -= offset;

    float peak = max(color.r, max(color.g, color.b));
    if (peak < start_compression) return color;

    const float d = 1.0 - start_compression;
    float new_peak = 1.0 - d * d / (peak + d - start_compression);
    color *= new_peak / peak;

    float g = 1.0 - 1.0 / (desaturation * (peak - new_peak) + 1.0);
    return mix(color, vec3(new_peak), g);
}

// https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
// sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT (the vectors are the rows: v * M)
const mat3 ACES_INPUT = mat3(
    vec3(0.59719, 0.35458, 0.04823),
    vec3(0.07600, 0.90834, 0.01566),
    vec3(0.02840, 0.13383, 0.83777));
// ODT_SAT => XYZ => D60_2_D65 => sRGB
const mat3 ACES_OUTPUT = mat3(
    vec3( 1.60475, -0.53108, -0.07367),
    vec3(-0.10208,  1.10813, -0.00605),
    vec3(-0.00327, -0.07276,  1.07602));

vec3 aces(vec3 color) {
    vec3 v = color * ACES_INPUT;
    vec3 a = v * (v + 0.0245786) - 0.000090537;
    vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return (a / b) * ACES_OUTPUT;
}

vec3 hable(vec3 x) {
    const float A = 0.15;
    const float B = 0.50;
    const float C = 0.10;
    const float D = 0.20;
    const float E = 0.02;
    const float F = 0.30;
    return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

vec3 linear_to_srgb(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(vec3(0.0031308), c));
}

vec3 curve(vec3 x) {
    if (u_curve == 1) return neutral(x) / neutral(vec3(WHITE_POINT)).r;
    if (u_curve == 2) return aces(x) / aces(vec3(WHITE_POINT)).r;
    if (u_curve == 3) return hable(HABLE_SCALE * x) / hable(vec3(HABLE_WHITE_POINT)).r;
    return x;
}

void main() {
    ivec2 px = ivec2(gl_FragCoord.xy);
    vec3  L  = texelFetch(u_texture, px, 0).rgb;
    float a_geom = clamp(texelFetch(u_texture_coverage, px, 0).r, 0.0, 1.0);
    float a_iso  = (u_use_iso != 0) ? clamp(texelFetch(u_texture_iso, px, 0).a, 0.0, 1.0) : 0.0;
    float a = 1.0 - (1.0 - a_geom) * (1.0 - a_iso);

    vec3 scene = (a > 1.0e-4) ? curve(max(L * (u_exposure / a), vec3(0.0))) : vec3(0.0);
    vec3 color = clamp(scene, 0.0, 1.0) * a + u_background * (1.0 - a);

    color = pow(linear_to_srgb(color), vec3(2.2 / u_gamma));
    out_frag = vec4(color, 1.0);
}
