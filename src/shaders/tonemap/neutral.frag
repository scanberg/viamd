#version 150 core

// Khronos PBR Neutral tone mapping (https://github.com/KhronosGroup/ToneMapping/tree/main/PBR_Neutral).
// Below the start of the highlight compression the colours pass through as they are, but for a small offset
// in the darkest ones (which takes off the 4% a dielectric reflects of a white environment), so the base
// colours of the geometry are reproduced. Only the highlights are compressed, and desaturated toward white
// as they approach it. No hue shifts, no toe that crushes dark colours.
// The curve only reaches 1 at infinity: as the ACES one does, it is scaled so the white point (24, the
// intensity of a white background) comes out white, which lifts everything else by about 1%.
// Then the sRGB encoding, with the gamma setting as an adjustment on top: at 2.2 it is the sRGB curve itself.

uniform sampler2D u_texture;
uniform float u_exposure = 1.0;
uniform float u_gamma = 2.2;

out vec4 out_frag;

vec3 pbr_neutral(vec3 color) {
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

vec3 linear_to_srgb(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(vec3(0.0031308), c));
}

void main() {
    vec4 in_frag = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);

    const float exposure_bias = 0.25;   // as the ACES tone mapper: 4 in the pipeline is 1 on the display
    const float white_point = 24.0;

    float white_scale = pbr_neutral(vec3(white_point * exposure_bias * u_exposure)).r;
    vec3 color = pbr_neutral(max(in_frag.rgb * exposure_bias * u_exposure, vec3(0.0))) / white_scale;
    color = pow(linear_to_srgb(color), vec3(2.2 / u_gamma));
    out_frag = vec4(color, in_frag.a);
}
