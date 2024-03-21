#version 150 core

#ifndef USE_TRANSPARENCY
#define USE_TRANSPARENCY 0
#endif

#ifndef USE_SSAO
#define USE_SSAO 0
#endif

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;
uniform sampler2D u_texture_transparent;
uniform sampler2D u_texture_ssao;

uniform mat4 u_inv_proj_mat;
uniform vec4 u_bg_color;
uniform float u_time;
uniform float u_gamma;

in vec2 tc;
out vec4 out_frag;

// Tonemapping operations (in case one writes to post tonemap)
float max3(float x, float y, float z) { return max(x, max(y, z)); }
float rcp(float x) { return 1.0 / x; }

const float exposure_bias = 0.5;

vec3 Tonemap(vec3 c) {
    c = c * exposure_bias;
    c = c * rcp(max3(c.r, c.g, c.b) + 1.0);
    c = pow(c, vec3(1.0 / u_gamma));
    return c;
}

vec3 TonemapInvert(vec3 c) {
    c = pow(c, vec3(u_gamma));
    c = c * rcp(1.0 - max3(c.r, c.g, c.b));
    c = c / exposure_bias;
    return c;
}

// TODO: Use linear depth instead and use uniform vec4 for unpacking to view coords.
vec4 depth_to_view_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

// https://aras-p.info/texts/CompactNormalStorage.html
vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec4 rand4(vec2 n) {
    return fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* vec4(43758.5453, 28001.8384, 50849.4141, 12996.89));
}

vec4 srand4(vec2 n) {
    return rand4(n) * 2.0 - 1.0;
}

const vec3 env_radiance = vec3(5.0);
const vec3 dir_radiance = vec3(10.0);
const vec3 L = normalize(vec3(1,1,1));
const float spec_exp = 100.0;

vec3 lambert(in vec3 radiance) {
    const float ONE_OVER_PI = 1.0 / 3.1415926535;
    return radiance * ONE_OVER_PI;
}

float fresnel(float H_dot_V) {
    const float n1 = 1.0;
    const float n2 = 1.5;
    const float R0 = pow((n1-n2)/(n1+n2), 2);
    return R0 + (1.0 - R0)*pow(1.0 - H_dot_V, 5);
}

vec3 shade(vec3 color, vec3 V, vec3 N) {
    vec3 H = normalize(L + V);
    float H_dot_V = clamp(dot(H, V), 0.0, 1.0);
    float N_dot_H = clamp(dot(N, H), 0.0, 1.0);
    float N_dot_L = clamp(dot(N, L), 0.0, 1.0);
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * lambert(env_radiance + N_dot_L * dir_radiance);
    vec3 specular = fr * (env_radiance + dir_radiance) * pow(N_dot_H, spec_exp);

    return diffuse + specular;
}

void main() {
    vec3 result = u_bg_color.rgb;
    
    vec4 color = texelFetch(u_texture_color, ivec2(gl_FragCoord.xy), 0);
    if (color.a != 0) {
        float depth = texelFetch(u_texture_depth, ivec2(gl_FragCoord.xy), 0).x;
        vec3 normal = decode_normal(texelFetch(u_texture_normal, ivec2(gl_FragCoord.xy), 0).xy);
        vec4 view_coord = depth_to_view_coord(tc, depth);

        // Add noise to reduce banding
        // Signed random to not affect overall luminance
        color.rgb = clamp(color.rgb + srand4(tc + u_time).rgb * 0.05, 0.0, 1.0);

        vec3 N = normal;
        vec3 V = -normalize(view_coord.xyz);
        color.rgb = shade(color.rgb, V, N);
        result = color.rgb;
    }

#if USE_SSAO
    float ssao = texelFetch(u_texture_ssao, ivec2(gl_FragCoord.xy), 0).x;
    result = result * ssao;
#endif

#if USE_TRANSPARENCY
    // Transparency layer is encoded as RGBA8 and thus not HDR
    // Therefore we tonemap, compose in LDR and invert tonemap back to HDR
    result = Tonemap(result);
    vec4 trans = texelFetch(u_texture_transparent, ivec2(gl_FragCoord.xy), 0);
    // The inverse Tonemap operation can result in NaNs if we have evil inputs which will mess up the temporal accumulation buffer and spread like wildfire
    // So we need to take some precautions here
    trans = clamp(trans, vec4(0.0), vec4(0.99, 0.99, 0.99, 1.0));
    result = trans.a * trans.rgb + result.rgb * (1.0 - trans.a);
    result = TonemapInvert(result);
#endif

    out_frag = vec4(result, 1);
}