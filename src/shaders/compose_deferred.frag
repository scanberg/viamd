#version 150 core

layout (std140) uniform UniformData {
    vec4  u_proj_info;
    vec3  u_bg_color;
    float u_time;
    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;
    vec3  u_light_dir;
    float u_basic;              // 1: shade_basic (postprocess_pipeline::ShadingMode_Basic), 0: shade
    float u_basic_ambient;
    float u_basic_specular;
};

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;

in vec2 tc;
out vec4 out_frag;

#ifndef COMPOSE_PROJECTION_PERSPECTIVE
#define COMPOSE_PROJECTION_PERSPECTIVE 1
#endif

vec3 linear_depth_to_view_coord(vec2 tex_coord, float linear_depth) {
    // proj_info uses tex_coord in [0,1]. For perspective, XY should be scaled
    // by positive eye-space distance while view-space Z stays negative in front
    // of the camera (OpenGL view convention).
    float eye_z = linear_depth;
    float view_z = -eye_z;
#if COMPOSE_PROJECTION_PERSPECTIVE
    return vec3((tex_coord * u_proj_info.xy + u_proj_info.zw) * eye_z, view_z);
#else
    return vec3((tex_coord * u_proj_info.xy + u_proj_info.zw), view_z);
#endif
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

const float PI = 3.1415926535;
const float ONE_OVER_PI = 1.0 / 3.1415926535;

float FresnelSchlick(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float FresnelSchlickRoughness(float cosTheta, float F0, float roughness) {
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float DistributionGGX(float NdotH, float roughness) {
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH2 = NdotH*NdotH;
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;
    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
    return num / denom;
}

float GeometrySmith(float NdotV, float NdotL, float roughness) {
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
    return ggx1 * ggx2;
}

vec3 shade(vec3 color, vec3 V, vec3 N) {
    vec3  L = normalize(u_light_dir);
    vec3  H = normalize(L + V);
    float H_dot_V = clamp(dot(H, V), 0.0, 1.0);
    float N_dot_H = clamp(dot(N, H), 0.0, 1.0);
    float N_dot_L = clamp(dot(N, L), 0.0, 1.0);
    float N_dot_V = clamp(dot(N, V), 0.0, 1.0);

    float F0 = clamp(u_F0, 0.0, 1.0);
    float roughness = clamp(u_roughness, 0.04, 1.0);
    vec3  albedo = color;

    vec3 Lo = vec3(0);

    // Directional light
    if (N_dot_L > 0.0) {
        float NDF = DistributionGGX(N_dot_H, roughness);
        float G   = GeometrySmith(N_dot_V, N_dot_L, roughness);
        float F   = FresnelSchlick(H_dot_V, F0);
        float kS = F;
        float kD = 1.0 - kS;
        float numerator   = NDF * G * F;
        float denominator = 4.0 * N_dot_V * N_dot_L + 0.0001;
        vec3 specular     = vec3(numerator / denominator);
        Lo += (kD * albedo * ONE_OVER_PI + specular) * u_dir_radiance * N_dot_L;
    }

    // Environment (still simplified; see note above about proper IBL)
    {
        float F = FresnelSchlickRoughness(N_dot_V, F0, roughness);
        float kS = F;
        float kD = 1.0 - kS;
        Lo += (kD * albedo * ONE_OVER_PI + vec3(kS)) * u_env_radiance;
    }

    return Lo;
}

// Basic shading: keeps the colours as given. Where a surface faces the light it shows its colour exactly;
// away from the light the colour darkens to u_basic_ambient of itself, along a wrapped ("half") Lambert
// falloff, (0.5 + 0.5 N.L)^2, which has no hard terminator and keeps some shape on the far side. The
// highlight (Blinn-Phong) blends toward white rather than adding to it, so no channel is pushed out of range
// and the hue stays. The light factor applies to the display (sRGB) values the colours were picked as, the
// way classic molecular viewers light them, and the result is returned as linear radiance: the display
// transform of the basic mode is the exact sRGB encoding, which takes it back to those values.
const float BASIC_SHININESS = 48.0;

vec3 srgb_to_linear(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c));
}

vec3 shade_basic(vec3 color, vec3 V, vec3 N) {
    vec3  L = normalize(u_light_dir);
    float w = 0.5 + 0.5 * dot(N, L);
    vec3  c = clamp(color, 0.0, 1.0) * mix(u_basic_ambient, 1.0, w * w);
    vec3  H = normalize(L + V);
    float s = u_basic_specular * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = mix(c, vec3(1.0), clamp(s, 0.0, 1.0));
    return srgb_to_linear(c);
}

void main() {
    vec3 result = u_bg_color.rgb;

    ivec2 pixel = ivec2(gl_FragCoord.xy);
    vec4 color = texelFetch(u_texture_color, pixel, 0);
    if (color.a > 0) {
        float linear_depth = texelFetch(u_texture_depth, pixel, 0).x;
        vec3 normal = decode_normal(texelFetch(u_texture_normal, pixel, 0).xy);

        vec2 tex_coord = (gl_FragCoord.xy + vec2(0.5)) / vec2(textureSize(u_texture_depth, 0));
        vec3 view_coord = linear_depth_to_view_coord(tex_coord, linear_depth);

        vec3 N = normal;
#if COMPOSE_PROJECTION_PERSPECTIVE
        vec3 V = -normalize(view_coord);
#else
        vec3 V = vec3(0, 0, 1);     // every ray runs along the view axis
#endif
        color.rgb = (u_basic > 0.5) ? shade_basic(color.rgb, V, N) : shade(color.rgb, V, N);
        result = color.rgb;
    }

    out_frag = vec4(result, 1);
}