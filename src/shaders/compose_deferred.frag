#version 150 core

layout (std140) uniform UniformData {
    vec4  u_proj_info;
    float u_linear_colors;      // 1: the colours are linear reflectances (ColorSpace_Linear), 0: sRGB
    float u_light_strength;     // of the directional light, relative to u_dir_radiance (and to the basic light)
    float u_use_shadow;         // 1: u_texture_shadow holds the visibility of the key light
    float u_time;
    vec3  u_sky_radiance;       // the environment above the horizon of the view (y > 0)
    float u_roughness;
    vec3  u_dir_radiance;       // the key light: irradiance at normal incidence
    float u_F0;
    vec3  u_light_dir;          // toward the key light
    float u_basic;              // 1: shade_basic (postprocess_pipeline::ShadingMode_Basic), 0: shade
    vec3  u_ground_radiance;    // the environment below the horizon
    float u_basic_ambient;
    float u_basic_specular;
    float u_use_ao;             // 1: u_texture_ao holds the visibility of the ambient occlusion
    vec2  u_inv_viewport;       // 1 / the size of the viewport, which can be smaller than the targets
};

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;
uniform sampler2D u_texture_ao;        // the visibility of the surroundings (ambient occlusion)
uniform sampler2D u_texture_shadow;    // the visibility of the key light (shadows)

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

// The colours are picked, and stored, as sRGB, or are linear reflectances already (u_linear_colors)
vec3 srgb_to_linear(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c));
}

// -----------------------------------------------------------------------------
// Physically based shading: a simplified Cook-Torrance model (GGX, Smith-Schlick, Schlick Fresnel) of a
// dielectric whose albedo is the colour, decoded from sRGB to linear. Lit by a key light and an environment
// of two uniform hemispheres split at the horizon of the view, the sky above and the ground below:
//  - the key light: diffuse (1 - F) albedo / pi * E N.L, specular D G F / (4 N.V N.L) * E N.L
//  - the environment, diffuse: the exact irradiance of the two hemispheres, pi * mix(ground, sky, (1 + N.y) / 2)
//  - the environment, specular: the radiance about the reflected direction, the horizon blurred by the
//    roughness, times the directional albedo of the specular lobe (the split sum, from Karis' analytic fit).
//    That albedo is what the surface reflects of its surroundings, so it also takes its share off the diffuse.
// The environment does not depend on the background, which is not part of the scene: a dark one does not leave
// the shadowed sides black, and a coloured one does not tint the geometry.
// -----------------------------------------------------------------------------

float FresnelSchlick(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
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

// The directional albedo of the specular lobe, F0 * x + y: Karis, "Physically Based Shading on Mobile" (2014)
vec2 EnvBRDFApprox(float NdotV, float roughness) {
    const vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
    const vec4 c1 = vec4( 1.0,  0.0425,  1.04, -0.04);
    vec4  r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    return vec2(-1.04, 1.04) * a004 + r.zw;
}

// The irradiance of the environment over pi
vec3 environment_irradiance(vec3 N) {
    return mix(u_ground_radiance, u_sky_radiance, 0.5 + 0.5 * N.y);
}

// The radiance of the environment about the direction R, as seen through a lobe of the roughness
vec3 environment_radiance(vec3 R, float roughness) {
    float w = 0.1 + roughness;
    return mix(u_ground_radiance, u_sky_radiance, smoothstep(-w, w, R.y));
}

// Each light has its own visibility. The ambient occlusion is the visibility of the environment, so it applies to the
// environment's light only: to its diffuse part as it is, and to its specular part as the visibility of the specular
// lobe (Lagarde and de Rousiers, "Moving Frostbite to PBR", 2014). The shadows are the visibility of the key light, and
// apply to it, diffuse and specular.
float specular_occlusion(float N_dot_V, float ao, float roughness) {
    return clamp(pow(N_dot_V + ao, exp2(-16.0 * roughness - 1.0)) - 1.0 + ao, 0.0, 1.0);    // roughness: alpha
}

vec3 shade(vec3 color, vec3 V, vec3 N, float ao, float shadow) {
    vec3  L = normalize(u_light_dir);
    vec3  H = normalize(L + V);
    float H_dot_V = clamp(dot(H, V), 0.0, 1.0);
    float N_dot_H = clamp(dot(N, H), 0.0, 1.0);
    float N_dot_L = clamp(dot(N, L), 0.0, 1.0);
    float N_dot_V = clamp(dot(N, V), 0.0, 1.0);

    float F0 = clamp(u_F0, 0.0, 1.0);
    float roughness = clamp(u_roughness, 0.04, 1.0);
    vec3  albedo = (u_linear_colors > 0.5) ? clamp(color, 0.0, 1.0) : srgb_to_linear(color);

    vec3 Lo = vec3(0);

    // Key light
    if (N_dot_L > 0.0) {
        float NDF = DistributionGGX(N_dot_H, roughness);
        float G   = GeometrySmith(N_dot_V, N_dot_L, roughness);
        float F   = FresnelSchlick(H_dot_V, F0);
        float numerator   = NDF * G * F;
        float denominator = 4.0 * N_dot_V * N_dot_L + 0.0001;
        vec3 specular     = vec3(numerator / denominator);
        Lo += ((1.0 - F) * albedo * ONE_OVER_PI + specular) * (u_dir_radiance * u_light_strength) * (N_dot_L * shadow);
    }

    // Environment
    {
        vec2  ab = EnvBRDFApprox(N_dot_V, roughness);
        float Fe = F0 * ab.x + ab.y;
        vec3  R  = reflect(-V, N);
        Lo += (1.0 - Fe) * albedo * environment_irradiance(N) * ao + Fe * environment_radiance(R, roughness) * specular_occlusion(N_dot_V, ao, roughness * roughness);
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

vec3 shade_basic(vec3 color, vec3 V, vec3 N, float shadow) {
    vec3  L = normalize(u_light_dir);
    float w = 0.5 + 0.5 * dot(N, L);
    // The light strength scales what the light adds over the ambient, and its highlight; the shadow takes both away
    float light = u_light_strength * shadow;
    vec3  c = clamp(color, 0.0, 1.0) * (u_basic_ambient + light * (1.0 - u_basic_ambient) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));     // brighter than the colour (light above 1): keep its hue, not clip per channel
    vec3  H = normalize(L + V);
    float s = light * u_basic_specular * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = mix(c, vec3(1.0), clamp(s, 0.0, 1.0));
    return (u_linear_colors > 0.5) ? c : srgb_to_linear(c);     // linear colours are lit as they are
}

void main() {
    // Nothing behind the geometry: the display transform puts the background there, by the coverage
    vec3 result = vec3(0.0);

    ivec2 pixel = ivec2(gl_FragCoord.xy);
    vec4 color = texelFetch(u_texture_color, pixel, 0);
    if (color.a > 0) {
        float linear_depth = texelFetch(u_texture_depth, pixel, 0).x;
        vec3 normal = decode_normal(texelFetch(u_texture_normal, pixel, 0).xy);

        vec2 tex_coord = gl_FragCoord.xy * u_inv_viewport;    // gl_FragCoord is at the pixel centre
        vec3 view_coord = linear_depth_to_view_coord(tex_coord, linear_depth);

        vec3 N = normal;
#if COMPOSE_PROJECTION_PERSPECTIVE
        vec3 V = -normalize(view_coord);
#else
        vec3 V = vec3(0, 0, 1);     // every ray runs along the view axis
#endif
        float ao = (u_use_ao > 0.5) ? texelFetch(u_texture_ao, pixel, 0).r : 1.0;
        float shadow = (u_use_shadow > 0.5) ? texelFetch(u_texture_shadow, pixel, 0).r : 1.0;
        // Basic: the occlusion darkens all of it, as the depth cue it is there
        color.rgb = (u_basic > 0.5) ? shade_basic(color.rgb, V, N, shadow) * ao : shade(color.rgb, V, N, ao, shadow);
        result = color.rgb;
    }

    out_frag = vec4(result, 1);
}