#version 410 core

#ifndef MAX_ISOVALUE_COUNT
#define MAX_ISOVALUE_COUNT 8
#endif

struct IsovalueParameters {
    float values[MAX_ISOVALUE_COUNT];
    vec4  colors[MAX_ISOVALUE_COUNT];
    int   count;
};

layout(std140) uniform UniformData {
    mat4 u_view_to_model_mat;
    mat4 u_model_to_view_mat;
    mat4 u_inv_proj_mat;
    mat4 u_model_view_proj_mat;

    vec2  u_inv_res;
    float u_time;
    float u_iso_optical_density_scale;

    vec3  u_clip_plane_min;
    float u_tf_min;
    vec3  u_clip_plane_max;
    float u_tf_inv_ext;

    vec3  u_gradient_spacing_world_space;
    float u_exposure;

    mat4  u_gradient_spacing_tex_space;

    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;

    float u_gamma;
};

uniform IsovalueParameters u_iso;

uniform sampler2D u_tex_entry;
uniform sampler2D u_tex_exit;
uniform sampler3D u_tex_density_volume;
uniform sampler3D u_tex_color_volume;
uniform sampler2D u_tex_tf;

layout(location = 0) out vec4 out_color;

const float REF_SAMPLING_RATE = 150.0;
const float ERT_THRESHOLD     = 0.995;
const float SAMPLING_RATE     = 2.0;
const float PI                = 3.1415926535;
const float INV_PI            = 1.0 / PI;
const float EPSILON           = 1e-6;

const vec3 LIGHT_DIR = normalize(vec3(1.0, 1.0, 1.0));

struct RayState {
    vec3 radiance;   // premultiplied accumulated radiance in output space
    float opacity;   // accumulated opacity = 1 - transmittance
    uint insideMask; // which isosurfaces currently enclose the ray
};

// -----------------------------------------------------------------------------
// Tone mapping
// -----------------------------------------------------------------------------

const mat3 ACES_INPUT_MAT = mat3(
    vec3(0.59719, 0.35458, 0.04823),
    vec3(0.07600, 0.90834, 0.01566),
    vec3(0.02840, 0.13383, 0.83777)
);

const mat3 ACES_OUTPUT_MAT = mat3(
    vec3( 1.60475, -0.53108, -0.07367),
    vec3(-0.10208,  1.10813, -0.00605),
    vec3(-0.00327, -0.07276,  1.07602)
);

vec3 RRTAndODTFit(vec3 v) {
    vec3 a = v * (v + 0.0245786) - 0.000090537;
    vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return a / b;
}

vec3 ACESFitted(vec3 color) {
    color = color * ACES_INPUT_MAT;
    color = RRTAndODTFit(color);
    color = color * ACES_OUTPUT_MAT;
    return clamp(color, 0.0, 1.0);
}

vec3 tonemap(vec3 hdr) {
    const float exposure_bias = 0.25;
    const float white_point   = 24.0;

    vec3 exposed = hdr * (exposure_bias * u_exposure);
    vec3 white_scale = ACESFitted(vec3(white_point * exposure_bias * u_exposure));

    vec3 color = ACESFitted(exposed) / white_scale;
    color = clamp(color, 0.0, 1.0);

    color = pow(color, vec3(1.0 / u_gamma)); // gamma correction
    return color;
}

// -----------------------------------------------------------------------------
// Sampling
// -----------------------------------------------------------------------------

float sampleDensity(vec3 texPos) {
    return texture(u_tex_density_volume, texPos).r;
}

vec4 sampleTransfer(float density) {
    float t = clamp((density - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
    return texture(u_tex_tf, vec2(t, 0.5));
}

vec4 sampleIsoColor(vec3 texPos, vec4 baseColor) {
#if defined(USE_COLOR_VOLUME)
    return baseColor * texture(u_tex_color_volume, texPos);
#else
    return baseColor;
#endif
}

vec3 sampleGradient(vec3 texPos) {
    vec3 dx = u_gradient_spacing_tex_space[0].xyz;
    vec3 dy = u_gradient_spacing_tex_space[1].xyz;
    vec3 dz = u_gradient_spacing_tex_space[2].xyz;

    vec3 g = vec3(
        sampleDensity(texPos + dx) - sampleDensity(texPos - dx),
        sampleDensity(texPos + dy) - sampleDensity(texPos - dy),
        sampleDensity(texPos + dz) - sampleDensity(texPos - dz)
    );

    return g / (2.0 * u_gradient_spacing_world_space);
}

float segmentLengthWorld(vec3 p0Tex, vec3 p1Tex) {
    vec3 deltaView = (u_model_to_view_mat * vec4(p1Tex - p0Tex, 0.0)).xyz;
    return length(deltaView);
}

// -----------------------------------------------------------------------------
// Ray accumulation
// -----------------------------------------------------------------------------

float remainingTransmittance(RayState state) {
    return 1.0 - state.opacity;
}

void accumulateEvent(inout RayState state, vec3 radiance, float opacity) {
    float T = remainingTransmittance(state);
    float a = clamp(opacity, 0.0, 1.0);

    state.radiance += T * radiance;
    state.opacity  += T * a;
}

float correctedOpacity(float alpha, float segmentLengthTex) {
    float a = clamp(alpha, 0.0, 1.0);
    float w = max(segmentLengthTex, 0.0) * REF_SAMPLING_RATE;
    return 1.0 - pow(max(1.0 - a, 1e-6), w);
}

void accumulateDvrSegment(inout RayState state, vec3 p0Tex, vec3 p1Tex) {
#if defined(INCLUDE_DVR)
    float lenTex = distance(p0Tex, p1Tex);
    if (lenTex <= EPSILON) {
        return;
    }

    vec3 midPos = 0.5 * (p0Tex + p1Tex);
    vec4 tf = sampleTransfer(sampleDensity(midPos));
    float a = correctedOpacity(tf.a, lenTex);

    accumulateEvent(state, tf.rgb * a, a);
#endif
}

void accumulateIsoInteriorAbsorption(inout RayState state, vec3 p0Tex, vec3 p1Tex) {
#if defined(INCLUDE_ISO)
    if (u_iso_optical_density_scale <= 0.0) {
        return;
    }

    int insideCount = bitCount(state.insideMask);
    if (insideCount <= 0) {
        return;
    }

    float lenWorld = segmentLengthWorld(p0Tex, p1Tex);
    if (lenWorld <= EPSILON) {
        return;
    }

    float tau = u_iso_optical_density_scale * float(insideCount) * lenWorld * REF_SAMPLING_RATE;
    float a = 1.0 - exp(-max(tau, 0.0));

    accumulateEvent(state, vec3(0.0), a);
#endif
}

void accumulateMediumSegment(inout RayState state, vec3 p0Tex, vec3 p1Tex) {
    accumulateDvrSegment(state, p0Tex, p1Tex);
    accumulateIsoInteriorAbsorption(state, p0Tex, p1Tex);
}

// -----------------------------------------------------------------------------
// Surface shading
// -----------------------------------------------------------------------------

float FresnelSchlick(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float FresnelSchlickRoughness(float cosTheta, float F0, float roughness) {
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float DistributionGGX(float NdotH, float roughness) {
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH2 = NdotH * NdotH;

    float denom = NdotH2 * (a2 - 1.0) + 1.0;
    return a2 / max(PI * denom * denom, 1e-6);
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = roughness + 1.0;
    float k = (r * r) / 8.0;
    return NdotV / max(NdotV * (1.0 - k) + k, 1e-6);
}

float GeometrySmith(float NdotV, float NdotL, float roughness) {
    return GeometrySchlickGGX(NdotV, roughness) * GeometrySchlickGGX(NdotL, roughness);
}

void accumulateSurfaceHit(inout RayState state, vec3 isoPosTex, vec4 baseColor) {
    vec3 V = -normalize((u_model_to_view_mat * vec4(isoPosTex, 1.0)).xyz);

    vec3 N = sampleGradient(isoPosTex);
    float nLen = length(N);
    N = (nLen > EPSILON) ? (N / nLen) : V;

    if (dot(N, V) < 0.0) {
        N = -N;
    }

    float F0 = clamp(u_F0, 0.0, 1.0);
    float roughness = clamp(u_roughness, 0.04, 1.0);
    vec3 albedo = baseColor.rgb;
    float baseAlpha = clamp(baseColor.a, 0.0, 1.0);

    float NdotL = clamp(dot(N, LIGHT_DIR), 0.0, 1.0);
    float NdotV = clamp(dot(N, V), 0.0, 1.0);

    vec3 hdr = vec3(0.0);

    if (NdotL > 0.0) {
        vec3 H = normalize(LIGHT_DIR + V);
        float NdotH = clamp(dot(N, H), 0.0, 1.0);
        float HdotV = clamp(dot(H, V), 0.0, 1.0);

        float F = FresnelSchlick(HdotV, F0);
        float D = DistributionGGX(NdotH, roughness);
        float G = GeometrySmith(NdotV, NdotL, roughness);

        float specular = (D * G * F) / max(4.0 * NdotV * NdotL, 1e-4);
        float kD = 1.0 - F;

        hdr += (kD * albedo * INV_PI + vec3(specular)) * u_dir_radiance * NdotL;
    }

    {
        float Fenv = FresnelSchlickRoughness(NdotV, F0, roughness);
        float kD = 1.0 - Fenv;
        hdr += (kD * albedo * INV_PI + vec3(Fenv)) * u_env_radiance;
    }

    float Fview = FresnelSchlick(NdotV, F0);
    float transmission = clamp((1.0 - baseAlpha) * (1.0 - Fview), 0.0, 1.0);
    float opacity = 1.0 - transmission;

    accumulateEvent(state, tonemap(hdr), opacity);
}

// -----------------------------------------------------------------------------
// Isosurface crossings
// -----------------------------------------------------------------------------

uint initialInsideMask(float density) {
#if !defined(INCLUDE_ISO)
    return 0u;
#else
    uint mask = 0u;
    float d = abs(density);

#if MAX_ISOVALUE_COUNT == 1
    if (d >= abs(u_iso.values[0])) {
        mask |= 1u;
    }
#else
    for (int i = 0; i < u_iso.count; ++i) {
        if (d >= abs(u_iso.values[i])) {
            mask |= (1u << uint(i));
        }
    }
#endif

    return mask;
#endif
}

int findIsoHits(
    float d0,
    float d1,
    out float hitFractions[MAX_ISOVALUE_COUNT],
    out int hitIndices[MAX_ISOVALUE_COUNT]
) {
#if !defined(INCLUDE_ISO)
    return 0;
#else
    float delta = d1 - d0;
    if (abs(delta) < EPSILON) {
        return 0;
    }

    int hitCount = 0;

    for (int i = 0; i < u_iso.count; ++i) {
        float frac = (u_iso.values[i] - d0) / delta;
        if (frac <= 0.0 || frac >= 1.0) {
            continue;
        }

        hitFractions[hitCount] = frac;
        hitIndices[hitCount] = i;
        ++hitCount;
    }

    for (int i = 1; i < hitCount; ++i) {
        float frac = hitFractions[i];
        int index = hitIndices[i];
        int j = i - 1;

        while (j >= 0 && hitFractions[j] > frac) {
            hitFractions[j + 1] = hitFractions[j];
            hitIndices[j + 1] = hitIndices[j];
            --j;
        }

        hitFractions[j + 1] = frac;
        hitIndices[j + 1] = index;
    }

    return hitCount;
#endif
}

void processSegment(
    inout RayState state,
    vec3 p0Tex,
    vec3 p1Tex,
    float d0,
    float d1
) {
    if (distance(p0Tex, p1Tex) <= EPSILON) {
        return;
    }

#if defined(INCLUDE_ISO)
    float hitFractions[MAX_ISOVALUE_COUNT];
    int hitIndices[MAX_ISOVALUE_COUNT];
    int hitCount = findIsoHits(d0, d1, hitFractions, hitIndices);

    float lastFrac = 0.0;

    for (int h = 0; h < hitCount; ++h) {
        float frac = hitFractions[h];
        int isoIdx = hitIndices[h];

        vec3 segStart = mix(p0Tex, p1Tex, lastFrac);
        vec3 segEnd   = mix(p0Tex, p1Tex, frac);

        accumulateMediumSegment(state, segStart, segEnd);

        vec3 isoPos = segEnd;
        vec4 isoColor = sampleIsoColor(isoPos, u_iso.colors[isoIdx]);
        accumulateSurfaceHit(state, isoPos, isoColor);

        state.insideMask ^= (1u << uint(isoIdx));
        lastFrac = frac;
    }

    accumulateMediumSegment(state, mix(p0Tex, p1Tex, lastFrac), p1Tex);
#else
    accumulateMediumSegment(state, p0Tex, p1Tex);
#endif
}

// -----------------------------------------------------------------------------
// Noise
// -----------------------------------------------------------------------------

float PDnrand(vec2 n) {
    return fract(sin(dot(n, vec2(12.9898, 78.233))) * 43758.5453);
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

void main() {
    vec3 entryPos = texelFetch(u_tex_entry, ivec2(gl_FragCoord.xy), 0).xyz;
    vec3 exitPos  = texelFetch(u_tex_exit,  ivec2(gl_FragCoord.xy), 0).xyz;

    float rayLength = distance(entryPos, exitPos);
    if (rayLength < 1e-5) {
        discard;
    }

    vec3 dir = (exitPos - entryPos) / rayLength;

    float voxelSamples = max(
        1.0,
        SAMPLING_RATE * length(dir * rayLength * vec3(textureSize(u_tex_density_volume, 0)))
    );
    float baseStep = max(rayLength / ceil(voxelSamples), 1e-4);

    float jitter = PDnrand(gl_FragCoord.xy + vec2(u_time));
    float t = jitter * baseStep;

    RayState state;
    state.radiance = vec3(0.0);
    state.opacity = 0.0;
    state.insideMask = 0u;

    vec3 prevPos = entryPos;
    float prevDensity = sampleDensity(prevPos);
    state.insideMask = initialInsideMask(prevDensity);

    while (t < rayLength && state.opacity < ERT_THRESHOLD) {
        float currT = min(t, rayLength);
        vec3 currPos = entryPos + currT * dir;
        float currDensity = sampleDensity(currPos);

        processSegment(state, prevPos, currPos, prevDensity, currDensity);

        prevPos = currPos;
        prevDensity = currDensity;
        t = currT + baseStep;
    }

    if (state.opacity < ERT_THRESHOLD) {
        float exitDensity = sampleDensity(exitPos);
        processSegment(state, prevPos, exitPos, prevDensity, exitDensity);
    }

    out_color = vec4(state.radiance, state.opacity);
}