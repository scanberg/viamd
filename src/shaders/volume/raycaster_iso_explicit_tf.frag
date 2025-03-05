#version 330 core

#define SHADING_ENABLED

#if !defined MAX_ISOVALUE_COUNT
#  define MAX_ISOVALUE_COUNT 8
#endif // MAX_ISOVALUE_COUNT

struct IsovalueParameters {
    float values[MAX_ISOVALUE_COUNT];
    vec4  colors[MAX_ISOVALUE_COUNT];
    int count;
};

layout (std140) uniform UniformData {
    mat4 u_view_to_model_mat;
    mat4 u_model_to_view_mat;
    mat4 u_inv_proj_mat;
    mat4 u_model_view_proj_mat;

    vec2  u_inv_res;
    float u_time;

    vec3  u_clip_plane_min;
    float u_tf_min;
    vec3  u_clip_plane_max;
    float u_tf_inv_ext;

    vec3 u_gradient_spacing_world_space;
    mat4 u_gradient_spacing_tex_space;

    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;
};

uniform IsovalueParameters u_iso;

uniform sampler2D u_tex_entry;
uniform sampler2D u_tex_exit;

uniform sampler3D u_tex_volume;
uniform sampler3D u_tex_tf_volume;
uniform sampler2D u_tex_tf;

layout(location = 0) out vec4  out_color;

const float REF_SAMPLING_RATE = 150.0;
const float ERT_THRESHOLD = 0.995;
const float samplingRate = 2.0;

float getVoxel(in vec3 samplePos) {
    return texture(u_tex_volume, samplePos).r;
}

float getTfVoxel(in vec3 samplePos) {
    return texture(u_tex_tf_volume, samplePos).r;
}

vec4 classify(in float density) {
    float t  = clamp((density - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
    vec4 col = texture(u_tex_tf, vec2(t, 0.5));
    return col;
}

vec4 compositing(in vec4 dstColor, in vec4 srcColor, in float tIncr) {
    srcColor.a = 1.0 - pow(1.0 - srcColor.a, tIncr * REF_SAMPLING_RATE);
    // pre-multiplied alpha
    srcColor.rgb *= srcColor.a;
    return dstColor + (1.0 - dstColor.a) * srcColor;
}

vec3 getGradient(in vec3 samplePos) {
    vec3 g = vec3(getVoxel(samplePos + u_gradient_spacing_tex_space[0].xyz), 
                  getVoxel(samplePos + u_gradient_spacing_tex_space[1].xyz),
                  getVoxel(samplePos + u_gradient_spacing_tex_space[2].xyz));
    g -= vec3(getVoxel(samplePos - u_gradient_spacing_tex_space[0].xyz), 
              getVoxel(samplePos - u_gradient_spacing_tex_space[1].xyz),
              getVoxel(samplePos - u_gradient_spacing_tex_space[2].xyz));
    return g / (2.0 * u_gradient_spacing_world_space);
}

const vec3 L = normalize(vec3(1,1,1));
const float spec_exp = 100.0;
const float PI = 3.1415926535;
const float ONE_OVER_PI = 1.0 / 3.1415926535;

// Tonemapping operations (in case one writes to post tonemap)
float max3(float x, float y, float z) { return max(x, max(y, z)); }
float rcp(float x) { return 1.0 / x; }

const float exposure_bias = 0.5;
const float gamma = 2.2;

vec3 Tonemap(vec3 c) {
    c = c * exposure_bias;
    c = c * rcp(max3(c.r, c.g, c.b) + 1.0);
    c = pow(c, vec3(1.0 / gamma));
    return c;
}

vec3 TonemapInvert(vec3 c) {
    c = pow(c, vec3(gamma));
    c = c * rcp(1.0 - max3(c.r, c.g, c.b));
    c = c / exposure_bias;
    return c;
}

// Cook-Torrance model From here:
// https://learnopengl.com/PBR/Theory

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

vec4 shade(vec4 color, vec3 V, vec3 N) {
    vec3  H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float N_dot_V = max(0.0, dot(N, V));

    float F0 = u_F0;
    float roughness = u_roughness;
    vec3  dir_radiance = u_dir_radiance;
    vec3  env_radiance = u_env_radiance;

    // Premultiply albedo with the alpha of the surface
    vec3  albedo = color.rgb * color.a;
    float alpha = color.a;

    vec3 Lo = vec3(0);

    {
        // Add contribution from directional light
        float NDF = DistributionGGX(N_dot_H, roughness);
        float G   = GeometrySmith(N_dot_V, N_dot_L, roughness);
        float F   = FresnelSchlick(H_dot_V, F0);
        float kS = F;
        float kD = 1.0 - kS;
        float numerator   = NDF * G * F;
        float denominator = 4.0 * N_dot_V * N_dot_L + 0.0001;
        float specular    = numerator / denominator;
        Lo += (kD * albedo * ONE_OVER_PI + vec3(specular)) * dir_radiance * N_dot_L;
    }

    {
        // Add contribution from environment
        float F = FresnelSchlickRoughness(N_dot_V, F0, roughness);
        float kS = F;
        float kD = 1.0 - kS;
        Lo += (kD * albedo * ONE_OVER_PI + vec3(kS)) * env_radiance; // Multiply with ao here

        // Modify the alpha term to 'approximate' the transmissive component in the reflection occuring at the surface
        // Only some portion of the energy enters the surface, the other portion reflects off.
        // We model this by an increase in opacity, which means the consequent contributions will be less
        alpha += kS;
    }

    return vec4(Lo, alpha);
}

/**
 * Draws an isosurface if the given isovalue is found along the ray in between the 
 * current and the previous volume sample. On return, tIncr refers to the distance
 * between the last valid isosurface and the current sampling position. That is
 * if no isosurface was found, tIncr is not modified.
 *
 * @param curResult        color accumulated so far during raycasting
 * @param isovalue         isovalue of isosurface to be drawn
 * @param isosurfaceColor  color of the isosurface used for blending
 * @param currentSample    scalar values of current sampling position
 * @param prevSample       scalar values of previous sample
 * @return in case of an isosurface, curResult is blended with the color of the isosurface. 
 *       Otherwise curResult is returned
 */
vec4 drawIsosurface(in vec4 curResult, in float isovalue, 
                    in float currentSample, in float prevSample,
                    in vec3 rayPosition, in vec3 rayDirection, 
                    in float t, in float raySegmentLen, inout float tIncr) {
    vec4 result = curResult;
    float sampleDelta = (currentSample - prevSample);

    // check if the isovalue is lying in between current and previous sample
    // found isosurface if differences between current/prev sample and isovalue have different signs
    if ((isovalue - currentSample) * (isovalue - prevSample) < 0) {

        // apply linear interpolation between current and previous sample to obtain location of isosurface
        float a = (currentSample - isovalue) / sampleDelta;
        // if a == 1, isosurface was already computed for previous sampling position
        if (a >= 1.0) {
            return result;
        }
        
        // adjust length of remaining ray segment
        tIncr = a * raySegmentLen;

        vec3 isopos = rayPosition - raySegmentLen * a * rayDirection;

        float isoTfValue = getTfVoxel(isopos);
        vec4 isocolor = classify(isoTfValue);
#if defined(SHADING_ENABLED)
        vec3 normal = getGradient(isopos);
        normal = normalize(normal);

        vec3 isoposView = -normalize((u_model_to_view_mat * vec4(isopos, 1.0)).xyz);

        // two-sided lighting
        if (dot(normal, isoposView) < 0.0) {
            normal = -normal;
        }

        isocolor = shade(isocolor, isoposView, normal);
#endif // SHADING_ENABLED

        // blend isosurface color with result accumulated so far
        result += (1.0 - result.a) * isocolor;
    }

    return result;
}

vec4 drawIsosurfaces(in vec4 curResult,
                     in float voxel, in float previousVoxel,
                     in vec3 rayPosition, in vec3 rayDirection,
                     inout float t, inout float tIncr) {
    // in case of zero isovalues return current color
    vec4 result = curResult;
    float raySegmentLen = tIncr;

#if MAX_ISOVALUE_COUNT == 1
    result = drawIsosurface(result, u_iso.values[0], voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
#else // MAX_ISOVALUE_COUNT
    // multiple isosurfaces, need to determine order of traversal
    if (voxel - previousVoxel > 0) {
        for (int i = 0; i < u_iso.count; ++i) {
            vec4 res = drawIsosurface(result, u_iso.values[i], voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
            result = res;
        }
    } else {
        for (int i = u_iso.count; i > 0; --i) {
            vec4 res = drawIsosurface(result, u_iso.values[i - 1], voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
            result = res;
        }
    }
#endif // MAX_ISOVALUE_COUNT
    return result;
}

float PDnrand( vec2 n ) {
    return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453 );
}

void main() {
    // Entry and exit positions are given in the volumes texture space
    vec3 entryPos = texelFetch(u_tex_entry, ivec2(gl_FragCoord.xy), 0).xyz;
    vec3 exitPos  = texelFetch(u_tex_exit,  ivec2(gl_FragCoord.xy), 0).xyz;

    if (entryPos == exitPos) discard;

    float dist = distance(entryPos, exitPos);

    vec3 ori = entryPos;
    vec3 dir = (exitPos - entryPos) / dist;

    float tEnd = dist;

    float jitter = PDnrand(gl_FragCoord.xy + vec2(u_time, u_time));

    float tIncr = min(tEnd, tEnd / (samplingRate * length(dir * tEnd * textureSize(u_tex_volume, 0))));
    float samples = ceil(tEnd / tIncr);
    float baseIncr = tEnd / samples;

	tIncr = baseIncr;

    float t = jitter * tIncr;

    vec4 result = vec4(0);
    float density = 0.0;
    vec3 samplePos = entryPos;

    while (t < tEnd) {
        samplePos = entryPos + t * dir;
        float prevDensity = density;
        density = getVoxel(samplePos);

        result = drawIsosurfaces(result, density, prevDensity, samplePos, dir, t, tIncr);

        if (result.a > ERT_THRESHOLD) {
            t = tEnd;
        } else {
            // make sure that tIncr has the correct length since drawIsoSurface will modify it
            tIncr = baseIncr;
            t += tIncr;
        }
    }

#ifdef SHADING_ENABLED
    result.rgb = Tonemap(result.rgb);
#endif

    out_color = result;
}