#version 330 core

#define SHADING_ENABLED

#if !defined MAX_ISOVALUE_COUNT
#  define MAX_ISOVALUE_COUNT 8
#endif // MAX_ISOVALUE_COUNT

#if !defined SHOW_ONLY_FIRST_ISO_HITS
#  define SHOW_ONLY_FIRST_ISO_HITS 0
#endif

struct IsovalueParameters {
    float values[MAX_ISOVALUE_COUNT];
    vec4  colors[MAX_ISOVALUE_COUNT];
    int count;
};

layout (std140) uniform UniformData
{
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
};

uniform IsovalueParameters u_iso;

uniform sampler2D u_tex_entry;
uniform sampler2D u_tex_exit;

uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_tf;

layout(location = 0) out vec4  out_color;
//layout(location = 1) out vec2  out_view_normal;
//out float gl_FragDepth;

const float REF_SAMPLING_RATE = 150.0;
const float ERT_THRESHOLD = 0.995;
const float samplingRate = 2.0;

float getVoxel(in vec3 samplePos) {
    return texture(u_tex_volume, samplePos).r;
}

vec4 classify(in float density) {
    float t = clamp((density - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
    return texture(u_tex_tf, vec2(t, 0.5));
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
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * lambert(env_radiance + N_dot_L * dir_radiance);
    vec3 specular = fr * (env_radiance + dir_radiance) * pow(N_dot_H, spec_exp);

    //return color.rgb * lambert(vec3(N_dot_L + 0.5) * 3.0) + specular;

    return diffuse + specular;
}

vec4 depth_to_view_coord(vec2 tc, float depth) {
    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

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
vec4 drawIsosurface(in vec4 curResult, in float isovalue, in vec4 isosurfaceColor, 
                    in float currentSample, in float prevSample,
                    in vec3 rayPosition, in vec3 rayDirection, 
                    in float t, in float raySegmentLen, inout float tIncr, out bool hit) {
    vec4 result = curResult;
    float sampleDelta = (currentSample - prevSample);
    hit = false;

    // check if the isovalue is lying in between current and previous sample
    // found isosurface if differences between current/prev sample and isovalue have different signs
    if ((isovalue - currentSample) * (isovalue - prevSample) < 0) {

        // apply linear interpolation between current and previous sample to obtain location of isosurface
        float a = (currentSample - isovalue) / sampleDelta;
        // if a == 1, isosurface was already computed for previous sampling position
        if (a >= 1.0) {
            hit = true;
            return result;
        }
        
        // adjust length of remaining ray segment
        tIncr = a * raySegmentLen;

        vec3 isopos = rayPosition - raySegmentLen * a * rayDirection;

        vec4 isocolor = isosurfaceColor;
#if defined(SHADING_ENABLED)
        vec3 normal = getGradient(isopos);
        normal = normalize(normal);

        vec3 isoposView = -normalize((u_model_to_view_mat * vec4(isopos, 1.0)).xyz);

        // two-sided lighting
        if (dot(normal, isoposView) < 0.0) {
            normal = -normal;
        }
        isocolor.rgb = shade(isocolor.rgb, isoposView, normal);

#endif // SHADING_ENABLED

#if defined(INCLUDE_DVR)
        // apply compositing of volumetric media from last sampling position up till isosurface
        vec4 voxelColor = classify(isovalue);
        if (voxelColor.a > 0) {
//#if defined(SHADING_ENABLED)
//            voxelColor.rgb = APPLY_LIGHTING(lighting, voxelColor.rgb, voxelColor.rgb, vec3(1.0),
//                                       isoposWorld, -gradient, toCameraDir);
//#endif // SHADING_ENABLED

            result = compositing(result, voxelColor, raySegmentLen - tIncr);
        }
#endif // INCLUDE_DVR

        isocolor.rgb *= isocolor.a; // use pre-multiplied alpha
        // blend isosurface color with result accumulated so far
        result += (1.0 - result.a) * isocolor;
        hit = true;
    }

    return result;
}

vec4 drawIsosurfaces(in vec4 curResult,
                     in float voxel, in float previousVoxel,
                     in vec3 rayPosition, in vec3 rayDirection,
                     inout float t, inout float tIncr, inout uint iso_surface_hit) {
    // in case of zero isovalues return current color
    vec4 result = curResult;
    float raySegmentLen = tIncr;
    bool hit;

#if MAX_ISOVALUE_COUNT == 1
#if SHOW_ONLY_FIRST_ISO_HITS
    if (iso_surface_hit & 1U) {
        //t = 1.0e19;
        return result;
    }
#endif
    result = drawIsosurface(result, u_iso.values[0], u_iso.colors[0],
                            voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr, hit);
#else // MAX_ISOVALUE_COUNT
    // multiple isosurfaces, need to determine order of traversal
    if (voxel - previousVoxel > 0) {
        for (int i = 0; i < u_iso.count; ++i) {
            vec4 res = drawIsosurface(result, u_iso.values[i], u_iso.colors[i],
                                      voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr, hit);
#if SHOW_ONLY_FIRST_ISO_HITS
            if (hit) {
                uint bit = 1U << uint(i);
                result = ((iso_surface_hit & bit) == 0U) ? res : result;
                iso_surface_hit = iso_surface_hit | bit;
            }
#else
            result = res;
#endif
        }
    } else {
        for (int i = u_iso.count; i > 0; --i) {
            vec4 res = drawIsosurface(result, u_iso.values[i - 1], u_iso.colors[i - 1],
                                  voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr, hit);
#if SHOW_ONLY_FIRST_ISO_HITS
            if (hit) {
                uint bit = 1U << uint(i);
                result = ((iso_surface_hit & bit) == 0U) ? res : result;
                iso_surface_hit = iso_surface_hit | bit;
            }
#else
            result = res;
#endif
        }
    }
#endif // MAX_ISOVALUE_COUNT
    return result;
}

float PDnrand( vec2 n ) {
    return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453 );
}

vec2 PDnrand2( vec2 n ) {
    return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* vec2(43758.5453, 28001.8384) );
}

vec2 encode_normal (vec3 n) {
   float p = sqrt(n.z * 8 + 8);
   return n.xy / p + 0.5;
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
    uint iso_surface_hit = 0U;
    vec3 samplePos = entryPos;

    while (t < tEnd) {
        samplePos = entryPos + t * dir;
        float prevDensity = density;
        density = getVoxel(samplePos);

#if defined(INCLUDE_ISO)
        result = drawIsosurfaces(result, density, prevDensity, samplePos, dir, t, tIncr, iso_surface_hit);
#endif 

#if defined(INCLUDE_DVR)
        vec4 srcColor = classify(density);
        if (srcColor.a > 0.0) {
            result = compositing(result, srcColor, tIncr);
        }
#endif

        if (result.a > ERT_THRESHOLD) {
            t = tEnd;
        } else {
            // make sure that tIncr has the correct length since drawIsoSurface will modify it
            tIncr = baseIncr;
            t += tIncr;
        }
    }

    result.rgb = Tonemap(result.rgb);
    out_color = result;
}