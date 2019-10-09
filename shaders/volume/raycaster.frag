#version 150 core

// #define INCLUDE_ISOSURFACES

#if (!defined(INCLUDE_DVR) && !defined(INCLUDE_ISOSURFACES))
#  define INCLUDE_DVR
#endif

#if !defined MAX_ISOVALUE_COUNT
#  define MAX_ISOVALUE_COUNT 6
#endif // MAX_ISOVALUE_COUNT

// need to ensure there is always at least one isovalue due to the use of the macro
// as array size in IsovalueParameters which will cause an error for size = 0.
#if MAX_ISOVALUE_COUNT < 1
#  undef MAX_ISOVALUE_COUNT
#  define MAX_ISOVALUE_COUNT 1
#endif


struct IsovalueParameters {
    float values[MAX_ISOVALUE_COUNT];
    vec4 colors[MAX_ISOVALUE_COUNT];
};

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_tf;
uniform float     u_scale = 1.0;
uniform float     u_alpha_scale = 1.0;
uniform vec2      u_inv_res;
uniform mat4      u_view_to_model_mat;
uniform mat4      u_model_to_tex_mat;
uniform mat4      u_inv_proj_mat;
uniform int       u_iso_enabled;
uniform IsovalueParameters u_isovalues;

in  vec3 model_pos;
in  vec3 model_eye;
in  vec3 color;

out vec4 out_frag;

const float REF_SAMPLING_RATE = 150.0;
const float ERT_THRESHOLD = 0.99;
const float samplingRate = 8.0;

vec4 classify(float density) {
    vec4 c = texture(u_tex_tf, vec2(density, 0.5));
    c.a = clamp(c.a * u_alpha_scale, 0.0, 1.0);
    return c;
}

vec4 compositing(in vec4 dstColor, in vec4 srcColor, in float tIncr) {
    srcColor.a = 1.0 - pow(1.0 - srcColor.a, tIncr * REF_SAMPLING_RATE);
    // pre-multiplied alpha
    srcColor.rgb *= srcColor.a;
    return dstColor + (1.0 - dstColor.a) * srcColor;
}

// Modified version from source found here:
// https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms#18459

bool ray_vs_aabb(out float t_entry, out float t_exit, in vec3 ori, in vec3 dir, in vec3 min_box, in vec3 max_box) {
    vec3 dir_frac = 1.0 / dir;

    vec3 tv_1 = (min_box - ori) * dir_frac;
    vec3 tv_2 = (max_box - ori) * dir_frac;

    vec3 tv_min = min(tv_1, tv_2);
    vec3 tv_max = max(tv_1, tv_2); 

    float t_min = max(max(tv_min.x, tv_min.y), tv_min.z);
    float t_max = min(min(tv_max.x, tv_max.y), tv_max.z);

    if (t_max < 0 || t_min > t_max) {
        return false;
    }

    t_entry = t_min;
    t_exit  = t_max;

    return true;
} 

vec4 depth_to_view_coord(vec2 tc, float depth) {
    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

/*
vec4 fetch_voxel(vec3 tc) {
    float a = min(texture(u_tex_volume, tc).x * u_scale, 1.0);
    return vec4(mix(vec3(0), u_color, min(a * 10.0, 1.0)), a);
}
*/

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
                    in float t, in float raySegmentLen, inout float tIncr) {
    vec4 result = curResult;
    float sampleDelta = (currentSample - prevSample);

    // check if the isovalue is lying in between current and previous sample
    // found isosurface if differences between current/prev sample and isovalue have different signs
    if ((isovalue - currentSample) * (isovalue - prevSample) <= 0) {
        // apply linear interpolation between current and previous sample to obtain location of isosurface
        float a = (currentSample - isovalue) / sampleDelta;
        // if a == 1, isosurface was already computed for previous sampling position
        if (a >= 1.0) return result;
        
        // adjust length of remaining ray segment
        tIncr = a * raySegmentLen;

        vec3 isopos = rayPosition - raySegmentLen * a * rayDirection;

        vec4 isocolor = isosurfaceColor;
#if defined(SHADING_ENABLED)
        vec3 gradient = COMPUTE_GRADIENT_FOR_CHANNEL(voxel, volume, volumeParameters, isopos, channel);
        gradient = normalize(gradient);

        // two-sided lighting
        if (dot(gradient, rayDirection) <= 0) {
            gradient = -gradient;
        }

        vec3 isoposWorld = (volumeParameters.textureToWorld * vec4(isopos, 1.0)).xyz;
        isocolor.rgb = APPLY_LIGHTING(lighting, isocolor.rgb, isocolor.rgb, vec3(1.0),
                           isoposWorld, -gradient, toCameraDir);
#endif // SHADING_ENABLED

#if defined(INCLUDE_DVR)
        // apply compositing of volumetric media from last sampling position up till isosurface
        vec4 voxelColor = classify(isovalue);
        if (voxelColor.a > 0) {
#if defined(SHADING_ENABLED)
            voxelColor.rgb = APPLY_LIGHTING(lighting, voxelColor.rgb, voxelColor.rgb, vec3(1.0),
                                       isoposWorld, -gradient, toCameraDir);
#endif // SHADING_ENABLED

            result = compositing(result, voxelColor, raySegmentLen - tIncr);
        }
#endif // INCLUDE_DVR

        isocolor.rgb *= isocolor.a; // use pre-multiplied alpha
        // blend isosurface color with result accumulated so far
        result += (1.0 - result.a) * isocolor;
    }

    return result;
}

vec4 drawIsosurfaces(in vec4 curResult,
                     in float voxel, in float previousVoxel,
                     in vec3 rayPosition, in vec3 rayDirection,
                     in float t, inout float tIncr) {
    // in case of zero isovalues return current color
    vec4 result = curResult;
    float raySegmentLen = tIncr;

#if MAX_ISOVALUE_COUNT == 1
    result = drawIsosurface(result, u_isovalues.values[0], u_isovalues.colors[0],
                            voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
#else // MAX_ISOVALUE_COUNT
    // multiple isosurfaces, need to determine order of traversal
    if (voxel - previousVoxel > 0) {
        for (int i = 0; i < MAX_ISOVALUE_COUNT; ++i) {
            result = drawIsosurface(result, u_isovalues.values[i], u_isovalues.colors[i],
                voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
        }
    } else {
        for (int i = MAX_ISOVALUE_COUNT; i > 0; --i) {
            result = drawIsosurface(result, u_isovalues.values[i - 1], u_isovalues.colors[i - 1],
                voxel, previousVoxel, rayPosition, rayDirection, t, raySegmentLen, tIncr);
        }
    }
#endif // MAX_ISOVALUE_COUNT
    return result;
}


float getVoxel(vec3 samplePos) {
    return texture(u_tex_volume, samplePos).r * u_scale;
}

void main() {
    out_frag = vec4(model_pos, 1.0);
    //out_frag = classify(model_pos.x);
    //return;

    // Do everything in model space
    vec3 ori = model_eye;
    vec3 dir = normalize(model_pos - model_eye);

    float t_entry, t_exit;
    if (!ray_vs_aabb(t_entry, t_exit, ori, dir, vec3(0), vec3(1))) discard;
    t_entry = max(0, t_entry);

    float depth = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
    if (depth < 1.0) {
        vec3 stop_pos = (u_view_to_model_mat * depth_to_view_coord(gl_FragCoord.xy * u_inv_res, depth)).xyz;
        t_exit = min(t_exit, dot(stop_pos - ori, dir));
    }

    if (t_entry == t_exit) discard;

    vec3 entryPos = ori + dir * t_entry;
    float tEnd = t_exit - t_entry;

    float tIncr = min(tEnd, tEnd / (samplingRate * length(dir * tEnd * textureSize(u_tex_volume, 0))));
    float samples = ceil(tEnd / tIncr);
    tIncr = tEnd / samples;

    vec4 result = vec4(0);
    float t = 0.5 * tIncr;
    float density = getVoxel(entryPos + t * dir); // need this for isosurface rendering
    while (t < tEnd) {
        vec3 samplePos = entryPos + t * dir;

        float prevDensity = density;
        float density = getVoxel(samplePos);

        if (u_iso_enabled > 0) {
            result = drawIsosurfaces(result, density, prevDensity, samplePos, dir, t, tIncr);
        } 
#if defined(INCLUDE_DVR)
        vec4 srcColor = classify(density);
        if (srcColor.a > 0.0) {
            result = compositing(result, srcColor, tIncr);
        }
#endif // INCLUDE_DVR
        if (result.a > ERT_THRESHOLD) {
            t = tEnd;
        } else {
            // make sure that tIncr has the correct length since drawIsoSurface will modify it
            tIncr = tEnd / samples;

            t += tIncr;
        }
    }

    out_frag = result;
}