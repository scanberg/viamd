#version 150 core

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_tf;
uniform float     u_scale = 1.0;
uniform float     u_alpha_scale = 1.0;
uniform vec2      u_inv_res;
uniform mat4      u_view_to_model_mat;
uniform mat4      u_model_to_tex_mat;
uniform mat4      u_inv_proj_mat;

in  vec3 model_pos;
in  vec3 model_eye;
in  vec3 color;

out vec4 out_frag;

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

vec4 classify(float density) {
    vec4 c = texture(u_tex_tf, vec2(density, 0.5));
    c.a = clamp(c.a * u_alpha_scale, 0.0, 1.0);
    return c;
}

const float REF_SAMPLING_RATE = 150.0;
const float ERT_THRESHOLD = 0.99;
const float samplingRate = 8.0;

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
    while (t < tEnd) {
        vec3 samplePos = entryPos + t * dir;

        if (samplePos.y > 0.5) {
            t += tIncr;
            continue;
        }

        float density = texture(u_tex_volume, samplePos).r;
        vec4 srcColor = classify(density * u_scale);

        if (srcColor.a > 0.0) {
            srcColor.a = 1.0 - pow(1.0 - srcColor.a, tIncr * REF_SAMPLING_RATE);
            // pre-multiplied alpha
            srcColor.rgb *= srcColor.a;
            result += (1.0 - result.a) * srcColor;
        }

        if (result.a > ERT_THRESHOLD) {
            t = tEnd;
        } else {
            t += tIncr;
        }
    }

    out_frag = result;
}