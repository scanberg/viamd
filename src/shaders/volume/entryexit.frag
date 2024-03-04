#version 330 core

layout (std140) uniform UniformData
{
    mat4 u_view_to_model_mat;
    mat4 u_model_to_view_mat;
    mat4 u_inv_proj_mat;
    mat4 u_model_view_proj_mat;

    vec2  u_inv_res;
    float u_density_scale;

    vec3  u_clip_plane_min;
    vec3  u_clip_plane_max;
    float u_time;

    vec3 u_gradient_spacing_world_space;
    mat4 u_gradient_spacing_tex_space;
};

#ifdef SAMPLE_DEPTH
uniform sampler2D u_tex_depth;
#endif

in  vec3 model_pos;
in  vec3 model_eye;

layout(location = 0) out vec3 out_entry;
layout(location = 1) out vec3 out_exit;

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

void main() {

    // Do everything in model space
    vec3 ori = model_eye;
    vec3 dir = normalize(model_pos - model_eye);

    float t_entry, t_exit;
    if (!ray_vs_aabb(t_entry, t_exit, ori, dir, u_clip_plane_min.xyz, u_clip_plane_max.xyz)) discard;
    t_entry = max(0, t_entry);

#ifdef SAMPLE_DEPTH
    float depth = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
    if (depth < 1.0) {
        vec3 stop_pos = (u_view_to_model_mat * depth_to_view_coord(gl_FragCoord.xy * u_inv_res, depth)).xyz;
        t_exit = min(t_exit, dot(stop_pos - ori, dir));
    }
#endif

    out_entry = ori + dir * t_entry;
    out_exit  = ori + dir * t_exit;
}