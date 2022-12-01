#version 330 core

uniform sampler2D u_tex_sdf;

in  vec3 model_pos;
in  vec3 model_eye;

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


//vec4 depth_to_view_coord(vec2 tc, float depth) {
//    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
//    vec4 view_coord = u_inv_proj_mat * clip_coord;
//    return view_coord / view_coord.w;
//}

void main() {
    vec3 ray_o = model_pos; // Start at point of entry not at eye...
    vec3 ray_d = normalize(model_pos - model_eye);

    float t_entry, t_exit;
    ray_vs_aabb(t_entry, t_exit, ray_o, ray_d, vec3(0,0,0), vec3(1,1,1));

    float voxel_ext = 1.0 / 256.0;
    int lod = 0;

    float t = 0;
    vec4 result = vec4(1,0,0,1);
    while (t < t_exit) {
        t += voxel_ext;
    }

    out_frag = result;
}