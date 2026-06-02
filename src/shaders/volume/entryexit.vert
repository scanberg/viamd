#version 330 core

layout (std140) uniform UniformData
{
    mat4 u_view_to_model_mat;
    mat4 u_model_to_view_mat;
    mat4 u_inv_proj_mat;
    mat4 u_model_view_proj_mat;

    vec2  u_inv_res;
    vec3  u_clip_plane_min;
    float u_time;
    vec3  u_clip_plane_max;
    float u_tf_min;

    vec3  u_gradient_spacing_world_space;
    float u_exposure;

    mat4 u_gradient_spacing_tex_space;

    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;
};

layout(location = 0) in vec3 in_pos;

out vec3 model_pos;

void main() {
    model_pos = u_clip_plane_min.xyz + in_pos * (u_clip_plane_max.xyz - u_clip_plane_min.xyz);
    gl_Position = u_model_view_proj_mat * vec4(model_pos, 1);
}