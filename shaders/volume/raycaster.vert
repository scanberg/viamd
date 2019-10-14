#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

layout(location = 0) in vec3 in_pos;

layout (std140) uniform UniformData
{
    mat4      u_view_to_model_mat;
    mat4      u_model_to_view_mat;
    mat4      u_inv_proj_mat;
    mat4      u_model_view_proj_mat;

    float     u_global_scale = 1.0;
    float     u_alpha_scale = 1.0;
    vec2      u_inv_res;

    float     u_isovalues[8];
    vec4      u_isocolors[8];
    int       u_iso_count;
    int       u_iso_enabled;
    float     _pad0[2];

    vec3      u_gradient_spacing_world_space;
    float     _pad1[1];

    mat3      u_gradient_spacing_tex_space;
};

out vec3 model_pos;
out vec3 model_eye;

void main() {
    model_pos = in_pos;
    model_eye = u_view_to_model_mat[3].xyz;
    gl_Position = u_model_view_proj_mat * vec4(in_pos, 1);
}