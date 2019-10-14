#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

#if !defined MAX_ISOVALUE_COUNT
#  define MAX_ISOVALUE_COUNT 8
#endif // MAX_ISOVALUE_COUNT

struct IsovalueParameters {
    float values[MAX_ISOVALUE_COUNT];
    vec4 colors[MAX_ISOVALUE_COUNT];
    int count;
};

layout (std140) uniform UniformData
{
    mat4 u_view_to_model_mat;
    mat4 u_model_to_view_mat;
    mat4 u_inv_proj_mat;
    mat4 u_model_view_proj_mat;

    vec2 u_inv_res;
    float u_density_scale;
    float u_alpha_scale;

    vec3 u_clip_plane_min;
    vec3 u_clip_plane_max;

    IsovalueParameters u_iso;

    vec3 u_gradient_spacing_world_space;
    mat3 u_gradient_spacing_tex_space;
};

layout(location = 0) in vec3 in_pos;

out vec3 model_pos;
out vec3 model_eye;

void main() {
    model_pos = u_clip_plane_min.xyz + in_pos * (u_clip_plane_max.xyz - u_clip_plane_min.xyz);
    model_eye = u_view_to_model_mat[3].xyz;
    gl_Position = u_model_view_proj_mat * vec4(model_pos, 1);
}