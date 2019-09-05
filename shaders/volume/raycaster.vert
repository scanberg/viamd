#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

layout(location = 0) in vec3 in_pos;

uniform mat4 u_view_to_model_mat;
uniform mat4 u_model_view_proj_mat;

out vec3 model_pos;
out vec3 model_eye;
out vec3 color;

void main() {
    color = mix(vec3(1,0,0), vec3(0,1,0), float(gl_VertexID) / 42.0);
    model_pos = in_pos;
    model_eye = u_view_to_model_mat[3].xyz;
    gl_Position = u_model_view_proj_mat * vec4(in_pos, 1);
}