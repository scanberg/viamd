#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform uint u_mask;

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;
layout(location = 2) in uint in_mask;

out Vertex {
    flat int discard_vert;
} out_vert;

void main() {
    gl_Position = u_view_mat * vec4(in_position, 1.0);
    out_vert.discard_vert = int((in_mask & u_mask) == 0U || in_color.a == 0.0);
}