#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;
layout(location = 2) in vec3 in_velocity;

out Vertex {
    flat vec4 view_velocity;
    flat vec4 color;
    flat vec4 picking_color;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

void main() {
    gl_Position = u_view_mat * vec4(in_position, 1.0);
    out_vert.view_velocity = u_view_mat * vec4(in_velocity, 0.0);
    out_vert.color = in_color;
    out_vert.picking_color = pack_u32(uint(gl_VertexID));
}