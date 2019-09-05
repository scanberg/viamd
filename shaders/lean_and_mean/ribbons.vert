#version 330 core

layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec3 v_support_tangent;
layout (location = 3) in uint v_atom_index;

uniform samplerBuffer u_atom_color_buffer;
uniform usamplerBuffer u_atom_mask_buffer;
uniform uint u_mask;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    int discard_vert;
} out_vert;

void main() {
    uint mask = texelFetch(u_atom_mask_buffer, int(v_atom_index)).x;
    vec4 color = texelFetch(u_atom_color_buffer, int(v_atom_index)).rgba;
    out_vert.control_point = vec4(v_control_point, 1);
    out_vert.support_vector = vec4(v_support_vector, 0);
    out_vert.support_tangent = vec4(v_support_tangent, 0);
    out_vert.discard_vert = int((mask & u_mask) == 0U || color.a == 0.0);
}