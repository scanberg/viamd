#version 330 core
uniform samplerBuffer u_atom_color_tex;

layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec3 v_support_tangent;
layout (location = 3) in vec4 v_classification;
layout (location = 4) in uint v_atom_index;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    vec4 weights;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

void main() {
    vec4 weights = vec4(0); // coil, sheet, helix weights
    weights.y = pow(v_classification.b, 1.0); // encoded as blue
    weights.z = pow(v_classification.r, 1.0); // encoded as red
    weights.x = clamp(1.0 - weights.y - weights.z, 0, 1);

    out_vert.control_point = vec4(v_control_point, 1);
    out_vert.support_vector = vec4(v_support_vector, 0);
    out_vert.support_tangent = vec4(v_support_tangent, 0);
    out_vert.color = texelFetch(u_atom_color_tex, int(v_atom_index));
    out_vert.picking_color = pack_u32(v_atom_index);
    out_vert.weights = weights;
}