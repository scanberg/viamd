#version 330 core
layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec4 v_classification;
layout (location = 3) in uint v_atom_index;

out Vertex {
    vec3 control_point;
    vec3 support_vector;
    vec4 classification;
    uint atom_index;
} out_vert;

void main() {
    out_vert.control_point = v_control_point;
    out_vert.support_vector = v_support_vector;
    out_vert.classification = v_classification;
    out_vert.atom_index = v_atom_index;
} 