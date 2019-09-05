#version 330 core

layout (location = 0) in vec3 v_position;

out uint atom_index;

void main() {
    atom_index = uint(gl_VertexID);
    gl_Position = vec4(v_position, 1);
} 
