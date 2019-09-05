#version 330 core
uniform mat4 u_view_proj_mat;
layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec3 v_tangent_vector;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 tangent_vector;
} out_vert;

void main() {
    out_vert.control_point = u_view_proj_mat * vec4(v_control_point, 1);
    out_vert.support_vector = u_view_proj_mat * vec4(v_support_vector, 0);
    out_vert.tangent_vector = u_view_proj_mat * vec4(v_tangent_vector, 0);
} 