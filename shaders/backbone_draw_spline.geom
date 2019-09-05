#version 330 core

layout(lines) in;
layout(line_strip, max_vertices = 6) out;

uniform vec4 u_s_color = vec4(0,1,0,1);
uniform vec4 u_v_color = vec4(1,0,0,1);
uniform vec4 u_t_color = vec4(0,0,1,1);

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 tangent_vector;
} in_vert[];

out vec4 color;

void main() {
    color = u_v_color;
    gl_Position = in_vert[0].control_point;
    EmitVertex();
    gl_Position = in_vert[0].control_point + in_vert[0].support_vector;
    EmitVertex();
    EndPrimitive();

    color = u_t_color;
    gl_Position = in_vert[0].control_point;
    EmitVertex();
    gl_Position = in_vert[0].control_point + in_vert[0].tangent_vector;
    EmitVertex();
    EndPrimitive();

    color = u_s_color;
    gl_Position = in_vert[0].control_point;
    EmitVertex();
    gl_Position = in_vert[1].control_point;
    EmitVertex();
    EndPrimitive();
}