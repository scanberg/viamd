#version 330 core

uniform mat4 u_mvp_mat;
uniform float u_radius_scale = 1.0;

uniform samplerBuffer u_buf_pos;
uniform samplerBuffer u_buf_rad;

layout(location = 0) in vec3 in_icosa_pos;

void main() {
    vec3  pos = texelFetch(u_buf_pos, gl_InstanceID).xyz;
    float rad = texelFetch(u_buf_rad, gl_InstanceID).r * u_radius_scale;
    gl_Position = u_mvp_mat * vec4(pos + in_icosa_pos * rad, 1);
}