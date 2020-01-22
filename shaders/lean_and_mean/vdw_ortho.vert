#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform float u_radius_scale = 1.0;
uniform uint u_mask;

uniform samplerBuffer u_buf_position;
uniform samplerBuffer u_buf_radius;
uniform samplerBuffer u_buf_color;
uniform usamplerBuffer u_buf_mask;

out VS_FS {
    smooth vec4 view_coord;
    flat   vec4 view_sphere;
    smooth vec2 uv;
} out_vert;

const float SQRT_3 = 1.7329598975688772;
const vec2 V_OFFSETS[3] = vec2[3](vec2(-SQRT_3,-1), vec2(SQRT_3,-1), vec2(0,2));

void main() {
    int g_idx = gl_VertexID / 3;
    int l_idx = gl_VertexID % 3;
    vec2 uv = V_OFFSETS[l_idx];

    vec3 position = texelFetch(u_buf_position, g_idx).xyz;
    float radius = texelFetch(u_buf_radius, g_idx).r;
    vec4 color = texelFetch(u_buf_color, g_idx).rgba;
    uint mask = texelFetch(u_buf_mask, g_idx).x;

    if ((mask & u_mask) == 0U || color.a == 0.0) {
        gl_Position = vec4(0,0,0,1);
        return;
    }

    float rad = radius * u_radius_scale;
    vec4 view_position = u_view_mat * vec4(position, 1.0);
    vec4 view_coord = view_position + vec4(uv * rad, 0, 0);
    vec4 view_sphere = vec4(view_position.xyz, rad);

    gl_Position = u_proj_mat * (view_coord + vec4(0, 0, rad, 0));

    out_vert.view_coord = view_coord;
    out_vert.view_sphere = view_sphere;
    out_vert.uv = uv;
}