#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform float u_radius_scale = 1.0;

uniform samplerBuffer u_buf_position;
uniform samplerBuffer u_buf_radius;
uniform samplerBuffer u_buf_color;
uniform samplerBuffer u_buf_view_velocity;

out VS_FS {
    smooth vec4 view_coord;
    flat vec4 view_sphere;
    flat vec4 view_velocity;
    flat vec4 picking_color;
    flat vec4 color;
    smooth vec2 uv;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

const float SQRT_3 = 1.7329598975688772;

uniform vec2 V_OFFSETS[3] = { vec2(-SQRT_3,-1), vec2(SQRT_3,-1), vec2(0,2) };

void main() {
    int g_idx = gl_VertexID / 3;
    int l_idx = gl_VertexID % 3;
    vec2 uv = V_OFFSETS[l_idx];

    vec3 position = texelFetch(u_buf_position, g_idx).xyz;
    float radius = texelFetch(u_buf_radius, g_idx).r;
    vec4 color = texelFetch(u_buf_color, g_idx).rgba;
    vec3 view_velocity = texelFetch(u_buf_view_velocity, g_idx).xyz;

    float rad = radius * u_radius_scale;
    vec4 view_position = u_view_mat * vec4(position, 1.0);
    vec4 view_coord = view_position + vec4(uv * rad, 0, 0);
    vec4 view_sphere = vec4(view_position.xyz, rad);
    vec4 view_vel = vec4(view_velocity, 0);

    gl_Position = u_proj_mat * (view_coord + vec4(0,0,rad,0));
    //gl_Position = vec4(offset * 1.0, 1.0);

    out_vert.view_coord = view_coord;
    out_vert.view_sphere = view_sphere;
    out_vert.view_velocity = view_vel;
    out_vert.picking_color = pack_u32(uint(gl_VertexID));
    out_vert.color = color;
    out_vert.uv = uv;
}