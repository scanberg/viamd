#version 410 core

// A sphere of a VdW Surface (vdw_surface.h), ray cast into the G-buffer as the opaque geometry writes it: the colour,
// the normal, the velocity of the camera's motion and the picking index of its atom. Coloured by its atom, or by the
// points on the surface: the values under their Gaussians about the hit, sum a v g / sum a g, looked up in the grid
// of cells as wide as the widest reaches.

in vec3 v_view;
flat in int  v_sphere;
flat in vec3 v_c;
flat in float v_r;

uniform samplerBuffer  u_sphere_color;  // the atom's colour (RGBA8)
uniform usamplerBuffer u_sphere_atom;   // the atom's index in the system
uniform samplerBuffer  u_point;         // xyz: the point (world), w: its value; sorted by cell
uniform samplerBuffer  u_point_w;       // x: its area, y: its Gaussian's exponent (1/Angstrom^2)
uniform usamplerBuffer u_cell;          // per cell: the first of its points, how many
uniform sampler2D      u_colormap;

uniform mat4  u_view_to_world;
uniform mat4  u_view_to_clip;
uniform mat4  u_curr_view_to_prev_clip;
uniform vec4  u_jitter_uv;
uniform int   u_ortho;
uniform int   u_by_points;
uniform vec2  u_range;
uniform vec3  u_grid_origin;
uniform float u_grid_cell;
uniform ivec3 u_grid_dim;
uniform uint  u_pick_base;

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_picking;

vec2 encode_normal(vec3 n) {
    float p = sqrt(n.z * 8.0 + 8.0);
    return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

// The points' value at x (world), and whether any reaches it
vec2 point_value(vec3 x) {
    ivec3 c = clamp(ivec3(floor((x - u_grid_origin) / u_grid_cell)), ivec3(0), u_grid_dim - 1);
    float sw = 0.0, swv = 0.0;
    for (int z = max(c.z - 1, 0); z <= min(c.z + 1, u_grid_dim.z - 1); ++z)
    for (int y = max(c.y - 1, 0); y <= min(c.y + 1, u_grid_dim.y - 1); ++y)
    for (int k = max(c.x - 1, 0); k <= min(c.x + 1, u_grid_dim.x - 1); ++k) {
        uvec2 cell = texelFetch(u_cell, (z * u_grid_dim.y + y) * u_grid_dim.x + k).xy;
        for (uint o = cell.x; o < cell.x + cell.y; ++o) {
            vec4 P = texelFetch(u_point, int(o));
            vec2 W = texelFetch(u_point_w, int(o)).xy;
            vec3 d = x - P.xyz;
            float g = W.x * exp(-W.y * dot(d, d));
            sw  += g;
            swv += g * P.w;
        }
    }
    return sw > 0.0 ? vec2(swv / sw, 1.0) : vec2(0.0);
}

void main() {
    if (v_sphere < 0) discard;
    vec3 ro = (u_ortho != 0) ? vec3(v_view.xy, 0.0) : vec3(0.0);
    vec3 rd = (u_ortho != 0) ? vec3(0.0, 0.0, -1.0) : normalize(v_view);
    vec3  oc = ro - v_c;
    float b  = dot(oc, rd);
    float h  = b * b - (dot(oc, oc) - v_r * v_r);
    if (h < 0.0) discard;
    float t  = -b - sqrt(h);
    vec3 x = ro + rd * t;
    vec3 n = (x - v_c) / v_r;

    vec4 clip = u_view_to_clip * vec4(x, 1.0);
    gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;

    vec3 color = texelFetch(u_sphere_color, v_sphere).rgb;
    if (u_by_points != 0) {
        vec2 fv = point_value((u_view_to_world * vec4(x, 1.0)).xyz);
        float ext = u_range.y - u_range.x;
        float param = abs(ext) > 0.0 ? clamp((fv.x - u_range.x) / ext, 0.0, 1.0) : 0.5;
        color = fv.y > 0.0 ? texture(u_colormap, vec2(param * (255.0 / 256.0) + 0.5 / 256.0, 0.5)).rgb : vec3(0.5);
    }

    vec4 prev = u_curr_view_to_prev_clip * vec4(x, 1.0);
    vec2 ss_vel = (clip.xy / clip.w - prev.xy / prev.w) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color    = vec4(color, 1.0);
    out_normal   = vec4(encode_normal(n), 0.0, 0.0);
    out_velocity = vec4(ss_vel, 0.0, 0.0);
    out_picking  = encode_index(u_pick_base + texelFetch(u_sphere_atom, v_sphere).r);
}
