#version 410 core

// A sphere of a VdW Surface (vdw_surface.h), ray cast into the G-buffer as the opaque geometry writes it: the colour,
// the normal, the velocity of the camera's motion and the picking index of its atom. Coloured by its atom; or, when
// coloured by the points on the surface, it marks its pixels with its index + 1 (out_index, the G-buffer's
// tex.index) for vdw_shade.frag to colour in a pass over them, once per pixel seen. Kept cheap: it runs for every
// sphere along a pixel's ray, as it writes its depth (gl_FragDepth) and so is depth tested only after it ran.

in vec3 v_view;
flat in int  v_sphere;
flat in vec3 v_c;
flat in float v_r;

uniform samplerBuffer  u_sphere_color;  // the atom's colour (RGBA8)
uniform usamplerBuffer u_sphere_atom;   // the atom's index in the system
uniform mat4  u_view_to_clip;
uniform mat4  u_curr_view_to_prev_clip;
uniform vec4  u_jitter_uv;
uniform int   u_ortho;
uniform uint  u_pick_base;

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_picking;
layout(location = 4) out uint out_index;    // to a draw buffer only when coloured by the points

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

    vec4 prev = u_curr_view_to_prev_clip * vec4(x, 1.0);
    vec2 ss_vel = (clip.xy / clip.w - prev.xy / prev.w) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color    = vec4(color, 1.0);
    out_normal   = vec4(encode_normal(n), 0.0, 0.0);
    out_velocity = vec4(ss_vel, 0.0, 0.0);
    out_picking  = encode_index(u_pick_base + texelFetch(u_sphere_atom, v_sphere).r);
    out_index    = uint(v_sphere + 1);
}
