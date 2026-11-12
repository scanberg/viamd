#version 410 core

// The spheres of a cavity into the G-buffer (cavity.h), each ray cast; the depth test makes their union. Writes what
// md_gl writes: the colour (grey; the full screen pass colours it), the normal, the velocity of the camera's motion
// and the picking index of the sphere.

in vec3 v_view;
flat in int  v_sphere;
flat in vec4 v_center;

uniform mat4  u_view_to_clip;
uniform mat4  u_curr_view_to_prev_clip;
uniform vec4  u_jitter_uv;
uniform int   u_ortho;
uniform uint  u_sphere_base;

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

void main() {
    vec3 c = v_center.xyz;
    float R = v_center.w;
    vec3 ro = (u_ortho != 0) ? vec3(v_view.xy, 0.0) : vec3(0.0);
    vec3 rd = (u_ortho != 0) ? vec3(0.0, 0.0, -1.0) : normalize(v_view);
    vec3  oc = ro - c;
    float b  = dot(oc, rd);
    float h  = b * b - (dot(oc, oc) - R * R);
    if (h < 0.0) discard;
    float t = -b - sqrt(h);
    vec3  x = ro + rd * t;
    vec3  n = (x - c) / R;

    vec4 clip = u_view_to_clip * vec4(x, 1.0);
    gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;

    vec4 prev = u_curr_view_to_prev_clip * vec4(x, 1.0);
    vec2 ss_vel = (clip.xy / clip.w - prev.xy / prev.w) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color    = vec4(0.5, 0.5, 0.5, 1.0);
    out_normal   = vec4(encode_normal(n), 0.0, 0.0);
    out_velocity = vec4(ss_vel, 0.0, 0.0);
    out_picking  = encode_index(u_sphere_base + uint(v_sphere));
}
