#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

#ifndef WRITE_DEPTH
#define WRITE_DEPTH 1
#endif

uniform vec4 u_color = vec4(1,1,1,1);
uniform mat4 u_proj_mat;

in GS_FS {
#if WRITE_DEPTH
    flat vec4 view_sphere;
    smooth vec4 view_coord;
#endif
    smooth vec2 uv;
} in_frag;

#if WRITE_DEPTH
#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
#endif
layout(location = 0) out vec4 out_color;

void main() {
    //if (dot(in_frag.uv, in_frag.uv) >= 1.0) discard;

#if WRITE_DEPTH
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;
    vec3 view_dir = -normalize(in_frag.view_coord.xyz);

    vec3 m = -center;
    vec3 d = -view_dir;
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);

    vec3 view_coord = d * t;
    vec4 clip_coord = u_proj_mat * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;
#endif

    out_color = u_color;
}