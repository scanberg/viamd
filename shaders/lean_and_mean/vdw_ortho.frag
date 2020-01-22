#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform vec4 u_color;

in VS_FS {
    smooth vec4 view_coord;
    flat vec4 view_sphere;
    smooth vec2 uv;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;

void main() {
    vec2 uv = in_frag.uv;
    float len2 = dot(uv,uv);
    if (len2 > 1.0) discard;

    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;

    vec3 view_coord = in_frag.view_coord.xyz + vec3(0, 0, radius * (sqrt(1.0 - len2)));
    vec4 clip_coord = u_proj_mat * vec4(view_coord, 1);

    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;
    out_color = u_color;
}