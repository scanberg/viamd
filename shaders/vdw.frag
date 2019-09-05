#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform sampler2D u_tex_noise;

uniform mat4 u_proj_mat;
uniform mat4 u_curr_view_to_prev_clip_mat;
uniform vec4 u_jitter_uv;
uniform uint u_frame;

in GS_FS {
    smooth vec4 view_coord;
    flat vec4 view_sphere;
    flat vec4 view_velocity;
    flat vec4 picking_color;
    flat vec4 color;
    flat uint atom_idx;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_ss_vel;
layout(location = 3) out vec4 out_emission;
layout(location = 4) out vec4 out_post_tonemap;
layout(location = 5) out vec4 out_picking_color;

const float GOLDEN_RATIO = 1.61803398875;

float PDnrand( vec2 n ) {
    return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453 );
}

// https://aras-p.info/texts/CompactNormalStorage.html
vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    
    //vec2 coord = gl_FragCoord.xy / 512;
    //vec4 noise = textureLod(u_tex_noise, coord, 0);
    //uint idx = in_frag.atom_idx;
    //noise = vec4(PDnrand(gl_FragCoord.xy * gl_FragCoord.z));

    //float val = mod(gl_FragCoord.x + gl_FragCoord.y, 2.0);
    //bool mask = val > 0.0;
    //bool mask = noise.x > 0.125;

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
    vec3 view_normal = (view_coord - center) / radius;
    vec3 view_vel = in_frag.view_velocity.xyz;
    vec4 clip_coord = u_proj_mat * vec4(view_coord, 1);

    vec3 prev_view_coord = view_coord - view_vel;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip_mat * vec4(prev_view_coord, 1);

    // Remove jitter from samples to provide the actual velocity
    // This is crucial for the temporal reprojection to work properly
    // Otherwise the velocity will push the samples outside of the "reprojection" region
    vec2 curr_ndc = clip_coord.xy / clip_coord.w;
    vec2 prev_ndc = prev_clip_coord.xy / prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    vec4 color = in_frag.color;
    vec4 picking_color = in_frag.picking_color;

    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;
    out_color_alpha = color;
    out_normal = encode_normal(view_normal);
    out_ss_vel = vec4(ss_vel, 0, 0);
    out_picking_color = picking_color;
}