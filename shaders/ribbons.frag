#version 330 core

uniform vec4 u_jitter_uv;

in Fragment {
    smooth vec4 color;
    smooth vec4 view_normal;
    smooth vec4 curr_clip_coord;
    smooth vec4 prev_clip_coord;
    flat vec4 picking_color;
} in_frag;

layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_ss_vel;
layout(location = 3) out vec4 out_emission;
layout(location = 4) out vec4 out_post_tonemap;
layout(location = 5) out vec4 out_picking_color;

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
	vec2 curr_ndc = in_frag.curr_clip_coord.xy / in_frag.curr_clip_coord.w;
	vec2 prev_ndc = in_frag.prev_clip_coord.xy / in_frag.prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color_alpha = in_frag.color;
    out_normal = encode_normal(normalize(in_frag.view_normal.xyz));
    out_ss_vel = vec4(ss_vel, 0, 0);
    out_picking_color = in_frag.picking_color;
}