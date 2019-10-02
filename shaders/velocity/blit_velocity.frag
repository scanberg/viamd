#version 150 core

uniform sampler2D u_tex_depth;
uniform mat4 u_curr_clip_to_prev_clip_mat;
uniform vec4 u_jitter_uv;

in  vec2 tc;
out vec4 out_ss_vel;

void main() {
    float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
    if (d == 1.0f) discard;
    //d = 1.0f;

    vec2 p_uv = tc;
    vec3 p_vs = vec3(tc, d);
    vec4 p_cs = vec4(p_vs * 2.0 - 1.0, 1.0); // [0, 1] -> [-1, 1]
    
	vec4 q_cs = u_curr_clip_to_prev_clip_mat * p_cs;
    vec2 q_uv = (q_cs.xy / q_cs.w) * 0.5 + 0.5; // [-1, 1] -> [0, 1]

    vec2 ss_vel = (p_uv - q_uv) + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_ss_vel = vec4(ss_vel, 0, 0);
}