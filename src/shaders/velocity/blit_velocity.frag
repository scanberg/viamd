#version 150 core

uniform sampler2D u_tex_depth;
uniform mat4 u_curr_clip_to_prev_clip_mat;
uniform vec4 u_jitter_uv;

in  vec2 tc;
out vec4 out_ss_vel;

void main() {
    float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
    // Background should output zero velocity so later passes (dilation, TAA, motion blur)
    // don't pick up undefined/previous-frame values.
    if (d == 1.0f) {
        out_ss_vel = vec4(0.0);
        return;
    }

    vec2 p_uv = tc;
    vec3 p_vs = vec3(tc, d);
    vec4 p_cs = vec4(p_vs * 2.0 - 1.0, 1.0); // [0, 1] -> [-1, 1]
    
	vec4 q_cs = u_curr_clip_to_prev_clip_mat * p_cs;
    vec2 q_uv = (q_cs.xy / max(q_cs.w, 1e-6)) * 0.5 + 0.5; // [-1, 1] -> [0, 1]
    q_uv = clamp(q_uv, vec2(0.0), vec2(1.0));

    vec2 ss_vel = (p_uv - q_uv) + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_ss_vel = vec4(ss_vel, 0, 0);
}