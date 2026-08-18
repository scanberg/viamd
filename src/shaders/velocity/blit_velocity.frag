#version 150 core

uniform sampler2D u_tex_depth;
uniform mat4 u_curr_clip_to_prev_clip_mat;
uniform vec4 u_jitter_uv;

out vec4 out_ss_vel;

void main() {
    float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
    // Background should output zero velocity so later passes (dilation, TAA, motion blur)
    // don't pick up undefined/previous-frame values.
    if (d == 1.0f) {
        out_ss_vel = vec4(0.0);
        return;
    }

    // gl_FragCoord.xy is already at the pixel centre -- the fragment covering
    // texel (i,j) reports (i+0.5, j+0.5) -- so dividing by the resolution gives
    // that texel's centre in UV directly. Adding another half pixel first put
    // the sample on the corner between texels instead.
    //
    // The offset propagates into p_cs and therefore into the reprojected q_uv,
    // so the velocity was not merely biased, it was evaluated half a texel away
    // from the fragment it was written for. Small, uniform, and permanent:
    // every TAA resolve fetched history half a texel off and blended it in,
    // which reads as the image never quite converging to sharp.
    vec2 p_uv = gl_FragCoord.xy / vec2(textureSize(u_tex_depth, 0));
    vec3 p_vs = vec3(p_uv, d);
    vec4 p_cs = vec4(p_vs * 2.0 - 1.0, 1.0); // [0, 1] -> [-1, 1]
    
	vec4 q_cs = u_curr_clip_to_prev_clip_mat * p_cs;
    vec2 q_uv = (q_cs.xy / max(q_cs.w, 1e-6)) * 0.5 + 0.5; // [-1, 1] -> [0, 1]
    q_uv = clamp(q_uv, vec2(0.0), vec2(1.0));

    vec2 ss_vel = (p_uv - q_uv) + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_ss_vel = vec4(ss_vel, 0, 0);
}