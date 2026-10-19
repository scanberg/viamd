#version 410 core

// One frame of the isosurface meshes (isomesh.frag) resolved and added to the running mean of the frames: frame
// n (from 0) enters with weight 1 / (n + 1), so after N frames the result is the plain average of N estimates.
//
// out_color: the radiance of the surfaces times their visibility, and their opacity times their visibility
//            (whose mean is the coverage, which the composite uses to correct the radiance)
// out_t:     the transmittance of the surfaces (exactly), and that times the absorption of what they enclose:
//            what reaches the eye of the scene behind them
// out_composite: the new mean as the composite (isomesh_utils.cpp) puts it over the caller's target, blended there
//            (ONE, ONE_MINUS_SRC_ALPHA): in the same pass, so the mean is not read again for it

#ifdef ISOMESH_MULTISAMPLE     // as it was (isomesh::Experiment_Multisampled)
uniform sampler2DMS u_tex_color;
uniform sampler2DMS u_tex_tod;
#else                           // a sample per pixel: plain textures, texelFetch of level s (0)
uniform sampler2D   u_tex_color;
uniform sampler2D   u_tex_tod;
#endif
uniform int         u_samples;
uniform sampler2D   u_tex_history_color;
uniform sampler2D   u_tex_history_t;
uniform float       u_weight;           // 1 / (n + 1); 1 discards the history

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_t;
layout(location = 2) out vec4 out_composite;

void main() {
    ivec2 px = ivec2(gl_FragCoord.xy);
    vec4 color = vec4(0.0);
    vec2 t = vec2(0.0);
    for (int s = 0; s < u_samples; ++s) {
        color += texelFetch(u_tex_color, px, s);
        vec2  tod = texelFetch(u_tex_tod, px, s).rg;
        float T = exp(tod.r);
        t += vec2(T, T * exp(-max(tod.g, 0.0)));
    }
    float inv = 1.0 / float(u_samples);
    color *= inv;
    t *= inv;

    if (u_weight >= 1.0) {
        out_color = color;
        out_t = vec4(t, 0.0, 0.0);
    } else {
        out_color = mix(texelFetch(u_tex_history_color, px, 0), color, u_weight);
        out_t = vec4(mix(texelFetch(u_tex_history_t, px, 0).rg, t, u_weight), 0.0, 0.0);
    }

    // As the composite: the radiance scaled to the exact coverage of the surfaces (Enderton's alpha correction), and
    // the coverage of the surfaces and what they enclose
    float coverage = 1.0 - out_t.r;
    vec3  L = (out_color.a > 1.0e-4) ? out_color.rgb * (coverage / out_color.a) : out_color.rgb;
    out_composite = vec4(L, 1.0 - out_t.g);
}
