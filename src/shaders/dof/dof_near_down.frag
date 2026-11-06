#version 410 core

// Part of the half-res depth of field: the source of near band L >= 1 at pyramid level L (see dof_near.frag).
// Box average over the 2^L x 2^L half-res texels of the band's share of every texel, so the separation into bands is
// made per half-res texel and a coarse texel never holds content of very different CoC.
uniform sampler2D u_tex_near_src;   // half res: (rgb, 1) * near energy
uniform sampler2D u_tex_coc_src;    // half res: .r = coc * near energy
uniform int u_level;                // band / level L
uniform int u_top;                  // coarsest band, it takes everything above it
layout(location = 0) out vec4 out_src;
layout(location = 1) out vec4 out_coc;

// Share of a near-field texel of CoC c (half-res px) in band L: tents over log2(c) centred on 2^(L + 1/2), so a band
// holds CoCs within a factor 2 of its centre and the shares of all bands sum to 1.
float band_weight(float c, int L, int top) {
    float x = log2(max(c, 1.0e-3)) - (float(L) + 0.5);
    if ((L == 0 && x < 0.0) || (L == top && x > 0.0)) return 1.0;
    return clamp(1.0 - abs(x), 0.0, 1.0);
}

void main() {
    ivec2 t = ivec2(gl_FragCoord.xy);
    int   n = 1 << u_level;
    ivec2 smax = textureSize(u_tex_near_src, 0) - 1;
    vec4  acc = vec4(0);
    float acc_coc = 0.0;
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            ivec2 q = t * n + ivec2(x, y);
            if (q.x > smax.x || q.y > smax.y) continue;
            vec4 s = texelFetch(u_tex_near_src, q, 0);
            if (s.a <= 1.0e-8) continue;
            float sc = texelFetch(u_tex_coc_src, q, 0).r;
            float w  = band_weight(sc / s.a, u_level, u_top);
            acc     += s * w;
            acc_coc += sc * w;
        }
    }
    float inv = 1.0 / float(n * n);
    out_src = acc * inv;
    out_coc = vec4(acc_coc * inv, 0.0, 0.0, 0.0);
}
