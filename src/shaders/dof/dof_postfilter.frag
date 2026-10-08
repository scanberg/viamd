#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// 3x3 tent over the half-res gather result: removes the sparse-sample stipple at large CoC (cheap: 9 taps at 1/4 res).
// Only applied where the gather actually blurred (alpha / CoC > 0), so in-focus texels stay sharp.
uniform sampler2D u_tex;       // gather output
uniform sampler2D u_tex_coc;   // half-res prepass (coc in .a)
out vec4 out_frag;
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    vec4 c = texelFetch(u_tex, p, 0);
    float cc = abs(texelFetch(u_tex_coc, p, 0).a);
    if (cc < 1.0 && c.a <= 0.0) { out_frag = c; return; }
    vec4 sum = vec4(0); float wsum = 0.0;
    for (int y = -1; y <= 1; ++y) for (int x = -1; x <= 1; ++x) {
        float w = (x == 0 ? 2.0 : 1.0) * (y == 0 ? 2.0 : 1.0);
        sum += texelFetch(u_tex, clamp(p + ivec2(x, y), ivec2(0), s), 0) * w;
        wsum += w;
    }
    out_frag = sum / wsum;
}
