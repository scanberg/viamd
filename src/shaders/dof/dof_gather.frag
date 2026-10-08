#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Half-res scatter-as-gather with separate near / far accumulation (in the spirit of Sousa 2013, Jimenez 2014).
// Fixed sample count spread over the (dilated) tile max CoC, so cost is independent of blur size.
#ifndef TILE
#define TILE 8
#endif
#ifndef NUM_SAMPLES
#define NUM_SAMPLES 32
#endif
uniform sampler2D u_tex;        // half-res: rgb, signed coc (half-res px)
uniform sampler2D u_tex_tile;   // dilated tile max |coc|
uniform int u_frame;
out vec4 out_frag;              // rgb = blurred colour, a = near-field coverage (far field blends on full-res CoC)

const float GOLDEN_ANGLE = 2.39996323;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    vec4 c = texelFetch(u_tex, p, 0);
    float R = texelFetch(u_tex_tile, min(p / TILE, textureSize(u_tex_tile, 0) - 1), 0).r;
    if (R < 0.5) { out_frag = vec4(c.rgb, 0.0); return; }

    float cc = abs(c.a);
    // interleaved gradient noise rotation (Jimenez 2014); the postfilter / TAA smooth what is left
    float rot = 6.2831853 * fract(52.9829189 * fract(dot(vec2(p) + float(u_frame) * 5.588238, vec2(0.06711056, 0.00583715))));
    vec2 cs = vec2(cos(rot), sin(rot));

    // far field (behind the focus plane): footprint limited by the centre's own CoC, so a blurred background never
    // bleeds over a sharper foreground.
    vec3  far_sum = c.rgb;
    float far_w = 1.0;
    // near field: samples in front of the focus plane, weighted by coverage of the centre (energy spread over pi*coc^2)
    vec3  near_sum = vec3(0);
    float near_w = 0.0;
    float near_cov = 0.0;
    float area = R * R / float(NUM_SAMPLES);   // disk area per sample / pi

    for (int i = 0; i < NUM_SAMPLES; ++i) {
        float r = sqrt((float(i) + 0.5) / float(NUM_SAMPLES)) * R;
        float a = float(i) * GOLDEN_ANGLE;
        vec2 d = vec2(cos(a), sin(a));   // constant after unrolling
        vec2 o = vec2(d.x * cs.x - d.y * cs.y, d.x * cs.y + d.y * cs.x) * r;
        vec4 sm = texelFetch(u_tex, clamp(p + ivec2(round(o)), ivec2(0), s), 0);
        float sc = abs(sm.a);
        if (sm.a < 0.0) {
            float w = clamp(sc - r + 1.0, 0.0, 1.0);
            near_sum += sm.rgb * w;
            near_w += w;
            near_cov += w * area / max(sc * sc, 1.0);
        } else {
            float w = clamp(min(sc, cc) - r + 1.0, 0.0, 1.0);
            far_sum += sm.rgb * w;
            far_w += w;
        }
    }
    vec3 far_col = far_sum / far_w;
    float near_a = clamp(near_cov, 0.0, 1.0);
    if (c.a < 0.0) near_a = max(near_a, smoothstep(0.5, 1.5, cc));
    vec3 col = near_w > 0.0 ? mix(far_col, near_sum / near_w, near_a) : far_col;
    out_frag = vec4(col, near_a);
}
