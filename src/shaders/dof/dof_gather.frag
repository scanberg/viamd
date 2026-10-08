#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Half-res scatter-as-gather with separate near / far layers (in the spirit of Sousa 2013, Jimenez 2014).
// The sample count grows with the gather radius (16..DOF_MAX_SAMPLES), so small blurs stay cheap and large ones are
// not reduced to a sparse stipple.
//
// Outputs (both consumed by dof_postfilter, then dof_composite):
//   out_far  rgb = far-field blur (what is at / behind the focus plane), a unused
//   out_near rgb = near-field colour * coverage (premultiplied), a = coverage of this pixel by the foreground blur
// Keeping the layers apart lets the composite apply the near coverage exactly once.
#ifndef TILE
#define TILE 8
#endif
#ifndef DOF_MAX_SAMPLES
#define DOF_MAX_SAMPLES 64
#endif
uniform sampler2D u_tex;        // half-res: rgb, signed coc (half-res px)
uniform sampler2D u_tex_tile;   // .r = gather radius per tile (half-res px)
uniform int u_frame;
layout(location = 0) out vec4 out_far;
layout(location = 1) out vec4 out_near;

const float GOLDEN_ANGLE = 2.39996323;

int sample_count(float R) { return clamp(int(ceil(R * 1.6)), 16, DOF_MAX_SAMPLES); }

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    vec4 c = texelFetch(u_tex, p, 0);
    float R = texelFetch(u_tex_tile, min(p / TILE, textureSize(u_tex_tile, 0) - 1), 0).r;
    if (R < 0.5) {
        out_far  = vec4(c.rgb, 0.0);
        out_near = vec4(0.0);
        return;
    }

    float cc = abs(c.a);
    int   n  = sample_count(R);
    // interleaved gradient noise rotation (Jimenez 2014); the postfilter / TAA smooth what is left
    float rot = 6.2831853 * fract(52.9829189 * fract(dot(vec2(p) + float(u_frame) * 5.588238, vec2(0.06711056, 0.00583715))));
    vec2  cs  = vec2(cos(rot), sin(rot));

    // far field (at / behind the focus plane): footprint limited by the centre's own CoC, so a blurred background never
    // bleeds over a sharper foreground.
    vec3  far_sum = c.rgb;
    float far_w   = 1.0;
    // near field: every sample in front of the focus plane whose CoC reaches the centre. Its energy is spread over
    // pi*coc^2, each sample stands for pi*R^2/n of area -> coverage += (R^2/n) / coc^2.
    vec3  near_sum = vec3(0);
    float near_w   = 0.0;
    float near_cov = 0.0;
    float area     = R * R / float(n);

    for (int i = 0; i < DOF_MAX_SAMPLES; ++i) {
        if (i >= n) break;
        float r = sqrt((float(i) + 0.5) / float(n)) * R;
        float a = float(i) * GOLDEN_ANGLE;
        vec2  d = vec2(cos(a), sin(a));
        vec2  o = vec2(d.x * cs.x - d.y * cs.y, d.x * cs.y + d.y * cs.x) * r;
        vec4  sm = texelFetch(u_tex, clamp(p + ivec2(round(o)), ivec2(0), s), 0);
        float sc = abs(sm.a);
        if (sm.a < 0.0) {
            float w = clamp(sc - r + 1.0, 0.0, 1.0);
            near_sum += sm.rgb * w;
            near_w   += w;
            near_cov += w * area / max(sc * sc, 1.0);
        } else {
            float w = clamp(min(sc, cc) - r + 1.0, 0.0, 1.0);
            far_sum += sm.rgb * w;
            far_w   += w;
        }
    }

    float near_a = clamp(near_cov, 0.0, 1.0);
    if (c.a < 0.0) near_a = max(near_a, smoothstep(0.5, 1.5, cc));   // the centre itself is foreground
    vec3 near_col = near_w > 0.0 ? near_sum / near_w : c.rgb;

    out_far  = vec4(far_sum / far_w, 0.0);
    out_near = vec4(near_col * near_a, near_a);
}
