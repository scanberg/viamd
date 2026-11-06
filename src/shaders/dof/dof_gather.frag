#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather (far) / near bands -> fill -> composite.
// Far field (what is at / behind the focus plane), scatter-as-gather evaluated as a quadrature rather than a point
// sampled estimate: the disk of the centre's own CoC is integrated with a fixed pattern (no per-pixel or per-frame
// rotation), every sample reading the mip level whose texels are as large as the spacing between samples, so it stands
// for the average over its share of the disk. The result is smooth in space and time on its own; a noisy estimate
// (rotated per frame for TAA to average) flickers instead, since TAA clips the history to the 3x3 neighbourhood of the
// current frame, which cannot hold the variance of noise correlated over several pixels.
//
// A sample contributes where its own CoC reaches the centre too, so a blurred background never bleeds over a sharper
// object in front of it.
//
// Where the centre belongs to the foreground (near field, dof_near), what is behind it is unknown here: it is left as
// a hole (alpha 0), which dof_fill fills from the far layer around.
//
// Output (consumed by dof_fill): rgb * a, a = 1 - how much the centre belongs to the near field.
#ifndef TILE
#define TILE 8
#endif
#ifndef FAR_SAMPLES
#define FAR_SAMPLES 32
#endif
uniform sampler2D u_tex;            // half-res: rgb, signed coc (half-res px)
uniform sampler2D u_tex_tile;       // .r = gather radius per tile (half-res px)
uniform sampler2D u_tex_far_src;    // mipmapped (rgb, 1) * far CoC
uniform sampler2D u_tex_coc_src;    // mipmapped (-, -, coc * far CoC, -)
out vec4 out_far;

const float GOLDEN_ANGLE = 2.39996323;
const float SQRT_PI      = 1.7724539;

// Sample i of n on the unit disk (Vogel spiral)
vec2 disk(int i, int n) {
    float a = float(i) * GOLDEN_ANGLE;
    return vec2(cos(a), sin(a)) * sqrt((float(i) + 0.5) / float(n));
}

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec4  c = texelFetch(u_tex, p, 0);
    float R = texelFetch(u_tex_tile, min(p / TILE, textureSize(u_tex_tile, 0) - 1), 0).r;
    if (R < 0.5) {
        out_far = vec4(c.rgb, 1.0);
        return;
    }

    vec3  far = c.rgb;
    float cc  = c.a;
    if (cc > 1.0) {
        vec2  inv_size = 1.0 / vec2(textureSize(u_tex_far_src, 0));
        vec2  uv0      = (vec2(p) + 0.5) * inv_size;
        const int NF = FAR_SAMPLES;
        float spacing = max(cc * SQRT_PI / sqrt(float(NF)), 1.0);   // cc * sqrt(pi / NF)
        float lod     = log2(spacing);
        vec4  acc     = vec4(0);
        for (int i = 0; i < NF; ++i) {
            vec2  o  = disk(i, NF) * cc;
            vec2  uv = uv0 + o * inv_size;
            vec4  s  = textureLod(u_tex_far_src, uv, lod);
            if (s.a <= 1.0e-4) continue;
            float sc    = textureLod(u_tex_coc_src, uv, lod).b / s.a;
            float reach = clamp((min(sc, cc) - length(o)) / spacing + 0.5, 0.0, 1.0);
            acc += s * reach;
        }
        // The centre counts as one sample of its own CoC: where little or nothing reaches it the result goes smoothly
        // to the centre's colour instead of switching to it.
        far = (acc.rgb + c.rgb * cc) / (acc.a + cc);
    }

    float valid = 1.0 - smoothstep(0.5, 1.5, -c.a);
    out_far = vec4(far * valid, valid);
}
