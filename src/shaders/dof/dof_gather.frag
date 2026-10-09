#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather -> composite.
// Scatter-as-gather with separate near / far layers (in the spirit of Sousa 2013, Jimenez 2014), evaluated as a
// quadrature rather than a point sampled estimate:
//
// Each layer is integrated over a disk with a fixed pattern (no per-pixel or per-frame rotation), every sample reading
// the mip level whose texels are as large as the spacing between samples, so it stands for the average over its share
// of the disk. The result is smooth in space and time on its own. A noisy estimate (rotated per frame for TAA to
// average) flickers instead: TAA clips the history to the 3x3 neighbourhood of the current frame, which cannot hold
// the variance of noise that is correlated over several pixels, so the noise passes through. That is what made the
// region around a blurred foreground flicker.
//
// The radii change continuously from pixel to pixel (the far field uses the centre's own CoC, the near field an
// interpolated reach) and the sample counts are fixed, otherwise the quadrature, which differs slightly with them,
// would show the tiles.
//
// The near-field coverage is the thin-lens one: 1/2 at the silhouette of a large foreground, rising to 1 within it and
// falling to 0 outside, over its CoC; a foreground thinner than its CoC never becomes opaque. Where it is partially
// transparent, the far layer of a foreground pixel holds what is behind it, filled in from the background around.
//
// Outputs (consumed by dof_composite):
//   out_far  rgb = far-field blur (what is at / behind the focus plane), a unused
//   out_near rgb = near-field colour * coverage (premultiplied), a = coverage of this pixel by the foreground blur
#ifndef TILE
#define TILE 8
#endif
#ifndef FAR_SAMPLES
#define FAR_SAMPLES 32
#endif
#ifndef NEAR_SAMPLES
#define NEAR_SAMPLES 32
#endif
uniform sampler2D u_tex;            // half-res: rgb, signed coc (half-res px)
uniform sampler2D u_tex_tile;       // .r = gather radius per tile (half-res px)
uniform sampler2D u_tex_tile_near;  // .r = 3x3 max of the near-field reach per tile (LINEAR)
uniform sampler2D u_tex_near_src;   // mipmapped (rgb, 1) * near energy density
uniform sampler2D u_tex_far_src;    // mipmapped (rgb, 1) * far CoC
uniform sampler2D u_tex_coc_src;    // mipmapped (coc * near energy, -, coc * far CoC, -)
layout(location = 0) out vec4 out_far;
layout(location = 1) out vec4 out_near;

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
        out_far  = vec4(c.rgb, 0.0);
        out_near = vec4(0.0);
        return;
    }

    vec2 inv_size = 1.0 / vec2(textureSize(u_tex_far_src, 0));
    vec2 uv0      = (vec2(p) + 0.5) * inv_size;

    // Far field (at / behind the focus plane): the disk of the centre's own CoC. A sample contributes where its own
    // CoC reaches the centre too, so a blurred background never bleeds over a sharper object in front of it.
    // At a foreground centre, the background around it within the foreground's CoC (what shows through its blur).
    vec3  far = c.rgb;
    float cc  = abs(c.a);
    if (cc > 1.0) {
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
        // The centre counts as one sample of its own CoC: where little or nothing reaches it (a foreground centre far
        // from the background) the result goes smoothly to the centre's colour instead of switching to it.
        far = (acc.rgb + c.rgb * cc) / (acc.a + cc);
    }

    // Near field (in front of the focus plane): every foreground texel whose CoC reaches the centre. Its energy is
    // spread over pi * coc^2, each sample stands for pi * RN^2 / NN of area -> coverage += (RN^2 / NN) * (m / coc^2).
    vec4  near_acc = vec4(0);
    vec2  tile_uv = (vec2(p) + 0.5) / (float(TILE) * vec2(textureSize(u_tex_tile_near, 0)));
    float RN      = texture(u_tex_tile_near, tile_uv).r;
    if (RN >= 0.5) {
        const int NN = NEAR_SAMPLES;
        float spacing = max(RN * SQRT_PI / sqrt(float(NN)), 1.0);   // RN * sqrt(pi / NN)
        float lod     = log2(spacing);
        for (int i = 0; i < NN; ++i) {
            vec2  o  = disk(i, NN) * RN;
            vec2  uv = uv0 + o * inv_size;
            vec4  s  = textureLod(u_tex_near_src, uv, lod);
            if (s.a <= 1.0e-8) continue;
            float sc    = textureLod(u_tex_coc_src, uv, lod).r / s.a;
            float reach = clamp((sc - length(o)) / spacing + 0.5, 0.0, 1.0);
            near_acc += s * reach;
        }
        near_acc *= RN * RN / float(NN);
    }

    float near_a   = clamp(near_acc.a, 0.0, 1.0);
    vec3  near_col = near_acc.a > 1.0e-6 ? near_acc.rgb / near_acc.a : c.rgb;

    out_far  = vec4(far, 0.0);
    out_near = vec4(near_col * near_a, near_a);
}
