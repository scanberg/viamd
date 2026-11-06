#version 410 core

// Part of the half-res depth of field: near field (in front of the focus plane), one pass per band L, coarsest first.
//
// Every foreground texel spreads its energy e = m / coc^2 over its CoC disk; the coverage of a pixel is the energy that
// reaches it (thin lens: 1/2 at the silhouette of a large foreground, 1 within it, 0 outside, a foreground thinner than
// its CoC never becomes opaque). The near field is split by CoC into bands an octave apart (band L: CoC around
// 2^(L + 1/2) half-res px, soft shares per half-res texel, see dof_near_down). Band L is gathered exactly at pyramid
// level L, where its CoC is 0.7 to 2.8 texels: every texel within reach, no sampling pattern. The bands are summed from
// the coarsest up, each level adding its band to the smoothly upsampled sum of the levels above it.
//
// A single gather over the largest CoC around (a tile radius) has to space its samples for that CoC. A less blurred
// structure next to a strongly blurred one is then sampled far too coarsely, and coarse texels that mix the two have a
// mean CoC that fits neither: the less blurred one comes out in blocks of the sampling grid.
//
// Output (premultiplied, unnormalised): rgb = sum of colour * coverage, a = sum of coverage (can exceed 1, the
// composite clamps it).
#ifndef TILE
#define TILE 8
#endif
uniform sampler2D u_tex_src;    // L = 0: half res (rgb, 1) * e; L >= 1: band L at level L (dof_near_down)
uniform sampler2D u_tex_coc;    // .r = coc * e, matching u_tex_src
uniform sampler2D u_tex_up;     // near field of the levels above (L + 1 and up), bilinear
uniform sampler2D u_tex_tile;   // tile pass 1: .b = smallest near CoC around the tile
uniform int u_level;
uniform int u_src_lod;          // level of u_tex_src / u_tex_coc holding this band
uniform int u_top;              // coarsest band, it takes everything above it
uniform int u_has_up;
out vec4 out_frag;

const float PI = 3.14159265;
// Band L holds CoCs up to 2^(L + 1.5) (2.83 texels at its level), the reach is softened by half a texel
const float REACH = 2.83 + 0.5;

float band_weight(float c, int L, int top) {
    float x = log2(max(c, 1.0e-3)) - (float(L) + 0.5);
    if ((L == 0 && x < 0.0) || (L == top && x > 0.0)) return 1.0;
    return clamp(1.0 - abs(x), 0.0, 1.0);
}

void main() {
    ivec2 t = ivec2(gl_FragCoord.xy);
    float scale = exp2(float(u_level));     // half-res px per texel of this level
    vec4  acc = vec4(0);

    // The two finest levels cover most pixels: skip them where nothing of their band is around (it reaches less than
    // a tile, the tile pass knows the smallest near CoC of the 3x3 tiles around)
    bool run = true;
    if (u_level <= 1 && u_level < u_top) {
        ivec2 tile = min((t << u_level) / TILE, textureSize(u_tex_tile, 0) - 1);
        run = texelFetch(u_tex_tile, tile, 0).b < exp2(float(u_level) + 1.5);
    }
    if (run) {
        ivec2 smax = textureSize(u_tex_src, u_src_lod) - 1;
        for (int y = -3; y <= 3; ++y) {
            for (int x = -3; x <= 3; ++x) {
                float d = length(vec2(x, y));
                if (d > REACH) continue;
                ivec2 q = t + ivec2(x, y);
                if (any(lessThan(q, ivec2(0))) || any(greaterThan(q, smax))) continue;
                vec4 s = texelFetch(u_tex_src, q, u_src_lod);
                if (s.a <= 1.0e-8) continue;
                float c = texelFetch(u_tex_coc, q, u_src_lod).r / s.a;
                if (u_level == 0) s *= band_weight(c, 0, u_top);
                acc += s * clamp(c / scale - d + 0.5, 0.0, 1.0);
            }
        }
        // Each texel holds the mean energy density of scale^2 half-res texels
        acc *= scale * scale / PI;
    }

    if (u_has_up == 1) {
        // Bilinear taps a quarter texel (of the level above) off centre: a slightly wider tent than plain bilinear
        // upsampling, so the coarse levels do not show their grid
        vec2 inv = 1.0 / vec2(textureSize(u_tex_up, 0));
        vec2 uv  = (vec2(t) + 0.5) * 0.5 * inv;
        vec2 d   = 0.25 * inv;
        acc += 0.25 * (textureLod(u_tex_up, uv + vec2(-d.x, -d.y), 0.0) +
                       textureLod(u_tex_up, uv + vec2( d.x, -d.y), 0.0) +
                       textureLod(u_tex_up, uv + vec2(-d.x,  d.y), 0.0) +
                       textureLod(u_tex_up, uv + vec2( d.x,  d.y), 0.0));
    }
    out_frag = acc;
}
