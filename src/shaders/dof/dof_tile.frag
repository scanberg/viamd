#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather -> composite.
//
// u_dilate = 0: per TILE x TILE block of half-res pixels, the largest near-field |CoC| (in front of the focus plane)
//               and the largest far-field CoC.
// u_dilate = 1: gather radius per tile (.r) and reach of the near field (.g). A near-field (foreground) blur spreads
//               over whatever is behind it, out to its own CoC, so every tile within that reach must gather at least
//               that far. The search covers u_reach_tiles tiles in each direction (enough for the CoC clamp) and only
//               takes a neighbour's near CoC if it actually reaches this tile. The far field only blurs within its own
//               CoC: no dilation needed. A fixed 3x3 dilation cut foreground blur off at tile edges as soon as the CoC
//               exceeded one tile, which shows up as axis aligned stair steps.
// u_dilate = 2: 3x3 max of the near-field reach, which the gather interpolates bilinearly (see there).
#ifndef TILE
#define TILE 8
#endif
uniform sampler2D u_tex;        // u_dilate = 0: half-res RGBA (signed coc in .a); 1: tile texture (near, far); 2: output of 1
uniform int u_dilate;
uniform int u_reach_tiles;      // ceil(max CoC / TILE), in tiles
out vec2 out_frag;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    if (u_dilate == 0) {
        float near = 0.0, far = 0.0;
        for (int y = 0; y < TILE; ++y) {
            for (int x = 0; x < TILE; ++x) {
                float coc = texelFetch(u_tex, min(p * TILE + ivec2(x, y), s), 0).a;
                near = max(near, -coc);
                far  = max(far,   coc);
            }
        }
        out_frag = vec2(near, far);
    } else if (u_dilate == 2) {
        // 3x3 max of the near-field reach. Interpolated bilinearly it is a continuous radius that is still at least the
        // reach of the tile it is sampled in (every value it interpolates includes that tile).
        float r = 0.0;
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                r = max(r, texelFetch(u_tex, clamp(p + ivec2(x, y), ivec2(0), s), 0).g);
            }
        }
        out_frag = vec2(r, 0.0);
    } else {
        vec2 own = texelFetch(u_tex, p, 0).rg;
        float rn = own.r;
        int k = u_reach_tiles;
        for (int y = -k; y <= k; ++y) {
            for (int x = -k; x <= k; ++x) {
                ivec2 q = p + ivec2(x, y);
                if (any(lessThan(q, ivec2(0))) || any(greaterThan(q, s))) continue;
                float near = texelFetch(u_tex, q, 0).r;
                // distance between the closest points of the two tiles, in half-res pixels
                vec2 gap = vec2(max(abs(x) - 1, 0), max(abs(y) - 1, 0)) * float(TILE);
                if (near > rn && dot(gap, gap) < near * near) rn = near;
            }
        }
        out_frag = vec2(max(rn, own.g), rn);   // .r: gather radius, .g: reach of the near field
    }
}
