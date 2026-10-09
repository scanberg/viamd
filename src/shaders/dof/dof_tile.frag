#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather (far) / near bands -> fill -> composite.
//
// u_dilate = 0: per TILE x TILE block of half-res pixels:
//               .r largest near-field |CoC| (in front of the focus plane), .g largest far-field CoC,
//               .b smallest CoC of the near-field content (as the near bands classify it, see dof_near)
// u_dilate = 1: .r gather radius of the far field, .g reach of the near field, .b 3x3 min of .b above
//               A near-field (foreground) blur spreads over whatever is behind it, out to its own CoC. The search
//               covers u_reach_tiles tiles in each direction (enough for the CoC clamp) and only takes a neighbour's
//               near CoC if it actually reaches this tile.
//               The smallest near CoC around a tile tells the finest near bands whether they have anything to gather
//               there: the content they hold reaches less than a tile.
#ifndef TILE
#define TILE 8
#endif
uniform sampler2D u_tex;            // u_dilate = 0: half-res RGBA (signed coc in .a); 1: output of pass 0
uniform sampler2D u_tex_near_src;   // u_dilate = 0: (rgb, 1) * near energy
uniform sampler2D u_tex_coc_src;    // u_dilate = 0: .r = coc * near energy
uniform int u_dilate;
uniform int u_reach_tiles;          // ceil(max CoC / TILE), in tiles
out vec4 out_frag;

const float NO_NEAR = 65000.0;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    if (u_dilate == 0) {
        float near = 0.0, far = 0.0, near_min = NO_NEAR;
        for (int y = 0; y < TILE; ++y) {
            for (int x = 0; x < TILE; ++x) {
                ivec2 q = min(p * TILE + ivec2(x, y), s);
                float coc = texelFetch(u_tex, q, 0).a;
                near = max(near, -coc);
                far  = max(far,   coc);
                float e = texelFetch(u_tex_near_src, q, 0).a;
                if (e > 1.0e-8) near_min = min(near_min, texelFetch(u_tex_coc_src, q, 0).r / e);
            }
        }
        out_frag = vec4(near, far, near_min, 0.0);
    } else {
        vec3 own = texelFetch(u_tex, p, 0).rgb;
        float rn = own.r;
        float near_min = NO_NEAR;
        int k = u_reach_tiles;
        for (int y = -k; y <= k; ++y) {
            for (int x = -k; x <= k; ++x) {
                ivec2 q = p + ivec2(x, y);
                if (any(lessThan(q, ivec2(0))) || any(greaterThan(q, s))) continue;
                vec3 t = texelFetch(u_tex, q, 0).rgb;
                // distance between the closest points of the two tiles, in half-res pixels
                vec2 gap = vec2(max(abs(x) - 1, 0), max(abs(y) - 1, 0)) * float(TILE);
                if (t.r > rn && dot(gap, gap) < t.r * t.r) rn = t.r;
                if (abs(x) <= 1 && abs(y) <= 1) near_min = min(near_min, t.b);
            }
        }
        out_frag = vec4(max(rn, own.g), rn, near_min, 0.0);
    }
}
