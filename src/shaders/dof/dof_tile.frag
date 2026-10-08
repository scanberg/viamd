#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Max |CoC| per tile (TILE x TILE half-res px). Second invocation with u_dilate = 1 takes the 3x3 neighbourhood max,
// so near-field blur can spread into neighbouring in-focus tiles.
#ifndef TILE
#define TILE 8
#endif
uniform sampler2D u_tex;   // u_dilate = 0: half-res RGBA (coc in .a); u_dilate = 1: tile texture (R)
uniform int u_dilate;
out float out_frag;
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 s = textureSize(u_tex, 0) - 1;
    float m = 0.0;
    if (u_dilate == 0) {
        for (int y = 0; y < TILE; ++y) for (int x = 0; x < TILE; ++x)
            m = max(m, abs(texelFetch(u_tex, min(p * TILE + ivec2(x, y), s), 0).a));
    } else {
        for (int y = -1; y <= 1; ++y) for (int x = -1; x <= 1; ++x)
            m = max(m, texelFetch(u_tex, clamp(p + ivec2(x, y), ivec2(0), s), 0).r);
    }
    out_frag = m;
}
