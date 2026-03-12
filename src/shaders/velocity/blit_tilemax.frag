#version 150 core

#pragma optionNV(unroll all)

uniform sampler2D u_tex_vel;
uniform sampler2D u_tex_linear_depth;

#ifndef TILE_SIZE
#define TILE_SIZE 20
#endif

#ifndef MAX_BLUR_PX
#define MAX_BLUR_PX 48
#endif

in vec2 tc;
out vec4 out_max_tile_vel;

void main() {
    ivec2 tileBase = ivec2(gl_FragCoord.xy) * TILE_SIZE; // output is width/TILE_SIZE

    float minD = 1e30;
    for (int j = 0; j < TILE_SIZE; ++j)
    for (int i = 0; i < TILE_SIZE; ++i) {
        ivec2 p = tileBase + ivec2(i,j);
        float d = texelFetch(u_tex_linear_depth, p, 0).x;
        if (d < minD) minD = d;
    }

    // allow a small thickness so you don’t drop the surface due to quantization
    float depthEps = 0.01 * minD; // tweak; could be constant too

    vec2 mv = vec2(0);
    float mv2 = 0;
    for (int j = 0; j < TILE_SIZE; ++j)
    for (int i = 0; i < TILE_SIZE; ++i) {
        ivec2 p = tileBase + ivec2(i,j);
        float d = texelFetch(u_tex_linear_depth, p, 0).x;
        if (d > minD + depthEps) continue; // reject farther stuff

        vec2 v = texelFetch(u_tex_vel, p, 0).xy;
        float v2 = dot(v,v);
        if (v2 > mv2) { mv = v; mv2 = v2; }
    }

    out_max_tile_vel = vec4(mv, 0.0, 0.0);
}