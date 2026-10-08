#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Small disk filter (centre + 8 bilinear taps) at the gather's own sample spacing, which is what the sparse sampling
// leaves behind as stipple. The near layer is premultiplied and smooth by nature, so it is filtered everywhere it is
// non-zero. The far layer only takes taps whose own CoC is at least the tap distance, so sharp foreground does not leak
// into a blurred background.
#ifndef TILE
#define TILE 8
#endif
#ifndef DOF_MAX_SAMPLES
#define DOF_MAX_SAMPLES 64
#endif
uniform sampler2D u_tex_far;    // gather output, location 0 (LINEAR)
uniform sampler2D u_tex_near;   // gather output, location 1 (LINEAR, premultiplied)
uniform sampler2D u_tex_coc;    // half-res prepass (signed coc in .a)
uniform sampler2D u_tex_tile;   // .r = gather radius per tile
layout(location = 0) out vec4 out_far;
layout(location = 1) out vec4 out_near;

int sample_count(float R) { return clamp(int(ceil(R * 1.6)), 16, DOF_MAX_SAMPLES); }

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec2  size = vec2(textureSize(u_tex_far, 0));
    vec4  far  = texelFetch(u_tex_far, p, 0);
    vec4  near = texelFetch(u_tex_near, p, 0);
    float R = texelFetch(u_tex_tile, min(p / TILE, textureSize(u_tex_tile, 0) - 1), 0).r;
    if (R < 0.5) {
        out_far = far;
        out_near = near;
        return;
    }

    float spacing = R * 1.7724539 / sqrt(float(sample_count(R)));   // R * sqrt(pi / n)
    float rho = clamp(0.5 * spacing, 1.0, 8.0);
    vec2  uv0 = (vec2(p) + 0.5) / size;

    vec4  near_sum = near * 2.0;
    float near_w = 2.0;
    vec3  far_sum = far.rgb * 2.0;
    float far_w = 2.0;
    for (int i = 0; i < 8; ++i) {
        float a = (float(i) + 0.5) * 0.78539816;
        vec2  o = vec2(cos(a), sin(a)) * rho;
        vec2  uv = uv0 + o / size;
        near_sum += texture(u_tex_near, uv);
        near_w += 1.0;
        float sc = abs(texelFetch(u_tex_coc, clamp(ivec2(vec2(p) + o + 0.5), ivec2(0), ivec2(size) - 1), 0).a);
        float w = clamp(sc - rho + 1.0, 0.0, 1.0);
        far_sum += texture(u_tex_far, uv).rgb * w;
        far_w += w;
    }
    out_far  = vec4(far_sum / far_w, 0.0);
    out_near = near_sum / near_w;
}
