#version 410 core

// Joint bilateral upsample of the half-res AO to full res; output is multiplied into the lit color (blend ZERO, SRC_COLOR).
// Two textureGather calls fetch the 2x2 (visibility, z) footprint.
#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
uniform sampler2D u_tex_ao;            // RG16F half res (visibility, linear depth), CLAMP_TO_EDGE
uniform sampler2D u_tex_linear_depth;  // full-res linear depth (level 0)
uniform float u_px_scale;
uniform float u_z_max;
out vec4 out_frag;
void main() {
    ivec2 q = ivec2(gl_FragCoord.xy);
    float z = texelFetch(u_tex_linear_depth, q, 0).r;
    if (z >= u_z_max) { out_frag = vec4(1); return; }
    vec2 half_size = vec2(textureSize(u_tex_ao, 0));
    vec2 hq = (vec2(q) + 0.5) * 0.5 - 0.5;
    vec2 base = floor(hq);
    vec2 f = hq - base;
    vec2 uv = (base + 1.0) / half_size;
    vec4 ao = textureGather(u_tex_ao, uv, 0);   // (0,1) (1,1) (1,0) (0,0)
    vec4 sz = textureGather(u_tex_ao, uv, 1);
    vec4 wb = vec4((1.0 - f.x) * f.y, f.x * f.y, f.x * (1.0 - f.y), (1.0 - f.x) * (1.0 - f.y));
#if AO_PERSPECTIVE
    float fp = u_px_scale * z;
#else
    float fp = u_px_scale;
#endif
    vec4 wz = 1.0 / (1.0 + abs(sz - z) / (2.0 * fp));
    vec4 w = wb * wz * wz + 1e-5;
    float vis = dot(ao, w) / dot(w, vec4(1));
    out_frag = vec4(vis);
}
