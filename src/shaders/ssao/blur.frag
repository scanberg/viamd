#version 410 core

// Separable bilateral blur at half resolution.
// Box kernel with half-weight end taps: with BLUR_RADIUS = 2*n it integrates exactly n periods of the 4x4 sample
// pattern, so the interleaving structure cancels on smooth surfaces. Edge-stopping uses the distance of the tap from
// the centre pixel's tangent plane (needs the centre normal only), measured in pixel footprints -> scale-free.
#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
#ifndef BLUR_RADIUS
#define BLUR_RADIUS 4
#endif
#define BLUR_TOLERANCE_PX 2.0     // tangent-plane distance tolerance, in full-res pixel footprints
uniform sampler2D u_tex;        // RG16F (visibility, z), half res
uniform sampler2D u_tex_normal; // full-res G-buffer normal
uniform vec4  u_proj_info;
uniform vec2  u_full_res;
uniform ivec2 u_dir;
uniform float u_px_scale;       // world size of a FULL-res pixel at depth 1 (persp) / absolute (ortho)
uniform float u_z_max;
out vec2 out_frag;

vec3 uv_to_view(vec2 uv, float z) {
#if AO_PERSPECTIVE
    return vec3((uv * u_proj_info.xy + u_proj_info.zw) * z, z);
#else
    return vec3((uv * u_proj_info.xy + u_proj_info.zw), z);
#endif
}
vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc * 4.0 - 2.0;
    float f = dot(fenc, fenc);
    float g = sqrt(1.0 - f / 4.0);
    return vec3(fenc * g, 1.0 - f / 2.0);
}
// full-res pixel centre represented by half-res texel p (matches the rotated-grid pick)
vec2 rep_px(ivec2 p) { return vec2(p * 2 + ivec2(p.y & 1, p.x & 1)) + 0.5; }

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 size = textureSize(u_tex, 0);
    vec2 c = texelFetch(u_tex, p, 0).rg;
    if (c.y >= u_z_max) { out_frag = c; return; }
    ivec2 fp = p * 2 + ivec2(p.y & 1, p.x & 1);
    vec3 N = decode_normal(texelFetch(u_tex_normal, fp, 0).xy) * vec3(1, 1, -1);
    vec3 P = uv_to_view(rep_px(p) / u_full_res, c.y);
#if AO_PERSPECTIVE
    float inv_tol = 1.0 / (BLUR_TOLERANCE_PX * u_px_scale * c.y);
#else
    float inv_tol = 1.0 / (BLUR_TOLERANCE_PX * u_px_scale);
#endif
    float sum = c.x, wsum = 1.0;
    for (int i = -BLUR_RADIUS; i <= BLUR_RADIUS; ++i) {
        if (i == 0) continue;
        ivec2 q = clamp(p + u_dir * i, ivec2(0), size - 1);
        vec2 s = texelFetch(u_tex, q, 0).rg;
        vec3 S = uv_to_view(rep_px(q) / u_full_res, s.y);
        float d = abs(dot(S - P, N)) * inv_tol;
        float w = (abs(i) == BLUR_RADIUS ? 0.5 : 1.0) * clamp(1.0 - d, 0.0, 1.0);
        sum += s.x * w;
        wsum += w;
    }
    out_frag = vec2(sum / wsum, c.y);
}
