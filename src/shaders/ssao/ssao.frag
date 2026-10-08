#version 410 core

// Scale-free screen-space ambient obscurance, evaluated at half resolution.
//
// There is no world-space radius. Samples are spread log-uniformly in screen space between u_r_min and u_r_max
// (full-res pixels; u_r_max is a fixed fraction of the viewport height), and every sample is judged at its *own* scale:
// its falloff radius is AO_FALLOFF_SCALE x its lateral distance. The occlusion test is therefore purely angular, which
// makes the result invariant to a uniform scaling of scene + camera: the same constants work for a C60 at 15 Å and a
// multi-million atom assembly at 3000 Å. Occluders that are much closer to the camera than their lateral offset (the
// usual SSAO halo around silhouettes) fall outside every sample's falloff and do not contribute.
//
// Cost is fixed per pixel: AO_NUM_SAMPLES point fetches from a rotated-grid depth mip chain, with the mip chosen so the
// fetch footprint grows with the sample distance (McGuire et al. 2012, "Scalable Ambient Obscurance"). This keeps the
// texture cache happy no matter how large the radius is on screen.

#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
#ifndef AO_NUM_SAMPLES
#define AO_NUM_SAMPLES 16
#endif
#define AO_LOG_Q            3       // a sample at distance s px reads mip floor(log2(s)) - AO_LOG_Q (but at least 1)
#define AO_MAX_MIP          5
#define AO_FALLOFF_SCALE    2.5     // occluders steeper than ~66 deg above a sample's own scale fade out
#define AO_BIAS             0.05    // n.v bias against self occlusion from depth quantisation
#define AO_SCALE_POWER      0.0     // sample weight ~ (s / r_max)^p: 0 weights every octave equally

uniform sampler2D u_tex_linear_depth;   // R32F linear view depth, mips 0..AO_MAX_MIP (rotated grid), texelFetch only
uniform sampler2D u_tex_normal;         // full-res G-buffer normal (RG16, spheremap encoded, view space)

uniform vec4  u_proj_info;
uniform vec2  u_full_res;       // full-res size in pixels
uniform float u_px_scale;       // world size of one full-res pixel at view depth 1 (persp) or absolute (ortho)
uniform float u_r_min;          // full-res px
uniform float u_r_max;          // full-res px
uniform float u_intensity;
uniform float u_z_max;
uniform int   u_frame;          // 0 for a stable pattern, frame index when TAA can integrate it

out vec2 out_frag;  // (visibility, linear depth): depth is carried along for the bilateral blur and upsample

const float GOLDEN_ANGLE = 2.39996323;

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

void main() {
    ivec2 hp = ivec2(gl_FragCoord.xy);
    // Full-res pixel that mip level 1 holds for this half-res texel (see depth_downsample.frag)
    ivec2 fp = hp * 2 + ivec2(hp.y & 1, hp.x & 1);

    float z = texelFetch(u_tex_linear_depth, hp, 1).r;
    if (z >= u_z_max) {
        out_frag = vec2(1.0, z);
        return;
    }

    vec2 frag_px = vec2(fp) + 0.5;
    vec3 P = uv_to_view(frag_px / u_full_res, z);
    vec3 N = decode_normal(texelFetch(u_tex_normal, fp, 0).xy) * vec3(1, 1, -1);
#if AO_PERSPECTIVE
    float px_world = u_px_scale * z;
#else
    float px_world = u_px_scale;
#endif

    // 4x4 interleaved pattern: 16 (rotation, radial phase) pairs from a Hammersley set, Bayer-ordered so that
    // neighbours differ. The blur integrates exactly one period of it. u_frame decorrelates frames for TAA.
    const int BAYER[16] = int[16](0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5);
    int   k   = BAYER[(hp.y & 3) * 4 + (hp.x & 3)];
    float rk  = float(bitfieldReverse(uint(k)) >> 28u) / 16.0;
    float rot = 6.2831853 * fract((float(k) + 0.5) / 16.0 + float(u_frame) * 0.618034);
    float jr  = fract(rk + 1.0 / 32.0 + float(u_frame) * 0.7548777);
    vec2  cs  = vec2(cos(rot), sin(rot));

    float log_ratio = log2(u_r_max / u_r_min);
    float ao = 0.0;
    float w_sum = 0.0;

    for (int i = 0; i < AO_NUM_SAMPLES; ++i) {
        float t = (float(i) + fract(jr + float(i) * 0.618034)) / float(AO_NUM_SAMPLES);
        float s = u_r_min * exp2(t * log_ratio);                  // lateral distance in full-res px
        float a = float(i) * GOLDEN_ANGLE;
        vec2  d = vec2(cos(a), sin(a));                            // constant after unrolling
        d = vec2(d.x * cs.x - d.y * cs.y, d.x * cs.y + d.y * cs.x);
        vec2  spx = frag_px + d * s;
        float ws = pow(s / u_r_max, AO_SCALE_POWER);
        w_sum += ws;
        if (any(lessThan(spx, vec2(0.0))) || any(greaterThanEqual(spx, u_full_res))) continue;  // off screen: unoccluded

        int   m  = clamp(int(log2(s)) - AO_LOG_Q, 1, AO_MAX_MIP);
        ivec2 tx = ivec2(spx) >> m;
        float sz = texelFetch(u_tex_linear_depth, tx, m).r;
        // Walk the rotated-grid picks back to the full-res pixel this depth came from, so the reconstructed point
        // lies exactly on the surface (reconstructing at spx instead gives false occlusion on curved atoms).
        for (int l = m; l > 0; --l) tx = tx * 2 + ivec2(tx.y & 1, tx.x & 1);
        vec3  v  = uv_to_view((vec2(tx) + 0.5) / u_full_res, sz) - P;
        float vv = dot(v, v);
        float r  = AO_FALLOFF_SCALE * s * px_world;               // falloff radius at this sample's own scale
        float fall = clamp(1.0 - vv / (r * r), 0.0, 1.0);
        float vn = dot(v, N) * inversesqrt(vv + 1e-12);
        ao += max(vn - AO_BIAS, 0.0) * fall * ws;
    }

    ao *= 1.0 / (w_sum * (1.0 - AO_BIAS));
    out_frag = vec2(clamp(1.0 - u_intensity * ao, 0.0, 1.0), z);
}
