#version 410 core

// Screen-space shadows of the key light, at half resolution: the visibility of the light, which the compose pass applies
// to its light (as the ambient occlusion is applied to the light of the surroundings). The output is the ambient
// occlusion's, (visibility, linear depth), so the same temporal filter, blur and upsample follow (shaders/ssao).
//
// The light is a disk of angular radius u_light_radius about u_light_dir. The view rays through the pixels on the line
// from a pixel toward the light, on screen, span a plane that holds the light's direction: walking that line finds the
// occluders of the light's chord in that plane, a slice along the light (Jimenez et al. 2016, "Practical Real-Time
// Strategies for Accurate Indirect Occlusion"). As in the Quality ambient occlusion (ssao_bitmask.frag), each depth
// sample occludes the arc between the directions to its front face and to the same point pushed back by an assumed
// thickness, and the union of the arcs is kept as a 32-sector mask (Therrien et al. 2023, "Screen Space Indirect Lighting
// with Visibility Bitmask"), here over the chord of the light: a thin occluder only shadows the part of the light it
// covers, and how much of the chord is covered makes the penumbra. The chord is taken at an offset from the centre of
// the disk across the slice, drawn per pixel and per frame with the density of the disk's chord lengths: the sectors
// are uniform along each chord, and the blur and the temporal filter average the chords into the whole disk.
//
// Scale-free like the ambient occlusion: the steps are spaced log-uniformly in screen space between u_r_min and u_r_max
// (shadows reach further than the occlusion), and the assumed thickness of an occluder is proportional to its lateral
// distance.
//
// Every step reads the half-res depth (mip 1), not a level of a footprint growing with the distance as the ambient
// occlusion does. A coarser level holds the depth of one pixel per texel (rotated grid): the far steps would then find
// an occluder's edge only to within the texel, the same texel for every frame and for all the pixels about it, and the
// error does not average out but shows as a fixed mottle in the soft shadows, the larger the closer the view. The steps
// of neighbouring pixels lie along the same line (the light's direction hardly changes over the screen), so they read
// close texels and stay in the texture cache, unlike the ambient occlusion's slices, which turn from pixel to pixel.
//
// What it cannot see: occluders outside the view, and what is hidden behind the nearest surface, whose thickness is a
// guess. The key light is on the viewer's side, so the occluders of what is seen are mostly in front of it, on screen;
// shadow maps would cover the rest.
//
// The number of steps, the thickness and the reach (ssao::SHADOW_R_MAX_FRACTION) are those closest to the ray traced
// visibility of the light over the atoms of proteins in space fill (tubulin, 4EB6, and the dengue envelope protein,
// 1K4R; a 10 degree light): a mean error of 0.05 and 0.07, unbiased. A thicker guess darkens what a thin part in front
// does not shadow (an atom or a loop standing out, with nothing behind it), a thinner one misses what is shadowed by the
// atoms hidden behind those in front. More steps, or a longer reach, find more of the former than of the latter.

#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
#define SH_NUM_STEPS    20      // along the line toward the light
#define SH_THICKNESS    1.0     // assumed thickness of an occluder, as a multiple of its lateral distance
#define SH_BIAS_PX      0.5     // height above the tangent plane, in pixel footprints, up to which a sample is taken
                                // for the receiving surface itself

uniform sampler2D u_tex_linear_depth;   // R32F linear view depth with mips (rotated grid): level 1 is read, texelFetch only
uniform sampler2D u_tex_normal;         // full-res G-buffer normal (RG16, spheremap encoded, view space)

uniform vec4  u_proj_info;
uniform vec2  u_full_res;
uniform float u_px_scale;
uniform float u_r_min;
uniform float u_r_max;
uniform float u_z_max;
uniform int   u_frame;
uniform vec3  u_light_dir;      // toward the light, normalized, in the view space of these passes (z into the screen)
uniform float u_light_radius;   // angular radius of the light (radians, > 0)

out vec2 out_frag;  // (visibility, linear depth)

const float PI = 3.14159265;

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

// acos, Abramowitz & Stegun 4.4.45: |error| <= 6.7e-5
float fast_acos(float x) {
    float a = abs(x);
    float r = sqrt(1.0 - a) * (1.5707288 + a * (-0.2121144 + a * (0.0742610 + a * -0.0187293)));
    return x >= 0.0 ? r : PI - r;
}

// The angle in the slice of the in-plane direction (a, b) = (along V, along o), from V toward o
float slice_angle(vec2 e) {
    float c = e.x * inversesqrt(max(dot(e, e), 1e-30));
    return (e.y < 0.0 ? -1.0 : 1.0) * fast_acos(clamp(c, -1.0, 1.0));
}

// Offsets of the chord from the centre of the disk, in units of its radius: the quantiles (k + 0.5) / 16 of the
// semicircle distribution, the density of the lengths of the chords
const float CHORD_OFFSET[16] = float[16](
    -0.85853, -0.70068, -0.57309, -0.45839, -0.35096, -0.24800, -0.14780, -0.04911,
     0.04911,  0.14780,  0.24800,  0.35096,  0.45839,  0.57309,  0.70068,  0.85853);

void main() {
    ivec2 hp = ivec2(gl_FragCoord.xy);
    ivec2 fp = hp * 2 + ivec2(hp.y & 1, hp.x & 1);   // full-res pixel that mip level 1 holds (depth_downsample.frag)

    float z = texelFetch(u_tex_linear_depth, hp, 1).r;
    if (z >= u_z_max) {
        out_frag = vec2(1.0, z);
        return;
    }

    vec2 frag_px = vec2(fp) + 0.5;
    vec3 P = uv_to_view(frag_px / u_full_res, z);
    vec3 N = decode_normal(texelFetch(u_tex_normal, fp, 0).xy) * vec3(1, 1, -1);
#if AO_PERSPECTIVE
    vec3  V = -normalize(P);
    float px_world = u_px_scale * z;
#else
    vec3  V = vec3(0, 0, -1);
    float px_world = u_px_scale;
#endif

    // The same 4x4 interleaved pattern, permuted per frame, as the ambient occlusion: k picks the chord, and its bit
    // reversal the phase of the steps (a Hammersley set over the 16 pixels and frames of a period)
    const int BAYER[16] = int[16](0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5);
    int   s_ofs = int(bitfieldReverse(uint(u_frame & 15)) >> 28u);
    int   k     = (BAYER[(hp.y & 3) * 4 + (hp.x & 3)] + s_ofs) & 15;
    float rk    = float(bitfieldReverse(uint(k)) >> 28u) / 16.0;

    // The chord of this pixel and frame: the light's direction turned across the plane of V and the light, and the
    // half width of the disk there
    vec3  L = u_light_dir;
    vec3  b = cross(V, L);
    float b_len = length(b);
    if (b_len < 1e-3) {
        // The light straight behind the viewer: what it shadows is hidden
        out_frag = vec2(1.0, z);
        return;
    }
    b /= b_len;
    float y  = CHORD_OFFSET[k] * u_light_radius;
    vec3  Lk = L * cos(y) + b * sin(y);
    float half_chord = sqrt(max(u_light_radius * u_light_radius - y * y, 1e-8));

    // The slice: the plane of V and Lk, o the direction in it perpendicular to V on the side of the light. On screen,
    // the line from the pixel along the projection of Lk.
    vec3  o = normalize(Lk - V * dot(Lk, V));
#if AO_PERSPECTIVE
    vec2  d = (Lk.xy * z - P.xy * Lk.z) / u_proj_info.xy * u_full_res;
#else
    vec2  d = Lk.xy / u_proj_info.xy * u_full_res;
#endif
    float d_len = length(d);
    if (d_len < 1e-6) {
        out_frag = vec2(1.0, z);
        return;
    }
    d /= d_len;

    // The chord, as angles in the slice from V
    float arc_lo  = fast_acos(clamp(dot(Lk, V), -1.0, 1.0)) - half_chord;
    float inv_arc = 1.0 / (2.0 * half_chord);
    float bias    = SH_BIAS_PX * px_world;

    // Steps at s_i = r_min * ratio^(i + phase)
    float log_ratio  = log2(u_r_max / u_r_min);
    float step_ratio = exp2(log_ratio / float(SH_NUM_STEPS));
    float s          = u_r_min * exp2(rk * log_ratio / float(SH_NUM_STEPS));

    // Stochastic rounding of the arcs onto the sectors keeps small arcs unbiased on average
    float bit_ofs = fract(rk * 7.0 + 0.5);
    uint  bits = 0u;

    for (int i = 0; i < SH_NUM_STEPS; ++i, s *= step_ratio) {
        vec2 spx = frag_px + d * s;                     // s: lateral distance in full-res px
        if (any(lessThan(spx, vec2(0.0))) || any(greaterThanEqual(spx, u_full_res))) break;

        ivec2 th = ivec2(spx) >> 1;
        float sz = texelFetch(u_tex_linear_depth, th, 1).r;
        if (sz >= u_z_max) continue;
        ivec2 tx = th * 2 + ivec2(th.y & 1, th.x & 1);  // the full-res pixel it holds
        vec3  S  = uv_to_view((vec2(tx) + 0.5) / u_full_res, sz);

        // Front face and the same point pushed back along its view ray by the assumed thickness
        float thick = SH_THICKNESS * s * px_world;
#if AO_PERSPECTIVE
        vec3 push = normalize(S) * thick;
#else
        vec3 push = vec3(0, 0, thick);
#endif
        vec3  Df = S - P;
        float hf = dot(Df, N);
        float hb = hf + dot(push, N);
        if (hf <= bias && hb <= bias) continue;         // under the tangent plane: the receiving surface, or behind it
        vec3  Db = Df + push;

        float uf = (slice_angle(vec2(dot(Df, V), dot(Df, o))) - arc_lo) * inv_arc;
        float ub = (slice_angle(vec2(dot(Db, V), dot(Db, o))) - arc_lo) * inv_arc;
        float lo = min(uf, ub);
        float hi = max(uf, ub);
        if (hi <= 0.0 || lo >= 1.0) continue;
        uint b0 = uint(clamp(floor(lo * 32.0 + bit_ofs), 0.0, 32.0));
        uint b1 = uint(clamp(floor(hi * 32.0 + bit_ofs), 0.0, 32.0));
        if (b1 > b0) {
            uint count = b1 - b0;
            bits |= (count >= 32u ? 0xFFFFFFFFu : ((1u << count) - 1u)) << b0;
            if (bits == 0xFFFFFFFFu) break;
        }
    }

    out_frag = vec2(1.0 - float(bitCount(bits)) / 32.0, z);
}
