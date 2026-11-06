#version 410 core

// Scale-free horizon-based ambient occlusion with a visibility bitmask, evaluated at half resolution ("Quality" mode;
// ssao.frag is the cheaper "Performance" mode). Same inputs and output as ssao.frag, so the temporal filter, blur and
// upsample are shared.
//
// Each pixel walks a few screen-space slices through itself (Jimenez et al. 2016, "Practical Real-Time Strategies for
// Accurate Indirect Occlusion"). Every depth sample along a slice occludes the arc between the directions to its front
// face and to the same point pushed back by an assumed thickness. The union of those arcs is kept as a 32-sector bitmask
// per slice (Therrien et al. 2023, "Screen Space Indirect Lighting with Visibility Bitmask"). A point sample either hits
// an occluder or misses it; a slice finds the occluder from every step that lands on it, so a crevasse is found by
// every slice that crosses it. Unlike a max-horizon (GTAO), a thin occluder (a small atom in front of a deep cavity)
// only hides the arc it covers, not everything below it.
//
// The sectors are spaced so that each carries the same share of the cosine-weighted hemisphere (the slice measure
// cos(theta - gamma) |sin theta| integrated in closed form), so visibility is simply 1 - occupied / 32.
//
// No atan, sin or cos per step: a step whose front and back points both lie behind the projected normal's tangent
// is skipped on a dot product (most steps around convex atoms), cos(2 theta - gamma) is expanded in the components of
// the normalized step direction, and the angle that remains in the measure comes from a polynomial acos (6.7e-5 rad,
// far below the 0.1 rad of a sector).
//
// Scale-free like ssao.frag: steps are spaced log-uniformly in screen space between u_r_min and u_r_max, the assumed
// thickness of a sample is proportional to its lateral distance, and the depth comes from the rotated-grid mip chain
// with the footprint growing with the step distance (McGuire et al. 2012).

#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
#ifndef AO_NUM_SLICES
#define AO_NUM_SLICES   2
#endif
#ifndef AO_NUM_STEPS
#define AO_NUM_STEPS    8       // per side of a slice
#endif
#define AO_LOG_Q            3   // a step at distance s px reads mip floor(log2(s)) - AO_LOG_Q (but at least 1)
#define AO_MAX_MIP          5
#define AO_THICKNESS        2.0 // assumed thickness of an occluder, as a multiple of its lateral distance
#define AO_POWER_SCALE      0.2 // output = visibility ^ (u_intensity * AO_POWER_SCALE): the default intensity (5) gives
                                // the plain cosine-weighted visibility, which sits close to ray-traced AO

uniform sampler2D u_tex_linear_depth;   // R32F linear view depth, mips 0..AO_MAX_MIP (rotated grid), texelFetch only
uniform sampler2D u_tex_normal;         // full-res G-buffer normal (RG16, spheremap encoded, view space)

uniform vec4  u_proj_info;
uniform vec2  u_full_res;
uniform float u_px_scale;
uniform float u_r_min;
uniform float u_r_max;
uniform float u_intensity;
uniform float u_z_max;
uniform int   u_frame;

out vec2 out_frag;  // (visibility, linear depth)

const float PI      = 3.14159265;
const float HALF_PI = 1.57079633;

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

// Cosine-weighted measure of the slice from the view direction (theta = 0) to the in-plane direction (a, b) =
// (cos theta, sin theta), for a projected normal at angle gamma: integral of cos(t - gamma) |sin t| dt
// = sign(theta) (cos gamma + 2 theta sin gamma - cos(2 theta - gamma)) / 4. Monotonic over the hemisphere.
float slice_measure(float a, float b, float cos_g, float sin_g) {
    float theta = (b < 0.0 ? -1.0 : 1.0) * fast_acos(a);
    float cos_2t_g = (a * a - b * b) * cos_g + 2.0 * a * b * sin_g;
    return (b < 0.0 ? -0.25 : 0.25) * (cos_g + 2.0 * theta * sin_g - cos_2t_g);
}

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

    // Same 4x4 interleaved (rotation, radial phase) pattern and per-frame permutation as ssao.frag
    const int BAYER[16] = int[16](0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5);
    int   s_ofs = int(bitfieldReverse(uint(u_frame & 15)) >> 28u);
    int   k     = (BAYER[(hp.y & 3) * 4 + (hp.x & 3)] + s_ofs) & 15;
    float rk    = float(bitfieldReverse(uint(k)) >> 28u) / 16.0;
    float phi0  = PI * (float(k) + 0.5) / (16.0 * float(AO_NUM_SLICES));
    vec2  d     = vec2(cos(phi0), sin(phi0));
    const float SLICE_STEP = PI / float(AO_NUM_SLICES);
    const vec2  slice_rot  = vec2(cos(SLICE_STEP), sin(SLICE_STEP));

    // Steps at s_i = r_min * ratio^(i + phase): one exp2 per side, then a multiply per step
    float log_ratio  = log2(u_r_max / u_r_min);
    float step_ratio = exp2(log_ratio / float(AO_NUM_STEPS));
    float s_first[2] = float[2](u_r_min * exp2(rk * log_ratio / float(AO_NUM_STEPS)),
                                u_r_min * exp2(fract(rk + 0.5) * log_ratio / float(AO_NUM_STEPS)));

    float vis_sum = 0.0;
    float w_sum   = 0.0;

    for (int sl = 0; sl < AO_NUM_SLICES; ++sl) {
        if (sl > 0) d = vec2(d.x * slice_rot.x - d.y * slice_rot.y, d.x * slice_rot.y + d.y * slice_rot.x);

        // Slice frame: o is the in-plane direction perpendicular to V on the +d side
        vec3  dir3 = vec3(d, 0.0);
        vec3  o    = normalize(dir3 - V * dot(dir3, V));
        vec3  axis = cross(o, V);
        vec3  n_p  = N - axis * dot(N, axis);
        float n_len = length(n_p);
        if (n_len < 1e-4) continue;
        float cos_g = clamp(dot(n_p, V) / n_len, -1.0, 1.0);
        float sgn_g = dot(n_p, o) < 0.0 ? -1.0 : 1.0;
        float gamma = sgn_g * fast_acos(cos_g);
        float sin_g = sgn_g * sqrt(1.0 - cos_g * cos_g);
        float total = cos_g + gamma * sin_g;                    // measure(hi) - measure(lo)
        float m_lo  = -0.5 * (cos_g + (gamma - HALF_PI) * sin_g); // measure(lo), lo = gamma - pi/2
        float inv_total = 1.0 / total;

        // Stochastic rounding of the arcs onto the sectors keeps small arcs unbiased on average
        float bit_ofs = fract(rk * 7.0 + float(sl) * 0.618034);
        uint  bits = 0u;

        for (int side = 0; side < 2; ++side) {
            vec2  sd = side == 0 ? d : -d;
            float s  = s_first[side];
            for (int i = 0; i < AO_NUM_STEPS; ++i, s *= step_ratio) {
                vec2  spx = frag_px + sd * s;                    // s: lateral distance in full-res px
                if (any(lessThan(spx, vec2(0.0))) || any(greaterThanEqual(spx, u_full_res))) break;

                int   m  = clamp(findMSB(int(s)) - AO_LOG_Q, 1, AO_MAX_MIP);
                ivec2 tx = ivec2(spx) >> m;
                float sz = texelFetch(u_tex_linear_depth, tx, m).r;
                if (sz >= u_z_max) continue;
                for (int l = m; l > 0; --l) tx = tx * 2 + ivec2(tx.y & 1, tx.x & 1);
                vec3 S = uv_to_view((vec2(tx) + 0.5) / u_full_res, sz);

                // Front face and the same point pushed back along its view ray by the assumed thickness
                float thick = AO_THICKNESS * s * px_world;
#if AO_PERSPECTIVE
                vec3 push = normalize(S) * thick;
#else
                vec3 push = vec3(0, 0, thick);
#endif
                vec3  Df = S - P;
                float hf = dot(Df, n_p);
                float hb = hf + dot(push, n_p);
                if (hf <= 0.0 && hb <= 0.0) continue;   // both behind the tangent: no part of the arc is in the hemisphere
                vec3  Db = Df + push;

                // In-plane directions of the two points; outside the hemisphere they clamp to its edge on their side
                vec2  ef = vec2(dot(Df, V), dot(Df, o));
                vec2  eb = vec2(dot(Db, V), dot(Db, o));
                ef *= inversesqrt(max(dot(ef, ef), 1e-30));
                eb *= inversesqrt(max(dot(eb, eb), 1e-30));
                float uf = hf > 0.0 ? (slice_measure(ef.x, ef.y, cos_g, sin_g) - m_lo) * inv_total : (ef.y > 0.0 ? 1.0 : 0.0);
                float ub = hb > 0.0 ? (slice_measure(eb.x, eb.y, cos_g, sin_g) - m_lo) * inv_total : (eb.y > 0.0 ? 1.0 : 0.0);
                uint  b0 = uint(clamp(floor(min(uf, ub) * 32.0 + bit_ofs), 0.0, 32.0));
                uint  b1 = uint(clamp(floor(max(uf, ub) * 32.0 + bit_ofs), 0.0, 32.0));
                if (b1 > b0) {
                    uint count = b1 - b0;
                    bits |= (count >= 32u ? 0xFFFFFFFFu : ((1u << count) - 1u)) << b0;
                }
            }
        }

        // Slice weight: projected normal length x its share of the cosine-weighted hemisphere
        float w = n_len * total;
        vis_sum += w * (1.0 - float(bitCount(bits)) / 32.0);
        w_sum   += w;
    }

    float vis = w_sum > 0.0 ? vis_sum / w_sum : 1.0;
    out_frag = vec2(pow(clamp(vis, 0.0, 1.0), u_intensity * AO_POWER_SCALE), z);
}
