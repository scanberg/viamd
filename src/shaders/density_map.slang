// density_map.slang -- density_map_composite's (density_map.h): the density layers of one plot into its canvas, a
// thread per texel, the first layer at the bottom. A layer is contours (fills of the regions above levels, lines along
// them of a width in pixels) or a heat map, of the weighted sum of the channels of its grid, sampled with a cubic
// B-spline (16 texels: wrapped where the domain is periodic, nothing outside it where it is bounded) whose gradient,
// in the canvas' pixels, makes the edges of the fills and the lines antialiased and of a constant width. The canvas may
// be clipped to a triangle, with antialiased edges. Out premultiplied. Mirrored in density_map.cpp.

#include "../../ext/mdlib/src/shaders/md_gpu.slang"

static const uint MAX_LEVELS = 8;       // DENSITY_MAP_MAX_LEVELS
static const uint GRID       = 512;     // DENSITY_MAP_DIM

// A draw's parameters: PARAMS float4 from params + PARAMS * its index (density_map.cpp's GpuDraw)
static const uint PARAMS   = 32;
static const uint P_WEIGHTS = 0;    // of the channels
static const uint P_LEVELS  = 1;    // 2: densities, the highest first
static const uint P_FILL    = 3;    // MAX_LEVELS, straight
static const uint P_LINE    = 11;   // the lines' colour
static const uint P_CMAP    = 12;   // 16
static const uint P_U0      = 28;   // mode (0 contours, 1 heat map), fills, lines, levels: their bits
static const uint P_U1      = 29;   // the colour map's colours, its log scale
static const uint P_F       = 30;   // line width (pixels), the colour map's maximum, opacity

struct Grid {
    DescriptorHandle<Texture2D<float4>> tex;
};

struct Args {
    float4  uv_rect;            // the uv of the grid at the canvas' lower left corner, and the extent
    float4  clip01;             // the triangle (uv): its first two corners
    float4  clip2;              // its third; z: 1 clipped to it
    uint    width, height, count, periodic;
    float4* params;
    Grid*   grids;              // a draw's
    DescriptorHandle<RWTexture2D<float4>> dst;
    uint    _p0, _p1;
};
MD_SHADER_ARGS(Args);

struct Draw {
    float4* p;
    Grid    grid;
    float4 weights;
};

float4 texel(Draw d, bool periodic, int2 p) {
    if (periodic) {
        p &= int(GRID - 1);
    } else if (any(p < 0) || any(p >= int(GRID))) {
        return float4(0.0);
    }
    return d.grid.tex.Load(int3(p, 0));
}

// The density at uv, and its derivatives along u and v
float3 sample_density(Draw d, bool periodic, float2 uv) {
    float2 st = uv * float(GRID) - 0.5;
    float2 i  = floor(st);
    float2 f  = st - i;
    float2 f2 = f * f, f3 = f2 * f;
    float2 w[4], dw[4];
    w[0] = (1.0 / 6.0) * (-f3 + 3.0 * f2 - 3.0 * f + 1.0);
    w[1] = (1.0 / 6.0) * (3.0 * f3 - 6.0 * f2 + 4.0);
    w[2] = (1.0 / 6.0) * (-3.0 * f3 + 3.0 * f2 + 3.0 * f + 1.0);
    w[3] = (1.0 / 6.0) * f3;
    dw[0] = -0.5 * (1.0 - f) * (1.0 - f);
    dw[1] = 1.5 * f2 - 2.0 * f;
    dw[2] = -1.5 * f2 + f + 0.5;
    dw[3] = 0.5 * f2;
    int2 base = int2(i) - 1;
    float v = 0.0, ds = 0.0, dt = 0.0;
    [ForceUnroll] for (int ky = 0; ky < 4; ++ky) {
        [ForceUnroll] for (int kx = 0; kx < 4; ++kx) {
            float t = dot(texel(d, periodic, base + int2(kx, ky)), d.weights);
            v  += w[kx].x  * w[ky].y  * t;
            ds += dw[kx].x * w[ky].y  * t;
            dt += w[kx].x  * dw[ky].y * t;
        }
    }
    return float3(v, ds * float(GRID), dt * float(GRID));
}

float4 over(float4 src, float4 dst) { return src + dst * (1.0 - src.a); }
float4 premultiply(float4 c) { return float4(c.rgb * c.a, c.a); }

float level_of(Draw d, uint i) { i &= MAX_LEVELS - 1u; return d.p[P_LEVELS + (i >> 2)][i & 3]; }

// How much of the pixel is within the triangle: its distance in pixels to the nearest edge, inside. duv: the uv a pixel
// spans along x and y.
float clip_coverage(Args a, float2 uv, float2 duv) {
    float2 t[3] = { a.clip01.xy, a.clip01.zw, a.clip2.xy };
    float s = sign((t[1].x - t[0].x) * (t[2].y - t[0].y) - (t[1].y - t[0].y) * (t[2].x - t[0].x));
    float cov = 1.0;
    [ForceUnroll] for (int k = 0; k < 3; ++k) {
        float2 p = t[k], q = t[(k + 1) % 3];
        float e = s * ((q.x - p.x) * (uv.y - p.y) - (q.y - p.y) * (uv.x - p.x));
        float2 de = s * float2(-(q.y - p.y), q.x - p.x) * duv;
        cov = min(cov, saturate(e / max(length(de), 1.0e-30) + 0.5));
    }
    return cov;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void composite(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    // The threads past the canvas run along (clamped) and store nothing: lavapipe, which makes the loads of the
    // levels' loop uniform, takes their address from any thread of the group, one that has returned too
    const bool inside = tid.x < a.width && tid.y < a.height;
    tid.xy = min(tid.xy, uint2(a.width, a.height) - 1u);
    // As GL had it: the fragment's coordinates, rows from the bottom
    float2 inv_res = 1.0 / float2(a.width, a.height);
    float2 fc  = float2(float(tid.x) + 0.5, float(a.height - tid.y) - 0.5);
    float2 uv  = a.uv_rect.xy + fc * inv_res * a.uv_rect.zw;
    float2 duv = inv_res * a.uv_rect.zw;
    float clip = a.clip2.z != 0.0 ? clip_coverage(a, uv, duv) : 1.0;
    bool periodic = a.periodic != 0u;

    float4 o = float4(0.0);
    for (uint n = 0u; n < a.count; ++n) {
        Draw d;
        d.p = a.params + n * PARAMS;
        d.grid = a.grids[n];
        d.weights = d.p[P_WEIGHTS];
        uint4  u0 = asuint(d.p[P_U0]);
        uint2  u1 = asuint(d.p[P_U1].xy);
        float4 f  = d.p[P_F];
        uint   num_levels = min(u0.w, MAX_LEVELS);
        float3 s = sample_density(d, periodic, uv);
        float dv = s.x;
        float g = max(length(s.yz * duv), 1.0e-30);
        float4 c = float4(0.0);
        if (u0.x == 0u) {
            if (u0.y != 0u) {
                for (uint k = 0u; k < num_levels; ++k) {
                    uint i = num_levels - 1u - k;
                    float cov = saturate((dv - level_of(d, i)) / g + 0.5);
                    c = over(premultiply(d.p[P_FILL + i]) * cov, c);
                }
            }
            if (u0.z != 0u) {
                float al = 0.0;
                for (uint i = 0u; i < num_levels; ++i) {
                    float dist = abs(dv - level_of(d, i)) / g;
                    al = max(al, saturate(0.5 * f.x + 0.5 - dist));
                }
                c = over(premultiply(d.p[P_LINE]) * al, c);
            }
        } else {
            float v = dv / max(f.y, 1.0e-30);
            if (u1.y != 0u) v = 1.0 + log10(max(v, 1.0e-4)) * 0.25;
            v = saturate(v);
            uint  count = clamp(u1.x, 1u, 16u);
            float sc = v * float(count - 1u);
            uint  i0 = uint(floor(sc));
            uint  i1 = min(i0 + 1u, count - 1u);
            float4 col = lerp(d.p[P_CMAP + i0], d.p[P_CMAP + i1], sc - float(i0));
            col.a *= smoothstep(0.0, 0.06, v);
            c = premultiply(col);
        }
        o = over(c * f.z * clip, o);
    }
    if (inside) a.dst[int2(tid.xy)] = o;
}
