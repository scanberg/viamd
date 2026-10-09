#version 410 core

// Isosurface raycaster, exact mode: exact intersection of the ray with the trilinear field. The fast mode,
// one filtered sample per voxel, is isosurface_fast.frag; both share the uniforms below.
//
// Inside a cell of the grid (the cube between eight voxel centres) the trilinearly interpolated field
// along a ray is a cubic in the ray parameter. The ray is walked cell by cell (3D DDA), and in a cell whose
// corner values bracket an isovalue the cubic is split where its derivative vanishes into pieces that are
// monotone, and each piece whose ends lie on different sides of the isovalue holds exactly one crossing,
// found by safeguarded Newton iteration on the cubic. No crossing of the trilinear surface is missed,
// however thin it is along the ray or however grazing the ray, and where a crossing is found does not
// depend on any sampling: no step length, no jitter, no temporal accumulation.
//
// Empty space: the block min/max grid (block_minmax.frag) is walked first. A block whose range holds no
// isovalue holds no surface and is crossed in one step; inside it the ray is on a known side of every
// surface, so its absorption is added analytically. Cells are only visited inside blocks that can hold a
// surface, and there the cell min/max grid (cell_minmax.frag) gives the range of each cell in one fetch:
// only the cells whose range holds an isovalue load their corners and get the cubic. When u_use_proxy is
// set, the block proxy (block_proxy.vert) has also narrowed every ray to the span that meets a block that
// can hold a surface, and discarded the pixels whose ray meets none.
//
// Membership: a point is INSIDE isosurface i when the density is on the far side of its value from zero,
// d >= v for v >= 0 and d <= v for v < 0. Crossing a surface toggles membership of exactly that surface,
// and the optical densities of the surfaces enclosing the ray attenuate it in between.
//
// Every crossing is shaded with the model and light of the deferred compose pass, and the result is
// written as premultiplied linear radiance + coverage, to be blended over the HDR scene before tone
// mapping.
//
// Coordinates: model space is the texture space [0,1]^3 of the volume. Voxel space puts voxel centres on
// the integers, x = p * dim - 0.5; cell c spans [c, c+1] and its corners are voxels c .. c+1, clamped to
// the volume like CLAMP_TO_EDGE filtering does (the half voxel border cells are constant across).

#ifndef MAX_ISO
#define MAX_ISO 8
#endif
#if MAX_ISO > 8
#error "the range test holds at most 8 isovalues (g_iso_a, g_iso_b)"
#endif

#define MAX_CELL_HITS 8

layout(std140) uniform IsoUniforms {
    mat4  u_clip_to_model;      // inverse(proj * view * model)
    mat4  u_model_to_view;
    mat4  u_grad_offsets;       // columns 0-2: model space offsets for one voxel along view x, y, z

    vec3  u_clip_min;
    float u_block_size;         // voxels per block side of u_tex_minmax
    vec3  u_clip_max;
    float u_use_depth;

    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;
    vec3  u_light_dir;          // view space, normalized
    float u_field_beg;

    vec2  u_inv_res;
    float u_field_inv_ext;
    float u_optical_scale;      // optical density -> extinction per world unit

    float u_use_proxy;          // u_tex_entry / u_tex_exit hold the block proxy depths
    float u_entry_from_near;    // the near plane cuts the volume: rays start there, not at the proxy entry
    float u_use_blocks;         // u_tex_minmax holds the block min/max grid
    float u_use_cells;          // u_tex_cells holds the cell min/max grid
};

uniform float u_iso_values[MAX_ISO];
uniform vec4  u_iso_colors[MAX_ISO];
uniform float u_iso_tau[MAX_ISO];
uniform int   u_iso_count;

uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_depth;
uniform sampler2D u_tex_entry;  // nearest depth of the block proxy, 1 where there is none
uniform sampler2D u_tex_exit;   // farthest depth of the block proxy, 0 where there is none
uniform sampler3D u_tex_minmax; // RG32F, min and max over every block and its one voxel apron
uniform sampler3D u_tex_cells;  // RG16F, the range of the corners of every cell (rounded outwards), cell c at texel c + 1

#if defined(USE_COLOR_VOLUME)
uniform sampler3D u_tex_color_volume;
#elif defined(USE_FIELD)
// A scalar field on the density volume's own grid and the colour map it is shown through
uniform sampler3D u_tex_field;
uniform sampler2D u_tex_field_colormap;
#endif

layout(location = 0) out vec4 out_color;

const float PI      = 3.1415926535;
const float INV_PI  = 1.0 / PI;
const float T_MIN   = 0.005;    // early ray termination: below this transmittance the ray is opaque
const float HUGE    = 1.0e30;
const int   MAX_STEPS = 4096;   // per loop, a guard only

// -----------------------------------------------------------------------------
// The ray, shared by the functions below
// -----------------------------------------------------------------------------

ivec3 g_dim_m1;     // volume dimensions - 1
vec3  g_p0;         // ray start, model space; p(t) = g_p0 + t * g_ray, t in [0, 1]
vec3  g_ray;
vec3  g_o;          // the same ray in voxel space: x(t) = g_o + t * g_d
vec3  g_d;
float g_dlen;       // |g_d|: voxels per unit t
vec3  g_dir;        // g_d / g_dlen, the direction the cubics are parameterized along (in voxels)
float g_ext_per_t;  // extinction per unit tau and unit t
vec3  g_V;          // view space direction towards the eye

// Shading state of the ray
vec3  g_L;          // premultiplied radiance
float g_T;          // transmittance
float g_t_abs;      // absorption has been applied up to here
uint  g_inside;     // membership per isosurface
float g_tau;        // optical density of the surfaces the ray is inside of
bool  g_known;      // membership determined yet
vec4  g_iso_a;      // isovalues 0-3 and 4-7, unused slots repeating the first, for a branch free range test
vec4  g_iso_b;

// -----------------------------------------------------------------------------
// Field
// -----------------------------------------------------------------------------

float sample_volume(vec3 p) {
    return texture(u_tex_volume, p).r;
}

float voxel(ivec3 c) {
    return texelFetch(u_tex_volume, clamp(c, ivec3(0), g_dim_m1), 0).r;
}

// Corners of cell c: lo = (f000, f100, f010, f110), hi = the same at z + 1
void fetch_corners(ivec3 c, out vec4 lo, out vec4 hi) {
    lo = vec4(voxel(c), voxel(c + ivec3(1,0,0)), voxel(c + ivec3(0,1,0)), voxel(c + ivec3(1,1,0)));
    hi = vec4(voxel(c + ivec3(0,0,1)), voxel(c + ivec3(1,0,1)), voxel(c + ivec3(0,1,1)), voxel(c + ivec3(1,1,1)));
}

// After stepping into cell c along one axis: the four corners on the shared face move over, the other
// four are fetched
void step_corners_x(int s, ivec3 c, inout vec4 lo, inout vec4 hi) {
    if (s > 0) {
        lo.xz = lo.yw; hi.xz = hi.yw;
        lo.y = voxel(c + ivec3(1,0,0)); lo.w = voxel(c + ivec3(1,1,0));
        hi.y = voxel(c + ivec3(1,0,1)); hi.w = voxel(c + ivec3(1,1,1));
    } else {
        lo.yw = lo.xz; hi.yw = hi.xz;
        lo.x = voxel(c);                lo.z = voxel(c + ivec3(0,1,0));
        hi.x = voxel(c + ivec3(0,0,1)); hi.z = voxel(c + ivec3(0,1,1));
    }
}

void step_corners_y(int s, ivec3 c, inout vec4 lo, inout vec4 hi) {
    if (s > 0) {
        lo.xy = lo.zw; hi.xy = hi.zw;
        lo.z = voxel(c + ivec3(0,1,0)); lo.w = voxel(c + ivec3(1,1,0));
        hi.z = voxel(c + ivec3(0,1,1)); hi.w = voxel(c + ivec3(1,1,1));
    } else {
        lo.zw = lo.xy; hi.zw = hi.xy;
        lo.x = voxel(c);                lo.y = voxel(c + ivec3(1,0,0));
        hi.x = voxel(c + ivec3(0,0,1)); hi.y = voxel(c + ivec3(1,0,1));
    }
}

void step_corners_z(int s, ivec3 c, inout vec4 lo, inout vec4 hi) {
    if (s > 0) {
        lo = hi;
        hi = vec4(voxel(c + ivec3(0,0,1)), voxel(c + ivec3(1,0,1)), voxel(c + ivec3(0,1,1)), voxel(c + ivec3(1,1,1)));
    } else {
        hi = lo;
        lo = vec4(voxel(c), voxel(c + ivec3(1,0,0)), voxel(c + ivec3(0,1,0)), voxel(c + ivec3(1,1,0)));
    }
}

// The trilinear interpolant of a cell along a + s * b (local coordinates, s in voxels): c0 + c1 s + c2 s^2 + c3 s^3
vec4 cubic_coeffs(vec4 lo, vec4 hi, vec3 a, vec3 b) {
    float k0 = lo.x;
    float k1 = lo.y - lo.x;
    float k2 = lo.z - lo.x;
    float k3 = hi.x - lo.x;
    float k4 = lo.w - lo.z - lo.y + lo.x;                               // xy
    float k5 = hi.z - hi.x - lo.z + lo.x;                               // yz
    float k6 = hi.y - hi.x - lo.y + lo.x;                               // xz
    float k7 = hi.w - hi.z - hi.y - lo.w + lo.y + hi.x + lo.z - lo.x;   // xyz
    float c0 = k0 + k1*a.x + k2*a.y + k3*a.z + k4*a.x*a.y + k5*a.y*a.z + k6*a.x*a.z + k7*a.x*a.y*a.z;
    float c1 = k1*b.x + k2*b.y + k3*b.z
             + k4*(a.x*b.y + b.x*a.y) + k5*(a.y*b.z + b.y*a.z) + k6*(a.x*b.z + b.x*a.z)
             + k7*(a.x*a.y*b.z + a.x*b.y*a.z + b.x*a.y*a.z);
    float c2 = k4*b.x*b.y + k5*b.y*b.z + k6*b.x*b.z + k7*(a.x*b.y*b.z + b.x*a.y*b.z + b.x*b.y*a.z);
    float c3 = k7*b.x*b.y*b.z;
    return vec4(c0, c1, c2, c3);
}

float cubic_eval(vec4 c, float s) {
    return c.x + s * (c.y + s * (c.z + s * c.w));
}

float cubic_deriv(vec4 c, float s) {
    return c.y + s * (2.0 * c.z + 3.0 * c.w * s);
}

// Where the cubic's derivative vanishes strictly inside (0, len), in order
int deriv_roots(vec4 c, float len, out float r0, out float r1) {
    float A = 3.0 * c.w;
    float B = 2.0 * c.z;
    float C = c.y;
    float scale = abs(A) + abs(B) + abs(C);
    float x0 = -1.0;
    float x1 = -1.0;
    if (scale > 0.0) {
        if (abs(A) <= 1e-7 * scale) {
            if (abs(B) > 1e-7 * scale) x0 = -C / B;
        } else {
            float disc = B * B - 4.0 * A * C;
            if (disc >= 0.0) {
                float sq = sqrt(disc);
                float q  = -0.5 * (B + (B >= 0.0 ? sq : -sq));
                x0 = q / A;
                x1 = (q != 0.0) ? C / q : x0;
                if (x1 < x0) { float tmp = x0; x0 = x1; x1 = tmp; }
            }
        }
    }
    int n = 0;
    r0 = len;
    r1 = len;
    if (x0 > 0.0 && x0 < len) { r0 = x0; n = 1; }
    if (x1 > 0.0 && x1 < len && x1 != x0) {
        if (n == 0) r0 = x1; else r1 = x1;
        n += 1;
    }
    return n;
}

// The crossing of g = cubic - v on [l, r], where g is monotone and its ends lie on different sides of v
// (or one is on it). Newton from the secant through the ends, kept inside the bracket by bisection.
float solve_crossing(vec4 c, float v, float l, float r, float gl, float gr) {
    if (gl == 0.0) return l;
    if (gr == 0.0) return r;
    if ((gl < 0.0) == (gr < 0.0)) return l;   // the change is at the start of the piece (it began on v)
    float s = clamp(l + (r - l) * gl / (gl - gr), l, r);
    for (int i = 0; i < 8; ++i) {
        float g = cubic_eval(c, s) - v;
        if (g == 0.0) return s;
        if ((g < 0.0) == (gl < 0.0)) { l = s; gl = g; } else { r = s; gr = g; }
        float dg = cubic_deriv(c, s);
        float sn = (dg != 0.0) ? s - g / dg : 0.5 * (l + r);
        sn = (sn > l && sn < r) ? sn : 0.5 * (l + r);
        if (abs(sn - s) < 1e-5) return sn;     // voxels
        s = sn;
    }
    return s;
}

bool is_inside(float d, float v) {
    return (v >= 0.0) ? (d >= v) : (d <= v);
}

uint membership(float d) {
    uint m = 0u;
    for (int i = 0; i < u_iso_count; ++i) {
        if (is_inside(d, u_iso_values[i])) m |= (1u << uint(i));
    }
    return m;
}

float tau_inside(uint mask) {
    float tau = 0.0;
    for (int i = 0; i < u_iso_count; ++i) {
        if ((mask & (1u << uint(i))) != 0u) tau += max(u_iso_tau[i], 0.0);
    }
    return tau;
}

// Whether some isovalue lies in [lo, hi] (which may have infinite ends)
bool range_holds_iso(float lo, float hi) {
    vec4 a = step(vec4(lo), g_iso_a) * step(g_iso_a, vec4(hi));
    vec4 b = step(vec4(lo), g_iso_b) * step(g_iso_b, vec4(hi));
    return any(greaterThan(a + b, vec4(0.0)));
}

// -----------------------------------------------------------------------------
// Shading: the model and light of compose_deferred.frag, so a surface with alpha 1 is lit like the
// opaque geometry around it (with the roughness the caller gives).
// -----------------------------------------------------------------------------

float FresnelSchlick(float cos_theta, float F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
}

float FresnelSchlickRoughness(float cos_theta, float F0, float roughness) {
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
}

float DistributionGGX(float NdotH, float roughness) {
    float a  = roughness * roughness;
    float a2 = a * a;
    float d  = NdotH * NdotH * (a2 - 1.0) + 1.0;
    return a2 / (PI * d * d);
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = roughness + 1.0;
    float k = (r * r) / 8.0;
    return NdotV / (NdotV * (1.0 - k) + k);
}

// The gradient along the view axes: its direction is the view space normal, no transform needed.
// Taken from the filtered field one voxel apart, which keeps the normals continuous across cells.
vec3 gradient_view(vec3 p) {
    vec3 dx = u_grad_offsets[0].xyz;
    vec3 dy = u_grad_offsets[1].xyz;
    vec3 dz = u_grad_offsets[2].xyz;
    return vec3(
        sample_volume(p + dx) - sample_volume(p - dx),
        sample_volume(p + dy) - sample_volume(p - dy),
        sample_volume(p + dz) - sample_volume(p - dz));
}

vec4 surface_color(vec3 p, int i) {
    vec4 base = u_iso_colors[i];
#if defined(USE_COLOR_VOLUME)
    return base * texture(u_tex_color_volume, p);
#elif defined(USE_FIELD)
    float f = texture(u_tex_field, p).r;
    float t = clamp((f - u_field_beg) * u_field_inv_ext, 0.0, 1.0);
    return base * vec4(texture(u_tex_field_colormap, vec2(t, 0.5)).rgb, 1.0);
#else
    return base;
#endif
}

// Adds one surface behind which the transmittance is g_T. The diffuse part is scaled by the surface's
// coverage (alpha); the specular part is not, so a clear surface still shows its reflections. The
// opacity is the coverage plus what Fresnel reflects away.
void shade_surface(vec3 albedo, float alpha, vec3 N, vec3 V) {
    float roughness = clamp(u_roughness, 0.04, 1.0);
    float F0 = clamp(u_F0, 0.0, 1.0);

    float NdotV = clamp(dot(N, V), 0.0, 1.0);
    float NdotL = clamp(dot(N, u_light_dir), 0.0, 1.0);

    vec3 diffuse  = vec3(0.0);
    vec3 specular = vec3(0.0);

    if (NdotL > 0.0) {
        vec3  H     = normalize(u_light_dir + V);
        float NdotH = clamp(dot(N, H), 0.0, 1.0);
        float HdotV = clamp(dot(H, V), 0.0, 1.0);
        float D = DistributionGGX(NdotH, roughness);
        float G = GeometrySchlickGGX(NdotV, roughness) * GeometrySchlickGGX(NdotL, roughness);
        float F = FresnelSchlick(HdotV, F0);
        specular += vec3(D * G * F / (4.0 * NdotV * NdotL + 0.0001)) * u_dir_radiance * NdotL;
        diffuse  += (1.0 - F) * albedo * INV_PI * u_dir_radiance * NdotL;
    }

    float Fe = FresnelSchlickRoughness(NdotV, F0, roughness);
    diffuse  += (1.0 - Fe) * albedo * INV_PI * u_env_radiance;
    specular += Fe * u_env_radiance;

    float Fv      = FresnelSchlick(NdotV, F0);
    float a       = clamp(alpha, 0.0, 1.0);
    float opacity = 1.0 - (1.0 - a) * (1.0 - Fv);

    g_L += g_T * (a * diffuse + specular);
    g_T *= 1.0 - opacity;
}

// -----------------------------------------------------------------------------
// Walking the ray
// -----------------------------------------------------------------------------

void absorb_to(float t) {
    if (g_tau > 0.0 && t > g_t_abs) {
        g_T *= exp(-g_tau * (t - g_t_abs) * g_ext_per_t);
    }
    g_t_abs = max(g_t_abs, t);
}

void hit(float t, int i) {
    absorb_to(t);
    vec3 p = g_p0 + t * g_ray;
    vec3 N = gradient_view(p);
    float nl = length(N);
    N = (nl > 1e-20) ? N / nl : g_V;
    if (dot(N, g_V) < 0.0) N = -N;
    vec4 c = surface_color(p, i);
    shade_surface(c.rgb, c.a, N, g_V);
    g_inside ^= (1u << uint(i));
    g_tau = tau_inside(g_inside);
}

// The part [tc, te] of the ray inside cell c, whose corners are lo / hi
void process_cell(ivec3 c, vec4 lo, vec4 hi, float tc, float te) {
    vec4  mn4  = min(lo, hi);
    vec4  mx4  = max(lo, hi);
    bool  cand = range_holds_iso(min(min(mn4.x, mn4.y), min(mn4.z, mn4.w)), max(max(mx4.x, mx4.y), max(mx4.z, mx4.w)));
    if (!cand && g_known) return;

    vec3  a   = (g_o + g_d * tc) - vec3(c);     // local coordinates where the ray enters the cell part
    float len = (te - tc) * g_dlen;            // its length in voxels

    if (!g_known) {
        // In a cell whose range holds no isovalue every corner is on the side of every surface the ray is on
        g_inside = membership(cand ? cubic_eval(cubic_coeffs(lo, hi, a, g_dir), 0.0) : lo.x);
        g_tau    = tau_inside(g_inside);
        g_known  = true;
        if (!cand) return;
    }

    vec4 cf = cubic_coeffs(lo, hi, a, g_dir);

    // Split into monotone pieces
    float r0, r1;
    int nr = deriv_roots(cf, len, r0, r1);
    float s[4];
    s[0] = 0.0;
    s[1] = (nr > 0) ? r0 : len;
    s[2] = (nr > 1) ? r1 : len;
    s[3] = len;
    int np = nr + 1;
    float f[4];
    f[0] = cubic_eval(cf, s[0]);
    f[1] = cubic_eval(cf, s[1]);
    f[2] = cubic_eval(cf, s[2]);
    f[3] = cubic_eval(cf, len);
    // f[np] must be the end of the last piece
    if (np == 1) f[1] = f[3];
    if (np == 2) f[2] = f[3];

    // Every crossing in the cell, in order along the ray. The side at the start of the cell is the one
    // carried along the ray, not re-evaluated, so a crossing on a cell face is counted exactly once.
    float hs[MAX_CELL_HITS];
    int   hi_[MAX_CELL_HITS];
    int   m = 0;
    for (int i = 0; i < u_iso_count; ++i) {
        float v = u_iso_values[i];
        bool flag = (g_inside & (1u << uint(i))) != 0u;
        for (int j = 0; j < np; ++j) {
            bool fe = is_inside(f[j + 1], v);
            if (fe != flag && m < MAX_CELL_HITS) {
                float h = solve_crossing(cf, v, s[j], s[j + 1], f[j] - v, f[j + 1] - v);
                int k = m;
                while (k > 0 && hs[k - 1] > h) {
                    hs[k]  = hs[k - 1];
                    hi_[k] = hi_[k - 1];
                    --k;
                }
                hs[k]  = h;
                hi_[k] = i;
                ++m;
            }
            flag = fe;
        }
    }

    for (int k = 0; k < m; ++k) {
        hit(tc + hs[k] / g_dlen, hi_[k]);
        if (g_T < T_MIN) return;
    }
}

// The corners of cell c: the four new ones after a step into it along axis prev (0-2) with the corners of
// the previous cell in lo / hi, all eight otherwise
void load_corners(int prev, ivec3 st, ivec3 c, inout vec4 lo, inout vec4 hi) {
    if      (prev == 0) step_corners_x(st.x, c, lo, hi);
    else if (prev == 1) step_corners_y(st.y, c, lo, hi);
    else if (prev == 2) step_corners_z(st.z, c, lo, hi);
    else                fetch_corners(c, lo, hi);
}

// Every cell the ray passes through in [t0, t1]. The range of a cell comes from the cell min/max grid, one
// fetch; only a cell whose range holds an isovalue loads its corners (reusing the face it shares with the
// previous cell when that one loaded them) and is intersected. Without the grid every cell loads them.
void walk_cells(float t0, float t1, float eps_t, vec3 inv, bvec3 flat_axis) {
    vec3  xs = g_o + g_d * (t0 + eps_t);
    ivec3 c  = clamp(ivec3(floor(xs)), ivec3(-1), g_dim_m1);
    ivec3 st = ivec3(mix(vec3(-1.0), vec3(1.0), greaterThanEqual(g_d, vec3(0.0))));
    vec3  nb = vec3(c) + vec3(greaterThanEqual(g_d, vec3(0.0)));        // the next boundary on each axis
    vec3  tn = mix((nb - g_o) * inv, vec3(HUGE), flat_axis);
    vec3  td = mix(abs(inv), vec3(HUGE), flat_axis);

    bool use_cells = u_use_cells > 0.5;
    vec4 lo = vec4(0.0);
    vec4 hi = vec4(0.0);
    int  prev = -1;     // lo / hi hold the corners of the cell before the last step, which was along this axis

    float tc = t0;
    for (int k = 0; k < MAX_STEPS; ++k) {
        float te = min(min(tn.x, tn.y), min(tn.z, t1));
        bool  loaded = false;
        if (te > tc) {
            bool cand = true;
            if (use_cells) {
                vec2 r = texelFetch(u_tex_cells, clamp(c + 1, ivec3(0), g_dim_m1 + 2), 0).xy;
                cand = range_holds_iso(r.x, r.y);
                if (!cand && !g_known) {
                    g_inside = membership(r.x);
                    g_tau    = tau_inside(g_inside);
                    g_known  = true;
                }
            }
            if (cand) {
                load_corners(prev, st, c, lo, hi);
                loaded = true;
                process_cell(c, lo, hi, tc, te);
                if (g_T < T_MIN) return;
            }
        }
        if (te >= t1) return;
        int axis;
        if (tn.x <= tn.y && tn.x <= tn.z) {
            c.x += st.x; tn.x += td.x; axis = 0;
        } else if (tn.y <= tn.z) {
            c.y += st.y; tn.y += td.y; axis = 1;
        } else {
            c.z += st.z; tn.z += td.z; axis = 2;
        }
        prev = loaded ? axis : -1;
        tc = max(tc, te);
    }
}

void main() {
    ivec2 px  = ivec2(gl_FragCoord.xy);
    vec2  ndc = (gl_FragCoord.xy * u_inv_res) * 2.0 - 1.0;

    float proxy_entry = 0.0;
    float proxy_exit  = 1.0;
    if (u_use_proxy > 0.5) {
        proxy_exit = texelFetch(u_tex_exit, px, 0).r;
        if (proxy_exit <= 0.0) discard;     // the ray meets no block that can hold a surface
        if (u_entry_from_near < 0.5) {
            proxy_entry = texelFetch(u_tex_entry, px, 0).r;
        }
    }

    // The ray from the near plane to the opaque scene (or the far plane), clipped to the clip box
    float depth = (u_use_depth > 0.5) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
    if (proxy_entry >= depth) discard;
    vec4 pn4 = u_clip_to_model * vec4(ndc, -1.0, 1.0);
    vec4 pf4 = u_clip_to_model * vec4(ndc, depth * 2.0 - 1.0, 1.0);
    vec3 pn  = pn4.xyz / pn4.w;
    vec3 pf  = pf4.xyz / pf4.w;
    vec3 seg = pf - pn;

    vec3 inv_seg = 1.0 / mix(seg, vec3(1e-20), equal(seg, vec3(0.0)));
    vec3 t0   = (u_clip_min - pn) * inv_seg;
    vec3 t1   = (u_clip_max - pn) * inv_seg;
    vec3 tmin = min(t0, t1);
    vec3 tmax = max(t0, t1);
    float s0 = max(max(tmin.x, tmin.y), max(tmin.z, 0.0));
    float s1 = min(min(tmax.x, tmax.y), min(tmax.z, 1.0));

    vec3 dim = vec3(textureSize(u_tex_volume, 0));
    if (u_use_proxy > 0.5) {
        // The proxy span as parameters on the same segment, widened by half a voxel so that depth
        // quantization cannot start a ray just past a surface that touches a block face
        float margin = 0.5 / max(length(seg * dim), 1e-6);
        float inv_ss = 1.0 / max(dot(seg, seg), 1e-30);
        if (proxy_entry > 0.0) {
            vec4 pe = u_clip_to_model * vec4(ndc, proxy_entry * 2.0 - 1.0, 1.0);
            s0 = max(s0, dot(pe.xyz / pe.w - pn, seg) * inv_ss - margin);
        }
        if (proxy_exit < depth) {
            vec4 px4 = u_clip_to_model * vec4(ndc, proxy_exit * 2.0 - 1.0, 1.0);
            s1 = min(s1, dot(px4.xyz / px4.w - pn, seg) * inv_ss + margin);
        }
    }
    if (s0 >= s1) discard;

    g_dim_m1 = ivec3(dim) - 1;
    g_p0  = pn + seg * s0;
    g_ray = seg * (s1 - s0);
    g_o   = g_p0 * dim - 0.5;
    g_d   = g_ray * dim;
    g_dlen = length(g_d);
    if (g_dlen < 1e-6) discard;
    g_dir = g_d / g_dlen;
    g_ext_per_t = length(mat3(u_model_to_view) * g_ray) * u_optical_scale;
    g_V = -normalize(mat3(u_model_to_view) * g_ray);

    g_L = vec3(0.0);
    g_T = 1.0;
    g_t_abs = 0.0;
    g_inside = 0u;
    g_tau = 0.0;
    g_known = false;
    for (int i = 0; i < 8; ++i) {
        float v = u_iso_values[clamp(i, 0, max(u_iso_count - 1, 0))];
        if (i < 4) g_iso_a[i] = v; else g_iso_b[i - 4] = v;
    }

    bvec3 flat_axis = lessThan(abs(g_d), vec3(1e-20));
    vec3  inv   = 1.0 / mix(g_d, vec3(1.0), flat_axis);
    float eps_t = 1e-3 / max(max(abs(g_d.x), abs(g_d.y)), abs(g_d.z));   // a thousandth of a voxel

    float B = u_block_size;
    ivec3 bdim = textureSize(u_tex_minmax, 0);

    float t = 0.0;
    for (int k = 0; k < MAX_STEPS && t < 1.0 && g_T >= T_MIN; ++k) {
        if (u_use_blocks < 0.5) {
            walk_cells(0.0, 1.0, eps_t, inv, flat_axis);
            break;
        }

        // The block just ahead of t, and where the ray leaves it
        vec3  xb  = (g_o + 0.5 + g_d * (t + eps_t)) / B;
        ivec3 b   = clamp(ivec3(floor(xb)), ivec3(0), bdim - 1);
        vec3  bnd = vec3(b) * B - 0.5 + B * vec3(greaterThanEqual(g_d, vec3(0.0)));
        vec3  tx  = mix((bnd - g_o) * inv, vec3(HUGE), flat_axis);
        float tb  = min(1.0, min(min(tx.x, tx.y), tx.z));
        tb = max(tb, t + eps_t);

        vec2 mm = texelFetch(u_tex_minmax, b, 0).xy;
        if (range_holds_iso(mm.x, mm.y)) {
            walk_cells(t, tb, eps_t, inv, flat_axis);
        } else {
            // No surface in the block: the ray is on one side of every surface throughout it
            absorb_to(t);
            g_inside = membership(mm.x);
            g_tau    = tau_inside(g_inside);
            g_known  = true;
        }
        t = tb;
    }
    absorb_to(1.0);

    out_color = vec4(g_L, 1.0 - g_T);
}
