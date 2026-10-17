#version 410 core

// Isosurface raycaster, fast mode. The exact mode, which intersects every cell the ray passes through, is
// isosurface.frag; both share the uniforms below (u_samples_per_voxel is u_block_size there).
//
// Marches the density volume at one sample per voxel along the ray, finds where the ray crosses each
// isovalue and refines the crossing on the trilinear field, so the position of a crossing does not depend
// on where the samples happen to fall (no jitter, no temporal accumulation needed). What is not found is
// a surface the ray enters and leaves within one step: thinner than a voxel along the ray, which happens
// at grazing silhouettes. Every crossing is shaded with the model and light of the deferred compose pass,
// and the result is written as premultiplied linear radiance + coverage, to be blended over the HDR
// scene before tone mapping.
//
// Membership: a point is INSIDE isosurface i when the density is on the far side of its value from
// zero, d >= v for v >= 0 and d <= v for v < 0. Crossing a surface toggles membership of exactly that
// surface, and the optical densities of the surfaces enclosing the ray attenuate it in between.
//
// Empty space: when u_use_proxy is set, the nearest and farthest depths of the blocks that can hold a
// surface (block_proxy.vert) narrow the ray to the part that can see one, and pixels whose ray meets no
// such block are discarded before any sampling.

#ifndef MAX_ISO
#define MAX_ISO 8
#endif

layout(std140) uniform IsoUniforms {
    mat4  u_clip_to_model;      // inverse(proj * view * model)
    mat4  u_model_to_view;
    mat4  u_grad_offsets;       // columns 0-2: model space offsets for one voxel along view x, y, z

    vec3  u_clip_min;
    float u_samples_per_voxel;
    vec3  u_clip_max;
    float u_use_depth;

    vec3  u_sky_radiance;       // the environment above the horizon of the view (y > 0)
    float u_roughness;
    vec3  u_dir_radiance;       // the key light: irradiance at normal incidence
    float u_F0;
    vec3  u_light_dir;          // toward the key light, view space, normalized
    float u_field_beg;

    vec2  u_inv_res;
    float u_field_inv_ext;
    float u_optical_scale;      // optical density -> extinction per world unit

    float u_use_proxy;          // u_tex_entry / u_tex_exit hold the block proxy depths
    float u_entry_from_near;    // the near plane cuts the volume: rays start there, not at the proxy entry
    float u_pad0;               // u_use_blocks of the exact mode
    float u_pad1;

    vec3  u_ground_radiance;    // the environment below the horizon of the view
    float u_basic;              // 1: shade_basic, the basic model of the compose pass; the radiances are unused

    float u_basic_ambient;
    float u_basic_specular;
    float u_linear_colors;      // 1: the colours are linear reflectances, 0: sRGB (as in the compose pass)
    float u_light_strength;     // of the directional light, relative to u_dir_radiance (and to the basic light)
};

uniform float u_iso_values[MAX_ISO];
uniform vec4  u_iso_colors[MAX_ISO];
uniform float u_iso_tau[MAX_ISO];
uniform int   u_iso_count;

uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_depth;
uniform sampler2D u_tex_entry;  // nearest depth of the block proxy, 1 where there is none
uniform sampler2D u_tex_exit;   // farthest depth of the block proxy, 0 where there is none

#if defined(USE_COLOR_VOLUME)
uniform sampler3D u_tex_color_volume;
#elif defined(USE_FIELD)
// A scalar field over the density volume's box (at its own resolution, sampled by the same
// normalised coordinates) and the colour map it is shown through
uniform sampler3D u_tex_field;
uniform sampler2D u_tex_field_colormap;
#endif

#if defined(USE_FIELD)
// Isolines (gfx/isolines.h): lines on the surfaces where the field takes these values, the persistent ones and
// the transient one, u_isoline_width pixels wide.
//
// A hit is found inside the loops along the ray, where screen space derivatives are undefined, so the field's
// change per pixel comes from ray differentials (Igehy 1999) instead: the rays through the next pixel along x and
// along y meet the tangent plane of the hit, which gives how far the surface point moves per pixel (dP/dx, dP/dy),
// and the field changes by its gradient along that, df/dx = grad f . dP/dx (grad f need not be projected onto
// the tangent plane: dP/dx lies in it). The gradient is taken from the field volume by central differences over
// the same offsets as the normals; the band the field is evaluated in is wide enough for them (surface_field.cpp,
// band_mask). The distance to the line of value v is then |f - v| / |(df/dx, df/dy)| pixels, and the lines keep
// their width in pixels however the surface is turned. The mesh shader does the same (isomesh.frag).
#define ISOLINES_MAX 8
uniform int   u_isoline_count;
uniform float u_isoline_values[ISOLINES_MAX];
uniform int   u_isoline_hover;
uniform float u_isoline_hover_value;
uniform float u_isoline_width;
uniform float u_field_grad_scale;   // 1 / (2 h), h the world length of the u_grad_offsets steps

bool g_lines = false;   // whether there are any lines to draw
vec3 g_rx_o;            // the rays through the next pixel along x and along y: origin and direction, view space
vec3 g_rx_d;
vec3 g_ry_o;
vec3 g_ry_d;

// The view space point at normalized device coordinates (ndc, z)
vec3 isoline_view_point(vec2 ndc, float z) {
    vec4 m = u_clip_to_model * vec4(ndc, z, 1.0);
    return (u_model_to_view * vec4(m.xyz / m.w, 1.0)).xyz;
}

void isoline_init(vec2 ndc) {
    g_lines = u_isoline_count > 0 || u_isoline_hover != 0;
    if (!g_lines) return;
    vec2 nx = ndc + vec2(2.0 * u_inv_res.x, 0.0);
    vec2 ny = ndc + vec2(0.0, 2.0 * u_inv_res.y);
    g_rx_o = isoline_view_point(nx, -1.0);
    g_rx_d = isoline_view_point(nx,  0.0) - g_rx_o;
    g_ry_o = isoline_view_point(ny, -1.0);
    g_ry_d = isoline_view_point(ny,  0.0) - g_ry_o;
}

// Where the ray o + s d meets the plane through P with normal N. At a grazing angle the point runs off far: the
// field changes a lot from one pixel to the next there, which is what it does on the screen.
vec3 isoline_plane_hit(vec3 o, vec3 d, vec3 P, vec3 N) {
    float dn = dot(N, d);
    float lo = 1.0e-6 * length(d);
    dn = (dn < 0.0) ? min(dn, -lo) : max(dn, lo);
    return o + d * (dot(N, P - o) / dn);
}

// The coverage of a line w pixels wide at distance d (pixels) from its middle, box filtered over a pixel
float isoline_coverage(float d, float w) {
    return clamp(0.5 * w + 0.5 - d, 0.0, 1.0);
}

// The isolines painted onto the colour of a surface: f the field there, g its change per pixel along the
// surface. The persistent lines dark with a light rim, the transient one light with a dark rim; a line is
// opaque paint, whatever the coverage of the surface under it. As in isomesh.frag.
vec4 apply_isolines(vec4 c, float f, float g) {
    const vec3  DARK  = vec3(0.02);
    const vec3  LIGHT = vec3(0.95);
    const float RIM   = 2.0;            // pixels, both sides together
    float inv_g = 1.0 / max(g, 1.0e-30);
    float w = u_isoline_width;

    float core = 0.0;
    float rim  = 0.0;
    for (int i = 0; i < u_isoline_count; ++i) {
        float d = abs(f - u_isoline_values[i]) * inv_g;
        core = max(core, isoline_coverage(d, w));
        rim  = max(rim,  isoline_coverage(d, w + RIM));
    }
    c.rgb = mix(c.rgb, LIGHT, rim);
    c.rgb = mix(c.rgb, DARK, core);
    float cover = rim;

    if (u_isoline_hover != 0) {
        float d  = abs(f - u_isoline_hover_value) * inv_g;
        float hc = isoline_coverage(d, w);
        float hr = isoline_coverage(d, w + RIM);
        c.rgb = mix(c.rgb, DARK, hr);
        c.rgb = mix(c.rgb, LIGHT, hc);
        cover = max(cover, hr);
    }
    c.a = mix(c.a, 1.0, cover);
    return c;
}

// The lines on the surface hit at p (model space), N its view space normal, c its colour
vec4 surface_isolines(vec4 c, vec3 p, vec3 N) {
    vec3 ox = u_grad_offsets[0].xyz;
    vec3 oy = u_grad_offsets[1].xyz;
    vec3 oz = u_grad_offsets[2].xyz;
    vec3 G = vec3(
        texture(u_tex_field, p + ox).r - texture(u_tex_field, p - ox).r,
        texture(u_tex_field, p + oy).r - texture(u_tex_field, p - oy).r,
        texture(u_tex_field, p + oz).r - texture(u_tex_field, p - oz).r) * u_field_grad_scale;
    vec3  P    = (u_model_to_view * vec4(p, 1.0)).xyz;
    vec3  dPdx = isoline_plane_hit(g_rx_o, g_rx_d, P, N) - P;
    vec3  dPdy = isoline_plane_hit(g_ry_o, g_ry_d, P, N) - P;
    float g    = length(vec2(dot(G, dPdx), dot(G, dPdy)));
    return apply_isolines(c, texture(u_tex_field, p).r, g);
}
#endif

layout(location = 0) out vec4 out_color;

const float PI      = 3.1415926535;
const float INV_PI  = 1.0 / PI;
const float T_MIN   = 0.005;    // early ray termination: below this transmittance the ray is opaque
const int   REFINE_STEPS = 2;   // regula falsi iterations per crossing

float sample_volume(vec3 p) {
    return texture(u_tex_volume, p).r;
}

bool is_inside(float d, float v) {
    return (v >= 0.0) ? (d >= v) : (d <= v);
}

vec3 unproject(vec2 ndc_xy, float depth) {
    vec4 p = u_clip_to_model * vec4(ndc_xy, depth * 2.0 - 1.0, 1.0);
    return p.xyz / p.w;
}

// Where along [pa, pb] (as a fraction) the field crosses v, given the values at the ends on opposite
// sides of it. Regula falsi on the trilinear field, then a final linear estimate within the bracket.
float refine_crossing(vec3 pa, vec3 pb, float da, float db, float v) {
    float a  = 0.0;
    float b  = 1.0;
    float fa = da - v;
    float fb = db - v;
    for (int i = 0; i < REFINE_STEPS; ++i) {
        float denom = fa - fb;
        if (denom == 0.0) break;
        float t  = a + (b - a) * fa / denom;
        float ft = sample_volume(mix(pa, pb, t)) - v;
        if (ft == 0.0) return t;        // on it; also keeps a bracket end from landing on the value
        if ((ft < 0.0) == (fa < 0.0)) {
            a  = t;
            fa = ft;
        } else {
            b  = t;
            fb = ft;
        }
    }
    float denom = fa - fb;
    return (denom != 0.0) ? a + (b - a) * fa / denom : 0.5 * (a + b);
}

// The gradient along the view axes: its direction is the view space normal, no transform needed
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

float tau_inside(uint mask) {
    float tau = 0.0;
    for (int i = 0; i < u_iso_count; ++i) {
        if ((mask & (1u << uint(i))) != 0u) {
            tau += max(u_iso_tau[i], 0.0);
        }
    }
    return tau;
}

// -----------------------------------------------------------------------------
// Shading: the model and lights of compose_deferred.frag, so a surface with alpha 1 is lit like the
// opaque geometry around it (with the roughness the caller gives): its colour decoded from sRGB to the
// linear albedo, a key light and an environment of sky and ground.
// -----------------------------------------------------------------------------

float FresnelSchlick(float cos_theta, float F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
}

// The directional albedo of the specular lobe, F0 * x + y: Karis, "Physically Based Shading on Mobile" (2014)
vec2 EnvBRDFApprox(float NdotV, float roughness) {
    const vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
    const vec4 c1 = vec4( 1.0,  0.0425,  1.04, -0.04);
    vec4  r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    return vec2(-1.04, 1.04) * a004 + r.zw;
}

// The irradiance of the sky and ground hemispheres over pi
vec3 environment_irradiance(vec3 N) {
    return mix(u_ground_radiance, u_sky_radiance, 0.5 + 0.5 * N.y);
}

// Their radiance about the direction R, as seen through a lobe of the roughness
vec3 environment_radiance(vec3 R, float roughness) {
    float w = 0.1 + roughness;
    return mix(u_ground_radiance, u_sky_radiance, smoothstep(-w, w, R.y));
}

// The colours are picked as sRGB
vec3 srgb_to_linear(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c));
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

// The basic model of compose_deferred.frag (shade_basic): the colour as given where the surface faces the
// light, darkened toward u_basic_ambient of itself away from it, a highlight blended toward white, all on
// the display (sRGB) values, returned as linear radiance.
const float BASIC_SHININESS = 48.0;

vec3 shade_basic(vec3 color, vec3 N, vec3 V) {
    float w = 0.5 + 0.5 * dot(N, u_light_dir);
    vec3  c = clamp(color, 0.0, 1.0) * (u_basic_ambient + u_light_strength * (1.0 - u_basic_ambient) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));
    vec3  H = normalize(u_light_dir + V);
    float s = u_light_strength * u_basic_specular * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = mix(c, vec3(1.0), clamp(s, 0.0, 1.0));
    return (u_linear_colors > 0.5) ? c : srgb_to_linear(c);
}

// Adds one surface to the premultiplied radiance L behind which the transmittance is T.
// The diffuse part is scaled by the surface's coverage (alpha); the specular part is not, so a clear
// surface still shows its reflections. The opacity is the coverage plus what the surface reflects away.
// In the basic model the surface is its shaded colour over its coverage.
void shade_surface(inout vec3 L, inout float T, vec3 color, float alpha, vec3 N, vec3 V) {
    if (u_basic > 0.5) {
        float a = clamp(alpha, 0.0, 1.0);
        L += T * a * shade_basic(color, N, V);
        T *= 1.0 - a;
        return;
    }

    float roughness = clamp(u_roughness, 0.04, 1.0);
    float F0 = clamp(u_F0, 0.0, 1.0);
    vec3  albedo = (u_linear_colors > 0.5) ? clamp(color, 0.0, 1.0) : srgb_to_linear(color);

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
        vec3  E = u_dir_radiance * u_light_strength;
        specular += vec3(D * G * F / (4.0 * NdotV * NdotL + 0.0001)) * E * NdotL;
        diffuse  += (1.0 - F) * albedo * INV_PI * E * NdotL;
    }

    // What the surface reflects of its surroundings is the directional albedo of its specular lobe:
    // that much of the environment is reflected and that much less passes through
    vec2  ab = EnvBRDFApprox(NdotV, roughness);
    float Fe = F0 * ab.x + ab.y;
    diffuse  += (1.0 - Fe) * albedo * environment_irradiance(N);
    specular += Fe * environment_radiance(reflect(-V, N), roughness);

    float a       = clamp(alpha, 0.0, 1.0);
    float opacity = 1.0 - (1.0 - a) * (1.0 - Fe);

    L += T * (a * diffuse + specular);
    T *= 1.0 - opacity;
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
    vec3 pn  = unproject(ndc, 0.0);
    vec3 pf  = unproject(ndc, depth);
    vec3 seg = pf - pn;

    vec3 inv_seg = 1.0 / mix(seg, vec3(1e-20), equal(seg, vec3(0.0)));
    vec3 t0   = (u_clip_min - pn) * inv_seg;
    vec3 t1   = (u_clip_max - pn) * inv_seg;
    vec3 tmin = min(t0, t1);
    vec3 tmax = max(t0, t1);
    float s0 = max(max(tmin.x, tmin.y), max(tmin.z, 0.0));
    float s1 = min(min(tmax.x, tmax.y), min(tmax.z, 1.0));

    if (u_use_proxy > 0.5) {
        // The proxy span as parameters on the same segment, widened by half a voxel so that depth
        // quantization cannot start a ray just past a surface that touches a block face
        vec3  dim_v  = vec3(textureSize(u_tex_volume, 0));
        float margin = 0.5 / max(length(seg * dim_v), 1e-6);
        float inv_ss = 1.0 / max(dot(seg, seg), 1e-30);
        if (proxy_entry > 0.0) {
            s0 = max(s0, dot(unproject(ndc, proxy_entry) - pn, seg) * inv_ss - margin);
        }
        if (proxy_exit < depth) {
            s1 = min(s1, dot(unproject(ndc, proxy_exit) - pn, seg) * inv_ss + margin);
        }
    }
    if (s0 >= s1) discard;

    vec3  p0  = pn + seg * s0;
    vec3  ray = seg * (s1 - s0);
    float len = length(ray);
    if (len < 1e-6) discard;

    vec3  dir   = ray / len;
    vec3  dim   = vec3(textureSize(u_tex_volume, 0));
    int   n     = max(1, int(ceil(len * length(dir * dim) * u_samples_per_voxel)));
    vec3  dp    = ray / float(n);
    float ext_step = length(mat3(u_model_to_view) * dp) * u_optical_scale;  // extinction per unit tau over one step
    vec3  V     = -normalize(mat3(u_model_to_view) * dir);

#if defined(USE_FIELD)
    isoline_init(ndc);
#endif

    float d0 = sample_volume(p0);
    uint inside = 0u;
    for (int i = 0; i < u_iso_count; ++i) {
        if (is_inside(d0, u_iso_values[i])) inside |= (1u << uint(i));
    }
    float tau = tau_inside(inside);

    vec3  L = vec3(0.0);
    float T = 1.0;
    vec3  pa = p0;

    for (int s = 1; s <= n; ++s) {
        vec3  pb = p0 + dp * float(s);
        float d1 = sample_volume(pb);

        uint crossed = 0u;
        for (int i = 0; i < u_iso_count; ++i) {
            if (is_inside(d0, u_iso_values[i]) != is_inside(d1, u_iso_values[i])) crossed |= (1u << uint(i));
        }

        if (crossed == 0u) {
            if (tau > 0.0) T *= exp(-tau * ext_step);
        } else {
            // Locate every crossing in this step and take them in order along the ray
            float frac[MAX_ISO];
            int   idx[MAX_ISO];
            int   m = 0;
            for (int i = 0; i < u_iso_count; ++i) {
                if ((crossed & (1u << uint(i))) == 0u) continue;
                float f = clamp(refine_crossing(pa, pb, d0, d1, u_iso_values[i]), 0.0, 1.0);
                int j = m;
                while (j > 0 && frac[j - 1] > f) {
                    frac[j] = frac[j - 1];
                    idx[j]  = idx[j - 1];
                    --j;
                }
                frac[j] = f;
                idx[j]  = i;
                ++m;
            }

            float last = 0.0;
            for (int h = 0; h < m; ++h) {
                if (tau > 0.0) T *= exp(-tau * ext_step * (frac[h] - last));
                last = frac[h];

                int  i  = idx[h];
                vec3 ph = mix(pa, pb, frac[h]);
                vec3 N  = gradient_view(ph);
                float nl = length(N);
                N = (nl > 1e-20) ? N / nl : V;
                if (dot(N, V) < 0.0) N = -N;

                vec4 c = surface_color(ph, i);
#if defined(USE_FIELD)
                if (g_lines) c = surface_isolines(c, ph, N);
#endif
                shade_surface(L, T, c.rgb, c.a, N, V);

                inside ^= (1u << uint(i));
                tau = tau_inside(inside);
            }
            if (tau > 0.0) T *= exp(-tau * ext_step * (1.0 - last));
        }

        if (T < T_MIN) break;
        pa = pb;
        d0 = d1;
    }

    out_color = vec4(L, 1.0 - T);
}
