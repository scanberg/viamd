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
    float u_pad0;               // u_use_blocks of the exact mode
    float u_pad1;
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
// A scalar field on the density volume's own grid and the colour map it is shown through
uniform sampler3D u_tex_field;
uniform sampler2D u_tex_field_colormap;
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

// Adds one surface to the premultiplied radiance L behind which the transmittance is T.
// The diffuse part is scaled by the surface's coverage (alpha); the specular part is not, so a clear
// surface still shows its reflections. The opacity is the coverage plus what Fresnel reflects away.
void shade_surface(inout vec3 L, inout float T, vec3 albedo, float alpha, vec3 N, vec3 V) {
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
