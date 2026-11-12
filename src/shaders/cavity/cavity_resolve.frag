#version 410 core

// The cavity's full screen pass over the G-buffer (cavity.h): where a sphere of the cavity was drawn, the exact point
// on it, the tessera nearest to it, and the field there through the colour map. Writes the colour (its alpha 1/255
// when unlit, which the composition reads as such and the coverage as geometry) and the picking index of the tessera.
//
// Of a switched cavity, the surface shows the weight GOSTSHYP gives its area there, the switching of VeloxChem's
// 'swig': the points of each sphere i are weighted by
//   F_i(x) = prod_j f((|x' - X_j| - (r_j - alpha gamma r_j)) / (gamma r_j)),   f(t) = 10t^3 - 15t^4 + 6t^5 on [0, 1]
// with x' the point on the sphere grown by r_ext, r_j = R_j + r_ext, gamma = sqrt(14 / N_lebedev) and
// alpha = 1/2 + 1/gamma - sqrt(1/gamma^2 - 1/28), and the area of a tessera is its Lebedev share times F. F is below 1
// over most of the spheres (the bands are about a van der Waals radius wide), so what is marked is where it is below a
// level, by default 1/2 (the bands along the seams that a cut at F = 1/2 would remove): screened by light dots whose
// area is the share missing below it, 1 - F / u_level, so that the surface washes out as it counts less, and
// optionally outlined by a light line where it is the level, an isoline of the weight.
// In screen space, a halftone on a grid at 45 degrees: on the surface its dots would stretch with the slope.

uniform sampler2D      u_tex_picking;   // a copy of the G-buffer's picking target
uniform samplerBuffer  u_tess_a;        // xyz: the centre of the tessera (the frame of the file), w: the kernel's sigma
uniform samplerBuffer  u_tess_b;        // xyz: its normal, w: its own value
uniform samplerBuffer  u_tess_c;        // x: what its kernel spreads per its normalization, y: its area per the same,
                                        // z: its weight in the power diagram (Å^2)
uniform samplerBuffer  u_spheres;       // xyz: the centre of the sphere, w: its radius
uniform usamplerBuffer u_nbr_range;     // x: the first of the sphere's neighbours in u_nbr_items, y: how many
uniform usamplerBuffer u_nbr_items;
uniform usamplerBuffer u_cells;         // x: the first of the cell's tesserae in u_items, y: how many
uniform usamplerBuffer u_items;
uniform sampler2D      u_colormap;

uniform mat4  u_clip_to_cavity;
uniform vec2  u_inv_size;
uniform uint  u_sphere_base;
uniform uint  u_sphere_count;
uniform uint  u_tessera_base;
uniform vec3  u_grid_min;
uniform float u_inv_cell;
uniform ivec3 u_grid_dim;
uniform int   u_mode;                   // 0: nearest, 1: kernel averaged by area, 2: kernel conservative
uniform vec2  u_range;
uniform int   u_unlit;
uniform int   u_edges;
uniform float u_power;                  // 1: nearest by power, |x - c_j|^2 - w_j; 0: by distance
uniform float u_power_reach2;           // a tessera further than this (squared) is never nearest by power
uniform int   u_switched;
uniform float u_gamma;
uniform float u_alpha;
uniform float u_r_ext;
uniform float u_screen_period;         // pixels, of the halftone
uniform float u_level;                  // the weight below which the surface is masked
uniform int   u_outline;                // and a line drawn where it is

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_picking;

const float KERNEL_CUT = 3.0;

uint decode_index(vec4 c) {
    uvec4 b = uvec4(round(c * 255.0));
    return b.x | (b.y << 8u) | (b.z << 16u) | (b.w << 24u);
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

float smoothstep5(float t) {
    t = clamp(t, 0.0, 1.0);
    return t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
}

// The switching of the point of sphere i in the direction n
float switching(int i, vec3 n) {
    vec4  s  = texelFetch(u_spheres, i);
    vec3  pe = s.xyz + n * (s.w + u_r_ext);
    uvec2 r  = texelFetch(u_nbr_range, i).xy;
    float F  = 1.0;
    for (uint k = r.x; k < r.x + r.y && F > 0.0; ++k) {
        vec4  sj = texelFetch(u_spheres, int(texelFetch(u_nbr_items, int(k)).r));
        float rj = sj.w + u_r_ext;
        F *= smoothstep5((length(pe - sj.xyz) - (rj - u_alpha * u_gamma * rj)) / (u_gamma * rj));
    }
    return F;
}

void ray(vec2 frag, out vec3 ro, out vec3 rd) {
    vec2 ndc = frag * u_inv_size * 2.0 - 1.0;
    vec4 a = u_clip_to_cavity * vec4(ndc, -1.0, 1.0);
    vec4 b = u_clip_to_cavity * vec4(ndc,  1.0, 1.0);
    ro = a.xyz / a.w;
    rd = normalize(b.xyz / b.w - ro);
}

void main() {
    uint idx = decode_index(texelFetch(u_tex_picking, ivec2(gl_FragCoord.xy), 0));
    if (idx < u_sphere_base || idx >= u_sphere_base + u_sphere_count) discard;
    uint sphere = idx - u_sphere_base;

    vec4 s = texelFetch(u_spheres, int(sphere));
    vec3 ro, rd;
    ray(gl_FragCoord.xy, ro, rd);
    // The ray hits the sphere it was drawn as; at its rim, by rounding, it may graze it: the nearest point then
    vec3  oc = ro - s.xyz;
    float b  = dot(oc, rd);
    float h  = b * b - (dot(oc, oc) - s.w * s.w);
    float t  = -b - sqrt(max(h, 0.0));
    vec3  x  = ro + rd * t;
    vec3  n  = normalize(x - s.xyz);
    x = s.xyz + n * s.w;

    float weight = (u_switched != 0) ? switching(int(sphere), n) : 1.0;
    // How much the weight changes over a pixel, for the line: here, before anything else can discard a pixel of
    // the quad. Across the border of two spheres it is not a change of the weight: no line there.
    float weight_fw = fwidth(weight);

    ivec3 c = ivec3(floor((x - u_grid_min) * u_inv_cell));
    if (any(lessThan(c, ivec3(0))) || any(greaterThanEqual(c, u_grid_dim))) discard;
    uvec2 range = texelFetch(u_cells, c.x + u_grid_dim.x * (c.y + u_grid_dim.y * c.z)).xy;

    float d1 = 1e30, d2 = 1e30;
    int   j1 = -1;
    vec3  c1 = vec3(0), c2 = vec3(0);
    float num = 0.0, den = 0.0;
    for (uint k = range.x; k < range.x + range.y; ++k) {
        int  j = int(texelFetch(u_items, int(k)).r);
        vec4 A = texelFetch(u_tess_a, j);
        vec3 nj = texelFetch(u_tess_b, j).xyz;
        if (dot(nj, n) <= 0.0) continue;
        vec4  C = texelFetch(u_tess_c, j);
        vec3  v = x - A.xyz;
        float d = dot(v, v);
        // The power of the point to the tessera, what the patches are split by: the squared distance less its weight
        float p = (u_power != 0.0) ? (d < u_power_reach2 ? d - C.z : 1e30) : d;
        if (p < d1) {
            d2 = d1; c2 = c1;
            d1 = p;  c1 = A.xyz; j1 = j;
        } else if (p < d2) {
            d2 = p; c2 = A.xyz;
        }
        if (u_mode != 0) {
            float s2 = A.w * A.w;
            if (d < KERNEL_CUT * KERNEL_CUT * s2) {
                float K = exp(-0.5 * d / s2);
                num += K * C.x;
                den += K * C.y;
            }
        }
    }
    if (j1 < 0) discard;    // no tessera faces this way here: the sphere's own picking and colour stay

    float own = texelFetch(u_tess_b, j1).w;
    float value = (u_mode == 0) ? own : (u_mode == 2 ? num : (den > 0.0 ? num / den : own));
    float ext = u_range.y - u_range.x;
    float param = abs(ext) > 0.0 ? clamp((value - u_range.x) / ext, 0.0, 1.0) : 0.5;
    vec3  color = texture(u_colormap, vec2(param * (255.0 / 256.0) + 0.5 / 256.0, 0.5)).rgb;

    if (u_edges != 0 && d2 < 1e29) {
        // How far the point is from the border of its patch, the plane between the two nearest tesserae where their
        // powers are equal ((p2 - p1) / 2 |c2 - c1| from it), against the size of a pixel there: the point the
        // neighbouring pixel's ray reaches at the same distance
        float dist = (d2 - d1) / (2.0 * max(length(c2 - c1), 1e-6));
        vec3 ro_x, rd_x;
        ray(gl_FragCoord.xy + vec2(1.0, 0.0), ro_x, rd_x);
        float pixel = max(length((ro_x + rd_x * t) - (ro + rd * t)), 1e-6);
        float edge = 1.0 - smoothstep(0.5, 1.5, dist / pixel);
        color *= mix(1.0, 0.35, edge);
    }

    if (u_switched != 0) {
        // Dots on a grid at 45 degrees, their area the share missing below the level, smoothed over a pixel
        float missing = clamp(1.0 - weight / u_level, 0.0, 1.0);
        if (missing > 0.0) {
            vec2  g  = mat2(0.70710678, -0.70710678, 0.70710678, 0.70710678) * gl_FragCoord.xy / u_screen_period;
            float d  = length(fract(g) - 0.5) * u_screen_period;            // pixels from the centre of its dot
            float r  = sqrt(missing / 3.14159265) * u_screen_period;        // pi r^2 = missing * period^2
            float dot_cover = 1.0 - smoothstep(r - 0.5, r + 0.5, d);
            color = mix(color, vec3(0.97), 0.8 * dot_cover);
        }
        // The line where it is the level, a pixel and a half wide
        if (u_outline != 0 && weight_fw > 0.0 && weight_fw < 0.25) {
            float px = abs(weight - u_level) / weight_fw;
            color = mix(color, vec3(1.0), 0.85 * (1.0 - smoothstep(0.5, 1.0, px)));
        }
    }

    out_color = vec4(color, u_unlit != 0 ? 1.0 / 255.0 : 1.0);
    out_picking = encode_index(u_tessera_base + uint(j1));
}
