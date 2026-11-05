#version 410 core

// Direct volume rendering: emission-absorption through a 1D transfer function, front to back.
//
// The ray is clipped analytically to the clip box, to the half-spaces of up to 8 planes (n . x + d >= 0 kept, in the
// volume's model space) and to the opaque depth. Samples are taken at
// SAMPLES_PER_VOXEL along it, each classified on its own (one fetch per sample), with the transfer
// function's alpha corrected for the step length: the alpha of the transfer function is the opacity
// of a slab 1/REF_SAMPLING_RATE of the volume thick. The output is premultiplied colour + coverage,
// in display space (the colours are the colour map's own).
//
// As slices, where the ray meets each plane within its rectangle and in front of the opaque depth, the transfer
// function's colour of the value there with the plane's opacity, composited front to back; the value from the volume
// (past the clip box too: its values at the faces carry on, the texture clamped to its edges), or from the plane's
// own image where it has one (evaluated on the rectangle, at a resolution of its own). In front of each, its edge: a
// line along the rectangle's border on a dark halo, as wide in pixels wherever it is (the rectangle's size in pixels
// by where the next pixels' rays meet the plane), so that it is behind the planes in front of it as they are.

layout(std140) uniform DvrUniforms {
    mat4  u_clip_to_model;
    vec3  u_clip_min;
    float u_tf_min;
    vec3  u_clip_max;
    float u_tf_inv_ext;
    vec2  u_inv_res;
    float u_time;           // 0 for a fixed jitter pattern
    float u_use_depth;
    vec4  u_planes[8];      // n . x + d >= 0 is kept, model space
    int   u_num_planes;
    float u_linear_out;     // 1: the colours decoded from sRGB to linear, for an HDR target
    float u_slice;          // 1: the volume on the planes (slices), each of its opacity, instead of through their cut
    float u_edge_px;        // half the width of the slices' edges' lines, pixels
    vec4  u_slice_opacity[2];   // slice i's at [i / 4][i % 4]: 0 draws its edge only
    vec4  u_slice_image[2];     // 1: slice i's values from its own image (u_tex_slices[i]), not the volume
    mat4  u_slice_rect_from_model[8];   // to its rectangle's unit square (x, y in [0, 1] on it), its image's too
    vec4  u_slice_edge[8];      // its edge's colour and alpha (0: none)
};

uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_depth;
uniform sampler2D u_tex_tf;
uniform sampler3D u_tex_slices[8];  // W x H x 1: each slice's image, evaluated on its plane

layout(location = 0) out vec4 out_color;

const float REF_SAMPLING_RATE = 150.0;
const float SAMPLES_PER_VOXEL = 2.0;
const float T_MIN = 0.005;

vec3 unproject(vec2 ndc_xy, float depth) {
    vec4 p = u_clip_to_model * vec4(ndc_xy, depth * 2.0 - 1.0, 1.0);
    return p.xyz / p.w;
}

vec3 srgb_to_linear(vec3 c) {
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c));
}

// The value of slice i at p (model space), q on its rectangle: from its image, or the volume. The images by constant
// indices: i is not dynamically uniform where the slices are sorted per pixel, so it is read here, where i is the loop's.
float slice_value(int i, vec3 p, vec2 q) {
    if (u_slice_image[i / 4][i % 4] < 0.5) return texture(u_tex_volume, p).r;
    vec3 t = vec3(q, 0.5);
    switch (i) {
        case 0: return texture(u_tex_slices[0], t).r;
        case 1: return texture(u_tex_slices[1], t).r;
        case 2: return texture(u_tex_slices[2], t).r;
        case 3: return texture(u_tex_slices[3], t).r;
        case 4: return texture(u_tex_slices[4], t).r;
        case 5: return texture(u_tex_slices[5], t).r;
        case 6: return texture(u_tex_slices[6], t).r;
        default: return texture(u_tex_slices[7], t).r;
    }
}

// Where the line o + d s meets slice i's plane, on its rectangle
vec2 slice_rect_hit(int i, vec3 o, vec3 d) {
    float a = dot(u_planes[i].xyz, o) + u_planes[i].w;
    float b = dot(u_planes[i].xyz, d);
    float s = -a / (abs(b) < 1e-20 ? 1e-20 : b);
    return (u_slice_rect_from_model[i] * vec4(o + d * s, 1.0)).xy;
}

float PDnrand(vec2 n) {
    return fract(sin(dot(n, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
    ivec2 px  = ivec2(gl_FragCoord.xy);
    vec2  ndc = (gl_FragCoord.xy * u_inv_res) * 2.0 - 1.0;

    float depth = (u_use_depth > 0.5) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
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
    if (u_slice > 0.5) {
        // the rays of the next pixels across and up: the rectangles' sizes in pixels where this one meets them
        vec2  ndc_x = ndc + vec2(2.0 * u_inv_res.x, 0.0);
        vec2  ndc_y = ndc + vec2(0.0, 2.0 * u_inv_res.y);
        vec3  pnx = unproject(ndc_x, 0.0), segx = unproject(ndc_x, 1.0) - pnx;
        vec3  pny = unproject(ndc_y, 0.0), segy = unproject(ndc_y, 1.0) - pny;
        float w = u_edge_px;
        float hit_s[8];
        int   hit_i[8];
        float hit_v[8];
        float hit_d[8];     // pixels in from the rectangle's border (out: negative)
        int   num_hits = 0;
        for (int i = 0; i < u_num_planes; ++i) {
            float a = dot(u_planes[i].xyz, pn) + u_planes[i].w;
            float b = dot(u_planes[i].xyz, seg);
            if (abs(b) < 1e-20) continue;
            float s = -a / b;
            if (s < 0.0 || s > 1.0) continue;
            vec3  p  = pn + seg * s;
            vec2  q  = (u_slice_rect_from_model[i] * vec4(p, 1.0)).xy;
            vec2  qx = slice_rect_hit(i, pnx, segx) - q;
            vec2  qy = slice_rect_hit(i, pny, segy) - q;
            // pixels per unit of each of q's coordinates, and how far in from each pair of sides
            vec2  px = 1.0 / max(vec2(length(vec2(qx.x, qy.x)), length(vec2(qx.y, qy.y))), vec2(1e-12));
            vec2  in2 = min(q, 1.0 - q) * px;
            float d = (in2.x >= 0.0 && in2.y >= 0.0) ? min(in2.x, in2.y) : -length(max(-in2, vec2(0.0)));
            float opacity = u_slice_opacity[i / 4][i % 4];
            bool  surface = opacity > 0.0 && d > -0.5;
            bool  edge = u_slice_edge[i].a > 0.0 && abs(d) < w + 1.5;
            if (!surface && !edge) continue;
            float v = surface ? slice_value(i, p, q) : 0.0;
            // in order along the ray
            int k = num_hits++;
            while (k > 0 && hit_s[k - 1] > s) {
                hit_s[k] = hit_s[k - 1];
                hit_i[k] = hit_i[k - 1];
                hit_v[k] = hit_v[k - 1];
                hit_d[k] = hit_d[k - 1];
                --k;
            }
            hit_s[k] = s;
            hit_i[k] = i;
            hit_v[k] = v;
            hit_d[k] = d;
        }
        if (num_hits == 0) discard;
        vec3  L = vec3(0.0);
        float T = 1.0;
        for (int h = 0; h < num_hits; ++h) {
            int   i = hit_i[h];
            float d = hit_d[h];
            // its edge, in front of it: the line over its halo (black), each feathered over a pixel
            vec4  e = u_slice_edge[i];
            if (e.a > 0.0) {
                float core = clamp(w + 0.5 - abs(d), 0.0, 1.0) * e.a;
                float halo = clamp(w + 1.5 - abs(d), 0.0, 1.0) * e.a * 0.45;
                vec3  ec = (u_linear_out > 0.5) ? srgb_to_linear(e.rgb) : e.rgb;
                L += T * core * ec;
                T *= 1.0 - (core + halo * (1.0 - core));
            }
            // the plane within its rectangle, its border feathered over a pixel
            float a = u_slice_opacity[i / 4][i % 4] * clamp(d + 0.5, 0.0, 1.0);
            if (a > 0.0) {
                float t = clamp((hit_v[h] - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
                vec3  c = texture(u_tex_tf, vec2(t, 0.5)).rgb;
                c = (u_linear_out > 0.5) ? srgb_to_linear(c) : c;
                L += T * a * c;
                T *= 1.0 - a;
            }
        }
        out_color = vec4(L, 1.0 - T);
        return;
    }
    for (int i = 0; i < u_num_planes; ++i) {
        float a = dot(u_planes[i].xyz, pn) + u_planes[i].w;   // at s = 0
        float b = dot(u_planes[i].xyz, seg);                  // per unit of s
        if (abs(b) < 1e-20) {
            if (a < 0.0) discard;
            continue;
        }
        float s = -a / b;
        if (b > 0.0) s0 = max(s0, s); else s1 = min(s1, s);
    }
    if (s0 >= s1) discard;

    vec3  p0  = pn + seg * s0;
    vec3  ray = seg * (s1 - s0);
    float len = length(ray);
    if (len < 1e-6) discard;

    vec3  dir = ray / len;
    vec3  dim = vec3(textureSize(u_tex_volume, 0));
    int   n   = max(1, int(ceil(len * length(dir * dim) * SAMPLES_PER_VOXEL)));
    vec3  dp  = ray / float(n);
    float w   = (len / float(n)) * REF_SAMPLING_RATE;   // opacity correction exponent, the same for every step
    float jitter = PDnrand(gl_FragCoord.xy + vec2(u_time));

    vec3  L = vec3(0.0);
    float T = 1.0;
    for (int i = 0; i < n; ++i) {
        float d  = texture(u_tex_volume, p0 + dp * (float(i) + jitter)).r;
        float t  = clamp((d - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
        vec4  tf = texture(u_tex_tf, vec2(t, 0.5));
        float a  = 1.0 - pow(max(1.0 - tf.a, 1e-6), w);
        vec3  c  = (u_linear_out > 0.5) ? srgb_to_linear(tf.rgb) : tf.rgb;
        L += T * a * c;
        T *= 1.0 - a;
        if (T < T_MIN) break;
    }

    out_color = vec4(L, 1.0 - T);
}
