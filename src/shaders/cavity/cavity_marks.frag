#version 410 core

// The patch of a tessera into the marks of the highlight (cavity.h): unjittered, at the size of the image shown. The
// view ray against the sphere of the tessera and those overlapping it; a point of their union whose nearest tessera
// is this one is of the patch. Marked as the marks of md_gl mark the highlight: b where it is in view (nothing of the
// scene in front of it), a where it is at all, and the group 0 (+8 in view) in the second target.

uniform samplerBuffer  u_tess_a;
uniform samplerBuffer  u_tess_b;
uniform samplerBuffer  u_tess_c;
uniform samplerBuffer  u_spheres;
uniform usamplerBuffer u_cells;
uniform usamplerBuffer u_items;
uniform sampler2D      u_tex_depth;     // the G-buffer's (jittered, at its render scale)

uniform mat4  u_clip_to_cavity;
uniform mat4  u_cavity_to_view;
uniform mat4  u_clip_to_view;
uniform float u_proj_11;
uniform int   u_ortho;
uniform vec2  u_inv_size;
uniform vec3  u_grid_min;
uniform float u_inv_cell;
uniform ivec3 u_grid_dim;
uniform vec4  u_local[32];              // the sphere of the tessera first, then those that overlap it
uniform int   u_num_local;
uniform uint  u_tessera;
uniform float u_power;
uniform float u_power_reach2;

layout(location = 0) out vec4 out_marks;
layout(location = 1) out vec4 out_group;

int nearest_tessera(vec3 x, vec3 n) {
    ivec3 c = ivec3(floor((x - u_grid_min) * u_inv_cell));
    if (any(lessThan(c, ivec3(0))) || any(greaterThanEqual(c, u_grid_dim))) return -1;
    uvec2 range = texelFetch(u_cells, c.x + u_grid_dim.x * (c.y + u_grid_dim.y * c.z)).xy;
    float best = 1e30;
    int   bj = -1;
    for (uint k = range.x; k < range.x + range.y; ++k) {
        int j = int(texelFetch(u_items, int(k)).r);
        if (dot(texelFetch(u_tess_b, j).xyz, n) <= 0.0) continue;
        vec3  v = x - texelFetch(u_tess_a, j).xyz;
        float d = dot(v, v);
        float p = (u_power != 0.0) ? (d < u_power_reach2 ? d - texelFetch(u_tess_c, j).z : 1e30) : d;
        if (p < best) { best = p; bj = j; }
    }
    return bj;
}

void main() {
    vec2 ndc = gl_FragCoord.xy * u_inv_size * 2.0 - 1.0;
    vec4 a = u_clip_to_cavity * vec4(ndc, -1.0, 1.0);
    vec4 b = u_clip_to_cavity * vec4(ndc,  1.0, 1.0);
    vec3 ro = a.xyz / a.w;
    vec3 rd = normalize(b.xyz / b.w - ro);

    // The nearest point along the ray of the union of the local spheres that is of the patch
    float best_t = 1e30;
    vec3  best_n = vec3(0);
    for (int i = 0; i < u_num_local; ++i) {
        vec4  s  = u_local[i];
        vec3  oc = ro - s.xyz;
        float bb = dot(oc, rd);
        float h  = bb * bb - (dot(oc, oc) - s.w * s.w);
        if (h < 0.0) continue;
        float t = -bb - sqrt(h);
        if (t <= 0.0 || t >= best_t) continue;
        vec3 x = ro + rd * t;
        bool buried = false;
        for (int k = 0; k < u_num_local; ++k) {
            if (k == i) continue;
            vec3 v = x - u_local[k].xyz;
            if (dot(v, v) < u_local[k].w * u_local[k].w * 0.9999) { buried = true; break; }
        }
        if (buried) continue;
        vec3 n = normalize(x - s.xyz);
        if (nearest_tessera(x, n) != int(u_tessera)) continue;
        best_t = t;
        best_n = n;
    }
    if (best_t >= 1e29) discard;

    // In view: no closer than the scene within a tolerance, as the marks of md_gl decide it
    vec3  x  = ro + rd * best_t;
    vec3  vc = (u_cavity_to_view * vec4(x, 1.0)).xyz;
    vec3  vn = normalize(mat3(u_cavity_to_view) * best_n);
    vec2  uv = gl_FragCoord.xy * u_inv_size;
    ivec2 size  = textureSize(u_tex_depth, 0);
    ivec2 texel = clamp(ivec2(uv * vec2(size)), ivec2(0), size - 1);
    float depth = texelFetch(u_tex_depth, texel, 0).r;
    vec4  scene = u_clip_to_view * vec4(uv * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
    float scene_z = -scene.z / scene.w;
    float z = -vc.z;
    float pixel  = (u_ortho != 0) ? 2.0 / (u_proj_11 * float(size.y)) : 2.0 * z / (u_proj_11 * float(size.y));
    float facing = (u_ortho != 0) ? abs(vn.z) : abs(dot(vn, normalize(vc)));
    float slope  = sqrt(max(1.0 - facing * facing, 0.0)) / max(facing, 0.125);
    bool  in_view = z <= scene_z + 2.0 * pixel * slope + 2e-3 * z;

    out_marks = vec4(0.0, 0.0, in_view ? 1.0 : 0.0, 1.0);
    out_group = vec4((1.0 + (in_view ? 8.0 : 0.0)) / 255.0, 0.0, 0.0, 0.0);
}
