#version 410 core

// A sphere of a VdW Surface (vdw_surface.h), ray cast into the G-buffer as the opaque geometry writes it: the colour,
// the normal, the velocity of the camera's motion and the picking index of its atom. Coloured by its atom, or by the
// points on the surface: the values under their Gaussians about the hit, sum a v g / sum a g, looked up in the grid
// of cells as wide as the widest reaches.
//
// ISOLINES (gfx/isolines.h; the variant built with it defined, drawn while there are lines to show): painted into the
// colour where the value takes given ones, |f - v| / (its change per pixel) within the width. The change per pixel is
// exact: the gradient of f, summed beside it (sum a g (v - f) (-2 w) (x - c) / sum a g), times how far the surface's
// point moves from this pixel's ray to the next one's, along the tangent plane of the sphere (ray differentials, as
// isomesh.frag's). No screen space derivatives, which are wrong where a quad straddles two spheres or the silhouette.

in vec3 v_view;
flat in int  v_sphere;
flat in vec3 v_c;
flat in float v_r;

uniform samplerBuffer  u_sphere_color;  // the atom's colour (RGBA8)
uniform usamplerBuffer u_sphere_atom;   // the atom's index in the system
uniform samplerBuffer  u_point;         // xyz: the point (world), w: its value; sorted by cell
uniform samplerBuffer  u_point_w;       // x: its area, y: its Gaussian's exponent (1/Angstrom^2)
uniform usamplerBuffer u_cell;          // per cell: the first of its points, how many
uniform sampler2D      u_colormap;

uniform mat4  u_view_to_world;
uniform mat4  u_view_to_clip;
uniform mat4  u_curr_view_to_prev_clip;
uniform vec4  u_jitter_uv;
uniform int   u_ortho;
uniform int   u_by_points;
uniform vec2  u_range;
uniform vec3  u_grid_origin;
uniform float u_grid_cell;
uniform ivec3 u_grid_dim;
uniform uint  u_pick_base;

#ifdef ISOLINES
#define ISOLINES_MAX 8
uniform int   u_isoline_count;
uniform float u_isoline_values[ISOLINES_MAX];
uniform int   u_isoline_hover;
uniform float u_isoline_hover_value;
uniform float u_isoline_width;      // pixels
uniform mat4  u_clip_to_view;
uniform vec2  u_inv_res;            // of the viewport
#endif

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_picking;

vec2 encode_normal(vec3 n) {
    float p = sqrt(n.z * 8.0 + 8.0);
    return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

// The points' value at x (world), and whether any reaches it; with ISOLINES also its gradient (world, per Angstrom)
vec2 point_value(vec3 x, out vec3 grad) {
    ivec3 c = clamp(ivec3(floor((x - u_grid_origin) / u_grid_cell)), ivec3(0), u_grid_dim - 1);
    float sw = 0.0, swv = 0.0;
    vec3  sg = vec3(0.0), sgv = vec3(0.0);   // the gradients of the two sums
    for (int z = max(c.z - 1, 0); z <= min(c.z + 1, u_grid_dim.z - 1); ++z)
    for (int y = max(c.y - 1, 0); y <= min(c.y + 1, u_grid_dim.y - 1); ++y)
    for (int k = max(c.x - 1, 0); k <= min(c.x + 1, u_grid_dim.x - 1); ++k) {
        uvec2 cell = texelFetch(u_cell, (z * u_grid_dim.y + y) * u_grid_dim.x + k).xy;
        for (uint o = cell.x; o < cell.x + cell.y; ++o) {
            vec4 P = texelFetch(u_point, int(o));
            vec2 W = texelFetch(u_point_w, int(o)).xy;
            vec3 d = x - P.xyz;
            float g = W.x * exp(-W.y * dot(d, d));
            sw  += g;
            swv += g * P.w;
#ifdef ISOLINES
            vec3 dg = (-2.0 * W.y * g) * d;
            sg  += dg;
            sgv += dg * P.w;
#endif
        }
    }
    grad = vec3(0.0);
    if (!(sw > 0.0)) return vec2(0.0);
    float f = swv / sw;
    grad = (sgv - f * sg) / sw;
    return vec2(f, 1.0);
}

#ifdef ISOLINES
// The view space ray through normalized device coordinates ndc: its origin on the near plane, its direction
void view_ray(vec2 ndc, out vec3 o, out vec3 d) {
    vec4 a = u_clip_to_view * vec4(ndc, -1.0, 1.0);
    vec4 b = u_clip_to_view * vec4(ndc,  0.0, 1.0);
    o = a.xyz / a.w;
    d = b.xyz / b.w - o;
}

// Where the ray o + s d meets the plane through P with normal N; far off at a grazing angle, as the value changes a
// lot from one pixel to the next there
vec3 plane_hit(vec3 o, vec3 d, vec3 P, vec3 N) {
    float dn = dot(N, d);
    float lo = 1.0e-6 * length(d);
    dn = (dn < 0.0) ? min(dn, -lo) : max(dn, lo);
    return o + d * (dot(N, P - o) / dn);
}

// The change of the value per pixel along the surface at P (view space), normal N, G its gradient (view space)
float change_per_pixel(vec3 P, vec3 N, vec3 G) {
    vec2 ndc  = gl_FragCoord.xy * u_inv_res * 2.0 - 1.0;
    vec2 step = 2.0 * u_inv_res;
    vec3 o, d;
    view_ray(ndc, o, d);
    vec3 H  = plane_hit(o, d, P, N);
    view_ray(ndc + vec2(step.x, 0.0), o, d);
    vec3 Hx = plane_hit(o, d, P, N);
    view_ray(ndc + vec2(0.0, step.y), o, d);
    vec3 Hy = plane_hit(o, d, P, N);
    return length(vec2(dot(G, Hx - H), dot(G, Hy - H)));
}

// The coverage of a line w pixels wide at distance d (pixels) from its middle, box filtered over a pixel
float coverage(float d, float w) {
    return clamp(0.5 * w + 0.5 - d, 0.0, 1.0);
}

// As isomesh.frag's apply_isolines: the persistent lines dark with a light rim, the transient one the other way around
vec3 apply_isolines(vec3 c, float f, float g) {
    const vec3  DARK  = vec3(0.02);
    const vec3  LIGHT = vec3(0.95);
    const float RIM   = 2.0;
    float inv_g = 1.0 / max(g, 1.0e-30);
    float w = u_isoline_width;
    float core = 0.0, rim = 0.0;
    for (int i = 0; i < u_isoline_count; ++i) {
        float d = abs(f - u_isoline_values[i]) * inv_g;
        core = max(core, coverage(d, w));
        rim  = max(rim,  coverage(d, w + RIM));
    }
    c = mix(c, LIGHT, rim);
    c = mix(c, DARK, core);
    if (u_isoline_hover != 0) {
        float d = abs(f - u_isoline_hover_value) * inv_g;
        c = mix(c, DARK,  coverage(d, w + RIM));
        c = mix(c, LIGHT, coverage(d, w));
    }
    return c;
}
#endif

void main() {
    if (v_sphere < 0) discard;
    vec3 ro = (u_ortho != 0) ? vec3(v_view.xy, 0.0) : vec3(0.0);
    vec3 rd = (u_ortho != 0) ? vec3(0.0, 0.0, -1.0) : normalize(v_view);
    vec3  oc = ro - v_c;
    float b  = dot(oc, rd);
    float h  = b * b - (dot(oc, oc) - v_r * v_r);
    if (h < 0.0) discard;
    float t  = -b - sqrt(h);
    vec3 x = ro + rd * t;
    vec3 n = (x - v_c) / v_r;

    vec4 clip = u_view_to_clip * vec4(x, 1.0);
    gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;

    vec3 color = texelFetch(u_sphere_color, v_sphere).rgb;
    if (u_by_points != 0) {
        vec3 grad;
        vec2 fv = point_value((u_view_to_world * vec4(x, 1.0)).xyz, grad);
        float ext = u_range.y - u_range.x;
        float param = abs(ext) > 0.0 ? clamp((fv.x - u_range.x) / ext, 0.0, 1.0) : 0.5;
        color = fv.y > 0.0 ? texture(u_colormap, vec2(param * (255.0 / 256.0) + 0.5 / 256.0, 0.5)).rgb : vec3(0.5);
#ifdef ISOLINES
        // the gradient into view space: the chain rule through x_world = M x_view, G_view = M^T G_world
        if (fv.y > 0.0) color = apply_isolines(color, fv.x, change_per_pixel(x, n, grad * mat3(u_view_to_world)));
#endif
    }

    vec4 prev = u_curr_view_to_prev_clip * vec4(x, 1.0);
    vec2 ss_vel = (clip.xy / clip.w - prev.xy / prev.w) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color    = vec4(color, 1.0);
    out_normal   = vec4(encode_normal(n), 0.0, 0.0);
    out_velocity = vec4(ss_vel, 0.0, 0.0);
    out_picking  = encode_index(u_pick_base + texelFetch(u_sphere_atom, v_sphere).r);
}
