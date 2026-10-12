#pragma once
#include <core/md_str.h>

namespace analysis_density {
inline constexpr str_t vertex_shader = STR_LIT(R"(
#version 150 core
void main() {
    uint idx = uint(gl_VertexID) % 3U;
    gl_Position = vec4(float(idx & 1U) * 4.0 - 1.0,
                       float((idx >> 1U) & 1U) * 4.0 - 1.0, 0, 1.0);
}
)");

// Shared with the periodic, four-channel Ramachandran renderer. The optional
// triangle mask uses rectangular UVs enclosing the nonperiodic shape simplex.
inline constexpr str_t map_shader = STR_LIT(R"(
#version 330 core
layout(location = 0) out vec4 out_frag0;
layout(location = 1) out vec4 out_frag1;
layout(location = 2) out vec4 out_frag2;
layout(location = 3) out vec4 out_frag3;
uniform sampler2D u_tex_den;
uniform vec4 u_viewport;
uniform vec2 u_inv_res;
uniform bool u_triangle_domain;
uniform vec4 u_map_colors[64];
uniform uint u_map_offset[4];
uniform uint u_map_length[4];
uniform vec2 u_map_range[4];
vec4 map_density(float val, uint map_idx) {
    uint offset = u_map_offset[map_idx];
    uint length = u_map_length[map_idx];
    if (length == 0U) return vec4(0);
    vec2 range = u_map_range[map_idx];
    val = clamp((val - range.x) / max(range.y - range.x, 1e-20), 0, 1);
    float s = val * float(length);
    float t = fract(s);
    uint i0 = clamp(uint(s), 0U, length - 1U);
    uint i1 = clamp(uint(s) + 1U, 0U, length - 1U);
    return mix(u_map_colors[offset + i0], u_map_colors[offset + i1], t);
}
void main() {
    vec2 uv = u_viewport.xy + gl_FragCoord.xy * u_inv_res * u_viewport.zw;
    if (u_triangle_domain && (uv.x < 0.5 * uv.y || uv.x > 1.0 - 0.5 * uv.y)) {
        out_frag0 = out_frag1 = out_frag2 = out_frag3 = vec4(0);
        return;
    }
    vec2 sample_uv = uv;
    if (u_triangle_domain) {
        vec2 size = vec2(textureSize(u_tex_den, 0));
        sample_uv = (0.5 + uv * (size - 1.0)) / size;
    }
    vec4 d = texture(u_tex_den, sample_uv);
    out_frag0 = map_density(d.x, 0U);
    out_frag1 = map_density(d.y, 1U);
    out_frag2 = map_density(d.z, 2U);
    out_frag3 = map_density(d.w, 3U);
}
)");

inline constexpr str_t iso_shader = STR_LIT(R"(
#version 330 core
layout(location = 0) out vec4 out_frag0;
layout(location = 1) out vec4 out_frag1;
layout(location = 2) out vec4 out_frag2;
layout(location = 3) out vec4 out_frag3;
uniform sampler2D u_tex_den;
uniform vec4 u_viewport;
uniform vec2 u_inv_res;
uniform bool u_triangle_domain;
uniform float u_iso_values[32];
uniform vec4 u_iso_level_colors[32];
uniform vec4 u_iso_contour_colors[32];
uniform uint u_iso_offset[4];
uniform uint u_iso_length[4];
uniform float u_iso_contour_line_scale;
vec4 map_density(float val, uint idx) {
    vec4 base = vec4(0);
    vec4 contour = vec4(0);
    uint offset = u_iso_offset[idx];
    uint length = u_iso_length[idx];
    if (length == 0U) return vec4(0);
    uint i = offset;
    for (; i < offset + length - 1U; ++i) {
        if (u_iso_values[i] <= val && val < u_iso_values[i + 1U]) break;
    }
    uint i0 = i;
    uint i1 = min(i + 1U, offset + length - 1U);
    float v1 = u_iso_values[i1];
    vec4 b0 = u_iso_level_colors[i0];
    vec4 b1 = u_iso_level_colors[i1];
    float band = max(fwidth(val) * 2.0 * u_iso_contour_line_scale, 1e-8);
    base = mix(b0, b1, smoothstep(v1 - band, v1 + band, val));
    for (uint j = offset; j < offset + length; ++j) {
        float v = u_iso_values[j];
        contour += u_iso_contour_colors[j] * smoothstep(v - band, v, val) *
                   (1.0 - smoothstep(v, v + band, val));
    }
    return contour + base * (1.0 - contour.a);
}
void main() {
    vec2 uv = u_viewport.xy + gl_FragCoord.xy * u_inv_res * u_viewport.zw;
    if (u_triangle_domain && (uv.x < 0.5 * uv.y || uv.x > 1.0 - 0.5 * uv.y)) {
        out_frag0 = out_frag1 = out_frag2 = out_frag3 = vec4(0);
        return;
    }
    vec2 sample_uv = uv;
    if (u_triangle_domain) {
        vec2 size = vec2(textureSize(u_tex_den, 0));
        sample_uv = (0.5 + uv * (size - 1.0)) / size;
    }
    vec4 d = texture(u_tex_den, sample_uv);
    out_frag0 = map_density(d.x, 0U);
    out_frag1 = map_density(d.y, 1U);
    out_frag2 = map_density(d.z, 2U);
    out_frag3 = map_density(d.w, 3U);
}
)");
}
