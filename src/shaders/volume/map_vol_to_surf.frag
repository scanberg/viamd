#version 330 core

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_tf;

uniform mat4 u_clip_to_tex_mat;
uniform vec2 u_inv_res;

uniform float u_tf_min;
uniform float u_tf_inv_ext;

layout(location = 0) out vec4  out_color;

float sample_vol(in vec3 tex_coord) {
    return texture(u_tex_volume, tex_coord).r;
}

vec4 classify(in float value) {
    float t = clamp((value - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
    return texture(u_tex_tf, vec2(t, 0.5));
}

vec3 depth_to_tex_coord(vec2 tc, float depth) {
    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
    vec4 tex_coord  = u_clip_to_tex_mat * clip_coord;
    return tex_coord.xyz / tex_coord.w;
}

void main() {
    vec2 uv = gl_FragCoord.xy * u_inv_res;
    float depth = texture(u_tex_depth, uv).x;
    if (depth == 1.0) discard;

    vec3 tex_coord = depth_to_tex_coord(gl_FragCoord.xy * u_inv_res, depth);
    float value    = sample_vol(tex_coord);
    vec4 color     = classify(value);

    out_color = color;
}