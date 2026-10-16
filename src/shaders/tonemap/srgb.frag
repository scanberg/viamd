#version 150 core

// The display transform of the basic shading mode: linear -> sRGB, clamped to [0,1]. The exact inverse of
// the decode compose_deferred.frag applies to the colours, so a colour in full light shows as it was picked.

uniform sampler2D u_texture;

out vec4 out_frag;

vec3 linear_to_srgb(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(vec3(0.0031308), c));
}

void main() {
    vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    out_frag = vec4(linear_to_srgb(color.rgb), color.a);
}
