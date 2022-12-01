#version 150 core

uniform sampler2D u_texture;
uniform float u_exposure = 1.0;
uniform float u_gamma = 2.2;

out vec4 out_frag;

void main() {
    vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    const float exposure_bias = 0.5;
    color.rgb = exposure_bias * u_exposure * pow(color.rgb, 1.0 / vec3(u_gamma));
    out_frag = color;
}