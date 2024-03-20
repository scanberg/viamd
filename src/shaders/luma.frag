#version 150 core

uniform sampler2D u_tex_rgba;
out vec4 out_rgbl;

void main() {
    vec4 color = texelFetch(u_tex_rgba, ivec2(gl_FragCoord.xy), 0);
    float luma = dot(color.rgb, vec3(0.299, 0.587, 0.114));
    out_rgbl = vec4(color.rgb, luma);
}