#version 150 core

uniform sampler2D u_texture;

out vec4 out_frag;

void main() {
    vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    out_frag = color;
}