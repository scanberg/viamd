#version 150 core

uniform sampler2D u_texture_atom_idx;
uniform usamplerBuffer u_buffer_selection;
uniform vec3 u_highlight = vec3(1,1,1);
uniform vec3 u_selection = vec3(1,1,1);
uniform vec3 u_outline = vec3(1,1,0);

in vec2 tc;
out vec4 out_frag;

uint unpackUnorm4x8(in vec4 v) {
    uvec4 uv = uvec4(v * 255.f);
    return uv.x | (uv.y << 8) | (uv.z << 16) | (uv.w << 24);
}

uint fetch_value(ivec2 coord) {
    vec4 c = texelFetch(u_texture_atom_idx, ivec2(gl_FragCoord.xy) + coord, 0);
    uint atom_idx = unpackUnorm4x8(c);
    return texelFetch(u_buffer_selection, int(atom_idx)).x;
}

void main() {
    uint c  = fetch_value(ivec2(0,0));
    uint xn = fetch_value(ivec2(-1,0));
    uint xp = fetch_value(ivec2(+1,0));
    uint yn = fetch_value(ivec2(0,-1));
    uint yp = fetch_value(ivec2(0,+1));

    float f_c = float(c != 0U);
    float f_n = float(xn != 0U) + float(xp != 0U) + float(yn != 0U) + float(yp != 0U);

    float line_t = 1.0 - step(0.0, 4.0 * f_c - f_n);
    float highlight_t = float((c & 1U) != 0U);
    float selection_t = float((c & 2U) != 0U);

    vec3 line_c = u_outline * line_t;
    vec3 fill_c = u_highlight * highlight_t + u_selection * selection_t;

    vec3 color = line_c + fill_c;
    out_frag = vec4(color, 0);
}