#version 150 core

#ifndef DILATE_RADIUS
#define DILATE_RADIUS 2
#endif

uniform sampler2D u_tex_vel;
uniform vec2 u_texel_size;

in vec2 tc;
out vec4 out_ss_vel;

void main() {
    vec2 du = vec2(u_texel_size.x, 0.0);
    vec2 dv = vec2(0.0, u_texel_size.y);

    vec2 mv = vec2(0.0);
    float rmv = 0.0;

    for (int y = -DILATE_RADIUS; y <= DILATE_RADIUS; ++y) {
        for (int x = -DILATE_RADIUS; x <= DILATE_RADIUS; ++x) {
            vec2 v = texture(u_tex_vel, tc + float(x) * du + float(y) * dv).xy;
            float rv = dot(v, v);
            if (rv > rmv) {
                mv = v;
                rmv = rv;
            }
        }
    }

    out_ss_vel = vec4(mv, 0.0, 0.0);
}
