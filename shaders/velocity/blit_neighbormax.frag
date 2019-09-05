#version 150 core

#pragma optionNV(unroll all)

#define EXTENT 1

uniform sampler2D u_tex_vel;
uniform vec2 u_tex_vel_texel_size;

in vec2 tc;
out vec4 out_max_neighbor_vel;

void main() {
	vec2 step = u_tex_vel_texel_size;
	vec2 base = tc;

	vec2 mv = vec2(0.0);
	float mv2 = 0.0;

	for (int i = -EXTENT; i <= EXTENT; i++) {
		for (int j = -EXTENT; j <= EXTENT; j++) {
			vec2 v = texture(u_tex_vel, base + vec2(i, j) * step).xy;
			float v2 = dot(v,v);
			if (v2 > mv2) {
				mv = v;
				mv2 = v2;
			}
		}
	}

	out_max_neighbor_vel = vec4(mv, 0.0, 0.0);
}