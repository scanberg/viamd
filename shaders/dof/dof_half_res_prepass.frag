#version 150 core

uniform sampler2D u_tex_depth; // Linear depth
uniform sampler2D u_tex_color;

uniform float u_focus_point;
uniform float u_focus_scale;

float getBlurSize(float d, float fp, float fs) {
	float coc = clamp((1.0 / fp - 1.0 / d) * fs, -1.0, 1.0);
	return abs(coc);
}

in vec2 tc;
out vec4 out_frag;

void main() {
	float depth = textureLod(u_tex_depth, tc, 1).r;
	vec3  color = textureLod(u_tex_color, tc, 1).rgb;
	float coc   = getBlurSize(depth, u_focus_point, u_focus_scale);

	out_frag = vec4(color, coc);
}