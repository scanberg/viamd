#version 150 core

uniform sampler2D u_tex_depth; // Linear depth
uniform sampler2D u_tex_color;

uniform float u_focus_point;
uniform float u_focus_scale;

float get_signed_coc(float d, float fp, float fs) {
	return clamp((1.0 / fp - 1.0 / d) * fs, -1.0, 1.0);
}

in vec2 tc;
out vec4 out_frag;

void main() {
	float depth = textureLod(u_tex_depth, tc, 1).r;
	vec3  color = texture(u_tex_color, tc).rgb;
	float coc_signed = get_signed_coc(depth, u_focus_point, u_focus_scale);

	// Store signed CoC in [0,1] so half-res path can preserve near/far classification.
	out_frag = vec4(color, coc_signed * 0.5 + 0.5);
}