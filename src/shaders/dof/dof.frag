#version 150 core
#pragma optionNV(unroll all)

// From http://blog.tuxedolabs.com/2018/05/04/bokeh-depth-of-field-in-single-pass.html

uniform sampler2D u_tex_half_res; // Half res color (rgb) and coc (a)
uniform sampler2D u_tex_color;    // Image to be processed 
uniform sampler2D u_tex_depth;    // Linear depth

uniform vec2  u_texel_size;    // The size of a pixel: vec2(1.0/width, 1.0/height) 
uniform float u_focus_depth; 
uniform float u_focus_scale;
uniform float u_time;

const float GOLDEN_ANGLE = 2.39996323; 
const float MAX_BLUR_SIZE = 15.0; 
const float RAD_SCALE = 1.5; // Smaller = nicer blur, larger = faster
const float PI = 3.1415926535;

//#define APPROX

float rand(vec2 n) {
	return fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453);
}

vec4 rand4( vec2 n ) {
    return fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* vec4(43758.5453, 28001.8384, 50849.4141, 12996.89));
}

float srand(vec2 n) {
	return rand(n) * 2.0 - 1.0;
}

vec4 srand4(vec2 n) {
    return rand4(n) * 2.0 - 1.0;
}

float get_signed_coc(float depth, float focus_point, float focus_scale) {
	return clamp((1.0 / focus_point - 1.0 / depth) * focus_scale, -1.0, 1.0);
}

float coc_to_radius(float coc_signed) {
	return abs(coc_signed) * MAX_BLUR_SIZE;
}

vec3 depth_of_field(vec2 tex_coord, float focus_point, float focus_scale) {
	float center_depth  = texture(u_tex_depth, tex_coord).r;
	vec3  center_color  = texture(u_tex_color, tex_coord).rgb;
	float center_coc_signed = get_signed_coc(center_depth, focus_point, focus_scale);
	float center_coc_radius = coc_to_radius(center_coc_signed);

	if (center_coc_radius < 0.5) {
		return center_color;
	}

	// Adaptive radius keeps quality near the focal plane and reduces cost when blur is small.
	float max_radius = clamp(center_coc_radius * 1.35 + 1.0, 2.0, MAX_BLUR_SIZE);

	vec3 near_sum = vec3(0.0);
	vec3 far_sum = vec3(0.0);
	float near_wsum = 0.0;
	float far_wsum = 0.0;
	float near_occ = 0.0;

	if (center_coc_signed < 0.0) {
		near_sum += center_color;
		near_wsum += 1.0;
	} else {
		far_sum += center_color;
		far_wsum += 1.0;
	}
	float radius        = RAD_SCALE;
	float ang           = rand(tex_coord + u_time) * 2.0 * PI;

	for (; radius < max_radius; ang += GOLDEN_ANGLE)
	{
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * u_texel_size * radius;
		float sample_depth = texture(u_tex_depth, tc).r;
		vec3  sample_color = texture(u_tex_color, tc).rgb;
		float sample_coc_signed = get_signed_coc(sample_depth, focus_point, focus_scale);
		float sample_coc_radius = coc_to_radius(sample_coc_signed);
		float coc_weight = smoothstep(radius - 0.5, radius + 0.5, sample_coc_radius);

		if (sample_coc_signed < 0.0) {
			// Foreground layer (near blur) should dominate at depth discontinuities.
			float depth_gate = smoothstep(0.0, 0.004, center_depth - sample_depth);
			float w = coc_weight * mix(0.2, 1.0, depth_gate);
			near_sum += sample_color * w;
			near_wsum += w;
			near_occ += w;
		} else {
			// Background layer contributes mostly behind/at center depth.
			float depth_gate = smoothstep(0.0, 0.004, sample_depth - center_depth);
			float w = coc_weight * mix(0.15, 1.0, depth_gate);
			far_sum += sample_color * w;
			far_wsum += w;
		}

		radius += RAD_SCALE/radius;

#ifdef APPROX
		if ((near_wsum + far_wsum) > 48.0) break;
		//if (radius > MAX_BLUR_SIZE * 0.15) break;
#endif
	}

#ifdef APPROX
const float HALF_RES_RAD_SCALE = 2.0;
	for (; radius < max_radius; ang += GOLDEN_ANGLE) {
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * u_texel_size * radius;
		vec4 sample_color_coc = texture(u_tex_half_res, tc);
		float sample_coc_signed = sample_color_coc.a * 2.0 - 1.0;
		float sample_coc_radius = coc_to_radius(sample_coc_signed);
		
		float sample_depth = textureLod(u_tex_depth, tc, 1).r;
		float coc_weight = smoothstep(radius - 0.5, radius + 0.5, sample_coc_radius);
		if (sample_coc_signed < 0.0) {
			float depth_gate = smoothstep(0.0, 0.004, center_depth - sample_depth);
			float w = coc_weight * mix(0.2, 1.0, depth_gate);
			near_sum += sample_color_coc.rgb * w;
			near_wsum += w;
			near_occ += w;
		} else {
			float depth_gate = smoothstep(0.0, 0.004, sample_depth - center_depth);
			float w = coc_weight * mix(0.15, 1.0, depth_gate);
			far_sum += sample_color_coc.rgb * w;
			far_wsum += w;
		}

		radius += HALF_RES_RAD_SCALE/radius;
	}
#endif

	vec3 near_color = near_wsum > 0.0 ? near_sum / near_wsum : center_color;
	vec3 far_color = far_wsum > 0.0 ? far_sum / far_wsum : center_color;

	// Near layer should be composited on top of far layer.
	float center_near = center_coc_signed < 0.0 ? 1.0 : 0.0;
	float near_presence = near_occ / (near_wsum + far_wsum + 1e-5);
	float near_alpha = clamp(max(center_near, near_presence), 0.0, 1.0);

	return mix(far_color, near_color, near_alpha);
}

in vec2 tc;
out vec4 out_frag;

void main() {
	vec3 dof = depth_of_field(tc, u_focus_depth, u_focus_scale);

	// To hide banding artifacts
	vec4 noise = rand4(tc + u_time) / 20.0;
	dof += noise.xyz;
	out_frag = vec4(dof, 1);
}