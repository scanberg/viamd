#version 150 core

// Credits to Brian Karis: http://graphicrants.blogspot.com/2013/12/tone-mapping.html
// And Timothy Lottes: https://gpuopen.com/optimized-reversible-tonemapper-for-resolve/

#ifndef USE_INVERSE
#define USE_INVERSE 0
#endif

#ifndef USE_LOWPASS
#define USE_LOWPASS 0
#endif

uniform sampler2D u_texture;

float max3(float x, float y, float z) { return max(x, max(y, z)); }
float rcp(float x) { return 1.0 / x; }

vec3 Tonemap(vec3 c) { return c * rcp(max3(c.r, c.g, c.b) + 1.0); }
vec3 TonemapWithWeight(vec3 c, float w) { return c * (w * rcp(max3(c.r, c.g, c.b) + 1.0)); }
vec3 TonemapInvert(vec3 c) { return c * rcp(1.0 - max3(c.r, c.g, c.b)); }

vec3 Fetch(ivec2 offset) {
	return texelFetch(u_texture, ivec2(gl_FragCoord.xy) + offset, 0).rgb;
}

out vec4 out_frag;

void main() {
	const float exposure_bias = 0.5;
	const float gamma = 2.4;
	vec3 color;
#if USE_INVERSE
	color = Fetch(ivec2(0, 0));
	color = pow(color, vec3(gamma));
	color = TonemapInvert(color);
    color = color / exposure_bias;
#else
	#if USE_LOWPASS
	color = 
	    TonemapWithWeight(exposure_bias * Fetch(ivec2(-1, -1)), 1.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 0, -1)), 2.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 1, -1)), 1.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2(-1,  0)), 2.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 0,  0)), 4.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 1,  0)), 2.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2(-1,  1)), 1.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 0,  1)), 2.0 / 16.0) +
	    TonemapWithWeight(exposure_bias * Fetch(ivec2( 1,  1)), 1.0 / 16.0);
	#else
	color = Fetch(ivec2(0, 0));
	color = color * exposure_bias;
	color = Tonemap(color);
	#endif
	color = pow(color, 1.0 / vec3(gamma));
#endif
	out_frag = vec4(color, 1);
}