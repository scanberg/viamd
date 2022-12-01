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
	vec3 color;
#if USE_INVERSE
    color = TonemapInvert(Fetch(ivec2(0, 0)));
#else
	#if USE_LOWPASS
	color = 
	    TonemapWithWeight(Fetch(ivec2(-1, -1)), 1.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 0, -1)), 2.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 1, -1)), 1.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2(-1,  0)), 2.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 0,  0)), 4.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 1,  0)), 2.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2(-1,  1)), 1.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 0,  1)), 2.0 / 16.0) +
	    TonemapWithWeight(Fetch(ivec2( 1,  1)), 1.0 / 16.0);
	#else
	color = Tonemap(Fetch(ivec2(0, 0)));
	#endif
#endif
	out_frag = vec4(color, 1);
}