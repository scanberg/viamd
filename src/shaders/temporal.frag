#version 150 core

#extension GL_ARB_explicit_attrib_location : enable
#pragma optionNV(unroll all)

// This is ported straight from the implementation given by playdeadgames (MIT License)
// https://github.com/playdeadgames/temporal/blob/master/Assets/Shaders/TemporalReprojection.shader

#ifndef USE_TONEMAP
#define USE_TONEMAP 0
#endif

#ifndef UNJITTER_COLORSAMPLES
#define UNJITTER_COLORSAMPLES 1
#endif
#ifndef UNJITTER_NEIGHBORHOOD
#define UNJITTER_NEIGHBORHOOD 0
#endif
#ifndef UNJITTER_REPROJECTION
#define UNJITTER_REPROJECTION 0
#endif
#ifndef UNJITTER_PREV_SAMPLES
#define UNJITTER_PREV_SAMPLES 0
#endif

#ifndef USE_YCOCG
#define USE_YCOCG 0
#endif
#ifndef USE_CLIPPING
#define USE_CLIPPING 1
#endif
#ifndef USE_DILATION
#define USE_DILATION 1
#endif
#ifndef USE_MOTION_BLUR
#define USE_MOTION_BLUR 1
#endif
#ifndef USE_MOTION_BLUR_NEIGHBORMAX
#define USE_MOTION_BLUR_NEIGHBORMAX 1
#endif
#ifndef USE_OPTIMIZATIONS
#define USE_OPTIMIZATIONS 1
#endif

#ifndef MINMAX_3X3
#define MINMAX_3X3 0
#endif
#ifndef MINMAX_3X3_ROUNDED
#define MINMAX_3X3_ROUNDED 1
#endif
#ifndef MINMAX_4TAP_VARYING
#define MINMAX_4TAP_VARYING 0
#endif
// Tints what the prefilter smooths, red by how much
#ifndef DEBUG_PREFILTER
#define DEBUG_PREFILTER 0
#endif

uniform sampler2D u_tex_main;
uniform sampler2D u_tex_prev;
uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_vel;
uniform sampler2D u_tex_vel_neighbormax;

uniform vec4 u_texel_size;
uniform vec4 u_jitter_uv;

uniform float u_time;
uniform float u_feedback_min = 0.88;
uniform float u_feedback_max = 0.97;
uniform float u_motion_scale = 0.1;

// The prefilter (1: on): the current frame smoothed where its detail is finer than a pixel, before the blend (MINMAX_3X3*)
uniform float u_prefilter = 0.0;
// A pixel's height in view space: x at linear depth 1 (y = 1, perspective) or at any (y = 0, orthographic)
uniform vec2  u_pixel_size = vec2(0.0, 1.0);

float g_prefilter_debug = 0.0;

const float FLT_EPS = 0.00000001f;

// https://software.intel.com/en-us/node/503873
vec3 rgb_to_ycocg(vec3 c)
{
	// Y = R/4 + G/2 + B/4
	// Co = R/2 - B/2
	// Cg = -R/4 + G/2 - B/4
	return vec3(
		 c.x/4.0 + c.y/2.0 + c.z/4.0,
		 c.x/2.0 - c.z/2.0,
		-c.x/4.0 + c.y/2.0 - c.z/4.0
	);
}

// https://software.intel.com/en-us/node/503873
vec3 ycocg_to_rgb(vec3 c)
{
	// R = Y + Co - Cg
	// G = Y + Cg
	// B = Y - Co - Cg
	return vec3(
		c.x + c.y - c.z,
		c.x + c.z,
		c.x - c.y - c.z
	);
}

float luminance(vec3 c) {
	const vec3 w = vec3(0.2125, 0.7154, 0.0721);
	return dot(c, w);
}

float PDnrand( vec2 n ) {
	return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453 );
}
vec2 PDnrand2( vec2 n ) {
	return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* vec2(43758.5453, 28001.8384) );
}
vec3 PDnrand3( vec2 n ) {
	return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* vec3(43758.5453, 28001.8384, 50849.4141 ) );
}
vec4 PDnrand4( vec2 n ) {
	return fract( sin(dot(n.xy, vec2(12.9898, 78.233)))* vec4(43758.5453, 28001.8384, 50849.4141, 12996.89) );
}


float PDsrand( vec2 n ) {
	return PDnrand( n ) * 2.0 - 1.0;
}
vec2 PDsrand2( vec2 n ) {
	return PDnrand2( n ) * 2.0 - 1.0;
}
vec3 PDsrand3( vec2 n ) {
	return PDnrand3( n ) * 2.0 - 1.0;
}
vec4 PDsrand4( vec2 n ) {
	return PDnrand4( n ) * 2.0 - 1.0;
}

// Tonemapper from http://gpuopen.com/optimized-reversible-tonemapper-for-resolve/
float max3(in vec3 value) {
    return max(max(value.x, value.y), value.z);
}

float rcp(float value) {
    return 1.0 / value;
}

vec3 tonemap(in vec3 color) {
    return color * rcp(max3(color) + 1.0);
}

vec3 tonemap(in vec3 color, in float weight) {
    return color * rcp(weight * max3(color) + 1.0);
}

vec3 tonemap_invert(in vec3 color) {
    return color * rcp(1.0 - max3(color));
}

float depth_sample_linear(vec2 uv) {
	return texture(u_tex_linear_depth, uv).x;
}

vec2 sample_velocity(vec2 uv) {
	return texture(u_tex_vel, uv).xy;
}

// The linear depths of the 3x3 around uv, row by row from (-1,-1): what the dilation of the velocity and the prefilter
// both read
void depth_neighbourhood(vec2 uv, out float d[9]) {
	vec2 dd = abs(u_texel_size.xy);
	vec2 du = vec2(dd.x, 0.0);
	vec2 dv = vec2(0.0, dd.y);
	d[0] = depth_sample_linear(uv - dv - du);
	d[1] = depth_sample_linear(uv - dv     );
	d[2] = depth_sample_linear(uv - dv + du);
	d[3] = depth_sample_linear(uv      - du);
	d[4] = depth_sample_linear(uv          );
	d[5] = depth_sample_linear(uv      + du);
	d[6] = depth_sample_linear(uv + dv - du);
	d[7] = depth_sample_linear(uv + dv     );
	d[8] = depth_sample_linear(uv + dv + du);
}

// The nearest of the 3x3: its uv and depth. The first of equals, in the order of d.
vec3 find_closest_fragment_3x3(vec2 uv, float d[9]) {
	int imin = 0;
	for (int i = 1; i < 9; ++i) {
		if (d[i] < d[imin]) imin = i;
	}
	vec2 offset = vec2(float(imin % 3) - 1.0, float(imin / 3) - 1.0);
	return vec3(uv + abs(u_texel_size.xy) * offset, d[imin]);
}

// PREFILTER
// Geometry finer than a pixel (atoms far away, thin lines, one pixel gaps) comes and goes from one frame's samples to
// the next as the jitter moves them: the history is clipped to whatever this frame caught, and it flickers. Where the
// 3x3 shows such detail, the current frame is taken from all nine samples, weighted by where each lies in the scene,
// rather than from one: it varies far less from frame to frame. Elsewhere it is left as it is, so resolved geometry
// keeps its edges. Measured, still view, against the average of the jittered frames: on a protein with atoms of about
// a pixel, the flicker a quarter lower, as much as FXAA before the blend lowered it, with less softening; on a peptide
// with atoms of 20 pixels, a quarter of FXAA's softening there. What it cannot reach: a feature no sample caught this
// frame, as nothing in the 3x3 shows it.

// How far c stands out from both a and b the same way (both above it or both below): a peak or a dip one sample wide.
// 0 on a step, where one of them is on c's side, and on a slope.
float ridge(float a, float b, float c) {
	float da = a - c;
	float db = b - c;
	return (da * db > 0.0) ? min(abs(da), abs(db)) : 0.0;
}

// Both at once, each in its unit (luminance 0.1; depth half a pixel's size): a ridge in luminance alone is often
// resolved, as the dark rim of a shaded cylinder against a light background (on a step in depth), and one in depth
// alone too, as the crease where two atoms meet (on a step in colour). Smoothing either would only soften it.
float ridge_both(int a, int b, int c, float l[9], float d[9], vec2 unit) {
	return min(ridge(l[a], l[b], l[c]) * unit.x, ridge(d[a], d[b], d[c]) * unit.y);
}

// How much the 3x3 holds detail finer than a pixel, 0 to 1: a ridge through the centre in any of four directions, or
// across one of the four samples beside it (a dot next to the centre)
float fine_detail(float l[9], float d[9]) {
	float px = u_pixel_size.x * (u_pixel_size.y > 0.0 ? d[4] : 1.0);
	vec2 unit = vec2(1.0 / 0.1, 1.0 / max(0.5 * px, 1e-8));
	float r = ridge_both(3, 5, 4, l, d, unit);
	r = max(r, ridge_both(1, 7, 4, l, d, unit));
	r = max(r, ridge_both(0, 8, 4, l, d, unit));
	r = max(r, ridge_both(2, 6, 4, l, d, unit));
	r = max(r, ridge_both(0, 2, 1, l, d, unit));
	r = max(r, ridge_both(6, 8, 7, l, d, unit));
	r = max(r, ridge_both(0, 6, 3, l, d, unit));
	r = max(r, ridge_both(2, 8, 5, l, d, unit));
	return smoothstep(0.5, 1.0, r);
}

// The current frame at the pixel's centre from the nine samples, each weighted by its distance from it in the scene:
// texel k holds the scene at k + j, j the jitter in pixels. A Gaussian, exp(-K d^2): K = 1.2 is wider than a fit of
// Blackman-Harris (2.29, Karis 2014); it gave the least flicker for the softness on atoms of a pixel's size.
vec4 reconstruct_3x3(vec4 c[9]) {
	const float K = 1.2;
	const vec3 offs = vec3(-1.0, 0.0, 1.0);
	vec2 j  = u_jitter_uv.xy * u_texel_size.zw;
	vec3 dx = offs + j.x;
	vec3 dy = offs + j.y;
	vec3 wx = exp(-K * dx * dx);
	vec3 wy = exp(-K * dy * dy);
	vec4 row0 = wx.x * c[0] + wx.y * c[1] + wx.z * c[2];
	vec4 row1 = wx.x * c[3] + wx.y * c[4] + wx.z * c[5];
	vec4 row2 = wx.x * c[6] + wx.y * c[7] + wx.z * c[8];
	return (wy.x * row0 + wy.y * row1 + wy.z * row2) / (dot(wx, vec3(1.0)) * dot(wy, vec3(1.0)));
}

vec4 sample_color(sampler2D tex, vec2 uv) {
	vec4 c = texture(tex, uv);
#if USE_YCOCG
	c.rgb = rgb_to_ycocg(c.rgb);
#endif
#if USE_TONEMAP
	c.rgb = tonemap(c.rgb);
#endif
	return c;
}

vec4 sample_color(sampler2D tex, vec2 uv, in float weight) {
	vec4 c = texture(tex, uv);
#if USE_YCOCG
	c.rgb = rgb_to_ycocg(c.rgb);
#endif
#if USE_TONEMAP
	c.rgb = tonemap(c.rgb, weight);
#endif
	return c;
}

vec4 resolve_color(vec4 c) {
#if USE_TONEMAP
	c.rgb = tonemap_invert(c.rgb);
#endif
#if USE_YCOCG
	c.rgb = ycocg_to_rgb(c.rgb);
#endif
	return c;
}

vec4 clip_aabb(vec3 aabb_min, vec3 aabb_max, vec4 p, vec4 q) {
#if USE_OPTIMIZATIONS
	// note: only clips towards aabb center (but fast!)
	vec3 p_clip = 0.5 * (aabb_max + aabb_min);
	vec3 e_clip = 0.5 * (aabb_max - aabb_min) + FLT_EPS;

	vec4 v_clip = q - vec4(p_clip, p.w);
	vec3 v_unit = v_clip.xyz / e_clip;
	vec3 a_unit = abs(v_unit);
	float ma_unit = max(a_unit.x, max(a_unit.y, a_unit.z));

	if (ma_unit > 1.0)
		return vec4(p_clip, p.w) + v_clip / ma_unit;
	else
		return q;// point inside aabb
#else
	vec4 r = q - p;
	vec3 rmax = aabb_max - p.xyz;
	vec3 rmin = aabb_min - p.xyz;

	const float eps = FLT_EPS;

	if (r.x > rmax.x + eps)
		r *= (rmax.x / r.x);
	if (r.y > rmax.y + eps)
		r *= (rmax.y / r.y);
	if (r.z > rmax.z + eps)
		r *= (rmax.z / r.z);

	if (r.x < rmin.x - eps)
		r *= (rmin.x / r.x);
	if (r.y < rmin.y - eps)
		r *= (rmin.y / r.y);
	if (r.z < rmin.z - eps)
		r *= (rmin.z / r.z);

	return p + r;
#endif
}

vec2 sample_velocity_dilated(sampler2D tex, vec2 uv, int support) {
	vec2 du = vec2(u_texel_size.x, 0.0);
	vec2 dv = vec2(0.0, u_texel_size.y);
	vec2 mv = vec2(0.0);
	float rmv = 0.0;

	int end = support + 1;
	for (int i = -support; i != end; i++)
	{
		for (int j = -support; j != end; j++)
		{
			vec2 v = texture(tex, uv + i * dv + j * du).xy;
			float rv = dot(v, v);
			if (rv > rmv)
			{
				mv = v;
				rmv = rv;
			}
		}
	}

	return mv;
}

vec4 sample_color_motion(sampler2D tex, vec2 uv, vec2 ss_vel) {
	const int taps = 5;// on either side!
	vec2 v = 0.5 * ss_vel;

	// Apply some noise offset
	// This breaks up some of the 'trailing shell' artifacts and makes it more noise like
	float srand = PDsrand(uv + vec2(u_time, u_time));
	vec2 vtap = v / taps;
	vec2 pos0 = uv + srand * vtap * 0.5;
	vec4 accu = vec4(0);
	float wsum = 0.0;

	for (int i = -taps; i <= taps; i++) {
        //float w = 1.0;// box
        float w = taps - abs(i) + 1;// triangle
        //float w = 1.0 / (1 + abs(i));// pointy triangle
        accu += w * sample_color(tex, pos0 + i * vtap);
        wsum += w;
	}

	return accu / wsum;
}

// McGuire2012Blur
float cone(float d, float t) {
	return clamp(1.0 - d/t, 0.0, 1.0);
}

float cylinder(float d, float t) {
	return 1.0 - smoothstep(0.95*t, 1.05*t, d);
}

float soft_depth_compare(float z_a, float z_b) {
	const float SOFT_Z_EXTENT = -0.05;
	return clamp(1.0 - (z_a - z_b) / SOFT_Z_EXTENT, 0.0, 1.0);
}

vec2 round_to_nearest(vec2 uv) {
	return ivec2(uv * u_texel_size.zw) * u_texel_size.xy;
}

vec4 reconstruct_color(sampler2D color_tex, vec2 uv0, vec2 ss_vel_max) {
    const int taps = 3;// on either side!

    if (length(ss_vel_max) < FLT_EPS) {
    	return sample_color(color_tex, uv0);
    }

    float srand = PDsrand(uv0 + vec2(u_time, u_time));
    vec2 vtap = 0.5 * ss_vel_max / taps;
    vec2 p0 = uv0 + vtap * (0.5 * srand);
    float d0 = depth_sample_linear(p0);
    float v0 = length(sample_velocity(p0));
    float w = 1.0;
    float wsum = w;
    vec4 accu = w * sample_color(color_tex, p0);

    for (int i = 1; i <= taps; i++)
    {
        vec2 uv = p0 - i * vtap; // Round to nearest?
        float d = depth_sample_linear(uv);
        float v = length(sample_velocity(uv));
        float f = soft_depth_compare(d0, d);
        float b = soft_depth_compare(d, d0);
        float duv = length(uv0 - uv);

        w = f * cone(duv, v) + b * cone(duv, v0) + cylinder(duv, v) * cylinder(duv, v0) * 2.0;
        wsum += w;
        accu += w * sample_color(color_tex, uv);
    }

    for (int i = 1; i <= taps; i++)
    {
        vec2 uv = p0 + i * vtap; // Round to nearest?
        float d = depth_sample_linear(uv);
        float v = length(sample_velocity(uv));
        float f = soft_depth_compare(d0, d);
        float b = soft_depth_compare(d, d0);
        float duv = length(uv0 - uv);

        w = f * cone(duv, v) + b * cone(duv, v0) + cylinder(duv, v) * cylinder(duv, v0) * 2.0;
        wsum += w;
        accu += w * sample_color(color_tex, uv);
    }

    return accu / wsum;
}

vec4 temporal_reprojection(vec2 ss_txc, vec2 ss_vel, float vs_dist, float d[9])
{
	// read texels
#if UNJITTER_COLORSAMPLES
	vec4 texel0 = sample_color(u_tex_main, ss_txc - u_jitter_uv.xy);
#else
	vec4 texel0 = sample_color(u_tex_main, ss_txc);
#endif

#if UNJITTER_PREV_SAMPLES
	vec4 texel1 = sample_color(u_tex_prev, ss_txc - u_jitter_uv.zw - ss_vel);
#else
	vec4 texel1 = sample_color(u_tex_prev, ss_txc - ss_vel);
#endif

	// calc min-max of current neighbourhood
#if UNJITTER_NEIGHBORHOOD
	vec2 uv = ss_txc - u_jitter_uv.xy;
#else
	vec2 uv = ss_txc;
#endif

#if MINMAX_3X3 || MINMAX_3X3_ROUNDED

	vec2 du = vec2(u_texel_size.x, 0.0);
	vec2 dv = vec2(0.0, u_texel_size.y);

	vec4 ctl = sample_color(u_tex_main, uv - dv - du);
	vec4 ctc = sample_color(u_tex_main, uv - dv		);
	vec4 ctr = sample_color(u_tex_main, uv - dv + du);
	vec4 cml = sample_color(u_tex_main, uv - du		);
	vec4 cmc = sample_color(u_tex_main, uv			);
	vec4 cmr = sample_color(u_tex_main, uv + du		);
	vec4 cbl = sample_color(u_tex_main, uv + dv - du);
	vec4 cbc = sample_color(u_tex_main, uv + dv		);
	vec4 cbr = sample_color(u_tex_main, uv + dv + du);

	vec4 cmin = min(ctl, min(ctc, min(ctr, min(cml, min(cmc, min(cmr, min(cbl, min(cbc, cbr))))))));
	vec4 cmax = max(ctl, max(ctc, max(ctr, max(cml, max(cmc, max(cmr, max(cbl, max(cbc, cbr))))))));

	#if MINMAX_3X3_ROUNDED || USE_YCOCG || USE_CLIPPING
		vec4 cavg = (ctl + ctc + ctr + cml + cmc + cmr + cbl + cbc + cbr) / 9.0;
	#endif

	#if MINMAX_3X3_ROUNDED
		vec4 cmin5 = min(ctc, min(cml, min(cmc, min(cmr, cbc))));
		vec4 cmax5 = max(ctc, max(cml, max(cmc, max(cmr, cbc))));
		vec4 cavg5 = (ctc + cml + cmc + cmr + cbc) / 5.0;
		cmin = 0.5 * (cmin + cmin5);
		cmax = 0.5 * (cmax + cmax5);
		cavg = 0.5 * (cavg + cavg5);
	#endif

	if (u_prefilter > 0.0) {
		vec4 c[9] = vec4[9](ctl, ctc, ctr, cml, cmc, cmr, cbl, cbc, cbr);
		float l[9];
		for (int i = 0; i < 9; ++i) l[i] = luminance(c[i].rgb);
		float k = u_prefilter * fine_detail(l, d);
		if (k > 0.0) {
			texel0 = mix(texel0, reconstruct_3x3(c), k);
		}
		g_prefilter_debug = k;
	}

#elif MINMAX_4TAP_VARYING// this is the method used in v2 (PDTemporalReprojection2)

	const float _SubpixelThreshold = 0.5;
	const float _GatherBase = 0.5;
	const float _GatherSubpixelMotion = 0.1666;

	vec2 texel_vel = ss_vel * u_texel_size.xy;
	float texel_vel_mag = length(texel_vel) * vs_dist;
	float k_subpixel_motion = clamp(_SubpixelThreshold / (FLT_EPS + texel_vel_mag), 0.0, 1.0);
	float k_min_max_support = _GatherBase + _GatherSubpixelMotion * k_subpixel_motion;

	vec2 ss_offset01 = k_min_max_support * vec2(-u_texel_size.x, u_texel_size.y);
	vec2 ss_offset11 = k_min_max_support * vec2(u_texel_size.x, u_texel_size.y);
	vec4 c00 = sample_color(u_tex_main, uv - ss_offset11);
	vec4 c10 = sample_color(u_tex_main, uv - ss_offset01);
	vec4 c01 = sample_color(u_tex_main, uv + ss_offset01);
	vec4 c11 = sample_color(u_tex_main, uv + ss_offset11);

	vec4 cmin = min(c00, min(c10, min(c01, c11)));
	vec4 cmax = max(c00, max(c10, max(c01, c11)));

	#if USE_YCOCG || USE_CLIPPING
		vec4 cavg = (c00 + c10 + c01 + c11) / 4.0;
	#endif
#endif

	// shrink chroma min-max
#if USE_YCOCG
	vec2 chroma_extent = vec2(0.25 * 0.5 * (cmax.r - cmin.r));
	vec2 chroma_center = texel0.gb;
	cmin.yz = chroma_center - chroma_extent;
	cmax.yz = chroma_center + chroma_extent;
	cavg.yz = chroma_center;
#endif

	// clamp to neighbourhood of current sample
#if USE_CLIPPING
	texel1 = clip_aabb(cmin.xyz, cmax.xyz, clamp(cavg, cmin, cmax), texel1);
#else
	texel1 = clamp(texel1, cmin, cmax);
#endif

	// feedback weight from unbiased luminance diff (t.lottes)
#if USE_YCOCG
	float lum0 = texel0.r;
	float lum1 = texel1.r;
#else
	float lum0 = luminance(texel0.rgb);
	float lum1 = luminance(texel1.rgb);
#endif
	float unbiased_diff = abs(lum0 - lum1) / max(lum0, max(lum1, 0.2));
	float unbiased_weight = 1.0 - unbiased_diff;
	float unbiased_weight_sqr = unbiased_weight * unbiased_weight;
	float k_feedback = mix(u_feedback_min, u_feedback_max, unbiased_weight_sqr);

	// output
	return mix(texel0, texel1, k_feedback);
}

in vec2 tc;
layout(location = 0) out vec4 out_buff;
layout(location = 1) out vec4 out_frag;

void main() {
	

#if UNJITTER_REPROJECTION
	vec2 uv = tc - u_jitter_uv.xy;
#else
	vec2 uv = tc;
#endif

#if USE_DILATION
	//--- 3x3 norm (sucks)
	//vec2 ss_vel = sample_velocity_dilated(_VelocityBuffer, uv, 1);
	//float vs_dist = depth_sample_linear(uv);

	//--- 5 tap nearest (decent)
	//vec3 c_frag = find_closest_fragment_5tap(uv);
	//vec2 ss_vel = tex2D(_VelocityBuffer, c_frag.xy).xy;
	//float vs_dist = depth_resolve_linear(c_frag.z);

	//--- 3x3 nearest (good)
	float d[9];
	depth_neighbourhood(uv, d);
	vec3 c_frag = find_closest_fragment_3x3(uv, d);
	vec2 ss_vel = texture(u_tex_vel, c_frag.xy).xy;
	float vs_dist = c_frag.z;
#else
	float d[9];
	depth_neighbourhood(uv, d);
	vec2 ss_vel = texture(u_tex_vel, uv).xy;
	float vs_dist = d[4];
#endif

	// temporal resolve
	vec4 color_temporal = temporal_reprojection(uv, ss_vel, vs_dist, d);

	// prepare outputs
	vec4 to_buffer = resolve_color(color_temporal);
	
#if USE_MOTION_BLUR
  #if USE_MOTION_BLUR_NEIGHBORMAX
    ss_vel = u_motion_scale * texture(u_tex_vel_neighbormax, uv).xy;
  #else
    ss_vel = u_motion_scale * ss_vel;
  #endif

	float vel_mag = length(ss_vel * u_texel_size.zw);
	const float vel_trust_full = 2.0;  // 2.0
	const float vel_trust_none = 15.0; // 15.0
	const float vel_trust_span = vel_trust_none - vel_trust_full;
	float trust = 1.0 - clamp(vel_mag - vel_trust_full, 0.0, vel_trust_span) / vel_trust_span;

	#if UNJITTER_COLORSAMPLES
	vec4 color_motion = sample_color_motion(u_tex_main, uv - u_jitter_uv.xy, ss_vel);
	//vec4 color_motion = reconstruct_color(u_tex_main, uv - u_jitter_uv.xy, ss_vel);
	#else
	vec4 color_motion = sample_color_motion(u_tex_main, uv, ss_vel);
	//vec4 color_motion = reconstruct_color(u_tex_main, uv, ss_vel);
	#endif

	vec4 to_screen = resolve_color(mix(color_motion, color_temporal, trust));
#else
	vec4 to_screen = resolve_color(color_temporal);
#endif

	//// NOTE: velocity debug
	//to_screen.g += 1000.0 * length(ss_vel);
	//to_screen = vec4(1000.0 * abs(ss_vel), 0.0, 0.0);
    //to_screen = vec4(ss_vel*1000,0,1);

	out_buff = to_buffer;
	out_frag = to_screen;
#if DEBUG_PREFILTER
	out_frag.rgb = mix(out_frag.rgb, vec3(1.0, 0.0, 0.0), 0.8 * g_prefilter_debug);
#endif

	//out_frag = texture(u_tex_vel_neighbormax, uv);

	//out_frag = texture(u_tex_prev, uv);
}