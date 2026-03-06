#version 330 core

#pragma optionNV(unroll all)

#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif

#ifndef AO_NUM_SAMPLES
#define AO_NUM_SAMPLES 16
#endif

struct HBAOData {
    float   radius_to_screen;
    float   neg_inv_r2;
    float   n_dot_v_bias;
    float   z_max;

    vec2    inv_full_res;
    float   ao_multiplier;
    float   pow_exponent;

    vec4    proj_info;

    vec4    sample_pattern[32];
};

layout(std140) uniform u_control_buffer {
    HBAOData control;
};

uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_normal;
uniform sampler2D u_tex_random;

in vec2 tc;
out vec4 out_frag;

vec3 uv_to_view(vec2 uv, float eye_z) {
#if AO_PERSPECTIVE
    return vec3((uv * control.proj_info.xy + control.proj_info.zw) * eye_z, eye_z);
#else
    return vec3((uv * control.proj_info.xy + control.proj_info.zw), eye_z);
#endif
}

vec3 fetch_view_pos(vec2 uv, float lod) {
    float view_depth = textureLod(u_tex_linear_depth, uv, lod).x;
    return uv_to_view(uv, view_depth);
}

vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec3 fetch_view_normal() {
    vec2 enc = texelFetch(u_tex_normal, ivec2(gl_FragCoord.xy), 0).xy;
    vec3 n = decode_normal(enc);
    return n * vec3(1,1,-1);
}

//----------------------------------------------------------------------------------
float falloff(float dist2) {
    // 1 scalar mad instruction
    return dist2 * control.neg_inv_r2 + 1.0;
}

//----------------------------------------------------------------------------------
// P = view-space position at the kernel center
// N = view-space normal at the kernel center
// S = view-space position of the current sample
//----------------------------------------------------------------------------------
float compute_pixel_obscurance(vec3 P, vec3 N, vec3 S) {
    vec3 V = S - P;
    float VdotV = dot(V, V);
    float NdotV = dot(N, V) * inversesqrt(VdotV);

    float falloff_mult = max(0.0, falloff(VdotV));
    return max(0.0, NdotV - control.n_dot_v_bias) * falloff_mult;
}

//----------------------------------------------------------------------------------
vec2 rotate_sample(vec2 sample, vec2 cos_sin) {
    return vec2(sample.x*cos_sin.x - sample.y*cos_sin.y, sample.x*cos_sin.y + sample.y*cos_sin.x);
}

//----------------------------------------------------------------------------------
vec4 get_jitter() {
    // (cos(Alpha),sin(Alpha),rand1,rand2)
    ivec2 coord = ivec2(gl_FragCoord.xy) & (AO_RANDOM_TEX_SIZE - 1);
    vec4 jitter = texelFetch(u_tex_random, coord, 0);

    return jitter;
}

//----------------------------------------------------------------------------------
float compute_ao(vec2 full_res_uv, float radius_pixels, vec4 jitter, vec3 view_position, vec3 view_normal) {
    const float global_mip_offset = -4.3; // -4.3 is recomended in the intel ASSAO implementation
    float mip_offset = log2(radius_pixels * 4) + global_mip_offset;

    float weight_sum = 0.0;
    float ao = 0.0;

    // Create a checkerboard mask offset to alternate the samples for neighboring pixels
    ivec2 coord = ivec2(gl_FragCoord.xy);

    // Use jitter.z to pick a different subset of the 32 pattern per pixel.
    // Map signed [-1,1] -> [0,31]
    int offset = int(floor((jitter.z * 0.5 + 0.5) * 32.0)) & 31;
    int stride = 1;

    float uv_scale = 0.5 + 0.5 * (0.5 + jitter.w * 0.5); // [0.5,1] to reduce the chance of sampling the same depth pixel multiple times for nearby geometry

    for (int i = 0; i < AO_NUM_SAMPLES; i++) {
        vec4 sample = control.sample_pattern[(offset + i * stride) & 31];
        vec2 uv = rotate_sample(sample.xy, jitter.xy) * uv_scale * radius_pixels;
        float weight_scale = sample.z;
        float mip_level = mip_offset + sample.w;
        
        // Skip snapping, it causes artifacts (noisy motion) in large scale scenes
        // Snapping is probably only relevant when rendering in lower resolution, which we don't do
        //vec2 snapped_uv = round(uv) * control.inv_full_res + full_res_uv;
        vec2 snapped_uv = uv * control.inv_full_res + full_res_uv;
        vec3 view_sample = fetch_view_pos(snapped_uv, mip_level);
        ao += compute_pixel_obscurance(view_position, view_normal, view_sample) * weight_scale;
        weight_sum += weight_scale;
    }
    ao *= control.ao_multiplier / weight_sum;

    return clamp(1.0 - ao, 0.0, 1.0);
}

//----------------------------------------------------------------------------------
void main() {
    float view_z = texelFetch(u_tex_linear_depth, ivec2(gl_FragCoord.xy), 0).x;
    if (view_z > control.z_max) discard;

    vec2 uv = tc;
    vec3 view_position = uv_to_view(uv, view_z);
    vec3 view_normal = fetch_view_normal();

  // Compute projection of disk of radius control.R into screen space
#if AO_PERSPECTIVE
    float radius_pixels = control.radius_to_screen / view_position.z;
#else 
    float radius_pixels = control.radius_to_screen;
#endif
    radius_pixels = max(radius_pixels, 3.0); // Avoid sampling the same pixel multiple times for nearby geometry

    // Get jitter vector for the current full-res pixel
    vec4 jitter = get_jitter();
    float ao = compute_ao(uv, radius_pixels, jitter, view_position, view_normal);

    out_frag = vec4(vec3(pow(ao, control.pow_exponent)), 1);
}