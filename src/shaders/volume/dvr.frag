#version 410 core

// Direct volume rendering: emission-absorption through a 1D transfer function, front to back.
//
// The ray is clipped analytically to the clip box and the opaque depth. Samples are taken at
// SAMPLES_PER_VOXEL along it, each classified on its own (one fetch per sample), with the transfer
// function's alpha corrected for the step length: the alpha of the transfer function is the opacity
// of a slab 1/REF_SAMPLING_RATE of the volume thick. The output is premultiplied colour + coverage,
// in display space (the colours are the colour map's own).

layout(std140) uniform DvrUniforms {
    mat4  u_clip_to_model;
    vec3  u_clip_min;
    float u_tf_min;
    vec3  u_clip_max;
    float u_tf_inv_ext;
    vec2  u_inv_res;
    float u_time;           // 0 for a fixed jitter pattern
    float u_use_depth;
};

uniform sampler3D u_tex_volume;
uniform sampler2D u_tex_depth;
uniform sampler2D u_tex_tf;

layout(location = 0) out vec4 out_color;

const float REF_SAMPLING_RATE = 150.0;
const float SAMPLES_PER_VOXEL = 2.0;
const float T_MIN = 0.005;

vec3 unproject(vec2 ndc_xy, float depth) {
    vec4 p = u_clip_to_model * vec4(ndc_xy, depth * 2.0 - 1.0, 1.0);
    return p.xyz / p.w;
}

float PDnrand(vec2 n) {
    return fract(sin(dot(n, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
    ivec2 px  = ivec2(gl_FragCoord.xy);
    vec2  ndc = (gl_FragCoord.xy * u_inv_res) * 2.0 - 1.0;

    float depth = (u_use_depth > 0.5) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
    vec3 pn  = unproject(ndc, 0.0);
    vec3 pf  = unproject(ndc, depth);
    vec3 seg = pf - pn;

    vec3 inv_seg = 1.0 / mix(seg, vec3(1e-20), equal(seg, vec3(0.0)));
    vec3 t0   = (u_clip_min - pn) * inv_seg;
    vec3 t1   = (u_clip_max - pn) * inv_seg;
    vec3 tmin = min(t0, t1);
    vec3 tmax = max(t0, t1);
    float s0 = max(max(tmin.x, tmin.y), max(tmin.z, 0.0));
    float s1 = min(min(tmax.x, tmax.y), min(tmax.z, 1.0));
    if (s0 >= s1) discard;

    vec3  p0  = pn + seg * s0;
    vec3  ray = seg * (s1 - s0);
    float len = length(ray);
    if (len < 1e-6) discard;

    vec3  dir = ray / len;
    vec3  dim = vec3(textureSize(u_tex_volume, 0));
    int   n   = max(1, int(ceil(len * length(dir * dim) * SAMPLES_PER_VOXEL)));
    vec3  dp  = ray / float(n);
    float w   = (len / float(n)) * REF_SAMPLING_RATE;   // opacity correction exponent, the same for every step
    float jitter = PDnrand(gl_FragCoord.xy + vec2(u_time));

    vec3  L = vec3(0.0);
    float T = 1.0;
    for (int i = 0; i < n; ++i) {
        float d  = texture(u_tex_volume, p0 + dp * (float(i) + jitter)).r;
        float t  = clamp((d - u_tf_min) * u_tf_inv_ext, 0.0, 1.0);
        vec4  tf = texture(u_tex_tf, vec2(t, 0.5));
        float a  = 1.0 - pow(max(1.0 - tf.a, 1e-6), w);
        L += T * a * tf.rgb;
        T *= 1.0 - a;
        if (T < T_MIN) break;
    }

    out_color = vec4(L, 1.0 - T);
}
