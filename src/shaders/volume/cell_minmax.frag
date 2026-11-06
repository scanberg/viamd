#version 410 core

// One texel per cell of the volume (the cube between eight voxel centres): the min and max of its eight
// corners. Cell c, from -1 to dim - 1 on each axis, is texel c + 1; its corners are voxels c and c + 1,
// clamped like CLAMP_TO_EDGE filtering does, so the half voxel border of the volume is covered too. The
// trilinear field anywhere inside a cell lies in this range. Stored in half precision, rounded outwards:
// the stored range always contains the true one. Rendered one layer (z) at a time.

uniform sampler3D u_tex_volume;
uniform int u_layer;

layout(location = 0) out vec2 out_range;

const float HALF_MAX        = 65504.0;
const float HALF_MIN_NORMAL = 6.103515625e-05;     // 2^-14

// Bounds that stay outside [mn, mx] whatever rounding the conversion to half uses: a margin of two half
// ulps, no subnormals (they may be flushed to zero), and overflow only to the infinity on the outside.
vec2 widen_to_half(float mn, float mx) {
    float lo = mn - abs(mn) * (1.0 / 512.0);
    float hi = mx + abs(mx) * (1.0 / 512.0);
    if (abs(lo) < HALF_MIN_NORMAL) lo = -HALF_MIN_NORMAL;
    if (abs(hi) < HALF_MIN_NORMAL) hi =  HALF_MIN_NORMAL;
    float inf = uintBitsToFloat(0x7F800000u);
    lo = (lo < -HALF_MAX) ? -inf : min(lo, HALF_MAX);
    hi = (hi >  HALF_MAX) ?  inf : max(hi, -HALF_MAX);
    return vec2(lo, hi);
}

void main() {
    ivec3 dim_m1 = textureSize(u_tex_volume, 0) - 1;
    ivec3 c = ivec3(ivec2(gl_FragCoord.xy), u_layer) - 1;

    float mn =  3.402823e38;
    float mx = -3.402823e38;
    for (int k = 0; k < 8; ++k) {
        ivec3 o = ivec3(k & 1, (k >> 1) & 1, (k >> 2) & 1);
        float v = texelFetch(u_tex_volume, clamp(c + o, ivec3(0), dim_m1), 0).r;
        mn = min(mn, v);
        mx = max(mx, v);
    }
    out_range = widen_to_half(mn, mx);
}
