#version 330

layout(location = 0) out vec4 outFrag;

uniform sampler2D uDepth;
uniform sampler2D uAO;

uniform vec2 u_inv_full_res;
uniform vec2 u_inv_half_res;

void main() {
    ivec2 fullPix = ivec2(gl_FragCoord.xy);
    float depthFull = texelFetch(uDepth, fullPix, 0).r;

    ivec2 halfSizeAO = textureSize(uAO, 0);

    // Convert full-res pixel position into half-res texel space (continuous)
    vec2 halfCoord = (gl_FragCoord.xy - vec2(0.5)) * 0.5;
    ivec2 base = ivec2(floor(halfCoord));

    float bestAO = 1.0;
    float minDepthDiff = 3.402823466e+38; // ~FLT_MAX

    for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
            ivec2 hp = clamp(base + ivec2(x, y), ivec2(0), halfSizeAO - ivec2(1));

            float depthHalf = texelFetch(uDepth, hp, 1).r; // mip 1 == half-res
            float aoHalf    = texelFetch(uAO,    hp, 0).r;

            float depthDiff = abs(depthHalf - depthFull);
            if (depthDiff < minDepthDiff) {
                minDepthDiff = depthDiff;
                bestAO = aoHalf;
            }
        }
    }

    outFrag = vec4(vec3(bestAO), 1.0);
}