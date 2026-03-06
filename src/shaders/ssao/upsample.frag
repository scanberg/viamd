#version 330

layout(location = 0) out vec4 outFrag;

uniform sampler2D uDepth; // full-res depth (no mip assumption)
uniform sampler2D uAO;    // half-res AO

void main() {
    ivec2 fullPix = ivec2(gl_FragCoord.xy);
    float depthFull = texelFetch(uDepth, fullPix, 0).r;

    ivec2 halfSize = textureSize(uAO, 0);
    ivec2 fullSize = textureSize(uDepth, 0);

    // Integer mapping: each 2x2 full-res block corresponds to one half-res texel
    ivec2 base = fullPix / 2;

    float bestAO = 1.0;
    float minDepthDiff = 3.402823466e+38;

    for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
            ivec2 hp = clamp(base + ivec2(x, y), ivec2(0), halfSize - ivec2(1));

            // Representative full-res depth for this half-res texel: center of its 2x2 footprint
            ivec2 repFull = hp * 2 + ivec2(1);
            repFull = clamp(repFull, ivec2(0), fullSize - ivec2(1));

            float depthRep = texelFetch(uDepth, repFull, 0).r;
            float aoHalf   = texelFetch(uAO, hp, 0).r;

            float depthDiff = abs(depthRep - depthFull);
            if (depthDiff < minDepthDiff) {
                minDepthDiff = depthDiff;
                bestAO = aoHalf;
            }
        }
    }

    outFrag = vec4(bestAO, bestAO, bestAO, 1.0);
}