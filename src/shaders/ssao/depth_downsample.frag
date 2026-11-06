#version 410 core

// Linear depth mip generation for SSAO.
// Rotated-grid subsampling (McGuire et al. 2012, "Scalable Ambient Obscurance"): every texel keeps the depth of one
// real full-res pixel instead of a min/max/average, so coarse levels never contain depths (and therefore view-space
// positions) that do not exist in the scene. The checkerboard pick spreads the chosen children evenly.
// ssao.frag relies on this exact pick when it reconstructs positions; keep the two in sync.

uniform sampler2D u_tex_linear_depth;
uniform int u_src_lod;
out vec4 out_frag;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 src_max = textureSize(u_tex_linear_depth, u_src_lod) - 1;
    ivec2 c = min(p * 2 + ivec2(p.y & 1, p.x & 1), src_max);
    out_frag = vec4(texelFetch(u_tex_linear_depth, c, u_src_lod).r);
}
