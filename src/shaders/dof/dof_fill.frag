#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather -> fill -> composite.
// Pull-push fill of the far layer: where the gather centre belongs to the foreground, what lies behind it is unknown
// (alpha 0 in the gather output). It is what shows through the semi-transparent edge of the foreground blur, so it has
// to continue the far layer around it smoothly. Gathering it from the background around (a quadrature over coarse
// texels that mix foreground and background, normalised by how much background it found) is noisy and blotchy, and
// its edge follows the half-res near / far classification.
//
// Pull: the mip chain of the premultiplied far layer (rgb * a, a) (box filtered).
// Push: from the coarsest level down, each level keeps what it has and takes the rest from the level above:
//   F_l = P_l.rgb + (1 - P_l.a) * up(F_l+1)
// so valid texels at level 0 are left exactly as they are, and a hole is filled from the smallest level that covers it.
uniform sampler2D u_tex_pull;    // premultiplied far layer, mipmapped
uniform sampler2D u_tex_coarse;  // F of the next coarser level (bilinear)
uniform int u_level;             // level of u_tex_pull written by this pass
uniform int u_top;               // 1: coarsest level, nothing above it
out vec4 out_frag;

void main() {
    ivec2 q = ivec2(gl_FragCoord.xy);
    vec4  P = texelFetch(u_tex_pull, q, u_level);
    if (u_top == 1) {
        out_frag = vec4(P.a > 1.0e-6 ? P.rgb / P.a : vec3(0.0), 1.0);
        return;
    }
    // Bilinear taps half a texel (of this level) off centre: a slightly wider tent than plain bilinear upsampling,
    // which keeps the coarse levels from showing as a diamond pattern in large holes
    vec2 inv  = 1.0 / vec2(textureSize(u_tex_pull, u_level));
    vec2 uv   = (vec2(q) + 0.5) * inv;
    vec2 d    = 0.5 * inv;
    vec3 up   = 0.25 * (textureLod(u_tex_coarse, uv + vec2(-d.x, -d.y), 0.0).rgb +
                        textureLod(u_tex_coarse, uv + vec2( d.x, -d.y), 0.0).rgb +
                        textureLod(u_tex_coarse, uv + vec2(-d.x,  d.y), 0.0).rgb +
                        textureLod(u_tex_coarse, uv + vec2( d.x,  d.y), 0.0).rgb);
    out_frag = vec4(P.rgb + (1.0 - clamp(P.a, 0.0, 1.0)) * up, 1.0);
}
