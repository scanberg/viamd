#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather (far) / near bands -> fill -> composite.
//
// out_frag: half-res colour + signed CoC. The CoC kept is the one of the nearest of the 4 children so foreground
//           silhouettes (sharp or blurred) win over what is behind them (tile classification, gather centre).
//
// The other outputs are the sources of the far gather (mipmapped after this pass) and of the near bands. Each is the
// box average of a per full-res pixel term, so a silhouette moving by a fraction of a pixel (TAA jitter) changes them
// by a fraction:
//   out_near = (rgb, 1) * e    e = m / coc^2: the energy density of a foreground pixel spread over its CoC disk,
//                              m = how much the pixel belongs to the blurred foreground (in front of the focus plane)
//   out_far  = (rgb, 1) * k    k = CoC behind the focus plane: how far the pixel spreads. Weighting by it keeps sharp
//                              pixels from spreading when a coarse mip level averages them with blurred ones.
//   out_coc  = (coc * e, 0, coc * k, 0)   -> mean CoC of what is in a texel, weighted like its colour
//
// Colours are blended in a compressed space, c / (1 + max(c)) (inverted in dof_composite). In linear HDR the backdrop
// (far brighter than lit geometry) dominates every blend it takes part in: a foreground blurred over it would vanish
// into it within a pixel or two of its silhouette instead of fading out over its CoC.
uniform sampler2D u_tex_color;   // full-res HDR
uniform sampler2D u_tex_linear_depth; // full-res linear depth
uniform float u_focus;
uniform float u_aperture_px;     // half-res px
uniform float u_max_coc_px;      // half-res px
layout(location = 0) out vec4 out_frag;
layout(location = 1) out vec4 out_near;
layout(location = 2) out vec4 out_far;
layout(location = 3) out vec4 out_coc;
// Thin-lens circle of confusion, normalised so it is independent of scene scale:
//   coc = aperture * (1 - z_focus / z)        (signed; < 0 in front of the focus plane)
// (1 - zf/z) = (1/zf - 1/z) * zf. Without the zf factor the blur would shrink as 1/zf when zooming out, so the same
// setting would look different for every scene scale.
// u_aperture_px is the CoC (in px of the target being written) of an object at infinity.
float signed_coc_px(float z, float z_focus, float aperture_px, float max_px) {
    return clamp(aperture_px * (1.0 - z_focus / z), -max_px, max_px);
}
vec3 compress(vec3 c) {
    return c / (1.0 + max(max(c.r, c.g), c.b));
}
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy) * 2;
    ivec2 fs = textureSize(u_tex_linear_depth, 0) - 1;
    vec3 c = vec3(0); float zmin = 1e30;
    vec4 near = vec4(0), far = vec4(0), coc = vec4(0);
    for (int i = 0; i < 4; ++i) {
        ivec2 q = min(p + ivec2(i & 1, i >> 1), fs);
        vec3  rgb = compress(max(texelFetch(u_tex_color, q, 0).rgb, vec3(0)));
        float z   = texelFetch(u_tex_linear_depth, q, 0).r;
        float s   = signed_coc_px(z, u_focus, u_aperture_px, u_max_coc_px);
        float m   = smoothstep(0.5, 1.5, -s);
        float e   = m / max(s * s, 1.0);
        float k   = max(s, 0.0);
        c += rgb;
        zmin = min(zmin, z);
        near += vec4(rgb, 1.0) * e;
        far  += vec4(rgb, 1.0) * k;
        coc  += vec4(-s * e, 0.0, s * k, 0.0);
    }
    out_frag = vec4(c * 0.25, signed_coc_px(zmin, u_focus, u_aperture_px, u_max_coc_px));
    out_near = near * 0.25;
    out_far  = far * 0.25;
    out_coc  = coc * 0.25;
}
