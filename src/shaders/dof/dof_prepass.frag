#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Half-res colour + signed CoC. The CoC kept is the one of the nearest of the 4 children so foreground
// silhouettes (sharp or blurred) win over what is behind them.
uniform sampler2D u_tex_color;   // full-res HDR
uniform sampler2D u_tex_linear_depth; // full-res linear depth
uniform float u_focus;
uniform float u_aperture_px;     // half-res px
uniform float u_max_coc_px;      // half-res px
out vec4 out_frag;
// Thin-lens circle of confusion, normalised so it is independent of scene scale:
//   coc = aperture * (1 - z_focus / z)        (signed; < 0 in front of the focus plane)
// (1 - zf/z) = (1/zf - 1/z) * zf. Without the zf factor the blur would shrink as 1/zf when zooming out, so the same
// setting would look different for every scene scale.
// u_aperture_px is the CoC (in px of the target being written) of an object at infinity.
float signed_coc_px(float z, float z_focus, float aperture_px, float max_px) {
    return clamp(aperture_px * (1.0 - z_focus / z), -max_px, max_px);
}
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy) * 2;
    ivec2 fs = textureSize(u_tex_linear_depth, 0) - 1;
    vec3 c = vec3(0); float zmin = 1e30;
    for (int i = 0; i < 4; ++i) {
        ivec2 q = min(p + ivec2(i & 1, i >> 1), fs);
        c += texelFetch(u_tex_color, q, 0).rgb;
        zmin = min(zmin, texelFetch(u_tex_linear_depth, q, 0).r);
    }
    out_frag = vec4(c * 0.25, signed_coc_px(zmin, u_focus, u_aperture_px, u_max_coc_px));
}
