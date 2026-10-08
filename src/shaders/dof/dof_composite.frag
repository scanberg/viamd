#version 410 core

// Part of the half-res depth of field: prepass -> tile max -> tile dilate -> gather -> postfilter -> composite.
// Full-res composite: sharp image where in focus, upsampled half-res blur elsewhere. The far-field blend uses the
// full-res pixel's own CoC so blurred background does not leak onto sharp edges at half-res granularity.
uniform sampler2D u_tex_color;   // full-res sharp HDR
uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_dof;     // half-res post-filtered gather result (bilinear)
uniform float u_focus;
uniform float u_aperture_px;     // full-res px
uniform float u_max_coc_px;      // full-res px
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
    ivec2 q = ivec2(gl_FragCoord.xy);
    vec3 sharp = texelFetch(u_tex_color, q, 0).rgb;
    float coc = signed_coc_px(texelFetch(u_tex_linear_depth, q, 0).r, u_focus, u_aperture_px, u_max_coc_px);
    vec4 d = texture(u_tex_dof, (vec2(q) + 0.5) * 0.5 / vec2(textureSize(u_tex_dof, 0)));
    float a = max(d.a, smoothstep(1.0, 3.0, abs(coc)));   // near coverage (half res) or own CoC (full res)
    out_frag = vec4(mix(sharp, d.rgb, a), 1.0);
}
