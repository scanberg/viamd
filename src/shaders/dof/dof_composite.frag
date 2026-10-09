#version 410 core

// Part of the half-res depth of field: prepass -> tiles -> gather (far) / near bands -> fill -> composite.
// Full-res composite. The far layer replaces the sharp image according to this pixel's own full-res CoC (so blurred
// background does not leak onto sharp edges at half-res granularity); the premultiplied near layer goes on top once.
uniform sampler2D u_tex_color;   // full-res sharp HDR
uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_far;     // half-res far layer, filled behind the foreground (bilinear), see dof_fill
uniform sampler2D u_tex_near;    // half-res near layer, premultiplied and unnormalised (bilinear), see dof_near
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
// The layers hold colours compressed by dof_prepass, the blend happens in that space
vec3 compress(vec3 c) {
    return c / (1.0 + max(max(c.r, c.g), c.b));
}
vec3 decompress(vec3 c) {
    return c / max(1.0 - max(max(c.r, c.g), c.b), 1.0e-4);
}
void main() {
    ivec2 q = ivec2(gl_FragCoord.xy);
    vec3 sharp_hdr = texelFetch(u_tex_color, q, 0).rgb;
    vec3 sharp = compress(max(sharp_hdr, vec3(0)));
    float coc = signed_coc_px(texelFetch(u_tex_linear_depth, q, 0).r, u_focus, u_aperture_px, u_max_coc_px);
    vec2 uv = (vec2(q) + 0.5) * 0.5 / vec2(textureSize(u_tex_far, 0));
    vec3 far  = texture(u_tex_far, uv).rgb;
    vec4 near = texture(u_tex_near, uv);
    // The near field holds the sum of colour * coverage and of coverage over all bands; overlapping foregrounds can
    // cover a pixel more than once
    float near_a = min(near.a, 1.0);
    near = vec4(near.rgb * (near_a / max(near.a, 1.0e-6)), near_a);
    float blur = smoothstep(1.0, 3.0, abs(coc));
    if (blur == 0.0 && near.a == 0.0) {
        out_frag = vec4(sharp_hdr, 1.0);   // exactly the input where nothing is blurred
        return;
    }
    vec3 base = mix(sharp, far, blur);
    out_frag = vec4(decompress(base * (1.0 - near.a) + near.rgb), 1.0);
}
