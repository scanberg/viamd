#version 410 core

// The glyphs of the tesserae of a cavity (cavity.h), one instance each: a quad facing the camera, as large as the
// silhouette of the sphere that bounds the glyph. A cylinder from the tessera's point outward along its normal, or a
// sphere centred on the point, of the radius and height the representation maps (cavity.cpp: upload_glyphs), never
// thinner than u_min_px pixels.

uniform samplerBuffer u_tess_a;         // xyz: the point of the tessera (the frame of the file)
uniform samplerBuffer u_tess_b;         // xyz: its normal, w: the value it is coloured by
uniform samplerBuffer u_glyph;          // x: the radius (Angstrom), y: the height of a cylinder (Angstrom)
uniform mat4  u_cavity_to_view;
uniform mat4  u_view_to_clip;
uniform int   u_ortho;
uniform int   u_only;                   // the one tessera to draw, or -1 for every instance
uniform int   u_sphere;                 // 1: spheres, 0: cylinders
uniform float u_min_px;
uniform float u_pixel;                  // the size of a pixel in view space: at a depth of 1 in perspective, else as it is

out vec3 v_view;
flat out int   v_tess;
flat out vec3  v_p0;
flat out vec3  v_p1;
flat out float v_r;

void cull() {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    v_view = vec3(0); v_tess = -1; v_p0 = v_p1 = vec3(0); v_r = 0.0;
}

void main() {
    int  j = (u_only >= 0) ? u_only : gl_InstanceID;
    vec4 A = texelFetch(u_tess_a, j);
    vec4 B = texelFetch(u_tess_b, j);
    vec4 G = texelFetch(u_glyph, j);
    float L = (u_sphere != 0) ? 0.0 : max(G.y, 1e-3);
    vec3 p0 = (u_cavity_to_view * vec4(A.xyz, 1.0)).xyz;
    vec3 ax = normalize(mat3(u_cavity_to_view) * B.xyz);
    vec3 p1 = p0 + ax * L;
    vec3 c  = 0.5 * (p0 + p1);
    float px = (u_ortho != 0) ? u_pixel : u_pixel * max(-c.z, 1e-3);
    float r  = max(G.x, u_min_px * px);
    float R  = 0.5 * L + r;

    float d = length(c);
    if (u_ortho == 0 && d <= R * 1.0001) { cull(); return; }
    v_tess = j; v_p0 = p0; v_p1 = p1; v_r = r;

    vec2 corner = vec2(float(gl_VertexID & 1) * 2.0 - 1.0, float((gl_VertexID >> 1) & 1) * 2.0 - 1.0);
    vec3 p;
    if (u_ortho != 0) {
        p = c + vec3(corner * R * 1.01, 0.0);
    } else {
        vec3 w = c / d;
        vec3 u = normalize(abs(w.y) < 0.99 ? cross(w, vec3(0, 1, 0)) : cross(w, vec3(1, 0, 0)));
        vec3 v = cross(u, w);
        float ext = R * d / sqrt(d * d - R * R) * 1.01;
        p = c + (corner.x * u + corner.y * v) * ext;
    }
    v_view = p;
    gl_Position = u_view_to_clip * vec4(p, 1.0);
}
