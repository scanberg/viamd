#version 330 core
#extension GL_ARB_shading_language_packing : enable

#ifndef NUM_SUBDIVISIONS
#define NUM_SUBDIVISIONS 16
#endif

#ifndef ORTHONORMALIZE
#define ORTHONORMALIZE 0
#endif

uniform float u_tension = 0.5;

layout(lines_adjacency) in;
layout(points, max_vertices = NUM_SUBDIVISIONS) out;

out vec3 out_control_point;
out uint out_support_vector_xy;
out uint out_support_vector_z_tangent_vector_x;
out uint out_tangent_vector_yz;
out uint out_classification;
out uint out_atom_index;

in Vertex {
    vec3 control_point;
    vec3 support_vector;
    vec4 classification;
    uint atom_index;
} in_vert[];

#ifndef GL_ARB_shading_language_packing
uint packSnorm2x16(in vec2 v) {
    ivec2 iv = ivec2(round(clamp(v, -1.0f, 1.0f) * 32767.0f));
    return uint(iv.y << 16) | uint(iv.x & 0xFFFF);
}

uint packUnorm4x8(in vec4 v) {
    ivec4 iv = ivec4(round(clamp(v, 0.0f, 1.0f) * 255.0f));
    return uint(iv.w << 24) | uint(iv.z << 16) | uint(iv.y << 8) | uint(iv.x);
}
#endif

vec3 catmull_rom(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s, float tension) {
    vec3 v0 = (p2 - p0) * tension;
    vec3 v1 = (p3 - p1) * tension;

    vec3 a = 2.0 * (p1 - p2) + (v0 + v1);
    vec3 b = 3.0 * (p2 - p1) - (2.0 * v0 + v1);

    vec3 res_0 = (a * s * s * s) + (b * s * s);
    vec3 res_1 = (v0 * s) + p1;
    return res_0 + res_1;
}

vec3 catmull_rom_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s, float tension) {
    vec3 v0 = (p2 - p0) * tension;
    vec3 v1 = (p3 - p1) * tension;

    return (2.0 * p1 - 2.0 * p2 + v0 + v1) * 3.0 * s * s + (-3.0 * p1 + 3.0 * p2 - 2.0 * v0 - v1) * 2.0 * s + v0;
}

// Inspired by this:
// https://www.it.uu.se/edu/course/homepage/grafik1/ht07/examples/curves.cpp
float b3(float t) {
    float at = abs(t);
    float t1 = pow(-at + 1.0, 3) * 2.0 / 3.0;
    float t2 = pow(-at + 2.0, 3) / 6.0;
    return mix(t2 - t1, t2, step(1.0, at));
}

float b3_t(float t) {
    float at = abs(t);
    float t1 = -sign(t) * pow(-at + 1.0, 2) * 2.0;
    float t2 = -sign(t) * pow(-at + 2.0, 2) * 0.5;
    return mix(t2 - t1, t2, step(1.0, at));
}

vec3 b_spline(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    return  p0 * b3(s + 1.0) +
            p1 * b3(s) +
            p2 * b3(s - 1.0) +
            p3 * b3(s - 2.0);
}

vec3 b_spline_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    return  p0 * b3_t(s + 1.0) +
            p1 * b3_t(s) +
            p2 * b3_t(s - 1.0) +
            p3 * b3_t(s - 2.0);
}

vec3 spline(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    //return catmull_rom(p0, p1, p2, p3, s, u_tension);
    return b_spline(p0, p1, p2, p3, s);
}

vec3 spline_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    //return catmull_rom_tangent(p0, p1, p2, p3, s, u_tension);
    return b_spline_tangent(p0, p1, p2, p3, s);
}

void main() {
    vec3 cp[4];
    vec3 sv[4];
    vec4 cl[2];
    uint ai[2];

    cp[0] = in_vert[0].control_point;
    cp[1] = in_vert[1].control_point;
    cp[2] = in_vert[2].control_point;
    cp[3] = in_vert[3].control_point;

    sv[0] = in_vert[0].support_vector;
    sv[1] = in_vert[1].support_vector;
    sv[2] = in_vert[2].support_vector;
    sv[3] = in_vert[3].support_vector;

    sv[0] *= sign(dot(sv[0], sv[1]));
    sv[2] *= sign(dot(sv[1], sv[2]));
    sv[3] *= sign(dot(sv[2], sv[3]));

    cl[0] = in_vert[1].classification;
    cl[1] = in_vert[2].classification;

    ai[0] = in_vert[1].atom_index;
    ai[1] = in_vert[2].atom_index;

    for (int i = 0; i < NUM_SUBDIVISIONS; i++) {
        float s = float(i) / float(NUM_SUBDIVISIONS);
        vec3 p = spline(cp[0], cp[1], cp[2], cp[3], s);
        vec3 v = normalize(spline(sv[0], sv[1], sv[2], sv[3], s));
        vec3 t = normalize(spline_tangent(cp[0], cp[1], cp[2], cp[3], s));

#if ORTHONORMALIZE
        v = normalize(v - t*dot(v,t)); // Othonormalize v with respect to t
#endif

        out_control_point = p;
        out_support_vector_xy = packSnorm2x16(v.xy);
        out_support_vector_z_tangent_vector_x = packSnorm2x16(vec2(v.z, t.x));
        out_tangent_vector_yz = packSnorm2x16(t.yz);
        out_classification = packUnorm4x8(mix(cl[0], cl[1], s));
        out_atom_index = s < 0.5 ? ai[0] : ai[1];
        
        EmitVertex();
        EndPrimitive();
    }
}