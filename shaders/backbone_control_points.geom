#version 330 core
#extension GL_ARB_shading_language_packing : enable

uniform sampler2D u_ramachandran_tex;

layout (triangles_adjacency) in;
layout (points, max_vertices = 1) out;

in uint atom_index[];

out vec3 out_control_point;
out uint out_support_vector_xy;
out uint out_support_vector_z_tangent_vector_x;
out uint out_tangent_vector_yz;
out uint out_classification;
out uint out_atom_index;

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

float dihedral_angle(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3) {
    vec3 b1 = p1 - p0;
    vec3 b2 = p2 - p1;
    vec3 b3 = p3 - p2;
    vec3 c1 = cross(b1, b2);
    vec3 c2 = cross(b2, b3);
    return atan(dot(cross(c1, c2), normalize(b2)), dot(c1, c2));
}

void main() {
    vec3 ca  = gl_in[0].gl_Position.xyz; // Ca[i]
    vec3 c   = gl_in[1].gl_Position.xyz; // C[i]
    vec3 o   = gl_in[2].gl_Position.xyz; // O[i]
    vec3 n   = gl_in[3].gl_Position.xyz; // N[i]
    vec3 c_p = gl_in[4].gl_Position.xyz; // C[i-1]
    vec3 n_n = gl_in[5].gl_Position.xyz; // N[i+1]

    vec3 p = ca;
    vec3 v = normalize(o - c);
    vec3 t = vec3(0);   // Placeholder, tangent is computed analytically in spline shader

    float phi = dihedral_angle(c_p, n, ca, c);
    float psi = dihedral_angle(n, ca, c, n_n);

    // [-PI, PI] -> [0, 1]... With a flip on y axis
    const float ONE_OVER_PI = 1.0 / 3.1515926535;
    vec2 uv = vec2(phi, -psi) * ONE_OVER_PI * 0.5 + 0.5;
    vec4 classification = texture(u_ramachandran_tex, uv);

    out_control_point = p;
    out_support_vector_xy = packSnorm2x16(v.xy);
    out_support_vector_z_tangent_vector_x = packSnorm2x16(vec2(v.z, t.x));
    out_tangent_vector_yz = packSnorm2x16(t.yz);
    out_classification = packUnorm4x8(classification);
    out_atom_index = atom_index[0];

    EmitVertex();
    EndPrimitive();
}