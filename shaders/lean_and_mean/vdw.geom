#version 150 core

#ifndef WRITE_DEPTH
#define WRITE_DEPTH 1
#endif

uniform mat4 u_inv_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float z;
} in_vert[];

out GS_FS {
#if WRITE_DEPTH
    flat vec4 view_sphere;
    smooth vec4 view_coord;
#endif
    smooth vec2 uv;
} out_frag;

void emit_vertex(vec2 uv) {
    vec2 axis_a = in_vert[0].axis_a;
    vec2 axis_b = in_vert[0].axis_b;
    vec2 center = in_vert[0].center;
    float z = in_vert[0].z;

    vec2 xy = (center + axis_a * uv.x + axis_b * uv.y);
    vec4 pc = vec4(xy, z, 1);

#if WRITE_DEPTH
    vec4 vc = u_inv_proj_mat * pc;
    out_frag.view_sphere = in_vert[0].view_sphere;
    out_frag.view_coord = vc / vc.w;
#endif
    out_frag.uv = uv;
    gl_Position = pc;
    EmitVertex();
}

void main()
{
    if (in_vert[0].view_sphere.w == 0.0) {
        EndPrimitive();
        return;
    }

    emit_vertex(vec2(-1,-1));
    emit_vertex(vec2( 1,-1));
    emit_vertex(vec2(-1, 1));
    emit_vertex(vec2( 1, 1));

    EndPrimitive();
}