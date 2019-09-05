#version 150 core

uniform mat4 u_inv_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
    flat vec4 view_velocity;
    flat vec4 color;
    flat vec4 picking_color;
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float z;
    flat uint atom_idx;
} in_vert[];

out GS_FS {
    smooth vec4 view_coord;
    flat vec4 view_sphere;
    flat vec4 view_velocity;
    flat vec4 picking_color;
    flat vec4 color;
    flat uint atom_idx;
} out_frag;

void emit_vertex(vec2 uv) {
    vec2 axis_a = in_vert[0].axis_a;
    vec2 axis_b = in_vert[0].axis_b;
    vec2 center = in_vert[0].center;
    float z = in_vert[0].z;

    vec2 xy = (center + axis_a * uv.x + axis_b * uv.y);
    vec4 pc = vec4(xy, z, 1);
    vec4 vc = u_inv_proj_mat * pc;

    out_frag.view_coord = vc / vc.w;
    gl_Position = pc;
    EmitVertex();
}

void main()
{
    if (in_vert[0].color.a == 0 || in_vert[0].view_sphere.w == 0) {
        EndPrimitive();
        return;
    }

    out_frag.view_sphere = in_vert[0].view_sphere;
    out_frag.view_velocity = in_vert[0].view_velocity;
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.atom_idx = in_vert[0].atom_idx;

    emit_vertex(vec2(-1,-1));
    emit_vertex(vec2( 1,-1));
    emit_vertex(vec2(-1, 1));
    emit_vertex(vec2( 1, 1));

    EndPrimitive();
}