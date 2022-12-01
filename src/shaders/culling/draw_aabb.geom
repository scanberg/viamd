#version 150 core

uniform mat4 u_view_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 14) out;

in Vertex {
    flat vec3 aabb_min;
    flat vec3 aabb_ext;
    flat uint id;
} in_vert[];

flat out uint id;

void main()
{
    id = in_vert[0].id;
    vec3 aabb_min = in_vert[0].aabb_min;
    vec3 aabb_ext = in_vert[0].aabb_ext;

    vec4 v[8];
    v[0] = u_view_proj_mat * vec4(aabb_min + vec3(0, 0, 0) * aabb_ext, 1);
    v[1] = u_view_proj_mat * vec4(aabb_min + vec3(1, 0, 0) * aabb_ext, 1);
    v[2] = u_view_proj_mat * vec4(aabb_min + vec3(0, 1, 0) * aabb_ext, 1);
    v[3] = u_view_proj_mat * vec4(aabb_min + vec3(1, 1, 0) * aabb_ext, 1);

    v[4] = u_view_proj_mat * vec4(aabb_min + vec3(0, 0, 1) * aabb_ext, 1);
    v[5] = u_view_proj_mat * vec4(aabb_min + vec3(1, 0, 1) * aabb_ext, 1);
    v[6] = u_view_proj_mat * vec4(aabb_min + vec3(0, 1, 1) * aabb_ext, 1);
    v[7] = u_view_proj_mat * vec4(aabb_min + vec3(1, 1, 1) * aabb_ext, 1);

    gl_Position = v[0]; EmitVertex();
    gl_Position = v[2]; EmitVertex();
    gl_Position = v[1]; EmitVertex();
    gl_Position = v[3]; EmitVertex();
    gl_Position = v[7]; EmitVertex();
    gl_Position = v[2]; EmitVertex();
    gl_Position = v[6]; EmitVertex();
    gl_Position = v[4]; EmitVertex();
    gl_Position = v[7]; EmitVertex();
    gl_Position = v[5]; EmitVertex();
    gl_Position = v[1]; EmitVertex();
    gl_Position = v[4]; EmitVertex();
    gl_Position = v[0]; EmitVertex();
    gl_Position = v[2]; EmitVertex();
    EndPrimitive(); 
}