#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

layout(location = 0) in vec3 in_center;
layout(location = 1) in vec3 in_extent;

out Vertex {
    flat vec3 aabb_min;
    flat vec3 aabb_ext;
    flat uint id;
} out_vert;

void main() {
    out_vert.aabb_min = in_center - in_extent;
    out_vert.aabb_ext = in_extent * 2.0;
    out_vert.id = uint(gl_VertexID);
}