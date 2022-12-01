#version 150 core

#extension GL_ARB_explicit_attrib_location : enable

layout (location = 0) in int in_res_offset;
layout (location = 1) in int in_res_count;

out vec3 out_aabb_center;
out vec3 out_aabb_extent;

uniform samplerBuffer buf_atom_pos;
uniform samplerBuffer buf_atom_rad;

void main() {
	vec3  pos = texelFetch(buf_atom_pos, in_res_offset).xyz;
	float rad = texelFetch(buf_atom_rad, in_res_offset).r; 
	vec3 min_box = pos - rad;
	vec3 max_box = pos + rad;

	for (int i = 1; i < in_res_count; i++) {
		pos = texelFetch(buf_atom_pos, in_res_offset + i).xyz;
		rad = texelFetch(buf_atom_rad, in_res_offset + i).r;
		min_box = min(min_box, pos - rad);
		max_box = max(max_box, pos + rad);
	}

	out_aabb_center = (min_box + max_box) * 0.5;
	out_aabb_extent = (max_box - min_box) * 0.5;
}