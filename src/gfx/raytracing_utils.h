#pragma once

#include <core/array_types.h>
#include <core/gl.h>
#include <mol/spatial_hash.h>
#include <mol/molecule_structure.h>

namespace render {

void initialize();
void shutdown();

void voxelize_scene(Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const uint32> atom_color, ivec3 resolution = {256, 256, 256},
                    vec3 min_box = {0, 0, 0}, vec3 max_box = {0, 0, 0});
void draw_voxelized_scene(const mat4& view_mat, const mat4& proj_mat);
void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, GLuint color_tex, const mat4& view_mat, const mat4& proj_mat);

}  // namespace render
