#pragma once

#include <core/array_types.h>
#include <core/gl.h>
#include <mol/spatial_hash.h>
#include <mol/molecule_structure.h>

namespace render {

void initialize();
void shutdown();

struct GPUVolume {
    GLuint texture_id;
    ivec3 resolution;
    vec3 min_box, max_box;
};

void voxelize_scene(Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const uint32> atom_color, ivec3 resolution = {256, 256, 256},
                    vec3 min_box = {0, 0, 0}, vec3 max_box = {0, 0, 0});

void voxelize_scene_gpu(Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const uint32> atom_color, ivec3 resolution = {256, 256, 256},
                        vec3 min_box = {0, 0, 0}, vec3 max_box = {0, 0, 0});

void illuminate_voxels_omnidirectional_constant(const vec3& intensity);

void draw_voxelized_scene(const mat4& view_mat, const mat4& proj_mat);
void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, GLuint color_alpha_tex, GLuint f0_smoothness_tex, const mat4& view_mat,
                      const mat4& proj_mat, float indirect_diffuse_scale, float indirect_specular_scale, float ambient_occlusion_scale);

}  // namespace render
