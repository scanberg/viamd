#pragma once

#include "gl.h"

#include <core/array_types.h>
#include <mol/molecule_structure.h>

namespace cone_trace {

struct GPUVolume {
    GLuint texture_id = 0;
	ivec3 resolution = { 0,0,0 };
	vec3 min_box = { 0,0,0 };
	vec3 max_box = { 0,0,0 };
	vec3 voxel_ext = { 0,0,0 };
};


void initialize(int gl_version_major = 3, int gl_version_minor = 3);
void shutdown();

void init_rgba_volume(GPUVolume* vol, ivec3 res, vec3 min_box, vec3 max_box);
void init_occlusion_volume(GPUVolume* vol, vec3 min_box, vec3 max_box, float voxel_ext_target = 4.0f);
void free_volume(GPUVolume* vol);

void compute_occupancy_volume(const GPUVolume& vol, const f32* atom_x, const f32* atom_y, const f32* atom_z, const f32* atom_radius, i32 num_atoms);

void voxelize_spheres_cpu(const GPUVolume& vol, Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const u32> atom_color);
void voxelize_spheres_gpu(const GPUVolume& vol, GLuint position_radius_buffer, GLuint color_buffer, i32 num_spheres);

void illuminate_voxels_omnidirectional_constant(const GPUVolume& vol, const vec3& intensity);

void draw_voxels_scene(const GPUVolume& vol, const mat4& view_mat, const mat4& proj_mat);

void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, GLuint color_alpha_tex, GLuint f0_smoothness_tex, const GPUVolume& vol, const mat4& view_mat,
                      const mat4& proj_mat, float indirect_diffuse_scale, float indirect_specular_scale, float ambient_occlusion_scale);

void render_directional_occlusion(GLuint depth_tex, GLuint normal_tex, const GPUVolume& vol, const mat4& view_mat, const mat4& proj_mat,
                                  float occlusion_scale);

}  // namespace render
