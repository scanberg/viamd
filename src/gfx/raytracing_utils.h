#pragma once

#include <core/array_types.h>
#include <core/gl.h>
#include <mol/spatial_hash.h>
#include <mol/molecule_structure.h>

namespace render {

struct GPUVolume {
    GLuint texture_id = 0;
	ivec3 resolution = { 0,0,0 };
	vec3 min_box = { 0,0,0 };
	vec3 max_box = { 0,0,0 };
	vec3 voxel_ext = { 0,0,0 };
};


void initialize(int gl_version_major = 3, int gl_version_minor = 3);
void shutdown();

void init_volume(GPUVolume* vol, ivec3 res, vec3 min_box, vec3 max_box);
void free_volume(GPUVolume* vol);

void voxelize_spheres_cpu(const GPUVolume& vol, Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const uint32> atom_color);

void voxelize_spheres_gpu(const GPUVolume& vol, GLuint position_radius_buffer, GLuint color_buffer, int32 num_spheres);

void illuminate_voxels_omnidirectional_constant(const GPUVolume& vol, const vec3& intensity);

void draw_voxels_scene(const GPUVolume& vol, const mat4& view_mat, const mat4& proj_mat);

void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, GLuint color_alpha_tex, GLuint f0_smoothness_tex, const GPUVolume& vol, const mat4& view_mat,
                      const mat4& proj_mat, float indirect_diffuse_scale, float indirect_specular_scale, float ambient_occlusion_scale);

}  // namespace render
