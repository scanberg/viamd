#pragma once

#include <core/types.h>
#include <core/string_types.h>

#include <volume.h>
#include <isosurface.h>
#include <gfx/gl.h>

namespace volume {

void initialize();
void shutdown();

void create_volume_texture(GLuint* texture, const ivec3& dim);
void free_volume_texture(GLuint texture);
void create_tf_texture(GLuint* texture, int* width, CStringView path_to_file);
void set_volume_texture_data(GLuint texture, ivec3 dim, void* data);
mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb);
mat4 compute_texture_to_model_matrix(const ivec3& dim);

void write_to_file(const Volume& volume, CStringView path_to_file);

/*
    Renders a volumetric texture using OpenGL.
    - volume_texture: An OpenGL 3D texture containing the data.
    - tf_texture:     An OpenGL 1D texture containing the transfer function
    - depth_texture:  An OpenGL 2D texture containing the depth data in the frame (for stopping ray traversal).
    - model_matrix:   Matrix containing model to world transformation of the volume, which is assumed to occupy a unit cube [0,1] in its model-space.
    - view_matrix:    Matrix containing world to view transformation of the camera.
    - proj_matrix:    Matrix containing view to clip transformation of the camera.
    - density_scale:  global scaling of density
    - alpha_scale:    global alpha scaling of the transfer function
    - isosurface:     information on isovalues and associated colors
    - voxel_spacing:  spacing of voxels in world space
    - clip_range_min: clip planes min x,y,z [0, 1]
    - clip_range_max: clip planes max x,y,z [0, 1]
*/

struct VolumeRenderDesc {
    struct {
        GLuint volume = 0;
        GLuint transfer_function = 0;
        GLuint depth = 0;
    } texture;

    struct {
        mat4 model = {};
        mat4 view = {};
        mat4 proj = {};
    } matrix;

    struct {
        vec3 min = {0, 0, 0};
        vec3 max = {1, 1, 1};
    } clip_planes;

    struct {
        float density = 1.0f;
        float alpha = 1.0f;
    } global_scaling;

    IsoSurfaces isosurface = {};
    bool isosurface_enabled = false;

    bool direct_volume_rendering_enabled = true;

    vec3 voxel_spacing = {};
};

void render_volume_texture(const VolumeRenderDesc& desc);


}  // namespace volume
