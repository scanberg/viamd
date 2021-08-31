#pragma once

#include <core/types.h>
#include <core/string_types.h>

#include <volume.h>
#include <isosurface.h>
#include <gfx/gl.h>

namespace volume {

void initialize();
void shutdown();

//bool init_tf_texture(GLuint* texture, int* width, CStringView path_to_file);
bool init_texture_2D(GLuint* texture, int width, int height, GLenum format);
bool init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format);

bool free_texture(GLuint* texture);

// We assume you set the entire data for the texture
bool set_texture_2D_data(GLuint texture, const void* data, GLenum format);
bool set_texture_3D_data(GLuint texture, const void* data, GLenum format);

mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb);
mat4 compute_world_to_model_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb);
mat4 compute_texture_to_model_matrix(const ivec3& dim);
mat4 compute_model_to_texture_matrix(const ivec3& dim);

void write_volume_to_file(const float* data, int64_t dim_x, int64_t dim_y, int64_t dim_z, CStringView path_to_file);

/*
    Renders a volumetric texture using OpenGL.
    - volume_texture: An OpenGL 3D texture containing the data
    - tf_texture:     An OpenGL 1D texture containing the transfer function
    - depth_texture:  An OpenGL 2D texture containing the depth data in the frame (for stopping ray traversal)
    - model_matrix:   Matrix containing model to world transformation of the volume, which is assumed to occupy a unit cube [0,1] in its model-space
    - view_matrix:    Matrix containing world to view transformation of the camera
    - proj_matrix:    Matrix containing view to clip transformation of the camera
    - density_scale:  global scaling of density
    - alpha_scale:    global alpha scaling of the transfer function
    - isosurface:     information on isovalues and associated colors
    - voxel_spacing:  spacing of voxels in world space
    - clip_volume:    define a subvolume (min, max)[0-1] which represents the visible portion of the volume
*/

struct RenderDesc {
    struct {
        GLuint texture;
        int width;
        int height;
    } render_target;

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
    } clip_volume;

    struct {
        float density = 1.0f;
        float alpha = 1.0f;
    } global_scaling;

    IsoSurfaces isosurface = {};
    bool isosurface_enabled = false;
    bool direct_volume_rendering_enabled = true;
    bool bounding_box_enabled = false;

    vec3 voxel_spacing = {};
};

void render_volume(const RenderDesc& desc);


}  // namespace volume
