#pragma once

#include <core/md_vec_math.h>
#include <core/md_str.h>

namespace volume {

void initialize();
void shutdown();

mat4_t compute_model_to_world_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb);
mat4_t compute_world_to_model_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb);
mat4_t compute_texture_to_model_matrix(int dim_x, int dim_y, int dim_z);
mat4_t compute_model_to_texture_matrix(int dim_x, int dim_y, int dim_z);

// Create a transfer function texture with an alpha ramp placed at an origin (0,1) with a specific steepness (scale)
void   compute_transfer_function_texture(uint32_t* texture, int implot_colormap, int resolution = 128, float alpha_origin = 0.0f, float alpha_scale = 1.0f);

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
    - clip_planes:    define a subvolume (min, max)[0-1] which represents the visible portion of the volume
*/

struct RenderDesc {
    struct {
        uint32_t depth = 0;
        uint32_t color = 0;
        uint32_t width = 0;
        uint32_t height = 0;
    } render_target;

    struct {
        uint32_t volume = 0;
        uint32_t transfer_function = 0;
    } texture;

    struct {
        mat4_t model = {};
        mat4_t view = {};
        mat4_t proj = {};
        mat4_t inv_proj = {};
    } matrix;

    struct {
        vec3_t min = {0, 0, 0};
        vec3_t max = {1, 1, 1};
    } clip_volume;

    struct {
        float density = 1.0f;
    } global_scaling;

    struct {
        // Enables temporal jittering of the ray-casting offset
        bool enabled = false;
    } temporal;

    struct {
        size_t count = 0;
        const float* values = NULL;
        const vec4_t* colors = NULL;
    } iso_surface;
    
    bool isosurface_enabled = false;
    bool direct_volume_rendering_enabled = true;

    vec3_t voxel_spacing = {};
};

void render_volume(const RenderDesc& desc);


}  // namespace volume
