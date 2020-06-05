#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/types.h>
#include <core/vector_types.h>

namespace postprocessing {

void initialize(int width, int height);
void shutdown();

typedef int Tonemapping;
enum Tonemapping_ { Tonemapping_Passthrough, Tonemapping_ExposureGamma, Tonemapping_Filmic };

struct Descriptor {
    struct {
        vec3 intensity = {20.f, 20.f, 20.f};
    } background;

    struct {
        bool enabled = true;
        float clip_point = 1.0f;
    } bloom;

    struct {
        bool enabled = true;
        Tonemapping mode = Tonemapping_Filmic;
        float exposure = 1.0f;
        float gamma = 2.2f;
    } tonemapping;

    struct {
        bool enabled = true;
        float radius = 6.0f;
        float intensity = 3.0f;
        float bias = 0.1f;
    } ambient_occlusion;

    struct {
        bool enabled = true;
        float focus_depth = 0.5f;
        float focus_scale = 10.f;
    } depth_of_field;

    struct {
        bool enabled = true;
        float feedback_min = 0.88f;
        float feedback_max = 0.97f;
        struct {
            bool enabled = true;
            float motion_scale = 1.f;
        } motion_blur;
    } temporal_reprojection;

    struct {
        GLuint depth = 0;
        GLuint color = 0;
        GLuint normal = 0;
        GLuint velocity = 0;
        GLuint emissive = 0;
        GLuint post_tonemap = 0;
    } input_textures;
};

void shade_and_postprocess(const Descriptor& desc, const ViewParam& view_param);

// void compute_linear_depth(GLuint depth_tex, float near_plane, float far_plane, bool orthographic = false);
void blit_static_velocity(GLuint tex_depth, const ViewParam& view_param);

// void shade_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4& inv_proj_matrix);
// void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer, const vec3& highlight, const vec3& selection, const vec3& outline);

void scale_hsv(GLuint color_tex, vec3 hsv_scale);

// void apply_ssao(GLuint depth_tex, GLuint normal_tex, const mat4& proj_mat, float intensity = 1.5f, float radius = 3.f, float bias = 0.1f);
// void apply_dof(GLuint depth_tex, GLuint color_tex, const mat4& proj_mat, float focus_point, float focus_scale);
// void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure = 1.0f, float gamma = 2.2f);
// void apply_temporal_aa(GLuint depth_tex, GLuint color_tex, GLuint velocity_tex, const ViewParam& param);
void blit_texture(GLuint tex);
void blit_color(vec4 color);

void blur_texture_gaussian(GLuint tex, int num_passes = 1);
void blur_texture_box(GLuint tex, int num_passes = 1);


}  // namespace postprocessing
