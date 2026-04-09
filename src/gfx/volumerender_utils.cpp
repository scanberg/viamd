#include "volumerender_utils.h"

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/postprocessing_utils.h>
#include <color_utils.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <implot.h>

#include <shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

static constexpr str_t v_shader_src_fs_quad = STR_LIT(
    R"(
#version 410 core

out vec2 tc;

uniform vec2 u_tc_scl = vec2(1,1);

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = (gl_Position.xy * 0.5 + 0.5) * u_tc_scl;
}
)");

namespace volume {

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ubo = 0;
    GLuint fbo = 0;
    GLuint ssbo = 0;

    GLuint tex_entry  = 0;
    GLuint tex_exit   = 0;
    GLuint tex_result = 0;

    uint32_t width  = 0;
    uint32_t height = 0;

    struct {
        GLuint entry_exit = 0;
        GLuint entry_exit_depth = 0;
        GLuint dvr_only = 0;
        GLuint iso_only = 0;
        GLuint iso_col_vol = 0;
        GLuint dvr_and_iso = 0;
        GLuint dvr_and_iso_col_vol = 0;
        GLuint splat_color = 0;
    } program;
} gl;

struct UniformData {
    mat4_t view_to_model_mat;
    mat4_t model_to_view_mat;
    mat4_t inv_proj_mat;
    mat4_t model_view_proj_mat;

    vec2_t inv_res;
    float  time;
    float  gamma;

    vec3_t clip_volume_min;
    float  tf_min;
    vec3_t clip_volume_max;
    float  tf_inv_ext;

    vec3_t gradient_spacing_world_space;
    float  exposure;
    mat4_t gradient_spacing_tex_space;

    vec3_t env_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
};

void initialize() {
    GLuint v_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_vert, entryexit_vert_size}, GL_VERTEX_SHADER);
    GLuint f_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER);
    GLuint f_shader_entry_exit_depth    = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define SAMPLE_DEPTH"));

    GLuint v_shader_vol                 = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
    GLuint f_shader_dvr_only            = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR"));
    GLuint f_shader_iso_only            = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_ISO"));
    GLuint f_shader_iso_only_col_vol    = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_ISO\n#define USE_COLOR_VOLUME"));
    GLuint f_shader_dvr_and_iso         = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR\n#define INCLUDE_ISO"));
    GLuint f_shader_dvr_and_iso_col_vol = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR\n#define INCLUDE_ISO\n#define USE_COLOR_VOLUME"));
    GLuint c_shader_splat_color         = gl::compile_shader_from_source({ (const char*)splat_color_comp, splat_color_comp_size }, GL_COMPUTE_SHADER);

    defer {
        glDeleteShader(v_shader_vol);
        glDeleteShader(f_shader_entry_exit);
        glDeleteShader(f_shader_dvr_only);
        glDeleteShader(f_shader_iso_only);
        glDeleteShader(f_shader_iso_only_col_vol);
        glDeleteShader(f_shader_dvr_and_iso);
        glDeleteShader(f_shader_dvr_and_iso_col_vol);
        glDeleteShader(c_shader_splat_color);
    };

    if (v_shader_entry_exit == 0 || v_shader_vol == 0 || f_shader_entry_exit == 0|| f_shader_dvr_only == 0 || f_shader_iso_only == 0 || f_shader_iso_only_col_vol == 0 || f_shader_dvr_and_iso == 0 || f_shader_dvr_and_iso_col_vol == 0) {
        MD_LOG_ERROR("shader compilation failed, shader program for raycasting will not be updated");
        return;
    }

    if (c_shader_splat_color == 0) {
        MD_LOG_ERROR("shader compilation failed, shader program for splat color computation will not be updated");
        return;
    }

    if (!gl.program.entry_exit) gl.program.entry_exit = glCreateProgram();
    if (!gl.program.entry_exit_depth) gl.program.entry_exit_depth = glCreateProgram();
    if (!gl.program.dvr_only) gl.program.dvr_only = glCreateProgram();
    if (!gl.program.iso_only) gl.program.iso_only = glCreateProgram();
    if (!gl.program.iso_col_vol) gl.program.iso_col_vol = glCreateProgram();
    if (!gl.program.dvr_and_iso) gl.program.dvr_and_iso = glCreateProgram();
    if (!gl.program.dvr_and_iso_col_vol) gl.program.dvr_and_iso_col_vol = glCreateProgram();

    if (!gl.program.splat_color) gl.program.splat_color = glCreateProgram();

    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit};
        gl::attach_link_detach(gl.program.entry_exit, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit_depth};
        gl::attach_link_detach(gl.program.entry_exit_depth, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_only};
        gl::attach_link_detach(gl.program.dvr_only, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_iso_only};
        gl::attach_link_detach(gl.program.iso_only, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_iso_only_col_vol};
        gl::attach_link_detach(gl.program.iso_col_vol, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_and_iso};
        gl::attach_link_detach(gl.program.dvr_and_iso, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_and_iso_col_vol};
        gl::attach_link_detach(gl.program.dvr_and_iso_col_vol, shaders, (int)ARRAY_SIZE(shaders));
    }

    gl::attach_link_detach(gl.program.splat_color, &c_shader_splat_color, 1);

    if (!gl.vbo) {
        // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
        constexpr uint8_t cube_strip[42] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                                            0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
        glGenBuffers(1, &gl.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube_strip), cube_strip, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!gl.vao) {
        glGenVertexArrays(1, &gl.vao);
        glBindVertexArray(gl.vao);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_UNSIGNED_BYTE, GL_FALSE, 0, (const GLvoid*)0);
        glBindVertexArray(0);
    }

    if (!gl.ubo) {
        glGenBuffers(1, &gl.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformData), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }

    if (!gl.tex_entry) {
        glGenTextures(1, &gl.tex_entry);
    }

    if (!gl.tex_exit) {
        glGenTextures(1, &gl.tex_exit);
    }

    if (!gl.tex_result) {
        glGenTextures(1, &gl.tex_result);
    }

    if (!gl.fbo) {
        glGenFramebuffers(1, &gl.fbo);
    }

    if (!gl.ssbo) {
        glGenBuffers(1, &gl.ssbo);
    }
}

void shutdown() {}

mat4_t compute_model_to_world_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            ext.x, 0, 0, 0,
            0, ext.y, 0, 0,
            0, 0, ext.z, 0,
            off.x, off.y, off.z, 1,
        };
}

mat4_t compute_world_to_model_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            1.0f / ext.x, 0, 0, 0,
            0, 1.0f / ext.y, 0, 0,
            0, 0, 1.0f / ext.z, 0,
            -off.x, -off.y, -off.z, 1,
        };
}

mat4_t compute_model_to_texture_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

mat4_t compute_texture_to_model_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

void compute_transfer_function_texture_simple(uint32_t* tex, int colormap, float alpha_scale, int res) {
    ASSERT(tex);
    if (res <= 0) {
        MD_LOG_ERROR("Bad input resolution");
        return;
    }

    size_t temp_pos = md_temp_get_pos();
    size_t bytes = sizeof(uint32_t) * res;
    uint32_t* pixel_data = (uint32_t*)md_temp_push(bytes);

    // Update colormap texture
    for (int i = 0; i < res; ++i) {
        float t = CLAMP((float)i / (float)(res - 1), 0.0f, 1.0f);
        ImVec4 col = ImPlot::SampleColormap(t, colormap);

        col.w = MIN(160 * t*t, 0.341176f) * alpha_scale;
        col.w = CLAMP(col.w, 0.0f, 1.0f);

        pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
    }

    gl::init_texture_2D(tex, res, 1, GL_RGBA8);
    gl::set_texture_2D_data(*tex, 0, pixel_data, GL_RGBA8);

    md_temp_set_pos_back(temp_pos);
}

void compute_transfer_function_texture(uint32_t* tex, int colormap, ramp_type_t type, float ramp_scale, float ramp_period, int res) {
    ASSERT(tex);
    if (res <= 0) {
        MD_LOG_ERROR("Bad input resolution");
        return;
    }

    size_t temp_pos = md_temp_get_pos();
    size_t bytes = sizeof(uint32_t) * res;
    uint32_t* pixel_data = (uint32_t*)md_temp_push(bytes);

    const float s = ramp_scale;
    const float p = ramp_period;

    // Update colormap texture
    for (int i = 0; i < res; ++i) {
        float t = (float)i / (float)(res - 1);
        ImVec4 col = ImPlot::SampleColormap(t, colormap);

        // Remap linear t which has a profile of '/' to '\/' to put the alpha ramp in the center (0.5)
        // y(x) = abs((1.0 - c) * (c - x)) Paste this into a grapher and you will see!alpha_origin;
        switch (type) {
        case RAMP_TYPE_SAWTOOTH:
            t = s * fmodf(t, p);
            break;
        case RAMP_TYPE_TRIANGLE:
            t = s * (2.0f/p) * fabsf(fmodf(t-p*0.5f,p)-p*0.5f);
            break;
        default:
            break;
        }

        t = CLAMP(t, 0.0f, 1.0f);
        col.w = t;

        pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
    }

    gl::init_texture_2D(tex, res, 1, GL_RGBA8);
    gl::set_texture_2D_data(*tex, 0, pixel_data, GL_RGBA8);

    md_temp_set_pos_back(temp_pos);
}

static void splat_point_color_volume_GPU(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    PUSH_GPU_SECTION("SPLAT COLOR VOLUME")
    glUseProgram(gl.program.splat_color);

    glUniformMatrix4fv(glGetUniformLocation(gl.program.splat_color, "u_world_to_model"), 1, GL_FALSE, (const float*)world_to_model);
    glUniformMatrix4fv(glGetUniformLocation(gl.program.splat_color, "u_voxel_to_world"), 1, GL_FALSE, (const float*)index_to_world);
    glUniform3iv(glGetUniformLocation(gl.program.splat_color, "u_volume_dim"),    1, volume_dim);
    glUniform3fv(glGetUniformLocation(gl.program.splat_color, "u_voxel_spacing"), 1, voxel_spacing);
    glUniform1i (glGetUniformLocation(gl.program.splat_color, "u_num_points"), (int)point_count);
    glUniform1f(glGetUniformLocation(gl.program.splat_color, "u_power"), power);

    size_t point_xyzw_offset    = 0;
	size_t point_xyzw_size      = ALIGN_TO(sizeof(vec4_t) * point_count, 256);
	size_t point_color_offset   = point_xyzw_offset + point_xyzw_size;
    size_t point_color_size     = ALIGN_TO(sizeof(uint32_t) * point_count, 256);
	size_t total_buffer_size    = ALIGN_TO(point_xyzw_size + point_color_size, 256);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, gl.ssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, total_buffer_size , NULL, GL_DYNAMIC_DRAW);

    glBufferSubData(GL_SHADER_STORAGE_BUFFER, point_xyzw_offset,  sizeof(vec4_t)   * point_count, point_xyzw);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, point_color_offset, sizeof(uint32_t) * point_count, point_color);

    // Bind two contiguous vec4 arrays as separate shader storage buffer bindings for position/radius and color
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, gl.ssbo, point_xyzw_offset,  point_xyzw_size);  // point_xyzw
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, gl.ssbo, point_color_offset, point_color_size); // point_color
    glBindImageTexture(0, vol_texture, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA8);

    const int num_wg[3] = {
        DIV_UP((int)volume_dim[0], 8),
        DIV_UP((int)volume_dim[1], 8),
        DIV_UP((int)volume_dim[2], 8),
    };

    glDispatchCompute(num_wg[0], num_wg[1], num_wg[2]);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    POP_GPU_SECTION()
}

static void splat_point_color_volume_CPU(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    ASSERT(volume_dim[0] % 4 == 0);
    ASSERT(volume_dim[1] % 4 == 0);
    ASSERT(volume_dim[2] % 4 == 0);

    (void)voxel_spacing;
    (void)world_to_model;

    md_allocator_i* temp_alloc = md_get_heap_allocator();

    size_t bytes = sizeof(vec4_t) * volume_dim[0] * volume_dim[1] * volume_dim[2];
    vec4_t* result = (vec4_t*)md_alloc(temp_alloc, bytes);
    ASSERT(result);

    mat4_t index_to_world_mat = mat4_load((const float*)index_to_world);

    for (int z = 0; z < volume_dim[2]; ++z) {
        for (int y = 0; y < volume_dim[1]; ++y) {
            for (int x = 0; x < volume_dim[0]; ++x) {
                vec4_t voxel_pos_index = {(float)x, (float)y, (float)z, 1.0f};
                vec4_t voxel_pos_world = mat4_mul_vec4(index_to_world_mat, voxel_pos_index);
                vec4_t acc = {0};
                float  sum = 0.0f;

                for (size_t i = 0; i < point_count; ++i) {
                    vec4_t point_world = point_xyzw[i];
                    float sigma = point_world.w;
                    point_world.w = 1.0f; // Ignore radius for distance calculation, we will use it as part of the influence factor instead
                    float dist2 = vec4_distance_squared(voxel_pos_world, point_world);

                    float inv2sig2  = 1.0 / (2.0 * sigma * sigma);
                    float w = exp(-dist2 * inv2sig2 * power);
                    acc += w * convert_color(point_color[i]);
                    sum += w;
                }

                vec4_t color = (sum > 0.0) ? acc / sum : vec4_set(1.0, 1.0, 1.0, 0.0);
                acc /= sum;
                int linear_idx = x + y * volume_dim[0] + z * volume_dim[0] * volume_dim[1];
                result[linear_idx] = color;
            }
        }
    }

    // Write result into volume texture
    glBindTexture(GL_TEXTURE_3D, vol_texture);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, volume_dim[0], volume_dim[1], volume_dim[2], GL_RGBA, GL_FLOAT, result);
    glBindTexture(GL_TEXTURE_3D, 0);

    md_free(temp_alloc, result, bytes);
}

// vol_origin is the origin of the volume in world space
// voxel_spacing is the spacing of voxels in world space
void compute_point_color_volume(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, double power) {
    if (!glIsTexture(vol_texture)) {
        MD_LOG_ERROR("Invalid volume texture");
        return;
    }

    int gl_major, gl_minor;
    glGetIntegerv(GL_MAJOR_VERSION, &gl_major);
    glGetIntegerv(GL_MINOR_VERSION, &gl_minor);

    if (gl_major > 4 || (gl_major == 4 && gl_minor >= 3)) {
        splat_point_color_volume_GPU(vol_texture, volume_dim, voxel_spacing, world_to_model, index_to_world, point_xyzw, point_color, point_count, (float)power);
    } else {
        splat_point_color_volume_CPU(vol_texture, volume_dim, voxel_spacing, world_to_model, index_to_world, point_xyzw, point_color, point_count, (float)power);
    }
}

void render_volume(const RenderDesc& desc) {
    if (!desc.dvr.enabled && !desc.iso.enabled) return;

    int    iso_count = CLAMP((int)desc.iso.count, 0, 8);
    float  iso_values[8];
    vec4_t iso_colors[8];
    float  iso_optical_densities[8] = { 0 };

    MEMCPY(iso_values, desc.iso.values, iso_count * sizeof(float));
    MEMCPY(iso_colors, desc.iso.colors, iso_count * sizeof(vec4_t));
    if (desc.iso.optical_densities) {
        MEMCPY(iso_optical_densities, desc.iso.optical_densities, iso_count * sizeof(float));
    }

    // Sort on iso value
    for (int i = 0; i < iso_count - 1; ++i) {
        for (int j = i + 1; j < iso_count; ++j) {
            if (iso_values[j] < iso_values[i]) {
                float  val_tmp = iso_values[i];
                vec4_t col_tmp = iso_colors[i];
                float  od_tmp = iso_optical_densities[i];
                iso_values[i] = iso_values[j];
                iso_colors[i] = iso_colors[j];
                iso_optical_densities[i] = iso_optical_densities[j];
                iso_values[j] = val_tmp;
                iso_colors[j] = col_tmp;
                iso_optical_densities[j] = od_tmp;
            }
        }
    }

    GLint bound_fbo;
    GLint bound_viewport[4];
    GLint bound_draw_buffer[8] = {0};
    GLint bound_draw_buffer_count = 0;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &bound_fbo);
    glGetIntegerv(GL_VIEWPORT, bound_viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &bound_draw_buffer[i]);
        // @NOTE: Assume that its tightly packed and if we stumple upon a zero draw buffer index, we enterpret that as the 'end'
        if (bound_draw_buffer[i] == GL_NONE) {
            bound_draw_buffer_count = i;
            break;
        }
    }

    if (gl.width < desc.render_target.width ||
        gl.height < desc.render_target.height)
    {
        gl.width = desc.render_target.width;
        gl.height = desc.render_target.height;
        gl::init_texture_2D(&gl.tex_entry,  gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_exit,   gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_result, gl.width, gl.height, GL_RGBA8);
    }

    const mat4_t model_to_view_matrix = mat4_mul(desc.matrix.view, desc.matrix.model);

    static float time = 0.0f;
    time += 1.0f / 100.0f;
    if (time > 100.0) time -= 100.0f;
    if (!desc.temporal.enabled) {
        time = 0.0f;
    }

    float tf_min = desc.dvr.min_tf_value;
    float tf_max = desc.dvr.max_tf_value;
    float tf_ext = tf_max - tf_min;
    float inv_tf_ext = tf_ext == 0 ? 1.0f : 1.0f / tf_ext;

    const float n1 = 1.0f;
    const float n2 = desc.shading.ior;
    const float F0 = powf((n1-n2)/(n1+n2), 2.0f);

    UniformData data;
    data.view_to_model_mat = mat4_inverse(model_to_view_matrix);
    data.model_to_view_mat = model_to_view_matrix;
    data.inv_proj_mat      = desc.matrix.inv_proj;
    data.model_view_proj_mat = desc.matrix.proj * model_to_view_matrix;
    data.inv_res = {1.f / (float)(desc.render_target.width), 1.f / (float)(desc.render_target.height)};
    data.time = time;
    data.gamma = desc.shading.gamma;
    data.clip_volume_min = desc.clip_volume.min;
    data.tf_min = tf_min;
    data.clip_volume_max = desc.clip_volume.max;
    data.tf_inv_ext = inv_tf_ext;
    data.gradient_spacing_world_space = desc.voxel_spacing;
    data.exposure = desc.shading.exposure;
    data.gradient_spacing_tex_space = data.view_to_model_mat * mat4_scale(desc.voxel_spacing.x, desc.voxel_spacing.y, desc.voxel_spacing.z);
    data.env_radiance = desc.shading.env_radiance;
    data.roughness = desc.shading.roughness;
    data.dir_radiance = desc.shading.dir_radiance;
    data.F0 = F0;

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformData), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    bool use_depth = desc.render_target.depth;

    if (use_depth) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, desc.render_target.depth);
    }

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);
    glBindVertexArray(gl.vao);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_entry, 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.tex_exit,  0);
    const GLuint draw_bufs[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, draw_bufs);

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(0, 0, desc.render_target.width, desc.render_target.height);

    glEnable(GL_CULL_FACE);
    {
        PUSH_GPU_SECTION("VOLUME ENTRY / EXIT");
        
        const GLuint prog = use_depth ? gl.program.entry_exit_depth : gl.program.entry_exit;
        const GLint uniform_block_index = glGetUniformBlockIndex(prog, "UniformData");
        const GLint uniform_loc_tex_depth = glGetUniformLocation(prog, "u_tex_depth");

        glUseProgram(prog);
        glUniform1i(uniform_loc_tex_depth, 0);
        glUniformBlockBinding(prog, uniform_block_index, 0);

        glCullFace(GL_FRONT);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
        POP_GPU_SECTION()
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.tex_entry);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gl.tex_exit);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, desc.texture.density_volume);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);

    if (desc.iso.use_color_volume) {
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_3D, desc.texture.color_volume);
    }

    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_result, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (desc.render_target.color) {
        ASSERT(glIsTexture(desc.render_target.color));
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc.render_target.color, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        if (desc.render_target.clear_color) {
            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);
        }
    } else {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
    }

    PUSH_GPU_SECTION("VOLUME RAYCASTING")
    {
        GLuint vol_prog = 0;
        if (desc.dvr.enabled) {
            if (desc.iso.enabled) {
                vol_prog = desc.iso.use_color_volume ? gl.program.dvr_and_iso_col_vol : gl.program.dvr_and_iso;
            } else {
                vol_prog = gl.program.dvr_only;
            }
        } else {
            if (desc.iso.use_color_volume) {
                vol_prog = gl.program.iso_col_vol;
            } else {
                vol_prog = desc.iso.enabled ? gl.program.iso_only : 0;
            }
        }

        if (vol_prog == 0) {
            MD_LOG_DEBUG("No raycasting shader program available for the current render description, skipping raycasting");
            return;
        }

        const GLint uniform_block_index             = glGetUniformBlockIndex(vol_prog, "UniformData");
        const GLint uniform_loc_tex_entry           = glGetUniformLocation(vol_prog, "u_tex_entry");
        const GLint uniform_loc_tex_exit            = glGetUniformLocation(vol_prog, "u_tex_exit");
        const GLint uniform_loc_tex_tf              = glGetUniformLocation(vol_prog, "u_tex_tf");
        const GLint uniform_loc_tex_density_volume  = glGetUniformLocation(vol_prog, "u_tex_density_volume");
        const GLint uniform_loc_tex_color_volume    = glGetUniformLocation(vol_prog, "u_tex_color_volume");
        const GLint uniform_loc_iso_values          = glGetUniformLocation(vol_prog, "u_iso.values");
        const GLint uniform_loc_iso_colors          = glGetUniformLocation(vol_prog, "u_iso.colors");
        const GLint uniform_loc_iso_optical_densities = glGetUniformLocation(vol_prog, "u_iso.optical_densities");
        const GLint uniform_loc_iso_count           = glGetUniformLocation(vol_prog, "u_iso.count");

        glUseProgram(vol_prog);

        glUniform1i(uniform_loc_tex_entry, 0);
        glUniform1i(uniform_loc_tex_exit,  1);
        glUniform1i(uniform_loc_tex_density_volume, 2);
        glUniform1i(uniform_loc_tex_tf, 3);
        glUniform1i(uniform_loc_tex_color_volume, 4);
        glUniform1fv(uniform_loc_iso_values, (GLsizei)iso_count, (const float*)iso_values);
        glUniform4fv(uniform_loc_iso_colors, (GLsizei)iso_count, (const float*)iso_colors);
        glUniform1fv(uniform_loc_iso_optical_densities, (GLsizei)iso_count, (const float*)iso_optical_densities);
        glUniform1i(uniform_loc_iso_count, (int)iso_count);
        glUniformBlockBinding(vol_prog, uniform_block_index, 0);

        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindVertexArray(0);
        glUseProgram(0);
    }
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glCullFace(GL_BACK);

    if (desc.render_target.color) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        if (bound_draw_buffer_count == 1)
            glDrawBuffer(bound_draw_buffer[0]);
        else
            glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
    }

}

}  // namespace volume
