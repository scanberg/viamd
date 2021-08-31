#include "volumerender_utils.h"

#include "image.h"
#include "gfx/gl_utils.h"
#include "gfx/immediate_draw_utils.h"
#include <glm/gtc/type_ptr.hpp>

#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/file.h>

#include <stdio.h>
#include <algorithm>

namespace volume {

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ubo = 0;
    GLuint fbo = 0;

    struct {
        GLuint dvr_only = 0;
        GLuint iso_only = 0;
        GLuint dvr_and_iso = 0;
    } program;
} gl;

struct UniformData {
    mat4 view_to_model_mat;
    mat4 model_to_view_mat;
    mat4 inv_proj_mat;
    mat4 model_view_proj_mat;

    vec2 inv_res;
    float density_scale = 1.0;
    float alpha_scale = 1.0;

    vec3 clip_volume_min;
    float _pad0;
    vec3 clip_volume_max;
    float time;

    vec3 gradient_spacing_world_space;
    float _pad1;
    mat4 gradient_spacing_tex_space;
};

void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.vert", GL_VERTEX_SHADER);
    GLuint f_shader_dvr_only = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_DVR");
    GLuint f_shader_iso_only = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_ISO");
    GLuint f_shader_dvr_and_iso =
        gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_DVR\n#define INCLUDE_ISO");
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(f_shader_dvr_only);
        glDeleteShader(f_shader_iso_only);
        glDeleteShader(f_shader_dvr_and_iso);
    };

    if (v_shader == 0u || f_shader_dvr_only == 0u || f_shader_iso_only == 0u || f_shader_dvr_and_iso == 0u) {
        LOG_WARNING("shader compilation failed, shader program for raycasting will not be updated");
        return;
    }

    if (!gl.program.dvr_only) gl.program.dvr_only = glCreateProgram();
    if (!gl.program.iso_only) gl.program.iso_only = glCreateProgram();
    if (!gl.program.dvr_and_iso) gl.program.dvr_and_iso = glCreateProgram();

    {
        const GLuint shaders[] = {v_shader, f_shader_dvr_only};
        gl::attach_link_detach(gl.program.dvr_only, shaders);
    }
    {
        const GLuint shaders[] = {v_shader, f_shader_iso_only};
        gl::attach_link_detach(gl.program.iso_only, shaders);
    }
    {
        const GLuint shaders[] = {v_shader, f_shader_dvr_and_iso};
        gl::attach_link_detach(gl.program.dvr_and_iso, shaders);
    }

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

    if (!gl.fbo) {
        glGenFramebuffers(1, &gl.fbo);
    }
}

void shutdown() {}

void create_tf_texture(GLuint* texture, int* width, CStringView path) {
    ASSERT(texture);
    // load transfer function
    if (*texture == 0 || !glIsTexture(*texture)) {
        glGenTextures(1, texture);
    }

    Image img;
    defer { free_image(&img); };
    if (read_image(&img, path)) {
        if (width) {
            *width = img.width;
        }
        glBindTexture(GL_TEXTURE_2D, *texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    } else {
        LOG_WARNING("could not read TF ('%.s')", path.length(), path.cstr());
    }
}

bool init_texture_2D(GLuint* texture, int width, int height, GLenum format) {
    ASSERT(texture);
    ASSERT(GL_RGBA8);

    if (glIsTexture(*texture)) {
        int x, y, fmt;
        glBindTexture(GL_TEXTURE_2D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &fmt);
        glBindTexture(GL_TEXTURE_2D, 0);
        if (width == x && width == y && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glBindTexture(GL_TEXTURE_2D, *texture);
    glTexStorage2D(GL_TEXTURE_2D, 1, format, width, height);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    return true;
}

bool init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format) {
    ASSERT(texture);
    ASSERT(format == GL_R32F || format == GL_R8);

    if (glIsTexture(*texture)) {
        int x, y, z, fmt;
        glBindTexture(GL_TEXTURE_3D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &z);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, &fmt);
        glBindTexture(GL_TEXTURE_3D, 0);
        if (width == x && width == y && width == z && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glBindTexture(GL_TEXTURE_3D, *texture);
    glTexStorage3D(GL_TEXTURE_3D, 1, format, width, height, depth);
    // glTexImage3D(GL_TEXTURE_3D, 0, GL_R8, dim.x, dim.y, dim.z, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);

    return true;
}

bool free_texture(GLuint* texture) {
    if (!glIsTexture(*texture)) return false;
    glDeleteTextures(1, texture);
    *texture = 0;
    return true;
}

bool set_texture_2D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (format == GL_RGBA8) {
        pixel_channel = GL_RGBA;
        pixel_type = GL_UNSIGNED_BYTE;
    }

    if (pixel_channel == 0 || pixel_type == 0) return false;

    glBindTexture(GL_TEXTURE_2D, texture);
    int w, h;
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_2D, 0);

    return true;
}


bool set_texture_3D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (format == GL_R32F) {
        pixel_channel = GL_RED;
        pixel_type = GL_FLOAT;
    }

    if (pixel_channel == 0 || pixel_type == 0) return false;

    glBindTexture(GL_TEXTURE_3D, texture);
    int w, h, d;
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH, &d);

    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, w, h, d, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_3D, 0);

    return true;
}

/*
void set_volume_texture_data(GLuint texture, ivec3 dim, const float* data, float max_value) {
    if (glIsTexture(texture)) {
        uint8_t* rescaled_data = (uint8_t*)TMP_MALLOC(dim.x * dim.y * dim.z);
        defer { TMP_FREE(rescaled_data); };

        if (max_value > 0) {
            const int size = dim.x * dim.y * dim.z;
            for (int i = 0; i < size; ++i) {
                rescaled_data[i] = (uint8_t)((double)data[i] / (double)max_value * 255);
            }
        }
        glBindTexture(GL_TEXTURE_3D, texture);
        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, dim.x, dim.y, dim.z, GL_RED, GL_UNSIGNED_BYTE, rescaled_data);
        glBindTexture(GL_TEXTURE_3D, 0);
    }
}
*/

mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb) {
    const vec3 ext = max_world_aabb - min_world_aabb;
    return mat4(vec4(ext.x, 0, 0, 0), vec4(0, ext.y, 0, 0), vec4(0, 0, ext.z, 0), vec4(min_world_aabb, 1));
}

mat4 compute_world_to_model_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb) {
    const vec3 ext = max_world_aabb - min_world_aabb;
    return mat4(vec4(1.0f / ext.x, 0, 0, 0), vec4(0, 1.0f / ext.y, 0, 0), vec4(0, 0, 1.0f / ext.z, 0), vec4(-min_world_aabb, 1));
}

mat4 compute_model_to_texture_matrix(const ivec3& dim) {
    (void)dim;
    return mat4(1);
}

mat4 compute_texture_to_model_matrix(const ivec3& dim) {
    (void)dim;
    return mat4(1);
}

void write_to_file(const Volume& volume, CStringView file) {
    FILE* f = fopen(file, "wb");
    defer { fclose(f); };

    if (!f) {
        LOG_ERROR("Could not open file %s", file);
        return;
    }
    fwrite(volume.voxel_data.data(), 1, volume.voxel_data.size_in_bytes(), f);
}

void render_volume(const RenderDesc& desc) {
    if (!desc.direct_volume_rendering_enabled && !desc.isosurface_enabled && !desc.bounding_box_enabled) return;

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

    const mat4 model_to_view_matrix = desc.matrix.view * desc.matrix.model;

    static float time = 0.0f;
    time += 1.0f / 60.0f;
    if (time > 100.0) time -= 100.0f;

    if (desc.render_target.texture) {
        ASSERT(glIsTexture(desc.render_target.texture));
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.fbo);
        // We don't need a depth buffer to render the volume
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc.render_target.texture, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glViewport(0, 0, desc.render_target.width, desc.render_target.height);
    }

    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT);

    if (desc.bounding_box_enabled) {
        const vec3 min_box = {0,0,0};
        const vec3 max_box = {1,1,1};

        immediate::set_model_view_matrix(model_to_view_matrix);
        immediate::set_proj_matrix(desc.matrix.proj);

        immediate::draw_box_wireframe(desc.clip_volume.min, desc.clip_volume.max, immediate::COLOR_RED);
        immediate::draw_box_wireframe(min_box, max_box, immediate::COLOR_BLACK);

        immediate::flush();
    }

    UniformData data;
    data.view_to_model_mat = math::inverse(model_to_view_matrix);
    data.model_to_view_mat = model_to_view_matrix;
    data.inv_proj_mat = math::inverse(desc.matrix.proj);
    data.model_view_proj_mat = desc.matrix.proj * model_to_view_matrix;
    data.inv_res = vec2(1.f / (float)(desc.render_target.width), 1.f / (float)(desc.render_target.height));
    data.density_scale = desc.global_scaling.density;
    data.alpha_scale = desc.global_scaling.alpha;
    data.clip_volume_min = desc.clip_volume.min;
    data.clip_volume_max = desc.clip_volume.max;
    data.time = time;
    data.gradient_spacing_world_space = desc.voxel_spacing;
    data.gradient_spacing_tex_space = mat4(glm::scale(data.view_to_model_mat, desc.voxel_spacing));

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformData), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, desc.texture.depth);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, desc.texture.volume);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);

    const GLuint program =
        desc.direct_volume_rendering_enabled ? (desc.isosurface_enabled ? gl.program.dvr_and_iso : gl.program.dvr_only) : gl.program.iso_only;

    const GLint uniform_loc_tex_depth = glGetUniformLocation(program, "u_tex_depth");
    const GLint uniform_loc_tex_volume = glGetUniformLocation(program, "u_tex_volume");
    const GLint uniform_loc_tex_tf = glGetUniformLocation(program, "u_tex_tf");
    const GLint uniform_block_index = glGetUniformBlockIndex(program, "UniformData");
    const GLint uniform_loc_iso_values = glGetUniformLocation(program, "u_iso.values");
    const GLint uniform_loc_iso_colors = glGetUniformLocation(program, "u_iso.colors");
    const GLint uniform_loc_iso_count = glGetUniformLocation(program, "u_iso.count");

    glUseProgram(program);

    glUniform1i(uniform_loc_tex_depth, 0);
    glUniform1i(uniform_loc_tex_volume, 1);
    glUniform1i(uniform_loc_tex_tf, 2);
    glUniform1fv(uniform_loc_iso_values, IsoSurfaces::MaxCount, desc.isosurface.values);
    glUniform4fv(uniform_loc_iso_colors, IsoSurfaces::MaxCount, &desc.isosurface.colors[0][0]);
    glUniform1i(uniform_loc_iso_count, desc.isosurface.count);
    glUniformBlockBinding(program, uniform_block_index, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glBindVertexArray(0);

    glUseProgram(0);

    glDisable(GL_CULL_FACE);
    glDisable(GL_BLEND);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
    glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
    glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);

}

}  // namespace volume
