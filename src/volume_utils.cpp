#include "volume_utils.h"

#include "image.h"
#include "gfx/gl_utils.h"
#include <glm/gtc/type_ptr.hpp>

#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>

#include <stdio.h>
#include <algorithm>

namespace volume {

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ubo = 0;

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

    vec3 clip_plane_min; float _pad0[1];
    vec3 clip_plane_max; float _pad1[1];

    IsoSurface isosurface; int _pad2[3];

    vec3 gradient_spacing_world_space; float _pad3[1];
    mat3 gradient_spacing_tex_space;
};

void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.vert", GL_VERTEX_SHADER);
    GLuint f_shader_dvr_only = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_DVR");
    GLuint f_shader_iso_only = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_ISO");
    GLuint f_shader_dvr_and_iso = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/volume/raycaster.frag", GL_FRAGMENT_SHADER, "#define INCLUDE_DVR\n#define INCLUDE_ISO");
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
        constexpr uint8_t cube_strip[42] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
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
}

void shutdown() {}

void create_volume_texture(GLuint* texture, const ivec3& dim) {
    ASSERT(texture);
    if (*texture != 0 && glIsTexture(*texture)) {
        glDeleteTextures(1, texture);
    }
    glGenTextures(1, texture);
    glBindTexture(GL_TEXTURE_3D, *texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R16F, dim.x, dim.y, dim.z, 0, GL_RED, GL_FLOAT, nullptr);
    glBindTexture(GL_TEXTURE_3D, 0);
}

void free_volume_texture(GLuint texture) {
    if (glIsTexture(texture)) glDeleteTextures(1, &texture);
}

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
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glBindTexture(GL_TEXTURE_2D, 0);
    } else {
        LOG_WARNING("could not read TF ('%.s')", path.length(), path.cstr());
    }
}

void set_volume_texture_data(GLuint texture, ivec3 dim, void* data) {
    if (glIsTexture(texture)) {
        glBindTexture(GL_TEXTURE_3D, texture);
        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, dim.x, dim.y, dim.z, GL_RED, GL_FLOAT, data);
        glBindTexture(GL_TEXTURE_3D, 0);
    }
}

mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb) {
    vec3 ext = max_world_aabb - min_world_aabb;
    return mat4(vec4(ext.x, 0, 0, 0), vec4(0, ext.y, 0, 0), vec4(0, 0, ext.z, 0), vec4(min_world_aabb, 1));
}

mat4 compute_texture_to_model_matrix(const ivec3& dim) {
    (void)dim;
    return mat4(1);
    /*
const vec3 cell_ext = 1.f / vec3(dim);
const vec3 scl = 1.f + cell_ext;
return math::inverse(mat4(vec4(scl.x, 0, 0, 0), vec4(0, scl.y, 0, 0), vec4(0, 0, scl.z, 0), vec4(-0.5f * cell_ext, 1.f)));
    */
}

void save_volume_to_file(const Volume& volume, CStringView file) {
    StringBuffer<512> zstr = file;
    FILE* f = fopen(zstr.cstr(), "wb");
    defer { fclose(f); };

    if (!f) {
        LOG_ERROR("Could not open file %s", file);
        return;
    }
    fwrite(volume.voxel_data.data(), 1, volume.voxel_data.size_in_bytes(), f);
}

void render_volume_texture(const VolumeRenderDesc& desc) {
    if (!desc.direct_volume_rendering_enabled && !desc.isosurface_enabled) return;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    const mat4 model_to_view_matrix = desc.matrix.view * desc.matrix.model;

    UniformData data;
    data.view_to_model_mat = math::inverse(model_to_view_matrix);
    data.model_to_view_mat = model_to_view_matrix;
    data.inv_proj_mat = math::inverse(desc.matrix.proj);
    data.model_view_proj_mat = desc.matrix.proj * model_to_view_matrix;
    data.inv_res = vec2(1.f / (float)(viewport[2]), 1.f / (float)(viewport[3]));
    data.density_scale = desc.global_scaling.density;
    data.alpha_scale = desc.global_scaling.alpha;
    data.clip_plane_min = desc.clip_planes.min;
    data.clip_plane_max = desc.clip_planes.max;
    memcpy(&data.isosurface, &desc.isosurface, sizeof(IsoSurface));
    data.gradient_spacing_world_space = desc.voxel_spacing;
    data.gradient_spacing_tex_space = mat3(glm::scale(data.view_to_model_mat, desc.voxel_spacing));

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformData), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, desc.texture.depth);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, desc.texture.volume);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);

    const GLuint program = desc.direct_volume_rendering_enabled ?
        (desc.isosurface_enabled ? gl.program.dvr_and_iso : gl.program.dvr_only) : gl.program.iso_only;

    const GLint uniform_loc_tex_depth = glGetUniformLocation(program, "u_tex_depth");
    const GLint uniform_loc_tex_volume = glGetUniformLocation(program, "u_tex_volume");
    const GLint uniform_loc_tex_tf = glGetUniformLocation(program, "u_tex_tf");
    const GLint uniform_block_index = glGetUniformBlockIndex(program, "UniformData");

    glUseProgram(program);

    glUniform1i(uniform_loc_tex_depth, 0);
    glUniform1i(uniform_loc_tex_volume, 1);
    glUniform1i(uniform_loc_tex_tf, 2);
    glUniformBlockBinding(program, uniform_block_index, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glBindVertexArray(0);

    glUseProgram(0);

    glDisable(GL_CULL_FACE);
}

}  // namespace volume
