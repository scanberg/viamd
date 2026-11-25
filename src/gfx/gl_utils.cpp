#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_file.h>
#include <core/md_path.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>

#include <string.h>
#include <stdlib.h>

#include "gl_utils.h"

bool gl::get_shader_compile_error(char* buffer, int max_length, GLuint shader) {
    GLint success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        int length;
        glGetShaderInfoLog(shader, max_length, &length, buffer);
        return true;
    } else {
        return false;
    }
}

bool gl::get_program_link_error(char* buffer, int max_length, GLuint program) {
    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
        int length;
        glGetProgramInfoLog(program, max_length, &length, buffer);
        return true;
    } else {
        return false;
    }
}

bool build_shader_src(md_writer_i* w, md_str_t src, md_str_t base_include_dir, md_allocator_i* alloc) {
    md_str_t line;
    while (md_parse_next_line(&line, &src)) {
        if (md_str_begins_with(line, MD_STR_LIT("#include "))) {
            md_str_t file = md_str_trim(md_str_substr(line, 9));
            if (!file || !(file.len > 2) || file[0] != '"' || file[file.len-1] != '"') {
                MD_LOG_ERROR("Failed to parse include file");
                return false;
            }
            file = md_str_substr(file, 1, file.len - 2);
            md_str_t path = md_alloc_printf(alloc, "%.*s%.*s", (int)base_include_dir.len, base_include_dir.ptr, (int)file.len, file.ptr);

            md_str_t inc_src = md_load_textfile(path, alloc);
            if (inc_src) {
                md_str_t base = md_path_folder(path);
                build_shader_src(w, inc_src, base, alloc);
            } else {
                MD_LOG_ERROR("Failed to open include file '%.*s'", (int)path.len, path.ptr);
                return false;
            }
        } else {
            md_write_str(w, line);
            md_write_char(w, '\n');
        }
    }

    return true;
}

GLuint gl::compile_shader_from_source(md_str_t src, GLenum type, md_str_t defines, md_str_t base_include_dir) {
    MD_ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
           type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    MD_DEFER { md_temp_end(temp); };

    GLuint shader = glCreateShader(type);
    md_str_writer_t sb = md_str_writer_create(alloc);
    
    if (defines) {
        md_str_t version_str = {};
        if (md_str_begins_with(src, MD_STR_LIT("#version "))) {
            if (!md_parse_next_line(&version_str, &src)) {
                MD_LOG_ERROR("Failed to extract version string!");
                return 0;
            }
            sb += version_str;
            sb += '\n';
            sb += defines;
            sb += '\n';
        }
        else {
            sb += defines;
            sb += '\n';
        }
    }

    build_shader_src(&sb.w, src, base_include_dir, alloc);

    md_str_t final_src = md_str_writer_to_str(&sb);
    glShaderSource(shader, 1, &final_src.ptr, 0);

    glCompileShader(shader);

    char buffer[1024];
    if (gl::get_shader_compile_error(buffer, sizeof(buffer), shader)) {
        MD_LOG_ERROR("%s\n", buffer);
        glDeleteShader(shader);
        shader = 0;
    }

    md_str_writer_free(&sb);

    return shader;
}

bool gl::attach_link_detach(GLuint program, const GLuint shaders[], int num_shaders) {
    MD_ASSERT(program);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];
    for (int i = 0; i < num_shaders; ++i) {
        MD_ASSERT(shaders[i]);
        glAttachShader(program, shaders[i]);
    }
    bool result = true;

    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, buffer_size, program)) {
        MD_LOG_ERROR("Linking program:\n%s\n", buffer);
        result = false;
    }

    for (int i = 0; i < num_shaders; ++i) {
        glDetachShader(program, shaders[i]);
    }
    return result;
}

bool gl::init_texture_2D(GLuint* texture, int width, int height, GLenum format) {
    MD_ASSERT(texture);    

    if (glIsTexture(*texture)) {
        int x, y;
        GLenum fmt;
        glBindTexture(GL_TEXTURE_2D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
        glBindTexture(GL_TEXTURE_2D, 0);
        if (width == x && height == y && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glBindTexture   (GL_TEXTURE_2D, *texture);
    glTexStorage2D  (GL_TEXTURE_2D, 1, format, width, height);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture   (GL_TEXTURE_2D, 0);

    return true;
}

bool gl::free_texture(GLuint* texture) {
    if (!glIsTexture(*texture)) return false;
    glDeleteTextures(1, texture);
    *texture = 0;
    return true;
}

static inline bool get_pixel_channel_type(GLenum& channel, GLenum& type, GLenum format) {
    switch (format) {
    case GL_RGBA8:
        channel = GL_RGBA;
        type = GL_UNSIGNED_BYTE;
        return true;
    case GL_RGBA32F:
        channel = GL_RGBA;
        type = GL_FLOAT;
        return true;
    case GL_R32F:
        channel = GL_RED;
        type = GL_FLOAT;
        return true;
    case GL_RGB10_A2:   // x in the low bits (md_gto_gpu_gradient_pack)
        channel = GL_RGBA;
        type = GL_UNSIGNED_INT_2_10_10_10_REV;
        return true;
    default:
        channel = 0;
        type = 0;
    }
    return false;
}

bool gl::set_texture_2D_data(GLuint texture, int level, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    int w, h;
    glBindTexture(GL_TEXTURE_2D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
    glTexSubImage2D(GL_TEXTURE_2D, level, 0, 0, w, h, pixel_channel, pixel_type, data);

    return true;
}
