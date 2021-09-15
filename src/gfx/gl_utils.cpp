#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>
#include <core/md_allocator.h>

#include <string.h>

#include "gl_utils.h"

static GLuint fbo = 0;

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

GLuint gl::compile_shader_from_source(const char* shader_src, GLenum type, const char* defines) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
           type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];

    GLuint shader = glCreateShader(type);
    if (defines) {
        str_t src = {shader_src, (int64_t)strlen(shader_src)};
        str_t version_str = {};
        if (compare_str_cstr_n(src, "#version", 8)) {
            if (!extract_line(&version_str, &src)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to extract version string!");
                return 0;
            }
            // We need a zero terminated cstr of version str.
            version_str = copy_str(version_str, default_temp_allocator);
            const char* sources[5] = {version_str.ptr, "\n", defines, "\n", src.ptr};
            glShaderSource(shader, 5, sources, 0);
            free_str(version_str, default_temp_allocator);

        }
        else {
            const char* sources[3] = {defines, "\n", src.ptr};
            glShaderSource(shader, 3, sources, 0);
        }
    } else {
        glShaderSource(shader, 1, &shader_src, 0);
    }

    glCompileShader(shader);
    if (gl::get_shader_compile_error(buffer, buffer_size, shader)) {
        md_printf(MD_LOG_TYPE_ERROR, "Compiling shader from source: \n%s\n", buffer);
        return 0;
    }

    return shader;
}

GLuint gl::compile_shader_from_file(const char* filename, GLenum type, const char* defines) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER || type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];

    str_t file = {filename, (int64_t)strlen(filename)};
    str_t src = load_textfile(file, default_temp_allocator);
    //defer { FREE(shader_src.cstr()); };

    if (!src.ptr) {
        md_printf(MD_LOG_TYPE_ERROR, "Could not load shader from file %.*s\n", (int)file.len, file.ptr);
        return 0;
    }

    GLuint shader = glCreateShader(type);
    if (defines) {
        str_t version_str = {};
        if (compare_str_cstr_n(src, "#version", 8)) {
            if (!extract_line(&version_str, &src)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to extract version string!");
                return 0;
            }
        }
        // We need a zero terminated cstr of version str.
        version_str = copy_str(version_str, default_temp_allocator);
        const char* sources[5] = {version_str.ptr, "\n", defines, "\n", src.ptr};
        glShaderSource(shader, 5, sources, 0);

        free_str(version_str, default_temp_allocator);
    } else {
        glShaderSource(shader, 1, &src.ptr, 0);
    }

    glCompileShader(shader);
    if (gl::get_shader_compile_error(buffer, buffer_size, shader)) {
        md_printf(MD_LOG_TYPE_ERROR, "Compiling shader (%.*s):\n%s\n", (int)file.len, file.ptr, buffer);
        return 0;
    }

    return shader;
}

bool gl::attach_link_detach(GLuint program, const GLuint shaders[], int num_shaders) {
    ASSERT(program);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];
    for (int i = 0; i < num_shaders; ++i) {
        ASSERT(shaders[i]);
        glAttachShader(program, shaders[i]);
    }
    bool result = true;

    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, buffer_size, program)) {
        md_printf(MD_LOG_TYPE_ERROR, "Linking program:\n%s\n", buffer);
        result = false;
    }

    for (int i = 0; i < num_shaders; ++i) {
        glDetachShader(program, shaders[i]);
    }
    return result;
}

bool gl::attach_link_detach_with_transform_feedback(GLuint program, const GLuint shaders[], int num_shaders, const char* varyings[], int num_varyings, GLenum buffer_capture_mode) {
    ASSERT(program);
    ASSERT(buffer_capture_mode == GL_INTERLEAVED_ATTRIBS || buffer_capture_mode == GL_SEPARATE_ATTRIBS);

    constexpr int buffer_size = 1024;
    char buffer[buffer_size];
    for (int i = 0; i < num_shaders; ++i) {
        ASSERT(shaders[i]);
        glAttachShader(program, shaders[i]);
    }

    glTransformFeedbackVaryings(program, num_varyings, varyings, buffer_capture_mode);

    bool result = true;

    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, buffer_size, program)) {
        md_printf(MD_LOG_TYPE_ERROR, "Linking program:\n%s\n", buffer);
        result = false;
    }

    for (int i = 0; i < num_shaders; ++i) {
        glDetachShader(program, shaders[i]);
    }
    return result;
}

bool gl::init_texture_2D(GLuint* texture, int width, int height, GLenum format) {
    ASSERT(texture);    

    if (glIsTexture(*texture)) {
        int x, y;
        GLenum fmt;
        glBindTexture(GL_TEXTURE_2D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
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

bool gl::init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format) {
    ASSERT(texture);
    ASSERT(format == GL_R32F || format == GL_R8);

    if (glIsTexture(*texture)) {
        int x, y, z;
        GLenum fmt;
        glBindTexture(GL_TEXTURE_3D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &z);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
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

bool gl::free_texture(GLuint* texture) {
    if (!glIsTexture(*texture)) return false;
    glDeleteTextures(1, texture);
    *texture = 0;
    return true;
}

bool gl::set_texture_2D_data(GLuint texture, const void* data, GLenum format) {
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

bool gl::set_texture_3D_data(GLuint texture, const void* data, GLenum format) {
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