#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_string_builder.h>

#include <string.h>

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

bool build_shader_src(md_string_builder_t* builder, str_t src, str_t base_include_dir) {
    str_t line;
    while (extract_line(&line, &src)) {
        if (str_equal_cstr_n(line, "#include ", 9)) {
            str_t file = str_trim_whitespace(substr(line, 9));
            if (!file || !(file.len > 2) || file[0] != '"' || file[file.len-1] != '"') {
                md_printf(MD_LOG_TYPE_ERROR, "Failed to parse include file");
                return false;
            }
            file = substr(file, 1, file.len - 2);
            str_t path = alloc_printf(default_temp_allocator, "%.*s%.*s", (int)base_include_dir.len, base_include_dir.ptr, (int)file.len, file.ptr);

            str_t inc_src = load_textfile(path, default_temp_allocator);
            if (inc_src) {
                build_shader_src(builder, inc_src, extract_path_without_file(path));
            } else {
                md_printf(MD_LOG_TYPE_ERROR, "Failed to open include file '%.*s'", (int)path.len, path.ptr);
                return false;
            }
        } else {
            md_string_builder_append_str(builder, line);
        }
    }

    return true;
}

GLuint gl::compile_shader_from_source(str_t src, GLenum type, str_t defines, str_t base_include_dir) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
           type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);

    GLuint shader = glCreateShader(type);
    md_string_builder_t builder = {0};
    md_string_builder_init(&builder, default_temp_allocator);
    
    if (defines) {
        str_t version_str = {};
        if (str_equal_cstr_n(src, "#version ", 9)) {
            if (!extract_line(&version_str, &src)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to extract version string!");
                return 0;
            }
            md_string_builder_append_str(&builder, version_str);
            md_string_builder_append_str(&builder, defines);
            md_string_builder_append_str(&builder, MAKE_STR("\n"));
        }
        else {
            md_string_builder_append_str(&builder, defines);
            md_string_builder_append_str(&builder, MAKE_STR("\n"));
        }
    }

    build_shader_src(&builder, src, base_include_dir);

    str_t final_src = md_string_builder_to_string(&builder);
    glShaderSource(shader, 1, &final_src.ptr, 0);
    md_string_builder_free(&builder);

    glCompileShader(shader);

    char buffer[1024];
    if (gl::get_shader_compile_error(buffer, sizeof(buffer), shader)) {
        md_printf(MD_LOG_TYPE_ERROR, "%s\n", buffer);
        return 0;
    }

    return shader;
}

GLuint gl::compile_shader_from_file(str_t filename, GLenum type, str_t defines) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
        type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);

    str_t src = load_textfile(filename, default_temp_allocator);
    if (!src) {
        md_printf(MD_LOG_TYPE_ERROR, "Failed to open source file for shader '%.*s'", (int)src.len, src.ptr);
        return 0;
    }

    GLuint shader = glCreateShader(type);
    md_string_builder_t builder = { 0 };
    md_string_builder_init(&builder, default_temp_allocator);

    if (defines) {
        str_t version_str = {};
        if (str_equal_cstr_n(src, "#version ", 9)) {
            if (!extract_line(&version_str, &src)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to extract version string!");
                return 0;
            }
            md_string_builder_append_str(&builder, version_str);
            md_string_builder_append_str(&builder, defines);
            md_string_builder_append_str(&builder, MAKE_STR("\n"));
        }
        else {
            md_string_builder_append_str(&builder, defines);
            md_string_builder_append_str(&builder, MAKE_STR("\n"));
        }
    }

    build_shader_src(&builder, src, extract_path_without_file(filename));

    str_t final_src = md_string_builder_to_string(&builder);
    glShaderSource(shader, 1, &final_src.ptr, 0);
    md_string_builder_free(&builder);

    glCompileShader(shader);

    char buffer[1024];
    if (gl::get_shader_compile_error(buffer, sizeof(buffer), shader)) {
        md_printf(MD_LOG_TYPE_ERROR, "%s\n", buffer);
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

bool gl::init_texture_1D(GLuint* texture, int width, GLenum format) {
    ASSERT(texture);    

    if (glIsTexture(*texture)) {
        int x, y;
        GLenum fmt;
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_WIDTH, &x);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
        glBindTexture(GL_TEXTURE_1D, 0);
        if (width == x && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glCreateTextures(GL_TEXTURE_1D, 1, texture);
    glTextureStorage1D(*texture, 1, format, width);
    glTextureParameteri(*texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    return true;
}

bool gl::init_texture_2D(GLuint* texture, int width, int height, GLenum format) {
    ASSERT(texture);    

    if (glIsTexture(*texture)) {
        int x, y;
        GLenum fmt;

        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);

        if (width == x && height == y && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glCreateTextures(GL_TEXTURE_2D, 1, texture);
    glTextureStorage2D(*texture, 1, format, width, height);
    glTextureParameteri(*texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    return true;
}

bool gl::init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format) {
    ASSERT(texture);
    ASSERT(format == GL_R32F || format == GL_R8);

    if (glIsTexture(*texture)) {
        int x, y, z;
        GLenum fmt;
        
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_DEPTH,  &z);
        glGetTextureLevelParameteriv(*texture, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
        
        if (width == x && height == y && depth == z && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glCreateTextures(GL_TEXTURE_3D, 1, texture);
    glTextureStorage3D(*texture, 1, format, width, height, depth);
    glTextureParameteri(*texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTextureParameteri(*texture, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return true;
}

bool gl::free_texture(GLuint* texture) {
    if (!glIsTexture(*texture)) return false;
    glDeleteTextures(1, texture);
    *texture = 0;
    return true;
}

static inline void get_pixel_channel_type(GLenum& channel, GLenum& type, GLenum format) {
    switch (format) {
    case GL_RGBA8:
        channel = GL_RGBA;
        type = GL_UNSIGNED_BYTE;
        break;
    case GL_RGBA32F:
        channel = GL_RGBA;
        type = GL_FLOAT;
        break;
    case GL_R32F:
        channel = GL_RED;
        type = GL_FLOAT;
        break;
    default:
        channel = 0;
        type = 0;
    }
}

bool gl::set_texture_1D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    get_pixel_channel_type(pixel_channel, pixel_type, format);
    if (pixel_channel == 0 || pixel_type == 0) return false;

    int w;
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_WIDTH, &w);
    glTextureSubImage1D(texture, 0, 0, w, pixel_channel, pixel_type, data);
    
    return true;
}

bool gl::set_texture_2D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    get_pixel_channel_type(pixel_channel, pixel_type, format);
    if (pixel_channel == 0 || pixel_type == 0) return false;

    int w, h;
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_WIDTH, &w);
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_HEIGHT, &h);

    glTextureSubImage2D(texture, 0, 0, 0, w, h, pixel_channel, pixel_type, data);

    return true;
}

bool gl::set_texture_3D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    get_pixel_channel_type(pixel_channel, pixel_type, format);
    if (pixel_channel == 0 || pixel_type == 0) return false;

    int w, h, d;
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_WIDTH, &w);
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTextureLevelParameteriv(texture, 0, GL_TEXTURE_DEPTH, &d);

    glTextureSubImage3D(texture, 0, 0, 0, 0, w, h, d, pixel_channel, pixel_type, data);

    return true;
}