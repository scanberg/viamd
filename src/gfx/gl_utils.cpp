#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_str_builder.h>

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

bool build_shader_src(md_strb_t* builder, str_t src, str_t base_include_dir) {
    str_t line;
    while (str_extract_line(&line, &src)) {
        if (str_eq_cstr_n(line, "#include ", 9)) {
            str_t file = str_trim(str_substr(line, 9));
            if (!file || !(file.len > 2) || file[0] != '"' || file[file.len-1] != '"') {
                MD_LOG_ERROR("Failed to parse include file");
                return false;
            }
            file = str_substr(file, 1, file.len - 2);
            str_t path = str_printf(md_get_temp_allocator(), "%.*s%.*s", (int)base_include_dir.len, base_include_dir.ptr, (int)file.len, file.ptr);

            str_t inc_src = load_textfile(path, md_get_temp_allocator());
            if (inc_src) {
                str_t base = {};
                extract_folder_path(&base, path);
                build_shader_src(builder, inc_src, base);
            } else {
                MD_LOG_ERROR("Failed to open include file '%.*s'", (int)path.len, path.ptr);
                return false;
            }
        } else {
            md_strb_push_str(builder, line);
            md_strb_push_char(builder, '\n');
        }
    }

    return true;
}

GLuint gl::compile_shader_from_source(str_t src, GLenum type, str_t defines, str_t base_include_dir) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
           type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);

    GLuint shader = glCreateShader(type);
    md_strb_t sb = md_strb_create(md_get_temp_allocator());
    
    if (defines) {
        str_t version_str = {};
        if (str_eq_cstr_n(src, "#version ", 9)) {
            if (!str_extract_line(&version_str, &src)) {
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

    build_shader_src(&sb, src, base_include_dir);

    str_t final_src = md_strb_to_str(sb);
    glShaderSource(shader, 1, &final_src.ptr, 0);

    glCompileShader(shader);

    char buffer[1024];
    if (gl::get_shader_compile_error(buffer, sizeof(buffer), shader)) {
        MD_LOG_ERROR("%s\n", buffer);
        glDeleteShader(shader);
        shader = 0;
    }

    md_strb_free(&sb);

    return shader;
}

GLuint gl::compile_shader_from_file(str_t filename, GLenum type, str_t defines) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER ||
        type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);

    str_t src = load_textfile(filename, md_get_temp_allocator());
    if (!src) {
        MD_LOG_ERROR("Failed to open source file for shader '%.*s'", (int)src.len, src.ptr);
        return 0;
    }

    GLuint shader = glCreateShader(type);
    md_strb_t sb = md_strb_create(md_get_temp_allocator());

    if (defines) {
        str_t version_str = {};
        if (str_eq_cstr_n(src, "#version ", 9)) {
            if (!str_extract_line(&version_str, &src)) {
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

    str_t folder_path;
    extract_folder_path(&folder_path, filename);
    build_shader_src(&sb, src, folder_path);

    str_t final_src = md_strb_to_str(sb);
    glShaderSource(shader, 1, &final_src.ptr, 0);
    md_strb_free(&sb);

    glCompileShader(shader);

    char buffer[1024];
    if (gl::get_shader_compile_error(buffer, sizeof(buffer), shader)) {
        MD_LOG_ERROR("%s\n", buffer);
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
        MD_LOG_ERROR("Linking program:\n%s\n", buffer);
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
        MD_LOG_ERROR("Linking program:\n%s\n", buffer);
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
        int x;
        GLenum fmt;
        glBindTexture(GL_TEXTURE_1D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_WIDTH, &x);
        glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
        glBindTexture(GL_TEXTURE_1D, 0);
        if (width == x && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glBindTexture  (GL_TEXTURE_1D, *texture);
    glTexStorage1D (GL_TEXTURE_1D, 1, format, width);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_1D, 0);

    return true;
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

bool gl::init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format) {
    ASSERT(texture);
    ASSERT(format == GL_R32F || format == GL_R16F || format == GL_R8);

    if (glIsTexture(*texture)) {
        int x, y, z;
        GLenum fmt;
        glBindTexture(GL_TEXTURE_3D, *texture);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &x);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &z);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&fmt);
        glBindTexture(GL_TEXTURE_3D, 0);

        if (width == x && height == y && depth == z && format == fmt)
            return true;
        else
            glDeleteTextures(1, texture);
    }

    glGenTextures(1, texture);
    glBindTexture  (GL_TEXTURE_3D, *texture);
    glTexStorage3D (GL_TEXTURE_3D, 1, format, width, height, depth);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glBindTexture  (GL_TEXTURE_3D, 0);

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
    glBindTexture(GL_TEXTURE_1D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_WIDTH, &w);
    glTexSubImage1D(GL_TEXTURE_1D, 0, 0, w, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_1D, 0);
    
    return true;
}

bool gl::set_texture_2D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    get_pixel_channel_type(pixel_channel, pixel_type, format);
    if (pixel_channel == 0 || pixel_type == 0) return false;

    int w, h;
    glBindTexture(GL_TEXTURE_2D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, pixel_channel, pixel_type, data);

    return true;
}

bool gl::set_texture_3D_data(GLuint texture, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    get_pixel_channel_type(pixel_channel, pixel_type, format);
    if (pixel_channel == 0 || pixel_type == 0) return false;

    int w, h, d;
    glBindTexture(GL_TEXTURE_3D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &d);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, w, h, d, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_3D, 0);

    return true;
}