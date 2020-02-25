#include <core/common.h>
#include <core/log.h>
#include <core/string_utils.h>
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

GLuint gl::compile_shader_from_file(CStringView filename, GLenum type, CStringView defines) {
    ASSERT(type == GL_VERTEX_SHADER || type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER || type == GL_COMPUTE_SHADER || type == GL_TESS_CONTROL_SHADER || type == GL_TESS_EVALUATION_SHADER);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];

    StringView shader_src = allocate_and_read_textfile(filename);
    defer { FREE(shader_src.cstr()); };

    if (!shader_src) {
        LOG_ERROR("Could not load shader from file %.*s\n", (i32)filename.size(), filename.cstr());
        return 0;
    }

    GLuint shader = glCreateShader(type);
    if (defines) {
        CStringView src = shader_src;
        StringBuffer<64> version_str;
        if (compare_n(shader_src, "#version", 8)) {
            version_str = extract_line(src);
        }
        const char* sources[5] = {version_str.cstr(), "\n", defines.cstr(), "\n", src.cstr()};
        glShaderSource(shader, 5, sources, 0);
    } else {
        const char* c_src = shader_src.cstr();
        glShaderSource(shader, 1, &c_src, 0);
    }

    glCompileShader(shader);
    if (gl::get_shader_compile_error(buffer, buffer_size, shader)) {
        LOG_ERROR("Compiling shader (%.*s):\n%s\n", (i32)filename.size(), filename.cstr(), buffer);
        return 0;
    }

    return shader;
}

bool gl::attach_link_detach(GLuint program, Array<const GLuint> shaders) {
    ASSERT(program);
    constexpr int buffer_size = 1024;
    char buffer[buffer_size];
    for (const auto shader : shaders) {
        ASSERT(shader);
        glAttachShader(program, shader);
    }
    defer {
        for (const auto shader : shaders) {
            glDetachShader(program, shader);
        }
    };

    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, buffer_size, program)) {
        LOG_ERROR("Linking program:\n%s\n", buffer);
        return false;
    }

    return true;
}

bool gl::attach_link_detach_with_transform_feedback(GLuint program, Array<const GLuint> shaders, Array<const char*> varyings, GLenum buffer_capture_mode) {
    ASSERT(program);
    ASSERT(buffer_capture_mode == GL_INTERLEAVED_ATTRIBS || buffer_capture_mode == GL_SEPARATE_ATTRIBS);

    constexpr int buffer_size = 1024;
    char buffer[buffer_size];
    for (const auto shader : shaders) {
        ASSERT(shader);
        glAttachShader(program, shader);
    }
    defer {
        for (const auto shader : shaders) {
            glDetachShader(program, shader);
        }
    };

    glTransformFeedbackVaryings(program, (GLsizei)varyings.size(), varyings.ptr, buffer_capture_mode);

    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, buffer_size, program)) {
        LOG_ERROR("Linking program:\n%s\n", buffer);
        return false;
    }

    return true;
}
