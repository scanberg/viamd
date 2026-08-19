#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>

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

bool build_shader_src(md_strb_t* builder, str_t src, str_t base_include_dir, md_allocator_i* alloc) {
    str_t line;
    while (str_extract_line(&line, &src)) {
        if (str_eq_cstr_n(line, "#include ", 9)) {
            str_t file = str_trim(str_substr(line, 9));
            if (!file || !(file.len > 2) || file[0] != '"' || file[file.len-1] != '"') {
                MD_LOG_ERROR("Failed to parse include file");
                return false;
            }
            file = str_substr(file, 1, file.len - 2);
            str_t path = str_printf(alloc, "%.*s%.*s", (int)base_include_dir.len, base_include_dir.ptr, (int)file.len, file.ptr);

            str_t inc_src = load_textfile(path, alloc);
            if (inc_src) {
                str_t base = {};
                extract_folder_path(&base, path);
                build_shader_src(builder, inc_src, base, alloc);
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
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    defer { md_temp_end(temp); };

    GLuint shader = glCreateShader(type);
    md_strb_t sb = md_strb_create(alloc);
    
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

    build_shader_src(&sb, src, base_include_dir, alloc);

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
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    defer { md_temp_end(temp); };

    str_t src = load_textfile(filename, alloc);
    if (!src) {
        MD_LOG_ERROR("Failed to open source file for shader '%.*s'", (int)src.len, src.ptr);
        return 0;
    }

    GLuint shader = glCreateShader(type);
    md_strb_t sb = md_strb_create(alloc);

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
    build_shader_src(&sb, src, folder_path, alloc);

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
    ASSERT(format == GL_R32F || format == GL_R16F || format == GL_R8 || format == GL_RGBA8);

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

bool gl::get_texture_dim(int dim[3], GLuint tex) {
    if (glIsTexture(tex)) {
        glBindTexture(GL_TEXTURE_3D, tex);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &dim[0]);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &dim[1]);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &dim[2]);
        glBindTexture(GL_TEXTURE_3D, 0);
        return true;
    }
    return false;
}

bool gl::get_texture_format(GLenum* format, GLuint tex) {
    if (glIsTexture(tex)) {
        glBindTexture(GL_TEXTURE_3D, tex);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);
        glBindTexture(GL_TEXTURE_3D, 0);
        return true;
    }
    return false;
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
    default:
        channel = 0;
        type = 0;
    }
    return false;
}

static inline size_t get_bytes_per_pixel(GLenum format) {
    switch (format) {
    case GL_RGBA8:
        return 4 * sizeof(uint8_t);
    case GL_RGBA32F:
        return 4 * sizeof(float);
    case GL_R32F:
        return sizeof(float);
    default:
        return 0;
    }
}

bool gl::set_texture_1D_data(GLuint texture, int level, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    int w;
    glBindTexture(GL_TEXTURE_1D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_WIDTH, &w);
    glTexSubImage1D(GL_TEXTURE_1D, level, 0, w, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_1D, 0);
    
    return true;
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

bool gl::set_texture_3D_data(GLuint texture, int level, const void* data, GLenum format) {
    if (!glIsTexture(texture)) return false;

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    int w, h, d;
    glBindTexture(GL_TEXTURE_3D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &d);

    glTexSubImage3D(GL_TEXTURE_3D, level, 0, 0, 0, w, h, d, pixel_channel, pixel_type, data);
    glBindTexture(GL_TEXTURE_3D, 0);

    return true;
}

// --- Streaming uploads through a pixel unpack buffer -------------------------
//
// A small cache of persistently mapped buffers, keyed by power-of-two size
// class. A buffer goes back into the cache carrying the fence from its last
// upload; it is only handed out again once the driver has signalled that fence,
// so we never overwrite bytes the GPU is still reading.

namespace {

struct PboSlot {
    GLuint  id      = 0;
    void*   mapped  = nullptr;
    size_t  size    = 0;
    GLsync  fence   = nullptr;
    bool    in_use  = false;
};

// Two per size class is enough to keep one upload in flight while the next is
// being filled; volumes here are large, so a deeper ring mostly wastes address
// space. Eight slots covers several size classes at once.
constexpr int    PBO_SLOT_COUNT = 8;
constexpr size_t PBO_MIN_SIZE   = 1u << 20;    // don't bother below 1 MB
// These are persistently mapped, so they are resident for as long as they are
// cached. A volume can be hundreds of MB, and eight unbounded slots would be
// several GB of address space and driver allocation; cap the total instead.
constexpr size_t PBO_TOTAL_BUDGET = 384u << 20;

PboSlot g_pbo[PBO_SLOT_COUNT];
int     g_pbo_active = -1;      // slot handed out by pbo_upload_begin()
size_t  g_pbo_bytes  = 0;       // total currently allocated

size_t pbo_size_class(size_t size) {
    size_t c = PBO_MIN_SIZE;
    while (c < size) c <<= 1;
    return c;
}

// True once the driver is done reading the slot's last upload.
bool pbo_slot_ready(PboSlot& s) {
    if (!s.fence) return true;
    GLenum r = glClientWaitSync(s.fence, 0, 0);
    if (r == GL_ALREADY_SIGNALED || r == GL_CONDITION_SATISFIED) {
        glDeleteSync(s.fence);
        s.fence = nullptr;
        return true;
    }
    return false;
}

void pbo_slot_release(PboSlot& s) {
    if (s.fence)  { glDeleteSync(s.fence); s.fence = nullptr; }
    if (s.mapped) {
        glBindBuffer(GL_PIXEL_UNPACK_BUFFER, s.id);
        glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);
        glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
        s.mapped = nullptr;
    }
    if (s.id) { glDeleteBuffers(1, &s.id); s.id = 0; }
    g_pbo_bytes -= s.size;
    s.size   = 0;
    s.in_use = false;
}

bool pbo_slot_create(PboSlot& s, size_t size) {
    // glBufferStorage is GL 4.4; viamd's GPU path already requires 4.3+, so
    // check rather than assume and let the caller fall back if it is absent.
    if (!glBufferStorage || !glMapBufferRange || !glFenceSync) return false;

    while (glGetError() != GL_NO_ERROR) {}   // so the check below is about us

    glGenBuffers(1, &s.id);
    if (!s.id) return false;
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, s.id);

    const GLbitfield storage = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
    glBufferStorage(GL_PIXEL_UNPACK_BUFFER, (GLsizeiptr)size, nullptr, storage);
    if (glGetError() != GL_NO_ERROR) {
        glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
        glDeleteBuffers(1, &s.id);
        s.id = 0;
        return false;
    }

    s.mapped = glMapBufferRange(GL_PIXEL_UNPACK_BUFFER, 0, (GLsizeiptr)size, storage);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    if (!s.mapped) {
        glDeleteBuffers(1, &s.id);
        s.id = 0;
        return false;
    }
    s.size = size;
    g_pbo_bytes += size;
    return true;
}

// Frees idle buffers, largest first, until `want` more bytes fit in the budget.
// Returns false if it cannot be made to fit, in which case the caller falls
// back to a plain client-pointer upload rather than forcing an eviction of
// something still in use.
bool pbo_make_room(size_t want) {
    if (want > PBO_TOTAL_BUDGET) return false;
    while (g_pbo_bytes + want > PBO_TOTAL_BUDGET) {
        int victim = -1;
        for (int i = 0; i < PBO_SLOT_COUNT; ++i) {
            PboSlot& s = g_pbo[i];
            if (!s.id || s.in_use || !pbo_slot_ready(s)) continue;
            if (victim < 0 || s.size > g_pbo[victim].size) victim = i;
        }
        if (victim < 0) return false;
        pbo_slot_release(g_pbo[victim]);
    }
    return true;
}

}  // namespace

void* gl::pbo_upload_begin(size_t size) {
    if (size == 0 || g_pbo_active >= 0) return nullptr;   // one at a time

    const size_t want = pbo_size_class(size);

    // Reuse an idle slot of the right class.
    for (int i = 0; i < PBO_SLOT_COUNT; ++i) {
        PboSlot& s = g_pbo[i];
        if (s.id && !s.in_use && s.size == want && pbo_slot_ready(s)) {
            s.in_use = true;
            g_pbo_active = i;
            return s.mapped;
        }
    }
    // Otherwise take a free slot, evicting an idle one of the wrong size if need be.
    int idx = -1;
    for (int i = 0; i < PBO_SLOT_COUNT && idx < 0; ++i) if (!g_pbo[i].id) idx = i;
    for (int i = 0; i < PBO_SLOT_COUNT && idx < 0; ++i) {
        if (!g_pbo[i].in_use && pbo_slot_ready(g_pbo[i])) { pbo_slot_release(g_pbo[i]); idx = i; }
    }
    if (idx < 0) return nullptr;   // everything busy; caller falls back

    if (!pbo_make_room(want))             return nullptr;
    if (!pbo_slot_create(g_pbo[idx], want)) return nullptr;
    g_pbo[idx].in_use = true;
    g_pbo_active = idx;
    return g_pbo[idx].mapped;
}

bool gl::pbo_upload_end_texture_3D(GLuint texture, int level, GLenum format) {
    if (g_pbo_active < 0) return false;
    PboSlot& s = g_pbo[g_pbo_active];
    g_pbo_active = -1;
    s.in_use = false;

    GLenum pixel_channel = 0, pixel_type = 0;
    if (!glIsTexture(texture) || !get_pixel_channel_type(pixel_channel, pixel_type, format)) {
        return false;
    }

    int w, h, d;
    glBindTexture(GL_TEXTURE_3D, texture);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &d);

    // With a buffer bound to GL_PIXEL_UNPACK_BUFFER the last argument is an
    // offset into it, not a pointer.
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, s.id);
    glTexSubImage3D(GL_TEXTURE_3D, level, 0, 0, 0, w, h, d, pixel_channel, pixel_type, (const void*)0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    glBindTexture(GL_TEXTURE_3D, 0);

    // Tells us when the buffer may be refilled.
    if (s.fence) glDeleteSync(s.fence);
    s.fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
    return true;
}

void gl::pbo_upload_abort(void) {
    if (g_pbo_active < 0) return;
    g_pbo[g_pbo_active].in_use = false;
    g_pbo_active = -1;
}

void gl::pbo_upload_shutdown(void) {
    g_pbo_active = -1;
    for (int i = 0; i < PBO_SLOT_COUNT; ++i) pbo_slot_release(g_pbo[i]);
}

bool gl::clear_texture_1D(GLuint texture, int level) {
    if (!glIsTexture(texture)) return false;

    glBindTexture(GL_TEXTURE_1D, texture);

    GLenum format = 0;
    glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;
    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    bool result = false;

    if (glClearTexImage) {
        float clear_value[4] = { 0 };
        glClearTexImage(texture, level, pixel_channel, pixel_type, clear_value);
        result = true;
    }
    else {
        int w;
        glGetTexLevelParameteriv(GL_TEXTURE_1D, 0, GL_TEXTURE_WIDTH,  &w);
        size_t bytes = w * get_bytes_per_pixel(format);
        if (bytes > 0) {
            md_temp_scope_t temp_scope = md_temp_begin();
            defer { md_temp_end(temp_scope); };
            void* data = md_temp_alloc_zero(temp_scope, bytes);
            glTexSubImage1D(GL_TEXTURE_1D, level, 0, w, pixel_channel, pixel_type, data);
            result = true;
        }
    }

    glBindTexture(GL_TEXTURE_1D, 0);
    return result;
}

bool gl::clear_texture_2D(GLuint texture, int level) {
    if (!glIsTexture(texture)) return false;
    
    glBindTexture(GL_TEXTURE_2D, texture);

    GLenum format = 0;
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;

    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    bool result = false;

    if (glClearTexImage) {
        float clear_value[4] = { 0 };
        glClearTexImage(texture, level, pixel_channel, pixel_type, clear_value);
        result = true;
    }
    else {
        int w, h;
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &w);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

        size_t bytes = w * h * get_bytes_per_pixel(format);
        if (bytes > 0) {
            md_temp_scope_t temp_scope = md_temp_begin();
            defer { md_temp_end(temp_scope); };
            void* data = md_temp_alloc_zero(temp_scope, bytes);
            glTexSubImage2D(GL_TEXTURE_2D, level, 0, 0, w, h, pixel_channel, pixel_type, data);
            result = true;
        }
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    return result;
}

bool gl::clear_texture_3D(GLuint texture, int level) {
    if (!glIsTexture(texture)) return false;
    
    glBindTexture(GL_TEXTURE_3D, texture);

    GLenum format = 0;
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);

    GLenum pixel_channel = 0;
    GLenum pixel_type = 0;
    if (!get_pixel_channel_type(pixel_channel, pixel_type, format)) return false;

    bool result = false;

    if (glClearTexImage) {
        float clear_value[4] = { 0 };
        glClearTexImage(texture, level, pixel_channel, pixel_type, clear_value);
        result = true;
    }
    else {
        int w, h, d;
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH,  &w);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &h);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_DEPTH,  &d);

        size_t bytes = w * h * d * get_bytes_per_pixel(format);
        if (bytes > 0) {
            md_temp_scope_t temp_scope = md_temp_begin();
            defer { md_temp_end(temp_scope); };
            void* data = md_temp_alloc_zero(temp_scope, bytes);
            glTexSubImage3D(GL_TEXTURE_3D, level, 0, 0, 0, w, h, d, pixel_channel, pixel_type, data);
            result = true;
        }
    }

    glBindTexture(GL_TEXTURE_3D, 0);
    return result;
}