#pragma once

#include "gl.h"
#include <core/md_str.h>

namespace gl {

bool get_shader_compile_error(char* buffer, int max_length, GLuint shader);
bool get_program_link_error(char* buffer, int max_length, GLuint program);
GLuint compile_shader_from_source(str_t source, GLenum shader_type, str_t defines = {}, str_t base_include_dir = {});
GLuint compile_shader_from_file(str_t filename, GLenum shader_type, str_t defines = {});
bool attach_link_detach(GLuint program, const GLuint shaders[], int num_shaders);
bool attach_link_detach_with_transform_feedback(GLuint program, const GLuint shaders[], int num_shaders, const char* varyings[], int num_varyings, GLenum buffer_capture_mode);

bool init_texture_1D(GLuint* texture, int width, GLenum format);
bool init_texture_2D(GLuint* texture, int width, int height, GLenum format);
bool init_texture_3D(GLuint* texture, int width, int height, int depth, GLenum format);

bool get_texture_dim(int dim[3], GLuint texture);
bool get_texture_format(GLenum* format, GLuint texture);

bool free_texture(GLuint* texture);

// Set the data for the entire texture
bool set_texture_1D_data(GLuint texture, int level, const void* data, GLenum format);
bool set_texture_2D_data(GLuint texture, int level, const void* data, GLenum format);
bool set_texture_3D_data(GLuint texture, int level, const void* data, GLenum format);

bool clear_texture_1D(GLuint texture, int level);
bool clear_texture_2D(GLuint texture, int level);
bool clear_texture_3D(GLuint texture, int level);

// Streaming uploads through a pixel unpack buffer.
//
// set_texture_3D_data() hands a client pointer to glTexSubImage3D, which makes
// the driver copy into its own pinned staging before it can start the DMA, and
// does not return until it has. Going through a persistently mapped PBO writes
// once into memory the GPU can already see and lets the transfer overlap.
//
// Intended shape:
//
//     void* dst = gl::pbo_upload_begin(size);       // NULL => no PBO available
//     if (dst) {
//         memcpy(dst, src, size);
//         gl::pbo_upload_end_texture_3D(texture, 0, GL_R32F);
//     } else {
//         gl::set_texture_3D_data(texture, 0, src, GL_R32F);   // fallback
//     }
//
// Buffers are cached by power-of-two size class and recycled once the driver
// signals it has consumed them, so a repeated upload of the same volume size
// reuses the same allocation. Everything here must run on the GL thread.

// Reserves a buffer of at least `size` bytes and returns its mapped range.
// Returns NULL if persistent mapping is unavailable or allocation failed, in
// which case the caller should use the plain client-pointer path.
void* pbo_upload_begin(size_t size);

// Uploads what was written into the pointer from pbo_upload_begin() into the
// whole of `texture`, then releases the buffer back to the cache. Must follow a
// successful pbo_upload_begin(); `size` must have covered the whole texture.
bool pbo_upload_end_texture_3D(GLuint texture, int level, GLenum format);

// Abandons an upload started with pbo_upload_begin() without touching a texture.
void pbo_upload_abort(void);

// Releases every cached buffer. Call before tearing down the GL context.
void pbo_upload_shutdown(void);

}  // namespace gl
