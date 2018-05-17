#pragma once

#include <core/gl.h>
#include <core/types.h>
#include <core/volume.h>

namespace volume {

void initialize();
void shutdown();

GLuint create_volume_texture();
void free_volume_texture(GLuint texture);
void set_volume_texture_data(GLuint texture, const Volume& volume);

void render_volume_texture(GLuint volume_texture, GLuint depth_texture, const mat4& basis, const mat4& view_matrix, const mat4& proj_matrix,
                           vec3 color = vec3(1, 0, 0), float opacity_scale = 1.f);

}  // namespace volume
