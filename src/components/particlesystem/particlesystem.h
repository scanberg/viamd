#pragma once

#include <core/md_vec_math.h>
#include <gfx/gl.h>

namespace particlesystem {

void initialize();
void shutdown();
void set_volume_texture(GLuint texture, int dim_x, int dim_y, int dim_z, vec3_t min, vec3_t max);
void draw_ui();

}  // namespace particlesystem
