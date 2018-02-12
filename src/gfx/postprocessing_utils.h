#pragma once

#include <core/gl.h>
#include <core/types.h>

namespace postprocessing {

void initialize(int width, int height);
void shutdown();

void apply_ssao(GLuint depth_tex, const mat4& proj_mat, float radius = 3.f, float strength = 1.5f, float bias = 0.1f);
void apply_tonemapping(GLuint color_tex);

}