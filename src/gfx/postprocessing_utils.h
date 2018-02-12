#pragma once

#include <core/gl.h>

namespace postprocessing {

void initialize(int width, int height);
void shutdown();

void apply_ssao(GLuint depth_tex, float strength);
void apply_tonemapping(GLuint color_tex);

}