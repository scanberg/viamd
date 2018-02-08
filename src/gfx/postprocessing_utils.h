#pragma once

#include <core/gl.h>

namespace postprocessing {

void initialize();
void shutdown();

void ssao(GLuint depth_tex, float strength);
void tonemapping(GLuint color_tex);

}