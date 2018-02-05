#pragma once

#include <mol/molecule.h>

void draw_vdw(const Array<vec3> atom_positions, const Array<float> atom_radii, const Array<uint32> atom_colors, const mat4 view_mat, const mat4 proj_mat, float radii_scale = 1.f);
