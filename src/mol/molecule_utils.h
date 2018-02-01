#pragma once

void transform_molecule(Array<vec3> atom_pos, const mat4 transformation);
void compute_bonds(Array<vec3> atom_bonds, const Array<vec3> atom_pos, const Array<Element> atom_elem);
