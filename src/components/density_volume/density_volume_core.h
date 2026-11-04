#pragma once

// Comparing a volume with another made elsewhere: another run, another system (a mutant against its wild type).
//
// The other volume comes as a Gaussian cube file, as viamd exports one (or any cube with atoms): its grid, its values,
// and the reference structure it was made around, placed in its frame. The two reference structures need not be the
// same: which atoms of the other play the part of which of this one is found by their maximum common substructure
// (md_mcs, by element: a cube has no names), and the other is superposed on this one over those atoms (Kabsch, by
// mass). The other volume is then sampled at this one's voxels, so the two can be shown side by side or subtracted.
//
// Free of ImGui, OpenGL and ApplicationState: tested in tests/.

#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <stddef.h>
#include <stdint.h>

struct md_allocator_i;
struct md_system_t;

namespace density_volume {

// Values on a grid: voxel [i][j][k] spans origin + (i..i+1) basis[0] + (j..j+1) basis[1] + (k..k+1) basis[2], in Å;
// values[i][j][k], k fastest. As VIAL publishes a volume (/origin, /basis).
struct Grid {
    uint32_t     dim[3] = {};
    vec3_t       origin = {};
    vec3_t       basis[3] = {};
    const float* values = nullptr;
};

// A cube file as read: its grid in Å (a cube's origin is the first voxel's centre: here it is that voxel's corner),
// its atoms in the grid's frame, and the bulk density viamd writes in its first line, NaN when it has none
struct Cube {
    Grid            grid;
    float*          values = nullptr;
    size_t          num_atoms = 0;
    int*            atomic_numbers = nullptr;
    vec3_t*         xyz = nullptr;
    double          bulk = 0;
    md_allocator_i* alloc = nullptr;
};

// False, with a reason in err, when text is not a cube with at least one voxel
bool cube_parse(Cube* out, str_t text, md_allocator_i* alloc, char* err = nullptr, size_t err_cap = 0);
bool cube_read(Cube* out, str_t path, md_allocator_i* alloc, char* err = nullptr, size_t err_cap = 0);
void cube_free(Cube* cube);

// The other's frame onto this one: x_this = R x_other + t
struct Alignment {
    mat3_t R = mat3_ident();
    vec3_t t = {};
    size_t matched = 0;         // Heavy atoms paired
    size_t this_heavy = 0;      // Of this structure's and the other's heavy atoms
    size_t other_heavy = 0;
    double rmsd = 0;            // Å, of the pairs superposed
    bool   incomplete = false;  // The search gave up: matched is the most it found
};

// This structure, atoms of sys, at xyz: the positions of every atom of sys, in this volume's frame. False when the
// cube has no atoms or fewer than three heavy atoms are common, which places nothing.
bool align(Alignment* out, const md_system_t* sys, const int32_t* atoms, size_t count, const vec3_t* xyz, const Cube& other,
           md_allocator_i* alloc, char* err = nullptr, size_t err_cap = 0);

// from's values at the centres of to's voxels, through al (to's frame is this one's), trilinearly; 0 where from does
// not reach. Returns the fraction of to's voxels it reaches.
double resample(float* out, const Grid& to, const Grid& from, const Alignment& al);

}  // namespace density_volume
