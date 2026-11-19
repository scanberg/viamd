#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_grid.h>
#include <core/md_vec_math.h>
#include <md_attributes.h>

struct md_system_t;
struct md_system_state_t;

// GRID FIELDS
//
// A field given by a file as its values on a grid: a cube file, which md_cube publishes as cube/<name> (md_cube.h):
// the values, {x, y, z} or {m, x, y, z} for m of them a point, x outermost; the corner of the first voxel and a step
// along each axis (origin, basis, Å); and the cube's atoms in its frame (atom/coordinate, Å), with which of the system's
// they are (atom/system_index). The grid is placed where those atoms are now: the rigid motion (Kabsch) of its atoms
// onto theirs, which is none for the frame the cube was made of, and follows the molecule through a trajectory. The
// values between the points are trilinear, and 0 outside the grid.

// The group the grids are published under, and whether an attribute is one of them (a values attribute, not one of its
// siblings)
#define GRID_FIELD_PREFIX "cube"
bool grid_field_is_grid(const md_attribute_t* attr);

// The first grid of the system, MD_ATTRIBUTE_INVALID for none
md_attribute_id_t grid_field_first(const md_system_t& sys);

// The values of a grid: how many a point (1 for {x, y, z}), 0 for an attribute that is not a grid
uint32_t grid_field_count(const md_system_t& sys, md_attribute_id_t key);

// The id of the member-th value of a point (an MO cube's orbital number), -1 where the cube gives none
int grid_field_member_id(const md_system_t& sys, md_attribute_id_t key, uint32_t member);

struct GridField {
    float*  values = nullptr;       // the member's, dim[0] * dim[1] * dim[2], z fastest; heap, grid_field_free
    int     dim[3] = {};
    // The grid where it is now, Å: the centre of voxel (i, j, k) at corner + axes * (i + 1/2, j + 1/2, k + 1/2)
    md_float3x3_t axes = {};        // columns: a step along x, y and z
    md_float3_t   corner = {};
    md_float3x3_t world_to_index = {};  // axes^-1
    double rmsd = 0.0;              // Å, of its atoms on the system's, as placed
    size_t num_placed_by = 0;       // the atoms it was placed by: 0 where it has none (placed as it is)
};

// The member-th values of the grid at key, placed on the atoms of state (without values: where it is alone). False for
// no such grid.
bool grid_field_get(GridField* out, const md_system_t& sys, const md_system_state_t& state, md_attribute_id_t key, uint32_t member, bool values = true);
void grid_field_free(GridField* f);

// The value at a point, bohr (where the fields of surface_field.h are evaluated)
float grid_field_sample_bohr(const GridField& f, const float xyz[3]);

// The grid itself as an md_grid_t (bohr), its voxels the field's; where its axes are skewed, the box of it along the
// world's axes at its finest step. False for a grid of no extent.
bool grid_field_md_grid(md_grid_t* grid, const GridField& f);
