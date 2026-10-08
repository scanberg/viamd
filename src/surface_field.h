#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_grid.h>
#include <core/md_unit.h>

#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;

// SURFACE FIELDS
//
// A FIELD is a scalar function over R^3 with an identity and a unit. A SURFACE is whatever a
// representation draws. A surface field colours the one by the other, and keeps three things apart
// that used to be one:
//
//   what      the field, evaluated by mdlib (md_gto_int for potentials)
//   where     a world space volume of field values, evaluated ONLY in the narrow band of voxels a
//             surface can sample: here the corners of the cells an isosurface passes through, of a
//             volume the field volume shares its grid with. Cost goes with the surface's area,
//             not the volume's, and a band that moves (a new isovalue) only evaluates what it has
//             not evaluated before
//   how       a colour map and a range, applied per pixel when the surface is shaded - changing
//             either costs nothing, and the values interpolate as values rather than as colours
//
// Values outside the band are never sampled by a surface of that band and are left at 0.
//
// UNITS: positions in bohr (the grid's), values in the field's unit (surface_field_unit).

enum class SurfaceFieldKind {
    // The electrostatic potential of the classical multipoles of the system: the charges
    // (atom/charge), and the dipoles and quadrupoles (atom/dipole, atom/quadrupole) where the
    // potential has them, with the atoms that carry none - the QM atoms of an embedding, NAN there -
    // left out. For a polarizable embedding this is its PERMANENT part only: the induced dipoles are
    // solved for during the calculation and not stored, so they are not in it.
    EmbeddingPotential,
    Count
};

inline const char* surface_field_kind_str[(int)SurfaceFieldKind::Count] = {
    "Embedding Potential",
};

// How a field is mapped to colour starts out as fits a potential: a range centred on zero, following
// the field's spread on the surface, on a diverging map - red for negative, blue for positive.
static inline ColorScale surface_field_default_scale() {
    ColorScale scale;
    scale.colormap   = COLOR_SCALE_COLORMAP_RDBU;
    scale.range_beg  = -0.05f;      // in the field's unit
    scale.range_end  =  0.05f;
    scale.symmetric  = true;
    scale.auto_range = true;
    return scale;
}

struct SurfaceFieldVolume {
    uint32_t  tex_id = 0;           // GL_R32F, the dimensions of the grid
    uint32_t  colormap_tex = 0;     // RGBA8, 256 x 1
    int       colormap = -1;        // which colour map colormap_tex holds
    int       dim[3] = {};

    float*    values = nullptr;     // [num_voxels] x fastest, the field where evaluated, 0 elsewhere
    uint64_t* done = nullptr;       // [num_voxels / 64] which voxels hold a value
    size_t    num_voxels = 0;
    size_t    num_evaluated = 0;

    uint64_t  source_hash = 0;      // what the values are OF; a change invalidates every one of them
    uint64_t  band_hash = 0;        // the isovalues the band was taken for

    // The spread of the field ON the surface: its values interpolated to where the isosurfaces
    // cross the edges of the grid. lo / hi are the 1st and 99th percentiles - a classical charge
    // can sit within half an Angstrom of a QM surface, and the spike it makes there would otherwise
    // be the range.
    size_t    num_surface_samples = 0;
    float     surface_min = 0.0f;
    float     surface_max = 0.0f;
    float     surface_lo = 0.0f;
    float     surface_hi = 0.0f;
};

md_unit_t surface_field_unit(SurfaceFieldKind kind);

// Whether the system has what the field is made from
bool surface_field_available(SurfaceFieldKind kind, const md_system_t& sys);

// Brings the field volume up to date for the isosurfaces at iso_values of `density`, sampled on
// `grid` (bohr, x fastest, voxel centres at (i + 0.5) * spacing). source_hash names everything the
// field depends on (geometry, charges, the grid); band_hash the isovalues. Returns false when the
// field cannot be made, and leaves the volume empty then.
bool surface_field_update(SurfaceFieldVolume* vol, SurfaceFieldKind kind, const md_system_t& sys, const md_system_state_t& state,
                          const md_grid_t& grid, const float* density, const float* iso_values, size_t num_iso,
                          uint64_t source_hash, uint64_t band_hash);

// The spread of the field on the surface, as a colour scale reads it: the 1st to 99th percentile for
// an automatic range, the extremes beside it. Empty while there are no samples.
ColorScaleSpan surface_field_span(const SurfaceFieldVolume& vol);

// The colour map texture for the mapping's colour map, rebuilt when it changed
uint32_t surface_field_colormap_texture(SurfaceFieldVolume* vol, int colormap);

void surface_field_free(SurfaceFieldVolume* vol);
