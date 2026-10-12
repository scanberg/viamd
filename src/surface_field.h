#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_grid.h>
#include <core/md_unit.h>
#include <core/md_gpu.h>

#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;
struct md_gto_int_charges_t;

// SURFACE FIELDS
//
// A FIELD is a scalar function over R^3 with an identity and a unit. A SURFACE is whatever a
// representation draws. A surface field colours the one by the other, and keeps three things apart
// that used to be one:
//
//   what      the field, evaluated by mdlib (md_gto_int for potentials)
//   where     a world space volume of field values over the box of the volume the surfaces are drawn
//             from, sampled by the same normalised texture coordinates. On the CPU it is evaluated
//             ONLY in the narrow band of voxels a surface can sample: the voxels around the cells an
//             isosurface passes through. Cost goes with the surface's area, not the volume's, and a
//             band that moves (a new isovalue, another density) only evaluates what it has not
//             evaluated before. On the GPU the whole volume is evaluated once, after which no surface
//             of it costs anything
//   how       a colour map and a range, applied per pixel when the surface is shaded - changing
//             either costs nothing, and the values interpolate as values rather than as colours
//
// The field volume has the density volume's grid, or a coarser one over the same box where the
// field is smooth and dear: the electrostatic potential of a QM density is a sum over tens of
// thousands of gaussians per point, and on an isodensity surface it varies on the scale of the
// distance to the nuclei. See surface_field_grid.
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

    // The electrostatic potential of the QM region itself, its nuclei and its ground state electron
    // density: V(r) = sum_A Z_A / |r - A| - sum_{mu,nu} D_{mu nu} (mu| 1/|r' - r| |nu), the molecular
    // electrostatic potential that is mapped onto an isodensity surface. D is the total SCF density
    // (orbital/total/density, or the alpha density where there is no beta channel to add), Z the
    // nuclear charges (qm/atom/nuclear_charge): the effective ones where the file states them, else
    // the atomic numbers - too much charge by the core electrons of an effective core potential the
    // file leaves unstated, which is checked for and reported.
    ElectrostaticPotential,
    Count
};

inline const char* surface_field_kind_str[(int)SurfaceFieldKind::Count] = {
    "Embedding Potential",
    "Electrostatic Potential",
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
    uint32_t  tex_id = 0;           // GL_R32F, the dimensions of grid
    uint32_t  colormap_tex = 0;     // RGBA8, 256 x 1
    int       colormap = -1;        // which colour map colormap_tex holds
    md_grid_t grid = {};            // what the values are sampled on: see surface_field_grid

    float*    values = nullptr;     // [num_voxels] x fastest, the field where evaluated, 0 elsewhere
    uint64_t* done = nullptr;       // [num_voxels / 64] which voxels hold a value
    size_t    num_voxels = 0;
    size_t    num_evaluated = 0;
    size_t    num_uploaded = 0;     // num_evaluated when the texture was last written
    bool      complete = false;     // every voxel holds a value (the GPU evaluates the whole volume)
    bool      on_gpu = false;       // and it was the GPU that evaluated it

    // The charge distribution the CPU evaluates, kept while the band grows over the same source
    md_gto_int_charges_t* charges = nullptr;

    uint64_t  source_hash = 0;      // what the values are OF; a change invalidates every one of them
    uint64_t  band_hash = 0;        // the surfaces the band was taken for

    // The spread of the field ON the surface: its values interpolated to where the isosurfaces
    // cross the edges of the density grid. lo / hi are the 1st and 99th percentiles - a classical
    // charge can sit within half an Angstrom of a QM surface, and the spike it makes there would
    // otherwise be the range.
    size_t    num_surface_samples = 0;
    float     surface_min = 0.0f;
    float     surface_max = 0.0f;
    float     surface_lo = 0.0f;
    float     surface_hi = 0.0f;
};

// Everything one update needs. The geometry is the caller's to choose, as everywhere an evaluation
// is made: the positions here are what the field is evaluated at, whatever frame they came from.
struct SurfaceFieldDesc {
    SurfaceFieldKind         kind = SurfaceFieldKind::EmbeddingPotential;
    const md_system_t*       sys = nullptr;
    const md_system_state_t* state = nullptr;           // the environment's sites (embedding potential)
    const vec3_t*            basis_atom_xyz = nullptr;  // bohr, one per basis atom: the QM nuclei and where
    size_t                   num_basis_atoms = 0;       // the basis sits (electrostatic potential)

    // The surfaces: the isosurfaces at iso_values of `density`, sampled on `grid` (bohr, x fastest,
    // voxel centres at (i + 0.5) * spacing)
    const md_grid_t*         grid = nullptr;
    const float*             density = nullptr;
    const float*             iso_values = nullptr;
    size_t                   num_iso = 0;

    // source_hash names everything the VALUES depend on (the kind, the geometry, the data);
    // band_hash everything the SURFACES depend on (the density and the isovalues). A new source
    // starts over; a new band keeps every value it already has.
    uint64_t                 source_hash = 0;
    uint64_t                 band_hash = 0;

    // Optional: evaluate on the device, into a 3D R32F storage texture of at least scratch_dim, and
    // read the result back. Used where it can serve the field (the whole field grid fits the scratch
    // and md_gto_int evaluates the distribution there), else the band is evaluated on the CPU.
    md_gpu_stream_t          gpu_stream = nullptr;
    md_gpu_texture_t         gpu_scratch = nullptr;
    uint32_t                 gpu_scratch_dim[3] = {};
};

md_unit_t surface_field_unit(SurfaceFieldKind kind);

// Whether the system has what the field is made from. Cheap: the attributes are looked up, not read.
bool surface_field_available(SurfaceFieldKind kind, const md_system_t& sys);

// The first kind the system has what to make of, Count for none
SurfaceFieldKind surface_field_first_available(const md_system_t& sys);

// The grid a field is sampled on, given the grid of the volume its surfaces are drawn from: that
// grid, or for the electrostatic potential the same box at about 0.4 bohr (never finer than the
// density's). Trilinear interpolation from 0.4 bohr is within 0.3-0.4% of the potential's spread
// over an isodensity surface (amide, 26 atom mol.h5), under one step of a 256 entry colour map,
// and the coarser grid is a fraction of the evaluation: for mol.h5 its band is ~3 times smaller than
// on the density's default resolution and ~12 times than on its high one, the whole grid 5 and 40
// times.
md_grid_t surface_field_grid(SurfaceFieldKind kind, const md_grid_t& density_grid);

// Brings the field volume up to date for the surfaces of desc. Returns false when the field cannot
// be made, and leaves the volume empty then.
bool surface_field_update(SurfaceFieldVolume* vol, const SurfaceFieldDesc& desc);

// Forgets what the values are of, keeping the allocations: the next update evaluates afresh. For
// when the data changes under the same geometry, a new system loaded over the old.
void surface_field_invalidate(SurfaceFieldVolume* vol);

// The spread of the field on the surface, as a colour scale reads it: the 1st to 99th percentile for
// an automatic range, the extremes beside it. Empty while there are no samples.
ColorScaleSpan surface_field_span(const SurfaceFieldVolume& vol);

// The colour map texture for the mapping's colour map, rebuilt when it changed
uint32_t surface_field_colormap_texture(SurfaceFieldVolume* vol, int colormap);

void surface_field_free(SurfaceFieldVolume* vol);
