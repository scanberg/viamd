#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_grid.h>
#include <core/md_unit.h>
#include <core/md_gpu.h>

#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;
struct SurfaceFieldWork;
struct SurfaceFieldPointsWork;

// SURFACE FIELDS
//
// A FIELD is a scalar function over R^3 with an identity and a unit. A SURFACE is whatever a
// representation draws. A surface field colours the one by the other, and keeps three things apart
// that used to be one:
//
//   what      the field, evaluated by mdlib (md_gto_int for potentials, md_gto_iso for what is made of
//             orbitals)
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
// NOTHING HERE BLOCKS. Making what the field is evaluated from (reading a QM density and building
// md_gto_int's distribution of it: seconds on a protein; the ionization energy's orbitals) runs on a worker thread, the evaluation on
// the GPU as an md_gto_int job advanced a little each frame, and its readback is waited for by
// looking. An update works out what is wanted and starts it; surface_field_poll, once a frame while
// surface_field_busy, takes what has finished. Until new values arrive a surface keeps the ones it
// has, and surface_field_progress names what is under way for a progress display.
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
    // left out. For a polarizable embedding whose run kept its induced dipoles (atom/induced_dipole)
    // those are part of it, added to the sites' dipoles: the total potential the QM region was
    // embedded in. Without them it is the permanent part only - see
    // surface_field_embedding_content for telling the two apart.
    EmbeddingPotential,

    // The electrostatic potential of the QM region itself, its nuclei and its ground state electron
    // density: V(r) = sum_A Z_A / |r - A| - sum_{mu,nu} D_{mu nu} (mu| 1/|r' - r| |nu), the molecular
    // electrostatic potential that is mapped onto an isodensity surface. D is the total SCF density
    // (orbital/total/density, or the alpha density where there is no beta channel to add), Z the
    // nuclear charges (qm/atom/nuclear_charge): the effective ones where the file states them, else
    // the atomic numbers - too much charge by the core electrons of an effective core potential the
    // file leaves unstated, which is checked for and reported.
    ElectrostaticPotential,

    // The average local ionization energy of the QM region (Sjoberg, Murray, Brinck and Politzer 1990),
    // I(r) = sum_k n_k |e_k| |phi_k(r)|^2 / sum_k n_k |phi_k(r)|^2 over its occupied orbitals: the energy
    // with which the electrons at r are held, on average. Mapped onto an isodensity surface, its lowest
    // values are where the least tightly held electrons are - the sites an electrophile or a radical
    // attacks - and it is never below |e_HOMO|. Hartree. Made of the canonical orbitals, their occupations
    // and their energies (md_qm_alie_orbitals); not for natural orbitals, and not for a file without
    // orbital energies. Kohn-Sham orbital energies give lower values than Hartree-Fock ones, as usual.
    AverageLocalIonizationEnergy,
    Count
};

inline const char* surface_field_kind_str[(int)SurfaceFieldKind::Count] = {
    "Embedding Potential",
    "Electrostatic Potential",
    "Local Ionization Energy",
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

// What a scale should be for a kind of field, when the kind is chosen: a potential's range centred on zero;
// the ionization energy's from its lowest values on the surface (red: the least tightly held electrons) to
// its highest (blue), which a range about zero would squeeze into its upper half. The colour map stays the
// user's.
static inline void surface_field_scale_for_kind(ColorScale* scale, SurfaceFieldKind kind) {
    if (!scale) return;
    scale->symmetric  = kind != SurfaceFieldKind::AverageLocalIonizationEnergy;
    scale->auto_range = true;
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
    bool      complete = false;     // every voxel holds a value
    bool      on_gpu = false;       // and it was the GPU that evaluated the last of them

    // What is under way (surface_field.cpp): the charge distribution, kept while the band grows over
    // the same source, and the evaluation
    SurfaceFieldWork* work = nullptr;

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
    size_t                   num_basis_atoms = 0;       // the basis sits (surface_field_uses_basis)

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

    // Optional: evaluate on the device, where mdlib can evaluate the field there (point charges and the QM
    // density by md_gto_int, the orbitals by md_gto_iso), else on the CPU's worker threads. The stream must
    // outlive the volume.
    md_gpu_stream_t          gpu_stream = nullptr;
};

md_unit_t surface_field_unit(SurfaceFieldKind kind);

// Whether the field is made of the QM region's basis, at basis_atom_xyz (SurfaceFieldDesc)
static inline bool surface_field_uses_basis(SurfaceFieldKind kind) {
    return kind == SurfaceFieldKind::ElectrostaticPotential || kind == SurfaceFieldKind::AverageLocalIonizationEnergy;
}

// Whether the system has what the field is made from. Cheap: the attributes are looked up, and for the
// ionization energy the orbitals' occupations read (whether they are canonical).
bool surface_field_available(SurfaceFieldKind kind, const md_system_t& sys);

// What the Embedding Potential of a dataset is made of: how many of the environment's sites carry
// each source (a finite value, not all zero). 'polarizable' counts the sites with a polarizability
// above zero; with none of them 'induced', the dataset is a polarizable embedding that did not keep
// its induced dipoles - VeloxChem writes them only from its fika embedding on - and the potential is
// its permanent part only, short of what the QM region was embedded in. Reads a column over every
// atom per source: for when it is shown, not every frame.
struct SurfaceFieldEmbeddingContent {
    size_t charges = 0;
    size_t dipoles = 0;       // permanent
    size_t quadrupoles = 0;
    size_t induced = 0;       // induced dipoles
    size_t polarizable = 0;
};
SurfaceFieldEmbeddingContent surface_field_embedding_content(const md_system_t& sys);

// The first kind the system has what to make of, Count for none
SurfaceFieldKind surface_field_first_available(const md_system_t& sys);

// The grid a field is sampled on, given the grid of the volume its surfaces are drawn from: that
// grid (the ionization energy: a point of it costs little), or for the electrostatic potential the same box at about 0.4 bohr (never finer than the
// density's). Trilinear interpolation from 0.4 bohr is within 0.3-0.4% of the potential's spread
// over an isodensity surface (amide, 26 atom mol.h5), under one step of a 256 entry colour map,
// and the coarser grid is a fraction of the evaluation: for mol.h5 its band is ~3 times smaller than
// on the density's default resolution and ~12 times than on its high one, the whole grid 5 and 40
// times.
md_grid_t surface_field_grid(SurfaceFieldKind kind, const md_grid_t& density_grid);

// Brings the field volume up to date for the surfaces of desc, without waiting: works out the band,
// and starts what it needs (or leaves it to the evaluation under way to finish first, then the rest).
// A new source drops what is under way; the texture keeps the old values until the new ones come.
// Returns false when the field cannot be made: no such field in the system, or its distribution
// failed to build (not tried again for the same source).
bool surface_field_update(SurfaceFieldVolume* vol, const SurfaceFieldDesc& desc);

// Takes what has finished and starts what is next; true when the texture and the span changed. Once a
// frame while surface_field_busy, from the GL thread.
bool surface_field_poll(SurfaceFieldVolume* vol);
bool surface_field_busy(const SurfaceFieldVolume& vol);

// What is under way for longer than a moment, every volume and point set, for a progress display:
// the field's name and how far it is (0 to 1; 0 while its distribution is being made). Returns how
// many (at most cap).
struct SurfaceFieldProgress {
    const char* label = "";
    float       fraction = 0.0f;
    bool        preparing = false;  // the distribution is being made
};
size_t surface_field_progress(SurfaceFieldProgress* out, size_t cap);

// Forgets what the values are of, keeping the allocations: the next update evaluates afresh. For
// when the data changes under the same geometry, a new system loaded over the old.
void surface_field_invalidate(SurfaceFieldVolume* vol);

// The spread of the field on the surface, as a colour scale reads it: the 1st to 99th percentile for
// an automatic range, the extremes beside it. Empty while there are no samples.
ColorScaleSpan surface_field_span(const SurfaceFieldVolume& vol);

// The colour map texture for the mapping's colour map, rebuilt when it changed
uint32_t surface_field_colormap_texture(SurfaceFieldVolume* vol, int colormap);

void surface_field_free(SurfaceFieldVolume* vol);

// FIELD AT POINTS
// The field evaluated exactly where it is wanted, at a list of points - the vertices of an isosurface
// mesh (iso_meshes.h) - rather than on a volume a surface samples: no grid, no band, no interpolation
// but across a triangle. The charge distribution is kept while desc.source_hash stays the same, so a
// new set of points (another isovalue) costs only their evaluation. On the GPU when desc.gpu_stream is
// set and md_gto_int can evaluate the distribution there (point charges, the QM density), else on the
// CPU's worker threads. Without waiting, as the volume: begin, then poll once a frame until done.
struct SurfaceFieldPoints {
    SurfaceFieldPointsWork* work = nullptr;
    bool     on_gpu = false;            // where the last evaluation ran
    double   seconds = 0.0;             // and what it took, from begin to the values
};

// The points, bohr. The GPU path reads gpu_xyz when there is one (points already on the device, say where
// an extraction left them: no upload), else uploads xyz; the CPU path reads xyz. Either alone will do,
// xyz is what lets the CPU take over from a GPU that cannot evaluate the field. On the GPU the ionization
// energy is evaluated where the points are (md_gto_iso_gpu_evaluate_launch), the potentials by md_gto_int.
struct SurfaceFieldPointSet {
    const float*  xyz = nullptr;        // host, packed float3
    md_gpu_addr_t gpu_xyz = 0;          // device, float xyz at gpu_xyz + i * gpu_stride
    uint32_t      gpu_stride = 0;       // bytes, a multiple of 4; 0 = packed
    size_t        count = 0;
};

// Starts evaluating the field at the points, which are copied (on the host, and on the device when
// they are there), replacing an evaluation under way. desc: kind, sys, state, basis_atom_xyz /
// num_basis_atoms, source_hash, gpu_stream (optional, must outlive fp); the rest is not used. Returns
// false when the field cannot be made, or cannot be evaluated where the points are.
bool surface_field_points_begin(SurfaceFieldPoints* fp, const SurfaceFieldPointSet& points, const SurfaceFieldDesc& desc);

enum class SurfaceFieldPointsStatus { None, Running, Done, Failed };

// Advances the evaluation begun last. Done: *values (count of them, owned by fp, valid until the next
// begin or free) holds the field at the points, in their order; Done is returned once.
SurfaceFieldPointsStatus surface_field_points_poll(SurfaceFieldPoints* fp, const float** values, size_t* count);

void surface_field_points_free(SurfaceFieldPoints* fp);

// The spread of values (1st to 99th percentile, and the extremes) as a colour scale reads it
ColorScaleSpan surface_field_span_of(const float* values, size_t count);
