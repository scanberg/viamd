#pragma once

// The SCF density matrices the CPU paths read, kept once made, in two forms while the consumers move from one to
// the other:
//
//   tiles    (density_cache_tiles) the one representation of a density: md_gto_dmat_t, in double, built from the
//            orbitals at DENSITY_TILES_CUTOFF and DENSITY_TILES_MARGIN, never the whole matrix. The ESP's charge
//            distribution (surface_field.cpp) reads it. 39 MB for a 1152 atom protein (9636 AOs), built in about a
//            second on 8 threads.
//   packed   (density_cache_packed) the packed upper triangle in double, row i holding columns i..N-1, the symmetric
//            part (md_qm_extract_packed_symmetric_f64), read from the file: what the critical point search
//            (quantum_chemistry.cpp) still reads. Seconds on a protein, 371 MB.
//
// The GPU paths do not read either: they build their tiles from the orbitals (md_gto_gpu_dmat).
//
// Read when a system is loaded (density_cache_preload, from init_system_data), on the task pool while the geometry is
// shown, so that the first of them finds it there.
//
// Of a system: cleared with it (free_system_data, once its tasks have ended). The cache itself lives as long as the
// application, so that a task holding a pointer to it is never left with a freed one across a load.

#include <core/md_str.h>
#include <md_gto_dmat.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>

struct md_system_t;
struct md_attribute_io_t;
struct DensityCache;

struct PackedDensity {
    size_t  dim = 0;            // N
    double* tri = nullptr;      // [N (N + 1) / 2]
    PackedDensity() = default;
    PackedDensity(const PackedDensity&) = delete;
    PackedDensity& operator=(const PackedDensity&) = delete;
    ~PackedDensity();
};

// The one tile representation's cutoff and margin (bohr), for every consumer: at 1e-6 ELF and LOL are within 1.1e-4
// of a negligible cutoff where rho > 1e-3, and the critical points of a 1152 atom protein and its ESP the same as of
// its whole matrix (md_dmat_cutoff_bench, md_topo_gto_bench --tiles, md_dmat_census_bench). The margin keeps the
// pattern while no atom has moved more than half of it.
#define DENSITY_TILES_CUTOFF 1.0e-6
#define DENSITY_TILES_MARGIN 1.0

struct DensityTiles {
    md_gto_dmat_t dm = {};
    DensityTiles() = default;
    DensityTiles(const DensityTiles&) = delete;
    DensityTiles& operator=(const DensityTiles&) = delete;
    ~DensityTiles();
};

DensityCache* density_cache_create();
void density_cache_destroy(DensityCache* cache);
void density_cache_clear(DensityCache* cache);

// The density at path, read the first time it is asked for and kept while the attributes it is made of keep their
// versions (a newer one replaces it; who holds the old one keeps it until they let it go). Read through io, which
// may be null (then every loop of the read on the calling thread). Any thread: one that asks while another reads
// waits for it. Null when the system has no such symmetric density, or it cannot be read (logged by the reader).
// cache null: read, not kept.
std::shared_ptr<const PackedDensity> density_cache_packed(DensityCache* cache, const md_system_t& sys, md_str_t path, md_attribute_io_t* io);

// The density at path as tiles over basis, made from its orbitals (md_qm_density_orbitals) the first time it is asked
// for, and again when the orbitals' versions change or an atom has moved more than DENSITY_TILES_MARGIN / 2 from where
// the pattern was made (atom_xyz: bohr, one per atom of the basis, stride bytes apart, 0 packed). Built through the
// task pool's executor. Any thread: one that asks while another builds waits for it. Null when the system has no
// orbitals for that density. cache null: built, not kept.
std::shared_ptr<const DensityTiles> density_cache_tiles(DensityCache* cache, const md_system_t& sys, md_str_t path, const md_gto_basis_t* basis,
                                                         const float* atom_xyz, size_t atom_xyz_stride);

// Reads the density at path into the cache on the task pool, now (a pool job of the cache's, listed while it runs):
// a consumer asking for it meanwhile waits for this read rather than making its own. A load cancels it before it has
// begun (task_system::cancel_all); one under way is waited for, as a read cannot stop halfway. sys must stay until
// it has ended: the application's, which a load frees only after waiting for the pool.
void density_cache_preload(DensityCache* cache, const md_system_t& sys, md_str_t path);

// The bytes kept (both forms), for a display or a log
size_t density_cache_bytes(DensityCache* cache);
