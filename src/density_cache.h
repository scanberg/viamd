#pragma once

// The SCF density matrices the CPU paths read - the ESP's charge distribution (surface_field.cpp), the critical
// point search (quantum_chemistry.cpp) - kept once read, each as its packed upper triangle in double: row i holding
// columns i..N-1, the symmetric part (md_qm_extract_packed_symmetric_f64; md_gto_int takes it as density_packed).
// Reading one from its file is seconds on a protein (a spherical matrix read whole and transformed), and each of
// them read it again whenever it was wanted; kept, it is read once, and is half the square's size (371 MB for 9636
// AOs). The GPU paths do not read it: they build their tiles from the orbitals (md_gto_gpu_dmat).
//
// Read when a system is loaded (density_cache_preload, from init_system_data), on the task pool while the geometry is
// shown, so that the first of them finds it there.
//
// Of a system: cleared with it (free_system_data, once its tasks have ended). The cache itself lives as long as the
// application, so that a task holding a pointer to it is never left with a freed one across a load.

#include <core/md_str.h>

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

DensityCache* density_cache_create();
void density_cache_destroy(DensityCache* cache);
void density_cache_clear(DensityCache* cache);

// The density at path, read the first time it is asked for and kept while the attributes it is made of keep their
// versions (a newer one replaces it; who holds the old one keeps it until they let it go). Read through io, which
// may be null (then every loop of the read on the calling thread). Any thread: one that asks while another reads
// waits for it. Null when the system has no such symmetric density, or it cannot be read (logged by the reader).
// cache null: read, not kept.
std::shared_ptr<const PackedDensity> density_cache_packed(DensityCache* cache, const md_system_t& sys, md_str_t path, md_attribute_io_t* io);

// Reads the density at path into the cache on the task pool, now (a pool job of the cache's, listed while it runs):
// a consumer asking for it meanwhile waits for this read rather than making its own. A load cancels it before it has
// begun (task_system::cancel_all); one under way is waited for, as a read cannot stop halfway. sys must stay until
// it has ended: the application's, which a load frees only after waiting for the pool.
void density_cache_preload(DensityCache* cache, const md_system_t& sys, md_str_t path);

// The bytes kept, for a display or a log
size_t density_cache_bytes(DensityCache* cache);
