#include "density_cache.h"
#include "task_system.h"

#include <md_system.h>
#include <md_attributes.h>
#include <md_qm.h>
#include <md_gto.h>
#include <md_gto_dmat.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <mutex>
#include <string>
#include <vector>

PackedDensity::~PackedDensity() {
    if (tri) md_free(md_get_heap_allocator(), tri, sizeof(double) * (dim * (dim + 1) / 2));
}

DensityTiles::~DensityTiles() {
    md_gto_dmat_free(&dm);
}

struct DensityCache {
    struct Entry {
        std::string path;
        uint64_t    key = 0;
        std::shared_ptr<const PackedDensity> density;
    };
    struct TilesEntry {
        std::string path;
        uint64_t    key = 0;
        std::shared_ptr<const DensityTiles> tiles;
    };
    std::mutex              mutex;
    std::vector<Entry>      entries;
    std::vector<TilesEntry> tiles;
    task_system::Job        preload;     // density_cache_preload's; ~Job cancels and waits for it
};

// What a density is made of: its attribute's version, and the orbitals' (an SCF density computed from them, as a
// file without the matrix publishes it, is a virtual attribute whose version need not follow theirs)
static uint64_t density_key(const md_system_t& sys, const md_attribute_t* attr) {
    static const md_str_t sources[4] = {
        MD_STR_INIT("orbital/alpha/coefficient"), MD_STR_INIT("orbital/alpha/occupation"),
        MD_STR_INIT("orbital/beta/coefficient"),  MD_STR_INIT("orbital/beta/occupation"),
    };
    uint64_t key = md_hash64(attr->path.ptr, attr->path.len, 0x7d1c4be3a6f50291ull);
    uint64_t v = md_attributes_version(&sys.attributes, attr->id);
    key = md_hash64(&v, sizeof(v), key);
    for (const md_str_t& p : sources) {
        const md_attribute_t* a = md_attributes_find(&sys.attributes, p);
        v = a ? md_attributes_version(&sys.attributes, a->id) : 0;
        key = md_hash64(&v, sizeof(v), key);
    }
    return key ? key : 1;
}

static std::shared_ptr<const PackedDensity> density_read(const md_system_t& sys, const md_attribute_t* attr, md_attribute_io_t* io) {
    const md_attribute_slice_t all = md_attribute_slice_all();
    const size_t n = md_qm_extract_packed_symmetric_f64(nullptr, 0, &sys.attributes, attr, all, io);
    if (n == 0) return {};
    const size_t len = n * (n + 1) / 2;
    auto d = std::make_shared<PackedDensity>();
    d->tri = (double*)md_alloc(md_get_heap_allocator(), sizeof(double) * len);
    if (!d->tri) {
        MD_LOG_ERROR("Density '" MD_STR_FMT "': no memory for its %.1f MB", MD_STR_ARG(attr->path), sizeof(double) * len / (1024.0 * 1024.0));
        return {};
    }
    d->dim = n;
    if (md_qm_extract_packed_symmetric_f64(d->tri, len, &sys.attributes, attr, all, io) != n) {
        return {};
    }
    return d;
}

DensityCache* density_cache_create() {
    return new DensityCache();
}

void density_cache_destroy(DensityCache* cache) {
    delete cache;
}

void density_cache_clear(DensityCache* cache) {
    if (!cache) return;
    std::lock_guard<std::mutex> lock(cache->mutex);
    cache->entries.clear();
    cache->tiles.clear();
}

// What the tiles of the density at path are made of: the orbitals' versions (md_qm_density_orbitals reads them)
static uint64_t tiles_key(const md_system_t& sys, md_str_t path) {
    static const md_str_t sources[4] = {
        MD_STR_INIT("orbital/alpha/coefficient"), MD_STR_INIT("orbital/alpha/occupation"),
        MD_STR_INIT("orbital/beta/coefficient"),  MD_STR_INIT("orbital/beta/occupation"),
    };
    uint64_t key = md_hash64(path.ptr, path.len, 0x3c6ef372fe94f82bull);
    for (const md_str_t& p : sources) {
        const md_attribute_t* a = md_attributes_find(&sys.attributes, p);
        const uint64_t v = a ? md_attributes_version(&sys.attributes, a->id) : 0;
        key = md_hash64(&v, sizeof(v), key);
    }
    return key ? key : 1;
}

static std::shared_ptr<const DensityTiles> tiles_build(const md_system_t& sys, md_str_t path, const md_gto_basis_t* basis, const float* atom_xyz,
                                                       size_t atom_xyz_stride) {
    md_allocator_i* heap = md_get_heap_allocator();
    double* rows = nullptr;
    double* w = nullptr;
    size_t num_ao = 0;
    const size_t num_mo = md_qm_density_orbitals(&rows, &w, &num_ao, &sys, path, heap);
    if (num_mo == 0) return {};
    MD_DEFER {
        md_free(heap, rows, sizeof(double) * num_mo * num_ao);
        md_free(heap, w, sizeof(double) * num_mo);
    };
    if (num_ao != md_gto_basis_num_ao(basis)) {
        MD_LOG_ERROR("Density '" MD_STR_FMT "': its orbitals span %zu atomic orbitals, the basis %zu", MD_STR_ARG(path), num_ao, md_gto_basis_num_ao(basis));
        return {};
    }
    auto t = std::make_shared<DensityTiles>();
    const md_gto_dmat_desc_t desc = { .basis = basis, .atom_xyz = atom_xyz, .atom_xyz_stride = atom_xyz_stride,
                                      .cutoff = DENSITY_TILES_CUTOFF, .margin = DENSITY_TILES_MARGIN };
    if (!md_gto_dmat_init(&t->dm, &desc, heap) || !md_gto_dmat_build(&t->dm, rows, num_ao, 1, w, num_mo, task_system::executor())) {
        return {};
    }
    return t;
}

std::shared_ptr<const DensityTiles> density_cache_tiles(DensityCache* cache, const md_system_t& sys, md_str_t path, const md_gto_basis_t* basis,
                                                         const float* atom_xyz, size_t atom_xyz_stride) {
    if (!basis || !atom_xyz) return {};
    if (!cache) return tiles_build(sys, path, basis, atom_xyz, atom_xyz_stride);

    const uint64_t key = tiles_key(sys, path);
    const std::string name(path.ptr, path.len);
    // held while it builds: another asking for it waits for this build rather than making its own
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t slot = cache->tiles.size();
    for (size_t i = 0; i < cache->tiles.size(); ++i) {
        if (cache->tiles[i].path != name) continue;
        const DensityTiles* t = cache->tiles[i].tiles.get();
        if (cache->tiles[i].key == key && t && t->dm.num_shells == basis->num_shells &&
            md_gto_dmat_pattern_valid(&t->dm, atom_xyz, atom_xyz_stride)) {
            return cache->tiles[i].tiles;
        }
        slot = i;   // of other orbitals, or of atoms that have moved too far: replaced
        break;
    }
    const md_tick_t t0 = md_tick_now();
    std::shared_ptr<const DensityTiles> t = tiles_build(sys, path, basis, atom_xyz, atom_xyz_stride);
    if (!t) return {};
    MD_LOG_DEBUG("Density '" MD_STR_FMT "': %llu of %llu tile pairs (%.1f MB in double) built from its orbitals in %.0f ms", MD_STR_ARG(path),
                 (unsigned long long)t->dm.num_tiles, (unsigned long long)t->dm.nt * (t->dm.nt + 1) / 2, 2048.0 * t->dm.num_tiles / (1024.0 * 1024.0),
                 md_tick_to_milliseconds(md_tick_now() - t0));
    if (slot == cache->tiles.size()) cache->tiles.push_back({});
    cache->tiles[slot].path = name;
    cache->tiles[slot].key = key;
    cache->tiles[slot].tiles = t;
    return t;
}

std::shared_ptr<const PackedDensity> density_cache_packed(DensityCache* cache, const md_system_t& sys, md_str_t path, md_attribute_io_t* io) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
    if (!attr) return {};
    if (!cache) return density_read(sys, attr, io);

    const uint64_t key = density_key(sys, attr);
    const std::string name(path.ptr, path.len);
    // held while it reads: another asking for it waits for this read rather than making its own
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t slot = cache->entries.size();
    for (size_t i = 0; i < cache->entries.size(); ++i) {
        if (cache->entries[i].path != name) continue;
        if (cache->entries[i].key == key) return cache->entries[i].density;
        slot = i;   // of other versions of what it is made of: replaced
        break;
    }
    const md_tick_t t0 = md_tick_now();
    std::shared_ptr<const PackedDensity> d = density_read(sys, attr, io);
    if (!d) return {};
    MD_LOG_DEBUG("Density '" MD_STR_FMT "': %zu atomic orbitals read in %.0f ms, kept as the packed triangle in double (%.1f MB)",
                 MD_STR_ARG(path), d->dim, md_tick_to_milliseconds(md_tick_now() - t0), sizeof(double) * (d->dim * (d->dim + 1) / 2) / (1024.0 * 1024.0));
    if (slot == cache->entries.size()) cache->entries.push_back({});
    cache->entries[slot].path = name;
    cache->entries[slot].key = key;
    cache->entries[slot].density = d;
    return d;
}

void density_cache_preload(DensityCache* cache, const md_system_t& sys, md_str_t path) {
    if (!cache || !md_attributes_find(&sys.attributes, path)) return;
    const md_system_t* s = &sys;
    // a static path or the caller's: copied, so that the job holds its own
    std::string name(path.ptr, path.len);
    cache->preload.start(MD_STR_LIT("Density matrix: reading"), [cache, s, name]() {
        // through the pool's executor: a density computed from its orbitals (a file without the matrix) is spread
        // over it (task_system: a job's loops run in it and stop with it)
        md_attribute_io_t* io = md_attribute_io_create(task_system::executor(), md_get_heap_allocator());
        density_cache_packed(cache, *s, md_str_t{name.data(), name.size()}, io);
        md_attribute_io_destroy(io);
        md_temp_trim();
    });
}

size_t density_cache_bytes(DensityCache* cache) {
    if (!cache) return 0;
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t bytes = 0;
    for (const DensityCache::Entry& e : cache->entries) {
        if (e.density) bytes += sizeof(double) * (e.density->dim * (e.density->dim + 1) / 2);
    }
    for (const DensityCache::TilesEntry& e : cache->tiles) {
        if (e.tiles) bytes += sizeof(double) * 256 * e.tiles->dm.num_tiles;
    }
    return bytes;
}
