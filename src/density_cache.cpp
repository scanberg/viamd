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
#include <core/md_vec_math.h>

#include <mutex>
#include <string>
#include <vector>

DensityTiles::~DensityTiles() {
    md_gto_dmat_free(&dm);
}

struct DensityCache {
    struct Entry {
        std::string path;
        uint64_t    key = 0;
        std::shared_ptr<const DensityTiles> tiles;
    };
    std::mutex         mutex;
    std::vector<Entry> entries;
    task_system::Job   preload;     // density_cache_preload's; ~Job cancels and waits for it
};

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
}

std::shared_ptr<const DensityTiles> density_cache_tiles(DensityCache* cache, const md_system_t& sys, md_str_t path, const md_gto_basis_t* basis,
                                                         const float* atom_xyz, size_t atom_xyz_stride) {
    if (!basis || !atom_xyz) return {};
    if (!cache) return tiles_build(sys, path, basis, atom_xyz, atom_xyz_stride);

    const uint64_t key = tiles_key(sys, path);
    const std::string name(path.ptr, path.len);
    // held while it builds: another asking for it waits for this build rather than making its own
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t slot = cache->entries.size();
    for (size_t i = 0; i < cache->entries.size(); ++i) {
        if (cache->entries[i].path != name) continue;
        const DensityTiles* t = cache->entries[i].tiles.get();
        if (cache->entries[i].key == key && t && t->dm.num_shells == basis->num_shells &&
            md_gto_dmat_pattern_valid(&t->dm, atom_xyz, atom_xyz_stride)) {
            return cache->entries[i].tiles;
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
    if (slot == cache->entries.size()) cache->entries.push_back({});
    cache->entries[slot].path = name;
    cache->entries[slot].key = key;
    cache->entries[slot].tiles = t;
    return t;
}

void density_cache_preload(DensityCache* cache, const md_system_t& sys, md_str_t path, const md_float3_t* atom_xyz, size_t num_atoms) {
    if (!cache || !atom_xyz || num_atoms == 0) return;
    const md_system_t* s = &sys;
    // a static path or the caller's, and the positions: copied, so that the job holds its own
    std::string name(path.ptr, path.len);
    std::vector<float> xyz(3 * num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
        xyz[3 * i + 0] = atom_xyz[i].x;
        xyz[3 * i + 1] = atom_xyz[i].y;
        xyz[3 * i + 2] = atom_xyz[i].z;
    }
    cache->preload.start(MD_STR_LIT("Density matrix: building its tiles"), [cache, s, name, xyz = std::move(xyz), num_atoms]() {
        md_allocator_i* heap = md_get_heap_allocator();
        md_gto_basis_t basis = {};
        if (!md_gto_basis_extract_attributes(&basis, &s->attributes, heap)) return;
        if (md_gto_basis_num_atoms(&basis) <= num_atoms) {
            density_cache_tiles(cache, *s, md_str_t{name.data(), name.size()}, &basis, xyz.data(), 0);
        }
        md_gto_basis_free(&basis, heap);
        md_temp_trim();
    });
}

size_t density_cache_bytes(DensityCache* cache) {
    if (!cache) return 0;
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t bytes = 0;
    for (const DensityCache::Entry& e : cache->entries) {
        if (e.tiles) bytes += sizeof(double) * 256 * e.tiles->dm.num_tiles;
    }
    return bytes;
}
