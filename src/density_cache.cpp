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

#if MD_ENABLE_GPU
DensityGpuTiles::~DensityGpuTiles() {
    if (dm) md_gto_gpu_dmat_destroy(stream, dm);
}
#endif

struct DensityCache {
    struct Entry {
        std::string path;
        md_attribute_slice_t slice = {};
        uint64_t    key = 0;
        uint64_t    used = 0;       // when it was last handed out (DensityCache::clock)
        std::shared_ptr<const DensityTiles> tiles;
    };
    std::mutex         mutex;
    std::vector<Entry> entries;
    uint64_t           clock = 0;
    task_system::Job   preload;     // density_cache_preload's; ~Job cancels and waits for it
};

static md_attribute_slice_t slice_of(const md_attribute_slice_t* slice) {
    return slice ? *slice : md_attribute_slice_all();
}

static bool slice_eq(const md_attribute_slice_t& a, const md_attribute_slice_t& b) {
    if (a.num_idx != b.num_idx) return false;
    for (uint32_t i = 0; i < a.num_idx; ++i) if (a.idx[i] != b.idx[i]) return false;
    return true;
}

uint64_t density_tiles_key(const md_system_t& sys, md_str_t path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_of(slice_ptr);
    uint64_t key = md_hash64(path.ptr, path.len, 0x3c6ef372fe94f82bull);
    key = md_hash64(&slice.num_idx, sizeof(slice.num_idx), key);
    if (slice.num_idx) key = md_hash64(slice.idx, sizeof(uint32_t) * slice.num_idx, key);
    // what md_qm_density_factors reads: the factors published beside the density, or an SCF density's orbitals
    std::string sources[4];
    if (md_qm_density_factors_available(&sys, path)) {
        const std::string p(path.ptr, path.len);
        sources[0] = p + "/factor";
        sources[1] = p + "/weight";
    } else {
        sources[0] = "orbital/alpha/coefficient";
        sources[1] = "orbital/alpha/occupation";
        sources[2] = "orbital/beta/coefficient";
        sources[3] = "orbital/beta/occupation";
    }
    for (const std::string& p : sources) {
        if (p.empty()) continue;
        const md_attribute_t* a = md_attributes_find(&sys.attributes, md_str_t{p.data(), p.size()});
        const uint64_t v = a ? md_attributes_version(&sys.attributes, a->id) : 0;
        key = md_hash64(&v, sizeof(v), key);
    }
    return key ? key : 1;
}

// The density's tiles, built from its factors a block at a time where they are orbitals read from a file
// (md_qm_density_dmat_build): the whole occupied block is never held
static std::shared_ptr<const DensityTiles> tiles_build(const md_system_t& sys, md_str_t path, const md_attribute_slice_t* slice,
                                                       const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride) {
    if (!md_qm_density_factors_available(&sys, path)) return {};
    auto t = std::make_shared<DensityTiles>();
    const md_gto_dmat_desc_t desc = { .basis = basis, .atom_xyz = atom_xyz, .atom_xyz_stride = atom_xyz_stride,
                                      .cutoff = DENSITY_TILES_CUTOFF, .margin = DENSITY_TILES_MARGIN };
    if (!md_gto_dmat_init(&t->dm, &desc, md_get_heap_allocator()) ||
        md_qm_density_dmat_build(&t->dm, &sys, path, slice, 0, task_system::executor()) == 0) {
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

std::shared_ptr<const DensityTiles> density_cache_tiles(DensityCache* cache, const md_system_t& sys, md_str_t path, const md_attribute_slice_t* slice_ptr,
                                                         const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride) {
    if (!basis || !atom_xyz) return {};
    const md_attribute_slice_t slice = slice_of(slice_ptr);
    const md_attribute_slice_t* sl = slice.num_idx ? &slice : nullptr;
    if (!cache) return tiles_build(sys, path, sl, basis, atom_xyz, atom_xyz_stride);

    const uint64_t key = density_tiles_key(sys, path, sl);
    const std::string name(path.ptr, path.len);
    // held while it builds: another asking for it waits for this build rather than making its own
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t slot = cache->entries.size();
    for (size_t i = 0; i < cache->entries.size(); ++i) {
        DensityCache::Entry& e = cache->entries[i];
        if (e.path != name || !slice_eq(e.slice, slice)) continue;
        const DensityTiles* t = e.tiles.get();
        if (e.key == key && t && t->dm.num_shells == basis->num_shells && md_gto_dmat_pattern_valid(&t->dm, atom_xyz, atom_xyz_stride)) {
            e.used = ++cache->clock;
            return e.tiles;
        }
        slot = i;   // of other factors, or of atoms that have moved too far: replaced
        break;
    }
    const md_tick_t t0 = md_tick_now();
    std::shared_ptr<const DensityTiles> t = tiles_build(sys, path, sl, basis, atom_xyz, atom_xyz_stride);
    if (!t) return {};
    MD_LOG_DEBUG("Density '" MD_STR_FMT "'%s: %llu of %llu tile pairs (%.1f MB in double) built from its factors in %.0f ms", MD_STR_ARG(path),
                 slice.num_idx ? " (a slice)" : "", (unsigned long long)t->dm.num_tiles, (unsigned long long)t->dm.nt * (t->dm.nt + 1) / 2,
                 2048.0 * t->dm.num_tiles / (1024.0 * 1024.0), md_tick_to_milliseconds(md_tick_now() - t0));
    if (slot == cache->entries.size()) {
        // a new sliced one: the least recently used of the others goes once there are as many as are kept
        if (slice.num_idx) {
            size_t n = 0, oldest = cache->entries.size();
            for (size_t i = 0; i < cache->entries.size(); ++i) {
                if (!cache->entries[i].slice.num_idx) continue;
                ++n;
                if (oldest == cache->entries.size() || cache->entries[i].used < cache->entries[oldest].used) oldest = i;
            }
            if (n >= DENSITY_CACHE_MAX_SLICED) slot = oldest;
        }
        if (slot == cache->entries.size()) cache->entries.push_back({});
    }
    DensityCache::Entry& e = cache->entries[slot];
    e.path = name;
    e.slice = slice;
    e.key = key;
    e.used = ++cache->clock;
    e.tiles = t;
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
            density_cache_tiles(cache, *s, md_str_t{name.data(), name.size()}, nullptr, &basis, xyz.data(), 0);
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
