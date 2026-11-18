#include "density_cache.h"

#include <md_system.h>
#include <md_qm.h>
#include <core/md_allocator.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <mutex>
#include <string>
#include <vector>

PackedDensity::~PackedDensity() {
    if (tri) md_free(md_get_heap_allocator(), tri, sizeof(double) * (dim * (dim + 1) / 2));
}

struct DensityCache {
    struct Entry {
        std::string path;
        uint64_t    key = 0;
        std::shared_ptr<const PackedDensity> density;
    };
    std::mutex         mutex;
    std::vector<Entry> entries;
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

size_t density_cache_bytes(DensityCache* cache) {
    if (!cache) return 0;
    std::lock_guard<std::mutex> lock(cache->mutex);
    size_t bytes = 0;
    for (const DensityCache::Entry& e : cache->entries) {
        if (e.density) bytes += sizeof(double) * (e.density->dim * (e.density->dim + 1) / 2);
    }
    return bytes;
}
