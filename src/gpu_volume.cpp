#include "gpu_volume.h"

#if MD_ENABLE_GPU

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <float.h>
#include <stdlib.h>

namespace {

struct Entry {
    uint32_t  id = 0;                 // 0: a free slot
    GpuVolume v = {};
};

struct {
    md_gpu_stream_t   stream = nullptr;
    md_array(Entry)   entries = 0;    // the slot of id is id - GPU_VOLUME_ID_BASE - 1
    md_array(uint32_t) free_slots = 0;
    uint64_t          version = 0;    // of the last write of any

    // The host copy the range and the histogram are taken from: of one volume's version
    uint32_t copy_id = 0;
    uint64_t copy_version = 0;
    float*   copy = nullptr;
    size_t   copy_count = 0;
    bool     copy_range = false;
    float    range_lo = 0.0f, range_hi = 0.0f;
    float*   counts = nullptr;
    int      bins = 0;
    float    lo = 0.0f, hi = 0.0f;
    uint64_t counts_version = 0;
    uint32_t counts_id = 0;
} g;

md_allocator_i* heap() { return md_get_heap_allocator(); }

Entry* find(uint32_t id) {
    if (!gpu_volume_is(id)) return nullptr;
    const size_t slot = (size_t)(id - GPU_VOLUME_ID_BASE - 1);
    if (slot >= md_array_size(g.entries) || g.entries[slot].id != id) return nullptr;
    return &g.entries[slot];
}

void release(GpuVolume* v) {
    if (v->tex) md_gpu_texture_destroy(v->tex);
    if (v->normals) md_gpu_texture_destroy(v->normals);
    v->tex = nullptr;
    v->normals = nullptr;
}

void copy_drop() {
    free(g.copy);
    g.copy = nullptr;
    g.copy_count = 0;
    g.copy_id = 0;
    g.copy_version = 0;
    g.copy_range = false;
}

// The host copy of the volume as it is (read back when it is of another volume or version)
const float* copy_of(uint32_t id, size_t* count) {
    const Entry* e = find(id);
    if (!e || !e->v.tex) return nullptr;
    if (g.copy && g.copy_id == id && g.copy_version == e->v.version) {
        *count = g.copy_count;
        return g.copy;
    }
    copy_drop();
    const size_t n = (size_t)e->v.dim[0] * e->v.dim[1] * e->v.dim[2];
    float* dst = (float*)malloc(sizeof(float) * MD_MAX(n, (size_t)1));
    if (!dst || !gpu_volume_read(id, dst, n)) {
        free(dst);
        return nullptr;
    }
    g.copy = dst;
    g.copy_count = n;
    g.copy_id = id;
    g.copy_version = e->v.version;
    *count = n;
    return dst;
}

}  // namespace

void gpu_volumes_init(md_gpu_stream_t stream) {
    g.stream = stream;
}

void gpu_volumes_shutdown() {
    for (size_t i = 0; i < md_array_size(g.entries); ++i) {
        if (g.entries[i].id) release(&g.entries[i].v);
    }
    md_array_free(g.entries, heap());
    md_array_free(g.free_slots, heap());
    g.entries = 0;
    g.free_slots = 0;
    copy_drop();
    free(g.counts);
    g.counts = nullptr;
    g.stream = nullptr;
}

uint32_t gpu_volume_create() {
    size_t slot = 0;
    if (md_array_size(g.free_slots) > 0) {
        slot = g.free_slots[md_array_size(g.free_slots) - 1];
        md_array_pop(g.free_slots);
    } else {
        slot = md_array_size(g.entries);
        md_array_push(g.entries, Entry{}, heap());
    }
    Entry& e = g.entries[slot];
    e = Entry{};
    e.id = GPU_VOLUME_ID_BASE + 1 + (uint32_t)slot;
    return e.id;
}

void gpu_volume_destroy(uint32_t* id) {
    if (!id) return;
    if (Entry* e = find(*id)) {
        release(&e->v);
        if (g.copy_id == e->id) copy_drop();
        const size_t slot = (size_t)(e->id - GPU_VOLUME_ID_BASE - 1);
        *e = Entry{};
        md_array_push(g.free_slots, (uint32_t)slot, heap());
        g.version += 1;
    }
    *id = 0;
}

uint32_t gpu_volume_ensure(uint32_t* id) {
    if (!find(*id)) *id = gpu_volume_create();
    return *id;
}

GpuVolume gpu_volume_get(uint32_t id) {
    const Entry* e = find(id);
    return e ? e->v : GpuVolume{};
}

void gpu_volume_set(uint32_t id, md_gpu_texture_t tex, md_gpu_texture_t normals, const md_float4x4_t& texture_to_world,
                    md_gpu_sync_t ready) {
    Entry* e = find(id);
    if (!e) {
        // of a volume no more: what was made for it goes
        if (tex) md_gpu_texture_destroy(tex);
        if (normals) md_gpu_texture_destroy(normals);
        return;
    }
    release(&e->v);
    e->v.tex = tex;
    e->v.normals = normals;
    for (int i = 0; i < 3; ++i) e->v.dim[i] = 0;
    if (tex) {
        const md_gpu_texture_desc_t* d = md_gpu_texture_desc(tex);
        e->v.dim[0] = d->width;
        e->v.dim[1] = d->height;
        e->v.dim[2] = d->depth_or_layers;
    }
    e->v.texture_to_world = texture_to_world;
    e->v.ready = ready;
    e->v.version = ++g.version;
}

bool gpu_volume_upload(uint32_t id, md_gpu_format_t format, const int dim[3], const void* data, const md_float4x4_t& texture_to_world) {
    if (!find(id) || !g.stream || !data || dim[0] <= 0 || dim[1] <= 0 || dim[2] <= 0) return false;
    md_gpu_texture_desc_t desc = {};
    desc.type = MD_GPU_TEX_3D;
    desc.format = format;
    desc.usage = MD_GPU_TEX_SAMPLED;
    desc.width = (uint32_t)dim[0];
    desc.height = (uint32_t)dim[1];
    desc.depth_or_layers = (uint32_t)dim[2];
    desc.label = "Representation volume";
    md_gpu_texture_t tex = md_gpu_texture_create(g.stream, &desc);
    if (!tex) {
        MD_LOG_ERROR("A volume of %dx%dx%d could not be made: %s", dim[0], dim[1], dim[2], md_gpu_last_error());
        return false;
    }
    const size_t size = md_gpu_texture_region_size(tex, nullptr);
    if (!md_gpu_upload_texture(g.stream, tex, nullptr, data, size)) {
        md_gpu_texture_destroy(tex);
        return false;
    }
    gpu_volume_set(id, tex, nullptr, texture_to_world, md_gpu_stream_record(g.stream));
    return true;
}

void gpu_volume_clear(uint32_t id) {
    if (Entry* e = find(id)) {
        if (e->v.tex || e->v.normals) gpu_volume_set(id, nullptr, nullptr, e->v.texture_to_world, md_gpu_sync_none());
    }
}

uint64_t gpu_volume_version(uint32_t id) {
    const Entry* e = find(id);
    return e ? e->v.version : 0;
}

uint64_t gpu_volume_version_any() {
    return g.version;
}

bool gpu_volume_read(uint32_t id, float* dst, size_t count) {
    const Entry* e = find(id);
    if (!e || !e->v.tex || !dst || !g.stream) return false;
    const md_gpu_texture_desc_t* d = md_gpu_texture_desc(e->v.tex);
    if (d->format != MD_GPU_FORMAT_R32_FLOAT) return false;
    const size_t size = md_gpu_texture_region_size(e->v.tex, nullptr);
    if (size != count * sizeof(float)) return false;
    md_gpu_mem_t rb = md_gpu_malloc(g.stream, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    md_gpu_stream_wait(g.stream, e->v.ready);
    const bool ok = md_gpu_copy_from_texture(g.stream, rb.gpu, e->v.tex, nullptr);
    if (ok) {
        md_gpu_stream_sync(g.stream);
        md_memcpy(dst, rb.cpu, size);
    }
    md_gpu_free(g.stream, rb.gpu);
    return ok;
}

bool gpu_volume_value_range(uint32_t id, float* lo, float* hi) {
    size_t n = 0;
    const float* v = copy_of(id, &n);
    if (!v || n == 0) return false;
    if (!g.copy_range) {
        float a = FLT_MAX, b = -FLT_MAX;
        for (size_t i = 0; i < n; ++i) {
            a = MD_MIN(a, v[i]);
            b = MD_MAX(b, v[i]);
        }
        g.range_lo = a;
        g.range_hi = b;
        g.copy_range = true;
    }
    *lo = g.range_lo;
    *hi = g.range_hi;
    return g.range_lo <= g.range_hi;
}

const float* gpu_volume_histogram(uint32_t id, float lo, float hi, int bins) {
    if (bins <= 0 || bins > 4096 || !(hi > lo)) return nullptr;
    const uint64_t version = gpu_volume_version(id);
    if (g.counts && g.counts_id == id && g.counts_version == version && g.lo == lo && g.hi == hi && g.bins == bins) return g.counts;
    size_t n = 0;
    const float* v = copy_of(id, &n);
    if (!v) return nullptr;
    if (g.bins != bins || !g.counts) {
        free(g.counts);
        g.counts = (float*)malloc(sizeof(float) * bins);
        if (!g.counts) return nullptr;
    }
    for (int i = 0; i < bins; ++i) g.counts[i] = 0.0f;
    const float scale = (float)bins / (hi - lo);
    for (size_t i = 0; i < n; ++i) {
        const float t = (v[i] - lo) * scale;
        if (t >= 0.0f && t < (float)bins) g.counts[(int)t] += 1.0f;
    }
    g.counts_id = id;
    g.counts_version = version;
    g.lo = lo;
    g.hi = hi;
    g.bins = bins;
    return g.counts;
}

#endif
