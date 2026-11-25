#pragma once

// The volumes of the representations, which the main view draws (md_render_bridge.cpp: md_compose's fields, colour
// volumes and normals, the planes' images): textures on md_gpu, each named by an id. What is written into a volume
// replaces what it held, whole: a new texture, the old one released once the GPU is done with it, so what is drawn is
// never half written and a frame never waits for an evaluation. Its version counts the writes (the view's key).
//
// The ids are of a range of their own (gpu_volume_is), apart from OpenGL's texture names: the evaluations (viamd.cpp)
// take either as the name of what they fill, a component window's GL texture or a representation's volume.

#include <stddef.h>
#include <stdint.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#include <core/md_vec_math.h>

constexpr uint32_t GPU_VOLUME_ID_BASE = 0x40000000u;
static inline bool gpu_volume_is(uint32_t id) { return id >= GPU_VOLUME_ID_BASE; }

// The stream its uploads and readbacks go through: the evaluations' (state->gpu_stream)
void gpu_volumes_init(md_gpu_stream_t stream);
void gpu_volumes_shutdown();

// An id for a volume of nothing yet. Destroyed with what it holds, and set to 0.
uint32_t gpu_volume_create();
void     gpu_volume_destroy(uint32_t* id);
// The id in *id, made if there is none
uint32_t gpu_volume_ensure(uint32_t* id);

// What it holds: its texture (null: nothing yet), the normals of its field (null: none), their dims, where they are
// (their texture space [0, 1]^3 to the world, Angstrom), the version of the write that made them, and the point of the
// evaluations' stream that wrote them (to wait for on another stream)
struct GpuVolume {
    md_gpu_texture_t tex = nullptr;
    md_gpu_texture_t normals = nullptr;
    uint32_t         dim[3] = {};
    md_float4x4_t    texture_to_world = {};
    uint64_t         version = 0;
    md_gpu_sync_t    ready = {};
};
GpuVolume gpu_volume_get(uint32_t id);

// What it holds replaced: by textures written by the work in the stream up to `ready` (taken over), the old released
void gpu_volume_set(uint32_t id, md_gpu_texture_t tex, md_gpu_texture_t normals, const md_float4x4_t& texture_to_world,
                    md_gpu_sync_t ready);
// ... by a texture of the format and dims (3D; one layer: an image) filled from the host, texels tightly packed
bool gpu_volume_upload(uint32_t id, md_gpu_format_t format, const int dim[3], const void* data, const md_float4x4_t& texture_to_world);
// ... by nothing
void gpu_volume_clear(uint32_t id);

uint64_t gpu_volume_version(uint32_t id);
// Counted up with every write of any volume
uint64_t gpu_volume_version_any();

// Its values (of an R32_FLOAT volume), read back and waited for: count floats, its voxels
bool gpu_volume_read(uint32_t id, float* dst, size_t count);
// Its values' range, and their histogram over [lo, hi] in bins (the texels in each; those outside not counted), from a
// copy read back once per version, of one volume at a time (the transfer function editor's). The counts until the
// next call, or null.
bool         gpu_volume_value_range(uint32_t id, float* lo, float* hi);
const float* gpu_volume_histogram(uint32_t id, float lo, float hi, int bins);
#endif
