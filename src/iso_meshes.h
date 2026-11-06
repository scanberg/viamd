#pragma once

// Isosurface meshes of the electronic structure, extracted from the fields themselves (md_gto_iso) rather than
// from a volume, for the mesh rendering mode (Settings::isosurface_rendering, isomesh_utils.h for the drawing).
//
// A set of meshes is up to ISO_MESHES_MAX_SLOTS slots, each a field at an isovalue (one or both sheets), all in
// one GL mesh. On the GPU (md_gto_iso_gpu) the slots are extracted together - the corner values are shared,
// a field in several slots costs its functions once - and the result is read back once; without a GPU device
// the double precision CPU path (md_gto_iso cache + mesh) does the same, slower. On the GPU the extraction is a
// job (md_gto_iso_gpu_job_t) advanced by every update for a few milliseconds at most: the mesh shown stays the
// previous one until the job is done, so a long extraction never stalls the frame (iso_meshes_progress tells how
// far it is). The CPU path blocks.
//
// Before a job, the fields go to the device (md_gto_iso_gpu_set_t), which for a large density is seconds of CPU
// (its density matrix): a pool task does it, on a GPU stream of its own, while the mesh shown stays. The fields
// of an electronic structure representation are made by pool tasks too. A field that a pool task may be reading
// is let go of with iso_field_release, never md_gto_iso_field_destroy.
//
// Positions are in bohr, at the atoms of the current frame: the meshes go to world space by a scale of
// BOHR_TO_ANGSTROM (iso_meshes_model_matrix).

#include <gfx/isomesh_utils.h>
#include <color_scale.h>
#include <core/md_vec_math.h>
#include <core/md_str.h>
#include <core/md_gpu.h>

#include <stdint.h>
#include <stddef.h>

struct ApplicationState;
struct ElectronicStructureRepresentation;
struct Representation;
struct IsoDesc;
struct md_gto_iso_field_t;
struct md_gto_iso_cache_t;
struct md_attribute_slice_t;
struct md_system_t;
enum class VolumeResolution;

static constexpr uint32_t ISO_MESHES_MAX_SLOTS = 16;

struct IsoMeshSlot {
    const md_gto_iso_field_t* field = nullptr;
    double   iso    = 0.0;      // > 0
    uint32_t sheets = 0;        // MD_GTO_ISO_SHEET_*
};

// Where a field's atoms are: the system whose attributes it is made of, and the positions (bohr) of its basis atoms
struct IsoFieldGeometry {
    const md_system_t* sys = nullptr;
    const vec3_t* atom_xyz = nullptr;
    size_t num_atoms = 0;
};

struct IsoMeshes {
    isomesh::Mesh mesh = {};
    uint32_t num_slots = 0;
    uint32_t index_offset[ISO_MESHES_MAX_SLOTS][2] = {};   // [slot][0] the positive sheet, [1] the negative
    uint32_t index_count [ISO_MESHES_MAX_SLOTS][2] = {};
    uint64_t key = 0;           // of what the mesh is of
    bool     ok  = false;       // whether the last extraction succeeded (a failed one is not retried until the key changes)
    bool     on_gpu = false;    // which path made it
    double   seconds = 0.0;     // what the last extraction took
    // The vertices, bohr, where a field is evaluated (surface_field_evaluate_points): on the host, packed
    // float3, and after a GPU extraction also where the extraction left them on the device (float4 per
    // vertex), so the evaluation reads them there instead of uploading them again
    float*   xyz = nullptr;
    uint32_t xyz_count = 0;
    uint32_t* tri_index = nullptr;  // the triangles on the host too, three per triangle, for the field's gradient
    uint32_t tri_index_count = 0;   // along the surface (isomesh::mesh_field_gradient)
    md_gpu_addr_t gpu_xyz = 0;

    // The GPU extraction under way (md_gto_iso_gpu_job_t), of the slots whose key is job_key: until it is done
    // the mesh above is the previous one
    void*    job = nullptr;
    uint64_t job_key = 0;
    uint32_t job_slots = 0;
    int64_t  job_start = 0;     // md_tick_t: when it began

    // Kept between extractions of the same fields
    uint64_t fields_key = 0;
    void*    gpu_set = nullptr;     // md_gto_iso_gpu_set_t
    void*    set_build = nullptr;   // the set being made on a worker thread (of these fields, or of the ones before)
    uint64_t set_failed_key = 0;    // the fields no set could be made of: extracted on the CPU
    struct Cache {
        const md_gto_iso_field_t* field = nullptr;
        md_gto_iso_cache_t* cache = nullptr;
        double lo = 0.0, hi = 0.0;
    } cache[ISO_MESHES_MAX_SLOTS] = {};
    double cache_cell_size = 0.0;
};

// Brings the mesh up to date with the slots: extracts unless it already holds them. fields_key names the
// contents of the fields (a field can be freed and another made at the same address). geometry: where the
// fields' atoms are, the current frame's when NULL. Returns m->ok.
bool iso_meshes_update(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t num_slots, double cell_size, uint64_t fields_key,
                       const IsoFieldGeometry* geometry = nullptr);
void iso_meshes_free(IsoMeshes* m, ApplicationState* state);

// A field no longer used: destroyed now, or when the pool task that reads it ends
void iso_field_release(md_gto_iso_field_t* field);

// Frees what the pool tasks that are no longer wanted left behind, once they end: once a frame, whether or not any
// meshes are drawn (a representation's field let go of in the volume rendering mode is held until then).
void iso_meshes_collect();

// Waits for the pool tasks under way and frees what they made. Before the GPU device and the task system go.
void iso_meshes_shutdown(ApplicationState* state);

// The GPU extractions under way for longer than a moment, for a progress display: how far each is (0 to 1), or
// negative for what is being prepared on a worker thread (a field, a set), which has no fraction. Returns how many
// (at most cap).
size_t iso_meshes_progress(float* fraction, size_t cap);

mat4_t iso_meshes_model_matrix();

// Lattice spacing in bohr for a volume resolution: the same three steps of detail as the volumes
double iso_meshes_cell_size(VolumeResolution resolution);

// Fields at the atoms of the current frame. NULL on failure. Free with iso_field_release once given to meshes.
md_gto_iso_field_t* iso_field_create_orbital(const ApplicationState* state, str_t coefficient_path, const md_attribute_slice_t* slice);
md_gto_iso_field_t* iso_field_create_density(const ApplicationState* state, str_t density_path, const md_attribute_slice_t* slice);
md_gto_iso_field_t* iso_field_create_density_matrix(const ApplicationState* state, const double* density_matrix, size_t dim);

// The surfaces of an electronic structure representation as meshes: the field its volume is evaluated from,
// cut at the values of 'iso' (electronic_structure_iso_desc_init). Fills 'out' (mesh, range, model matrix,
// colour, optical density; the caller sets the colouring) and returns how many there are; 0 when there is
// nothing to draw or the extraction failed.
//
// Coloured by a field (SurfaceColoring::Field), the field is evaluated at the vertices of the meshes
// (surface_field_evaluate_points: on the GPU where it can be) and the surfaces come with
// Coloring_VertexField; the caller sets the colour map and the range. The colour scale's range follows
// the values on the vertices when it is automatic.
size_t electronic_structure_mesh_surfaces(isomesh::Surface* out, size_t cap, ApplicationState* state, Representation* rep, const IsoDesc& iso);

// The spread of the field over the vertices of the representation's meshes; not valid without any
ColorScaleSpan electronic_structure_mesh_field_span(const ElectronicStructureRepresentation& es);
void   electronic_structure_mesh_free(ApplicationState* state, ElectronicStructureRepresentation* es);
