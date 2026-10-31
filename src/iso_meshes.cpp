#include <iso_meshes.h>

#include <viamd.h>
#include <task_system.h>
#include <surface_field.h>

#include <md_gto_iso.h>
#include <md_qm.h>
#include <core/md_allocator.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <math.h>
#include <stdlib.h>

// Newton steps onto the surface and tangential relaxation passes: display quality (md_gto_iso_mesh_desc_t)
static constexpr uint32_t PROJECT_ITERATIONS = 8;
static constexpr uint32_t RELAX_ITERATIONS   = 3;

mat4_t iso_meshes_model_matrix() {
    const float s = (float)BOHR_TO_ANGSTROM;
    return mat4_scale(s, s, s);
}

double iso_meshes_cell_size(VolumeResolution resolution) {
    switch (resolution) {
    case VolumeResolution::Low:  return 0.30;
    case VolumeResolution::High: return 0.12;
    case VolumeResolution::Mid:
    default:                     return 0.20;
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Fields
// ---------------------------------------------------------------------------------------------------------------

struct FieldContext {
    md_gto_basis_t basis = {};
    const vec3_t*  atom_xyz = nullptr;   // bohr
    size_t         num_atoms = 0;
};

// The current frame's: what the public functions below make fields of
static IsoFieldGeometry current_geometry(md_temp_scope_t temp, const ApplicationState* state) {
    IsoFieldGeometry g;
    g.sys = &state->mold.sys;
    g.atom_xyz = basis_atom_positions_extract(&g.num_atoms, temp, state->mold.sys, state->mold.state);
    return g;
}

static bool field_context(FieldContext* ctx, md_temp_scope_t temp, const IsoFieldGeometry& g) {
    if (!g.sys || !md_gto_basis_extract_attributes(&ctx->basis, &g.sys->attributes, md_temp_allocator(temp))) {
        return false;
    }
    ctx->atom_xyz = g.atom_xyz;
    ctx->num_atoms = g.num_atoms;
    if (!ctx->atom_xyz || ctx->num_atoms != md_gto_basis_num_atoms(&ctx->basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis and the atom positions do not agree");
        return false;
    }
    return true;
}

static md_gto_iso_field_t* field_create_orbital(const IsoFieldGeometry& g, str_t coefficient_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    FieldContext ctx;
    if (!field_context(&ctx, temp, g)) return nullptr;
    size_t num_ao = 0;
    const double* coeffs = orbital_coefficients_extract(&num_ao, temp, *g.sys, coefficient_path, slice);
    if (!coeffs) return nullptr;
    if (num_ao != md_gto_basis_num_ao(&ctx.basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis spans %zu atomic orbitals and '" STR_FMT "' %zu", md_gto_basis_num_ao(&ctx.basis), STR_ARG(coefficient_path), num_ao);
        return nullptr;
    }
    return md_gto_iso_field_create(&ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), coeffs, md_get_heap_allocator());
}

md_gto_iso_field_t* iso_field_create_orbital(const ApplicationState* state, str_t coefficient_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    return field_create_orbital(current_geometry(temp, state), coefficient_path, slice);
}

static md_gto_iso_field_t* field_create_density_matrix(const IsoFieldGeometry& g, const double* density_matrix, size_t dim) {
    if (!density_matrix || dim == 0) return nullptr;
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    FieldContext ctx;
    if (!field_context(&ctx, temp, g)) return nullptr;
    if (dim != md_gto_basis_num_ao(&ctx.basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis spans %zu atomic orbitals and the density matrix %zu", md_gto_basis_num_ao(&ctx.basis), dim);
        return nullptr;
    }
    return md_gto_iso_field_create_density(&ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), density_matrix, md_get_heap_allocator());
}

md_gto_iso_field_t* iso_field_create_density_matrix(const ApplicationState* state, const double* density_matrix, size_t dim) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    return field_create_density_matrix(current_geometry(temp, state), density_matrix, dim);
}

static md_gto_iso_field_t* field_create_density(const IsoFieldGeometry& g, str_t density_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    size_t dim = 0;
    const double* D = density_matrix_extract(&dim, temp, *g.sys, density_path, slice);
    return D ? field_create_density_matrix(g, D, dim) : nullptr;
}

md_gto_iso_field_t* iso_field_create_density(const ApplicationState* state, str_t density_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    return field_create_density(current_geometry(temp, state), density_path, slice);
}

// An SCF density (orbital/{alpha,beta,total,difference}/density) from the orbitals it is the sum of
// (md_qm_density_orbitals): the field takes them as they are, where from the matrix it would factorize
// it again, which for a protein's density is minutes rather than a second. NULL when the system does
// not publish the orbitals, or the density is zero.
static md_gto_iso_field_t* field_create_scf_density(const IsoFieldGeometry& g, str_t density_path) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    FieldContext ctx;
    if (!field_context(&ctx, temp, g)) return nullptr;
    double* rows = nullptr;
    double* weights = nullptr;
    size_t num_ao = 0;
    const size_t count = md_qm_density_orbitals(&rows, &weights, &num_ao, g.sys, density_path, md_temp_allocator(temp));
    if (count == 0) return nullptr;
    if (num_ao != md_gto_basis_num_ao(&ctx.basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis spans %zu atomic orbitals and the orbitals of '" STR_FMT "' %zu", md_gto_basis_num_ao(&ctx.basis), STR_ARG(density_path), num_ao);
        return nullptr;
    }
    return md_gto_iso_field_create_products(&ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), rows, nullptr, weights, count, md_get_heap_allocator());
}

// ---------------------------------------------------------------------------------------------------------------
// Meshes
// ---------------------------------------------------------------------------------------------------------------

#if MD_ENABLE_GPU
// The jobs under way, for iso_meshes_progress: a job and when it began
static constexpr uint32_t MAX_JOBS = 64;
static md_gto_iso_gpu_job_t* active_job[MAX_JOBS];
static md_tick_t             active_t0[MAX_JOBS];
static uint32_t              num_active = 0;

static void job_register(md_gto_iso_gpu_job_t* job) {
    if (num_active < MAX_JOBS) {
        active_job[num_active] = job;
        active_t0[num_active] = md_tick_now();
        num_active++;
    }
}

static void job_unregister(md_gto_iso_gpu_job_t* job) {
    for (uint32_t i = 0; i < num_active; ++i) {
        if (active_job[i] == job) {
            num_active--;
            active_job[i] = active_job[num_active];
            active_t0[i] = active_t0[num_active];
            return;
        }
    }
}

// Drops the job under way, if any: what it made is lost, the mesh shown stays
static void job_cancel(IsoMeshes* m) {
    if (!m->job) return;
    md_gto_iso_gpu_job_t* job = (md_gto_iso_gpu_job_t*)m->job;
    job_unregister(job);
    md_gto_iso_gpu_job_end(job, nullptr);
    m->job = nullptr;
    m->job_key = 0;
    m->job_slots = 0;
}
#endif

// ---------------------------------------------------------------------------------------------------------------
// Made on worker threads
// ---------------------------------------------------------------------------------------------------------------
// A large density's field and its set on the device take seconds of CPU: the orbitals copied, the field made
// of them, and for a set of more functions than the factored kernels take the density matrix D = C W C^T and
// its upload (md_gto_iso_gpu_set_create). They are made by pool tasks, the set on a stream of its own that the
// task destroys once its uploads are done, while the meshes shown stay what they were.
//
// What a task reads must outlive it: the system's attributes do (a load waits for the pool), the positions are
// its own copy, and the fields a set is made of are held for it - a field let go of while a set of it is being
// made (iso_field_release) goes when that task ends. A task whose result is no longer wanted is orphaned and
// cleaned up once it ends (iso_meshes_collect); nothing here waits for one, but iso_meshes_shutdown.

// What an electronic structure representation's field is of, copied for the task that makes it
struct IsoFieldSpec {
    ElectronicStructureSource source = ElectronicStructureSource::MolecularOrbital;
    ElectronicStructureSpin spin = ElectronicStructureSpin::Alpha;
    ElectronicStructureNtoComponent nto_component = ElectronicStructureNtoComponent::Particle;
    ElectronicStructureTransitionDensityComponent transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;
    int orbital_idx = 0;
    int excited_state_idx = 0;
    int nto_lambda_idx = 0;
    md_attribute_id_t density_property_key = MD_ATTRIBUTE_INVALID;
};

struct IsoBuild {
    enum Kind { Field, Set } kind = Field;
    task_system::ID task = task_system::INVALID_ID;
    bool      ran = false;       // an interrupted pool task never runs its function
    bool      orphan = false;    // nobody wants what it makes
    md_tick_t t0 = 0;
    uint64_t  key = 0;
    const md_system_t* sys = nullptr;
    vec3_t*   xyz = nullptr;     // bohr, its copy
    size_t    num_atoms = 0;
    // Field
    IsoFieldSpec spec = {};
    md_gto_iso_field_t* field = nullptr;
    // Set
    const md_gto_iso_field_t* fields[ISO_MESHES_MAX_SLOTS] = {};
    uint32_t  num_fields = 0;
#if MD_ENABLE_GPU
    // The stream the set is made on, the task's; destroyed by the task once its work is done (giving back its
    // staging, the density matrix's upload: tens of ms on the GL thread). The GL thread's stream frees a set
    // nobody took.
    md_gpu_stream_t stream = nullptr;
    md_gpu_stream_t gl_stream = nullptr;
    void* set = nullptr;              // md_gto_iso_gpu_set_t
#endif
};

static IsoBuild* g_builds[256];
static uint32_t  g_num_builds = 0;
static md_gto_iso_field_t* g_dead_fields[256];   // let go of while a set of them was being made
static uint32_t  g_num_dead = 0;

static md_gto_iso_field_t* field_make(const IsoFieldGeometry& g, const IsoFieldSpec& spec);

static IsoBuild* build_new(IsoBuild::Kind kind, uint64_t key, const md_system_t* sys, const vec3_t* xyz, size_t num_atoms) {
    if (g_num_builds == ARRAY_SIZE(g_builds)) return nullptr;
    IsoBuild* b = new IsoBuild();
    b->kind = kind;
    b->key = key;
    b->sys = sys;
    b->t0 = md_tick_now();
    if (xyz && num_atoms) {
        b->xyz = (vec3_t*)malloc(sizeof(vec3_t) * num_atoms);
        if (b->xyz) {
            MEMCPY(b->xyz, xyz, sizeof(vec3_t) * num_atoms);
            b->num_atoms = num_atoms;
        }
    }
    g_builds[g_num_builds++] = b;
    return b;
}

static bool build_running(const IsoBuild* b) {
    return b && task_system::task_is_running(b->task);
}

// On a worker thread
static void build_run(IsoBuild* b) {
    b->ran = true;
    IsoFieldGeometry g;
    g.sys = b->sys;
    g.atom_xyz = b->xyz;
    g.num_atoms = b->num_atoms;
    if (b->kind == IsoBuild::Field) {
        b->field = field_make(g, b->spec);
        return;
    }
#if MD_ENABLE_GPU
    md_temp_scope_t temp = md_temp_begin();
    FieldContext ctx;
    if (field_context(&ctx, temp, g)) {
        b->set = md_gto_iso_gpu_set_create(b->stream, &ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), b->fields, b->num_fields, md_get_heap_allocator());
    }
    md_temp_end(temp);
    md_gpu_stream_destroy(b->stream);   // once its uploads are done: the set is then ready for any stream
    b->stream = nullptr;
#endif
}

static void build_start(IsoBuild* b) {
    // and the worker's temp arenas given back after: a density's orbitals are read there, hundreds of MB that a
    // pool thread would otherwise hold until it exits
    b->task = task_system::create_pool_task(STR_LIT("##Isosurface mesh preparation"), [b]() { build_run(b); md_temp_trim(); });
    task_system::enqueue_task(b->task);
}

static bool field_held(const md_gto_iso_field_t* f) {
    for (uint32_t i = 0; i < g_num_builds; ++i) {
        const IsoBuild* b = g_builds[i];
        if (b->kind != IsoBuild::Set || !build_running(b)) continue;
        for (uint32_t k = 0; k < b->num_fields; ++k) if (b->fields[k] == f) return true;
    }
    return false;
}

// Fields destroyed by pool tasks: a protein's density field is a few hundred MB in large blocks, tens of ms to
// give back, at the frame that takes the next one
static task_system::ID g_destroys[64];
static uint32_t        g_num_destroys = 0;

static void field_destroy(md_gto_iso_field_t* field) {
    for (uint32_t i = 0; i < g_num_destroys;) {
        if (!task_system::task_is_running(g_destroys[i])) {
            g_destroys[i] = g_destroys[--g_num_destroys];
            continue;
        }
        ++i;
    }
    if (g_num_destroys == ARRAY_SIZE(g_destroys)) {
        md_gto_iso_field_destroy(field);
        return;
    }
    const task_system::ID task = task_system::create_pool_task(STR_LIT("##Isosurface field release"), [field]() { md_gto_iso_field_destroy(field); });
    task_system::enqueue_task(task);
    g_destroys[g_num_destroys++] = task;
}

void iso_field_release(md_gto_iso_field_t* field) {
    if (!field) return;
    if (field_held(field) && g_num_dead < ARRAY_SIZE(g_dead_fields)) {
        g_dead_fields[g_num_dead++] = field;
        return;
    }
    if (field_held(field)) {
        // no room to defer it: the task is waited for, which is rare and only costs the time it has left
        for (uint32_t i = 0; i < g_num_builds; ++i) task_system::task_wait_for(g_builds[i]->task);
    }
    field_destroy(field);
}

// What a finished build made, dropped
static void build_destroy(IsoBuild* b) {
    if (b->field) field_destroy(b->field);
#if MD_ENABLE_GPU
    if (b->set) md_gto_iso_gpu_set_destroy(b->gl_stream, (md_gto_iso_gpu_set_t*)b->set);
    if (b->stream) md_gpu_stream_destroy(b->stream);   // of a task that never ran: nothing on it
#endif
    free(b->xyz);
    delete b;
}

static void build_remove(IsoBuild* b) {
    for (uint32_t i = 0; i < g_num_builds; ++i) {
        if (g_builds[i] == b) {
            g_builds[i] = g_builds[--g_num_builds];
            return;
        }
    }
}

// The orphans that have ended, the fields let go of that no task holds any more
static void builds_collect() {
    for (uint32_t i = 0; i < g_num_builds;) {
        IsoBuild* b = g_builds[i];
        if (b->orphan && !build_running(b)) {
            g_builds[i] = g_builds[--g_num_builds];
            build_destroy(b);
            continue;
        }
        ++i;
    }
    for (uint32_t i = 0; i < g_num_dead;) {
        if (!field_held(g_dead_fields[i])) {
            field_destroy(g_dead_fields[i]);
            g_dead_fields[i] = g_dead_fields[--g_num_dead];
            continue;
        }
        ++i;
    }
}

static void build_orphan(IsoBuild* b) {
    if (b) b->orphan = true;
}

void iso_meshes_collect() {
    builds_collect();
}

void iso_meshes_shutdown(ApplicationState* state) {
    for (uint32_t i = 0; i < g_num_builds; ++i) task_system::task_wait_for(g_builds[i]->task);
    for (uint32_t i = 0; i < g_num_builds; ++i) g_builds[i]->orphan = true;
    builds_collect();
    for (uint32_t i = 0; i < g_num_destroys; ++i) task_system::task_wait_for(g_destroys[i]);
    g_num_destroys = 0;
    (void)state;
}

size_t iso_meshes_progress(float* fraction, size_t cap) {
    size_t n = 0;
    const md_tick_t now = md_tick_now();
    // what is being made on worker threads: no fraction (negative)
    for (uint32_t i = 0; i < g_num_builds && n < cap; ++i) {
        const IsoBuild* b = g_builds[i];
        if (b->orphan || !build_running(b) || md_tick_to_seconds(now - b->t0) < 0.25) continue;
        fraction[n++] = -1.0f;
    }
#if MD_ENABLE_GPU
    for (uint32_t i = 0; i < num_active && n < cap; ++i) {
        // a job done within a few frames shows nothing
        if (md_tick_to_seconds(now - active_t0[i]) < 0.25) continue;
        fraction[n++] = md_gto_iso_gpu_job_progress(active_job[i]);
    }
#endif
    return n;
}

static void release_extraction_state(IsoMeshes* m, ApplicationState* state) {
#if MD_ENABLE_GPU
    job_cancel(m);      // it holds the set
    if (m->gpu_set) {
        md_gto_iso_gpu_set_destroy(state->gpu_stream, (md_gto_iso_gpu_set_t*)m->gpu_set);
    }
#else
    (void)state;
#endif
    // a set under way of the fields that were stays m->set_build: set_ensure lets it end before it begins another
    m->set_failed_key = 0;
    m->gpu_set = nullptr;
    for (auto& c : m->cache) {
        md_gto_iso_cache_destroy(c.cache);
        c = {};
    }
    m->cache_cell_size = 0.0;
    m->fields_key = 0;
}

// The vertices on the host, for what is evaluated at them
static void keep_xyz(IsoMeshes* m, const float* pos, size_t stride, uint32_t count) {
    free(m->xyz);
    m->xyz = count ? (float*)malloc(sizeof(float) * 3 * count) : nullptr;
    m->xyz_count = m->xyz ? count : 0;
    for (uint32_t i = 0; i < m->xyz_count; ++i) {
        for (int k = 0; k < 3; ++k) m->xyz[3 * i + k] = pos[stride * i + k];
    }
}

// The triangles on the host, for the gradient of what is evaluated at the vertices
static void keep_index(IsoMeshes* m, const uint32_t* idx, uint32_t count) {
    free(m->tri_index);
    m->tri_index = count ? (uint32_t*)malloc(sizeof(uint32_t) * count) : nullptr;
    m->tri_index_count = m->tri_index ? count : 0;
    if (m->tri_index) MEMCPY(m->tri_index, idx, sizeof(uint32_t) * count);
}

// The device copy of the vertices, which goes with the mesh it is of
static void drop_gpu_xyz(IsoMeshes* m, ApplicationState* state) {
#if MD_ENABLE_GPU
    if (m->gpu_xyz && state->gpu_stream) md_gpu_free(state->gpu_stream, m->gpu_xyz);
#else
    (void)state;
#endif
    m->gpu_xyz = 0;
}

void iso_meshes_free(IsoMeshes* m, ApplicationState* state) {
    ASSERT(m);
    release_extraction_state(m, state);
    build_orphan((IsoBuild*)m->set_build);
    drop_gpu_xyz(m, state);
    isomesh::mesh_free(&m->mesh);
    free(m->xyz);
    free(m->tri_index);
    const uint64_t version = m->mesh.version;
    *m = {};
    m->mesh.version = version;
}

static void set_ranges(IsoMeshes* m, uint32_t num_slots, const uint32_t offset[][2], const uint32_t count[][2]) {
    m->num_slots = num_slots;
    for (uint32_t k = 0; k < ISO_MESHES_MAX_SLOTS; ++k) {
        for (int s = 0; s < 2; ++s) {
            m->index_offset[k][s] = k < num_slots ? offset[k][s] : 0;
            m->index_count[k][s]  = k < num_slots ? count[k][s]  : 0;
        }
    }
}

#if MD_ENABLE_GPU
// How long an update waits for a job's GPU work at most: a job of a few milliseconds is done within the frame
// that began it or the next, a longer one goes on over frames
static constexpr double JOB_POLL_MS = 3.0;

// Starts the extraction of the slots as a job on the GPU; false when it cannot be (no device, no set for these
// fields)
// The slots' fields on the device (m->gpu_set), made by a pool task on a stream of its own (IsoBuild): 1 when
// they are there, 0 while they are being made, -1 when they cannot be (no device, fields the GPU path does
// not take). g: where the fields' atoms are.
static int set_ensure(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t n, const IsoFieldGeometry& g) {
    if (!state->gpu_device || !state->gpu_stream) return -1;
    if (m->gpu_set) return 1;
    if (m->set_failed_key == m->fields_key) return -1;
    IsoBuild* b = (IsoBuild*)m->set_build;
    if (b) {
        // One at a time, also of fields that have since gone: for a large density each holds the density
        // matrix, hundreds of MB on the host and on the device, and a frame played every second would pile them up
        if (build_running(b)) return 0;
        m->set_build = nullptr;
        build_remove(b);
        if (b->key == m->fields_key && b->ran) {
            if (b->set) {
                // its uploads are done: the task waited for them before it ended (destroying its stream)
                m->gpu_set = b->set;
                b->set = nullptr;
                MD_LOG_DEBUG("Isosurface meshes: %u fields on the device in %.1f ms (on a worker thread)", n, md_tick_to_milliseconds(md_tick_now() - b->t0));
            } else {
                m->set_failed_key = m->fields_key;
                MD_LOG_DEBUG("Isosurface meshes: no GPU set for these fields, extracting on the CPU");
            }
        }
        build_destroy(b);
        if (m->gpu_set) return 1;
        if (m->set_failed_key == m->fields_key) return -1;
        // of other fields, or interrupted before it ran: one for these
    }
    IsoBuild* nb = build_new(IsoBuild::Set, m->fields_key, g.sys, g.atom_xyz, g.num_atoms);
    if (!nb) return -1;
    nb->num_fields = n;
    for (uint32_t k = 0; k < n; ++k) nb->fields[k] = slots[k].field;     // held: iso_field_release
    nb->stream = md_gpu_stream_create(state->gpu_device, MD_GPU_STREAM_COMPUTE, "isosurface fields");
    nb->gl_stream = state->gpu_stream;
    if (!nb->stream) {
        build_remove(nb);
        build_destroy(nb);
        return -1;
    }
    build_start(nb);
    m->set_build = nb;
    return 0;
}

static bool job_start(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t n, double cell_size, uint64_t key) {
    if (!state->gpu_device || !state->gpu_stream || !m->gpu_set) return false;
    md_gto_iso_gpu_desc_t desc = {};
    desc.cell_size = cell_size;
    desc.project_iterations = PROJECT_ITERATIONS;
    desc.relax_iterations = RELAX_ITERATIONS;
    desc.read_back = true;      // the mesh comes to the host with the job, for GL
    for (uint32_t k = 0; k < n; ++k) {
        desc.iso[k] = slots[k].iso;
        desc.sheets[k] = slots[k].sheets;
    }
    md_gto_iso_gpu_job_t* job = md_gto_iso_gpu_job_begin(state->gpu_stream, (const md_gto_iso_gpu_set_t*)m->gpu_set, &desc);
    if (!job) return false;
    m->job = job;
    m->job_key = key;
    m->job_slots = n;
    m->job_start = md_tick_now();
    job_register(job);
    return true;
}

// Advances the job. Done: its mesh replaces the one shown (m->key becomes the job's). Returns 1 while it runs,
// 2 when it is done, 0 when it failed (the mesh shown is kept)
static int job_poll(IsoMeshes* m, ApplicationState* state) {
    md_gto_iso_gpu_job_t* job = (md_gto_iso_gpu_job_t*)m->job;
    // waited for in its first frames only: a short job then ends within them, and a long one costs the frame
    // nothing more than a look
    const double age_ms = md_tick_to_milliseconds(md_tick_now() - m->job_start);
    if (md_gto_iso_gpu_job_poll(job, age_ms < 4.0 * JOB_POLL_MS ? JOB_POLL_MS : 0.0) == MD_GTO_ISO_GPU_JOB_RUNNING) return 1;
    job_unregister(job);
    md_gto_iso_gpu_mesh_t gm = {};
    const bool done = md_gto_iso_gpu_job_end(job, &gm);
    const uint64_t key = m->job_key;
    const uint32_t n = m->job_slots;
    m->job = nullptr;
    m->job_key = 0;
    m->job_slots = 0;
    if (!done) return 0;
    drop_gpu_xyz(m, state);
    const bool ok = isomesh::mesh_upload(&m->mesh, gm.host_position, gm.host_normal, 4, gm.vertex_count, gm.host_index, gm.index_count);
    if (ok) {
        set_ranges(m, n, gm.index_offset, gm.index_count_of);
        keep_xyz(m, gm.host_position, 4, gm.vertex_count);
        keep_index(m, gm.host_index, gm.index_count);
        // the positions stay on the device for what is evaluated at them
        m->gpu_xyz = gm.position;
        gm.position = 0;
    }
    m->seconds = gm.seconds_total;
    md_gto_iso_gpu_mesh_free(state->gpu_stream, &gm);
    m->key = key;
    m->ok = ok;
    m->on_gpu = true;
    if (ok) {
        MD_LOG_DEBUG("Isosurface meshes: %u slots, %u vertices, %u triangles in %.1f ms (GPU, over frames)", n, m->mesh.vertex_count,
                     m->mesh.index_count / 3, m->seconds * 1.0e3);
    } else {
        MD_LOG_ERROR("Isosurface meshes: the upload failed");
    }
    return 2;
}

#endif

// The cache of a field that serves 'iso' at this spacing, made if there is none: over [iso/2, 4 iso], so a
// value dragged within that range costs only the mesh
static md_gto_iso_cache_t* cpu_cache(IsoMeshes* m, const md_gto_iso_field_t* field, double iso, double cell_size) {
    if (m->cache_cell_size != cell_size) {
        for (auto& c : m->cache) {
            md_gto_iso_cache_destroy(c.cache);
            c = {};
        }
        m->cache_cell_size = cell_size;
    }
    IsoMeshes::Cache* slot = nullptr;
    for (auto& c : m->cache) {
        if (c.field == field && c.cache) {
            if (iso >= c.lo && iso <= c.hi) return c.cache;
            slot = &c;  // the field's, rebuilt for the new range
            break;
        }
    }
    if (!slot) {
        for (auto& c : m->cache) {
            if (!c.cache) { slot = &c; break; }
        }
    }
    if (!slot) return nullptr;
    md_gto_iso_cache_destroy(slot->cache);
    *slot = {};

    md_gto_iso_cache_desc_t desc = {};
    desc.cell_size = cell_size;
    desc.iso_min = iso * 0.5;
    desc.iso_max = iso * 4.0;
    desc.defer_corners = true;
    md_gto_iso_cache_t* cache = md_gto_iso_cache_create(field, &desc, md_get_heap_allocator());
    if (!cache) return nullptr;

    // The corner values are the expensive part, and the bricks are independent
    md_gto_iso_cache_info_t info = {};
    md_gto_iso_cache_info(&info, cache);
    if (info.num_bricks > 0) {
        task_system::ID task = task_system::create_pool_task(STR_LIT("Isosurface mesh corners"), info.num_bricks,
            [cache](uint32_t beg, uint32_t end, uint32_t) { md_gto_iso_cache_evaluate(cache, beg, end - beg); }, 16);
        task_system::enqueue_task(task);
        task_system::task_wait_for(task);
    }

    slot->field = field;
    slot->cache = cache;
    slot->lo = desc.iso_min;
    slot->hi = desc.iso_max;
    return cache;
}

static bool extract_cpu(IsoMeshes* m, const IsoMeshSlot* slots, uint32_t n, double cell_size) {
    md_gto_iso_mesh_t meshes[ISO_MESHES_MAX_SLOTS] = {};
    uint32_t offset[ISO_MESHES_MAX_SLOTS][2] = {}, count[ISO_MESHES_MAX_SLOTS][2] = {};
    size_t nv = 0, ni = 0;
    bool ok = true;
    for (uint32_t k = 0; k < n && ok; ++k) {
        md_gto_iso_cache_t* cache = cpu_cache(m, slots[k].field, slots[k].iso, cell_size);
        md_gto_iso_mesh_desc_t desc = {};
        desc.iso = slots[k].iso;
        desc.sheets = slots[k].sheets;
        desc.project_iterations = PROJECT_ITERATIONS;
        desc.relax_iterations = RELAX_ITERATIONS;
        ok = cache && md_gto_iso_mesh_extract(&meshes[k], cache, &desc, md_get_heap_allocator());
        if (ok) {
            for (int s = 0; s < 2; ++s) {
                offset[k][s] = (uint32_t)ni + meshes[k].sheet_index_offset[s];
                count[k][s]  = meshes[k].sheet_index_count[s];
            }
            nv += meshes[k].vertex_count;
            ni += meshes[k].index_count;
        }
    }
    if (ok) {
        float*    pos = (float*)calloc(3 * (nv + 1), sizeof(float));
        float*    nrm = (float*)calloc(3 * (nv + 1), sizeof(float));
        uint32_t* idx = (uint32_t*)calloc(ni + 1, sizeof(uint32_t));
        size_t v0 = 0, i0 = 0;
        for (uint32_t k = 0; k < n && pos && nrm && idx; ++k) {
            const md_gto_iso_mesh_t& mk = meshes[k];
            MEMCPY(pos + 3 * v0, mk.position, sizeof(float) * 3 * mk.vertex_count);
            MEMCPY(nrm + 3 * v0, mk.normal,   sizeof(float) * 3 * mk.vertex_count);
            for (uint32_t i = 0; i < mk.index_count; ++i) idx[i0 + i] = mk.index[i] + (uint32_t)v0;
            v0 += mk.vertex_count;
            i0 += mk.index_count;
        }
        ok = pos && nrm && idx && isomesh::mesh_upload(&m->mesh, pos, nrm, 3, (uint32_t)nv, idx, (uint32_t)ni);
        if (ok) {
            set_ranges(m, n, offset, count);
            keep_xyz(m, pos, 3, (uint32_t)nv);
            keep_index(m, idx, (uint32_t)ni);
        }
        free(pos);
        free(nrm);
        free(idx);
    }
    for (uint32_t k = 0; k < n; ++k) {
        if (meshes[k].alloc) md_gto_iso_mesh_free(&meshes[k]);
    }
    return ok;
}

static bool extract_now(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t n, double cell_size, uint64_t key);

bool iso_meshes_update(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t num_slots, double cell_size, uint64_t fields_key,
                       const IsoFieldGeometry* geometry) {
    ASSERT(m);
    ASSERT(state);
    const uint32_t n = MIN(num_slots, ISO_MESHES_MAX_SLOTS);
    builds_collect();

    uint64_t fkey = md_hash64(&n, sizeof(n), fields_key);
    uint64_t key  = md_hash64(&cell_size, sizeof(cell_size), fields_key);
    for (uint32_t k = 0; k < n; ++k) {
        fkey = md_hash64(&slots[k].field, sizeof(slots[k].field), fkey);
        key  = md_hash64(&slots[k].field, sizeof(slots[k].field), key);
        key  = md_hash64(&slots[k].iso, sizeof(slots[k].iso), key);
        key  = md_hash64(&slots[k].sheets, sizeof(slots[k].sheets), key);
    }
    key = md_hash64(&n, sizeof(n), key);
    fkey = fkey ? fkey : 1;
    key  = key ? key : 1;
#if MD_ENABLE_GPU
    // A job under way: on with it if it makes what is wanted, else dropped (the mesh shown stays until the
    // next one is done)
    if (m->job && m->job_key != key) job_cancel(m);
    if (m->job) {
        if (job_poll(m, state) != 0) return m->ok;
        // it failed: extracted now (the CPU path)
        return extract_now(m, state, slots, n, cell_size, key);
    }
#endif
    if (key == m->key) return m->ok;

    for (uint32_t k = 0; k < n; ++k) {
        if (!slots[k].field || !(slots[k].iso > 0.0)) {
            m->key = key;
            m->ok  = false;
            drop_gpu_xyz(m, state);
            return false;
        }
    }
    if (fkey != m->fields_key) {
        release_extraction_state(m, state);
        m->fields_key = fkey;
    }
#if MD_ENABLE_GPU
    // on the GPU as a job, over frames if it takes long, once the fields are on the device
    if (n > 0) {
        md_temp_scope_t temp = md_temp_begin();
        const IsoFieldGeometry g = geometry ? *geometry : current_geometry(temp, state);
        const int set = set_ensure(m, state, slots, n, g);
        md_temp_end(temp);
        if (set == 0) return m->ok;     // being made: the meshes shown stay
        if (set > 0 && job_start(m, state, slots, n, cell_size, key)) {
            if (job_poll(m, state) != 0) return m->ok;
        }
        // no set or no job: on the CPU
        return extract_now(m, state, slots, n, cell_size, key);
    }
#else
    (void)geometry;
#endif
    return extract_now(m, state, slots, n, cell_size, key);
}

// Extracts the slots on the CPU before it returns: without a device, or for fields the GPU path does not take.
// The mesh becomes the slots' (m->key = key) whatever comes of it.
static bool extract_now(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t n, double cell_size, uint64_t key) {
    m->key = key;
    m->ok  = false;
    drop_gpu_xyz(m, state);

    const md_tick_t t0 = md_tick_now();
    bool ok = false;
    if (n == 0) {
        ok = isomesh::mesh_upload(&m->mesh, nullptr, nullptr, 3, 0, nullptr, 0);
        set_ranges(m, 0, nullptr, nullptr);
        keep_xyz(m, nullptr, 3, 0);
        keep_index(m, nullptr, 0);
    } else {
        ok = extract_cpu(m, slots, n, cell_size);
    }
    m->seconds = md_tick_to_milliseconds(md_tick_now() - t0) * 1.0e-3;
    m->ok = ok;
    m->on_gpu = false;
    if (ok) {
        MD_LOG_DEBUG("Isosurface meshes: %u slots, %u vertices, %u triangles in %.1f ms (CPU)", n, m->mesh.vertex_count, m->mesh.index_count / 3,
                     m->seconds * 1.0e3);
    } else {
        MD_LOG_ERROR("Isosurface meshes: the extraction failed");
    }
    return ok;
}

// ---------------------------------------------------------------------------------------------------------------
// Electronic structure representations
// ---------------------------------------------------------------------------------------------------------------

// What the meshes hold of their fields (the set, the caches, a job), let go of before the fields are: the mesh
// shown stays until one of the new fields is extracted
static void iso_meshes_drop_fields(IsoMeshes* m, ApplicationState* state) {
    release_extraction_state(m, state);
}

// Whether a mesh is being made: its fields' set on a worker thread, or a job on the GPU
static bool iso_meshes_busy(const IsoMeshes* m) {
    return m->job || build_running((const IsoBuild*)m->set_build);
}

struct EsIsoMesh {
    uint64_t field_key = 0;
    md_gto_iso_field_t* field = nullptr;
    IsoFieldSpec field_spec = {};   // what it is of
    bool field_failed = false;
    vec3_t*  field_xyz = nullptr;   // where the field's atoms are (bohr): what its set on the device is made at
    size_t   field_num_atoms = 0;
    IsoBuild* field_build = nullptr;    // the field of another key, being made on a worker thread
    IsoMeshes meshes = {};

    // A field at the vertices (SurfaceColoring::Field)
    SurfaceFieldPoints field_points = {};
    uint64_t values_key = 0;        // of the mesh and the field the values are of
    ColorScaleSpan span = {};
};

static bool field_spec_equal(const IsoFieldSpec& a, const IsoFieldSpec& b) {
    return a.source == b.source && a.spin == b.spin && a.nto_component == b.nto_component &&
           a.transition_density_component == b.transition_density_component && a.orbital_idx == b.orbital_idx &&
           a.excited_state_idx == b.excited_state_idx && a.nto_lambda_idx == b.nto_lambda_idx &&
           a.density_property_key == b.density_property_key;
}

static IsoFieldSpec field_spec_of(const ElectronicStructureRepresentation& es) {
    IsoFieldSpec f;
    f.source = es.source;
    f.spin = es.spin;
    f.nto_component = es.nto_component;
    f.transition_density_component = es.transition_density_component;
    f.orbital_idx = es.orbital_idx;
    f.excited_state_idx = es.excited_state_idx;
    f.nto_lambda_idx = es.nto_lambda_idx;
    f.density_property_key = es.density_property_key;
    return f;
}

// The field of a representation at the positions of g. Reads only the system's attributes and g: a pool task
// makes it (IsoBuild).
static md_gto_iso_field_t* field_make(const IsoFieldGeometry& g, const IsoFieldSpec& es) {
    switch (es.source) {
    case ElectronicStructureSource::MolecularOrbital: {
        const str_t path = (es.spin == ElectronicStructureSpin::Beta) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.orbital_idx);
        return field_create_orbital(g, path, &slice);
    }
    case ElectronicStructureSource::NaturalTransitionOrbital: {
        const str_t path = (es.nto_component == ElectronicStructureNtoComponent::Particle) ? es_path::nto_particle : es_path::nto_hole;
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)es.excited_state_idx, (uint32_t)es.nto_lambda_idx);
        return field_create_orbital(g, path, &slice);
    }
    case ElectronicStructureSource::TransitionDensity: {
        str_t path = es_path::attachment_density;
        switch (es.transition_density_component) {
        case ElectronicStructureTransitionDensityComponent::Detachment: path = es_path::detachment_density; break;
        case ElectronicStructureTransitionDensityComponent::Difference: path = es_path::transition_diff;    break;
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default: break;
        }
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.excited_state_idx);
        return field_create_density(g, path, &slice);
    }
    case ElectronicStructureSource::ElectronDensity: {
        // from its orbitals where the system has them, from its matrix otherwise. Its maxima are at the
        // nuclei: the GPU extraction evaluates it only around the surface (the flood)
        const str_t path = es_electron_density_path(es.spin);
        md_gto_iso_field_t* f = field_create_scf_density(g, path);
        if (!f) f = field_create_density(g, path, nullptr);
        md_gto_iso_field_set_maxima_at_atoms(f, true);
        return f;
    }
    case ElectronicStructureSource::DensityProperty: {
        const md_attribute_t* attr = md_attributes_get(&g.sys->attributes, es.density_property_key);
        return attr ? field_create_density(g, attr->path, nullptr) : nullptr;
    }
    default:
        return nullptr;
    }
}

// What names the field: what it is of and where the atoms are
static uint64_t electronic_structure_field_key(const ApplicationState* state, const ElectronicStructureRepresentation& es) {
    uint64_t h = 0x150;
    const md_system_state_t& st = state->mold.state;
    if (st.xyz && st.num_atoms) h = md_hash64(st.xyz, sizeof(st.xyz[0]) * st.num_atoms, h);
    h = md_hash64(&es.source, sizeof(es.source), h);
    h = md_hash64(&es.spin, sizeof(es.spin), h);
    h = md_hash64(&es.nto_component, sizeof(es.nto_component), h);
    h = md_hash64(&es.transition_density_component, sizeof(es.transition_density_component), h);
    h = md_hash64(&es.orbital_idx, sizeof(es.orbital_idx), h);
    h = md_hash64(&es.excited_state_idx, sizeof(es.excited_state_idx), h);
    h = md_hash64(&es.nto_lambda_idx, sizeof(es.nto_lambda_idx), h);
    h = md_hash64(&es.density_property_key, sizeof(es.density_property_key), h);
    return h ? h : 1;
}

void electronic_structure_mesh_free(ApplicationState* state, ElectronicStructureRepresentation* es) {
    ASSERT(es);
    if (!es->iso_mesh) return;
    iso_meshes_free(&es->iso_mesh->meshes, state);
    iso_field_release(es->iso_mesh->field);     // when no set of it is being made
    build_orphan(es->iso_mesh->field_build);
    free(es->iso_mesh->field_xyz);
    surface_field_points_free(&es->iso_mesh->field_points);
    delete es->iso_mesh;
    es->iso_mesh = nullptr;
}

// The field of what is wanted now, begun on a worker thread when none is being made, and not while the set of the
// one shown is (so after the meshes' update, which begins that): for a large density each is seconds of the CPU and
// a few hundred MB (the field's orbitals, the set's density matrix), and one at a time keeps both to what one
// extraction needs; the next field is made while the GPU extracts the last
static void field_build_begin(EsIsoMesh& em, ApplicationState* state, uint64_t field_key, const IsoFieldSpec& spec) {
    if (field_key != em.field_key && !em.field_build && !build_running((const IsoBuild*)em.meshes.set_build) &&
        // a natural transition orbital: not until its state has been read (on a worker thread, es_nto_ready)
        (spec.source != ElectronicStructureSource::NaturalTransitionOrbital || es_nto_ready(state, (size_t)spec.excited_state_idx))) {
        md_temp_scope_t temp = md_temp_begin();
        const IsoFieldGeometry g = current_geometry(temp, state);
        IsoBuild* b = build_new(IsoBuild::Field, field_key, g.sys, g.atom_xyz, g.num_atoms);
        md_temp_end(temp);
        if (b) {
            b->spec = spec;
            build_start(b);
            em.field_build = b;
        }
    }
}

size_t electronic_structure_mesh_surfaces(isomesh::Surface* out, size_t cap, ApplicationState* state, Representation* rep, const IsoDesc& iso) {
    ASSERT(state);
    ASSERT(rep);
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    if (!es.iso_mesh) es.iso_mesh = new EsIsoMesh();
    EsIsoMesh& em = *es.iso_mesh;

    // The field, made on a worker thread: the old one and its meshes stay until the new one is there. A field
    // of what is wanted at other positions (a frame since played) is taken all the same, as the newest there is,
    // but not while the mesh of the one shown is being made, if that is of the same: dropping its work for a
    // field that will itself be out of date by the time its mesh is would, on a trajectory played faster than
    // a mesh is made, leave the surface where it was. A field of something else (another orbital) is not taken.
    const uint64_t field_key = electronic_structure_field_key(state, es);
    const IsoFieldSpec spec = field_spec_of(es);
    const bool shown_is_wanted = em.field && field_spec_equal(em.field_spec, spec);
    if (em.field_build && !build_running(em.field_build) &&
        !(em.field_build->ran && field_spec_equal(em.field_build->spec, spec) && shown_is_wanted && iso_meshes_busy(&em.meshes))) {
        IsoBuild* b = em.field_build;
        em.field_build = nullptr;
        build_remove(b);
        if (b->ran && field_spec_equal(b->spec, spec)) {
            // in: the meshes of the old field go with their state on the device (their set holds it)
            iso_meshes_drop_fields(&em.meshes, state);
            iso_field_release(em.field);
            em.field = b->field;
            b->field = nullptr;
            free(em.field_xyz);
            em.field_xyz = b->xyz;
            em.field_num_atoms = b->num_atoms;
            b->xyz = nullptr;
            b->num_atoms = 0;
            em.field_key = b->key;
            em.field_spec = b->spec;
            em.field_failed = em.field == nullptr;
            if (em.field_failed) {
                MD_LOG_ERROR("Isosurface meshes: no field for the representation '%s'", rep->name);
            } else {
                MD_LOG_DEBUG("Isosurface meshes: the field of '%s' made in %.1f ms (on a worker thread)", rep->name, md_tick_to_milliseconds(md_tick_now() - b->t0));
            }
        }
        build_destroy(b);   // of another field, or interrupted before it ran: made again below
    }
    if (!em.field) {
        field_build_begin(em, state, field_key, spec);
        return 0;
    }

    // One slot per distinct isovalue, its sheets those of the entries that share it: the + and - surfaces of
    // an orbital are the two sheets of one slot. A magnitude |f| = v is both sheets of f at v.
    const bool magnitude = electronic_structure_uses_magnitude_toggle(es) && es.use_magnitude;
    IsoMeshSlot slots[ISO_MESHES_MAX_SLOTS];
    uint32_t num_slots = 0;
    int      entry_slot[8];
    uint32_t entry_sheets[8];
    const size_t count = MIN(iso.count, (size_t)8);
    for (size_t i = 0; i < count; ++i) {
        entry_slot[i] = -1;
        const double v = (double)iso.values[i];
        const double a = fabs(v);
        if (!(a > 0.0)) continue;
        const uint32_t sheets = magnitude ? MD_GTO_ISO_SHEET_BOTH : (v > 0.0 ? MD_GTO_ISO_SHEET_POSITIVE : MD_GTO_ISO_SHEET_NEGATIVE);
        int slot = -1;
        for (uint32_t k = 0; k < num_slots; ++k) {
            if (slots[k].iso == a) { slot = (int)k; break; }
        }
        if (slot < 0) {
            slot = (int)num_slots++;
            slots[slot] = {em.field, a, 0};
        }
        slots[slot].sheets |= sheets;
        entry_slot[i] = slot;
        entry_sheets[i] = sheets;
    }

    IsoFieldGeometry geometry;
    geometry.sys = &state->mold.sys;
    geometry.atom_xyz = em.field_xyz;
    geometry.num_atoms = em.field_num_atoms;
    const bool ok = iso_meshes_update(&em.meshes, state, slots, num_slots, iso_meshes_cell_size(es.resolution), em.field_key, &geometry);
    field_build_begin(em, state, field_key, spec);
    if (!ok) return 0;

    // Coloured by a field: the field at the vertices
    bool vertex_field = false;
    if (es.coloring == SurfaceColoring::Field && surface_field_available(es.field_kind, state->mold.sys)) {
        // What the values are of: the field (its kind, the frame) and the vertices
        uint64_t source = md_hash64(&es.field_kind, sizeof(es.field_kind), 0x5f1e1d);
        source = md_hash64(&state->animation.frame, sizeof(state->animation.frame), source);
        const md_system_state_t& st = state->mold.state;
        if (st.xyz && st.num_atoms) source = md_hash64(st.xyz, sizeof(st.xyz[0]) * st.num_atoms, source);
        const uint64_t values_key = md_hash64_combine(source, em.meshes.key);
        if (values_key != em.values_key) {
            // Begun now, taken below once it is done: the surfaces go uncoloured until then
            em.values_key = values_key;
            em.span = {};
            md_temp_scope_t temp = md_temp_begin();
            SurfaceFieldDesc desc = {};
            desc.kind        = es.field_kind;
            desc.sys         = &state->mold.sys;
            desc.state       = &state->mold.state;
            desc.source_hash = source;
            if (surface_field_uses_basis(es.field_kind)) {
                desc.basis_atom_xyz = basis_atom_positions_extract(&desc.num_basis_atoms, temp, state->mold.sys, state->mold.state);
            }
#if MD_ENABLE_GPU
            desc.gpu_stream = state->gpu_stream;
#endif
            SurfaceFieldPointSet points = {};
            points.xyz        = em.meshes.xyz;
            points.gpu_xyz    = em.meshes.gpu_xyz;
            points.gpu_stride = sizeof(float) * 4;
            points.count      = em.meshes.xyz_count;
            if (!surface_field_points_begin(&em.field_points, points, desc)) {
                isomesh::mesh_upload_values(&em.meshes.mesh, nullptr, nullptr, 0);
            }
            md_temp_end(temp);
        }
        const float* values = nullptr;
        size_t nv = 0;
        if (surface_field_points_poll(&em.field_points, &values, &nv) == SurfaceFieldPointsStatus::Done && values && nv == em.meshes.xyz_count) {
            // and its gradient along the surface, which the isolines take their width from
            md_temp_scope_t temp = md_temp_begin();
            float* gradients = (float*)md_temp_alloc(temp, sizeof(float) * 3 * MAX(nv, (size_t)1));
            isomesh::mesh_field_gradient(gradients, em.meshes.xyz, values, (uint32_t)nv, em.meshes.tri_index, em.meshes.tri_index_count);
            if (isomesh::mesh_upload_values(&em.meshes.mesh, values, gradients, (uint32_t)nv)) {
                em.span = surface_field_span_of(values, nv);
                color_scale_update_range(&es.field_map, em.span);
                MD_LOG_DEBUG("Isosurface meshes: %s at %zu vertices in %.1f ms (%s, over frames)", surface_field_kind_info(es.field_kind).label, nv,
                             em.field_points.seconds * 1.0e3, em.field_points.on_gpu ? "GPU" : "CPU");
            } else {
                isomesh::mesh_upload_values(&em.meshes.mesh, nullptr, nullptr, 0);
            }
            md_temp_end(temp);
        }
        vertex_field = em.span.valid;
    }

    const mat4_t model = iso_meshes_model_matrix();
    size_t n = 0;
    for (size_t i = 0; i < count; ++i) {
        if (entry_slot[i] < 0) continue;
        const int k = entry_slot[i];
        for (int s = 0; s < 2; ++s) {
            if (!(entry_sheets[i] & (1u << s)) || em.meshes.index_count[k][s] == 0 || n >= cap) continue;
            isomesh::Surface& surf = out[n++];
            surf = {};
            surf.mesh = &em.meshes.mesh;
            surf.index_offset = em.meshes.index_offset[k][s];
            surf.index_count  = em.meshes.index_count[k][s];
            surf.model = model;
            surf.color = iso.colors[i];
            surf.optical_density = iso.optical_densities[i];
            if (vertex_field) surf.coloring = isomesh::Coloring_VertexField;
        }
    }
    return n;
}

ColorScaleSpan electronic_structure_mesh_field_span(const ElectronicStructureRepresentation& es) {
    return es.iso_mesh ? es.iso_mesh->span : ColorScaleSpan{};
}
