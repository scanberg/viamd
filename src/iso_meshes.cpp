#include <iso_meshes.h>

#include <viamd.h>
#include <task_system.h>
#include <surface_field.h>

#include <md_gto_iso.h>
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

static bool field_context(FieldContext* ctx, md_temp_scope_t temp, const ApplicationState* state) {
    if (!md_gto_basis_extract_attributes(&ctx->basis, &state->mold.sys.attributes, md_temp_allocator(temp))) {
        return false;
    }
    ctx->atom_xyz = basis_atom_positions_extract(&ctx->num_atoms, temp, state->mold.sys, state->mold.state);
    if (!ctx->atom_xyz || ctx->num_atoms != md_gto_basis_num_atoms(&ctx->basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis and the atom positions do not agree");
        return false;
    }
    return true;
}

md_gto_iso_field_t* iso_field_create_orbital(const ApplicationState* state, str_t coefficient_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    FieldContext ctx;
    if (!field_context(&ctx, temp, state)) return nullptr;
    size_t num_ao = 0;
    const double* coeffs = orbital_coefficients_extract(&num_ao, temp, state->mold.sys, coefficient_path, slice);
    if (!coeffs) return nullptr;
    if (num_ao != md_gto_basis_num_ao(&ctx.basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis spans %zu atomic orbitals and '" STR_FMT "' %zu", md_gto_basis_num_ao(&ctx.basis), STR_ARG(coefficient_path), num_ao);
        return nullptr;
    }
    return md_gto_iso_field_create(&ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), coeffs, md_get_heap_allocator());
}

md_gto_iso_field_t* iso_field_create_density_matrix(const ApplicationState* state, const double* density_matrix, size_t dim) {
    if (!density_matrix || dim == 0) return nullptr;
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    FieldContext ctx;
    if (!field_context(&ctx, temp, state)) return nullptr;
    if (dim != md_gto_basis_num_ao(&ctx.basis)) {
        MD_LOG_ERROR("Isosurface meshes: the basis spans %zu atomic orbitals and the density matrix %zu", md_gto_basis_num_ao(&ctx.basis), dim);
        return nullptr;
    }
    return md_gto_iso_field_create_density(&ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), density_matrix, md_get_heap_allocator());
}

md_gto_iso_field_t* iso_field_create_density(const ApplicationState* state, str_t density_path, const md_attribute_slice_t* slice) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    size_t dim = 0;
    const double* D = density_matrix_extract(&dim, temp, state->mold.sys, density_path, slice);
    return D ? iso_field_create_density_matrix(state, D, dim) : nullptr;
}

// ---------------------------------------------------------------------------------------------------------------
// Meshes
// ---------------------------------------------------------------------------------------------------------------

static void release_extraction_state(IsoMeshes* m, ApplicationState* state) {
#if MD_ENABLE_GPU
    if (m->gpu_set) {
        md_gto_iso_gpu_set_destroy(state->gpu_stream, (md_gto_iso_gpu_set_t*)m->gpu_set);
    }
#else
    (void)state;
#endif
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
static bool extract_gpu(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t n, double cell_size) {
    if (!state->gpu_device || !state->gpu_stream) return false;
    md_gpu_stream_t stream = state->gpu_stream;

    if (!m->gpu_set) {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        FieldContext ctx;
        if (!field_context(&ctx, temp, state)) return false;
        const md_gto_iso_field_t* fields[ISO_MESHES_MAX_SLOTS];
        for (uint32_t k = 0; k < n; ++k) fields[k] = slots[k].field;
        m->gpu_set = md_gto_iso_gpu_set_create(stream, &ctx.basis, (const float*)ctx.atom_xyz, sizeof(vec3_t), fields, n, md_get_heap_allocator());
        if (!m->gpu_set) {
            MD_LOG_DEBUG("Isosurface meshes: no GPU set for these fields, extracting on the CPU");
            return false;
        }
    }

    md_gto_iso_gpu_desc_t desc = {};
    desc.cell_size = cell_size;
    desc.project_iterations = PROJECT_ITERATIONS;
    desc.relax_iterations = RELAX_ITERATIONS;
    for (uint32_t k = 0; k < n; ++k) {
        desc.iso[k] = slots[k].iso;
        desc.sheets[k] = slots[k].sheets;
    }

    md_gto_iso_gpu_mesh_t gm = {};
    if (!md_gto_iso_gpu_extract(stream, (const md_gto_iso_gpu_set_t*)m->gpu_set, &desc, &gm)) {
        return false;
    }

    bool ok = false;
    const size_t nv = gm.vertex_count, ni = gm.index_count;
    float*    pos = (float*)malloc(sizeof(float) * 4 * (nv + 1));
    float*    nrm = (float*)malloc(sizeof(float) * 4 * (nv + 1));
    uint32_t* idx = (uint32_t*)malloc(sizeof(uint32_t) * (ni + 1));
    if (pos && nrm && idx && md_gto_iso_gpu_mesh_read_all(stream, &gm, pos, nrm, idx)) {
        ok = isomesh::mesh_upload(&m->mesh, pos, nrm, 4, gm.vertex_count, idx, gm.index_count);
        if (ok) {
            set_ranges(m, n, gm.index_offset, gm.index_count_of);
            keep_xyz(m, pos, 4, gm.vertex_count);
            keep_index(m, idx, gm.index_count);
            // The positions stay on the device for what is evaluated at them; the normals and indices
            // are in GL now
            m->gpu_xyz = gm.position;
            gm.position = 0;
        }
    }
    free(pos);
    free(nrm);
    free(idx);
    md_gto_iso_gpu_mesh_free(stream, &gm);
    return ok;
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

bool iso_meshes_update(IsoMeshes* m, ApplicationState* state, const IsoMeshSlot* slots, uint32_t num_slots, double cell_size, uint64_t fields_key) {
    ASSERT(m);
    ASSERT(state);
    const uint32_t n = MIN(num_slots, ISO_MESHES_MAX_SLOTS);

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
    if (key == m->key) return m->ok;
    m->key = key;
    m->ok  = false;
    drop_gpu_xyz(m, state);

    for (uint32_t k = 0; k < n; ++k) {
        if (!slots[k].field || !(slots[k].iso > 0.0)) return false;
    }
    if (fkey != m->fields_key) {
        release_extraction_state(m, state);
        m->fields_key = fkey;
    }

    const md_tick_t t0 = md_tick_now();
    bool ok = false;
    bool on_gpu = false;
    if (n == 0) {
        ok = isomesh::mesh_upload(&m->mesh, nullptr, nullptr, 3, 0, nullptr, 0);
        set_ranges(m, 0, nullptr, nullptr);
        keep_xyz(m, nullptr, 3, 0);
        keep_index(m, nullptr, 0);
    } else {
#if MD_ENABLE_GPU
        ok = on_gpu = extract_gpu(m, state, slots, n, cell_size);
#endif
        if (!ok) {
            ok = extract_cpu(m, slots, n, cell_size);
        }
    }
    m->seconds = md_tick_to_milliseconds(md_tick_now() - t0) * 1.0e-3;
    m->ok = ok;
    m->on_gpu = on_gpu;
    if (ok) {
        MD_LOG_DEBUG("Isosurface meshes: %u slots, %u vertices, %u triangles in %.1f ms (%s)", n, m->mesh.vertex_count, m->mesh.index_count / 3,
                     m->seconds * 1.0e3, on_gpu ? "GPU" : "CPU");
    } else {
        MD_LOG_ERROR("Isosurface meshes: the extraction failed");
    }
    return ok;
}

// ---------------------------------------------------------------------------------------------------------------
// Electronic structure representations
// ---------------------------------------------------------------------------------------------------------------

struct EsIsoMesh {
    uint64_t field_key = 0;
    md_gto_iso_field_t* field = nullptr;
    bool field_failed = false;
    IsoMeshes meshes = {};

    // A field at the vertices (SurfaceColoring::Field)
    SurfaceFieldPoints field_points = {};
    uint64_t values_key = 0;        // of the mesh and the field the values are of
    ColorScaleSpan span = {};
};

static md_gto_iso_field_t* electronic_structure_field_create(const ApplicationState* state, const ElectronicStructureRepresentation& es) {
    switch (es.source) {
    case ElectronicStructureSource::MolecularOrbital: {
        const str_t path = (es.spin == ElectronicStructureSpin::Beta) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.orbital_idx);
        return iso_field_create_orbital(state, path, &slice);
    }
    case ElectronicStructureSource::NaturalTransitionOrbital: {
        const str_t path = (es.nto_component == ElectronicStructureNtoComponent::Particle) ? es_path::nto_particle : es_path::nto_hole;
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)es.excited_state_idx, (uint32_t)es.nto_lambda_idx);
        return iso_field_create_orbital(state, path, &slice);
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
        return iso_field_create_density(state, path, &slice);
    }
    case ElectronicStructureSource::ElectronDensity:
        return iso_field_create_density(state, es_electron_density_path(es.spin), nullptr);
    case ElectronicStructureSource::DensityProperty: {
        const md_attribute_t* attr = md_attributes_get(&state->mold.sys.attributes, es.density_property_key);
        return attr ? iso_field_create_density(state, attr->path, nullptr) : nullptr;
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
    md_gto_iso_field_destroy(es->iso_mesh->field);
#if MD_ENABLE_GPU
    surface_field_points_free(&es->iso_mesh->field_points, state->gpu_stream);
#else
    surface_field_points_free(&es->iso_mesh->field_points, nullptr);
#endif
    delete es->iso_mesh;
    es->iso_mesh = nullptr;
}

size_t electronic_structure_mesh_surfaces(isomesh::Surface* out, size_t cap, ApplicationState* state, Representation* rep, const IsoDesc& iso) {
    ASSERT(state);
    ASSERT(rep);
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    if (!es.iso_mesh) es.iso_mesh = new EsIsoMesh();
    EsIsoMesh& em = *es.iso_mesh;

    const uint64_t field_key = electronic_structure_field_key(state, es);
    if (field_key != em.field_key) {
        // the set holds the old field: dropped before the field goes
        iso_meshes_free(&em.meshes, state);
        md_gto_iso_field_destroy(em.field);
        em.field = electronic_structure_field_create(state, es);
        em.field_key = field_key;
        em.field_failed = em.field == nullptr;
        if (em.field_failed) {
            MD_LOG_ERROR("Isosurface meshes: no field for the representation '%s'", rep->name);
        }
    }
    if (!em.field) return 0;

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

    if (!iso_meshes_update(&em.meshes, state, slots, num_slots, iso_meshes_cell_size(es.resolution), field_key)) {
        return 0;
    }

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
            em.values_key = values_key;
            em.span = {};
            md_temp_scope_t temp = md_temp_begin();
            SurfaceFieldDesc desc = {};
            desc.kind        = es.field_kind;
            desc.sys         = &state->mold.sys;
            desc.state       = &state->mold.state;
            desc.source_hash = source;
            if (es.field_kind == SurfaceFieldKind::ElectrostaticPotential) {
                desc.basis_atom_xyz = basis_atom_positions_extract(&desc.num_basis_atoms, temp, state->mold.sys, state->mold.state);
            }
#if MD_ENABLE_GPU
            desc.gpu_stream = state->gpu_stream;
#endif
            const uint32_t nv = em.meshes.xyz_count;
            float* values = (float*)md_temp_alloc(temp, sizeof(float) * MAX(nv, 1u));
            SurfaceFieldPointSet points = {};
            points.xyz        = em.meshes.xyz;
            points.gpu_xyz    = em.meshes.gpu_xyz;
            points.gpu_stride = sizeof(float) * 4;
            points.count      = nv;
            // and its gradient along the surface, which the isolines take their width from
            float* gradients = (float*)md_temp_alloc(temp, sizeof(float) * 3 * MAX(nv, 1u));
            bool evaluated = surface_field_evaluate_points(&em.field_points, values, points, desc);
            if (evaluated) {
                isomesh::mesh_field_gradient(gradients, em.meshes.xyz, values, nv, em.meshes.tri_index, em.meshes.tri_index_count);
            }
            if (evaluated && isomesh::mesh_upload_values(&em.meshes.mesh, values, gradients, nv)) {
                em.span = surface_field_span_of(values, nv);
                color_scale_update_range(&es.field_map, em.span);
                MD_LOG_DEBUG("Isosurface meshes: %s at %u vertices in %.1f ms (%s)", surface_field_kind_str[(int)es.field_kind], nv,
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
