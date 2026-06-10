#if VIAMD_TREXIO

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>

#define IMGUI_DEFINE_MATH_OPERATORS

#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/volumerender_utils.h>

#include <md_gto.h>
#include <md_util.h>

#include <viamd.h>
#include <viamd_event.h>

#include <color_utils.h>
#include <imgui_widgets.h>
#include <loader.h>
#include <qm_ui.h>
#include <trexio_data.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

#include <imgui_internal.h>

#define TREXIO_ANGSTROM_TO_BOHR 1.8897261246257702
#define TREXIO_BOHR_TO_ANGSTROM 0.529177210903
#define TREXIO_DEFAULT_GTO_CUTOFF_VALUE 1.0e-6

static const double trexio_volume_resolution_samples_per_angstrom[3] = {4.0, 8.0, 16.0};
static const double trexio_orbital_grid_samples_per_angstrom = 8.0;
constexpr uint64_t interaction_surface_trexio_orb = HASH_STR_LIT64("interaction surface trexio orb");

struct TrexioAABB {
    vec3_t min_ext = {0};
    vec3_t max_ext = {0};
};

struct TrexioShellEntry {
    const trexio_data::Shell* shell = nullptr;
    uint32_t atom_idx = 0;
    uint32_t l = 0;
    uint32_t source_ao_offset = 0;
    uint32_t source_func_count = 0;
    uint32_t source_order = 0;
};

enum class TrexioAoOrderingMode : uint8_t {
    TrexioSpherical,
    SourceOrder,
    SourcePxyz,
};

static const char* trexio_ao_ordering_mode_str(TrexioAoOrderingMode mode) {
    switch (mode) {
    case TrexioAoOrderingMode::TrexioSpherical: return "TREXIO spherical";
    case TrexioAoOrderingMode::SourceOrder: return "Source / mdlib order";
    case TrexioAoOrderingMode::SourcePxyz: return "Source p: px,py,pz";
    default: return "Unknown";
    }
}

static inline md_gto_op_t trexio_gto_op_from_use_magnitude(bool use_magnitude) {
    return MD_GTO_OP_SET | (use_magnitude ? MD_GTO_OP_ABS : 0);
}

static inline mat4_t trexio_compute_texture_to_world_mat(const mat3_t& orientation, const vec3_t& origin, const vec3_t& extent) {
    mat4_t T = mat4_translate_vec3(origin);
    mat4_t R = mat4_from_mat3(orientation);
    mat4_t S = mat4_scale_vec3(extent);
    return T * R * S;
}

static inline mat4_t trexio_compute_world_to_model_mat(const mat3_t& orientation, const vec3_t& origin) {
    return mat4_from_mat3(mat3_transpose(orientation)) * mat4_translate_vec3(-origin);
}

static inline void trexio_compute_dim(int out_dim[3], const vec3_t& extent, double samples_per_unit_length) {
    out_dim[0] = CLAMP(ALIGN_TO((int)(extent.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = CLAMP(ALIGN_TO((int)(extent.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = CLAMP(ALIGN_TO((int)(extent.z * samples_per_unit_length), 8), 8, 512);
}

static void trexio_init_grid(md_grid_t* grid, const vec3_t& min_ext, const vec3_t& max_ext, double samples_per_unit_length) {
    vec3_t extent = max_ext - min_ext;
    trexio_compute_dim(grid->dim, extent, samples_per_unit_length);
    grid->orientation = mat3_ident();
    grid->origin = min_ext;
    grid->spacing = vec3_div(extent, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
}

static void trexio_init_volume(Volume* vol, const md_grid_t& grid, GLenum format = GL_R32F) {
    ASSERT(vol);
    MEMCPY(vol->dim, grid.dim, sizeof(vol->dim));
    const float scl = TREXIO_BOHR_TO_ANGSTROM;
    vec3_t extent = md_grid_extent(&grid);
    vol->world_to_model = trexio_compute_world_to_model_mat(grid.orientation, grid.origin * scl);
    vol->texture_to_world = trexio_compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
    vol->voxel_size = grid.spacing * scl;
    gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
}

static void trexio_calculate_bounds(TrexioAABB* aabb, const vec4_t* atom_xyzw, size_t count) {
    ASSERT(aabb);
    ASSERT(atom_xyzw);

    vec4_t min_v = vec4_set1(FLT_MAX);
    vec4_t max_v = vec4_set1(-FLT_MAX);
    for (size_t i = 0; i < count; ++i) {
        min_v = vec4_min(min_v, atom_xyzw[i]);
        max_v = vec4_max(max_v, atom_xyzw[i]);
    }

    const float pad = 6.0f;
    min_v -= pad;
    max_v += pad;
    aabb->min_ext = vec3_from_vec4(min_v);
    aabb->max_ext = vec3_from_vec4(max_v);
}

static bool trexio_orbital_matches_spin(const trexio_data::Orbital& orbital, ElectronicStructureSpin spin) {
    if (spin == ElectronicStructureSpin::Beta) {
        return orbital.spin == 1;
    }
    return orbital.spin != 1;
}

static size_t trexio_shell_num_spherical_functions(int l) {
    return (size_t)(2 * l + 1);
}

static size_t trexio_shell_component_source_index(int l, TrexioAoOrderingMode mode, size_t component_idx) {
    if (l > 0) {
        ASSERT(component_idx < trexio_shell_num_spherical_functions(l));
        switch (mode) {
        case TrexioAoOrderingMode::TrexioSpherical: {
            // TREXIO spherical AO order is m = 0,+1,-1,+2,-2,...
            // mdlib/VIAMD shell order is m = -l,...,-1,0,+1,...,+l.
            const int m = (int)component_idx - l;
            if (m == 0) return 0;
            return m > 0 ? (size_t)(2 * m - 1) : (size_t)(-2 * m);
        }
        case TrexioAoOrderingMode::SourcePxyz:
            if (l == 1) {
                // Source p order is px,py,pz; mdlib/VIAMD p order is py,pz,px.
                static const size_t p_xyz_to_mdlib[] = {1, 2, 0};
                return p_xyz_to_mdlib[component_idx];
            }
            return component_idx;
        case TrexioAoOrderingMode::SourceOrder:
        default:
            return component_idx;
        }
    }
    return component_idx;
}

struct TrexioComponent : viamd::EventHandler {
    struct Settings {
        TrexioAoOrderingMode ao_ordering_mode = TrexioAoOrderingMode::TrexioSpherical;
    } settings;

    struct Summary {
        bool show_window = false;
    } summary;

    struct OrbitalGrid {
        bool show_window = false;
        Volume vol[16] = {};
        int vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        ElectronicStructureSpin vol_mo_spin[16] = {};
        uint32_t iso_tex[16] = {};
        GBuffer gbuf = {};
        PickingSurface picking_surface = {};
        Camera camera = {};
        ViewTransform target = {};
        ViewTransform default_view = {};
        bool show_coordinate_system_widget = true;
        int num_x = 3;
        int num_y = 3;
        int mo_idx = 0;
        int scroll_to_idx = -1;
        IsoDesc iso = {
            .count = 2,
            .values = {0.05f, -0.05f},
            .colors = {{0.f/255.f,75.f/255.f,135.f/255.f,0.75f}, {255.f/255.f,205.f/255.f,0.f/255.f,0.75f}},
        };
        md_gl_rep_t gl_rep = {};
    } orb;

    md_allocator_i* arena = nullptr;
    trexio_data::TrexioData data = {};
    md_gto_basis_t basis = {};
    md_array(uint32_t) ordered_to_source_ao = nullptr;
    md_array(int32_t) alpha_orbital_map = nullptr;
    md_array(int32_t) beta_orbital_map = nullptr;
    md_array(double) alpha_occupation_cache = nullptr;
    md_array(double) beta_occupation_cache = nullptr;
    md_array(double) alpha_energy_cache = nullptr;
    md_array(double) beta_energy_cache = nullptr;
    bool loaded = false;
    bool basis_ready = false;

    TrexioComponent() {
        viamd::event_system_register_handler(*this);
    }

    void rebuild_basis_mapping() {
        basis_ready = build_basis();
        loaded = basis_ready;
        invalidate_orbital_caches();
    }

    void on_ao_ordering_mode_changed(ApplicationState& state) {
        rebuild_basis_mapping();
        update_representation_info(&state);
        flag_all_representations_as_dirty(&state);
    }

    void invalidate_orbital_caches() {
        for (int i = 0; i < (int)ARRAY_SIZE(orb.vol_mo_idx); ++i) {
            orb.vol_mo_idx[i] = -1;
            orb.vol_mo_spin[i] = {};
        }
    }

    void reset() {
        for (Volume& vol : orb.vol) {
            if (vol.tex_id) {
                glDeleteTextures(1, &vol.tex_id);
                vol.tex_id = 0;
            }
        }
        for (uint32_t& tex : orb.iso_tex) {
            if (tex) {
                glDeleteTextures(1, &tex);
                tex = 0;
            }
        }
        destroy_gbuffer(&orb.gbuf);
        picking_surface_free(&orb.picking_surface);
        md_gl_rep_destroy(orb.gl_rep);
        settings = {};
        orb = {};
        data = {};
        basis = {};
        ordered_to_source_ao = nullptr;
        alpha_orbital_map = nullptr;
        beta_orbital_map = nullptr;
        alpha_occupation_cache = nullptr;
        beta_occupation_cache = nullptr;
        alpha_energy_cache = nullptr;
        beta_energy_cache = nullptr;
        loaded = false;
        basis_ready = false;
        if (arena) {
            md_arena_allocator_reset(arena);
        }
    }

    void ensure_arena() {
        if (!arena) {
            arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
        }
    }

    size_t num_orbitals(ElectronicStructureSpin spin) const {
        return spin == ElectronicStructureSpin::Beta ? md_array_size(beta_orbital_map) : md_array_size(alpha_orbital_map);
    }

    bool has_beta_spin() const {
        return md_array_size(beta_orbital_map) > 0;
    }

    const double* orbital_occupancy(ElectronicStructureSpin spin) const {
        return spin == ElectronicStructureSpin::Beta ? beta_occupation_cache : alpha_occupation_cache;
    }

    const double* orbital_energy(ElectronicStructureSpin spin) const {
        return spin == ElectronicStructureSpin::Beta ? beta_energy_cache : alpha_energy_cache;
    }

    size_t find_homo_idx(ElectronicStructureSpin spin) const {
        size_t homo_idx = 0;
        bool found = false;
        size_t local_idx = 0;
        for (const auto& orbital : data.orbitals) {
            if (!trexio_orbital_matches_spin(orbital, spin)) continue;
            if (orbital.occupation > 0.0) {
                homo_idx = local_idx;
                found = true;
            }
            ++local_idx;
        }
        return found ? homo_idx : 0;
    }

    int homo_idx(ElectronicStructureSpin spin) const {
        return (int)find_homo_idx(spin);
    }

    int lumo_idx(ElectronicStructureSpin spin) const {
        const size_t count = num_orbitals(spin);
        return count ? MIN(homo_idx(spin) + 1, (int)count - 1) : -1;
    }

    void build_orbital_maps() {
        alpha_orbital_map = nullptr;
        beta_orbital_map = nullptr;
        alpha_occupation_cache = nullptr;
        beta_occupation_cache = nullptr;
        alpha_energy_cache = nullptr;
        beta_energy_cache = nullptr;

        for (size_t i = 0; i < data.orbitals.size(); ++i) {
            if (data.orbitals[i].spin == 1) {
                md_array_push(beta_orbital_map, (int32_t)i, arena);
                md_array_push(beta_occupation_cache, data.orbitals[i].occupation, arena);
                md_array_push(beta_energy_cache, data.orbitals[i].energy, arena);
            } else {
                md_array_push(alpha_orbital_map, (int32_t)i, arena);
                md_array_push(alpha_occupation_cache, data.orbitals[i].occupation, arena);
                md_array_push(alpha_energy_cache, data.orbitals[i].energy, arena);
            }
        }
    }

    bool build_basis() {
        basis = {};
        ordered_to_source_ao = nullptr;

        if (data.shells.empty() || data.orbitals.empty() || data.ao_num <= 0) {
            return false;
        }

        std::vector<uint32_t> shell_source_offset(data.shells.size(), UINT32_MAX);
        std::vector<uint32_t> shell_source_count(data.shells.size(), 0);

        if (data.ao_shell.size() == (size_t)data.ao_num) {
            for (int32_t ao = 0; ao < data.ao_num; ++ao) {
                const int32_t shell_idx = data.ao_shell[(size_t)ao];
                if (shell_idx < 0 || shell_idx >= (int32_t)data.shells.size()) continue;
                uint32_t& offset = shell_source_offset[(size_t)shell_idx];
                if (offset == UINT32_MAX) offset = (uint32_t)ao;
                shell_source_count[(size_t)shell_idx] += 1;
            }
        } else {
            uint32_t offset = 0;
            for (size_t i = 0; i < data.shells.size(); ++i) {
                const int l = data.shells[i].angular_momentum;
                shell_source_offset[i] = offset;
                shell_source_count[i] = (uint32_t)trexio_shell_num_spherical_functions(l);
                offset += shell_source_count[i];
            }
        }

        md_array(TrexioShellEntry) shell_entries = nullptr;
        uint32_t source_order = 0;
        for (size_t i = 0; i < data.shells.size(); ++i) {
            const trexio_data::Shell& src_shell = data.shells[i];
            const int l = src_shell.angular_momentum;
            if (l < 0) {
                MD_LOG_ERROR("Unsupported TREXIO shell angular momentum %d", l);
                return false;
            }
            if (l >= 2 && data.cartesian != 0) {
                MD_LOG_ERROR("Cartesian TREXIO shells with l >= 2 are not supported yet");
                return false;
            }
            if (src_shell.atom_index < 0 || src_shell.atom_index >= (int32_t)data.atoms.size()) {
                MD_LOG_ERROR("TREXIO basis references unknown nucleus index %d", src_shell.atom_index);
                return false;
            }
            if (src_shell.primitives.empty()) {
                continue;
            }

            const size_t expected_count = trexio_shell_num_spherical_functions(l);
            const uint32_t source_count = shell_source_count[i];
            const uint32_t source_offset = shell_source_offset[i];
            if (source_offset == UINT32_MAX || source_count < expected_count) {
                MD_LOG_ERROR("TREXIO AO metadata mismatch for shell %zu", i);
                return false;
            }

            TrexioShellEntry entry = {
                .shell = &src_shell,
                .atom_idx = (uint32_t)src_shell.atom_index,
                .l = (uint32_t)l,
                .source_ao_offset = source_offset,
                .source_func_count = source_count,
                .source_order = source_order++,
            };
            md_array_push(shell_entries, entry, arena);
        }

        std::sort(shell_entries, shell_entries + md_array_size(shell_entries), [](const TrexioShellEntry& a, const TrexioShellEntry& b) {
            if (a.l != b.l) return a.l < b.l;
            if (a.atom_idx != b.atom_idx) return a.atom_idx < b.atom_idx;
            return a.source_order < b.source_order;
        });

        for (size_t i = 0; i < md_array_size(shell_entries); ++i) {
            const TrexioShellEntry& entry = shell_entries[i];
            md_gto_shell_t shell = {
                .atom_idx = entry.atom_idx,
                .primitive_offset = basis.num_primitives,
                .num_primitives = (uint32_t)entry.shell->primitives.size(),
                .l = entry.l,
            };
            md_array_push(basis.shells, shell, arena);
            basis.num_shells += 1;

            for (const auto& primitive : entry.shell->primitives) {
                md_array_push(basis.alpha, (float)primitive.exponent, arena);
                md_array_push(basis.coeff, (float)primitive.coefficient, arena);
                basis.num_primitives += 1;
            }

            const size_t num_funcs = trexio_shell_num_spherical_functions((int)entry.l);
            for (size_t func_idx = 0; func_idx < num_funcs; ++func_idx) {
                const size_t source_component_idx = trexio_shell_component_source_index((int)entry.l, settings.ao_ordering_mode, func_idx);
                if (source_component_idx >= entry.source_func_count) {
                    MD_LOG_ERROR("TREXIO AO component index out of range");
                    return false;
                }
                md_array_push(ordered_to_source_ao, entry.source_ao_offset + (uint32_t)source_component_idx, arena);
            }
        }

        for (const auto& orbital : data.orbitals) {
            if (orbital.coefficients.size() != (size_t)data.ao_num) {
                MD_LOG_ERROR("TREXIO MO coefficient count mismatch: expected %d, got %zu", data.ao_num, orbital.coefficients.size());
                return false;
            }
        }

        if (md_array_size(ordered_to_source_ao) != (size_t)data.ao_num) {
            MD_LOG_ERROR("TREXIO AO remap size mismatch: expected %d, got %zu", data.ao_num, md_array_size(ordered_to_source_ao));
            return false;
        }

        return basis.num_shells > 0;
    }

    bool load_file(str_t path) {
        ensure_arena();
        reset();

        std::string error;
        data = trexio_data::parse_trexio_file(std::string(path.ptr, path.len), &error);
        if (data.atoms.empty()) {
            if (!error.empty()) {
                MD_LOG_ERROR("Failed to load TREXIO data from '%.*s': %s", (int)path.len, path.ptr, error.c_str());
            }
            return false;
        }

        build_orbital_maps();
        basis_ready = build_basis();
        loaded = basis_ready;
        if (!basis_ready) {
            MD_LOG_ERROR("Failed to initialize TREXIO orbital basis from '%.*s'", (int)path.len, path.ptr);
        }
        return loaded;
    }

    const int32_t* orbital_map_for_spin(ElectronicStructureSpin spin, size_t* count) const {
        const int32_t* map = spin == ElectronicStructureSpin::Beta ? beta_orbital_map : alpha_orbital_map;
        if (count) {
            *count = spin == ElectronicStructureSpin::Beta ? md_array_size(beta_orbital_map) : md_array_size(alpha_orbital_map);
        }
        return map;
    }

    const trexio_data::Orbital* orbital_from_local_idx(ElectronicStructureSpin spin, int orbital_idx) const {
        size_t map_count = 0;
        const int32_t* orbital_map = orbital_map_for_spin(spin, &map_count);
        if (!orbital_map || orbital_idx < 0 || orbital_idx >= (int)map_count) return nullptr;
        const int32_t global_idx = orbital_map[orbital_idx];
        if (global_idx < 0 || global_idx >= (int32_t)data.orbitals.size()) return nullptr;
        return &data.orbitals[(size_t)global_idx];
    }

    const trexio_data::Orbital* selected_orbital(const ElectronicStructureRepresentation& es) const {
        return orbital_from_local_idx(es.spin, es.orbital_idx);
    }

    void reorder_coefficients(double* dst_coeffs, const trexio_data::Orbital& orbital) const {
        const size_t count = md_array_size(ordered_to_source_ao);
        for (size_t i = 0; i < count; ++i) {
            dst_coeffs[i] = orbital.coefficients[ordered_to_source_ao[i]];
        }
    }

    static vec4_t* build_atom_positions_bohr(const md_system_t* sys, md_allocator_i* alloc) {
        ASSERT(sys);
        ASSERT(alloc);

        vec4_t* atom_xyzw = nullptr;
        md_array_resize(atom_xyzw, sys->atom.count, alloc);
        for (size_t i = 0; i < sys->atom.count; ++i) {
            atom_xyzw[i] = vec4_set(
                (float)(sys->atom.x[i] * TREXIO_ANGSTROM_TO_BOHR),
                (float)(sys->atom.y[i] * TREXIO_ANGSTROM_TO_BOHR),
                (float)(sys->atom.z[i] * TREXIO_ANGSTROM_TO_BOHR),
                1.0f);
        }
        return atom_xyzw;
    }

    void ensure_visual_state(const ApplicationState& state) {
        if (orb.gl_rep.id) return;
        if (!state.mold.gl_mol.id || state.mold.sys.atom.count == 0) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        uint32_t* colors = md_temp_alloc_array(temp, uint32_t, state.mold.sys.atom.count);
        color_atoms_type(colors, state.mold.sys.atom.count, state.mold.sys);

        reset_view(&orb.default_view, state.mold.sys);
        orb.target = orb.default_view;
        orb.camera = orb.target;
        orb.mo_idx = homo_idx(ElectronicStructureSpin::Alpha);
        orb.scroll_to_idx = orb.mo_idx;

        orb.gl_rep = md_gl_rep_create(state.mold.gl_mol);
        md_gl_rep_set_atom_colors(orb.gl_rep, 0, (uint32_t)state.mold.sys.atom.count, colors, 0);
        picking_surface_init(&orb.picking_surface, interaction_surface_trexio_orb);
    }

    bool evaluate_orbital_texture(uint32_t tex_id, const md_grid_t& grid, const vec4_t* atom_xyzw, ElectronicStructureSpin spin, int orbital_idx) const {
        const trexio_data::Orbital* orbital = orbital_from_local_idx(spin, orbital_idx);
        if (!orbital) return false;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        double* ao_coeffs = md_temp_alloc_array(temp, double, md_array_size(ordered_to_source_ao));
        reorder_coefficients(ao_coeffs, *orbital);

        md_gto_grid_evaluate_mo_GL(tex_id, &grid, &basis, (const float*)atom_xyzw, sizeof(vec4_t), ao_coeffs, TREXIO_DEFAULT_GTO_CUTOFF_VALUE, MD_GTO_EVAL_MODE_PSI, MD_GTO_OP_SET);
        return true;
    }

    static void fill_orbitals(::MolecularOrbital* dst, const TrexioComponent& component, ElectronicStructureSpin spin, md_allocator_i* alloc) {
        ASSERT(dst);
        ASSERT(alloc);

        const size_t num_orbitals = component.num_orbitals(spin);
        if (num_orbitals == 0) {
            *dst = {};
            return;
        }

        dst->num_orbitals = num_orbitals;
        dst->homo_idx = component.find_homo_idx(spin);
        dst->lumo_idx = MIN(dst->homo_idx + 1, num_orbitals - 1);

        md_array_resize(dst->label, num_orbitals, alloc);
        md_array_resize(dst->occupation, num_orbitals, alloc);
        md_array_resize(dst->energy, num_orbitals, alloc);

        size_t count = 0;
        const int32_t* map = component.orbital_map_for_spin(spin, &count);
        for (size_t i = 0; i < count; ++i) {
            const trexio_data::Orbital& orbital = component.data.orbitals[(size_t)map[i]];
            if (!orbital.symmetry.empty()) {
                dst->label[i] = str_printf(alloc, "%zu %s", i + 1, orbital.symmetry.c_str());
            } else {
                dst->label[i] = str_printf(alloc, "%zu", i + 1);
            }
            dst->occupation[i] = orbital.occupation;
            dst->energy[i] = orbital.energy;
        }
    }

    void draw_summary_window(ApplicationState& state) {
        if (!summary.show_window || !loaded) return;

        ImGui::SetNextWindowSize({380, 420}, ImGuiCond_FirstUseEver);
        if (!qm_ui::begin_summary_window("TREXIO", "TrexioQMSummary", &summary.show_window)) {
            ImGui::End();
            return;
        }

        if (ImGui::TreeNode("File Information")) {
            qm_ui::draw_source_label("TREXIO");
            ImGui::Text("AO Ordering:          %s", trexio_ao_ordering_mode_str(settings.ao_ordering_mode));
            if (ImGui::TreeNode("Source Details / Advanced")) {
                ImGui::TextWrapped("Path: %s", data.path.c_str());
                ImGui::Text("Basis Type:           %s", data.basis_type.empty() ? "<none>" : data.basis_type.c_str());
                ImGui::Text("MO Type:              %s", data.mo_type.empty() ? "<none>" : data.mo_type.c_str());
                ImGui::Text("AO Representation:    %s", data.cartesian ? "Cartesian" : "Spherical");
                TrexioAoOrderingMode ao_mode = settings.ao_ordering_mode;
                if (ImGui::BeginCombo("AO Ordering", trexio_ao_ordering_mode_str(ao_mode))) {
                    const TrexioAoOrderingMode modes[] = {
                        TrexioAoOrderingMode::TrexioSpherical,
                        TrexioAoOrderingMode::SourceOrder,
                        TrexioAoOrderingMode::SourcePxyz,
                    };
                    for (TrexioAoOrderingMode mode : modes) {
                        const bool selected = mode == ao_mode;
                        if (ImGui::Selectable(trexio_ao_ordering_mode_str(mode), selected)) {
                            settings.ao_ordering_mode = mode;
                            on_ao_ordering_mode_changed(state);
                        }
                        if (selected) ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }
                ImGui::SetItemTooltip("Manual debug override for unusual TREXIO-like sources. Leave this on TREXIO spherical for standard files.");
                ImGui::TreePop();
            }
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("System Information")) {
            ImGui::Text("Num Nuclei:           %-6zu", data.atoms.size());
            ImGui::Text("Num Electrons:        %-6d", data.electron_num);
            ImGui::Text("Num Alpha Electrons:  %-6d", data.electron_up_num);
            ImGui::Text("Num Beta Electrons:   %-6d", data.electron_dn_num);
            ImGui::Text("Num Shells:           %-6zu", data.shells.size());
            ImGui::Text("Num AOs:              %-6d", data.ao_num);
            ImGui::Text("Num MOs:              %-6d", data.mo_num);
            ImGui::TreePop();
        }

        if (!data.orbitals.empty() && ImGui::TreeNode("Orbital Information")) {
            const int alpha_homo = homo_idx(ElectronicStructureSpin::Alpha);
            const int alpha_lumo = lumo_idx(ElectronicStructureSpin::Alpha);
            ImGui::Text("Alpha HOMO:           %d", alpha_homo + 1);
            ImGui::Text("Alpha LUMO:           %d", alpha_lumo + 1);
            if (has_beta_spin()) {
                const int beta_homo = homo_idx(ElectronicStructureSpin::Beta);
                const int beta_lumo = lumo_idx(ElectronicStructureSpin::Beta);
                ImGui::Text("Beta HOMO:            %d", beta_homo + 1);
                ImGui::Text("Beta LUMO:            %d", beta_lumo + 1);
            }
            ImGui::TreePop();
        }

        if (!data.atoms.empty() && ImGui::TreeNode("Geometry")) {
            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;
            if (ImGui::BeginTable("TREXIO Geometry Table", 6, flags, ImVec2(520, -1), 0)) {
                ImGui::TableSetupColumn("Atom");
                ImGui::TableSetupColumn("Symbol");
                ImGui::TableSetupColumn("Z");
                ImGui::TableSetupColumn("X (Bohr)");
                ImGui::TableSetupColumn("Y (Bohr)");
                ImGui::TableSetupColumn("Z (Bohr)");
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                for (size_t row_n = 0; row_n < data.atoms.size(); ++row_n) {
                    const trexio_data::Atom& atom = data.atoms[row_n];
                    const bool is_sel = md_bitfield_test_bit(&state.selection.selection_mask, row_n);
                    const bool is_hov = md_bitfield_test_bit(&state.selection.highlight_mask, row_n);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char label[16];
                    snprintf(label, sizeof(label), "%zu", row_n + 1);
                    ImGui::Selectable(label, is_sel || is_hov, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap);
                    ImGui::TableNextColumn();
                    str_t symbol = md_atomic_number_symbol((md_atomic_number_t)CLAMP(atom.atomic_number, 0, 118));
                    ImGui::Text("%.*s", (int)symbol.len, symbol.ptr ? symbol.ptr : "");
                    ImGui::TableNextColumn();
                    ImGui::Text("%d", atom.atomic_number);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.x);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.y);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.z);
                }
                ImGui::EndTable();
            }
            ImGui::TreePop();
        }

        ImGui::End();
    }

    void draw_orb_window(const ApplicationState& state) {
        if (!orb.show_window || !loaded || !basis_ready) return;

        ensure_visual_state(state);
        if (!orb.gl_rep.id) return;

        const bool unrestricted = has_beta_spin();
        const size_t alpha_count = num_orbitals(ElectronicStructureSpin::Alpha);
        const size_t beta_count = num_orbitals(ElectronicStructureSpin::Beta);
        const int num_list_orbitals = (int)(unrestricted ? MAX(alpha_count, beta_count) : alpha_count);
        if (num_list_orbitals <= 0) return;

        ImGui::SetNextWindowSize({600, 300}, ImGuiCond_FirstUseEver);
        if (!qm_ui::begin_orbital_grid_window("TREXIO", "TrexioQMOrbitalGrid", &orb.show_window)) {
            ImGui::End();
            return;
        }

        const double* occ_alpha = orbital_occupancy(ElectronicStructureSpin::Alpha);
        const double* occ_beta = orbital_occupancy(ElectronicStructureSpin::Beta);
        const double* ene_alpha = orbital_energy(ElectronicStructureSpin::Alpha);
        const double* ene_beta = orbital_energy(ElectronicStructureSpin::Beta);

        const float text_base_height = ImGui::GetTextLineHeightWithSpacing();
        int num_x = unrestricted ? 2 : orb.num_x;
        int num_y = orb.num_y;
        int num_mos = num_x * num_y;
        int beg_mo_idx = orb.mo_idx - num_mos / 2 + (num_mos % 2 == 0 ? 1 : 0);
        int window_size = num_mos;
        if (unrestricted) {
            beg_mo_idx = orb.mo_idx - num_y / 2 + (num_y % 2 == 0 ? 1 : 0);
            window_size = num_y;
        }

        {
            ImGui::BeginChild("left pane", ImVec2(300, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
            ImGui::SameLine();
            ImGui::PushItemWidth(-1);
            ImGui::BeginGroup();

            ImGui::SliderInt("##Rows", &orb.num_y, 1, 4);
            if (unrestricted) ImGui::PushDisabled();
            ImGui::SliderInt("##Cols", &orb.num_x, 1, 4);
            if (unrestricted) ImGui::PopDisabled();

            const double iso_min = 1.0e-4;
            const double iso_max = 5.0;
            double iso_val = orb.iso.values[0];
            ImGui::TextUnformatted(electronic_structure_iso_value_label());
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            orb.iso.values[0] = (float)iso_val;
            orb.iso.values[1] = -(float)iso_val;
            orb.iso.count = 2;
            ImGui::ColorEdit4("##Color Positive", orb.iso.colors[0].elem);
            ImGui::SetItemTooltip("Color Positive");
            ImGui::ColorEdit4("##Color Negative", orb.iso.colors[1].elem);
            ImGui::SetItemTooltip("Color Negative");

            if (ImGui::IsWindowAppearing()) orb.scroll_to_idx = orb.mo_idx;
            if (ImGui::Button("Goto HOMO", ImVec2(-1, 0))) {
                orb.scroll_to_idx = homo_idx(ElectronicStructureSpin::Alpha);
                orb.mo_idx = homo_idx(ElectronicStructureSpin::Alpha);
            }

            const ImGuiTableFlags flags = ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY;
            const int num_cols = unrestricted ? 5 : 3;
            if (ImGui::BeginTable("TREXIO Molecular Orbitals", num_cols, flags)) {
                ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed);
                ImGui::TableSetupColumn(unrestricted ? (const char*)u8"Occ. α" : "Occupancy", ImGuiTableColumnFlags_WidthFixed);
                if (unrestricted) ImGui::TableSetupColumn((const char*)u8"Occ. β", ImGuiTableColumnFlags_WidthFixed);
                ImGui::TableSetupColumn(unrestricted ? (const char*)u8"Ene. α" : "Energy", ImGuiTableColumnFlags_WidthFixed);
                if (unrestricted) ImGui::TableSetupColumn((const char*)u8"Ene. β", ImGuiTableColumnFlags_WidthFixed);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                for (int n = num_list_orbitals - 1; n >= 0; --n) {
                    ImGui::PushID(n + 1);
                    ImGui::TableNextRow();
                    const bool is_selected = (beg_mo_idx <= n && n < beg_mo_idx + window_size);
                    ImGui::TableNextColumn();
                    if (orb.scroll_to_idx != -1 && n == orb.scroll_to_idx) {
                        orb.scroll_to_idx = -1;
                        ImGui::SetScrollHereY();
                    }
                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i", n + 1);
                    if (ImGui::Selectable(buf, is_selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) {
                        orb.mo_idx = n;
                    }

                    ImGui::TableNextColumn();
                    if (occ_alpha && n < (int)alpha_count) ImGui::Text("%.1f %s", occ_alpha[n], n == homo_idx(ElectronicStructureSpin::Alpha) ? "HOMO" : (n == lumo_idx(ElectronicStructureSpin::Alpha) ? "LUMO" : ""));
                    else ImGui::TextUnformatted("-");
                    if (unrestricted) {
                        ImGui::TableNextColumn();
                        if (occ_beta && n < (int)beta_count) ImGui::Text("%.1f %s", occ_beta[n], n == homo_idx(ElectronicStructureSpin::Beta) ? "HOMO" : (n == lumo_idx(ElectronicStructureSpin::Beta) ? "LUMO" : ""));
                        else ImGui::TextUnformatted("-");
                    }
                    ImGui::TableNextColumn();
                    if (ene_alpha && n < (int)alpha_count) ImGui::Text("%.4f", ene_alpha[n]); else ImGui::TextUnformatted("-");
                    if (unrestricted) {
                        ImGui::TableNextColumn();
                        if (ene_beta && n < (int)beta_count) ImGui::Text("%.4f", ene_beta[n]); else ImGui::TextUnformatted("-");
                    }
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }

            ImGui::EndGroup();
            ImGui::PopItemWidth();
            ImGui::EndChild();
        }

        ImGui::SameLine();

        int vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        ElectronicStructureSpin vol_mo_spin[16] = {};
        for (int i = 0; i < num_mos; ++i) {
            if (unrestricted) {
                vol_mo_idx[i] = beg_mo_idx + i / 2;
                vol_mo_spin[i] = (i & 1) ? ElectronicStructureSpin::Alpha : ElectronicStructureSpin::Beta;
            } else {
                vol_mo_idx[i] = beg_mo_idx + i;
                vol_mo_spin[i] = ElectronicStructureSpin::Alpha;
            }
        }

        int job_queue[16];
        int num_jobs = 0;
        for (int i = 0; i < num_mos; ++i) {
            if (orb.vol_mo_idx[i] == vol_mo_idx[i] && orb.vol_mo_spin[i] == vol_mo_spin[i]) continue;
            bool found = false;
            for (int j = 0; j < num_mos; ++j) {
                if (i == j) continue;
                if (vol_mo_idx[i] == orb.vol_mo_idx[j] && vol_mo_spin[i] == orb.vol_mo_spin[j]) {
                    ImSwap(orb.vol[i], orb.vol[j]);
                    ImSwap(orb.vol_mo_idx[i], orb.vol_mo_idx[j]);
                    ImSwap(orb.vol_mo_spin[i], orb.vol_mo_spin[j]);
                    found = true;
                    break;
                }
            }
            if (!found) job_queue[num_jobs++] = i;
        }

        if (num_jobs > 0) {
            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };
            vec4_t* atom_xyzw = build_atom_positions_bohr(&state.mold.sys, md_temp_allocator(temp));
            TrexioAABB aabb = {};
            trexio_calculate_bounds(&aabb, atom_xyzw, md_array_size(atom_xyzw));
            md_grid_t grid = {};
            trexio_init_grid(&grid, aabb.min_ext, aabb.max_ext, trexio_orbital_grid_samples_per_angstrom * TREXIO_BOHR_TO_ANGSTROM);

            for (int i = 0; i < num_jobs; ++i) {
                const int slot_idx = job_queue[i];
                const int mo_idx = vol_mo_idx[slot_idx];
                const ElectronicStructureSpin spin = vol_mo_spin[slot_idx];
                orb.vol_mo_idx[slot_idx] = mo_idx;
                orb.vol_mo_spin[slot_idx] = spin;
                if (orbital_from_local_idx(spin, mo_idx)) {
                    trexio_init_volume(&orb.vol[slot_idx], grid);
                    evaluate_orbital_texture(orb.vol[slot_idx].tex_id, grid, atom_xyzw, spin, mo_idx);
                }
            }
        }

        ImVec2 canvas_sz = ImMax(ImVec2(50.0f, 50.0f), ImGui::GetContentRegionAvail());
        camera_animate(&orb.camera, orb.target, state.app.timing.delta_s);
        ImGui::Dummy(canvas_sz);

        ImGuiIO& io = ImGui::GetIO();
        ImVec2 canvas_min = ImGui::GetItemRectMin();
        ImVec2 canvas_max = ImGui::GetItemRectMax();
        ImVec2 orb_win_sz = ImFloor(canvas_sz / ImVec2((float)num_x, (float)num_y));
        const float aspect = orb_win_sz.x / orb_win_sz.y;

        const mat4_t view_to_world = camera_view_to_world_matrix(orb.camera);
        const mat4_t clip_to_view = camera_clip_to_view_matrix_persp(orb.camera, aspect);
        (void)view_to_world;
        (void)clip_to_view;

        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddRectFilled(canvas_min, canvas_max, IM_COL32(255, 255, 255, 255));

        for (int i = 0; i < num_mos; ++i) {
            ImGui::PushID(i);
            defer { ImGui::PopID(); };
            const int mo_idx = unrestricted ? beg_mo_idx + i / 2 : beg_mo_idx + i;
            const ElectronicStructureSpin spin = unrestricted ? ((i & 1) ? ElectronicStructureSpin::Alpha : ElectronicStructureSpin::Beta) : ElectronicStructureSpin::Alpha;
            const bool valid = orbital_from_local_idx(spin, mo_idx) != nullptr;
            const int x = num_x - i % num_x - 1;
            const int y = num_y - i / num_x - 1;

            if (valid) {
                const ImVec2 p0 = canvas_min + orb_win_sz * ImVec2((float)(x + 0), (float)(y + 0));
                const ImVec2 p1 = canvas_min + orb_win_sz * ImVec2((float)(x + 1), (float)(y + 1));
                const ImVec2 text_pos_bl = ImVec2(p0.x + text_base_height * 0.5f, p1.y - text_base_height);
                const ImVec2 text_pos_tl = ImVec2(p0.x + text_base_height * 0.5f, p0.y + text_base_height * 0.25f);
                const ImVec2 text_pos_br = ImVec2(p1.x - text_base_height * 0.5f, p1.y - text_base_height);
                const ImVec2 sz = p1 - p0;

                ImGui::SetCursorScreenPos(p0);
                InteractionSurfaceState surface_state = interaction_surface(interaction_surface_trexio_orb, vec_cast(sz), InteractionSurfaceFlags_NoRegionSelect);
                InteractionSurfaceViewTransformArgs view_args = {.camera = orb.camera, .trackball_param = state.view.trackball_param};
                InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&orb.target, surface_state, view_args);
                if (view_result.reset_requested) orb.target = orb.default_view;

                const char* lbl = "";
                if (mo_idx == homo_idx(spin)) lbl = "(HOMO)";
                else if (mo_idx == lumo_idx(spin)) lbl = "(LUMO)";

                char buf[32];
                snprintf(buf, sizeof(buf), "%i %s", mo_idx + 1, lbl);
                draw_list->AddImage((ImTextureID)(intptr_t)orb.gbuf.tex.transparency, p0, p1, {0, 1}, {1, 0});
                draw_list->AddImage((ImTextureID)(intptr_t)orb.iso_tex[i], p0, p1, {0, 1}, {1, 0});
                draw_list->AddText(text_pos_bl, ImColor(0, 0, 0), buf);

                if (unrestricted) {
                    draw_list->AddText(text_pos_tl, ImColor(0, 0, 0), spin == ElectronicStructureSpin::Alpha ? (const char*)u8"α" : (const char*)u8"β");
                    snprintf(buf, sizeof(buf), "%.4f", spin == ElectronicStructureSpin::Alpha ? ene_alpha[mo_idx] : ene_beta[mo_idx]);
                } else {
                    snprintf(buf, sizeof(buf), "%.4f", ene_alpha[mo_idx]);
                }
                float width = ImGui::CalcTextSize(buf).x;
                draw_list->AddText(text_pos_br - ImVec2(width, 0), ImColor(0, 0, 0), buf);
            }
        }

        for (int x = 1; x < num_x; ++x) draw_list->AddLine({canvas_min.x + orb_win_sz.x * x, canvas_min.y}, {canvas_min.x + orb_win_sz.x * x, canvas_max.y}, IM_COL32(0, 0, 0, 255));
        for (int y = 1; y < num_y; ++y) draw_list->AddLine({canvas_min.x, canvas_min.y + orb_win_sz.y * y}, {canvas_max.x, canvas_min.y + orb_win_sz.y * y}, IM_COL32(0, 0, 0, 255));

        int width = MAX(1, (int)(orb_win_sz.x * io.DisplayFramebufferScale.x));
        int height = MAX(1, (int)(orb_win_sz.y * io.DisplayFramebufferScale.y));
        if ((int)orb.gbuf.width != width || (int)orb.gbuf.height != height) {
            init_gbuffer(&orb.gbuf, width, height);
            for (int i = 0; i < num_mos; ++i) gl::init_texture_2D(orb.iso_tex + i, width, height, GL_RGBA8);
        }

        if (orb.show_coordinate_system_widget) {
            float ext = MIN(orb_win_sz.x, orb_win_sz.y) * 0.4f;
            float pad = 20.0f;
            ImGui::SetCursorScreenPos(ImVec2(canvas_min.x + pad, canvas_max.y - ext - pad));
            quat_t out_orientation = orb.target.orientation;
            if (ImGui::CoordinateSystemWidget(&out_orientation, orb.camera.orientation, ImVec2(ext, ext))) {
                const vec3_t look_at = camera_get_look_at(orb.target);
                orb.target.orientation = quat_normalize(out_orientation);
                orb.target.position = camera_position_from_look_at(look_at, orb.target.orientation, orb.target.distance);
            }
        }

        const mat4_t view_mat = camera_world_to_view_matrix(orb.camera);
        const mat4_t proj_mat = camera_view_to_clip_matrix_persp(orb.camera, aspect);
        const mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(orb.camera, aspect);

        clear_gbuffer(&orb.gbuf);
        const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY, GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY};

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDepthMask(GL_TRUE);
        glEnable(GL_SCISSOR_TEST);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, orb.gbuf.fbo);
        glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
        glViewport(0, 0, orb.gbuf.width, orb.gbuf.height);
        glScissor(0, 0, orb.gbuf.width, orb.gbuf.height);

        md_gl_draw_op_t draw_op = {};
        draw_op.type = MD_GL_REP_BALL_AND_STICK;
        draw_op.args.ball_and_stick.ball_scale = 1.0f;
        draw_op.args.ball_and_stick.stick_radius = 1.0f;
        draw_op.rep = orb.gl_rep;
        md_gl_draw_args_t draw_args = {
            .shaders = state.gl.shaders,
            .draw_operations = {.count = 1, .ops = &draw_op},
            .view_transform = {.view_matrix = (const float*)view_mat.elem, .proj_matrix = (const float*)proj_mat.elem},
            .picking_offset = {.atom_base = state.picking_range_atom.beg, .bond_base = state.picking_range_bond.beg},
        };
        md_gl_draw(&draw_args);

        glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
        glClearColor(1, 1, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT);

        postprocessing::Descriptor postprocess_desc = {
            .background = {.color = {24.f, 24.f, 24.f}},
            .tonemapping = {.enabled = state.visuals.tonemapping.enabled, .mode = state.visuals.tonemapping.tonemapper, .exposure = state.visuals.tonemapping.exposure, .gamma = state.visuals.tonemapping.gamma},
            .ambient_occlusion = {.enabled = false},
            .depth_of_field = {.enabled = false},
            .fxaa = {.enabled = true},
            .temporal_aa = {.enabled = false},
            .sharpen = {.enabled = false},
            .input_textures = {.depth = orb.gbuf.tex.depth, .color = orb.gbuf.tex.color, .normal = orb.gbuf.tex.normal, .velocity = orb.gbuf.tex.velocity},
        };
        ViewParam view_param = {
            .matrix = {.curr = {.view = view_mat, .proj = proj_mat, .norm = view_mat}, .inv = {.proj = inv_proj_mat}},
            .clip_planes = {.near = orb.camera.near_plane, .far = orb.camera.far_plane},
            .resolution = {orb_win_sz.x, orb_win_sz.y},
            .fov_y = orb.camera.fov_y,
        };
        postprocessing::shade_and_postprocess(postprocess_desc, view_param);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);
        glDisable(GL_SCISSOR_TEST);

        for (int i = 0; i < num_mos; ++i) {
            if (!orb.vol[i].tex_id || !orb.iso_tex[i]) continue;
            volume::RenderDesc vol_desc = {
                .render_target = {.depth = orb.gbuf.tex.depth, .color = orb.iso_tex[i], .width = orb.gbuf.width, .height = orb.gbuf.height, .clear_color = true},
                .texture = {.density_volume = orb.vol[i].tex_id},
                .matrix = {.model = orb.vol[i].texture_to_world, .view = view_mat, .proj = proj_mat, .inv_proj = inv_proj_mat},
                .iso = {.enabled = true, .count = orb.iso.count, .values = orb.iso.values, .colors = orb.iso.colors},
                .shading = {.env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f, .roughness = 0.3f, .dir_radiance = {10, 10, 10}, .ior = 1.5f},
                .voxel_spacing = orb.vol[i].voxel_size,
            };
            volume::render_volume(vol_desc);
        }

        ImGui::End();
    }

    bool evaluate_molecular_orbital(EvalElectronicStructure& eval) const {
        Representation* rep = eval.rep;
        const ElectronicStructureRepresentation& es = rep->electronic_structure;
        const trexio_data::Orbital* orbital = selected_orbital(es);
        if (!orbital) {
            MD_LOG_ERROR("Invalid TREXIO orbital selection: spin=%d idx=%d", (int)es.spin, es.orbital_idx);
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_allocator_i* temp_alloc = md_temp_allocator(temp);
        vec4_t* atom_xyzw = build_atom_positions_bohr(eval.sys, temp_alloc);
        TrexioAABB aabb = {};
        trexio_calculate_bounds(&aabb, atom_xyzw, md_array_size(atom_xyzw));

        const double samples_per_angstrom = trexio_volume_resolution_samples_per_angstrom[(int)es.resolution];
        const double samples_per_unit_length = samples_per_angstrom * TREXIO_BOHR_TO_ANGSTROM;
        md_grid_t grid = {};
        trexio_init_grid(&grid, aabb.min_ext, aabb.max_ext, samples_per_unit_length);
        trexio_init_volume(&rep->electronic_structure.density_vol, grid, GL_R32F);

        double* ao_coeffs = md_temp_alloc_array(temp, double, md_array_size(ordered_to_source_ao));
        reorder_coefficients(ao_coeffs, *orbital);
        md_gto_grid_evaluate_mo_GL(rep->electronic_structure.density_vol.tex_id, &grid, &basis, (const float*)atom_xyzw, sizeof(vec4_t), ao_coeffs, TREXIO_DEFAULT_GTO_CUTOFF_VALUE, MD_GTO_EVAL_MODE_PSI, trexio_gto_op_from_use_magnitude(es.use_magnitude));
        return true;
    }

    void process_events(const viamd::Event* events, size_t num_events) override {
        for (size_t i = 0; i < num_events; ++i) {
            const auto& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdFrameTick: {
                if (!loaded) break;
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw_summary_window(state);
                draw_orb_window(state);
            } break;

            case viamd::EventType_ViamdWindowDrawMenu:
                if (loaded && ImGui::BeginMenu("TREXIO")) {
                    ImGui::Checkbox(qm_ui::summary_label, &summary.show_window);
                    ImGui::Checkbox(qm_ui::orbital_grid_label, &orb.show_window);
                    ImGui::EndMenu();
                }
                break;

            case viamd::EventType_ViamdLoadData: {
                auto& payload = *(LoadDataPayload*)e.payload;
                if (payload.loader_state.type == LoaderType_TREXIO) {
                    load_file(payload.path_to_file);
                }
            } break;

            case viamd::EventType_ViamdSystemFree:
                reset();
                break;

            case viamd::EventType_ViamdRepresentationInfoFill: {
                if (!loaded) break;
                RepresentationInfo& info = *(RepresentationInfo*)e.payload;
                if (!info.alloc) break;
                fill_orbitals(&info.alpha, *this, ElectronicStructureSpin::Alpha, info.alloc);
                fill_orbitals(&info.beta, *this, ElectronicStructureSpin::Beta, info.alloc);
                if (info.alpha.num_orbitals > 0 || info.beta.num_orbitals > 0) {
                    info.electronic_structure_source_mask |= ElectronicStructureSourceFlag_MolecularOrbital;
                }
            } break;

            case viamd::EventType_ViamdRepresentationEvalElectronicStructure: {
                if (!loaded || !basis_ready) break;
                EvalElectronicStructure& eval = *(EvalElectronicStructure*)e.payload;
                if (!eval.rep || !eval.sys) break;
                Representation* rep = eval.rep;
                const ElectronicStructureRepresentation& es = rep->electronic_structure;
                if (es.source != ElectronicStructureSource::MolecularOrbital) break;

                uint64_t frame_hash = md_hash64(&eval.frame, sizeof(eval.frame), 0);
                uint64_t vol_hash = 0;
                vol_hash = md_hash64(&es.source, sizeof(es.source), vol_hash);
                vol_hash = md_hash64(&es.use_magnitude, sizeof(es.use_magnitude), vol_hash);
                vol_hash = md_hash64(&es.spin, sizeof(es.spin), vol_hash);
                vol_hash = md_hash64(&es.resolution, sizeof(es.resolution), vol_hash);
                vol_hash = md_hash64(&es.orbital_idx, sizeof(es.orbital_idx), vol_hash);
                vol_hash = md_hash64_combine(vol_hash, frame_hash);
                if (vol_hash != rep->electronic_structure.vol_hash && evaluate_molecular_orbital(eval)) {
                    rep->electronic_structure.vol_hash = vol_hash;
                }
            } break;
            }
        }
    }
};

static TrexioComponent instance = {};

#endif // VIAMD_TREXIO
