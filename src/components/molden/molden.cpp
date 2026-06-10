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
#include <molden.h>
#include <qm_ui.h>

#include <algorithm>
#include <app/IconsFontAwesome6.h>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <string>

#include <imgui_internal.h>

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903
#define DEFAULT_GTO_CUTOFF_VALUE 1.0e-6

static const double PI_D = 3.14159265358979323846;

static const double volume_resolution_samples_per_angstrom[3] = {
    4.0,
    8.0,
    16.0,
};

static const double orbital_grid_samples_per_angstrom = 8.0;
constexpr uint64_t interaction_surface_molden_orb = HASH_STR_LIT64("interaction surface molden orb");

struct AABB {
    vec3_t min_ext = {0};
    vec3_t max_ext = {0};
};

struct MoldenShellEntry {
    const molden::ContractedShell* shell = nullptr;
    uint32_t atom_idx = 0;
    uint32_t l = 0;
    uint32_t source_ao_offset = 0;
    uint32_t source_order = 0;
    bool use_sp_coeff = false;
};

enum class MoldenAoOrderingMode : uint8_t {
    Auto,
    Identity,
    PStandard,
    PSpSwapped,
};

static inline md_gto_op_t gto_op_from_use_magnitude(bool use_magnitude) {
    return MD_GTO_OP_SET | (use_magnitude ? MD_GTO_OP_ABS : 0);
}

static inline mat4_t compute_texture_to_world_mat(const mat3_t& orientation, const vec3_t& origin, const vec3_t& extent) {
    mat4_t T = mat4_translate_vec3(origin);
    mat4_t R = mat4_from_mat3(orientation);
    mat4_t S = mat4_scale_vec3(extent);
    return T * R * S;
}

static inline mat4_t compute_world_to_model_mat(const mat3_t& orientation, const vec3_t& origin) {
    return mat4_from_mat3(mat3_transpose(orientation)) * mat4_translate_vec3(-origin);
}

static inline void compute_dim(int out_dim[3], const vec3_t& extent, double samples_per_unit_length) {
    out_dim[0] = CLAMP(ALIGN_TO((int)(extent.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = CLAMP(ALIGN_TO((int)(extent.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = CLAMP(ALIGN_TO((int)(extent.z * samples_per_unit_length), 8), 8, 512);
}

static void init_grid(md_grid_t* grid, const vec3_t& min_ext, const vec3_t& max_ext, double samples_per_unit_length) {
    vec3_t extent = max_ext - min_ext;
    compute_dim(grid->dim, extent, samples_per_unit_length);
    grid->orientation = mat3_ident();
    grid->origin = min_ext;
    grid->spacing = vec3_div(extent, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
}

static void init_volume(Volume* vol, const md_grid_t& grid, GLenum format = GL_R32F) {
    ASSERT(vol);
    MEMCPY(vol->dim, grid.dim, sizeof(vol->dim));
    const float scl = BOHR_TO_ANGSTROM;
    vec3_t extent = md_grid_extent(&grid);
    vol->world_to_model = compute_world_to_model_mat(grid.orientation, grid.origin * scl);
    vol->texture_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
    vol->voxel_size = grid.spacing * scl;
    gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
}

static void calculate_bounds(AABB* aabb, const vec4_t* atom_xyzw, size_t count) {
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

static bool orbital_matches_spin(const molden::MolecularOrbital& orbital, molden::SpinType spin) {
    if (spin == molden::SpinType::Beta) {
        return orbital.spin == molden::SpinType::Beta;
    }
    return orbital.spin != molden::SpinType::Beta;
}

static int shell_l_value(molden::ShellType shell_type) {
    switch (shell_type) {
    case molden::ShellType::S: return 0;
    case molden::ShellType::P: return 1;
    case molden::ShellType::D: return 2;
    case molden::ShellType::F: return 3;
    case molden::ShellType::G: return 4;
    default: return -1;
    }
}

static size_t shell_num_spherical_functions(int l) {
    return (size_t)(2 * l + 1);
}

static const char* ao_ordering_mode_str(MoldenAoOrderingMode mode) {
    switch (mode) {
    case MoldenAoOrderingMode::Auto: return "Auto";
    case MoldenAoOrderingMode::Identity: return "Source Order";
    case MoldenAoOrderingMode::PStandard: return "Standard P Order";
    case MoldenAoOrderingMode::PSpSwapped: return "SP-Derived P Order";
    default: return "Unknown";
    }
}

static bool contains_case_insensitive(const std::string& text, const char* needle) {
    if (!needle || !needle[0]) return true;
    const size_t needle_len = strlen(needle);
    if (needle_len > text.size()) return false;

    for (size_t i = 0; i + needle_len <= text.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < needle_len; ++j) {
            const unsigned char lhs = (unsigned char)text[i + j];
            const unsigned char rhs = (unsigned char)needle[j];
            if (std::tolower(lhs) != std::tolower(rhs)) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}

static size_t shell_component_source_index(int l, MoldenAoOrderingMode mode, size_t component_idx) {
    // md_gto expands p shells in the internal order [py, pz, px]. Molden files
    // in the wild are not consistent, so the remap is exposed explicitly.
    if (l == 1) {
        static const size_t p_order[] = {1, 2, 0};
        static const size_t sp_order[] = {0, 2, 1};
        static const size_t identity_order[] = {0, 1, 2};
        ASSERT(component_idx < ARRAY_SIZE(p_order));
        switch (mode) {
        case MoldenAoOrderingMode::Identity:
            return identity_order[component_idx];
        case MoldenAoOrderingMode::PStandard:
            return p_order[component_idx];
        case MoldenAoOrderingMode::PSpSwapped:
            return sp_order[component_idx];
        case MoldenAoOrderingMode::Auto:
        default:
            return identity_order[component_idx];
        }
    }
    return component_idx;
}

static int find_atom_vector_index(const molden::MoldenData& data, int32_t atom_index) {
    for (size_t i = 0; i < data.atoms.size(); ++i) {
        if (data.atoms[i].atom_index == atom_index) {
            return (int)i;
        }
    }
    return -1;
}

static double primitive_overlap(int l, const float* exponents, const float* coeffs, int i, int j) {
    const double fab = 1.0 / (exponents[i] + exponents[j]);
    const double fab2 = fab * fab;
    const double ovl = coeffs[i] * coeffs[j] * pow(PI_D * fab, 1.5);

    switch (l) {
    case 0: return ovl;
    case 1: return 0.5 * fab * ovl;
    case 2: return 3.0 * fab2 * ovl;
    case 3: return 7.5 * fab2 * fab * ovl;
    case 4: return 420.0 * fab2 * fab2 * ovl;
    default:
        ASSERT(false);
        return 0.0;
    }
}

static void normalize_shell_coefficients(int l, float* exponents, float* coeffs, int count) {
    if (count <= 0) return;

    if (count == 1) {
        coeffs[0] = 1.0f;
    }

    const double fpi = 2.0 / PI_D;
    static const double f_table[] = {
        0.0,
        2.0,
        1.1547005383792515,
        1.0327955589886445,
        0.19518001458970664,
    };

    for (int i = 0; i < count; ++i) {
        coeffs[i] = (float)(coeffs[i] * pow(exponents[i] * fpi, 0.75));
        if (l > 0) {
            coeffs[i] = (float)(coeffs[i] * pow(f_table[l] * exponents[i], l * 0.5));
        }
    }

    double ovl = 0.0;
    for (int i = 0; i < count; ++i) {
        ovl += primitive_overlap(l, exponents, coeffs, i, i);
        for (int j = i + 1; j < count; ++j) {
            ovl += 2.0 * primitive_overlap(l, exponents, coeffs, i, j);
        }
    }

    const double scale = ovl > 0.0 ? 1.0 / sqrt(ovl) : 1.0;
    for (int i = 0; i < count; ++i) {
        coeffs[i] = (float)(coeffs[i] * scale);
    }
}

struct MoldenComponent : viamd::EventHandler {
    struct Settings {
        MoldenAoOrderingMode ao_ordering_mode = MoldenAoOrderingMode::Auto;
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
    molden::MoldenData data = {};
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

    MoldenComponent() {
        viamd::event_system_register_handler(*this);
    }

    MoldenAoOrderingMode resolved_ao_ordering_mode() const {
        if (settings.ao_ordering_mode != MoldenAoOrderingMode::Auto) {
            return settings.ao_ordering_mode;
        }

        if (contains_case_insensitive(data.title, "gansu")) {
            return MoldenAoOrderingMode::PStandard;
        }

        bool has_sp_shell = false;
        bool has_standalone_p_shell = false;
        for (const auto& atom_basis : data.basis_sets) {
            for (const auto& shell : atom_basis.shells) {
                has_sp_shell |= shell.shell_type == molden::ShellType::SP;
                has_standalone_p_shell |= shell.shell_type == molden::ShellType::P;
            }
        }

        if (has_sp_shell && !has_standalone_p_shell) {
            return MoldenAoOrderingMode::PSpSwapped;
        }

        return MoldenAoOrderingMode::PStandard;
    }

    void invalidate_orbital_caches() {
        for (int i = 0; i < (int)ARRAY_SIZE(orb.vol_mo_idx); ++i) {
            orb.vol_mo_idx[i] = -1;
            orb.vol_mo_spin[i] = {};
        }
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

    static bool is_molden_path(str_t path) {
        str_t ext = {};
        return extract_ext(&ext, path) &&
               (str_eq_ignore_case(ext, STR_LIT("molden")) ||
                str_eq_ignore_case(ext, STR_LIT("mold")));
    }

    static size_t find_homo_idx(const molden::MoldenData& molden_data, molden::SpinType spin) {
        size_t homo_idx = 0;
        bool found = false;
        size_t local_idx = 0;
        for (const auto& orbital : molden_data.orbitals) {
            if (!orbital_matches_spin(orbital, spin)) continue;
            if (orbital.occupation > 0.0) {
                homo_idx = local_idx;
                found = true;
            }
            ++local_idx;
        }
        return found ? homo_idx : 0;
    }

    static const char* coordinate_unit_str(molden::CoordinateUnit unit) {
        switch (unit) {
        case molden::CoordinateUnit::Angstrom: return "Angstrom";
        case molden::CoordinateUnit::AtomicUnit: return "Bohr";
        default: return "Unknown";
        }
    }

    static const char* basis_format_str(molden::BasisFormat format) {
        switch (format) {
        case molden::BasisFormat::Cartesian: return "Cartesian";
        case molden::BasisFormat::Spherical: return "Spherical";
        default: return "Unknown";
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

    static void fill_orbitals(::MolecularOrbital* dst, const molden::MoldenData& molden_data, molden::SpinType spin, md_allocator_i* alloc) {
        ASSERT(dst);
        ASSERT(alloc);

        const size_t num_orbitals = spin == molden::SpinType::Beta ? molden_data.num_beta_orbitals : molden_data.num_alpha_orbitals;
        if (num_orbitals == 0) {
            *dst = {};
            return;
        }

        dst->num_orbitals = num_orbitals;
        dst->homo_idx = find_homo_idx(molden_data, spin);
        dst->lumo_idx = MIN(dst->homo_idx + 1, num_orbitals - 1);

        md_array_resize(dst->label, num_orbitals, alloc);
        md_array_resize(dst->occupation, num_orbitals, alloc);
        md_array_resize(dst->energy, num_orbitals, alloc);

        size_t local_idx = 0;
        for (const auto& orbital : molden_data.orbitals) {
            if (!orbital_matches_spin(orbital, spin)) continue;

            if (!orbital.symmetry.empty()) {
                dst->label[local_idx] = str_printf(alloc, "%zu %s", local_idx + 1, orbital.symmetry.c_str());
            } else {
                dst->label[local_idx] = str_printf(alloc, "%zu", local_idx + 1);
            }
            dst->occupation[local_idx] = orbital.occupation;
            dst->energy[local_idx] = orbital.energy;
            ++local_idx;
        }
    }

    void build_orbital_maps() {
        alpha_orbital_map = nullptr;
        beta_orbital_map = nullptr;
        alpha_occupation_cache = nullptr;
        beta_occupation_cache = nullptr;
        alpha_energy_cache = nullptr;
        beta_energy_cache = nullptr;

        for (size_t i = 0; i < data.orbitals.size(); ++i) {
            if (orbital_matches_spin(data.orbitals[i], molden::SpinType::Beta)) {
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

        if (data.basis_sets.empty() || data.orbitals.empty()) {
            return false;
        }

        md_array(MoldenShellEntry) shell_entries = nullptr;
        uint32_t source_ao_offset = 0;
        uint32_t source_order = 0;

        for (const auto& atom_basis : data.basis_sets) {
            const int atom_idx = find_atom_vector_index(data, atom_basis.atom_index);
            if (atom_idx < 0) {
                MD_LOG_ERROR("Molden basis references unknown atom index %d", atom_basis.atom_index);
                return false;
            }

            for (const auto& shell : atom_basis.shells) {
                const size_t source_func_count = molden::util::get_num_basis_functions(shell.shell_type, data.basis_format);
                if (shell.shell_type == molden::ShellType::SP) {
                    MoldenShellEntry s_entry = { &shell, (uint32_t)atom_idx, 0u, source_ao_offset, source_order++, false };
                    MoldenShellEntry p_entry = { &shell, (uint32_t)atom_idx, 1u, source_ao_offset + 1, source_order++, true };
                    md_array_push(shell_entries, s_entry, arena);
                    md_array_push(shell_entries, p_entry, arena);
                } else {
                    const int l = shell_l_value(shell.shell_type);
                    if (l < 0) {
                        MD_LOG_ERROR("Unsupported Molden shell type encountered");
                        return false;
                    }
                    if (l >= 2 && data.basis_format == molden::BasisFormat::Cartesian) {
                        MD_LOG_ERROR("Cartesian Molden shells with l >= 2 are not supported yet");
                        return false;
                    }

                    MoldenShellEntry entry = { &shell, (uint32_t)atom_idx, (uint32_t)l, source_ao_offset, source_order++, false };
                    md_array_push(shell_entries, entry, arena);
                }

                source_ao_offset += (uint32_t)source_func_count;
            }
        }

        std::sort(shell_entries, shell_entries + md_array_size(shell_entries), [](const MoldenShellEntry& a, const MoldenShellEntry& b) {
            if (a.l != b.l) return a.l < b.l;
            if (a.atom_idx != b.atom_idx) return a.atom_idx < b.atom_idx;
            return a.source_order < b.source_order;
        });

        for (size_t i = 0; i < md_array_size(shell_entries); ++i) {
            const MoldenShellEntry& entry = shell_entries[i];

            md_gto_shell_t shell = {
                .atom_idx = entry.atom_idx,
                .primitive_offset = basis.num_primitives,
                .num_primitives = (uint32_t)entry.shell->primitives.size(),
                .l = entry.l,
            };
            md_array_push(basis.shells, shell, arena);
            basis.num_shells += 1;

            const size_t primitive_offset = md_array_size(basis.alpha);
            for (const auto& primitive : entry.shell->primitives) {
                md_array_push(basis.alpha, (float)primitive.exponent, arena);
                md_array_push(basis.coeff, (float)((entry.use_sp_coeff ? primitive.coefficient_sp : primitive.coefficient) * entry.shell->scale_factor), arena);
                basis.num_primitives += 1;
            }

            normalize_shell_coefficients((int)entry.l, basis.alpha + primitive_offset, basis.coeff + primitive_offset, (int)entry.shell->primitives.size());

            const size_t num_funcs = shell_num_spherical_functions((int)entry.l);
            for (size_t func_idx = 0; func_idx < num_funcs; ++func_idx) {
                (void)entry.use_sp_coeff;
                const size_t source_component_idx = shell_component_source_index((int)entry.l, resolved_ao_ordering_mode(), func_idx);
                md_array_push(ordered_to_source_ao, entry.source_ao_offset + (uint32_t)source_component_idx, arena);
            }
        }

        for (const auto& orbital : data.orbitals) {
            if (orbital.coefficients.size() != source_ao_offset) {
                MD_LOG_ERROR("Molden AO coefficient count mismatch: expected %u, got %zu", source_ao_offset, orbital.coefficients.size());
                return false;
            }
        }

        if (md_array_size(ordered_to_source_ao) != source_ao_offset) {
            MD_LOG_ERROR("Molden AO remap size mismatch: expected %u, got %zu", source_ao_offset, md_array_size(ordered_to_source_ao));
            return false;
        }

        return basis.num_shells > 0;
    }

    bool load_file(str_t path) {
        ensure_arena();
        reset();

        std::string error;
        data = molden::parse_molden_file(std::string(path.ptr, path.len), &error);
        if (data.atoms.empty() || data.orbitals.empty()) {
            if (!error.empty()) {
                MD_LOG_ERROR("Failed to load Molden orbital data from '%.*s': %s", (int)path.len, path.ptr, error.c_str());
            }
            return false;
        }

        build_orbital_maps();
        basis_ready = build_basis();
        loaded = basis_ready;
        if (!basis_ready) {
            MD_LOG_ERROR("Failed to initialize Molden orbital basis from '%.*s'", (int)path.len, path.ptr);
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

    const molden::MolecularOrbital* selected_orbital(const ElectronicStructureRepresentation& es) const {
        size_t map_count = 0;
        const int32_t* orbital_map = orbital_map_for_spin(es.spin, &map_count);
        if (!orbital_map || es.orbital_idx < 0 || es.orbital_idx >= (int)map_count) {
            return nullptr;
        }

        const int32_t global_idx = orbital_map[es.orbital_idx];
        if (global_idx < 0 || global_idx >= (int32_t)data.orbitals.size()) {
            return nullptr;
        }
        return &data.orbitals[(size_t)global_idx];
    }

    const molden::MolecularOrbital* orbital_from_local_idx(ElectronicStructureSpin spin, int orbital_idx) const {
        size_t map_count = 0;
        const int32_t* orbital_map = orbital_map_for_spin(spin, &map_count);
        if (!orbital_map || orbital_idx < 0 || orbital_idx >= (int)map_count) {
            return nullptr;
        }

        const int32_t global_idx = orbital_map[orbital_idx];
        if (global_idx < 0 || global_idx >= (int32_t)data.orbitals.size()) {
            return nullptr;
        }
        return &data.orbitals[(size_t)global_idx];
    }

    int homo_idx(ElectronicStructureSpin spin) const {
        return (int)find_homo_idx(data, spin == ElectronicStructureSpin::Beta ? molden::SpinType::Beta : molden::SpinType::Alpha);
    }

    int lumo_idx(ElectronicStructureSpin spin) const {
        const size_t count = num_orbitals(spin);
        return count ? MIN(homo_idx(spin) + 1, (int)count - 1) : -1;
    }

    void reorder_coefficients(double* dst_coeffs, const molden::MolecularOrbital& orbital) const {
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
                (float)(sys->atom.x[i] * ANGSTROM_TO_BOHR),
                (float)(sys->atom.y[i] * ANGSTROM_TO_BOHR),
                (float)(sys->atom.z[i] * ANGSTROM_TO_BOHR),
                1.0f
            );
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
        picking_surface_init(&orb.picking_surface, interaction_surface_molden_orb);
    }

    bool evaluate_orbital_texture(uint32_t tex_id, const md_grid_t& grid, const vec4_t* atom_xyzw, ElectronicStructureSpin spin, int orbital_idx) const {
        const molden::MolecularOrbital* orbital = orbital_from_local_idx(spin, orbital_idx);
        if (!orbital) return false;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        double* ao_coeffs = md_temp_alloc_array(temp, double, md_array_size(ordered_to_source_ao));
        reorder_coefficients(ao_coeffs, *orbital);

        md_gto_grid_evaluate_mo_GL(
            tex_id,
            &grid,
            &basis,
            (const float*)atom_xyzw,
            sizeof(vec4_t),
            ao_coeffs,
            DEFAULT_GTO_CUTOFF_VALUE,
            MD_GTO_EVAL_MODE_PSI,
            MD_GTO_OP_SET
        );
        return true;
    }

    void draw_summary_window(ApplicationState& state) {
        if (!summary.show_window || !loaded) return;

        ImGui::SetNextWindowSize({380, 420}, ImGuiCond_FirstUseEver);
        if (!qm_ui::begin_summary_window("Molden", "MoldenQMSummary", &summary.show_window)) {
            ImGui::End();
            return;
        }

        if (ImGui::TreeNode("File Information")) {
            qm_ui::draw_source_label("Molden");
            ImGui::Text("Title:                 %s", data.title.empty() ? "<none>" : data.title.c_str());
            const MoldenAoOrderingMode resolved_ordering_mode = resolved_ao_ordering_mode();
            if (settings.ao_ordering_mode == MoldenAoOrderingMode::Auto) {
                ImGui::Text("AO Ordering:           Auto (%s)", ao_ordering_mode_str(resolved_ordering_mode));
            } else {
                ImGui::Text("AO Ordering:           %s", ao_ordering_mode_str(settings.ao_ordering_mode));
            }
            if (ImGui::TreeNode("Source Details / Advanced")) {
                ImGui::Text("Coordinate Units:      %s", coordinate_unit_str(data.coord_unit));
                ImGui::Text("Basis Format:          %s", basis_format_str(data.basis_format));
                MoldenAoOrderingMode ao_mode = settings.ao_ordering_mode;
                if (ImGui::BeginCombo("AO Ordering", ao_ordering_mode_str(ao_mode))) {
                    const MoldenAoOrderingMode modes[] = {
                        MoldenAoOrderingMode::Auto,
                        MoldenAoOrderingMode::Identity,
                        MoldenAoOrderingMode::PStandard,
                        MoldenAoOrderingMode::PSpSwapped,
                    };
                    for (MoldenAoOrderingMode mode : modes) {
                        const bool selected = mode == ao_mode;
                        if (ImGui::Selectable(ao_ordering_mode_str(mode), selected)) {
                            settings.ao_ordering_mode = mode;
                            on_ao_ordering_mode_changed(state);
                        }
                        if (selected) ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }
                ImGui::SetItemTooltip("Manual debug override for unusual Molden producers. Leave this on Auto for normal files.");
                ImGui::TreePop();
            }
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("System Information")) {
            ImGui::Text("Num Atoms:             %-6zu", data.atoms.size());
            ImGui::Text("Num Basis Sets:        %-6zu", data.basis_sets.size());
            ImGui::Text("Total Basis Functions: %-6zu", data.total_basis_functions);
            ImGui::Text("Num Alpha Orbitals:    %-6zu", data.num_alpha_orbitals);
            ImGui::Text("Num Beta Orbitals:     %-6zu", data.num_beta_orbitals);
            ImGui::Text("Num Shells:            %-6u", basis.num_shells);
            ImGui::TreePop();
        }

        if (!data.orbitals.empty() && ImGui::TreeNode("Orbital Information")) {
            const int alpha_homo = homo_idx(ElectronicStructureSpin::Alpha);
            const int alpha_lumo = lumo_idx(ElectronicStructureSpin::Alpha);
            ImGui::Text("Alpha HOMO:            %d", alpha_homo + 1);
            ImGui::Text("Alpha LUMO:            %d", alpha_lumo + 1);
            if (has_beta_spin()) {
                const int beta_homo = homo_idx(ElectronicStructureSpin::Beta);
                const int beta_lumo = lumo_idx(ElectronicStructureSpin::Beta);
                ImGui::Text("Beta HOMO:             %d", beta_homo + 1);
                ImGui::Text("Beta LUMO:             %d", beta_lumo + 1);
            }
            ImGui::TreePop();
        }

        if (!data.atoms.empty() && ImGui::TreeNode("Geometry")) {
            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX |
                ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;
            static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;

            if (ImGui::BeginTable("Geometry Table", 6, flags, ImVec2(520, -1), 0)) {
                const ImVec4 hover_color = ImVec4(1.0f, 1.0f, 0.5f, 0.3f);
                const ImVec4 base_color = ImVec4(0.5f, 0.5f, 1.0f, 0.3f);
                ImGui::TableSetupColumn("Atom", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Symbol", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Z", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Coord X", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, hover_color);
                ImGui::PushStyleColor(ImGuiCol_Header, base_color);
                bool item_hovered = false;

                for (size_t row_n = 0; row_n < data.atoms.size(); ++row_n) {
                    const molden::Atom& atom = data.atoms[row_n];
                    const bool is_sel = md_bitfield_test_bit(&state.selection.selection_mask, row_n);
                    const bool is_hov = md_bitfield_test_bit(&state.selection.highlight_mask, row_n);

                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                    ImGui::TableNextColumn();
                    ImGui::PushStyleColor(ImGuiCol_Header, is_hov ? hover_color : base_color);

                    char label[16];
                    snprintf(label, sizeof(label), "%zu", row_n + 1);
                    ImGui::Selectable(label, is_sel || is_hov, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap);
                    if (ImGui::TableGetHoveredRow() == (int)(row_n + 1) && state.mold.sys.atom.count > row_n) {
                        md_bitfield_clear(&state.selection.highlight_mask);
                        md_bitfield_set_bit(&state.selection.highlight_mask, row_n);
                        item_hovered = true;

                        if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                            md_bitfield_set_bit(&state.selection.selection_mask, row_n);
                        } else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                            md_bitfield_clear_bit(&state.selection.selection_mask, row_n);
                        }
                    }

                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(atom.element_symbol.c_str());
                    ImGui::TableNextColumn();
                    ImGui::Text("%d", atom.atomic_number);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.x);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.y);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", atom.z);

                    ImGui::PopStyleColor(1);
                }

                if (!item_hovered && ImGui::IsWindowHovered()) {
                    md_bitfield_clear(&state.selection.highlight_mask);
                }

                ImGui::PopStyleColor(2);
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
        if (!qm_ui::begin_orbital_grid_window("Molden", "MoldenQMOrbitalGrid", &orb.show_window)) {
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

        int orb_homo_idx = unrestricted ? MAX(homo_idx(ElectronicStructureSpin::Alpha), homo_idx(ElectronicStructureSpin::Beta)) : homo_idx(ElectronicStructureSpin::Alpha);
        int orb_lumo_idx = unrestricted ? MIN(lumo_idx(ElectronicStructureSpin::Alpha), lumo_idx(ElectronicStructureSpin::Beta)) : lumo_idx(ElectronicStructureSpin::Alpha);

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
            ImGui::SetItemTooltip("Visual surface threshold. Orbital grids render both +isovalue and -isovalue surfaces.");

            orb.iso.values[0] = (float)iso_val;
            orb.iso.values[1] = -(float)iso_val;
            orb.iso.count = 2;

            ImGui::ColorEdit4("##Color Positive", orb.iso.colors[0].elem);
            ImGui::SetItemTooltip("Color Positive");
            ImGui::ColorEdit4("##Color Negative", orb.iso.colors[1].elem);
            ImGui::SetItemTooltip("Color Negative");

            if (ImGui::IsWindowAppearing()) {
                orb.scroll_to_idx = orb.mo_idx;
            }
            if (ImGui::Button("Goto HOMO", ImVec2(-1, 0))) {
                orb.scroll_to_idx = homo_idx(ElectronicStructureSpin::Alpha);
                orb.mo_idx = homo_idx(ElectronicStructureSpin::Alpha);
            }

            const ImGuiTableFlags flags =
                ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY;
            const int num_cols = unrestricted ? 5 : 3;
            if (ImGui::BeginTable("Molecular Orbitals", num_cols, flags)) {
                if (unrestricted) {
                    ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn((const char*)u8"Occ. α", ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn((const char*)u8"Occ. β", ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn((const char*)u8"Ene. α", ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn((const char*)u8"Ene. β", ImGuiTableColumnFlags_WidthFixed);
                } else {
                    ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn("Occupancy", ImGuiTableColumnFlags_WidthFixed);
                    ImGui::TableSetupColumn("Energy", ImGuiTableColumnFlags_WidthFixed);
                }
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
                    if (occ_alpha && n < (int)alpha_count) {
                        const char* lbl = "";
                        if (n == orb_homo_idx && n == homo_idx(ElectronicStructureSpin::Alpha)) lbl = "HOMO";
                        else if (n == orb_lumo_idx && n == lumo_idx(ElectronicStructureSpin::Alpha)) lbl = "LUMO";
                        ImGui::Text("%.1f %s", occ_alpha[n], lbl);
                    } else {
                        ImGui::TextUnformatted("-");
                    }

                    if (unrestricted) {
                        ImGui::TableNextColumn();
                        if (occ_beta && n < (int)beta_count) {
                            const char* lbl = "";
                            if (n == orb_homo_idx && n == homo_idx(ElectronicStructureSpin::Beta)) lbl = "HOMO";
                            else if (n == orb_lumo_idx && n == lumo_idx(ElectronicStructureSpin::Beta)) lbl = "LUMO";
                            ImGui::Text("%.1f %s", occ_beta[n], lbl);
                        } else {
                            ImGui::TextUnformatted("-");
                        }
                    }

                    ImGui::TableNextColumn();
                    if (ene_alpha && n < (int)alpha_count) {
                        ImGui::Text("%.4f", ene_alpha[n]);
                    } else {
                        ImGui::TextUnformatted("-");
                    }

                    if (unrestricted) {
                        ImGui::TableNextColumn();
                        if (ene_beta && n < (int)beta_count) {
                            ImGui::Text("%.4f", ene_beta[n]);
                        } else {
                            ImGui::TextUnformatted("-");
                        }
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
            if (!found) {
                job_queue[num_jobs++] = i;
            }
        }

        if (num_jobs > 0) {
            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };

            vec4_t* atom_xyzw = build_atom_positions_bohr(&state.mold.sys, md_temp_allocator(temp));
            AABB aabb = {};
            calculate_bounds(&aabb, atom_xyzw, md_array_size(atom_xyzw));

            md_grid_t grid = {};
            init_grid(&grid, aabb.min_ext, aabb.max_ext, orbital_grid_samples_per_angstrom * BOHR_TO_ANGSTROM);

            for (int i = 0; i < num_jobs; ++i) {
                const int slot_idx = job_queue[i];
                const int mo_idx = vol_mo_idx[slot_idx];
                const ElectronicStructureSpin spin = vol_mo_spin[slot_idx];
                orb.vol_mo_idx[slot_idx] = mo_idx;
                orb.vol_mo_spin[slot_idx] = spin;

                if (orbital_from_local_idx(spin, mo_idx)) {
                    init_volume(&orb.vol[slot_idx], grid);
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
        const mat4_t clip_to_world = mat4_mul(view_to_world, clip_to_view);

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
                InteractionSurfaceState surface_state = interaction_surface(interaction_surface_molden_orb, vec_cast(sz), InteractionSurfaceFlags_NoRegionSelect);
                InteractionSurfaceViewTransformArgs view_args = {
                    .camera = orb.camera,
                    .trackball_param = state.view.trackball_param,
                };
                InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&orb.target, surface_state, view_args);
                if (view_result.reset_requested) {
                    orb.target = orb.default_view;
                }

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

        for (int x = 1; x < num_x; ++x) {
            ImVec2 p0 = {canvas_min.x + orb_win_sz.x * x, canvas_min.y};
            ImVec2 p1 = {canvas_min.x + orb_win_sz.x * x, canvas_max.y};
            draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
        }
        for (int y = 1; y < num_y; ++y) {
            ImVec2 p0 = {canvas_min.x, canvas_min.y + orb_win_sz.y * y};
            ImVec2 p1 = {canvas_max.x, canvas_min.y + orb_win_sz.y * y};
            draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
        }

        int width = MAX(1, (int)(orb_win_sz.x * io.DisplayFramebufferScale.x));
        int height = MAX(1, (int)(orb_win_sz.y * io.DisplayFramebufferScale.y));
        if ((int)orb.gbuf.width != width || (int)orb.gbuf.height != height) {
            init_gbuffer(&orb.gbuf, width, height);
            for (int i = 0; i < num_mos; ++i) {
                gl::init_texture_2D(orb.iso_tex + i, width, height, GL_RGBA8);
            }
        }

        if (orb.show_coordinate_system_widget) {
            float ext = MIN(orb_win_sz.x, orb_win_sz.y) * 0.4f;
            float pad = 20.0f;
            ImVec2 size = ImVec2(ext, ext);
            ImGui::SetCursorScreenPos(ImVec2(canvas_min.x + pad, canvas_max.y - ext - pad));
            quat_t out_orientation = orb.target.orientation;
            if (ImGui::CoordinateSystemWidget(&out_orientation, orb.camera.orientation, size)) {
                const vec3_t look_at = camera_get_look_at(orb.target);
                orb.target.orientation = quat_normalize(out_orientation);
                orb.target.position = camera_position_from_look_at(look_at, orb.target.orientation, orb.target.distance);
            }
        }

        const mat4_t view_mat = camera_world_to_view_matrix(orb.camera);
        const mat4_t proj_mat = camera_view_to_clip_matrix_persp(orb.camera, aspect);
        const mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(orb.camera, aspect);

        clear_gbuffer(&orb.gbuf);
        const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
            GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY};

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
            .draw_operations = {
                .count = 1,
                .ops = &draw_op,
            },
            .view_transform = {
                .view_matrix = (const float*)view_mat.elem,
                .proj_matrix = (const float*)proj_mat.elem,
            },
            .picking_offset = {
                .atom_base = state.picking_range_atom.beg,
                .bond_base = state.picking_range_bond.beg,
            },
        };
        md_gl_draw(&draw_args);

        glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
        glClearColor(1, 1, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT);

        postprocessing::Descriptor postprocess_desc = {
            .background = {
                .color = {24.f, 24.f, 24.f},
            },
            .tonemapping = {
                .enabled = state.visuals.tonemapping.enabled,
                .mode = state.visuals.tonemapping.tonemapper,
                .exposure = state.visuals.tonemapping.exposure,
                .gamma = state.visuals.tonemapping.gamma,
            },
            .ambient_occlusion = {
                .enabled = false,
            },
            .depth_of_field = {
                .enabled = false,
            },
            .fxaa = {
                .enabled = true,
            },
            .temporal_aa = {
                .enabled = false,
            },
            .sharpen = {
                .enabled = false,
            },
            .input_textures = {
                .depth = orb.gbuf.tex.depth,
                .color = orb.gbuf.tex.color,
                .normal = orb.gbuf.tex.normal,
                .velocity = orb.gbuf.tex.velocity,
            },
        };

        ViewParam view_param = {
            .matrix = {
                .curr = {
                    .view = view_mat,
                    .proj = proj_mat,
                    .norm = view_mat,
                },
                .inv = {
                    .proj = inv_proj_mat,
                },
            },
            .clip_planes = {
                .near = orb.camera.near_plane,
                .far = orb.camera.far_plane,
            },
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
                .render_target = {
                    .depth = orb.gbuf.tex.depth,
                    .color = orb.iso_tex[i],
                    .width = orb.gbuf.width,
                    .height = orb.gbuf.height,
                    .clear_color = true,
                },
                .texture = {
                    .density_volume = orb.vol[i].tex_id,
                },
                .matrix = {
                    .model = orb.vol[i].texture_to_world,
                    .view = view_mat,
                    .proj = proj_mat,
                    .inv_proj = inv_proj_mat,
                },
                .iso = {
                    .enabled = true,
                    .count = orb.iso.count,
                    .values = orb.iso.values,
                    .colors = orb.iso.colors,
                },
                .shading = {
                    .env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f,
                    .roughness = 0.3f,
                    .dir_radiance = {10, 10, 10},
                    .ior = 1.5f,
                },
                .voxel_spacing = orb.vol[i].voxel_size,
            };
            volume::render_volume(vol_desc);
        }

        ImGui::End();
    }

    bool evaluate_molecular_orbital(EvalElectronicStructure& eval) const {
        Representation* rep = eval.rep;
        const ElectronicStructureRepresentation& es = rep->electronic_structure;
        const molden::MolecularOrbital* orbital = selected_orbital(es);
        if (!orbital) {
            MD_LOG_ERROR("Invalid Molden orbital selection: spin=%d idx=%d", (int)es.spin, es.orbital_idx);
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        md_allocator_i* temp_alloc = md_temp_allocator(temp);
        vec4_t* atom_xyzw = build_atom_positions_bohr(eval.sys, temp_alloc);
        AABB aabb = {};
        calculate_bounds(&aabb, atom_xyzw, md_array_size(atom_xyzw));

        const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)es.resolution];
        const double samples_per_unit_length = samples_per_angstrom * BOHR_TO_ANGSTROM;

        md_grid_t grid = {};
        init_grid(&grid, aabb.min_ext, aabb.max_ext, samples_per_unit_length);
        init_volume(&rep->electronic_structure.density_vol, grid, GL_R32F);

        double* ao_coeffs = md_temp_alloc_array(temp, double, md_array_size(ordered_to_source_ao));
        reorder_coefficients(ao_coeffs, *orbital);

        md_gto_grid_evaluate_mo_GL(
            rep->electronic_structure.density_vol.tex_id,
            &grid,
            &basis,
            (const float*)atom_xyzw,
            sizeof(vec4_t),
            ao_coeffs,
            DEFAULT_GTO_CUTOFF_VALUE,
            MD_GTO_EVAL_MODE_PSI,
            gto_op_from_use_magnitude(es.use_magnitude)
        );

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
                if (loaded && ImGui::BeginMenu("Molden")) {
                    ImGui::Checkbox(qm_ui::summary_label, &summary.show_window);
                    ImGui::Checkbox(qm_ui::orbital_grid_label, &orb.show_window);
                    ImGui::EndMenu();
                }
                break;

            case viamd::EventType_ViamdLoadData: {
                auto& payload = *(LoadDataPayload*)e.payload;
                if (is_molden_path(payload.path_to_file)) {
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

                fill_orbitals(&info.alpha, data, molden::SpinType::Alpha, info.alloc);
                fill_orbitals(&info.beta, data, molden::SpinType::Beta, info.alloc);

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
                if (es.source != ElectronicStructureSource::MolecularOrbital) {
                    break;
                }

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

static MoldenComponent instance = {};
