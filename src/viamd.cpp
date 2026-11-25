
#include <md_infer.h>
#include <md_element.h>
#include <core/md_geom.h>
#include <md_run.h>
#include <md_dssp.h>
#include <md_qm.h>
#include "density_cache.h"
#include <md_gto_iso.h>
#include <md_vial.h>
#include <md_pucker.h>
#include <vector>
#include <algorithm>

#include <atomic>
#include <stdarg.h>

#include <script_engine.h>
#include <attribute_access.h>
#include <script_suggestions.h>

#include <core/md_log.h>
#include <core/md_path.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_file.h>
#include <core/md_writer.h>
#include <core/md_arena_allocator.h>
#include <core/md_vec_math.h>

#include <viamd.h>
#include <drop_files.h>
#include <representation.h>
#include <system_scene.h>
#include <script_editor.h>
#include <viamd_event.h>
#include <event.h>
#include <loader.h>
#include <color_utils.h>
#include <app_settings.h>
#include <display_units.h>
#include <frame_axis.h>
#include <serialization_utils.h>
#include <plot_periodic.h>

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>

#include <imgui.h>
#include <imgui_notify.h>
#include <es_planes.h>
#include <implot.h>
#include <implot_internal.h>

md_float3x3_t mat3_PCA(const md_float4_t* xyzw, size_t count) {
    md_float3_t acc = md_float3_zero();
    for (size_t i = 0; i < count; ++i) {
        acc = md_float3_add(acc, md_float3_from_float4(xyzw[i]));
    }
    md_float3_t mean = acc / (float)count;

    md_float3x3_t cov = md_float3x3_covariance_matrix_float4(xyzw, nullptr, count, mean);
    md_float3x3_eigen_t eigen = md_float3x3_eigen(cov);
    md_float3x3_t PCA = md_float3x3_orthonormalize(md_float3x3_extract_rotation(eigen.vectors));
    return PCA;
}

void calculate_bounds(float out_min[3], float out_max[3], const md_float4_t* xyzw, size_t count, const md_float3x3_t& orientation) {
    md_float4_t min_v = md_float4_set1(FLT_MAX);
    md_float4_t max_v = md_float4_set1(-FLT_MAX);

    md_float4x4_t rot = md_float4x4_from_float3x3(md_float3x3_transpose(orientation));

    for (size_t i = 0; i < count; ++i) {
        md_float4_t v = md_float4x4_mul_float4(rot, xyzw[i]);
        min_v = md_float4_min(min_v, v);
        max_v = md_float4_max(max_v, v);
    }

    // Padding
    const float pad = 6.0f;
    min_v -= pad;
    max_v += pad;

    md_memcpy(out_min, &min_v, sizeof(float) * 3);
    md_memcpy(out_max, &max_v, sizeof(float) * 3);
}

// Construct texture to world transformation matrix for Volume
// extent is the extent of the volume (dim * voxel_size)
md_float4x4_t compute_texture_to_world_mat(const md_float3x3_t& orientation, const md_float3_t& origin, const md_float3_t& extent) {
    md_float4x4_t T = md_float4x4_translate_float3(origin);
    md_float4x4_t R = md_float4x4_from_float3x3(orientation);
    md_float4x4_t S = md_float4x4_scale_float3(extent);
    return T * R * S;
}

md_float4x4_t compute_world_to_model_mat(const md_float3x3_t& orientation, const md_float3_t& origin) {
    md_float4x4_t world_to_model = md_float4x4_from_float3x3(md_float3x3_transpose(orientation)) * md_float4x4_translate_float3(-origin);
    return world_to_model;
}

md_float4x4_t compute_index_to_world_mat(const md_float3x3_t& orientation, const md_float3_t& in_origin, const md_float3_t& stepsize) {
    md_float3_t step_x = orientation.col[0] * stepsize.x;
    md_float3_t step_y = orientation.col[1] * stepsize.y;
    md_float3_t step_z = orientation.col[2] * stepsize.z;
    // Shift origin by half voxel
    md_float3_t origin = in_origin + orientation * (stepsize * 0.5f);

    md_float4x4_t index_to_world = {
        step_x.x, step_x.y, step_x.z, 0.0f,
        step_y.x, step_y.y, step_y.z, 0.0f,
        step_z.x, step_z.y, step_z.z, 0.0f,
        origin.x, origin.y, origin.z, 1.0f,
    };

    return index_to_world;
}

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
void compute_dim(int out_dim[3], const md_float3_t& in_ext, double samples_per_unit_length) {
    out_dim[0] = MD_CLAMP(MD_ALIGN_TO((int)(in_ext.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = MD_CLAMP(MD_ALIGN_TO((int)(in_ext.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = MD_CLAMP(MD_ALIGN_TO((int)(in_ext.z * samples_per_unit_length), 8), 8, 512);
}

// Grid units are BOHR, matching what md_gto evaluates in; a Volume's transforms are Angstrom,
// matching the world the camera lives in. init_volume is where the two meet, and it is the only
// place that conversion belongs.
void init_grid(md_grid_t* grid, const md_float3x3_t& orientation, const md_float3_t& min_ext, const md_float3_t& max_ext, double samples_per_unit_length) {
    MD_ASSERT(grid);
    md_float3_t extent = max_ext - min_ext;
    compute_dim(grid->dim, extent, samples_per_unit_length);
    md_float3_t voxel_size = md_float3_div(extent, md_float3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
    grid->orientation = orientation;
    grid->origin = orientation * min_ext;
    grid->spacing = voxel_size;
}

void init_volume(Volume* vol, const md_grid_t& grid, GLenum format) {
    MD_ASSERT(vol);
    md_memcpy(vol->dim, grid.dim, sizeof(vol->dim));

    // A grid is in BOHR and a Volume's transforms are in Angstrom; this is where the two meet.
    const float scl = (float)BOHR_TO_ANGSTROM;

    md_float3_t extent = md_grid_extent(&grid);
    vol->world_to_model   = compute_world_to_model_mat(grid.orientation, grid.origin * scl);
    vol->texture_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
    vol->voxel_size       = grid.spacing * scl;
    gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
    volume::notify_data_changed(vol->tex_id);
}

// PICKING TOOLTIP
//
// The tooltip describes the object under the cursor: an atom, a bond, or a backbone segment of the cartoon or ribbons,
// which stands for its residue. What a click selects is the highlight's to show, as that follows the selection
// granularity; the tooltip does not change with it, so the same atom always reads the same.
//
// The title starts with what the object is, the word VIAL selects it by (particle, bond, ring, residue), then names it
// and where it sits (residue, chain): with a wider granularity, a residue's or a chain's atoms lit up, the title still
// says that it is a particle that is described. Each row gives one property, grouped by what it means
// rather than by the bitfield it is stored in: kind, role, chemistry, region, source. Keys, and values that are our own
// words, are lower case; names from the file (atoms, residues, chains, entities) and element symbols are printed as
// given. A row with nothing to say is left out. The index row spells the object in the script's terms, from 1.
//
// The font carries Latin-1, Greek and U+2010-2027 (application.cpp): the middle dot, the en dash, Å, °, ², φ and ψ
// print; an arrow, the minus sign U+2212 and the prime U+2032 do not, so '-', '/' and "'" stand in for them.

// The coordinates of the atoms of mask with their masses (w), packed in the order of the atoms
static void extract_xyzw_from_mask(md_float4_t* out_xyzw, const md_bitset_t* mask, const md_system_t* sys, const md_system_state_t* state) {
    MD_ASSERT(out_xyzw);
    MD_ASSERT(mask);
    MD_ASSERT(sys);
    md_bitset_iter_t it = md_bitset_iter(*mask);
    size_t count = 0;
    while (md_bitset_iter_next(&it)) {
        const size_t i = it.idx;
        out_xyzw[count++] = md_float4_from_float3(md_state_coord(state, i), md_atom_mass(&sys->atom, i));
    }
}

static md_str_t tooltip_vprintf(md_allocator_i* alloc, const char* format, va_list args) {
    va_list args_len;
    va_copy(args_len, args);
    const int len = vsnprintf(NULL, 0, format, args_len);
    va_end(args_len);
    if (len <= 0) return {};

    char* buf = (char*)md_alloc(alloc, (size_t)len + 1);
    vsnprintf(buf, (size_t)len + 1, format, args);
    return {buf, (size_t)len};
}

void tooltip_title(PickingTooltipTextRequest* req, const char* format, ...) {
    MD_ASSERT(req && req->alloc);
    va_list args;
    va_start(args, format);
    const TooltipLine line = {.kind = TooltipLineKind::Title, .text = tooltip_vprintf(req->alloc, format, args)};
    va_end(args);
    md_array_push(req->lines, line, req->alloc);
}

void tooltip_row(PickingTooltipTextRequest* req, md_str_t key, const char* format, ...) {
    MD_ASSERT(req && req->alloc);
    va_list args;
    va_start(args, format);
    const TooltipLine line = {.kind = TooltipLineKind::Row, .key = key, .text = tooltip_vprintf(req->alloc, format, args)};
    va_end(args);
    md_array_push(req->lines, line, req->alloc);
}

// A row of a list built in sb, left out when the list is empty
static void tooltip_row_list(PickingTooltipTextRequest* req, md_str_t key, const md_str_writer_t& sb) {
    if (!md_str_writer_empty(&sb)) {
        tooltip_row(req, key, "%s", md_str_writer_to_cstr(&sb));
    }
}

static void tooltip_list_add(md_str_writer_t* sb, const char* item) {
    if (!md_str_writer_empty(sb)) *sb += ", ";
    *sb += item;
}

// A lower case copy, for names the system spells capitalised ("Carbon") where the tooltip uses them as words
static md_str_t tooltip_lower(md_allocator_i* alloc, md_str_t str) {
    if (md_str_empty(str)) return {};
    char* buf = (char*)md_alloc(alloc, str.len + 1);
    md_memcpy(buf, str.ptr, str.len);
    buf[str.len] = '\0';
    md_chars_to_lower(buf, str.len);
    return {buf, str.len};
}

// " e", " Å", or nothing for a dimensionless value
static void tooltip_unit_suffix(char* buf, size_t cap, md_unit_t unit) {
    MD_ASSERT(buf && cap > 1);
    buf[0] = '\0';
    if (md_unit_print(buf + 1, cap - 1, unit) > 0) {
        buf[0] = ' ';
    }
}

// Where an atom sits: its residue and its part, -1 for either it is not in
struct TooltipPlace {
    int res_idx = -1;
    int part_idx = -1;
    int atom_idx = -1;      // An atom of it, which tells its molecule
};

static TooltipPlace tooltip_place_of_atom(const md_system_t& sys, size_t atom_idx) {
    return {
        .res_idx = md_system_residue_find_by_atom_idx(&sys, atom_idx),
        .part_idx = md_system_part_find_by_atom_idx(&sys, atom_idx),
        .atom_idx = (int)atom_idx,
    };
}

static TooltipPlace tooltip_place_of_residue(const md_system_t& sys, int res_idx) {
    return {
        .res_idx = res_idx,
        .part_idx = res_idx >= 0 ? md_system_part_find_by_residue_idx(&sys, (size_t)res_idx) : -1,
    };
}

// "ALA 42", or the residue's index when it has no name
static void tooltip_append_residue(md_str_writer_t* sb, const md_system_t& sys, int res_idx) {
    const md_str_t name = md_system_residue_name(&sys, res_idx);
    if (!md_str_empty(name)) {
        md_str_writer_fmt(sb, MD_STR_FMT " %d", MD_STR_ARG(name), md_system_residue_seq_id(&sys, res_idx));
    } else {
        md_str_writer_fmt(sb, "residue %d", res_idx + 1);
    }
}

// The chain of a place (md_chain_data_t), -1 for one in none
static int tooltip_place_chain(const md_system_t& sys, const TooltipPlace& place) {
    return place.res_idx >= 0 ? md_system_chain_find_by_residue_idx(&sys, (size_t)place.res_idx) : -1;
}

// The position of the part among the strands of the polymers, as polymer(i) counts them: -1 for a part which is no strand
static int tooltip_polymer_position(const md_system_t& sys, int part_idx) {
    if (part_idx < 0 || !md_part_kind_is_polymer(md_system_part_kind(&sys, (size_t)part_idx))) return -1;
    int pos = 0;
    for (int p = 0; p < part_idx; ++p) pos += md_part_kind_is_polymer(md_system_part_kind(&sys, (size_t)p)) ? 1 : 0;
    return pos;
}

// "chain A", the id the data gives the chain (md_chain_data_t: a chain is only where the data names one)
static void tooltip_append_chain(md_str_writer_t* sb, const md_system_t& sys, int chain) {
    md_str_writer_fmt(sb, "chain " MD_STR_FMT, MD_STR_ARG(md_system_chain_id(&sys, (size_t)chain)));
}

// The rest of a title: " · ALA 42 · chain A" for one place, and for the two ends of a bond that crosses over
// " · ALA 42 / GLY 43 · chain A" or " · CYS 10 · chain A / CYS 50 · chain B". A place in no chain is named by its
// residue alone; where it is in a script's terms (polymer(2), molecule(7)) is the index row's.
static void tooltip_append_place(md_str_writer_t* sb, const md_system_t& sys, const TooltipPlace& a, const TooltipPlace* b = nullptr) {
    auto sep = [sb]() { if (!md_str_writer_empty(sb)) *sb += TOOLTIP_SEP; };
    const int chain_a = tooltip_place_chain(sys, a);
    const int chain_b = b ? tooltip_place_chain(sys, *b) : -1;
    const bool two_chains = b && chain_b != chain_a;
    const bool two_residues = b && b->res_idx != a.res_idx;

    if (a.res_idx != -1) {
        sep();
        tooltip_append_residue(sb, sys, a.res_idx);
        if (two_residues && !two_chains && b->res_idx != -1) {
            *sb += " / ";
            tooltip_append_residue(sb, sys, b->res_idx);
        }
    }
    if (chain_a != -1) {
        sep();
        tooltip_append_chain(sb, sys, chain_a);
    }
    if (two_chains && (b->res_idx != -1 || chain_b != -1)) {
        *sb += " / ";
        if (b->res_idx != -1) {
            tooltip_append_residue(sb, sys, b->res_idx);
            if (chain_b != -1) *sb += TOOLTIP_SEP;
        }
        if (chain_b != -1) {
            tooltip_append_chain(sb, sys, chain_b);
        }
    }
}

// Two units: atom(3) atom(9) in v1, and particle[[2, 8]] in VIAL, which measurements take as their one argument
static void tooltip_index_pair(md_str_writer_t* sb, script::Language lang, const char* kind, int64_t a, int64_t b) {
    if (lang == script::Language::V2) {
        const char* sep = md_str_writer_empty(sb) ? "" : " ";
        md_str_writer_fmt(sb, "%s%s[[%lld, %lld]]", sep, strcmp(kind, "atom") ? kind : "particle", (long long)a, (long long)b);
    } else {
        script_write_position(sb, lang, kind, a);
        script_write_position(sb, lang, kind, b);
    }
}

// "residue(56) part(3) chain(1)", after whatever sb holds: the place in the script's terms, residue[55] part[2] chain[0]
// in VIAL. The part is given where it is
// less than its molecule (a glycan, a strand bonded to another). Then what it is filed under: its chain, else the strand
// of a polymer it is by its position (polymer(2)), else its molecule by its position among the molecules.
// Positions are written as the script's language writes them (script_write_position): atom(5) or particle[4]
static void tooltip_append_index(md_str_writer_t* sb, const md_system_t& sys, const TooltipPlace& place, script::Language lang) {
    const int chain = tooltip_place_chain(sys, place);
    int atom = place.atom_idx;
    if (atom == -1 && place.res_idx != -1) atom = (int)md_system_residue_atom_range(&sys, place.res_idx).beg;
    const md_molecule_idx_t mol = atom != -1 ? md_system_molecule_find_by_atom_idx(&sys, (size_t)atom) : -1;

    if (place.res_idx != -1) script_write_position(sb, lang, "residue", place.res_idx);
    const int32_t part = place.part_idx;
    if (part != -1 && mol != -1) {
        const md_urange_t part_atoms = md_system_part_atom_range(&sys, (size_t)part);
        const size_t mol_size = md_system_molecule_atom_count(&sys, (size_t)mol);
        if (part_atoms.end - part_atoms.beg < mol_size) {
            script_write_position(sb, lang, "part", part);
        }
    }
    if (chain != -1) {
        script_write_position(sb, lang, "chain", chain);
        return;
    }
    const int pos = tooltip_polymer_position(sys, part);
    if (pos >= 0) {
        script_write_position(sb, lang, "polymer", pos);
        return;
    }
    if (mol != -1) script_write_position(sb, lang, "molecule", mol);
}

// What a place is of: its part, by its name and its kind as part(name=, kind=) take them (the name its type's
// description, an mmCIF entity's, or inferred: the kind of a polymer, the residue names of anything else), and its
// molecule, by its type, where that says more than the part (a glycoprotein, a molecule of a topology, named as the
// topology names it)
static void tooltip_part_and_molecule_rows(PickingTooltipTextRequest* req, const md_system_t& sys, const TooltipPlace& place) {
    md_str_t part_name = {};
    if (place.part_idx != -1) {
        const md_str_t name = md_vial_part_name(&sys, (size_t)place.part_idx);
        const char* kind = md_vial_part_kind(&sys, (size_t)place.part_idx);
        if (md_str_empty(name) || md_str_eq_cstr(name, kind)) {
            tooltip_row(req, MD_STR_LIT("part"), "%s", kind);
        } else {
            tooltip_row(req, MD_STR_LIT("part"), MD_STR_FMT " (%s)", MD_STR_ARG(name), kind);
        }
        part_name = name;
    }
    int atom = place.atom_idx;
    if (atom == -1 && place.res_idx != -1) atom = (int)md_system_residue_atom_range(&sys, place.res_idx).beg;
    const md_molecule_idx_t mol = atom != -1 ? md_system_molecule_find_by_atom_idx(&sys, (size_t)atom) : -1;
    if (mol != -1) {
        const md_molecule_type_idx_t type = md_system_molecule_type_idx(&sys, (size_t)mol);
        const md_str_t name = type != -1 ? md_system_molecule_type_name(&sys, (size_t)type) : md_str_t{};
        if (!md_str_empty(name) && !md_str_eq(name, part_name)) {
            tooltip_row(req, MD_STR_LIT("molecule"), MD_STR_FMT, MD_STR_ARG(name));
        }
    }
}

// A fragment (md_fragment.h) by its name, as fragment(name) selects it: "phenyl", "peptide", "C3 linker". Its kind,
// which a script does not select by, only where it has no name
static void tooltip_append_fragment(md_str_writer_t* sb, const md_system_t& sys, md_fragment_idx_t frag) {
    const md_str_t name = md_system_fragment_name(&sys, frag);
    if (!md_str_empty(name)) {
        *sb += name;
    } else {
        *sb += md_str_from_cstr(md_fragment_kind_str(md_system_fragment_kind(&sys, frag)));
    }
}

// The fragment row of an atom: what chemical piece of its molecule it is in, and, as fragments and residues do not
// nest, the other residues that piece is of: "ester, with POPC 12". The fragments are of the
// system's bonds, which the frame shown may have changed.
static void tooltip_fragment_row(PickingTooltipTextRequest* req, const md_system_t& sys, size_t atom_idx) {
    const md_fragment_idx_t frag = md_system_fragment_find_by_atom_idx(&sys, atom_idx);
    if (frag == MD_FRAGMENT_NONE) return;

    md_str_writer_t sb = md_str_writer_create(req->alloc);
    tooltip_append_fragment(&sb, sys, frag);

    const md_residue_idx_t own = md_system_residue_find_by_atom_idx(&sys, atom_idx);
    md_residue_idx_t res[4];
    const size_t num_res = md_system_fragment_residues(res, MD_COUNTOF(res), &sys, frag);
    size_t shown = 0;
    for (size_t i = 0; i < MD_MIN(num_res, MD_COUNTOF(res)); ++i) {
        if (res[i] == own) continue;
        sb += shown == 0 ? MD_STR_LIT(", with ") : MD_STR_LIT(", ");
        tooltip_append_residue(&sb, sys, res[i]);
        shown += 1;
    }
    // Of more residues than were looked up (a residue of its own is among them)
    if (num_res > MD_COUNTOF(res)) {
        md_str_writer_fmt(&sb, " and %zu more", num_res - MD_COUNTOF(res));
    }
    tooltip_row_list(req, MD_STR_LIT("fragment"), sb);
}

// The atom's name, else its element's symbol, else its index
static void tooltip_append_atom_name(md_str_writer_t* sb, const md_system_t& sys, size_t atom_idx) {
    md_str_t name = md_atom_name(&sys.atom, atom_idx);
    if (md_str_empty(name)) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom_idx);
        if (z) name = md_atomic_number_symbol(z);
    }
    if (!md_str_empty(name)) {
        *sb += name;
    } else {
        md_str_writer_fmt(sb, "atom %zu", atom_idx + 1);
    }
}

// An attribute's row as a script reads it: bfactor(s) for the B-factors, attr(name, of=s) with the name below atom/ for
// the others; the label, which is the representation's, only for one outside atom/
static md_str_t tooltip_attribute_key(const md_attribute_t* attr) {
    if (md_str_eq(attr->path, MD_STR_LIT("atom/b_factor"))) return MD_STR_LIT("bfactor");
    const md_str_t atom = MD_STR_LIT("atom/");
    if (attr->path.len > atom.len && md_str_begins_with(attr->path, atom)) return md_str_substr(attr->path, atom.len);
    return atom_attribute_label(attr);
}

// Values of the attribute table for one atom: the field a visible representation is coloured by, at the variant it
// shows, which is the value behind the colour; then every field with one value per atom (a partial charge, a B-factor,
// an occupancy). Fields the atom's rows already give from the system itself are skipped, as are absent values.
static void tooltip_atom_attribute_rows(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t atom_idx) {
    const md_system_t& sys = sc.data.sys;

    struct Field {
        const md_attribute_t* attr;
        int variant;    // -1 for a field without a variant axis
    };
    Field fields[6];
    size_t num_fields = 0;

    auto add = [&](const md_attribute_t* attr, int variant) {
        if (!attr || num_fields == MD_COUNTOF(fields)) return;
        for (size_t i = 0; i < num_fields; ++i) {
            if (fields[i].attr == attr && fields[i].variant == variant) return;
        }
        fields[num_fields++] = {attr, variant};
    };

    for (size_t i = 0; i < md_array_size(sc.representation.reps); ++i) {
        const Representation& rep = sc.representation.reps[i];
        if (!representation_shown(sc, rep) || rep.color_mapping != ColorMapping::Attribute) continue;
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, rep.atom_attribute.key);
        if (!attr) continue;
        const int variant = attr->format.rank > 1 ? MD_CLAMP(rep.atom_attribute.variant_idx, 0, atom_attribute_variant_count(attr) - 1) : -1;
        add(attr, variant);
    }

    md_attribute_id_t ids[32];
    const size_t num_ids = MD_MIN(atom_attribute_query(ids, MD_COUNTOF(ids), sys), MD_COUNTOF(ids));
    for (size_t i = 0; i < num_ids; ++i) {
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, ids[i]);
        if (!attr || attr->format.rank != 1) continue;
        if (md_str_eq(attr->path, MD_STR_LIT("atom/formal_charge")) || md_str_eq(attr->path, MD_STR_LIT("atom/mass"))) continue;
        add(attr, -1);
    }

    for (size_t i = 0; i < num_fields; ++i) {
        const md_attribute_t* attr = fields[i].attr;
        const int variant = fields[i].variant;
        const md_attribute_slice_t slice = variant >= 0 ? md_attribute_slice_2((uint32_t)variant, atom_idx) : md_attribute_slice_1(atom_idx);

        float value = 0.0f;
        if (md_attribute_extract_f32(&value, 1, &sys.attributes, attr, slice, md_unit_none()) != 1 || atom_attribute_value_absent(value)) continue;

        md_str_t key = tooltip_attribute_key(attr);
        if (variant >= 0) {
            key = md_alloc_printf(req->alloc, MD_STR_FMT " (%d)", MD_STR_ARG(key), variant + 1);
        }
        char unit[32];
        tooltip_unit_suffix(unit, sizeof(unit), attr->unit);
        // A charge reads with its sign and to a thousandth of e; anything else by its significant digits
        const bool charge = md_unit_base_equal(attr->unit, md_unit_elementary_charge());
        tooltip_row(req, key, charge ? "%+.3f%s" : "%.4g%s", value, unit);
    }
}

static void tooltip_atom(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t atom_idx) {
    const md_system_t& sys = sc.data.sys;
    const TooltipPlace place = tooltip_place_of_atom(sys, atom_idx);

    // Titled by what it is, as every tooltip is: "particle · CA · ALA5 · A". The kind row says which: atom, bead, virtual
    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += MD_STR_LIT("particle" TOOLTIP_SEP);
    tooltip_append_atom_name(&title, sys, atom_idx);
    tooltip_append_place(&title, sys, place);
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    const md_atomic_number_t z        = md_atom_atomic_number(&sys.atom, atom_idx);
    const md_particle_kind_t particle = md_atom_particle_kind(&sys.atom, atom_idx);
    // Its chemistry in the frame shown where the bonds tracked changed it, that of the system elsewhere
    int formal_charge = 0;
    int num_h = 0;
    const md_atom_flags_t    flags    = frame_atom_flags(sc, atom_idx, &formal_charge, &num_h);
    const bool nucleotide = place.res_idx != -1 && md_system_residue_kind(&sys, place.res_idx) == MD_RESIDUE_KIND_NUCLEOTIDE;

    // What it is, its kind as particle(kind=...) takes it: atom, bead or virtual, as the kind row of a bond or a ring
    tooltip_row(req, MD_STR_LIT("kind"), "%s", particle == MD_PARTICLE_VIRTUAL_SITE ? "virtual" : md_particle_kind_name(particle));
    if (z) {
        const md_str_t symbol = md_atomic_number_symbol(z);
        const md_str_t name   = tooltip_lower(req->alloc, md_atomic_number_name(z));
        tooltip_row(req, MD_STR_LIT("element"), MD_STR_FMT " (" MD_STR_FMT ")", MD_STR_ARG(symbol), MD_STR_ARG(name));
    }
    const md_str_t ff_type = md_particle_type_ff_type(&sys.atom.type, (size_t)md_atom_particle_type_idx(&sys.atom, atom_idx));
    if (!md_str_empty(ff_type)) {
        tooltip_row(req, MD_STR_LIT("type"), MD_STR_FMT, MD_STR_ARG(ff_type));
    }
    if (particle != MD_PARTICLE_ATOM) {
        // An atom weighs what its element does, a bead what its force field says
        const float mass = md_atom_mass(&sys.atom, atom_idx);
        if (mass > 0.0f) tooltip_row(req, MD_STR_LIT("mass"), "%.3f Da", mass);
    }

    // Its roles, as particle(role=...) takes them: what the strand runs through (backbone), it and what of its residue
    // hangs off it that is no side (main), the side (an amino acid's side chain, a nucleotide's nucleoside) and a
    // nucleotide's finer base and sugar, its end of the strand, and the QM region. terminus, either end, is left out, as
    // the end it is says it.
    md_str_writer_t role = md_str_writer_create(req->alloc);
    if (flags & MD_ATOM_FLAG_BACKBONE)   tooltip_list_add(&role, "backbone");
    if (flags & MD_ATOM_FLAG_MAIN_CHAIN) tooltip_list_add(&role, "main");
    if (flags & (MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_NUCLEOSIDE)) tooltip_list_add(&role, "side");
    if (flags & MD_ATOM_FLAG_NUCLEOBASE) {
        tooltip_list_add(&role, "base");
    } else if (flags & MD_ATOM_FLAG_NUCLEOSIDE) {
        tooltip_list_add(&role, "sugar");
    }
    if (flags & MD_ATOM_FLAG_TERMINAL_BEG) tooltip_list_add(&role, nucleotide ? "5_terminus" : "n_terminus");
    if (flags & MD_ATOM_FLAG_TERMINAL_END) tooltip_list_add(&role, nucleotide ? "3_terminus" : "c_terminus");
    if (md_atom_is_qm(&sys.atom, atom_idx)) tooltip_list_add(&role, "qm");
    tooltip_row_list(req, MD_STR_LIT("role"), role);

    // The chemical piece of its molecule it is in
    tooltip_fragment_row(req, sys, atom_idx);

    // Its chemistry: a row per keyword of particle() which selects it, keyed as the keyword and with the value as a
    // script writes it, unquoted (hybridization='sp2', aromatic=True, hydrogens=1, formal_charge=-1). A row says what
    // the atom has, so none where it has nothing to say.
    const md_hybridization_t hyb = md_atom_flags_hybridization(flags);
    if (hyb != MD_HYBRIDIZATION_UNKNOWN) tooltip_row(req, MD_STR_LIT("hybridization"), "%s", md_hybridization_name(hyb));
    if (flags & MD_ATOM_FLAG_AROMATIC)   tooltip_row(req, MD_STR_LIT("aromatic"), "True");
    // Of a metallic phase: an atom of a metallic bond (md_bond_metallic_atoms), which only a metal is. A metal takes no
    // part in the perception, which gives it no hybridization, so only an atom without one is looked for in the bonds.
    if (hyb == MD_HYBRIDIZATION_UNKNOWN && z > 1) {
        for (size_t i = 0; i < sys.bond.count; ++i) {
            const md_atom_pair_t p = sys.bond.pairs[i];
            if (md_bond_kind(sys.bond.flags[i]) == MD_BOND_KIND_METALLIC && (p.idx[0] == (md_atom_idx_t)atom_idx || p.idx[1] == (md_atom_idx_t)atom_idx)) {
                tooltip_row(req, MD_STR_LIT("metallic"), "True");
                break;
            }
        }
    }
    if (num_h > 0)     tooltip_row(req, MD_STR_LIT("hydrogens"), "%d", num_h);
    if (formal_charge) tooltip_row(req, MD_STR_LIT("formal_charge"), "%+d", formal_charge);

    tooltip_atom_attribute_rows(req, sc, atom_idx);

    const md_float3_t pos = md_state_coord(&picking_hit_instance(sc, req->hit).state, atom_idx);
    tooltip_row(req, MD_STR_LIT("position"), "%.3f, %.3f, %.3f " TOOLTIP_ANGSTROM, pos.x, pos.y, pos.z);

    md_str_writer_t index = md_str_writer_create(req->alloc);
    script_write_position(&index, sc.app->script.language, "atom", atom_idx);
    tooltip_append_index(&index, sys, place, sc.app->script.language);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// The index of the system's bond of this pair and kind (md_bond_data_t), as VIAL's bond[i] numbers it; -1 when it stores
// none. The graph's and the interactions' are looked up; the elastic and metallic bonds, which the graph leaves out, are
// looked for along the first range.
static md_bond_idx_t tooltip_stored_bond(const md_system_t& sys, md_atom_pair_t pair, md_bond_flags_t flags) {
    const md_bond_data_t* bd = &sys.bond;
    const md_bond_kind_t kind = md_bond_kind(flags);
    md_bond_idx_t b = md_bond_find(bd, pair.idx[0], pair.idx[1]);
    if (b >= 0 && md_bond_kind(bd->flags[b]) == kind) return b;
    b = md_bond_find_interaction(bd, pair.idx[0], pair.idx[1]);
    if (b >= 0 && md_bond_kind(bd->flags[b]) == kind) return b;
    for (size_t i = 0; i < bd->count; ++i) {
        const md_atom_pair_t p = bd->pairs[i];
        const bool same = (p.idx[0] == pair.idx[0] && p.idx[1] == pair.idx[1]) || (p.idx[0] == pair.idx[1] && p.idx[1] == pair.idx[0]);
        if (same && md_bond_kind(bd->flags[i]) == kind) return (md_bond_idx_t)i;
    }
    return -1;
}

// A bond of the system, one formed in the frame shown, or two atoms a representation draws a bond between which none
// joins (picked_bond)
static void tooltip_bond(PickingTooltipTextRequest* req, const SystemContext& sc, md_atom_pair_t pair, md_bond_flags_t flags, PickedBondSource source) {
    const md_system_t& sys = sc.data.sys;
    const md_system_state_t& sys_state = picking_hit_instance(sc, req->hit).state;
    const int order = md_bond_order(flags);

    char symbol = '-';
    if (order == MD_BOND_ORDER_DOUBLE)    symbol = '=';
    if (order == MD_BOND_ORDER_TRIPLE)    symbol = '#';
    if (order == MD_BOND_ORDER_QUADRUPLE) symbol = '$';
    if (flags & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) symbol = ':';

    const TooltipPlace place[2] = {
        tooltip_place_of_atom(sys, pair.idx[0]),
        tooltip_place_of_atom(sys, pair.idx[1]),
    };

    // "bond · C1-C2 · LIG1"; two atoms a representation joins with no bond are a pair
    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += source == PickedBondSource::Unbonded ? MD_STR_LIT("pair" TOOLTIP_SEP) : MD_STR_LIT("bond" TOOLTIP_SEP);
    tooltip_append_atom_name(&title, sys, pair.idx[0]);
    title += symbol;
    tooltip_append_atom_name(&title, sys, pair.idx[1]);
    tooltip_append_place(&title, sys, place[0], &place[1]);
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    // What it is, as bond() selects it: its kind, its order and, below, its origin, each as the keyword writes it
    if (source == PickedBondSource::Unbonded) {
        tooltip_row(req, MD_STR_LIT("kind"), "no bond (a pair of the representation)");
    } else {
        tooltip_row(req, MD_STR_LIT("kind"), "%s", bond_kind_str(flags));
        // A resonance's word, else 1 to 4; none where it is unknown, as a coordinate bond's
        md_str_writer_t ord = md_str_writer_create(req->alloc);
        if (flags & MD_BOND_FLAG_AROMATIC) {
            ord += MD_STR_LIT("aromatic");
        } else if (flags & MD_BOND_FLAG_DELOCALIZED) {
            ord += MD_STR_LIT("delocalized");
        } else if (order > 0) {
            md_str_writer_fmt(&ord, "%d", order);
        }
        if (!md_str_empty(md_str_writer_to_str(&ord))) {
            if (flags & MD_BOND_FLAG_ORDER_PERCEIVED) ord += MD_STR_LIT(" (perceived)");
            tooltip_row_list(req, MD_STR_LIT("order"), ord);
        }
    }

    // The fragment it is within, or the two it joins: "amine / phenyl"
    if (source == PickedBondSource::System) {
        const md_fragment_idx_t f0 = md_system_fragment_find_by_atom_idx(&sys, pair.idx[0]);
        const md_fragment_idx_t f1 = md_system_fragment_find_by_atom_idx(&sys, pair.idx[1]);
        if (f0 != MD_FRAGMENT_NONE && f0 == f1) {
            tooltip_fragment_row(req, sys, pair.idx[0]);
        } else if (f0 != MD_FRAGMENT_NONE && f1 != MD_FRAGMENT_NONE) {
            md_str_writer_t frags = md_str_writer_create(req->alloc);
            tooltip_append_fragment(&frags, sys, f0);
            frags += " / ";
            tooltip_append_fragment(&frags, sys, f1);
            tooltip_row_list(req, MD_STR_LIT("fragments"), frags);
        }
    }

    // Across a periodic boundary the bond is the short way round
    md_float3_t d = md_float3_sub(md_state_coord(&sys_state, pair.idx[1]), md_state_coord(&sys_state, pair.idx[0]));
    md_geom_min_image(d.elem, 0, 1, &sys_state.unitcell);
    tooltip_row(req, source == PickedBondSource::Unbonded ? MD_STR_LIT("distance") : MD_STR_LIT("length"), "%.3f " TOOLTIP_ANGSTROM, md_float3_length(d));

    if (source == PickedBondSource::Formed) {
        // Not a bond of the system: tracked, of this frame alone, inferred from its distances
        tooltip_row(req, MD_STR_LIT("origin"), "%s (formed in this frame)", md_bond_origin_name(md_bond_origin(flags)));
    } else if (source == PickedBondSource::Found) {
        // Not a bond of the system: what the filter found for this frame, a hydrogen, ionic or halogen bond or a bond formed
        tooltip_row(req, MD_STR_LIT("origin"), "%s (found in this frame)", md_bond_origin_name(md_bond_origin(flags)));
    } else if (source == PickedBondSource::System) {
        tooltip_row(req, MD_STR_LIT("origin"), "%s", md_bond_origin_name(md_bond_origin(flags)));
    }

    // A bond the system stores has its number in VIAL, bond[i], beside its particles; one found or formed in a frame is
    // stored nowhere, and is named by its particles alone, which is also what stays the same one whatever is inferred
    md_str_writer_t index = md_str_writer_create(req->alloc);
    if (source == PickedBondSource::System && sc.app->script.language == script::Language::V2) {
        const md_bond_idx_t b = tooltip_stored_bond(sys, pair, flags);
        if (b >= 0) md_str_writer_fmt(&index, "bond[%lld]", (long long)b);
    }
    tooltip_index_pair(&index, sc.app->script.language, "atom", pair.idx[0], pair.idx[1]);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// A backbone segment stands for its residue: the cartoon and ribbons have no atoms to speak of. What the segment adds is
// its secondary structure, as residue(ss=) selects it: DSSP of the frame shown (the state's), and the state the cartoon
// draws where that differs, as the cartoon smooths the states over frames; and its dihedrals, as VIAL's
// backbone_angles() gives them of the coordinates shown (md_backbone_segment_dihedrals): phi needs a neighbour before
// it in the chain, psi and omega one after it, and a trace has none.
static void tooltip_backbone_segment_geometry(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t seg_idx) {
    const md_system_t& sys = sc.data.sys;
    const md_protein_backbone_data_t& bb = sys.protein_backbone;

    // Indexed by md_secondary_structure_t, as residue(ss=) spells them; unknown has none
    static const char* ss_words[] = { "", "coil", "turn", "bend", "helix_310", "helix_alpha", "helix_pi", "sheet", "bridge" };
    auto ss_word = [](const md_secondary_structure_t* ss, uint32_t i) -> const char* {
        return ss && (size_t)ss[i] < MD_COUNTOF(ss_words) ? ss_words[ss[i]] : "";
    };
    const SystemInstance& inst = picking_hit_instance(sc, req->hit);
    const md_secondary_structure_t* drawn = displayed_secondary_structure(inst);
    const md_secondary_structure_t* frame = md_system_state_secondary_structure(&inst.state, &sys);
    const char* ss = frame ? ss_word(frame, seg_idx) : ss_word(drawn, seg_idx);
    const char* ss_drawn = ss_word(drawn, seg_idx);
    if (ss[0] && ss_drawn[0] && strcmp(ss, ss_drawn) != 0) {
        tooltip_row(req, MD_STR_LIT("ss"), "%s (drawn as %s)", ss, ss_drawn);
    } else if (ss[0]) {
        tooltip_row(req, MD_STR_LIT("ss"), "%s", ss);
    }

    const md_system_state_t& sys_state = inst.state;
    float a[3];
    const uint32_t has = sys_state.xyz ? md_backbone_segment_dihedrals(a, seg_idx, sys_state.xyz, &sys_state.unitcell, &bb) : 0;
    if (has) {
        // phi, psi, omega and the degree sign in UTF-8 escapes
        static const char* names[3] = { "\xCF\x86", "\xCF\x88", "\xCF\x89" };
        md_str_writer_t angles = md_str_writer_create(req->alloc);
        for (int c = 0; c < 3; ++c) {
            if (!(has & (1u << c))) continue;
            if (!md_str_writer_empty(&angles)) angles += ", ";
            md_str_writer_fmt(&angles, "%s %.1f\xC2\xB0", names[c], MD_RAD_TO_DEG(a[c]));
        }
        tooltip_row(req, MD_STR_LIT("backbone_angles"), "%s", md_str_writer_to_cstr(&angles));
    }
}

// What a residue is, a row per keyword of residue() which selects it, keyed as the keyword and with the value as a script
// writes it, unquoted: its kind (kind='amino_acid'), its ends (role='n_terminus'; terminus, either end, is left out, as
// the end it is says it), a break of its chain beside it (gap='before') and, of an amino acid or a nucleotide whose
// backbone was not found, resolved=False
static void tooltip_residue_rows(PickingTooltipTextRequest* req, const md_system_t& sys, int res_idx) {
    const md_residue_kind_t  kind  = md_system_residue_kind(&sys, res_idx);
    const md_residue_flags_t flags = md_system_residue_flags(&sys, res_idx);
    const bool nucleotide = kind == MD_RESIDUE_KIND_NUCLEOTIDE;

    // Indexed by md_residue_kind_t, as residue(kind=) spells them
    static const char* kinds[] = { "other", "amino_acid", "nucleotide", "water", "ion", "lipid" };
    if ((size_t)kind < MD_COUNTOF(kinds)) tooltip_row(req, MD_STR_LIT("kind"), "%s", kinds[kind]);

    md_str_writer_t role = md_str_writer_create(req->alloc);
    if (flags & MD_RESIDUE_FLAG_TERMINAL_BEG) tooltip_list_add(&role, nucleotide ? "5_terminus" : "n_terminus");
    if (flags & MD_RESIDUE_FLAG_TERMINAL_END) tooltip_list_add(&role, nucleotide ? "3_terminus" : "c_terminus");
    tooltip_row_list(req, MD_STR_LIT("role"), role);

    md_str_writer_t gap = md_str_writer_create(req->alloc);
    if (flags & MD_RESIDUE_FLAG_BREAK_BEFORE) tooltip_list_add(&gap, "before");
    if (flags & MD_RESIDUE_FLAG_BREAK_AFTER)  tooltip_list_add(&gap, "after");
    tooltip_row_list(req, MD_STR_LIT("gap"), gap);

    if ((kind == MD_RESIDUE_KIND_AMINO_ACID || nucleotide) && !(flags & MD_RESIDUE_FLAG_RESOLVED)) {
        tooltip_row(req, MD_STR_LIT("resolved"), "False");
    }
}

static void tooltip_backbone_segment(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t seg_idx) {
    const md_system_t& sys = sc.data.sys;
    const md_protein_backbone_data_t& bb = sys.protein_backbone;
    MD_ASSERT(seg_idx < bb.segment.count && bb.segment.res_idx);

    const int res_idx = bb.segment.res_idx[seg_idx];
    if (res_idx < 0 || (size_t)res_idx >= sys.residue.count) {
        // No residue (a backbone found from the bonds alone, as for an XYZ file): no residue to name, only what the
        // segment adds, and the atom of its CA
        tooltip_title(req, "amino acid");
        tooltip_backbone_segment_geometry(req, sc, seg_idx);
        md_str_writer_t index = md_str_writer_create(req->alloc);
        script_write_position(&index, sc.app->script.language, "atom", bb.segment.atoms[seg_idx].ca);
        tooltip_row_list(req, MD_STR_LIT("index"), index);
        return;
    }
    const md_urange_t range = md_system_residue_atom_range(&sys, res_idx);
    const TooltipPlace place = tooltip_place_of_residue(sys, res_idx);

    // "residue · ALA5 · A": the segment stands for its residue
    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += MD_STR_LIT("residue");
    tooltip_append_place(&title, sys, place);
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    tooltip_residue_rows(req, sys, res_idx);

    tooltip_backbone_segment_geometry(req, sc, seg_idx);

    tooltip_row(req, MD_STR_LIT("atoms"), "%u", range.end - range.beg);

    tooltip_part_and_molecule_rows(req, sys, place);

    md_str_writer_t index = md_str_writer_create(req->alloc);
    tooltip_append_index(&index, sys, place, sc.app->script.language);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// The gap of a chain over a break after the residue res_idx (md_gl_mol_set_chain_gaps): the residues on either side,
// how many residues are missing between them by their numbering, and how far apart the atoms a link would join are in
// the current frame. The numbering is that of the residues of the sequence in an mmCIF file, which counts the missing
// ones exactly, but the author's in a PDB file, which may skip where nothing is missing; a model built in pieces may
// number on over a break where something is.
static void tooltip_chain_gap(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t res_idx) {
    const md_system_t& sys = sc.data.sys;
    MD_ASSERT((size_t)res_idx + 1 < sys.residue.count);
    if (!(md_system_residue_flags(&sys, res_idx) & MD_RESIDUE_FLAG_BREAK_AFTER)) return;

    const md_urange_t range = md_system_residue_atom_range(&sys, res_idx);
    const TooltipPlace place[2] = {
        tooltip_place_of_residue(sys, (int)res_idx),
        tooltip_place_of_residue(sys, (int)res_idx + 1),
    };

    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += "chain break";
    tooltip_append_place(&title, sys, place[0], &place[1]);
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    const int missing = md_system_residue_seq_id(&sys, res_idx + 1) - md_system_residue_seq_id(&sys, res_idx) - 1;
    if (missing > 0) {
        tooltip_row(req, MD_STR_LIT("missing"), "%d residue%s by the numbering", missing, missing == 1 ? "" : "s");
    } else {
        tooltip_row(req, MD_STR_LIT("missing"), "none by the numbering");
    }

    md_atom_idx_t atoms[2];
    if (md_system_chain_break_atoms(atoms, &sys, res_idx)) {
        const md_system_state_t& sys_state = picking_hit_instance(sc, req->hit).state;
        md_float3_t d = md_float3_sub(md_state_coord(&sys_state, atoms[1]), md_state_coord(&sys_state, atoms[0]));
        md_geom_min_image(d.elem, 0, 1, &sys_state.unitcell);
        md_str_writer_t span = md_str_writer_create(req->alloc);
        md_str_writer_fmt(&span, "%.1f " TOOLTIP_ANGSTROM " (", md_float3_length(d));
        tooltip_append_atom_name(&span, sys, atoms[0]);
        span += " - ";
        tooltip_append_atom_name(&span, sys, atoms[1]);
        span += ")";
        tooltip_row(req, MD_STR_LIT("span"), "%s", md_str_writer_to_cstr(&span));
    }

    md_str_writer_t index = md_str_writer_create(req->alloc);
    tooltip_index_pair(&index, sc.app->script.language, "residue", res_idx, (int64_t)res_idx + 1);
    const md_chain_idx_t chain = md_system_chain_find_by_residue_idx(&sys, res_idx);
    const int strand = tooltip_polymer_position(sys, place[0].part_idx);
    const md_molecule_idx_t mol = range.beg < range.end ? md_system_molecule_find_by_atom_idx(&sys, range.beg) : -1;
    if (chain != -1)       script_write_position(&index, sc.app->script.language, "chain", chain);
    else if (strand != -1) script_write_position(&index, sc.app->script.language, "polymer", strand);
    else if (mol != -1)    script_write_position(&index, sc.app->script.language, "molecule", mol);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// A ring (the circle of an aromatic ring): what it is, its atoms around it and their elements, how flat it is in the
// current frame, the RMS distance of its atoms from their plane (the plane the circle is drawn in, ring.vert), and its
// pucker
static void tooltip_ring(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t ring_idx) {
    const md_system_t& sys = sc.data.sys;
    const md_system_state_t& sys_state = picking_hit_instance(sc, req->hit).state;
    // A ring of the system, or one of the frame where its bonds changed (picked_ring)
    const int32_t* atoms = nullptr;
    size_t count = 0;
    md_ring_flags_t flags = MD_RING_FLAG_NONE;
    if (!picked_ring(sc, ring_idx, &atoms, &count, &flags) || count == 0) return;
    for (size_t i = 0; i < count; ++i) {
        if (atoms[i] < 0 || (size_t)atoms[i] >= sys.atom.count) return;
    }
    const bool frame_ring = ring_idx >= md_ring_data_count(&sys.ring);

    // Where it is: the place of its first atom, and of another where the ring reaches over (a ring is mostly in one)
    const TooltipPlace place = tooltip_place_of_atom(sys, (size_t)atoms[0]);
    TooltipPlace other = place;
    for (size_t i = 1; i < count; ++i) {
        const TooltipPlace p = tooltip_place_of_atom(sys, (size_t)atoms[i]);
        if (p.res_idx != place.res_idx || p.part_idx != place.part_idx) {
            other = p;
            break;
        }
    }

    // "ring", whatever it is: the kind row below says aromatic, as ring(kind=) selects it
    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += MD_STR_LIT("ring");
    tooltip_append_place(&title, sys, place, &other);
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    // Its kinds, each as ring(kind=...) selects it: every tag it has. fused_aromatic is a ring of a fused aromatic system
    // (its circle is dashed), whose sextets are not assigned to particular rings.
    md_str_writer_t kind = md_str_writer_create(req->alloc);
    if (flags & MD_RING_FLAG_AROMATIC) tooltip_list_add(&kind, "aromatic");
    if (flags & MD_RING_FLAG_HETERO)   tooltip_list_add(&kind, "hetero");
    if (flags & MD_RING_FLAG_FUSED)    tooltip_list_add(&kind, "fused");
    if ((flags & MD_RING_FLAG_AROMATIC) && (flags & MD_RING_FLAG_FUSED_AROMATIC)) tooltip_list_add(&kind, "fused_aromatic");
    if (flags & MD_RING_FLAG_SUGAR)    tooltip_list_add(&kind, "sugar");
    tooltip_row_list(req, MD_STR_LIT("kind"), kind);
    tooltip_row(req, MD_STR_LIT("size"), "%zu", count);

    // The ring system it is of, as a fragment: "indolyl"
    if (!frame_ring) {
        tooltip_fragment_row(req, sys, (size_t)atoms[0]);
    }

    // The atoms in order around the ring, then its composition: C first, the other elements by atomic number
    md_str_writer_t names = md_str_writer_create(req->alloc);
    int element_count[128] = {0};
    for (size_t i = 0; i < count; ++i) {
        if (i > 0) names += ' ';
        tooltip_append_atom_name(&names, sys, (size_t)atoms[i]);
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, (size_t)atoms[i]);
        if (0 < z && z < (int)MD_COUNTOF(element_count)) element_count[z] += 1;
    }
    md_str_writer_t formula = md_str_writer_create(req->alloc);
    for (int pass = 0; pass < 2; ++pass) {
        for (int z = 1; z < (int)MD_COUNTOF(element_count); ++z) {
            if (!element_count[z] || (pass == 0) != (z == 6)) continue;
            formula += md_atomic_number_symbol((md_atomic_number_t)z);
            if (element_count[z] > 1) md_str_writer_fmt(&formula, "%d", element_count[z]);
        }
    }
    if (!md_str_writer_empty(&formula)) {
        tooltip_row(req, MD_STR_LIT("atoms"), "%s (%s)", md_str_writer_to_cstr(&names), md_str_writer_to_cstr(&formula));
    } else {
        tooltip_row(req, MD_STR_LIT("atoms"), "%s", md_str_writer_to_cstr(&names));
    }

    // Planarity, as planarity() of VIAL has it: the RMS distance of the atoms from their plane of least spread (the
    // eigenvector of the least eigenvalue of their scatter), the atoms the short way round from the first across a
    // periodic boundary
    md_float3_t p[MD_RING_MAX_SIZE];
    if (count <= MD_COUNTOF(p)) {
        const md_float3_t p0 = md_state_coord(&sys_state, (size_t)atoms[0]);
        md_float3_t c = {0, 0, 0};
        for (size_t i = 0; i < count; ++i) {
            md_float3_t d = md_float3_sub(md_state_coord(&sys_state, (size_t)atoms[i]), p0);
            md_geom_min_image(d.elem, 0, 1, &sys_state.unitcell);
            p[i] = md_float3_add(p0, d);
            c = md_float3_add(c, p[i]);
        }
        c = md_float3_mul1(c, 1.0f / (float)count);
        const md_float3x3_eigen_t eig = md_float3x3_eigen(md_float3x3_covariance_matrix(p, nullptr, nullptr, count, c));
        md_float3_t n = eig.vectors.col[2];     // of the least eigenvalue: they are in descending order
        const float len = md_float3_length(n);
        if (len > 0.0f) {
            n = md_float3_mul1(n, 1.0f / len);
            float sum2 = 0.0f;
            for (size_t i = 0; i < count; ++i) {
                const float h = md_float3_dot(md_float3_sub(p[i], c), n);
                sum2 += h * h;
            }
            tooltip_row(req, MD_STR_LIT("planarity"), "%.3f " TOOLTIP_ANGSTROM, sqrtf(sum2 / (float)count));
        }
    }

    // Its pucker, as VIAL's ring(pucker=) selects it and pucker() and pseudorotation() give it (md_pucker_ring): the kind
    // of its conformation, Q, and theta and phi where the ring's order is fixed by its chemistry, theta of a ring of six;
    // the pseudorotation of a furanose. A ring of the frame is none of the system's, which VIAL knows
    md_pucker_ring_t pr;
    if (!frame_ring && md_pucker_ring(&pr, &sys, &sys_state, ring_idx)) {
        // Indexed by md_pucker_kind_t, as ring(pucker=) spells them
        static const char* kinds[] = { "", "planar", "chair", "boat", "skew", "half_chair", "envelope", "twist" };
        const md_pucker_kind_t kind = md_pucker_conformer(&pr.cp).kind;
        md_str_writer_t pucker = md_str_writer_create(req->alloc);
        if ((size_t)kind < MD_COUNTOF(kinds) && kinds[kind][0]) md_str_writer_fmt(&pucker, "%s, ", kinds[kind]);
        md_str_writer_fmt(&pucker, "Q %.2f " TOOLTIP_ANGSTROM, pr.cp.amplitude);
        if (pr.anchored && pr.size == 6) md_str_writer_fmt(&pucker, ", \xCE\xB8 %.1f\xC2\xB0", pr.cp.theta);
        if (pr.anchored) md_str_writer_fmt(&pucker, ", \xCF\x86 %.1f\xC2\xB0", pr.cp.phi2);
        tooltip_row(req, MD_STR_LIT("pucker"), "%s", md_str_writer_to_cstr(&pucker));
        if (pr.has_pseudo) {
            tooltip_row(req, MD_STR_LIT("pseudorotation"), "P %.1f\xC2\xB0, \xCF\x84m %.1f\xC2\xB0", pr.pseudo_phase, pr.pseudo_amplitude);
        }
    }

    // ring(k) of v1, ring[k - 1] of VIAL, is ring k - 1 of the system, as the circle picks it. A ring of the frame is
    // none of the system's: the bonds formed in the frame close it, and the script knows the rings of the system
    md_str_writer_t index = md_str_writer_create(req->alloc);
    if (frame_ring) {
        tooltip_list_add(&index, "of this frame, its bonds as they are tracked");
    } else {
        script_write_position(&index, sc.app->script.language, "ring", ring_idx);
    }
    tooltip_append_index(&index, sys, place, sc.app->script.language);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// A point of the cartoon of the lipids and nucleic acids stands for its residue: the residue, and what part of it
// the point is
static void tooltip_cartoon_glyph(PickingTooltipTextRequest* req, const SystemContext& sc, uint32_t point_idx) {
    const md_system_t& sys = sc.data.sys;
    const CartoonGlyphs& glyphs = system_instance(sc).cartoon_glyphs;
    // Indexed by CartoonGlyphKind
    static const char* parts[] = { "head", "backbone", "tail", "sterol skeleton", "backbone", "sugar", "base" };
    static_assert(MD_COUNTOF(parts) == CartoonGlyph_Count, "A name for every kind of glyph");

    const int32_t res_idx = cartoon_glyphs_point_residue(glyphs, point_idx);
    if (res_idx < 0 || (size_t)res_idx >= sys.residue.count) {
        // No residue (a nucleic backbone found from the bonds alone, as for an XYZ file): only what part it is
        if (glyphs.point_kind[point_idx] < MD_COUNTOF(parts)) tooltip_title(req, "nucleotide" TOOLTIP_SEP "%s", parts[glyphs.point_kind[point_idx]]);
        return;
    }
    const md_urange_t range = md_system_residue_atom_range(&sys, res_idx);
    const TooltipPlace place = tooltip_place_of_residue(sys, res_idx);

    // The title says what part of the residue the point is, which no keyword selects: "residue · DPPC 12 · head"
    md_str_writer_t title = md_str_writer_create(req->alloc);
    title += MD_STR_LIT("residue");
    tooltip_append_place(&title, sys, place);
    if (glyphs.point_kind[point_idx] < MD_COUNTOF(parts)) {
        title += TOOLTIP_SEP;
        title += md_str_from_cstr(parts[glyphs.point_kind[point_idx]]);
    }
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    tooltip_residue_rows(req, sys, res_idx);

    tooltip_row(req, MD_STR_LIT("atoms"), "%u", range.end - range.beg);

    tooltip_part_and_molecule_rows(req, sys, place);

    md_str_writer_t index = md_str_writer_create(req->alloc);
    tooltip_append_index(&index, sys, place, sc.app->script.language);
    tooltip_row_list(req, MD_STR_LIT("index"), index);
}

// The hit names the group and the element outright - its range was reserved against that group's attribute - so there
// is no list to rebuild and nothing to look the index up in.
static void tooltip_dipole(PickingTooltipTextRequest* req, const SystemContext& sc, const PickingHit& hit) {
    const md_system_t& sys = sc.data.sys;
    DipoleGroup group = {};
    if (!dipole_group_from_key(&group, sys, hit.key) || hit.local_idx >= group.count) return;

    char label[64];
    const int label_len = dipole_entry_label(label, sizeof(label), group, hit.local_idx);
    tooltip_title(req, "dipole" TOOLTIP_SEP "%.*s", label_len, label);

    md_float3_t vec = {0, 0, 0};
    if (dipole_moment_read(&vec, nullptr, sys, group.key, hit.local_idx)) {
        char unit[32];
        if (md_unit_is_atomic(group.unit)) {
            snprintf(unit, sizeof(unit), " a.u.");
        } else {
            tooltip_unit_suffix(unit, sizeof(unit), group.unit);
        }
        tooltip_row(req, MD_STR_LIT("magnitude"), "%.3f%s", md_float3_length(vec), unit);
        tooltip_row(req, MD_STR_LIT("vector"), "%.3f, %.3f, %.3f%s", vec.x, vec.y, vec.z, unit);
    }
}

// A point of a point set drawn as glyphs: what it is, every quantity of its set at it
static void tooltip_glyph(PickingTooltipTextRequest* req, const SystemContext& sc, const PickingHit& hit) {
    const Representation* rep = nullptr;
    for (size_t i = 0; i < md_array_size(sc.representation.reps); ++i) {
        if (sc.representation.reps[i].node.id == (uint32_t)hit.key) rep = &sc.representation.reps[i];
    }
    if (!rep) return;
    const GlyphRepresentation& g = rep->glyphs;
    const md_system_t& sys = sc.data.sys;
    if (hit.local_idx >= glyphs_count(g)) return;

    // titled by its group's leaf: "tessera 12"
    const md_str_t group = glyphs_group(sys, g.points);
    const size_t slash = md_str_rfind_char(group, '/');
    const md_str_t name = slash == MD_STR_NPOS ? group : md_str_substr(group, slash + 1, group.len);
    md_str_writer_t title = md_str_writer_create(req->alloc);
    md_writef(&title.w, MD_STR_FMT " %u", MD_STR_ARG(name), hit.local_idx + 1);
    const uint32_t atom = glyphs_point_atom(sys, g.points, hit.local_idx);
    if (atom != UINT32_MAX) {
        title += TOOLTIP_SEP;
        tooltip_append_atom_name(&title, sys, atom);
    }
    tooltip_title(req, "%s", md_str_writer_to_cstr(&title));

    md_attribute_id_t q[64];
    const size_t nq = MD_MIN(glyphs_quantities(q, MD_COUNTOF(q), sys, g.points), MD_COUNTOF(q));
    for (size_t i = 0; i < nq; ++i) {
        const md_attribute_t* a = md_attributes_get(&sys.attributes, q[i]);
        const md_str_t leaf = md_attribute_leaf(a);
        if (md_str_eq(leaf, MD_STR_LIT("atom_index")) || md_str_eq(leaf, MD_STR_LIT("system_index"))) continue;
        const double v = glyphs_value(sys, q[i], hit.local_idx);
        if (!isfinite(v)) continue;
        char unit[64];
        double f = 1.0;
        glyphs_unit_print(unit, sizeof(unit), &f, a->unit);
        const md_str_t label = a->label.len ? a->label : leaf;
        tooltip_row(req, label, unit[0] ? "%.4g %s" : "%.4g", v * f, unit);
    }
    if (atom != UINT32_MAX) {
        md_str_writer_t index = md_str_writer_create(req->alloc);
        script_write_position(&index, sc.app->script.language, "atom", atom);
        tooltip_row_list(req, MD_STR_LIT("index"), index);
    }
}

static void fill_picking_tooltip(PickingTooltipTextRequest* req, const SystemContext& sc, const PickingHit& hit) {
    MD_ASSERT(req);
    const md_system_t& sys = sc.data.sys;

    if (hit.domain == PickingDomain_Atom) {
        if (hit.local_idx < sys.atom.count) tooltip_atom(req, sc, hit.local_idx);
    } else if (hit.domain == PickingDomain_Bond) {
        md_atom_pair_t pair;
        md_bond_flags_t flags;
        PickedBondSource source;
        if (picked_bond(*sc.app, hit, &pair, &flags, &source)) tooltip_bond(req, sc, pair, flags, source);
    } else if (hit.domain == PickingDomain_BackboneSegment) {
        if (hit.local_idx < sys.protein_backbone.segment.count && sys.protein_backbone.segment.res_idx) tooltip_backbone_segment(req, sc, hit.local_idx);
    } else if (hit.domain == PickingDomain_Ring) {
        tooltip_ring(req, sc, hit.local_idx);
    } else if (hit.domain == PickingDomain_CartoonGlyph) {
        const uint32_t point_idx = cartoon_glyphs_picked_point(system_instance(sc).cartoon_glyphs, hit.local_idx);
        if (point_idx < system_instance(sc).cartoon_glyphs.num_points) tooltip_cartoon_glyph(req, sc, point_idx);
    } else if (hit.domain == PickingDomain_ChainGap) {
        if ((size_t)hit.local_idx + 1 < sys.residue.count) tooltip_chain_gap(req, sc, hit.local_idx);
    } else if (hit.domain == PickingDomain_Dipole) {
        tooltip_dipole(req, sc, hit);
    } else if (hit.domain == PickingDomain_Glyph) {
        tooltip_glyph(req, sc, hit);
    }
}

static void tooltip_text(md_str_t str) {
    if (md_str_empty(str)) {
        ImGui::TextUnformatted("");
    } else {
        ImGui::TextUnformatted(str.ptr, str.ptr + str.len);
    }
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state) {
    if (hit.raw_idx == INVALID_PICKING_IDX) return;
    
	md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
    MD_DEFER { md_temp_end(temp); };

    PickingTooltipTextRequest tooltip_request = {
        .app = state,
        .hit = hit,
        .alloc = state.allocator.frame,
    };

    viamd::event_system_broadcast_event(viamd::EventType_ViamdPickingTooltipTextRequest, viamd::EventPayloadType_PickingTooltipTextRequest, &tooltip_request);

    const size_t num_lines = md_array_size(tooltip_request.lines);
    if (num_lines == 0) return;

    const ImVec2 offset = { 10.f, 18.f };
    const ImVec2 new_pos = {ImGui::GetMousePos().x + offset.x, ImGui::GetMousePos().y + offset.y};
    ImGui::SetNextWindowPos(new_pos);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Picking Tooltip Window", 0,
        ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);

    // A title begins a section, and the rows under it are one table so that their keys and values line up
    bool table_open = false;
    int  num_tables = 0;
    for (size_t i = 0; i < num_lines; ++i) {
        const TooltipLine& line = tooltip_request.lines[i];
        if (line.kind == TooltipLineKind::Title) {
            if (table_open) {
                ImGui::EndTable();
                table_open = false;
            }
            if (i > 0) ImGui::Separator();
            tooltip_text(line.text);
        } else {
            if (!table_open) {
                char id[16];
                snprintf(id, sizeof(id), "##rows%d", num_tables++);
                table_open = ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingFixedFit);
                if (!table_open) continue;
            }
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
            tooltip_text(line.key);
            ImGui::PopStyleColor();
            ImGui::TableSetColumnIndex(1);
            tooltip_text(line.text);
        }
    }
    if (table_open) ImGui::EndTable();

    ImGui::End();
    ImGui::PopStyleColor();
}

void interrupt_async_tasks(ApplicationState* state) {
    task_system::cancel_all();

    script_evals_interrupt(state);

    task_system::wait_all();
}

static inline void clear_frame_cache(FrameCache* cache) {
#if FRAME_CACHE_SIZE == 4
    md_mm_storeu_epi32(cache->frame_idx, md_mm_set1_epi32(-1));
    md_lru_cache4_init(&cache->lru);
#elif FRAME_CACHE_SIZE == 8
    md_mm256_storeu_epi32(cache->frame_idx, md_mm256_set1_epi32(-1));
    md_lru_cache8_init(&cache->lru);
#endif
}

static inline void init_frame_cache(FrameCache* cache, size_t num_atoms, md_allocator_i* alloc) {
    clear_frame_cache(cache);
    size_t capacity = MD_ALIGN_TO(num_atoms, 16);
    for (size_t i = 0; i < FRAME_CACHE_SIZE; ++i) {
        md_array_resize(cache->states[i].xyz, capacity, alloc);
    }
}

static inline void free_frame_cache(FrameCache* cache, md_allocator_i* alloc) {
    for (size_t i = 0; i < FRAME_CACHE_SIZE; ++i) {
        md_array_free(cache->states[i].xyz, alloc);
    }
    clear_frame_cache(cache);
}

static inline bool find_frame_in_cache(int* out_slot_idx, int64_t frame_idx, const FrameCache* cache) {
    MD_ASSERT(out_slot_idx);
#if FRAME_CACHE_SIZE == 4
    md_128i frame_indices = md_mm_loadu_epi32(cache->frame_idx);
    md_128i cmp_mask = md_mm_cmpeq_epi32(frame_indices, md_mm_set1_epi32((int32_t)frame_idx));
    int mask = md_mm_movemask_epi8(cmp_mask);
    *out_slot_idx = md_ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
    return mask != 0;
#elif FRAME_CACHE_SIZE == 8
    md_256i frame_indices = md_mm256_loadu_epi32(cache->frame_idx);
    md_256i cmp_mask = md_mm256_cmpeq_epi32(frame_indices, md_mm256_set1_epi32((int32_t)frame_idx));
    int mask = md_mm256_movemask_epi8(cmp_mask);
    *out_slot_idx = md_ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
    return mask != 0;
#endif
}

static inline int find_lru_cache_slot(const FrameCache* cache) {
#if FRAME_CACHE_SIZE == 4
    return md_lru_cache4_get_lru(cache->lru);
#elif FRAME_CACHE_SIZE == 8
    return md_lru_cache8_get_lru(cache->lru);
#endif
}

static inline void set_mru_cache_slot(FrameCache* cache, int slot_idx) {
#if FRAME_CACHE_SIZE == 4
    md_lru_cache4_set_mru(&cache->lru, slot_idx);
#elif FRAME_CACHE_SIZE == 8
    md_lru_cache8_set_mru(&cache->lru, slot_idx);
#endif
}

void clear_system_frame_cache(SystemContext& sys) {
    for (size_t i = 0; i < sys.num_instances; ++i) clear_frame_cache(&sys.instances[i]->frame_cache);
}

enum class SecondaryStructureRenderClass : uint8_t {
    Coil,
    Helix,
    Sheet,
};

static SecondaryStructureRenderClass secondary_structure_render_class(md_secondary_structure_t ss) {
    switch (ss) {
    case MD_SECONDARY_STRUCTURE_HELIX_310:
    case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
    case MD_SECONDARY_STRUCTURE_HELIX_PI:
        return SecondaryStructureRenderClass::Helix;
    case MD_SECONDARY_STRUCTURE_BETA_SHEET:
    case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:
        return SecondaryStructureRenderClass::Sheet;
    case MD_SECONDARY_STRUCTURE_COIL:
    case MD_SECONDARY_STRUCTURE_TURN:
    case MD_SECONDARY_STRUCTURE_BEND:
    case MD_SECONDARY_STRUCTURE_UNKNOWN:
    default:
        return SecondaryStructureRenderClass::Coil;
    }
}

static md_secondary_structure_t secondary_structure_from_render_class(SecondaryStructureRenderClass cls) {
    switch (cls) {
    case SecondaryStructureRenderClass::Helix:
        return MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
    case SecondaryStructureRenderClass::Sheet:
        return MD_SECONDARY_STRUCTURE_BETA_SHEET;
    case SecondaryStructureRenderClass::Coil:
    default:
        return MD_SECONDARY_STRUCTURE_COIL;
    }
}

// Fills an isolated coil between two segments of the same structure in with that structure, within each chain, to
// keep the cartoon from flickering at the noise of the assignment. Presentation only.
static void secondary_structure_weights_fill_isolated_coils(md_gl_secondary_structure_t* weights, const md_protein_backbone_data_t* backbone) {
    auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
        return a.helix == b.helix && a.sheet == b.sheet;
    };
    const md_gl_secondary_structure_t ss_coil  = { 0, 0 };
    const md_gl_secondary_structure_t ss_helix = { .helix = 1.0f };
    const md_gl_secondary_structure_t ss_sheet = { .sheet = 1.0f };
    for (size_t r = 0; r < backbone->range.count; ++r) {
        for (size_t j = backbone->range.offset[r] + 1; j + 1 < backbone->range.offset[r + 1]; ++j) {
            if (!is_eq(weights[j], ss_coil)) continue;
            if (is_eq(weights[j - 1], ss_helix) && is_eq(weights[j + 1], ss_helix)) weights[j] = ss_helix;
            if (is_eq(weights[j - 1], ss_sheet) && is_eq(weights[j + 1], ss_sheet)) weights[j] = ss_sheet;
        }
    }
}

// The cartoon's secondary structure for the displayed state's own labels, uploaded directly: the weights are renderer
// input and are kept nowhere but on the GPU. Without a run there is nothing to blend between; interpolate_system_state
// uploads its own blend of the frames around the displayed one.
static void upload_secondary_structure_weights(SystemInstance& inst) {
    const md_system_t& sys = inst.sys->data.sys;
    const size_t num_segments = sys.protein_backbone.segment.count;
    const md_secondary_structure_t* ss = md_system_state_secondary_structure(&inst.state, &sys);
    if (!ss || num_segments == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_gl_secondary_structure_t* weights = md_temp_alloc_array(temp, md_gl_secondary_structure_t, num_segments);
    for (size_t i = 0; i < num_segments; ++i) {
        weights[i] = md_gl_secondary_structure_convert(ss[i]);
    }
    secondary_structure_weights_fill_isolated_coils(weights, &sys.protein_backbone);
    md_gl_mol_set_backbone_secondary_structure(inst.gl_mol, 0, (uint32_t)num_segments, weights, 0);
}

static void secondary_structure_render_denoise(md_secondary_structure_t* dst, const md_secondary_structure_t* src, size_t num_frames, size_t stride) {
    MD_ASSERT(dst);
    MD_ASSERT(src);

    if (num_frames == 0 || stride == 0) {
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    MD_DEFER { md_temp_end(temp_scope); };

    SecondaryStructureRenderClass* classes  = md_temp_alloc_array(temp_scope, SecondaryStructureRenderClass, num_frames);
    SecondaryStructureRenderClass* filtered = md_temp_alloc_array(temp_scope, SecondaryStructureRenderClass, num_frames);

    constexpr int window_radius = 2;

    for (size_t seg_idx = 0; seg_idx < stride; ++seg_idx) {
        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            classes[frame_idx] = secondary_structure_render_class(src[frame_idx * stride + seg_idx]);
        }

        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            int counts[3] = {0, 0, 0};
            size_t beg = frame_idx > window_radius ? frame_idx - window_radius : 0;
            size_t end = MD_MIN(frame_idx + window_radius + 1, num_frames);
            for (size_t j = beg; j < end; ++j) {
                counts[(int)classes[j]] += 1;
            }

            SecondaryStructureRenderClass majority = classes[frame_idx];
            int majority_count = 0;
            for (int j = 0; j < 3; ++j) {
                if (counts[j] > majority_count) {
                    majority = (SecondaryStructureRenderClass)j;
                    majority_count = counts[j];
                }
            }

            filtered[frame_idx] = majority_count > (int)((end - beg) / 2) ? majority : classes[frame_idx];
        }

        size_t run_beg = 0;
        while (run_beg < num_frames) {
            size_t run_end = run_beg + 1;
            while (run_end < num_frames && filtered[run_end] == filtered[run_beg]) {
                ++run_end;
            }

            size_t run_len = run_end - run_beg;
            if (run_beg > 0 && run_end < num_frames && filtered[run_beg - 1] == filtered[run_end]) {
                SecondaryStructureRenderClass replacement = filtered[run_beg - 1];
                bool replace_single = run_len <= 1;
                bool replace_short_structured = run_len <= 2 && replacement != SecondaryStructureRenderClass::Coil;
                if (replace_single || replace_short_structured) {
                    for (size_t j = run_beg; j < run_end; ++j) {
                        filtered[j] = replacement;
                    }
                }
            }

            run_beg = run_end;
        }

        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            dst[frame_idx * stride + seg_idx] = secondary_structure_from_render_class(filtered[frame_idx]);
        }
    }
}

// #trajectorydata

SystemRun* instance_run(SystemInstance& inst) {
    return (inst.run >= 0 && (size_t)inst.run < inst.sys->num_runs) ? inst.sys->runs[inst.run] : nullptr;
}

const SystemRun* instance_run(const SystemInstance& inst) {
    return (inst.run >= 0 && (size_t)inst.run < inst.sys->num_runs) ? inst.sys->runs[inst.run] : nullptr;
}

SystemRun* system_run(SystemContext& sys) {
    return sys.num_instances ? instance_run(system_instance(sys)) : nullptr;
}

const SystemRun* system_run(const SystemContext& sys) {
    return sys.num_instances ? instance_run(system_instance(sys)) : nullptr;
}

md_str_t run_path(const SystemContext& sys) {
    const SystemRun* run = system_run(sys);
    return run ? md_str_from_cstr(run->path) : md_str_t{};
}

md_str_t run_attribute_path(char* buf, size_t cap, const SystemRun* run, md_str_t leaf) {
    MD_ASSERT(buf && cap > 0);
    buf[0] = '\0';
    if (!run || run->path[0] == '\0') {
        return {};
    }
    const int len = snprintf(buf, cap, "%s/" MD_STR_FMT, run->path, MD_STR_ARG(leaf));
    if (len <= 0 || (size_t)len >= cap) {
        buf[0] = '\0';
        return {};
    }
    return {buf, (size_t)len};
}

const md_attribute_t* run_attribute(const SystemContext& sys, const SystemRun* run, md_str_t leaf) {
    if (!run || run->path[0] == '\0') return nullptr;
    return md_attributes_find_in(&sys.data.sys.attributes, md_str_from_cstr(run->path), leaf);
}

const md_attribute_t* run_time_axis(const SystemContext& sys, const SystemRun* run) {
    const md_attribute_t* axis = run_attribute(sys, run, MD_STR_LIT("time"));
    return md_attribute_view(axis, MD_ATTRIBUTE_TYPE_F64, 1, 1) ? axis : nullptr;
}

size_t run_num_frames(const SystemContext& sys, const SystemRun* run) {
    const md_attribute_t* axis = run_time_axis(sys, run);
    return axis ? axis->format.shape[0] : 0;
}

const double* run_frame_times(const SystemContext& sys, const SystemRun* run) {
    return (const double*)md_attribute_view(run_time_axis(sys, run), MD_ATTRIBUTE_TYPE_F64, 1, 1);
}

md_unit_t run_time_unit(const SystemContext& sys, const SystemRun* run) {
    const md_attribute_t* axis = run_time_axis(sys, run);
    return axis ? axis->unit : md_unit_none();
}

double run_time_at_frame(const SystemContext& sys, double frame) {
    return time_at_frame(run_frame_times(sys), run_num_frames(sys), frame);
}

// ---------------------------------------------------------------------------------------------------------
// The clock
// ---------------------------------------------------------------------------------------------------------

const SystemInstance* clock_reference(const ApplicationState* app) {
    for (size_t i = 0; i < app->num_systems; ++i) {
        const SystemContext& sc = *app->systems[i];
        for (size_t k = 0; k < sc.num_instances; ++k) {
            if (run_num_frames(sc, instance_run(*sc.instances[k])) > 0) return sc.instances[k];
        }
    }
    return nullptr;
}

md_unit_t clock_unit(const ApplicationState* app) {
    const SystemInstance* ref = clock_reference(app);
    return ref ? run_time_unit(*ref->sys, instance_run(*ref)) : md_unit_none();
}

// The run's index among the system's, which is its index among the runs of its files (RunFiles) too. -1 for none.
static int run_index(const SystemContext& sys, const SystemRun* run) {
    for (size_t i = 0; i < sys.num_runs; ++i) {
        if (sys.runs[i] == run) return (int)i;
    }
    return -1;
}

double run_time_offset(const SystemContext& sys, const SystemRun* run) {
    const int i = run_index(sys, run);
    return i >= 0 && (size_t)i < md_array_size(sys.files.runs) ? sys.files.runs[i].time_offset : 0.0;
}

void run_set_time_offset(SystemContext& sys, const SystemRun* run, double offset) {
    const int i = run_index(sys, run);
    if (i < 0 || (size_t)i >= md_array_size(sys.files.runs) || sys.files.runs[i].time_offset == offset) return;
    sys.files.runs[i].time_offset = offset;
    for (size_t k = 0; k < sys.num_instances; ++k) {
        if (instance_run(*sys.instances[k]) == run) sys.instances[k]->interpolate_system_state = true;
    }
}

bool run_align_start(SystemContext& sys, const SystemRun* run) {
    const ApplicationState* app = sys.app;
    const SystemInstance* ref = clock_reference(app);
    if (!ref || instance_run(*ref) == run) return false;
    const ClockAxis a = run_clock_axis(app, sys, run);
    if (!clock_axis_timed(a)) return false;
    const double start = clock_axis_time(run_clock_axis(app, *ref->sys, instance_run(*ref)), 0.0);
    run_set_time_offset(sys, run, clock_axis_offset_to_start_at(a, start));
    return true;
}

ClockAxis run_clock_axis(const ApplicationState* app, const SystemContext& sys, const SystemRun* run) {
    ClockAxis a;
    a.times = run_frame_times(sys, run);
    a.count = a.times ? run_num_frames(sys, run) : 0;
    a.offset = run_time_offset(sys, run);
    // Frames numbered beside frames numbered as they are, times taken into the clock's unit; the one beside the other
    // not by its times at all
    const md_unit_t from = run_time_unit(sys, run);
    const md_unit_t to   = clock_unit(app);
    double f = 1.0;
    if (md_unit_is_none(from) != md_unit_is_none(to)) f = 0.0;
    else if (!md_unit_is_none(from) && !md_unit_conversion_factor(&f, from, to)) f = 0.0;
    a.scale = f;
    return a;
}

double run_frame_at_clock(const ApplicationState* app, const SystemContext& sys, const SystemRun* run, double time) {
    const ClockAxis a = run_clock_axis(app, sys, run);
    if (a.count == 0) return 0.0;
    if (clock_axis_timed(a)) return clock_axis_frame(a, time);
    // Frame by frame as the reference
    return MD_CLAMP(clock_frame(app), 0.0, (double)(a.count - 1));
}

double run_clock_time(const ApplicationState* app, const SystemContext& sys, const SystemRun* run, double frame) {
    const ClockAxis a = run_clock_axis(app, sys, run);
    if (clock_axis_timed(a)) return clock_axis_time(a, frame);
    // Frame by frame as the reference: the reference's time at the same frame
    const SystemInstance* ref = clock_reference(app);
    if (!ref) return 0.0;
    return clock_axis_time(run_clock_axis(app, *ref->sys, instance_run(*ref)), frame);
}

double instance_current_frame(const SystemInstance& inst) {
    const SystemRun* run = instance_run(inst);
    const size_t num_frames = run_num_frames(*inst.sys, run);
    if (num_frames == 0) return -1.0;
    if (inst.frame >= 0) return (double)MD_MIN((size_t)inst.frame, num_frames - 1);
    return run_frame_at_clock(inst.sys->app, *inst.sys, run, inst.sys->app->animation.time);
}

double clock_frame(const ApplicationState* app) {
    const SystemInstance* ref = clock_reference(app);
    if (!ref) return 0.0;
    const ClockAxis a = run_clock_axis(app, *ref->sys, instance_run(*ref));
    return clock_axis_timed(a) ? clock_axis_frame(a, app->animation.time) : 0.0;
}

void clock_set_frame(ApplicationState* app, double frame) {
    const SystemInstance* ref = clock_reference(app);
    if (!ref) {
        app->animation.time = 0.0;
        return;
    }
    const ClockAxis a = run_clock_axis(app, *ref->sys, instance_run(*ref));
    if (clock_axis_timed(a)) app->animation.time = clock_axis_time(a, frame);
}

void clock_span(const ApplicationState* app, double* beg, double* end) {
    bool any = false;
    *beg = *end = 0.0;
    for (size_t i = 0; i < app->num_systems; ++i) {
        const SystemContext& sc = *app->systems[i];
        for (size_t k = 0; k < sc.num_instances; ++k) {
            const SystemRun* run = instance_run(*sc.instances[k]);
            const size_t n = run_num_frames(sc, run);
            if (n == 0) continue;
            const double b = run_clock_time(app, sc, run, 0.0);
            const double e = run_clock_time(app, sc, run, (double)(n - 1));
            *beg = any ? MD_MIN(*beg, b) : b;
            *end = any ? MD_MAX(*end, e) : e;
            any = true;
        }
    }
}

double clock_frame_duration(const ApplicationState* app) {
    const SystemInstance* ref = clock_reference(app);
    return ref ? clock_axis_frame_duration(run_clock_axis(app, *ref->sys, instance_run(*ref))) : 1.0;
}

void clock_step(ApplicationState* app, double frames) {
    const SystemInstance* ref = clock_reference(app);
    if (!ref) return;
    double beg, end;
    clock_span(app, &beg, &end);
    app->animation.time = clock_step_time(run_clock_axis(app, *ref->sys, instance_run(*ref)), app->animation.time, frames, beg, end);
}

void clock_update(ApplicationState* app) {
    double beg, end;
    clock_span(app, &beg, &end);
    const md_unit_t unit = clock_unit(app);
    if (!md_unit_equal(unit, app->animation.time_unit)) {
        // Another reference, in another unit: the same time in it, else its first
        double f = 1.0;
        const bool convert = !md_unit_is_none(unit) && !md_unit_is_none(app->animation.time_unit) && md_unit_conversion_factor(&f, app->animation.time_unit, unit);
        if (convert) {
            app->animation.time *= f;
            app->timeline.view_range.beg_x *= f;
            app->timeline.view_range.end_x *= f;
            TimelinePeriods& ps = app->timeline.periods;
            for (int i = 0; i < ps.count; ++i) {
                ps.items[i].beg *= f;
                ps.items[i].end *= f;
            }
        } else {
            app->animation.time = beg;
            app->timeline.view_range = {beg, end > beg ? end : beg + 1.0};
        }
        app->animation.time_unit = unit;
    }
    app->animation.time = MD_CLAMP(app->animation.time, beg, end);
}

void playback_range(const ApplicationState* app, double* beg, double* end) {
    clock_span(app, beg, end);
    const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, app->animation.loop_period);
    if (!p || p->follow) return;
    double pb, pe;
    timeline_period_span(&pb, &pe, *p, app->animation.time);
    pb = MD_MAX(pb, *beg);
    pe = MD_MIN(pe, *end);
    if (pe > pb) {
        *beg = pb;
        *end = pe;
    }
}

void playback_toggle(ApplicationState* app) {
    if (app->animation.mode == PlaybackMode::Playing) {
        app->animation.mode = PlaybackMode::Stopped;
        for (size_t i = 0; i < app->num_systems; ++i) system_instances_flag(*app->systems[i], MolBit_ClearVelocity);
        return;
    }
    app->animation.mode = PlaybackMode::Playing;
    // A bounce's turn back lasts as long as it bounces
    if (app->animation.loop != PlaybackLoop_Bounce) app->animation.direction = 1;
    // From where it begins in the direction played: when it stopped at the end, or is outside the span it goes over
    double beg, end;
    playback_range(app, &beg, &end);
    const bool forward = app->animation.fps * (float)app->animation.direction >= 0.0f;
    double& t = app->animation.time;
    if (forward && (t >= end || t < beg)) t = beg;
    else if (!forward && (t <= beg || t > end)) t = end;
}

void playback_update(ApplicationState* app, double dt) {
    double beg, end;
    playback_range(app, &beg, &end);
    const double delta = dt * (double)app->animation.fps * clock_frame_duration(app);
    if (app->animation.direction != 1 && app->animation.direction != -1) app->animation.direction = 1;
    if (app->animation.loop != PlaybackLoop_Bounce) app->animation.direction = 1;
    // Outside the span it goes over (another period chosen to loop over while playing): from where that begins
    const bool forward = delta * (double)app->animation.direction >= 0.0;
    double& t = app->animation.time;
    if (t < beg || t > end) t = forward ? beg : end;
    if (!clock_advance(&app->animation.time, delta, beg, end, app->animation.loop, &app->animation.direction)) {
        app->animation.mode = PlaybackMode::Stopped;
    }
}

bool run_period_frames(int* beg, int* end, const ApplicationState* app, const SystemContext& sys, const SystemRun* run, const TimelinePeriod& period) {
    const int n = (int)run_num_frames(sys, run);
    if (n == 0) return false;
    const ClockAxis a = run_clock_axis(app, sys, run);
    if (clock_axis_timed(a)) return timeline_period_frames(beg, end, period, app->animation.time, a);
    // Frame by frame as the reference
    const SystemInstance* ref = clock_reference(app);
    int b, e;
    if (!ref || !timeline_period_frames(&b, &e, period, app->animation.time, run_clock_axis(app, *ref->sys, instance_run(*ref))) || b > n - 1) return false;
    *beg = b;
    *end = MD_MIN(e, n - 1);
    return true;
}

const SystemInstance& picking_hit_instance(const SystemContext& sc, const PickingHit& hit) {
    const int k = picking_key_instance(hit.key);
    if (k >= 0 && (size_t)k < sc.num_instances) return *sc.instances[k];
    for (size_t i = 0; i < md_array_size(sc.representation.reps); ++i) {
        if (sc.representation.reps[i].node.id == hit.key) return representation_instance(sc, sc.representation.reps[i]);
    }
    return system_instance(sc);
}

const md_secondary_structure_t* displayed_secondary_structure(const SystemInstance& inst) {
    const SystemContext& sys = *inst.sys;
    const md_system_state_t* state = &inst.state;
    const SystemRun* run = instance_run(inst);
    const SystemTrajectoryData::Render empty = {};
    const auto& render = run ? run->trajectory_data.secondary_structure_render : empty;
    if (render.data && render.stride == sys.data.sys.protein_backbone.segment.count && md_state_has_frame(state)) {
        const size_t frame = (size_t)md_state_frame_nearest(state);
        if ((frame + 1) * render.stride <= render.count) {
            return render.data + frame * render.stride;
        }
    }
    return md_system_state_secondary_structure(state, &sys.data.sys);
}

static const md_str_t frame_extract_paths[] = { MD_STR_INIT("atom/position"), MD_STR_INIT("unitcell") };

bool extract_frame(const SystemContext& sys, const SystemRun* run, int64_t frame, md_system_state_t* out) {
    MD_ASSERT(out);
    if (!run) return false;
    md_system_extract_t* ex = md_system_extract_begin(&sys.data.sys, md_str_from_cstr(run->path),
        frame_extract_paths, MD_COUNTOF(frame_extract_paths), md_get_heap_allocator());
    if (!ex) {
        return false;
    }
    const bool ok = md_system_extract_frame(ex, frame, out);
    md_system_extract_end(ex);
    return ok;
}


static void end_frame_extracts(SystemRun& run) {
    for (size_t i = 0; i < MD_COUNTOF(run.frame_extract); ++i) {
        md_system_extract_end(run.frame_extract[i]);
        run.frame_extract[i] = nullptr;
    }
}

// "run/<stem>" for a trajectory file. The stem is folded to letters, digits, '_' and '-' so that it
// is one path segment whatever the file is called, and it is taken from the file name alone so the
// same file gives the same run - and the same attribute ids - every time it is loaded.
static void run_path_from_file(char* buf, size_t cap, md_str_t path) {
    md_str_t file = md_path_stem(path);
    size_t len = (size_t)snprintf(buf, cap, "run/");
    for (size_t i = 0; i < file.len && len + 1 < cap; ++i) {
        const char c = file.ptr[i];
        const bool keep = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
        buf[len++] = keep ? c : '_';
    }
    if (len == 4) {
        len += (size_t)snprintf(buf + len, cap - len, "trajectory");
    }
    buf[MD_MIN(len, cap - 1)] = '\0';
}

// The path made unique among the runs of the system: another run of a file of the same name (two md.xtc of two
// replicas) gets "_2", "_3" after it, as its attributes would otherwise be the other's
static void run_path_make_unique(char* buf, size_t cap, const md_attributes_t& attributes) {
    const size_t len = strnlen(buf, cap);
    for (int n = 2; n < 1000 && md_attributes_find_in(&attributes, md_str_from_cstr(buf), MD_STR_LIT("time")); ++n) {
        snprintf(buf + len, cap - len, "_%d", n);
    }
}

// What a run holds, released: its task stopped, its readers ended (they hold its files open), what was derived along it
// and its attributes (its frame axis, the quantities per frame, the series loaded along it). What reads it was stopped
// by the caller; its files in the list are the caller's.
static void free_run(ApplicationState* state, SystemContext& sys, SystemRun& run) {
    run.backbone_computations.cancel_and_wait();
    end_frame_extracts(run);

    // backbone_angles.data and secondary_structure.data are VIEWS into sys.attributes and are not md_arrays - freeing
    // them here would hand the wrong pointer to the array allocator. The table owns them; dropping the views is all
    // that is owed. The denoised render copy is an ordinary array owned here.
    run.trajectory_data.backbone_angles = {};
    run.trajectory_data.secondary_structure = {};
    md_array_free(run.trajectory_data.secondary_structure_render.data, state->allocator.persistent);
    run.trajectory_data.secondary_structure_render = {};

    if (run.path[0]) {
        md_attributes_remove_prefix(&sys.data.sys.attributes, md_str_from_cstr(run.path));
        run.path[0] = '\0';
    }
}

// What was derived of the runs of the system, by whoever holds on to it, goes: before the runs it is of change
static void runs_changing(ApplicationState* state, SystemContext& sys) {
    viamd::event_system_broadcast_event(viamd::EventType_ViamdTrajectoryFree, viamd::EventPayloadType_System, &sys);
    // The orientation reference was read from the first frame of the run shown
    sys.operations.initial_frame.valid = false;
    // Converted axes and histograms, and the evaluation of the script, of what is about to go or no longer be shown.
    // The attribute versions they are keyed on only mean something within one table, and the next run may be
    // allocated where this was; taking over from an evaluation of another run (script::eval_create) would carry its
    // values over to this one. The tasks which use them were stopped by the caller.
    series_cache_free(state);
    script_eval_free(state);
}

void free_trajectory_data(ApplicationState* state, SystemContext& sys) {
    MD_ASSERT(state);
    if (sys.num_runs > 0) runs_changing(state, sys);
    for (size_t i = 0; i < sys.num_runs; ++i) {
        free_run(state, sys, *sys.runs[i]);
        delete sys.runs[i];
        sys.runs[i] = nullptr;
    }
    sys.num_runs = 0;
    for (size_t k = 0; k < sys.num_instances; ++k) {
        sys.instances[k]->run = -1;
        sys.instances[k]->frame = -1;
        clear_frame_cache(&sys.instances[k]->frame_cache);
    }
    system_instances_update(sys, MolBit_ClearVelocity);
    // The runs, and with them the files they were opened from and the series loaded along them
    file_list_drop_runs(&sys.files);
}

// ### INSTANCES ###

// What an instance has on the GPU, and its coordinates, made for the system as loaded
static void instance_init_data(SystemInstance& inst) {
    SystemContext& sys = *inst.sys;
    md_allocator_i* alloc = sys.app->allocator.persistent;
    inst.gl_mol = md_gl_mol_create(&sys.data.sys);
    // A new system is shown as it is: the subsets of the parts of an assembly are made when one is shown
    inst.scene.assembly = -1;
    md_array_shrink(inst.scene.assembly_transforms, 0);
    md_array_shrink(inst.scene.assembly_copy_offset, 0);
    // Before the representations, which each make a representation of the glyphs
    cartoon_glyphs_init(&inst.cartoon_glyphs, sys.data.sys, inst.gl_mol, alloc);
    init_frame_cache(&inst.frame_cache, sys.data.sys.atom.count, alloc);

    // The system's own coordinates until they are made from its source (interpolate_system_state), with the backbone
    // they carry
    md_system_state_copy(&inst.state, &sys.data.state);
    if (md_system_state_backbone_compute(&inst.state, &sys.data.sys)) {
        upload_secondary_structure_weights(inst);
    }
    md_float3x3_t A;
    md_unitcell_A_extract_float(A.elem, &inst.state.unitcell);
    const md_float3_t c = md_float3x3_mul_float3(A, md_float3_set(0.5f, 0.5f, 0.5f));
    inst.unitcell_transform = md_float4x4_translate(-c.x, -c.y, -c.z);
    inst.state_rotation = md_float4x4_ident();
    inst.last_interpolated_nearest_frame = -1;
    inst.last_interpolated_frame = -1.0;
    inst.interpolate_system_state = true;
    inst.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
}

static void instance_free_data(SystemInstance& inst) {
    md_allocator_i* alloc = inst.sys->app->allocator.persistent;
    assembly_subsets_free(inst);
    md_gl_mol_destroy(inst.gl_mol);
    // Not left naming the destroyed one: what runs before the next system is loaded (a workspace creating its
    // representations, a frame after a load that failed) sees no molecule rather than a stale one
    inst.gl_mol = {};
    cartoon_glyphs_free(&inst.cartoon_glyphs);
    free_frame_cache(&inst.frame_cache, alloc);
    md_system_state_free(&inst.state);
    inst.flags_hash = 0;
}

SystemInstance* instance_add(SystemContext& sys) {
    if (sys.num_instances >= INSTANCE_MAX) return nullptr;
    SystemInstance* inst = new SystemInstance();
    inst->sys = &sys;
    inst->id = sys.next_instance_id++;
    inst->state.alloc = sys.app->allocator.persistent;
    sys.instances[sys.num_instances++] = inst;
    if (sys.data.sys.atom.count > 0) instance_init_data(*inst);
    return inst;
}

// Everything it holds, and itself
static void instance_destroy(SystemInstance* inst) {
    md_allocator_i* alloc = inst->sys->app->allocator.persistent;
    instance_free_data(*inst);
    md_array_free(inst->scene.periodic_transforms, alloc);
    md_array_free(inst->scene.periodic_transforms_prev, alloc);
    md_array_free(inst->scene.assembly_transforms, alloc);
    md_array_free(inst->scene.assembly_copy_offset, alloc);
    delete inst;
}

void instance_remove(SystemContext& sys, size_t i) {
    if (i >= sys.num_instances || sys.num_instances == 1) return;
    interrupt_async_tasks(sys.app);
    SystemInstance* inst = sys.instances[i];
    remove_instance_representations(sys, *inst);
    instance_destroy(inst);
    for (size_t k = i; k + 1 < sys.num_instances; ++k) sys.instances[k] = sys.instances[k + 1];
    sys.instances[--sys.num_instances] = nullptr;
    // The first may be another now: the timeline, the script and the tools follow it
    if (i == 0) instance_show(system_instance(sys), system_instance(sys).run, system_instance(sys).frame);
}

SystemInstance* instance_find(SystemContext& sys, uint32_t id) {
    for (size_t k = 0; k < sys.num_instances; ++k) if (sys.instances[k]->id == id) return sys.instances[k];
    return nullptr;
}

const SystemInstance* instance_find(const SystemContext& sys, uint32_t id) {
    for (size_t k = 0; k < sys.num_instances; ++k) if (sys.instances[k]->id == id) return sys.instances[k];
    return nullptr;
}

void system_instances_update(SystemContext& sys, uint32_t dirty_bits) {
    for (size_t k = 0; k < sys.num_instances; ++k) {
        sys.instances[k]->interpolate_system_state = true;
        sys.instances[k]->dirty_gpu_buffers |= dirty_bits;
    }
}

void instance_show(SystemInstance& inst, int run, int64_t frame) {
    SystemContext& sys = *inst.sys;
    ApplicationState* app = sys.app;
    if (run < 0 || (size_t)run >= sys.num_runs) run = -1;
    // The timeline shows all of the clock's span as long as it showed all of it before: a longer run widens it
    double span_beg, span_end;
    clock_span(app, &span_beg, &span_end);
    const bool view_all = app->timeline.view_range.beg_x <= span_beg && app->timeline.view_range.end_x >= span_end;
    const bool first = &inst == &system_instance(sys);
    // What was derived along the run the first instance shows is of another one now
    if (first && run != inst.run) runs_changing(app, sys);
    inst.run = run;
    inst.frame = run >= 0 ? frame : -1;
    inst.last_interpolated_nearest_frame = -1;
    inst.last_interpolated_frame = -1.0;
    clear_frame_cache(&inst.frame_cache);
    inst.interpolate_system_state = true;
    inst.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
    if (view_all) {
        clock_update(app);
        clock_span(app, &span_beg, &span_end);
        app->timeline.view_range = {span_beg, span_end};
    }
    if (!first) return;

    // The timeline shows the clock (clock.h), whose reference this is likely to be: all of it, in its unit, which the
    // display unit only changes the printing of (frame_axis.h)
    if (&sys == &app->system_context && run_num_frames(sys, instance_run(inst)) > 0) {
        clock_update(app);
        double beg, end;
        clock_span(app, &beg, &end);
        app->timeline.view_range = {beg, end};
    }
    // The script reads along the run shown, and what it measured was of the one before
    app->script.compile_ir = true;
    viamd::event_system_broadcast_event(viamd::EventType_ViamdTrajectoryInit, viamd::EventPayloadType_System, &sys);
}

void close_run(ApplicationState* app, SystemContext& sys, size_t i) {
    MD_ASSERT(app);
    if (i >= sys.num_runs) return;
    interrupt_async_tasks(app);
    runs_changing(app, sys);
    free_run(app, sys, *sys.runs[i]);
    delete sys.runs[i];
    for (size_t j = i; j + 1 < sys.num_runs; ++j) sys.runs[j] = sys.runs[j + 1];
    sys.runs[--sys.num_runs] = nullptr;
    file_list_remove_run(&sys.files, i);

    // An instance of it goes with it, as it came with it, but for the first, which shows the system's own coordinates;
    // those of the runs after it show the same run in its new place
    for (size_t k = sys.num_instances; k-- > 0;) {
        SystemInstance& inst = *sys.instances[k];
        if (inst.run == (int)i) {
            if (k > 0) instance_remove(sys, k);
            else       instance_show(inst, -1);
        } else if (inst.run > (int)i) {
            inst.run -= 1;
        }
    }
}

// The backbone angles and secondary structure of every frame of the run, computed on the pool. Started when the run is
// opened, and again once the backbones are inferred (they are not there before: SystemInference)
static void trajectory_backbone_start(SystemContext& sys, SystemRun& run) {
    ApplicationState* data = sys.app;
    const size_t num_frames = run_num_frames(sys, &run);
    if (num_frames == 0) return;
    char path_buf[256];
    if (sys.data.sys.protein_backbone.segment.count > 0) {
        // The angles and the per frame secondary structure are TEMPORAL attributes: the table
        // owns the storage and the trajectory_data fields below are views onto it. That is the
        // point of the move - 'stride' and 'count' were a hand rolled shape {F,S}, and the one
        // place they could disagree with the buffer was here, in three lines repeated per
        // quantity. Now the shape IS the declaration, and md_attributes_iter answers "what
        // varies over time in this dataset" without anybody maintaining a second list.
        //
        // Both are declared together so the frame axis and the segment axis are stated once for
        // the pair. Resident rather than computed on demand, because the ramachandran density
        // reads EVERY frame at once - a per plane provider would recompute the whole trajectory
        // each time the window redraws. Reserved-and-filled-lazily is the interesting third
        // option, and it only starts paying at trajectory lengths where this allocation hurts.
        const size_t num_segments = sys.data.sys.protein_backbone.segment.count;
        md_attributes_t* attributes = &sys.data.sys.attributes;

        // One still under way writes through the pointers replaced below: stopped first
        run.backbone_computations.cancel_and_wait();

        // Both are checked against "<run>/time", published above, so a declaration whose
        // outermost extent disagrees is refused at the point the mistake is made rather than
        // surviving until something reads past the end of it.

        // md_secondary_structure_t is a 4 byte enum, so I32 is the storage and the PATH is what
        // tells a consumer these are labels rather than numbers to average - which is the rule
        // md_system.h already states for integral attributes.
        const md_attribute_desc_t ss_desc = {
            .path   = run_attribute_path(path_buf, sizeof(path_buf), &run, MD_STR_LIT("backbone/secondary_structure")),
            .format = {
                .type = MD_ATTRIBUTE_TYPE_I32, .components = 1,
                .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
            },
            .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit   = md_unit_none(),
            .label  = MD_STR_INIT("Secondary Structure"),
        };
        md_attribute_id_t ss_id = md_attributes_replace(attributes, &ss_desc);

        run.trajectory_data.secondary_structure.data   = (md_secondary_structure_t*)md_attributes_data(attributes, ss_id, MD_ATTRIBUTE_TYPE_I32);
        run.trajectory_data.secondary_structure.stride = num_segments;
        run.trajectory_data.secondary_structure.count  = num_segments * num_frames;

        const md_attribute_format_t angle_format = {
            .type = MD_ATTRIBUTE_TYPE_F32, .components = 2,   // phi, psi - one value, two parts
            .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
        };
        const md_attribute_desc_t angle_desc = {
            .path   = run_attribute_path(path_buf, sizeof(path_buf), &run, MD_STR_LIT("backbone/angle")),
            .format = angle_format,
            .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit   = md_unit_radian(),
            .label  = MD_STR_INIT("Backbone Angles"),
        };
        md_attribute_id_t angle_id = md_attributes_replace(attributes, &angle_desc);

        // md_backbone_angles_t is exactly {float phi; float psi;}, so the table's storage IS a
        // md_backbone_angles_t array - no repacking, and the existing consumers keep their type.
        run.trajectory_data.backbone_angles.data   = (md_backbone_angles_t*)md_attributes_data(attributes, angle_id, MD_ATTRIBUTE_TYPE_F32);
        run.trajectory_data.backbone_angles.stride = num_segments;
        run.trajectory_data.backbone_angles.count  = num_segments * num_frames;

        // The denoised copy stays an ordinary array: it is a PRESENTATION smoothing of the one
        // above, read only by the renderer, and publishing it would put two answers to "what is
        // the secondary structure at frame f" in the same table.
        run.trajectory_data.secondary_structure_render.stride = sys.data.sys.protein_backbone.segment.count;
        run.trajectory_data.secondary_structure_render.count = sys.data.sys.protein_backbone.segment.count * num_frames;
        md_array_resize(run.trajectory_data.secondary_structure_render.data, sys.data.sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
        md_memset(run.trajectory_data.secondary_structure_render.data, 0, md_array_bytes(run.trajectory_data.secondary_structure_render.data));

        // Launch work to compute the values
        uint64_t time = (uint64_t)md_tick_now();
        run.backbone_computations.start_for(MD_STR_LIT("Backbone Operations"), (uint32_t)num_frames, 1, [s = &sys, r = &run](uint32_t range_beg, uint32_t range_end) {
            const md_system_t* sys = &s->data.sys;

            md_temp_scope_t temp = md_temp_begin();
            MD_DEFER { md_temp_end(temp); };

            md_system_state_t frame_state = { .alloc = temp.arena };
				md_system_state_init(&frame_state, sys->atom.count);

            // One context for the range, so the run's files stay open across its frames.
            md_system_extract_t* ex = md_system_extract_begin(sys, md_str_from_cstr(r->path),
                frame_extract_paths, MD_COUNTOF(frame_extract_paths), md_get_heap_allocator());
            if (!ex) {
                return;
            }
            for (uint32_t frame_idx = range_beg; frame_idx < range_end && !task_system::cancelled(); ++frame_idx) {
                md_backbone_angles_t* bb_dst = r->trajectory_data.backbone_angles.data + r->trajectory_data.backbone_angles.stride * frame_idx;
                md_secondary_structure_t* ss_dst = r->trajectory_data.secondary_structure.data + r->trajectory_data.secondary_structure.stride * frame_idx;

                if (!md_system_extract_frame(ex, frame_idx, &frame_state)) {
                    continue;
                }
                md_backbone_angles_compute(bb_dst, r->trajectory_data.backbone_angles.stride, frame_state.xyz, &frame_state.unitcell, &sys->protein_backbone);
                md_dssp_compute(ss_dst, r->trajectory_data.secondary_structure.stride, frame_state.xyz, &frame_state.unitcell, &sys->protein_backbone);
            }
            md_system_extract_end(ex);
        }, [s = &sys, r = &run, t0 = time, num_frames](task_system::Status status) {
            // Stopped part of the way (a load, or the next trajectory): the data is not for this one, or not all there
            if (status != task_system::Status::Completed) return;
            secondary_structure_render_denoise(
                r->trajectory_data.secondary_structure_render.data,
                r->trajectory_data.secondary_structure.data,
                num_frames,
                r->trajectory_data.secondary_structure.stride);

            uint64_t t1 = (uint64_t)md_tick_now();
            double elapsed = md_tick_to_seconds(t1 - t0);
            MD_LOG_INFO("Finished computing trajectory data (%.2fs)", elapsed);
            // The fill is finished, so the two temporal attributes have changed - which is the
            // one moment a consumer caching something derived from them needs to hear about.
            // The table carries that now; the hand rolled fingerprints beside the data are gone,
            // because a stamp that lives next to storage the table owns is a second answer
            // waiting to disagree with it.
            md_attributes_t* attributes = &s->data.sys.attributes;
            char buf[256];
            md_attributes_touch(attributes, md_attributes_id_from_path(run_attribute_path(buf, sizeof(buf), r, MD_STR_LIT("backbone/angle"))));
            md_attributes_touch(attributes, md_attributes_id_from_path(run_attribute_path(buf, sizeof(buf), r, MD_STR_LIT("backbone/secondary_structure"))));
            r->trajectory_data.secondary_structure_render.fingerprint = generate_fingerprint();

				// The data just went from placeholder to real values for every frame - force the
				// next interpolate_system_state to actually push it, even if playback is paused on
				// the same nearest frame it was on before this task started.
				for (size_t k = 0; k < s->num_instances; ++k) {
                    if (instance_run(*s->instances[k]) == r) s->instances[k]->last_interpolated_nearest_frame = -1;
                }
                system_instances_update(*s, MolBit_ClearVelocity);
            flag_all_representations_as_dirty(*s);
        });
    }
}

void assembly_transforms(md_array(md_float4x4_t)* transforms, md_array(uint32_t)* copy_offset, const md_system_t& sys, int assembly, md_allocator_i* alloc) {
    MD_ASSERT(transforms);
    md_array_shrink(*transforms, 0);
    if (copy_offset) md_array_shrink(*copy_offset, 0);
    const md_assembly_data_t& a = sys.assembly;
    if (assembly < 0 || (size_t)assembly >= a.count) return;
    if (copy_offset) md_array_push(*copy_offset, 0, alloc);
    // The groups and their operators from the offsets as they are (assembly_group)
    for (size_t p = a.group_offset[assembly]; p < a.group_offset[assembly + 1]; ++p) {
        const AssemblyGroup group = assembly_group(a, p);
        for (const uint32_t* o = group.oper; o < group.oper_end; ++o) {
            md_array_push(*transforms, a.oper_transform[*o], alloc);
        }
        if (copy_offset) md_array_push(*copy_offset, (uint32_t)md_array_size(*transforms), alloc);
    }
}

void assembly_subsets_free(SystemInstance& inst) {
    for (size_t i = 0; i < md_array_size(inst.scene.assembly_subsets); ++i) {
        if (inst.scene.assembly_subsets[i].id) md_gl_subset_destroy(inst.scene.assembly_subsets[i]);
    }
    for (size_t i = 0; i < md_array_size(inst.scene.assembly_glyph_subsets); ++i) {
        if (inst.scene.assembly_glyph_subsets[i].id) md_gl_subset_destroy(inst.scene.assembly_glyph_subsets[i]);
    }
    md_array_shrink(inst.scene.assembly_subsets, 0);
    md_array_shrink(inst.scene.assembly_glyph_subsets, 0);
    inst.scene.assembly_subsets_of  = -1;
    inst.scene.assembly_subsets_mol = 0;
}

void assembly_subsets_update(SystemInstance& inst) {
    SystemContext& sc = *inst.sys;
    ApplicationState* app = sc.app;
    MD_ASSERT(app);
    const int assembly = inst.scene.assembly;
    const uint32_t mol_id = inst.gl_mol.id;
    if (assembly == inst.scene.assembly_subsets_of && mol_id == inst.scene.assembly_subsets_mol) return;
    assembly_subsets_free(inst);
    const md_system_t& sys = sc.data.sys;
    if (assembly < 0 || (size_t)assembly >= sys.assembly.count || !mol_id) return;
    inst.scene.assembly_subsets_of  = assembly;
    inst.scene.assembly_subsets_mol = mol_id;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_allocator_i* alloc = app->allocator.persistent;
    const md_assembly_data_t& a = sys.assembly;
    const CartoonGlyphs& glyphs = inst.cartoon_glyphs;
    const bool has_glyphs = cartoon_glyphs_valid(glyphs);
    uint32_t* atoms = md_temp_alloc_array(temp, uint32_t, MD_MAX(sys.atom.count, (size_t)1));
    // The groups and their instances from the offsets as they are (assembly_group), not from a range of them
    for (size_t p = a.group_offset[assembly]; p < a.group_offset[assembly + 1]; ++p) {
        const AssemblyGroup group = assembly_group(a, p);
        md_gl_subset_t subset = {};
        md_gl_subset_t glyph_subset = {};
        // A group of all instances (each once) shows the whole system
        if ((size_t)(group.inst_end - group.inst) < sys.instance.count) {
            uint32_t count = 0;
            for (const md_instance_idx_t* k = group.inst; k < group.inst_end; ++k) {
                const md_urange_t r = md_system_instance_atom_range(&sys, (size_t)*k);
                const uint32_t num = r.end > r.beg ? r.end - r.beg : 0;
                for (uint32_t i = 0; i < num; ++i) atoms[count++] = r.beg + i;
            }
            subset = md_gl_subset_create(inst.gl_mol);
            md_gl_subset_set_atoms(subset, atoms, count);
            if (has_glyphs) {
                glyph_subset = md_gl_subset_create(glyphs.gl_mol);
                md_gl_subset_set_derived(glyph_subset, subset);
            }
        }
        md_array_push(inst.scene.assembly_subsets, subset, alloc);
        md_array_push(inst.scene.assembly_glyph_subsets, glyph_subset, alloc);
    }
}

void init_system_data(ApplicationState* data, SystemContext& sys) {
    if (sys.data.sys.atom.count) {
        md_bitset_free(&sys.operations.recenter_query.mask, data->allocator.persistent);
        md_bitset_free(&sys.operations.selection_mask, data->allocator.persistent);
        sys.operations.recenter_query.filter.result = {};
        sys.operations.recenter_query.evaluated_version = 0;
        sys.operations.recenter_query.ir_fingerprint = 0;
        sys.operations.initial_frame.valid = false;
        recenter_mark_query_dirty(sys);
        {
            // Of the bonds of another system
            auto& bt = sys.operations.bond_tracking;
            md_bond_delta_free(&bt.delta, data->allocator.persistent);
            md_chem_delta_free(&bt.chem, data->allocator.persistent);
            bt.delta = {};
            bt.key = 0;
            bt.formed_capacity = 0;
            bt.ring_capacity = 0;
            md_bitset_free(&bt.scope_mask, data->allocator.persistent);
            bt.scope.result = {};
            bt.scope_dirty = true;
        }

        // Each instance its molecule on the GPU and its coordinates (before the representations, which each make a
        // representation of the glyphs)
        for (size_t k = 0; k < sys.num_instances; ++k) instance_init_data(*sys.instances[k]);

        md_float3x3_t A;
		md_unitcell_A_extract_float(A.elem, &sys.data.state.unitcell);

        md_float3_t aabb_min = {};
        md_float3_t aabb_max = {};
        if (md_system_state_has_coords(&sys.data.state)) {
            md_geom_aabb(aabb_min.elem, aabb_max.elem, (const float*)sys.data.state.xyz, 0, nullptr, sys.data.state.num_atoms);
        }

        const md_float3_t cell_ext = md_float3x3_mul_float3(A, md_float3_set1(1.0f));
        const float max_cell_ext = md_float3_reduce_max(cell_ext);
        const float max_aabb_ext = md_float3_reduce_max(md_float3_sub(aabb_max, aabb_min));

        data->view.camera.near_plane = 1.0f;
        data->view.camera.far_plane = 100000.0f;
        data->view.trackball_param.max_distance = MD_MAX(data->view.trackball_param.max_distance, MD_MAX(max_cell_ext, max_aabb_ext) * 10.0f);

        // The view goes to the system when it is the only one: one beside others is framed with them where it was
        // loaded from (load_file_into)
        if (system_count(*data) <= 1 && md_system_state_has_coords(&sys.data.state)) {
            ViewTransform default_view = {};
            reset_view(&default_view, sys.data.state, (const md_bitset_t*)nullptr, data->view.camera.fov_y, data->view.region);
            data->view.trackball_param.max_distance = MD_MAX(max_cell_ext, max_aabb_ext) * 10.0f;
            data->view.target = default_view;
            data->view.camera = default_view;
        }
        

    }
    // The components are told once the system is inferred, when it is inferred in the background: what they derive from
    // it (the particle and residue types of the dataset window and its overrides, the pucker rows) is of the inferred one
    if (!sys.data.inference.pending) {
        viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemInit, viamd::EventPayloadType_System, &sys);
    }

    init_all_representations(sys);
    data->script.compile_ir = true;

    // The first excited state's natural transition orbitals, read on a worker thread while the geometry is
    // shown: the state a representation of them starts at. Of the first system, whose quantum chemistry is evaluated.
    if (&sys == &data->system_context) {
        es_nto_ready(data, 0);

        // The SCF density the CPU paths read (the ESP, the critical points: the total density, the alpha one without
        // it), into the density cache on the task pool while the geometry is shown: the first of them finds it there
        // rather than reading it then (seconds on a protein). Only for a system with a basis to make a field of it with.
        const md_system_t& msys = sys.data.sys;
        const md_str_t path = md_attributes_find(&msys.attributes, es_path::total_density) ? es_path::total_density : es_path::alpha_density;
        const bool basis = md_attributes_find(&msys.attributes, MD_STR_LIT("basis/shell/atom_index")) != nullptr &&
                           md_attributes_find(&msys.attributes, MD_STR_LIT("basis/primitive/exponent")) != nullptr;
        if (basis) density_cache_preload(density_cache_of(data), msys, path);
    }
}

static void system_inference_free(SystemContext& sys);

void free_system_data(ApplicationState* data, SystemContext& sys) {
    MD_ASSERT(data);
    interrupt_async_tasks(data);
    // Before the arena is rewound: the inferred arrays live in arenas of their own, which go with the system
    system_inference_free(sys);
    es_nto_warm_free(sys);     // its tasks were waited for above
    density_cache_clear(sys.data.density_cache);   // the densities of the system: nothing reads them any more

    // The arena OUTLIVES the system it backs, and the zeroing below is why there used to be a
    // second handle to it: md_memset over md_system_t clears sys.alloc along with everything else, so
    // the only pointer to the arena has to be carried across by hand. Rewound, not destroyed - the
    // next load reuses it.
    md_allocator_i* sys_alloc = sys.data.sys.alloc;
    // What was inferred is in an arena of the persistent allocator's pages, which the arena's rewind does not give back
    md_system_derived_free(&sys.data.sys);
    md_arena_allocator_reset(sys_alloc);
    md_memset(&sys.data.sys, 0, sizeof(sys.data.sys));
    md_memset(&sys.data.state, 0, sizeof(sys.data.state));
    sys.data.sys.alloc   = sys_alloc;
    sys.data.state.alloc = sys_alloc;

    md_array_free(sys.operations.initial_frame.rel_xyzw, data->allocator.persistent);
    sys.operations.initial_frame.rel_xyzw = nullptr;
    sys.operations.initial_frame.valid = false;

    // What each instance has on the GPU and its coordinates are of the system
    for (size_t k = 0; k < sys.num_instances; ++k) instance_free_data(*sys.instances[k]);

    // The dataset's GPU data goes with the dataset, for the same reason gl_mol does.
    system_gpu_data_free(sys);

    // Everything loaded was loaded onto the system
    file_list_clear(&sys.files);



    md_bitmask_clear(&sys.selection.selection_mask);
    md_bitmask_clear(&sys.selection.highlight_mask);
    if (data->selection.single_selection_sequence.system == sys.id) {
        single_selection_sequence_clear(&data->selection.single_selection_sequence);
    }

    // The script is compiled against the first system and evaluated in it (until it is per system): what it computed
    // goes with that one alone
    if (&sys == &data->system_context) {
        // The frames selected are of its run, which the timeline shows
        frame_selection_clear(data);
        md_bitmask_clear(&data->selection.frame_highlight_mask);
        script_vis_reset(data);

        // The tasks were waited for above: with the fields cleared, nothing uses any IR
        data->script.ir = nullptr;
        data->script.eval_ir = nullptr;
        script_ir_collect(data);
        script_eval_free(data);
    }

    sys.name[0] = '\0';

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemFree, viamd::EventPayloadType_System, &sys);
}


// ## The systems

void system_context_init(SystemContext* sys, ApplicationState* app) {
    MD_ASSERT(sys && app);
    md_allocator_i* alloc = app->allocator.persistent;
    sys->app = app;
    // One arena per system, held by the system and the state that share it. Both are rewound rather than destroyed
    // between loads (free_system_data), so this handle is set once and never reassigned.
    sys->data.sys.alloc   = md_arena_allocator_create(alloc, MD_MEGABYTES(1));
    sys->data.state.alloc = sys->data.sys.alloc;
    // What is inferred of the system (its bonds, molecules, rings, chemistry, ...) is in an arena of its own, replaced
    // whole whenever it is inferred again (an edit of its bonds, Recompute Bonds), whose pages come from where this one
    // takes its own: the persistent allocator, which gives them back (md_system_t, ALLOCATORS). Destroyed with the
    // system in free_system_data.
    // The masks are edited in place. Over any atom index: they hold only the words they write, and a mask over a
    // system's atoms needs no resizing when another is loaded.
    md_bitmask_init(&sys->selection.selection_mask,       0, MD_BITSET_MAX_INDEX, alloc);
    md_bitmask_init(&sys->selection.highlight_mask,       0, MD_BITSET_MAX_INDEX, alloc);
    md_bitmask_init(&sys->representation.visibility_mask, 0, MD_BITSET_MAX_INDEX, alloc);
    // Sets, replaced whole
    sys->operations.selection_mask = md_bitset_empty_set();
    sys->operations.recenter_query.mask = md_bitset_empty_set();
    // A system has an instance at least: the first, which the timeline, the script and the tools read
    instance_add(*sys);
}

void system_context_destroy(SystemContext* sys) {
    MD_ASSERT(sys && sys->app);
    md_allocator_i* alloc = sys->app->allocator.persistent;
    md_bitmask_free(&sys->selection.selection_mask);
    md_bitmask_free(&sys->selection.highlight_mask);
    md_bitmask_free(&sys->representation.visibility_mask);
    md_bitset_free(&sys->operations.selection_mask, alloc);
    md_bitset_free(&sys->operations.recenter_query.mask, alloc);
    md_bitset_free(&sys->operations.bond_tracking.scope_mask, alloc);
    file_list_free(&sys->files, alloc);
    md_array_free(sys->representation.reps, alloc);
    md_array_free(sys->representation.groups, alloc);
    for (size_t k = 0; k < sys->num_instances; ++k) {
        instance_destroy(sys->instances[k]);
        sys->instances[k] = nullptr;
    }
    sys->num_instances = 0;
    if (sys->data.sys.alloc) md_arena_allocator_destroy(sys->data.sys.alloc);
    sys->data.sys.alloc = sys->data.state.alloc = nullptr;
}

bool system_drawn(const SystemContext& sys) {
    return sys.data.sys.atom.count > 0 && md_system_state_has_coords(&sys.data.state);
}

bool system_loaded(const SystemContext& sys) {
    return sys.data.sys.atom.count > 0 || file_list_count(&sys.files) > 0;
}

SystemContext* system_add(ApplicationState* app) {
    MD_ASSERT(app);
    SystemContext* first = &app->system_context;
    if (!system_loaded(*first)) {
        // Listed last, as an added one is
        size_t i = 0;
        while (i < app->num_systems && app->systems[i] != first) ++i;
        for (; i + 1 < app->num_systems; ++i) app->systems[i] = app->systems[i + 1];
        app->systems[app->num_systems - 1] = first;
        // Of whatever was loaded into it before, nothing refers to it by its id now, nor by its directory
        if (app->num_systems > 1) {
            first->id = app->next_system_id++;
            first->dir[0] = '\0';
            system_dir_number(*first);
        }
        return first;
    }
    return system_alloc(app);
}

SystemContext* system_alloc(ApplicationState* app) {
    MD_ASSERT(app);
    if (app->num_systems >= SYSTEM_MAX) return nullptr;
    SystemContext* sys = new SystemContext();
    system_context_init(sys, app);
    sys->id = app->next_system_id++;
    system_dir_number(*sys);
    app->systems[app->num_systems++] = sys;
    return sys;
}

void system_release(ApplicationState* app, SystemContext* sys) {
    MD_ASSERT(app && sys);
    if (sys == &app->system_context) return;
    size_t i = 0;
    while (i < app->num_systems && app->systems[i] != sys) ++i;
    if (i == app->num_systems) return;
    for (; i + 1 < app->num_systems; ++i) app->systems[i] = app->systems[i + 1];
    app->systems[--app->num_systems] = nullptr;
    system_context_destroy(sys);
    delete sys;
}

SystemContext* system_find(ApplicationState* app, uint32_t id) {
    MD_ASSERT(app);
    if (id == SYSTEM_ID_NONE) return nullptr;
    for (size_t i = 0; i < app->num_systems; ++i) {
        if (app->systems[i]->id == id) return app->systems[i];
    }
    return nullptr;
}

const SystemContext* system_find(const ApplicationState& app, uint32_t id) {
    return system_find(const_cast<ApplicationState*>(&app), id);
}

// --- THE PALETTE ---

void palette_builtin(Palette* p) {
    MD_ASSERT(p);
    for (int z = 0; z < MD_Z_Count; ++z) p->element[z] = md_float4_from_u32(md_atomic_number_cpk_color((md_atomic_number_t)z));
    for (int m = 0; m < MD_MONOMER_COUNT; ++m) p->monomer[m] = md_float4_from_u32(md_monomer_color((md_monomer_t)m));
    // As the representations had them: an unknown structure the darker grey
    p->ss_unknown = md_float4_from_u32(0xFF555555);
    p->ss_coil    = {0.86f, 0.86f, 0.86f, 1.0f};
    p->ss_helix   = {0.12f, 0.86f, 0.12f, 1.0f};
    p->ss_sheet   = {0.12f, 0.12f, 0.86f, 1.0f};
}

// Each entry of a palette by the name it is kept by, and where it is in it
struct PaletteEntry {
    char         name[16];
    md_float4_t* color;
};

static size_t palette_entries(PaletteEntry* out, size_t cap, Palette* p) {
    size_t n = 0;
    auto add = [&](md_str_t prefix, md_str_t id, md_float4_t* c) {
        if (n >= cap || id.len == 0) return;
        snprintf(out[n].name, sizeof(out[n].name), MD_STR_FMT ":" MD_STR_FMT, MD_STR_ARG(prefix), MD_STR_ARG(id));
        out[n].color = c;
        n += 1;
    };
    for (int z = 1; z < MD_Z_Count; ++z) add(MD_STR_LIT("e"), md_atomic_number_symbol((md_atomic_number_t)z), &p->element[z]);
    for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) add(MD_STR_LIT("m"), md_monomer_id((md_monomer_t)m), &p->monomer[m]);
    add(MD_STR_LIT("ss"), MD_STR_LIT("unknown"), &p->ss_unknown);
    add(MD_STR_LIT("ss"), MD_STR_LIT("coil"),    &p->ss_coil);
    add(MD_STR_LIT("ss"), MD_STR_LIT("helix"),   &p->ss_helix);
    add(MD_STR_LIT("ss"), MD_STR_LIT("sheet"),   &p->ss_sheet);
    return n;
}

#define PALETTE_MAX_ENTRIES (MD_Z_Count + MD_MONOMER_COUNT + 4)

void palette_write(char* buf, size_t cap, const Palette& palette) {
    MD_ASSERT(buf && cap > 0);
    Palette builtin;
    palette_builtin(&builtin);
    PaletteEntry cur[PALETTE_MAX_ENTRIES], ref[PALETTE_MAX_ENTRIES];
    const size_t n = palette_entries(cur, PALETTE_MAX_ENTRIES, const_cast<Palette*>(&palette));
    palette_entries(ref, PALETTE_MAX_ENTRIES, &builtin);
    size_t len = 0;
    buf[0] = '\0';
    for (size_t i = 0; i < n; ++i) {
        const uint32_t c = md_u32_from_float4(*cur[i].color);
        if (c == md_u32_from_float4(*ref[i].color)) continue;
        // RGBA, as written: the bytes from the lowest
        const int w = snprintf(buf + len, cap - len, "%s%s=%02x%02x%02x%02x", len ? " " : "", cur[i].name,
                               c & 0xFF, (c >> 8) & 0xFF, (c >> 16) & 0xFF, (c >> 24) & 0xFF);
        if (w < 0 || (size_t)w >= cap - len) {
            VIAMD_LOG_ERROR("Palette: more of it changed than the settings keep; the rest is not kept");
            buf[len] = '\0';
            return;
        }
        len += (size_t)w;
    }
}

static int hex_digit(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

void palette_read(Palette* palette, md_str_t text) {
    MD_ASSERT(palette);
    palette_builtin(palette);
    PaletteEntry entries[PALETTE_MAX_ENTRIES];
    const size_t n = palette_entries(entries, PALETTE_MAX_ENTRIES, palette);
    while (text.len > 0) {
        text = md_str_trim_beg(text);
        size_t end = 0;
        while (end < text.len && text.ptr[end] != ' ') ++end;
        const md_str_t item = md_str_substr(text, 0, end);
        text = md_str_substr(text, end, SIZE_MAX);
        md_str_t name, value;
        if (!md_str_cut(&name, &value, item, md_str_find_char(item, '=')) || value.len != 8) continue;
        uint32_t c = 0;
        bool ok = true;
        for (int k = 0; k < 4 && ok; ++k) {
            const int hi = hex_digit(value.ptr[2 * k]), lo = hex_digit(value.ptr[2 * k + 1]);
            ok = hi >= 0 && lo >= 0;
            c |= (uint32_t)((hi << 4) | lo) << (8 * k);
        }
        if (!ok) continue;
        for (size_t i = 0; i < n; ++i) {
            if (md_str_eq_cstr(name, entries[i].name)) {
                *entries[i].color = md_float4_from_u32(c);
                break;
            }
        }
    }
}

void palette_changed(ApplicationState* app) {
    MD_ASSERT(app);
    app->palette_version += 1;
    palette_write(app->settings.palette, sizeof(app->settings.palette), app->palette);
    app_settings::mark_dirty();
}

// --- SYSTEM DIRECTORIES ---

static bool dir_char_valid(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
}

static bool dir_is_number(md_str_t dir) {
    if (dir.len == 0) return false;
    for (size_t i = 0; i < dir.len; ++i) {
        if (dir.ptr[i] < '0' || dir.ptr[i] > '9') return false;
    }
    return true;
}

// The system other than sys with the directory, NULL when none has it
static SystemContext* dir_holder(const ApplicationState& app, const SystemContext& sys, md_str_t dir) {
    for (size_t i = 0; i < app.num_systems; ++i) {
        SystemContext* sc = app.systems[i];
        if (sc != &sys && md_str_eq_cstr(dir, sc->dir)) return sc;
    }
    return nullptr;
}

void system_dir_number(SystemContext& sys) {
    MD_ASSERT(sys.app);
    char buf[sizeof(sys.dir)];
    for (uint32_t n = 0;; ++n) {
        const int len = snprintf(buf, sizeof(buf), "%u", n);
        if (!dir_holder(*sys.app, sys, md_str_t{buf, (size_t)len})) break;
    }
    md_memcpy(sys.dir, buf, sizeof(buf));
}

const char* system_dir_invalid(const ApplicationState& app, const SystemContext& sys, md_str_t dir) {
    if (dir.len == 0) return "it is empty";
    if (dir.len >= sizeof(sys.dir)) return "it is too long";
    for (size_t i = 0; i < dir.len; ++i) {
        if (!dir_char_valid(dir.ptr[i])) return "only letters, digits, '_', '-' and '.' go in it";
    }
    if (dir_is_number(dir)) return "the numbers are given to the systems";
    if (dir_holder(app, sys, dir)) return "another system has it";
    return nullptr;
}

bool system_set_dir(SystemContext& sys, md_str_t dir) {
    MD_ASSERT(sys.app);
    if (system_dir_invalid(*sys.app, sys, dir)) return false;
    md_str_copy_to_char_buf(sys.dir, sizeof(sys.dir), dir);
    return true;
}

void system_dir_take(SystemContext& sys, md_str_t dir) {
    MD_ASSERT(sys.app);
    if (dir.len == 0 || dir.len >= sizeof(sys.dir)) return;
    for (size_t i = 0; i < dir.len; ++i) {
        if (!dir_char_valid(dir.ptr[i])) return;
    }
    if (SystemContext* other = dir_holder(*sys.app, sys, dir)) {
        if (!dir_is_number(md_str_from_cstr(other->dir))) return;
        md_str_copy_to_char_buf(sys.dir, sizeof(sys.dir), dir);
        system_dir_number(*other);
        return;
    }
    md_str_copy_to_char_buf(sys.dir, sizeof(sys.dir), dir);
}

size_t system_count(const ApplicationState& app) {
    size_t n = 0;
    for (size_t i = 0; i < app.num_systems; ++i) n += system_loaded(*app.systems[i]);
    return n;
}

SystemContext* picking_hit_system(ApplicationState* app, const PickingHit& hit) {
    return system_find(app, hit.system);
}

const SystemContext* picking_hit_system(const ApplicationState& app, const PickingHit& hit) {
    return system_find(app, hit.system);
}

// Opens the run in path for the loaded system, beside those it has, and shows it. False when the file holds no frames to
// make a run of, which is not an error for a structure file, or the system has as many runs as it can.
static bool open_run(ApplicationState* state, SystemContext& sys, md_str_t path, const loader::LoaderState& load_state) {
    if (sys.num_runs >= RUN_MAX) {
        VIAMD_LOG_ERROR("'" MD_STR_FMT "' is not opened: at most %zu runs of a system are opened at once", MD_STR_ARG(path), RUN_MAX);
        return false;
    }
    interrupt_async_tasks(state);
    // A system without coordinates (a topology alone) takes those of the first frame as its own: room for them
    const bool coordinate_less = !md_system_state_has_coords(&sys.data.state);
    if (coordinate_less) {
        const md_unitcell_t cell = sys.data.state.unitcell;
        md_system_state_init(&sys.data.state, sys.data.sys.atom.count);
        sys.data.state.unitcell = cell;
    }

    // Publishing the run IS opening the trajectory: the run is named after the file, and its positions are read from it
    // frame by frame from here on. Everything sampled along it lives below the same prefix - the positions, the per
    // frame backbone data, an energy file loaded later - so that closing it is removing the prefix.
    SystemRun* run = new SystemRun();
    run_path_from_file(run->path, sizeof(run->path), path);
    run_path_make_unique(run->path, sizeof(run->path), sys.data.sys.attributes);
    const uint32_t flags = (load_state.flags & LoaderFlag_DisableCacheWrite) ? MD_RUN_FLAG_DISABLE_CACHE_WRITE : 0;
    if (!loader::publish_run(&sys.data.sys, path, md_str_from_cstr(run->path), flags) || run_num_frames(sys, run) == 0) {
        md_attributes_remove_prefix(&sys.data.sys.attributes, md_str_from_cstr(run->path));
        delete run;
        return false;
    }
    const size_t num_frames = run_num_frames(sys, run);
    const bool first_run = sys.num_runs == 0;
    sys.runs[sys.num_runs++] = run;
    file_list_record(&sys.files, file_entry_from(path, load_state), FileRole_Trajectory, state->allocator.persistent);
    md_array_last(sys.files.runs)->frame_count = num_frames;
    trajectory_backbone_start(sys, *run);

    // Its first frame is the system's own coordinates when it had none: what it is inferred of
    // (system_inference_defer), and what its instances without a run show
    if (coordinate_less) {
        extract_frame(sys, run, 0, &sys.data.state);
        system_instances_update(sys, MolBit_ClearVelocity);
    }

    // Shown by the first instance when it is the system's first run, else by an instance of its own beside the others,
    // drawn as the first is - unless what opens it places the runs in the instances itself (instances_held)
    if (first_run) {
        instance_show(system_instance(sys), 0);
    } else if (!sys.instances_held) {
        if (SystemInstance* inst = instance_add(sys)) {
            instance_show(*inst, (int)sys.num_runs - 1);
            representations_copy_to_instance(sys, system_instance(sys), *inst);
        } else {
            VIAMD_LOG_INFO("'" MD_STR_FMT "' is opened as a run of '%s', which has as many instances as it can show (%zu): it is shown in none", MD_STR_ARG(path), sys.name, INSTANCE_MAX);
        }
    }
    return true;
}


// ### INFERENCE IN THE BACKGROUND ###
//
// The work (md_system_infer_prepare) shares the core of the system, which nothing writes until it is published or
// discarded, and has of its own a derived arena (its pages on the heap, the system's derived backing), with the
// coordinates it is inferred from, which become the system's reference state. It is inferred on a worker with its loops
// on the pool. Published, its arena becomes the system's and the one it replaces is destroyed; the system keeps its
// core, its attribute table and its force field. What is left of the work in its arena (its copies of the types of the
// atoms) goes with it. The work of a topology which changed meanwhile (a supplement, an edit) is dropped and made again,
// from the same coordinates.

static void inference_discard_work(SystemInference& inf) {
    if (inf.arena) {
        md_system_infer_discard(&inf.work);
        md_arena_allocator_destroy(inf.arena);
        inf.arena = nullptr;
    }
    inf.work = {};
}

static void system_inference_on_main(SystemContext& sys, task_system::Status status);

// What the system is inferred from when nothing says: the coordinates of the work under way, should it be made again
// (those Recompute Bonds gave it), else those the system was inferred from (an edit of its bonds keeps the bonds inferred
// as they were), else those shown (as loaded)
static const md_system_state_t* inference_coordinates(const SystemContext& sc) {
    const SystemInference& inf = sc.data.inference;
    if (inf.arena && inf.work.sys.reference.xyz) return &inf.work.sys.reference;
    if (sc.data.sys.reference.xyz) return &sc.data.sys.reference;
    return &sc.data.state;
}

void system_inference_start(SystemContext& sc, const md_system_state_t* state) {
    SystemInference& inf = sc.data.inference;
    inf.job.cancel_and_wait();
    const md_system_t& sys = sc.data.sys;
    if (sys.atom.count == 0) {
        inference_discard_work(inf);
        inf.pending = false;
        return;
    }

    // Made before the work it replaces is dropped: the coordinates may be that work's
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_infer_work_t work = {};
    const bool prepared = md_system_infer_prepare(&work, &sys, state ? state : inference_coordinates(sc), arena);
    inference_discard_work(inf);
    if (!prepared) {
        VIAMD_LOG_ERROR("Could not prepare the inference of the structure of the system");
        md_arena_allocator_destroy(arena);
        inf.pending = false;
        return;
    }
    inf.arena = arena;
    inf.work = work;
    inf.done = false;
    inf.pending = true;
    inf.restart = false;

    SystemInference* p = &inf;
    inf.job.start(MD_STR_LIT("Inferring structure"),
        [p]() {
            md_system_infer_run(&p->work, task_system::executor());
            // Bond inference leaves its scratch committed in this thread's arenas (it rewinds them, which keeps the
            // pages): given back, as this was large and is done
            md_temp_trim();
            if (!task_system::cancelled()) p->done = true;
        },
        [s = &sc](task_system::Status status) { system_inference_on_main(*s, status); });
}

bool system_inference_pending(const SystemContext& sys) {
    return sys.data.inference.pending;
}

// What the system shown is built of, made again of the inferred one: the molecule on the GPU (its bonds, backbones and
// chain gaps), the glyphs of the cartoon, the secondary structure of the coordinates shown, the representations of
// them, the backbone angles and secondary structure of the runs' frames, and what the components derived from the
// system (ViamdSystemInit). However inference ended: published, or done on the main thread (system_inference_finish).
static void system_inferred_refresh(SystemContext& sys) {
    ApplicationState* app = sys.app;
    // The runs' backbones, which there were none of before
    for (size_t i = 0; i < sys.num_runs; ++i) trajectory_backbone_start(sys, *sys.runs[i]);
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        md_gl_rep_destroy(sys.representation.reps[i].md_rep);
        sys.representation.reps[i].md_rep = {};
    }
    for (size_t k = 0; k < sys.num_instances; ++k) {
        SystemInstance& inst = *sys.instances[k];
        assembly_subsets_free(inst);
        cartoon_glyphs_free(&inst.cartoon_glyphs);
        // Made again in place: the handle stays, so what components made of the molecule when the system was loaded
        // (the representations of the QM windows) stays valid without their knowing it was. Destroyed and created, it
        // had a new id, and those were left drawing a molecule which no longer existed. The atoms are those it was made
        // of; should they not be, it is a new molecule after all.
        if (!md_gl_mol_rebuild(inst.gl_mol, &sys.data.sys)) {
            md_gl_mol_destroy(inst.gl_mol);
            inst.gl_mol = md_gl_mol_create(&sys.data.sys);
        }
        cartoon_glyphs_init(&inst.cartoon_glyphs, sys.data.sys, inst.gl_mol, app->allocator.persistent);
        // Its backbone, of the backbones inferred, when its coordinates are made again
        inst.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_DirtyRadius | MolBit_DirtyFlags | MolBit_DirtyBonds | MolBit_ClearVelocity;
        inst.interpolate_system_state = true;
        inst.last_interpolated_nearest_frame = -1;
        inst.gpu_buffers_version += 1;
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemInit, viamd::EventPayloadType_System, &sys);
    init_all_representations(sys);
    flag_all_representations_as_dirty(sys);
    app->script.compile_ir = true;
    app->script.time_since_last_change = 1.0e9;
}

// The default representations, decided from what is now known of the system, in place of the placeholder shown while it
// was inferred, unless it was changed or others were made meanwhile. However inference ended: published, or done again
// on the main thread (system_inference_finish, as when a supplement or a run comes along before it is done).
static void inference_defaults(SystemContext& sys) {
    SystemInference& inf = sys.data.inference;
    if (!inf.defaults_pending) return;
    inf.defaults_pending = false;
    const size_t num_reps = md_array_size(sys.representation.reps);
    const Representation* rep = num_reps == 1 ? &sys.representation.reps[0] : nullptr;
    const bool untouched = rep && strcmp(rep->name, inf.placeholder_name) == 0 && rep->type == RepresentationType::SpaceFill &&
        rep->color_mapping == ColorMapping::Type && rep->enabled && rep->filter.text[0] == '\0' && rep->scale.x == 1.0f;
    if (untouched) {
        remove_all_representations(sys);
        create_default_representations(sys);
        recompute_atom_visibility_mask(sys);
    }
}

// The work published: the main thread, nothing else running which reads the system. False when the topology changed
// since the work was made of it, and then nothing is done.
static bool inference_publish(SystemContext& sc) {
    SystemInference& inf = sc.data.inference;
    md_system_t& sys = sc.data.sys;

    // A version of its own: what was derived from the core is stale
    if (!md_system_infer_publish(&sys, &inf.work)) return false;

    // Its derived arena is the system's now, the one it replaced destroyed: the arena has nothing of the system
    md_arena_allocator_destroy(inf.arena);
    inf.arena = nullptr;
    inf.work = {};
    inf.pending = false;
    inf.restart = false;

    system_inferred_refresh(sc);

    inference_defaults(sc);
    md_temp_trim();
    VIAMD_LOG_INFO("The structure of the system is inferred: %zu bonds, %zu molecules", sys.bond.count, md_system_molecule_count(&sys));
    return true;
}

static void system_inference_on_main(SystemContext& sys, task_system::Status status) {
    ApplicationState* app = sys.app;
    SystemInference& inf = sys.data.inference;
    if (!inf.pending) return;
    if (status == task_system::Status::Completed && inf.done && sys.data.sys.topology_version == inf.work.base_version) {
        script_evals_interrupt(app);
        task_system::wait_all();
        inference_publish(sys);
    } else {
        // Interrupted, or of a topology which changed meanwhile: made again of the system as it is, next frame
        inf.restart = true;
    }
}

void system_inference_update(SystemContext& sys) {
    SystemInference& inf = sys.data.inference;
    if (inf.pending && inf.restart && !inf.job.running()) {
        system_inference_start(sys, nullptr);
    }
    // The runs' backbones, whose computation was stopped by what changed another system (interrupt_async_tasks stops
    // the work of all of them), and of this one once it is inferred again (system_inferred_refresh starts them)
    if (!inf.pending) {
        for (size_t i = 0; i < sys.num_runs; ++i) {
            SystemRun& run = *sys.runs[i];
            if (run.backbone_computations.status() == task_system::Status::Cancelled) trajectory_backbone_start(sys, run);
        }
    }
}

void system_inference_finish(SystemContext& sys) {
    ApplicationState* app = sys.app;
    SystemInference& inf = sys.data.inference;
    if (!inf.pending) return;
    inf.job.wait();
    script_evals_interrupt(app);
    task_system::wait_all();
    if (!(inf.done && inference_publish(sys))) {
        // Inferred on the system itself, here and now, from the coordinates the work was to be inferred from (copied by
        // the inference before the work is dropped)
        inf.job.cancel_and_wait();
        md_system_infer(&sys.data.sys, inference_coordinates(sys), task_system::executor());
        inference_discard_work(inf);
        inf.pending = false;
        inf.restart = false;
        system_inferred_refresh(sys);
        inference_defaults(sys);
    }
}

// Inferred again once what is being loaded is in: the files of a drop, of a workspace, are each loaded before it is
// started, the next time the systems are updated (system_inference_update) - once, of all of them. What is under way
// stops, as it shares the core the files write. A workspace waits for it instead (system_inference_finish).
static void system_inference_defer(SystemContext& sc) {
    SystemInference& inf = sc.data.inference;
    inf.job.cancel_and_wait();
    inf.pending = sc.data.sys.atom.count > 0;
    inf.restart = inf.pending;
    inf.done = false;
}

static void system_inference_free(SystemContext& sys) {
    SystemInference& inf = sys.data.inference;
    inf.job.cancel_and_wait();
    inference_discard_work(inf);
    inf.pending = false;
    inf.restart = false;
    inf.defaults_pending = false;
}

// The user's edits of the bonds (md_system_t::bond_edit) given to the system again, of the same atoms: a workspace's, or
// those of the system before it was built again. The removals first: a bond the user made where one was removed is the
// user's, and stays.
static void apply_bond_edits(SystemContext& sc, const md_bond_list_t& edits) {
    const size_t num_atoms = sc.data.sys.atom.count;
    const size_t num_edits = md_bond_list_count(&edits);
    if (num_edits == 0) return;
    system_core_edit_begin(sc);
    for (int add = 0; add < 2; ++add) {
        for (size_t i = 0; i < num_edits; ++i) {
            const md_atom_pair_t pair = edits.pairs[i];
            const md_bond_flags_t flags = edits.flags[i];
            if (((flags & MD_BOND_FLAG_INACTIVE) == 0) != (add == 1)) continue;
            if ((size_t)pair.idx[0] >= num_atoms || (size_t)pair.idx[1] >= num_atoms || pair.idx[0] == pair.idx[1]) continue;
            if (add) md_system_bond_insert(&sc.data.sys, pair.idx[0], pair.idx[1], flags);
            else     md_system_bond_remove(&sc.data.sys, pair.idx[0], pair.idx[1], flags);
        }
    }
    system_inference_defer(sc);
}

// A system of its own joined to the atoms of sys (md_system_merge): its coordinates, its names, its bonds, ... as the
// kinds go (md_io.h, JOINING). Inferred again after, once what is loaded with it is in (system_inference_defer).
static LoadResult join_system(ApplicationState* state, SystemContext& sys, md_str_t path, const loader::LoaderState& load_state, const md_system_t* loaded, const md_system_state_t* loaded_state) {
    md_allocator_i* arena = nullptr;
    md_system_t file_sys = {};
    md_system_state_t file_state = {};
    if (!loaded) {
        arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
        file_sys.alloc = arena;
        file_state.alloc = arena;
        if (loader::load(&file_sys, &file_state, path, load_state) != MD_LOAD_OK) {
            VIAMD_LOG_ERROR("Failed to load molecular data from file '" MD_STR_FMT "'", MD_STR_ARG(path));
            md_system_free(&file_sys);
            md_arena_allocator_destroy(arena);
            md_temp_trim();
            return LoadResult_Failed;
        }
        loaded = &file_sys;
        loaded_state = &file_state;
    }

    // The core is edited: what reads it stops, and the system is inferred again of it once what is loaded is in
    system_core_edit_begin(sys);
    const bool had_coords = md_system_state_has_coords(&sys.data.state);
    const bool joined = md_system_merge(&sys.data.sys, &sys.data.state, loaded, loaded_state);
    if (arena) {
        md_system_free(&file_sys);
        md_arena_allocator_destroy(arena);
    }
    md_temp_trim();
    system_inference_defer(sys);
    if (!joined) {
        VIAMD_LOG_ERROR("'" MD_STR_FMT "' is not of the atoms of '%s'", MD_STR_ARG(path), sys.name);
        return LoadResult_Failed;
    }

    file_list_record(&sys.files, file_entry_from(path, load_state), FileRole_Supplement, state->allocator.persistent);
    system_instances_update(sys, MolBit_DirtyPosition | MolBit_DirtyRadius | MolBit_DirtyFlags | MolBit_DirtyBonds | MolBit_ClearVelocity);
    for (size_t k = 0; k < sys.num_instances; ++k) sys.instances[k]->last_interpolated_nearest_frame = -1;
    state->script.compile_ir = true;
    if (!had_coords && md_system_state_has_coords(&sys.data.state)) {
        // Shown now that it has coordinates: framed as a system loaded is
        reset_view_scene(&state->view.target, *state, state->view.camera.fov_y, state->view.region);
    }
    VIAMD_LOG_SUCCESS("Joined '" MD_STR_FMT "' to '%s'", MD_STR_ARG(path), sys.name);
    return LoadResult_Joined;
}

// A system of its own built onto sys without atoms: a structure, a topology alone (without coordinates), coordinates
// alone or a trajectory alone (of points)
static LoadResult build_system(ApplicationState* state, SystemContext& sys, md_str_t path_to_file, const loader::LoaderState& load_state, bool opens_own_run, const md_system_t* loaded, const md_system_state_t* loaded_state) {
    interrupt_async_tasks(state);
    free_trajectory_data(state, sys);
    free_system_data(state, sys);

    // sys.alloc and state.alloc are already the dataset arena - set at startup and put
    // back by free_system_data above. The loaded coordinates are the initial current state;
    // inference copies them into sys.reference, so the two stay independent.
    MD_ASSERT(sys.data.sys.alloc && sys.data.state.alloc == sys.data.sys.alloc);
    const bool ok = loaded ? md_system_merge(&sys.data.sys, &sys.data.state, loaded, loaded_state)
                           : loader::load(&sys.data.sys, &sys.data.state, path_to_file, load_state) == MD_LOAD_OK;
    // A loader's scratch is in this thread's temp arenas (its read buffer, what it holds per block), which keep
    // their pages once rewound: given back, as a load is large and is done
    md_temp_trim();
    if (!ok || sys.data.sys.atom.count == 0) {
        VIAMD_LOG_ERROR("Failed to load molecular data from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));
        return LoadResult_Failed;
    }
    VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));

    // Named after its file, until the user names it
    {
        const md_str_t stem = md_path_stem(path_to_file);
        md_str_copy_to_char_buf(sys.name, sizeof(sys.name), stem);
    }

    // In the list before the system is initialized, so that what handles that sees what it was built from. A trajectory
    // alone is the system's and its run's in one: it is recorded as the run it opens.
    const bool trajectory = (load_state.flags & LoaderFlag_Trajectory) && !(load_state.flags & LoaderFlag_System);
    if (!trajectory) file_list_record(&sys.files, file_entry_from(path_to_file, load_state), FileRole_System, state->allocator.persistent);
    sys.files.atom_count = sys.data.sys.atom.count;

    // Inferred whatever the data is, coarse grained or not: mdlib tells beads from atoms (and Martini from
    // other models) itself, and infers no bond between atoms a topology gives all the bonds of
    // (md_atom_data_t::bonds_complete). The colours of the particle types now, which the core is shown with
    // (they are of the core); the rest in the background once what is loaded with it is in, put in place when it is
    // done (system_inference_defer)
    md_system_infer_type_colors(&sys.data.sys);
    sys.data.inference.pending = true;
    init_system_data(state, sys);

    // A structure file of several frames (a multi model PDB, an XYZ trajectory) is its own run as well, and so is a
    // trajectory alone. Here rather than by the caller, so that the run is there for what handles ViamdLoadData.
    if (opens_own_run || trajectory) {
        open_run(state, sys, path_to_file, load_state);
    }
    system_inference_defer(sys);
    return LoadResult_System;
}

LoadResult load_data(ApplicationState* state, SystemContext& sys, md_str_t filepath, const loader::LoaderState& load_state, FileRole role, bool opens_own_run,
                     const md_system_t* loaded, const md_system_state_t* loaded_state, int along) {
    MD_ASSERT(state);

    LoadResult result = LoadResult_Failed;
    md_str_t path_to_file = md_path_make_canonical(filepath, state->allocator.frame);
    if (path_to_file) {
        // What the file is to the system: told by what each is, the role telling a series apart, and the frames of a
        // structure of several (a multi-model PDB) where it opens the run of its own system (file_list.h)
        const LoaderFlags flags = load_state.flags;
        const bool has_atoms  = sys.data.sys.atom.count > 0;
        const bool trajectory = role == FileRole_Trajectory || ((flags & LoaderFlag_Trajectory) && !(flags & LoaderFlag_System));
        const bool own_system = (flags & LoaderFlag_System) && !trajectory;
        const bool series     = role == FileRole_Series || (flags & LoaderFlag_Temporal);
        // Data of the system's atoms: a file which is no system of its own, or one which is that as well (an itp's
        // molecules, the QM data of a region of them) when its check says it is of these atoms
        const bool adds = has_atoms && !series && !trajectory &&
                          (!own_system || ((flags & LoaderFlag_Supplemental) && loader::check(&sys.data.sys, path_to_file, load_state) == MD_LOAD_OK));

        if (series || adds) {
            const FileRole data_role = series ? FileRole_Series : FileRole_Supplement;
            // A series is along the run asked for, else along the one shown
            const size_t run_index = (along >= 0 && (size_t)along < sys.num_runs) ? (size_t)along : (size_t)MD_MAX(system_instance(sys).run, 0);
            const SystemRun* run_along = run_index < sys.num_runs ? sys.runs[run_index] : nullptr;
            const md_str_t run = run_along ? md_str_from_cstr(run_along->path) : md_str_t{};
            if (!has_atoms) {
                VIAMD_LOG_ERROR("'" MD_STR_FMT "' adds to a system; load the system first", MD_STR_ARG(path_to_file));
                return LoadResult_Failed;
            }
            if (series && md_str_empty(run)) {
                VIAMD_LOG_ERROR("'" MD_STR_FMT "' holds data along a trajectory; load the trajectory first", MD_STR_ARG(path_to_file));
                return LoadResult_Failed;
            }
            // The table is about to be written to, and worker threads read it: script evaluation
            // (attr() and its frames), the backbone task and playback's frame loads, all through
            // extraction contexts. The evaluations are interrupted, since the script is recompiled
            // below and they would start over anyway; the rest are left to FINISH rather than
            // interrupted - an interrupted backbone task leaves a just loaded trajectory without its
            // backbone data.
            script_evals_interrupt(state);
            // Loaded before the system is inferred again, which is then of it (system_inference_defer): an inference
            // under way stops, as it shares the core the file may write
            const bool inferring = system_inference_pending(sys);
            if (inferring) system_core_edit_begin(sys);
            task_system::wait_all();

            const size_t num_atoms = sys.data.sys.atom.count;
            const uint64_t topology_version = sys.data.sys.topology_version;
            if (loader::load(&sys.data.sys, &sys.data.state, path_to_file, load_state, run) == MD_LOAD_OK) {
                // The contract of a supplement: it adds to the atoms, and never changes how many there are.
                // What is kept by atom index (selections, user bonds, particle type edits, and all of it in a
                // workspace) would silently apply to other atoms. The loaders check their file against the
                // system before they change it; this is where one that does not would show.
                MD_ASSERT(sys.data.sys.atom.count == num_atoms);
                if (sys.data.sys.atom.count != num_atoms) {
                    VIAMD_LOG_ERROR("'" MD_STR_FMT "' changed the number of atoms from %zu to %zu; what is kept by atom index no longer applies", MD_STR_ARG(path_to_file), num_atoms, sys.data.sys.atom.count);
                }
                result = LoadResult_Supplement;
                FileEntry entry = file_entry_from(path_to_file, load_state);
                entry.run = series ? (int)run_index : -1;
                file_list_record(&sys.files, entry, data_role, state->allocator.persistent);

                // A topology replaces bonds (and with them structures), which the GPU holds a copy of: the core
                // changed, and the system is inferred again of it in the background (a load never infers)
                system_instances_flag(sys, MolBit_DirtyBonds);
                if (inferring || sys.data.sys.topology_version != topology_version) system_inference_defer(sys);
                // attr() in the script resolves against the table when it is compiled.
                state->script.compile_ir = true;
                if (data_role == FileRole_Series) {
                    const char* hint = (load_state.type == LoaderType_EDR) ? "edr/<term>" :
                                       (load_state.type == LoaderType_AMBER_ENERGY) ? "amber_energy/<term>" :
                                       (load_state.type == LoaderType_XVG) ? "xvg/<file>/<legend>" : "csv/<file>/<column>";
                    VIAMD_LOG_SUCCESS("Loaded '" MD_STR_FMT "' into '" MD_STR_FMT "'; plot it from System > Series, or read it in the script with attr(\"%s\")", MD_STR_ARG(path_to_file), MD_STR_ARG(run), hint);
                } else {
                    VIAMD_LOG_SUCCESS("Successfully loaded supplemental data from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));
                }
            } else {
                VIAMD_LOG_ERROR("Failed to load supplemental data from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));
                if (inferring) system_inference_defer(sys);
            }
        } else if (trajectory && has_atoms) {
            if (open_run(state, sys, path_to_file, load_state)) {
                result = LoadResult_Trajectory;
                VIAMD_LOG_SUCCESS("Successfully opened trajectory from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));
            } else {
                VIAMD_LOG_ERROR("Failed to open trajectory from file '" MD_STR_FMT "'", MD_STR_ARG(path_to_file));
            }
        } else if (own_system && has_atoms) {
            result = join_system(state, sys, path_to_file, load_state, loaded, loaded_state);
        } else if (own_system || trajectory) {
            result = build_system(state, sys, path_to_file, load_state, opens_own_run, loaded, loaded_state);
        } else {
            VIAMD_LOG_ERROR("'" MD_STR_FMT "' adds to a system; load the system first", MD_STR_ARG(path_to_file));
            return LoadResult_Failed;
        }
        if (result == LoadResult_Failed) return result;

#if MD_ENABLE_GPU
        // The dataset's uploaded GTO basis, built here rather than by whichever component happens to
        // want it first: it is derived from the system's own basis/ attributes, every consumer of
        // that system wants the same one, and a system's GPU basis must not depend on which UI
        // component is compiled in. Returns false when the system publishes no basis, which is the
        // normal case rather than a failure.
        system_gpu_data_update(state, DEFAULT_GTO_CUTOFF_VALUE);
#endif

        // The frames selected were of the run there was
        if (load_state.flags & (LoaderFlag_System | LoaderFlag_Trajectory)) {
            frame_selection_clear(state);
            md_bitmask_clear(&state->selection.frame_highlight_mask);
        }

        LoadDataPayload data = {
            .app_state = state,
            .system = &sys,
            .loader_state = load_state,
            .path_to_file = path_to_file,
        };
        viamd::event_system_broadcast_event(viamd::EventType_ViamdLoadData, viamd::EventPayloadType_LoadData, &data);
    }

    return result;
}

// The next name of a list of them separated by commas, trimmed. False at the end.
static bool next_list_name(md_str_t* tok, md_str_t* list) {
    if (md_str_empty(*list)) return false;
    md_str_cut(tok, list, *list, md_str_find_char(*list, ','));
    *tok = md_str_trim(*tok);
    return true;
}

const char* bond_kind_str(md_bond_flags_t flags) {
    return md_bond_kind_name(md_bond_kind(flags));
}

bool bond_kind_parse(md_bond_flags_t* out, md_str_t str) {
    MD_ASSERT(out);
    str = md_str_trim(str);
    // 'interaction' was the word of a contact before it named every interaction
    if (md_str_eq_cstr(str, "interaction")) str = MD_STR_LIT("contact");
    for (uint32_t v = 0; v <= (uint32_t)MD_BOND_FLAG_KIND_MASK; v += 0x20u) {
        if (md_str_eq_cstr(str, md_bond_kind_name((md_bond_kind_t)v))) {
            *out = (md_bond_flags_t)v;
            return true;
        }
    }
    return false;
}

// Origins (a mask of 1 << md_bond_origin_t) by name, "file,inferred"
const char* bond_origins_str(char* buf, size_t cap, uint32_t origins) {
    MD_ASSERT(buf && cap > 0);
    buf[0] = '\0';
    size_t len = 0;
    for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED; ++o) {
        if (!(origins & (1u << o))) continue;
        const int n = snprintf(buf + len, cap - len, "%s%s", len ? "," : "", md_bond_origin_name((md_bond_origin_t)o));
        if (n < 0 || (size_t)n >= cap - len) break;
        len += (size_t)n;
    }
    return buf;
}

bool bond_origins_parse(uint32_t* out, md_str_t str) {
    MD_ASSERT(out);
    uint32_t mask = 0;
    md_str_t tok;
    while (next_list_name(&tok, &str)) {
        if (md_str_empty(tok)) continue;
        bool found = false;
        for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED && !found; ++o) {
            if (md_str_eq_cstr(tok, md_bond_origin_name((md_bond_origin_t)o))) {
                mask |= 1u << o;
                found = true;
            }
        }
        if (!found) return false;
    }
    *out = mask;
    return true;
}

// #workspace
//
// A workspace is everything that makes up one analysis, as opposed to the settings that follow the
// user between analyses (those are in the ImGui .ini, see app_settings): which files are open, how
// they are shown, the script, the selections, the plots and the windows showing them.
//
// LOADING goes in three steps, and the order is the design:
//   1. Reset. Everything a workspace holds goes back to its default first - here, and in every
//      component through EventType_ViamdDeserializeBegin - so a value the file does not mention
//      ends up as the default, never as whatever the previous workspace left behind.
//   2. Read. Every section is parsed. What does not need the data is applied as it is read; what
//      does (atom indices, frames, masks) is kept aside.
//   3. Load and apply. The files are loaded, then what was kept aside is applied against them, and
//      EventType_ViamdDeserializeEnd tells the components the data they asked for is there.
//
// SAVING builds the whole text first and only then touches the file, so a failure half way through
// leaves the previous workspace where it was.
//
// The files are stored per system and run, with the ids every other section refers to them by (file_list.h).
//
// VERSION, in [Workspace]: 1 (no [Workspace] section) stored the files in slots ([Files] MoleculeFile,
// TrajectoryFile, EnergyFile, SeriesFile); 2 stores systems, their runs and their files ([System], [Run],
// [Input]); 3 keys what is within an object by its path (ElectronicStructure.Field.Range, Slice.0.Opacity), see
// serialization_utils.h; 4 has several systems, and what is of one names it by its id (System=, the first entry of
// [Representation], [RepresentationGroup], [SystemScene], [UserBonds], [RemovedBonds], before each Mask of
// [ActiveSelection]): a section that names none is the first's. The operations, the periodic images and the assembly
// are the system's ([SystemScene]), no longer in [Operations] and [RenderSettings], which are still read. A version
// newer than the reader's is read for what the reader knows: a version 2 reader keeps a representation's name, filter,
// type and colours, and leaves what moved into blocks at its defaults. 5 has several instances of a system ([Instance]:
// the run each shows, its placement, periodic images and assembly, which were the system's in [SystemScene]), and the
// representations and groups name theirs (Instance=); one which names none is of the first.
static const int WORKSPACE_VERSION = 5;

struct WorkspaceWindow {
    char  name[32];
    bool* show;
};
static WorkspaceWindow workspace_windows[64];
static size_t num_workspace_windows = 0;

bool* workspace_window_flag(const char* name) {
    for (size_t i = 0; i < num_workspace_windows; ++i) {
        if (strcmp(workspace_windows[i].name, name) == 0) return workspace_windows[i].show;
    }
    return nullptr;
}

void workspace_register_window(const char* name, bool* show) {
    MD_ASSERT(name && show);
    for (size_t i = 0; i < num_workspace_windows; ++i) {
        if (strcmp(workspace_windows[i].name, name) == 0) {
            workspace_windows[i].show = show;
            return;
        }
    }
    if (num_workspace_windows < MD_COUNTOF(workspace_windows)) {
        WorkspaceWindow& w = workspace_windows[num_workspace_windows++];
        snprintf(w.name, sizeof(w.name), "%s", name);
        w.show = show;
    } else {
        MD_LOG_ERROR("Workspace: too many windows registered, '%s' is not stored", name);
    }
}

// What is read before the data it refers to is loaded, applied once it is
struct WorkspacePending {
    // The systems, their runs and the files of each ([System], [Run], [Input], or a [Files] section of a
    // workspace written before version 2)
    WorkspaceFiles files;
    int version;

    // The time shown, on the clock in its unit ([Animation] Time, TimeUnit); a workspace written before the clock has
    // the frame of the first system's run (Frame) only
    bool      has_frame;
    double    frame;
    bool      has_time;
    double    time;
    md_unit_t time_unit;
    // Whether the first system's run has the frames it had when the workspace was saved: what is kept by frame (the
    // frame, a view range in frames) is of other frames otherwise
    bool      frames_as_saved;

    // Kept aside rather than written to the camera: loading the molecule puts the camera at the
    // default view of what was loaded (init_system_data), which would overwrite it
    bool          has_camera;
    ViewTransform camera;

    // The systems it loads, by the ids it names them by ([System] Id), and what refers to the atoms of each: kept until
    // its files are loaded. In the order the workspace names them, which the Scene window lists them in; the first is
    // loaded into the first system held (system_context), the others into systems of their own.
    struct System {
        int            id;
        SystemContext* sys;

        // The user's edits of the bonds (md_system_t::bond_edit): a bond added, of the kind given if one was, or removed
        // (MD_BOND_FLAG_INACTIVE), of the kind which tells its range and the origin it had
        md_bond_list_t bond_edits;

        bool has_selection_mask;
        md_bitset_t selection_mask;

        // The user's groups ([Selection], as the stored selections they were): of the system, which loading its files
        // makes
        struct UserGroup {
            md_str_t    name;
            md_bitset_t atoms;
        };
        md_array(UserGroup) user_groups;

        bool has_recenter_target;
        md_bitset_t recenter_target;

        // The id of the assembly shown (system.assembly), of a workspace written before the instances: the first's
        char assembly_id[8];

        // Its instances ([Instance]), made as they are read, the first being the system's first: what of them refers to
        // its files (the run each shows, the assembly) is given once they are loaded
        struct Instance {
            int     run;
            int64_t frame;
            char    assembly_id[8];
        } instances[INSTANCE_MAX];
        size_t num_instances;
    } systems[SYSTEM_MAX];
    size_t num_systems;
    bool   systems_mapped;      // Every system the files named before the first section that refers to one
    size_t skipped_systems;     // Named, but beyond SYSTEM_MAX
    ApplicationState* app;

    // Timeline filter and zoom, in frames: a time would depend on the unit it is shown in. The
    // filter is how a workspace from before the timeline's periods held it, and becomes one of them:
    // what was shown over the filter is shown over that period.
    bool   filter_enabled;
    bool   has_filter_range;
    double filter_beg_frame;
    double filter_end_frame;
    bool   temporal_window_enabled;
    double temporal_window_extent;
    bool   has_periods;
    // The ends of each period as read, by its id: times in a unit, or frames of the first system's run (offsets in its
    // frames when it follows) of a workspace from before the clock, and of the filter made a period. Taken into the
    // clock's time once every system is loaded.
    struct PeriodEnds {
        uint32_t  id;
        bool      in_frames;
        md_unit_t unit;
    };
    PeriodEnds period_ends[TIMELINE_MAX_PERIODS];
    int        num_period_ends;
    // The timeline's view, ViewTime: in the run's time and its unit. ViewFrames, as workspaces before it wrote it, in frames.
    bool   has_view_range;
    bool   view_range_in_frames;
    double view_range[2];

    // The frames selected, of the first system's run
    bool has_frame_mask;
    md_bitset_t frame_mask;

    // The simulation box of a workspace from before it was a representation ([RenderSettings]), which becomes one
    bool   box_enabled;
    md_float4_t box_color;
};

static void workspace_reset(ApplicationState* data) {
    data->editor.SetText("");
    data->script.other_text.clear();
    data->files.workspace[0] = '\0';

    // The workspace's files are loaded onto nothing, as they were when they were first loaded: what a file
    // is to the system can depend on what is there (a VeloxChem file of the loaded atoms supplements them,
    // and builds the system otherwise). A workspace that names no files leaves nothing loaded. Every system goes:
    // the first is emptied, which the workspace's first is loaded into, the others are taken out.
    interrupt_async_tasks(data);
    for (size_t i = data->num_systems; i-- > 0;) {
        SystemContext* sc = data->systems[i];
        system_remove(*sc);
        system_release(data, sc);
    }
    // The first's directory, as the workspace gives it ([SystemScene] Dir=), else the first number
    data->system_context.dir[0] = '\0';
    system_dir_number(data->system_context);

    // The fields rather than the struct: whether the window is open is decided by [Windows]. No unit to the time yet:
    // what the workspace holds on the clock is taken into the clock's unit once its runs are loaded, not converted from
    // the unit of what was loaded before (clock_update).
    data->animation.time = 0.0;
    data->animation.time_unit = md_unit_none();
    data->animation.fps = 10.0f;
    data->animation.tension = 0.0f;
    data->animation.stride = 1;
    data->animation.interpolation = InterpolationMode::CubicSpline;
    data->animation.mode = PlaybackMode::Stopped;
    data->animation.loop = PlaybackLoop_Once;
    data->animation.loop_period = 0;
    data->animation.direction = 1;

    plot_clear(data->timeline.subplots, PLOT_MAX_SUBPLOTS);
    plot_clear(data->distributions.subplots, PLOT_MAX_SUBPLOTS);
    periodic_views_clear();
    data->timeline.num_subplots = 1;
    data->distributions.num_subplots = 1;
    data->timeline.periods = TimelinePeriods{};

    // The look to the factory's, which a workspace written before one of its settings existed was seen in. The rest
    // of the visuals is the computer's (the .ini), which no workspace sets.
    {
        const decltype(data->visuals) factory = {};
        data->visuals.background  = factory.background;
        data->visuals.shading     = factory.shading;
        data->visuals.tonemapping = factory.tonemapping;
        data->visuals.dof         = factory.dof;
    }
    data->view.mode = CameraMode::Perspective;
    data->view.camera.fov_y = Camera{}.fov_y;

    data->selection.granularity = SelectionGranularity::Atom;
    single_selection_sequence_clear(&data->selection.single_selection_sequence);
}

// A mask by the atoms it holds. Base64 now; a workspace written before that holds the raw bytes,
// which is read when it survived the text format (a mask whose bytes held no line break).
static bool deserialize_mask(md_bitset_t* mask, md_str_t arg, md_allocator_i* alloc) {
    if (md_str_begins_with(arg, MD_STR_LIT("###"))) {
        return viamd::extract_bitset(mask, arg, alloc);
    }
    md_str_t raw;
    viamd::extract_str(raw, arg);
    return !md_str_empty(raw) && md_bitset_deserialize(mask, raw.ptr, raw.len, alloc);
}

// The system of the workspace with the id, made when first asked for; -1 for the first. NULL for one it does not name, or
// one beyond as many as are held. A workspace that names none (no files) has the first, which its sections are of.
static WorkspacePending::System* workspace_system(WorkspacePending* p, int id) {
    ApplicationState* app = p->app;
    const size_t named = md_array_size(p->files.systems);
    auto add = [p, app](int ws_id) -> WorkspacePending::System* {
        if (p->num_systems == SYSTEM_MAX) {
            p->skipped_systems += 1;
            return nullptr;
        }
        SystemContext* sys = p->num_systems == 0 ? &app->system_context : system_alloc(app);
        if (!sys) {
            p->skipped_systems += 1;
            return nullptr;
        }
        WorkspacePending::System& s = p->systems[p->num_systems++];
        s = {};
        s.id = ws_id;
        s.sys = sys;
        s.selection_mask = md_bitset_empty_set();
        s.recenter_target = md_bitset_empty_set();
        return &s;
    };
    // The systems named so far, in their order, the first time: the [System] sections come before what refers to them
    if (!p->systems_mapped) {
        p->systems_mapped = true;
        for (size_t i = 0; i < named; ++i) add(p->files.systems[i].id);
        if (named == 0) add(0);
    }
    if (id < 0) return p->num_systems > 0 ? &p->systems[0] : nullptr;
    for (size_t i = 0; i < p->num_systems; ++i) {
        if (p->systems[i].id == id) return &p->systems[i];
    }
    // Named after what refers to it (a workspace edited by hand), or not at all
    for (size_t i = 0; i < named; ++i) {
        if (p->files.systems[i].id == id) return add(id);
    }
    return nullptr;
}

static SystemContext* workspace_system_fn(void* user, int id) {
    WorkspacePending::System* s = workspace_system((WorkspacePending*)user, id);
    return s ? s->sys : nullptr;
}

SystemContext* workspace_system_of(const viamd::deserialization_state_t& state, int id) {
    return state.reader ? workspace_system_fn(state.reader, id) : nullptr;
}

void load_workspace(ApplicationState* data, md_str_t filename) {
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_str_t txt = md_load_textfile(filename, temp_alloc);

    if (md_str_empty(txt)) {
        VIAMD_LOG_ERROR("Could not open workspace file: '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return;
    }

    viamd::deserialization_state_t state = {
        .filename = filename,
        .text = txt,
    };

    // ## 1. Reset
    workspace_reset(data);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserializeBegin, viamd::EventPayloadType_DeserializationState, &state);

    // ## 2. Read
    WorkspacePending pending = {};
    pending.app = data;
    pending.frame_mask = md_bitset_empty_set();
    state.reader = &pending;

    md_str_t folder = md_path_folder(filename);

    md_str_t section;
    while (viamd::next_section_header(section, state)) {
        md_str_t ident, arg;
        if (md_str_eq(section, MD_STR_LIT("Workspace"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Version"))) {
                    viamd::extract_int(pending.version, arg);
                }
            }
            if (pending.version > WORKSPACE_VERSION) {
                VIAMD_LOG_INFO("Workspace: written by a newer viamd (version %d, this one writes %d); what this one does not know is left out", pending.version, WORKSPACE_VERSION);
            }
        } else if (workspace_files_read(&pending.files, state, folder, temp_alloc)) {
            // The systems and their files
        } else if (md_str_eq(section, MD_STR_LIT("Animation"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Frame"))) {
                    pending.has_frame = viamd::extract_dbl(pending.frame, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Time"))) {
                    pending.has_time = viamd::extract_dbl(pending.time, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("TimeUnit"))) {
                    md_str_t u = {};
                    if (viamd::extract_str(u, arg)) md_unit_parse(&pending.time_unit, u);
                } else if (md_str_eq(ident, MD_STR_LIT("Fps"))) {
                    viamd::extract_flt(data->animation.fps, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Tension"))) {
                    viamd::extract_flt(data->animation.tension, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Stride"))) {
                    if (viamd::extract_int(data->animation.stride, arg)) {
                        data->animation.stride = MD_MAX(1, data->animation.stride);
                    }
                } else if (md_str_eq(ident, MD_STR_LIT("Interpolation"))) {
                    viamd::extract_enum(data->animation.interpolation, arg, (int)InterpolationMode::Count);
                } else if (md_str_eq(ident, MD_STR_LIT("Loop"))) {
                    viamd::extract_enum(data->animation.loop, arg, (int)PlaybackLoop_Count);
                } else if (md_str_eq(ident, MD_STR_LIT("LoopPeriod"))) {
                    int id = 0;
                    if (viamd::extract_int(id, arg)) data->animation.loop_period = (uint32_t)MD_MAX(id, 0);
                } else if (md_str_eq(ident, MD_STR_LIT("Transport"))) {
                    viamd::extract_enum(data->animation.transport, arg, (int)TransportAnchor_Count);
                }
            }
        } else if (md_str_eq(section, MD_STR_LIT("TimelineView"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("FilterEnabled"))) {
                    viamd::extract_bool(pending.filter_enabled, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("FilterFrames"))) {
                    float r[2];
                    if (viamd::extract_flt_vec(r, 2, arg)) {
                        pending.has_filter_range = true;
                        pending.filter_beg_frame = r[0];
                        pending.filter_end_frame = r[1];
                    }
                } else if (md_str_eq(ident, MD_STR_LIT("TemporalWindowEnabled"))) {
                    viamd::extract_bool(pending.temporal_window_enabled, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("TemporalWindowExtent"))) {
                    viamd::extract_dbl(pending.temporal_window_extent, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("PeriodLayout"))) {
                    viamd::extract_enum(data->timeline.periods.layout, arg, (int)TimelinePeriodLayout_Count);
                } else if (md_str_eq(ident, MD_STR_LIT("PeriodShadeAll"))) {
                    viamd::extract_bool(data->timeline.periods.shade_all, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("ViewTime"))) {
                    if (viamd::extract_dbl_vec(pending.view_range, 2, arg)) {
                        pending.has_view_range = true;
                        pending.view_range_in_frames = false;
                    }
                } else if (md_str_eq(ident, MD_STR_LIT("ViewFrames"))) {
                    if (viamd::extract_dbl_vec(pending.view_range, 2, arg)) {
                        pending.has_view_range = true;
                        pending.view_range_in_frames = true;
                    }
                }
            }
        } else if (md_str_eq(section, MD_STR_LIT("TimelinePeriod"))) {
            pending.has_periods = true;
            char unit[64];
            bool in_frames = false;
            const TimelinePeriod* r = timeline_periods_deserialize_period(state, &data->timeline.periods, unit, sizeof(unit), &in_frames);
            if (!r) {
                VIAMD_LOG_ERROR("Workspace: more than %d timeline periods, the rest are left out", TIMELINE_MAX_PERIODS);
            } else if (pending.num_period_ends < TIMELINE_MAX_PERIODS) {
                WorkspacePending::PeriodEnds& e = pending.period_ends[pending.num_period_ends++];
                e.id = r->id;
                e.in_frames = in_frames;
                e.unit = md_unit_none();
                if (unit[0]) md_unit_parse(&e.unit, md_str_from_cstr(unit));
            }
        } else if (md_str_eq(section, MD_STR_LIT("RenderSettings"))) {
            auto& v = data->visuals;
            while (viamd::next_entry(ident, arg, state)) {
                if      (md_str_eq(ident, MD_STR_LIT("BackgroundColor")))      viamd::extract_vec3(v.background.color, arg);
                // BackgroundIntensity scaled the background radiance up to the white point of the tone mappers.
                // The background now shows as picked, so it is ignored.
                else if (md_str_eq(ident, MD_STR_LIT("ShadingMode")))          viamd::extract_enum(v.shading.mode, arg, (int)postprocess_pipeline::ShadingMode_Count);
                else if (md_str_eq(ident, MD_STR_LIT("ShadingColors")))        viamd::extract_enum(v.shading.colors, arg, (int)postprocess_pipeline::ColorSpace_Count);
                else if (md_str_eq(ident, MD_STR_LIT("ShadingLight")))         viamd::extract_flt(v.shading.light, arg);
                else if (md_str_eq(ident, MD_STR_LIT("ShadingEnvironment")))   viamd::extract_flt(v.shading.environment, arg);
                else if (md_str_eq(ident, MD_STR_LIT("ShadingAmbient")))       viamd::extract_flt(v.shading.ambient, arg);
                else if (md_str_eq(ident, MD_STR_LIT("ShadingSpecular")))      viamd::extract_flt(v.shading.specular, arg);
                // The ambient occlusion and the anti-aliasing (Ssao*, Fxaa*, Taa*, MotionBlur*, Sharpen*) were written by
                // earlier versions and are ignored: they are the computer's settings (the .ini), not the figure's. Without
                // a radius of the scene, nothing of the ambient occlusion is the figure's. DofFocusScale belonged to the
                // previous scale-dependent DOF and is ignored too.
                else if (md_str_eq(ident, MD_STR_LIT("TonemapEnabled")))       viamd::extract_bool(v.tonemapping.enabled, arg);
                else if (md_str_eq(ident, MD_STR_LIT("Tonemapper")))           viamd::extract_enum(v.tonemapping.tonemapper, arg, (int)postprocess_pipeline::Tonemapper_Count);
                else if (md_str_eq(ident, MD_STR_LIT("TonemapExposure")))      viamd::extract_flt(v.tonemapping.exposure, arg);
                else if (md_str_eq(ident, MD_STR_LIT("TonemapGamma")))         viamd::extract_flt(v.tonemapping.gamma, arg);
                else if (md_str_eq(ident, MD_STR_LIT("DofEnabled")))           viamd::extract_bool(v.dof.enabled, arg);
                else if (md_str_eq(ident, MD_STR_LIT("DofAperture")))          viamd::extract_flt(v.dof.aperture, arg);
                // DEPRECATED here: the first system's, before they were each system's ([SystemScene])
                else if (md_str_eq(ident, MD_STR_LIT("PeriodicImages"))) {
                    float n[3];
                    WorkspacePending::System* ws = workspace_system(&pending, -1);
                    if (ws && viamd::extract_flt_vec(n, 3, arg)) {
                        for (int k = 0; k < 3; ++k) system_instance(*ws->sys).scene.periodic_images[k] = MD_CLAMP((int)n[k], 0, PERIODIC_IMAGES_MAX);
                    }
                }
                else if (md_str_eq(ident, MD_STR_LIT("Assembly"))) {
                    md_str_t id = {};
                    WorkspacePending::System* ws = workspace_system(&pending, -1);
                    if (ws && viamd::extract_str(id, arg)) md_str_copy_to_char_buf(ws->assembly_id, sizeof(ws->assembly_id), id);
                }
                else if (md_str_eq(ident, MD_STR_LIT("SimulationBoxEnabled"))) viamd::extract_bool(pending.box_enabled, arg);
                else if (md_str_eq(ident, MD_STR_LIT("SimulationBoxColor")))   viamd::extract_vec4(pending.box_color, arg);
            }
        } else if (md_str_eq(section, MD_STR_LIT("Camera"))) {
            // A transform is only restored whole: half of one mixed with the default view is neither
            bool has_position = false, has_orientation = false, has_distance = false;
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Position"))) {
                    has_position = viamd::extract_vec3(pending.camera.position, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Orientation")) || md_str_eq(ident, MD_STR_LIT("Rotation"))) {
                    // Rotation: DEPRECATED name
                    has_orientation = viamd::extract_quat(pending.camera.orientation, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Distance"))) {
                    has_distance = viamd::extract_flt(pending.camera.distance, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Mode"))) {
                    viamd::extract_enum(data->view.mode, arg, (int)CameraMode::Count);
                } else if (md_str_eq(ident, MD_STR_LIT("FovY"))) {
                    viamd::extract_flt(data->view.camera.fov_y, arg);
                }
            }
            pending.has_camera = has_position && has_orientation && has_distance;
        } else if (md_str_eq(section, MD_STR_LIT("Operations"))) {
            // Of a workspace written before [SystemScene]: the operations, all of them the one system's
            WorkspacePending::System* ws = workspace_system(&pending, -1);
            while (viamd::next_entry(ident, arg, state)) {
                if (!ws) continue;
                if (md_str_eq(ident, MD_STR_LIT("RecenterTarget"))) ws->has_recenter_target = deserialize_mask(&ws->recenter_target, arg, temp_alloc);
                else system_scene_read_entry(*ws->sys, ident, arg);
            }
        } else if (md_str_eq(section, MD_STR_LIT("SystemScene"))) {
            // Of the system it names first (System=): of one the workspace does not load, left out
            WorkspacePending::System* ws = workspace_system(&pending, -1);
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("System"))) {
                    int id = -1;
                    ws = viamd::extract_int(id, arg) ? workspace_system(&pending, id) : nullptr;
                } else if (!ws) {
                    continue;
                } else if (md_str_eq(ident, MD_STR_LIT("RecenterTarget"))) {
                    ws->has_recenter_target = deserialize_mask(&ws->recenter_target, arg, temp_alloc);
                } else if (md_str_eq(ident, MD_STR_LIT("Assembly"))) {
                    md_str_t id = {};
                    if (viamd::extract_str(id, arg)) md_str_copy_to_char_buf(ws->assembly_id, sizeof(ws->assembly_id), id);
                } else {
                    system_scene_read_entry(*ws->sys, ident, arg);
                }
            }
        } else if (md_str_eq(section, MD_STR_LIT("Instance"))) {
            // Of the system it names first (System=), made when it is given its id (Id=), the entries after it are its
            WorkspacePending::System* ws = workspace_system(&pending, -1);
            SystemInstance* inst = nullptr;
            WorkspacePending::System::Instance* pi = nullptr;
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("System"))) {
                    int id = -1;
                    ws = viamd::extract_int(id, arg) ? workspace_system(&pending, id) : nullptr;
                } else if (md_str_eq(ident, MD_STR_LIT("Id"))) {
                    int id = 0;
                    if (!ws || inst || !viamd::extract_int(id, arg) || id <= 0) continue;
                    SystemContext& sc = *ws->sys;
                    if (ws->num_instances > 0 && instance_find(sc, (uint32_t)id)) {
                        VIAMD_LOG_ERROR("Workspace: names the instance %d of a system twice; the second is left out", id);
                        continue;
                    }
                    inst = ws->num_instances == 0 ? &system_instance(sc) : instance_add(sc);
                    if (!inst) {
                        VIAMD_LOG_ERROR("Workspace: names more instances of a system than it can show (%zu); the rest are left out", INSTANCE_MAX);
                        continue;
                    }
                    inst->id = (uint32_t)id;
                    sc.next_instance_id = MD_MAX(sc.next_instance_id, (uint32_t)id + 1);
                    pi = &ws->instances[ws->num_instances++];
                    *pi = {};
                    pi->run = -1;
                    pi->frame = -1;
                } else if (!inst) {
                    continue;
                } else if (md_str_eq(ident, MD_STR_LIT("Name"))) {
                    viamd::extract_to_char_buf(inst->name, sizeof(inst->name), arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Run"))) {
                    viamd::extract_int(pi->run, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Frame"))) {
                    int frame = -1;
                    if (viamd::extract_int(frame, arg)) pi->frame = frame;
                } else if (md_str_eq(ident, MD_STR_LIT("Assembly"))) {
                    md_str_t id = {};
                    if (viamd::extract_str(id, arg)) md_str_copy_to_char_buf(pi->assembly_id, sizeof(pi->assembly_id), id);
                } else {
                    system_scene_read_instance_entry(inst->scene, ident, arg);
                }
            }
        } else if (representations_deserialize(section, state, workspace_system_fn, &pending)) {
        } else if (md_str_eq(section, MD_STR_LIT("UserBonds")) || md_str_eq(section, MD_STR_LIT("RemovedBonds"))) {
            // Each bond its atoms, and the kind it was given after them, if one was, and of a bond removed its origin; of
            // the system it names first
            const bool removed = md_str_eq(section, MD_STR_LIT("RemovedBonds"));
            WorkspacePending::System* ws = workspace_system(&pending, -1);
            bool last_taken = false;
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("System"))) {
                    int id = -1;
                    ws = viamd::extract_int(id, arg) ? workspace_system(&pending, id) : nullptr;
                    last_taken = false;
                } else if (!ws) {
                    continue;
                } else if (md_str_eq(ident, MD_STR_LIT("atoms"))) {
                    int atom_indices[2] = {-1, -1};
                    last_taken = viamd::extract_int_vec(atom_indices, 2, arg) && atom_indices[0] >= 0 && atom_indices[1] >= 0;
                    if (last_taken) {
                        const md_bond_flags_t flags = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER) | (removed ? MD_BOND_FLAG_INACTIVE : MD_BOND_FLAG_NONE);
                        md_bond_list_push(&ws->bond_edits, atom_indices[0], atom_indices[1], flags, temp_alloc);
                    }
                } else if (md_str_eq(ident, MD_STR_LIT("kind")) && last_taken) {
                    md_str_t str;
                    md_bond_flags_t kind = MD_BOND_FLAG_NONE;
                    if (!viamd::extract_str(str, arg) || !bond_kind_parse(&kind, str)) {
                        // A contact, which is no particular kind and leaves the molecules as they are
                        VIAMD_LOG_INFO("A user bond of the workspace has a kind that is not known; it is kept as a contact");
                        kind = (md_bond_flags_t)MD_BOND_KIND_CONTACT;
                    }
                    md_bond_flags_t* last = md_array_last(ws->bond_edits.flags);
                    *last = md_bond_flags_set_origin(kind | MD_BOND_FLAG_KIND_GIVEN, md_bond_origin(*last)) | (*last & MD_BOND_FLAG_INACTIVE);
                } else if (md_str_eq(ident, MD_STR_LIT("origin")) && last_taken && removed) {
                    // Where the bond removed came from: shown, of no effect
                    md_str_t str;
                    if (viamd::extract_str(str, arg)) {
                        for (int o = 0; o < 4; ++o) {
                            if (md_str_eq_cstr(md_str_trim(str), md_bond_origin_name((md_bond_origin_t)o))) {
                                md_bond_flags_t* last = md_array_last(ws->bond_edits.flags);
                                *last = md_bond_flags_set_origin(*last, (md_bond_origin_t)o);
                            }
                        }
                    }
                }
            }
        } else if (md_str_eq(section, MD_STR_LIT("Script"))) {
            // A workspace without the key is v1, as every workspace was before v2. Text is the script in it, OtherText
            // the script in the other language, which a workspace from before there were two has none of
            script::Language lang = script::Language::V1;
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Language"))) {
                    md_str_t str;
                    if (viamd::extract_str(str, arg) && md_str_eq(str, MD_STR_LIT("VIAL"))) lang = script::Language::V2;
                } else if (md_str_eq(ident, MD_STR_LIT("Text"))) {
                    md_str_t str;
                    viamd::extract_str(str, arg);
                    data->editor.SetText(std::string(str.ptr, str.len));
                } else if (md_str_eq(ident, MD_STR_LIT("OtherText"))) {
                    md_str_t str;
                    viamd::extract_str(str, arg);
                    data->script.other_text.assign(str.ptr, str.len);
                }
            }
            script_set_language(data, lang);
        } else if (md_str_eq(section, MD_STR_LIT("Selection"))) {
            // A group of the user's of the system it names (System=), the first when it names none
            md_str_t label = {};
            md_bitset_t mask = {};
            bool has_mask = false;
            int id = -1;
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Label"))) {
                    viamd::extract_str(label, arg);
                } else if (md_str_eq(ident, MD_STR_LIT("Mask"))) {
                    has_mask = deserialize_mask(&mask, arg, temp_alloc);
                } else if (md_str_eq(ident, MD_STR_LIT("System"))) {
                    viamd::extract_int(id, arg);
                }
            }
            WorkspacePending::System* ws = workspace_system(&pending, id);
            if (ws && !md_str_empty(label) && has_mask) {
                md_array_push(ws->user_groups, (WorkspacePending::System::UserGroup{label, mask}), temp_alloc);
            }
        } else if (md_str_eq(section, MD_STR_LIT("ActiveSelection"))) {
            // How it is selected, and what is selected of the system it names (System=), the first when it names none
            WorkspacePending::System* ws = workspace_system(&pending, -1);
            while (viamd::next_entry(ident, arg, state)) {
                if (md_str_eq(ident, MD_STR_LIT("Granularity"))) {
                    viamd::extract_enum(data->selection.granularity, arg, (int)SelectionGranularity::Count);
                } else if (md_str_eq(ident, MD_STR_LIT("System"))) {
                    int id = -1;
                    ws = viamd::extract_int(id, arg) ? workspace_system(&pending, id) : nullptr;
                } else if (md_str_eq(ident, MD_STR_LIT("Mask")) && ws) {
                    ws->has_selection_mask = deserialize_mask(&ws->selection_mask, arg, temp_alloc);
                } else if (md_str_eq(ident, MD_STR_LIT("Frames"))) {
                    // Of the first system's run, which the timeline shows
                    pending.has_frame_mask = deserialize_mask(&pending.frame_mask, arg, temp_alloc);
                }
            }
        } else if (md_str_eq(section, MD_STR_LIT("Windows"))) {
            // Windows not named keep their state: a workspace older than this section says nothing about them
            while (viamd::next_entry(ident, arg, state)) {
                for (size_t i = 0; i < num_workspace_windows; ++i) {
                    if (md_str_eq_cstr(ident, workspace_windows[i].name)) {
                        viamd::extract_bool(*workspace_windows[i].show, arg);
                        break;
                    }
                }
            }
        } else if (plot_layout_deserialize(state, MD_STR_LIT("Timeline"), MD_STR_LIT("TimelineSeries"), data, data->timeline.subplots, &data->timeline.num_subplots)) {
        } else if (plot_layout_deserialize(state, MD_STR_LIT("Distributions"), MD_STR_LIT("DistributionSeries"), data, data->distributions.subplots, &data->distributions.num_subplots)) {
        } else if (periodic_views_deserialize(state)) {
        } else {
            const char* before = state.text.ptr;
            viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserialize, viamd::EventPayloadType_DeserializationState, &state);
            if (state.text.ptr == before && viamd::next_entry(ident, arg, state)) {
                // Written by a component this build does not have, or by a newer version
                MD_LOG_DEBUG("Workspace: nothing reads the section [" MD_STR_FMT "], it is skipped", MD_STR_ARG(section));
            }
        }
    }

    // The box as it was drawn: lines of a pixel, of its colour
    WorkspacePending::System* first_ws = workspace_system(&pending, -1);
    if (pending.box_enabled && first_ws) {
        if (Representation* box = create_representation(*first_ws->sys, RepresentationType::SimulationBox)) {
            snprintf(box->name, sizeof(box->name), "box");
            if (pending.box_color.w > 0.0f) box->box.color = pending.box_color;
            box->box.thickness = 0.0f;
        }
    }

    // ## 3. Load the files, then apply what refers to them
    md_str_copy_to_char_buf(data->files.workspace, sizeof(data->files.workspace), filename);

    // Every system the workspace names, also those no section refers to
    workspace_system(&pending, -1);
    if (pending.skipped_systems) {
        VIAMD_LOG_ERROR("Workspace: names %zu systems; at most %zu are loaded at once, the rest are left out", pending.num_systems + pending.skipped_systems, SYSTEM_MAX);
    }

    // Each system: its own files in order, then those of each of its runs, a run being opened on the complete system.
    // This build loads RUN_MAX runs per system; what the workspace names beyond that is reported and left out.
    for (size_t k = 0; k < pending.num_systems; ++k) {
        WorkspacePending::System& ws = pending.systems[k];
        SystemContext& sc = *ws.sys;
        const bool first = k == 0;

        SystemPlan plan;
        workspace_files_plan(&plan, &pending.files, ws.id, RUN_MAX, temp_alloc);
        // The runs it opens are shown by the instances it names, below; one which names none shows its first run
        sc.instances_held = ws.num_instances > 0;
        if (plan.skipped_runs)    VIAMD_LOG_ERROR("Workspace: names %zu runs of a system; at most %zu are opened at once, the rest are left out", plan.skipped_runs + RUN_MAX, RUN_MAX);
        if (plan.skipped_inputs)  VIAMD_LOG_ERROR("Workspace: %zu files belong to no system or run it names, and are left out", plan.skipped_inputs);

        for (size_t i = 0; i < md_array_size(plan.inputs); ++i) {
            const SystemPlan::Input& in = plan.inputs[i];
            const md_str_t path = md_str_from_cstr(in.entry.path);

            // The role is the file's place; the loader is the one stored, or the one the file tells
            const loader::LoaderState loader_state = file_entry_loader_state(in.entry);

            LoadResult result = LoadResult_Failed;
            if (!md_path_is_valid(path)) {
                VIAMD_LOG_ERROR("Workspace: the file '%s' is not there", in.entry.path);
            } else if ((loader_state.flags & LoaderFlag_RequiresDialogue) && !loader_state.arg) {
                // What the dialog would have asked (the LAMMPS atom format) is not stored with the file
                VIAMD_LOG_ERROR("Workspace: '%s' has to be loaded through File > Open File..., the workspace does not say how to read it", in.entry.path);
            } else {
                // A series along its own run, which was opened before it (or is held in its place, below)
                const int along = in.role == FileRole_Series ? in.run : -1;
                result = load_data(data, sc, path, loader_state, in.role, in.opens_own_run, nullptr, nullptr, along);
            }

            if (result == LoadResult_Failed) {
                // Kept in its place, so that saving the workspace does not drop it: a file that was moved is found
                // again by putting it back, not by building the workspace again
                FileEntry missing = in.entry;
                missing.type    = loader_state.type;
                missing.missing = true;
                missing.run     = in.run;
                // A run not opened holds its place among the system's runs, without frames: run i of the list stays run
                // i of the system
                if (in.role == FileRole_Trajectory) {
                    if (sc.num_runs >= RUN_MAX) continue;
                    sc.runs[sc.num_runs++] = new SystemRun();
                }
                file_list_record(&sc.files, missing, in.role, data->allocator.persistent);
            }
        }

        // What each instance shows: a run not opened (missing, or beyond what this build opens), the system's own
        // coordinates. Before the time and the timeline, which showing the first instance's run starts over.
        sc.instances_held = false;
        for (size_t i = 0; i < ws.num_instances && i < sc.num_instances; ++i) {
            const WorkspacePending::System::Instance& pi = ws.instances[i];
            const bool opened = pi.run >= 0 && (size_t)pi.run < sc.num_runs && run_num_frames(sc, sc.runs[pi.run]) > 0;
            if (pi.run >= 0 && !opened) {
                VIAMD_LOG_ERROR("Workspace: an instance of '%s' shows a run which is not opened; it shows the system's own coordinates", sc.name);
            }
            instance_show(*sc.instances[i], opened ? pi.run : -1, pi.frame);
        }

        // The user's edits of the bonds: of the core, so given before the system is inferred, which they are part of. As
        // many atoms as when the workspace was saved, or they would join others.
        {
            const size_t loaded_atoms = sc.data.sys.atom.count;
            if (loaded_atoms > 0 && (plan.atom_count == 0 || plan.atom_count == loaded_atoms)) {
                apply_bond_edits(sc, ws.bond_edits);
            }
        }

        // What follows refers to the system as inferred (its bonds, residue types, structures): waited for here, as the
        // workspace's steps were applied after inference when it was not in the background
        system_inference_finish(sc);

        // What the workspace was saved with stays the reference while the files are missing, so that a save keeps it
        FileList& list = sc.files;
        if (list.atom_count == 0) {
            list.atom_count = plan.atom_count;
        }
        for (size_t r = 0; r < md_array_size(list.runs) && r < md_array_size(plan.runs); ++r) {
            if (list.runs[r].frame_count == 0) list.runs[r].frame_count = plan.runs[r].frame_count;
            list.runs[r].time_offset = plan.runs[r].time_offset;
        }

        const size_t num_atoms  = sc.data.sys.atom.count;
        const size_t num_frames = run_num_frames(sc);

        // The counts the workspace was saved with. When the files no longer have them, what it keeps by atom index
        // or by frame would land on other atoms and other frames, and is left out instead.
        const bool   atoms_as_saved  = plan.atom_count == 0 || plan.atom_count == num_atoms;
        const size_t saved_frames    = md_array_size(plan.runs) ? plan.runs[0].frame_count : 0;
        const bool   frames_as_saved = saved_frames == 0 || saved_frames == num_frames;
        if (num_atoms > 0 && !atoms_as_saved) {
            VIAMD_LOG_ERROR("Workspace: the system has %zu atoms, but had %zu when the workspace was saved; its selections, groups and bonds are left out", num_atoms, plan.atom_count);
        }
        if (num_frames > 0 && !frames_as_saved) {
            VIAMD_LOG_ERROR("Workspace: the trajectory has %zu frames, but had %zu when the workspace was saved; its frame and timeline ranges are left out", num_frames, saved_frames);
        }

        if (first) pending.frames_as_saved = num_frames > 0 && frames_as_saved;

        // Loading a trajectory sets the filter to all of it; the workspace's comes after. The time and the view of the
        // timeline are of the clock, once every system is loaded (below).
        if (first && num_frames > 0 && frames_as_saved) {
            const double last = (double)(num_frames - 1);
            // A workspace from before the periods: its filter becomes a period, over which what was shown over the
            // filter is shown
            if (!pending.has_periods && pending.filter_enabled) {
                TimelinePeriod* r = nullptr;
                if (pending.temporal_window_enabled) {
                    const double  extent = pending.temporal_window_extent > 0.0 ? pending.temporal_window_extent : 10.0;
                    const int32_t half = (int32_t)round(extent * 0.5);
                    r = timeline_periods_add(&data->timeline.periods, "Filter", -half, half, true);
                } else {
                    const double beg = pending.has_filter_range ? MD_CLAMP(pending.filter_beg_frame, 0.0, last) : 0.0;
                    const double end = pending.has_filter_range ? MD_CLAMP(pending.filter_end_frame, beg, last) : last;
                    r = timeline_periods_add(&data->timeline.periods, "Filter", (int32_t)round(beg), (int32_t)round(end), false);
                }
                if (r) {
                    data->timeline.periods.legacy_filter_id = r->id;
                    if (pending.num_period_ends < TIMELINE_MAX_PERIODS) pending.period_ends[pending.num_period_ends++] = { r->id, true, md_unit_none() };
                }
            }
        }

        // Masks of atoms that are not there are dropped rather than applied to whatever holds those indices now
        auto mask_fits = [num_atoms, atoms_as_saved](const md_bitset_t* mask) {
            return num_atoms > 0 && atoms_as_saved && (md_bitset_empty(*mask) || mask->end <= num_atoms);
        };
        if (ws.has_selection_mask && mask_fits(&ws.selection_mask)) {
            md_bitmask_assign(&sc.selection.selection_mask, ws.selection_mask);
        }
        // The user's groups, after those the files gave. One of a name a group of the files has (a stored selection named
        // SOL, from before they were groups) is given a name of its own rather than merged with that group by group('SOL').
        for (size_t i = 0; i < md_array_size(ws.user_groups); ++i) {
            const WorkspacePending::System::UserGroup& g = ws.user_groups[i];
            if (!mask_fits(&g.atoms)) continue;
            md_str_t name = g.name;
            char buf[96];
            for (int k = 2; !user_group_name_free(sc.data.sys, name); ++k) {
                snprintf(buf, sizeof(buf), MD_STR_FMT "%d", MD_STR_ARG(g.name), k);
                name = md_str_from_cstr(buf);
            }
            if (!md_str_eq(name, g.name)) {
                VIAMD_LOG_INFO("Workspace: the system has a group '" MD_STR_FMT "' of its own; the workspace's is named '" MD_STR_FMT "'", MD_STR_ARG(g.name), MD_STR_ARG(name));
            }
            md_system_group_add(&sc.data.sys, name, g.atoms, MD_GROUP_ORIGIN_USER);
        }
        // Frames that are not there (another run, or a shorter one) are dropped likewise
        if (first && pending.has_frame_mask && frames_as_saved && (md_bitset_empty(pending.frame_mask) || pending.frame_mask.end <= (uint32_t)num_frames)) {
            frame_selection_modify(data, pending.frame_mask, SelectionOperator::Set);
        }
        if (ws.has_recenter_target && mask_fits(&ws.recenter_target)) {
            md_bitset_assign(&sc.operations.selection_mask, ws.recenter_target, data->allocator.persistent);
            recenter_update_target_data(sc);
        }

        // The assembly each instance shows, by its id in the files loaded
        for (size_t i = 0; i < sc.num_instances; ++i) {
            const char* id = i < ws.num_instances ? ws.instances[i].assembly_id : "";
            if (!id[0] && i == 0) id = ws.assembly_id;
            if (!id[0]) continue;
            const md_assembly_data_t& a = sc.data.sys.assembly;
            for (size_t j = 0; j < a.count; ++j) {
                if (md_str_eq(md_assembly_id(&a, j), md_str_from_cstr(id))) sc.instances[i]->scene.assembly = (int)j;
            }
            if (sc.instances[i]->scene.assembly < 0) {
                VIAMD_LOG_ERROR("Workspace: the system has no assembly '%s'; it is shown as it is", id);
            }
        }
    }

    // A system it names with no files at all is not one: taken out (the first stays, empty)
    for (size_t k = pending.num_systems; k-- > 0;) {
        SystemContext* sc = pending.systems[k].sys;
        if (!system_loaded(*sc)) {
            remove_all_representations(*sc);
            system_release(data, sc);
        }
    }

    // The time shown, and the view of the timeline, on the clock of every run loaded: a time as it was saved (in the
    // clock's unit, converted when the reference is in another), a frame of the first system's run as it was before
    // the clock, while that run has its frames as saved
    clock_update(data);
    if (clock_reference(data)) {
        double beg, end;
        clock_span(data, &beg, &end);
        const md_unit_t unit = clock_unit(data);
        double f = 1.0;
        const bool time_fits = pending.has_time && (md_unit_is_none(pending.time_unit) ? md_unit_is_none(unit) : (!md_unit_is_none(unit) && md_unit_conversion_factor(&f, pending.time_unit, unit)));
        if (time_fits) {
            data->animation.time = MD_CLAMP(pending.time * f, beg, end);
        } else if (pending.has_frame && pending.frames_as_saved) {
            clock_set_frame(data, pending.frame);
        }
        if (pending.has_view_range && (time_fits || pending.frames_as_saved)) {
            double b = pending.view_range[0];
            double e = pending.view_range[1];
            if (pending.view_range_in_frames) {
                const SystemInstance* ref = clock_reference(data);
                b = run_clock_time(data, *ref->sys, instance_run(*ref), b);
                e = run_clock_time(data, *ref->sys, instance_run(*ref), e);
            } else {
                b *= f;
                e *= f;
            }
            b = MD_CLAMP(b, beg, end);
            e = MD_CLAMP(e, beg, end);
            if (e > b) data->timeline.view_range = {b, e};
        }

        // The periods on the clock: frames of the reference's run (the first system's, which the frames were of) at
        // their times, times converted into the clock's unit. Left as they are where they cannot be.
        const SystemInstance* ref = clock_reference(data);
        const ClockAxis ref_axis = run_clock_axis(data, *ref->sys, instance_run(*ref));
        for (int i = 0; i < pending.num_period_ends; ++i) {
            const WorkspacePending::PeriodEnds& pe = pending.period_ends[i];
            TimelinePeriod* r = timeline_periods_find(&data->timeline.periods, pe.id);
            if (!r) continue;
            double pf = 1.0;
            if (pe.in_frames) {
                timeline_period_frames_to_time(r, ref_axis);
            } else if (md_unit_is_none(pe.unit) ? md_unit_is_none(unit) : (!md_unit_is_none(unit) && md_unit_conversion_factor(&pf, pe.unit, unit))) {
                r->beg *= pf;
                r->end *= pf;
            }
        }
    }

    // The camera flies in: from every system to where the workspace looked from. Without a stored camera, both are
    // every system (or the assembly shown, of one).
    reset_view_scene(&data->view.camera, *data, data->view.camera.fov_y, data->view.region);
    data->view.trackball_param.max_distance = MD_MAX(data->view.trackball_param.max_distance, data->view.camera.distance * 10.0f);
    data->view.target = pending.has_camera ? pending.camera : (ViewTransform)data->view.camera;

    // What was shown over the timeline filter, before there were periods, is shown over the period the filter followed.
    // The residues take theirs on the end of the load, after which there is no filter to speak of.
    plot_resolve_legacy_periods(data, data->timeline.subplots, PLOT_MAX_SUBPLOTS);
    plot_resolve_legacy_periods(data, data->distributions.subplots, PLOT_MAX_SUBPLOTS);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserializeEnd, viamd::EventPayloadType_DeserializationState, &state);
    data->timeline.periods.legacy_filter_id = 0;
}

bool save_workspace(ApplicationState* app_state, md_str_t filename) {
    md_allocator_i* temp_alloc = app_state->allocator.frame;

    viamd::serialization_state_t state {
        .filename = filename,
        .sw = md_str_writer_create(temp_alloc),
    };

    constexpr md_str_t header_snippet = MD_STR_LIT(
        R"(
        #01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#
        #                                                                                #
        #            VIAMD — Visual Interactive Analysis of Molecular Dynamics           #
        #                                                                                #
        #                    github: https://github.com/scanberg/viamd                   #
        #                 manual: https://github.com/scanberg/viamd/wiki                 #
        #                    youtube playlist: https://bit.ly/4aRsPrh                    #
        #                                twitter: @VIAMD_                                #
        #                                                                                #
        #                If you use VIAMD in your research, please cite:                 #
        #   "VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics"    #
        #       Robin Skånberg, Ingrid Hotz, Anders Ynnerman, and Mathieu Linares        #
        #                 J. Chem. Inf. Model. 2023, 63, 23, 7382–7391                   #
        #                   https://doi.org/10.1021/acs.jcim.3c01033                     #
        #                                                                                #
        #01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#
        )");

    // Write big ass header
    state.sw += header_snippet;
    state.sw += '\n';

    viamd::write_section_header(state, MD_STR_LIT("Workspace"));
    viamd::write_int(state, MD_STR_LIT("Version"), WORKSPACE_VERSION);

    // The systems, their runs and their files, missing ones too: what could not be loaded is still what the
    // workspace is of. A system's id is the one every section which refers to it uses.
    {
        int next_run_id = 0;
        for (size_t i = 0; i < app_state->num_systems; ++i) {
            const SystemContext& sc = *app_state->systems[i];
            if (file_list_count(&sc.files) > 0) {
                next_run_id = file_list_serialize(state, &sc.files, (int)sc.id, next_run_id);
            }
        }
    }

    viamd::write_section_header(state, MD_STR_LIT("Animation"));
    // The time on the clock, in its unit, and the frame of the reference's run at it, which is all a workspace written
    // before the clock knows
    viamd::write_dbl(state, MD_STR_LIT("Time"), app_state->animation.time);
    {
        char unit[64] = "";
        md_unit_print(unit, sizeof(unit), clock_unit(app_state));
        if (unit[0]) viamd::write_str(state, MD_STR_LIT("TimeUnit"), md_str_from_cstr(unit));
    }
    viamd::write_dbl(state, MD_STR_LIT("Frame"), clock_frame(app_state));
    viamd::write_flt(state, MD_STR_LIT("Fps"), app_state->animation.fps);
    viamd::write_flt(state, MD_STR_LIT("Tension"), app_state->animation.tension);
    viamd::write_int(state, MD_STR_LIT("Stride"), app_state->animation.stride);
    viamd::write_int(state, MD_STR_LIT("Interpolation"), (int)app_state->animation.interpolation);
    viamd::write_int(state, MD_STR_LIT("Loop"), (int)app_state->animation.loop);
    if (timeline_periods_find(&app_state->timeline.periods, app_state->animation.loop_period)) {
        viamd::write_int(state, MD_STR_LIT("LoopPeriod"), (int)app_state->animation.loop_period);
    }
    viamd::write_int(state, MD_STR_LIT("Transport"), (int)app_state->animation.transport);

    plot_layout_serialize(state, MD_STR_LIT("Timeline"), MD_STR_LIT("TimelineSeries"), app_state->timeline.subplots, app_state->timeline.num_subplots);
    plot_layout_serialize(state, MD_STR_LIT("Distributions"), MD_STR_LIT("DistributionSeries"), app_state->distributions.subplots, app_state->distributions.num_subplots);
    periodic_views_serialize(state);

    {
        // The view in the run's time and its unit, as the timeline holds it: the same trajectory, which its frame count
        // checks on loading, has the same time. Of the first system's run, which the timeline shows.
        const auto& tl = app_state->timeline;
        viamd::write_section_header(state, MD_STR_LIT("TimelineView"));
        viamd::write_int(state, MD_STR_LIT("PeriodLayout"), (int)tl.periods.layout);
        viamd::write_bool(state, MD_STR_LIT("PeriodShadeAll"), tl.periods.shade_all);
        if (run_num_frames(app_state->system_context) > 0) {
            const double view[2] = { tl.view_range.beg_x, tl.view_range.end_x };
            viamd::write_dbl_vec(state, MD_STR_LIT("ViewTime"), view, 2);
        }
        // On the clock, in its unit
        char unit[64] = "";
        md_unit_print(unit, sizeof(unit), clock_unit(app_state));
        timeline_periods_serialize(state, &tl.periods, unit);
    }

    {
        const auto& v = app_state->visuals;
        viamd::write_section_header(state, MD_STR_LIT("RenderSettings"));
        viamd::write_vec3(state, MD_STR_LIT("BackgroundColor"), v.background.color);
        viamd::write_int(state,  MD_STR_LIT("ShadingMode"), (int)v.shading.mode);
        viamd::write_int(state,  MD_STR_LIT("ShadingColors"), (int)v.shading.colors);
        viamd::write_flt(state,  MD_STR_LIT("ShadingLight"), v.shading.light);
        viamd::write_flt(state,  MD_STR_LIT("ShadingEnvironment"), v.shading.environment);
        viamd::write_flt(state,  MD_STR_LIT("ShadingAmbient"), v.shading.ambient);
        viamd::write_flt(state,  MD_STR_LIT("ShadingSpecular"), v.shading.specular);
        viamd::write_bool(state, MD_STR_LIT("TonemapEnabled"), v.tonemapping.enabled);
        viamd::write_int(state,  MD_STR_LIT("Tonemapper"), (int)v.tonemapping.tonemapper);
        viamd::write_flt(state,  MD_STR_LIT("TonemapExposure"), v.tonemapping.exposure);
        viamd::write_flt(state,  MD_STR_LIT("TonemapGamma"), v.tonemapping.gamma);
        viamd::write_bool(state, MD_STR_LIT("DofEnabled"), v.dof.enabled);
        viamd::write_flt(state,  MD_STR_LIT("DofAperture"), v.dof.aperture);
        // The ambient occlusion and the anti-aliasing are the computer's settings (the .ini), the simulation box a
        // representation, the periodic images and the assembly the system's ([SystemScene])
    }

    // The target, not the camera: the camera eases towards it over several frames (camera_animate),
    // so saving during a fly-in or right after a drag would store a point along the way
    viamd::write_section_header(state, MD_STR_LIT("Camera"));
    viamd::write_vec3(state, MD_STR_LIT("Position"), app_state->view.target.position);
    viamd::write_quat(state, MD_STR_LIT("Orientation"), app_state->view.target.orientation);
    viamd::write_flt(state,  MD_STR_LIT("Distance"), app_state->view.target.distance);
    viamd::write_int(state,  MD_STR_LIT("Mode"), (int)app_state->view.mode);
    viamd::write_flt(state,  MD_STR_LIT("FovY"), app_state->view.camera.fov_y);

    // Each system as the scene shows it, and its representations
    for (size_t i = 0; i < app_state->num_systems; ++i) {
        const SystemContext& sc = *app_state->systems[i];
        if (!system_loaded(sc)) continue;
        system_scene_serialize(state, sc);
        representations_serialize(state, sc);
    }

    {
        std::string text = app_state->editor.GetText();
        viamd::write_section_header(state, MD_STR_LIT("Script"));
        // Written for v2 alone, so a v1 workspace is as it was
        if (app_state->script.language == script::Language::V2) viamd::write_str(state, MD_STR_LIT("Language"), MD_STR_LIT("VIAL"));
        viamd::write_str(state, MD_STR_LIT("Text"), md_str_t{text.c_str(), text.size()});
        const std::string& other = app_state->script.other_text;
        if (!other.empty()) viamd::write_str(state, MD_STR_LIT("OtherText"), md_str_t{other.c_str(), other.size()});
    }

    // The user's groups of each system, in the section the stored selections were written in, which they were: of the
    // system named (System=), the first's naming none, as they were. The files give the others.
    for (size_t k = 0; k < app_state->num_systems; ++k) {
        const SystemContext& sc = *app_state->systems[k];
        if (!system_loaded(sc)) continue;
        const md_group_data_t& groups = sc.data.sys.group;
        for (size_t i = 0; i < groups.count; ++i) {
            if (groups.origin[i] != MD_GROUP_ORIGIN_USER) continue;
            viamd::write_section_header(state, MD_STR_LIT("Selection"));
            if (k > 0) viamd::write_int(state, MD_STR_LIT("System"), (int)sc.id);
            viamd::write_str(state, MD_STR_LIT("Label"), groups.name[i]);
            viamd::write_bitset(state, MD_STR_LIT("Mask"), &groups.atoms[i]);
        }
    }

    // How it is selected, then what is selected of each system, and the frames selected of the first system's run
    viamd::write_section_header(state, MD_STR_LIT("ActiveSelection"));
    viamd::write_int(state, MD_STR_LIT("Granularity"), (int)app_state->selection.granularity);
    const md_bitset_t frames = md_bitmask_view(&app_state->selection.frame_mask);
    if (!md_bitset_empty(frames)) {
        viamd::write_bitset(state, MD_STR_LIT("Frames"), &frames);
    }
    for (size_t k = 0; k < app_state->num_systems; ++k) {
        const SystemContext& sc = *app_state->systems[k];
        const md_bitset_t selection = md_bitmask_view(&sc.selection.selection_mask);
        if (system_loaded(sc) && !md_bitset_empty(selection)) {
            viamd::write_int(state, MD_STR_LIT("System"), (int)sc.id);
            viamd::write_bitset(state, MD_STR_LIT("Mask"), &selection);
        }
    }

    // The user's edits of the bonds of each system (md_system_t::bond_edit), whether the system is inferred of them yet
    // or not: the bonds added, of the graph and interactions, with the kind they were given, then the bonds removed
    // (given by the data or inferred), each with its kind, which tells its range (the graph, the elastic bonds or the
    // interactions), and where it came from
    for (size_t k = 0; k < app_state->num_systems; ++k) {
        const SystemContext& sc = *app_state->systems[k];
        if (!system_loaded(sc)) continue;
        const md_bond_list_t& edits = sc.data.sys.bond_edit;
        const size_t num_edits = md_bond_list_count(&edits);
        for (int removed = 0; removed < 2; ++removed) {
            bool has_section = false;
            for (size_t i = 0; i < num_edits; ++i) {
                const md_bond_flags_t flags = edits.flags[i];
                if (((flags & MD_BOND_FLAG_INACTIVE) != 0) != (removed == 1)) continue;
                if (!has_section) {
                    viamd::write_section_header(state, removed ? MD_STR_LIT("RemovedBonds") : MD_STR_LIT("UserBonds"));
                    viamd::write_int(state, MD_STR_LIT("System"), (int)sc.id);
                    has_section = true;
                }
                viamd::write_int_vec(state, MD_STR_LIT("atoms"), edits.pairs[i].idx, 2);
                if (removed || (flags & MD_BOND_FLAG_KIND_GIVEN)) {
                    viamd::write_str(state, MD_STR_LIT("kind"), md_str_from_cstr(bond_kind_str(flags)));
                }
                if (removed) {
                    viamd::write_str(state, MD_STR_LIT("origin"), md_str_from_cstr(md_bond_origin_name(md_bond_origin(flags))));
                }
            }
        }
    }

    viamd::write_section_header(state, MD_STR_LIT("Windows"));
    for (size_t i = 0; i < num_workspace_windows; ++i) {
        viamd::write_bool(state, md_str_from_cstr(workspace_windows[i].name), *workspace_windows[i].show);
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSerialize, viamd::EventPayloadType_SerializationState, &state);

    // Only now, with all of it in hand, is the file touched
    const md_str_t text = md_str_writer_to_str(&state.sw);
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Could not open workspace file for writing: '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }
    const size_t written = md_file_write(file, md_str_ptr(text), md_str_len(text));
    md_file_close(&file);
    if (written != md_str_len(text)) {
        VIAMD_LOG_ERROR("Could not write all of the workspace file: '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

    // Saved there, it is the workspace from here on: the next plain save goes to the same file
    md_str_copy_to_char_buf(app_state->files.workspace, sizeof(app_state->files.workspace), filename);
    VIAMD_LOG_SUCCESS("Saved workspace '" MD_STR_FMT "'", MD_STR_ARG(filename));
    return true;
}

script::FilterCtx filter_context(const SystemContext& sys) {
    // On the coordinates of the first instance: those the selection and the script are of
    return filter_context(system_instance(sys));
}

script::FilterCtx filter_context(const SystemInstance& inst) {
    const SystemContext& sys = *inst.sys;
    const ApplicationState* app = sys.app;
    const bool names = app && &sys == &app->system_context;
    const SystemRun* run = instance_run(inst);
    return script::FilterCtx{&sys.data.sys, &inst.state, names ? app->script.ir : nullptr, nullptr, run ? md_str_from_cstr(run->path) : md_str_t{}};
}

// --- THE USER'S GROUPS ---
bool user_group_name_free(const md_system_t& sys, md_str_t name, size_t self) {
    if (md_str_empty(name)) return false;
    for (size_t i = 0; i < sys.group.count; ++i) {
        if (i != self && md_str_eq(sys.group.name[i], name)) return false;
    }
    return true;
}

int store_selection_as_group(SystemContext& sys) {
    const md_bitset_t atoms = md_bitmask_view(&sys.selection.selection_mask);
    if (md_bitset_empty(atoms)) return -1;
    char name[32];
    for (int k = 1;; ++k) {
        snprintf(name, sizeof(name), "sel%d", k);
        if (user_group_name_free(sys.data.sys, md_str_from_cstr(name))) break;
    }
    const size_t idx = md_system_group_add(&sys.data.sys, md_str_from_cstr(name), atoms, MD_GROUP_ORIGIN_USER);
    VIAMD_LOG_INFO("The selection of '%s' is the group '%s' (System > Groups), group('%s') in VIAL", sys.name, name, name);
    return (int)idx;
}

// "ground_state" reads as "Ground State" in a menu. The path segment is the identity; this is
// presentation only, which is why it lives here and not in mdlib. Writes into a caller buffer so
// gathering stays allocation free.
int dipole_label_pretty(char* buf, size_t cap, md_str_t group) {
    if (!buf || cap == 0) return 0;
    size_t n = MD_MIN(group.len, cap - 1);
    bool boundary = true;
    for (size_t i = 0; i < n; ++i) {
        char c = group.ptr[i];
        if (c == '_' || c == '-') {
            buf[i] = ' ';
            boundary = true;
            continue;
        }
        buf[i] = (boundary && c >= 'a' && c <= 'z') ? (char)(c - 'a' + 'A') : c;
        boundary = false;
    }
    buf[n] = '\0';
    return (int)n;
}

int dipole_entry_label(char* buf, size_t cap, const DipoleGroup& group, uint32_t index) {
    if (!buf || cap == 0) return 0;

    int len = dipole_label_pretty(buf, cap, group.label);
    if (group.count > 1 && (size_t)len + 1 < cap) {
        len += snprintf(buf + len, cap - (size_t)len, " %u", index + 1);
    }
    return len;
}

// The origin lives beside the vector: same group, last path segment swapped. Derived rather than
// carried around, so a key is the only thing anyone has to hold on to.
static const md_attribute_t* dipole_origin_of(const md_system_t& sys, const md_attribute_t* vec) {
    return md_attributes_sibling(&sys.attributes, vec, MD_STR_LIT("origin"));
}

// A vector and an origin form a group only if both are there, both are 3 residue, and the origin
// is addressed by the vector's index space - either one per element, or a single anchor shared over
// all of them, which is how one centre of charge serves every excited state.
static bool dipole_group_qualifies(const md_attribute_t* vec, const md_attribute_t* org) {
    if (!vec || !org) return false;
    if (vec->format.components != 3 || org->format.components != 3) return false;

    const size_t num_elem = md_attribute_value_count(&vec->format);
    const size_t num_org  = md_attribute_value_count(&org->format);
    return num_elem > 0 && (num_org == num_elem || num_org == 1);
}

size_t dipole_groups_gather(DipoleGroup out[], size_t cap, const md_system_t& sys) {
    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter_children(&sys.attributes, MD_STR_LIT("dipole")); md_attributes_next(&it);) {
        const md_attribute_t* vec = md_attributes_find_in(&sys.attributes, it.child_path, MD_STR_LIT("vector"));
        if (!vec || !dipole_group_qualifies(vec, dipole_origin_of(sys, vec))) continue;

        if (out && count < cap) {
            out[count] = {
                .key   = vec->id,
                .label = it.child,
                .count = (uint32_t)md_attribute_value_count(&vec->format),
                .unit  = vec->unit,
            };
        }
        count += 1;
    }

    return count;
}

bool dipole_group_from_key(DipoleGroup* out, const md_system_t& sys, md_attribute_id_t key) {
    const md_attribute_t* vec = md_attributes_get(&sys.attributes, key);
    if (!vec || !dipole_group_qualifies(vec, dipole_origin_of(sys, vec))) return false;

    if (out) {
        // The group name is the path segment before the leaf, which is the same string the gather
        // above hands back - so a label built from either route reads identically.
        md_str_t label = vec->path;
        label = md_str_substr(label, 0, md_str_rfind_char(label, '/'));
        label = md_str_substr(label, md_str_rfind_char(label, '/') + 1);

        *out = {
            .key   = key,
            .label = label,
            .count = (uint32_t)md_attribute_value_count(&vec->format),
            .unit  = vec->unit,
        };
    }
    return true;
}

bool dipole_moment_read(md_float3_t* out_vec, md_float3_t* out_origin, const md_system_t& sys, md_attribute_id_t key, uint32_t index) {
    const md_attribute_t* vec = md_attributes_get(&sys.attributes, key);
    if (!vec) return false;
    const md_attribute_t* org = dipole_origin_of(sys, vec);
    if (!dipole_group_qualifies(vec, org)) return false;
    if (index >= md_attribute_value_count(&vec->format)) return false;

    // One slice, clamped to each half's own rank. A shared origin is rank 0 and takes no index; a
    // per state one is rank 1 and takes the same index as the vector. That clamp is the whole cost
    // of not requiring the two to have equal shapes.
    const md_attribute_slice_t vec_slice = vec->format.rank > 0 ? md_attribute_slice_1(index) : md_attribute_slice_all();
    const md_attribute_slice_t org_slice = org->format.rank > 0 ? vec_slice : md_attribute_slice_all();

    // The vector comes back as stored, whatever the producer chose; the origin must be a length in
    // system space, and extraction refuses if it is not.
    float v[3] = {0, 0, 0};
    float o[3] = {0, 0, 0};
    md_attribute_extract_f32(v, MD_COUNTOF(v), &sys.attributes, vec, vec_slice, md_unit_none());
    md_attribute_extract_f32(o, MD_COUNTOF(o), &sys.attributes, org, org_slice, md_unit_angstrom());

    if (out_vec)    *out_vec    = md_float3_set(v[0], v[1], v[2]);
    if (out_origin) *out_origin = md_float3_set(o[0], o[1], o[2]);
    return true;
}


// A per atom scalar field: values one residue wide, the atom axis LAST, and at most one axis of
// variants ahead of it. mdlib deliberately does not know that an "atom/..." path is over this
// system's atoms - categories are not predeclared - so this is where that convention is checked.
static bool atom_attribute_qualifies(const md_attribute_t* attr, const md_system_t& sys) {
    const md_attribute_format_t& fmt = attr->format;
    if (fmt.rank < 1 || fmt.rank > 2) return false;
    if (fmt.components != 1) return false;
    return fmt.shape[fmt.rank - 1] == (uint32_t)sys.atom.count;
}

size_t atom_attribute_query(md_attribute_id_t out_ids[], size_t cap, const md_system_t& sys) {
    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, MD_STR_LIT("atom")); md_attributes_next(&it);) {
        if (!atom_attribute_qualifies(it.attr, sys)) continue;

        if (out_ids && count < cap) {
            out_ids[count] = it.attr->id;
        }
        count += 1;
    }

    return count;
}

md_str_t atom_attribute_label(const md_attribute_t* attr) {
    if (!attr) return md_str_t{};
    // An empty label is a valid state, and the leaf is what the path spells for itself.
    return md_str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
}

int atom_attribute_variant_count(const md_attribute_t* attr) {
    if (!attr) return 0;
    return attr->format.rank > 1 ? (int)attr->format.shape[0] : 1;
}

ColorScaleSpan atom_attribute_span(const md_attributes_t* table, const md_attribute_t* attr, const md_bitset_t* mask) {
    ColorScaleSpan span;
    if (!attr) return span;

    const size_t num_values = md_attribute_element_count(&attr->format);
    const size_t num_atoms  = attr->format.rank > 0 ? attr->format.shape[attr->format.rank - 1] : 0;
    if (num_values == 0 || num_atoms == 0) return span;

    md_temp_scope_t temp = md_temp_begin();
    float* values = (float*)md_temp_alloc(temp, sizeof(float) * num_values);

    // Deliberately the whole attribute and not one variant: a span measured per variant would make
    // the colours shift as the index slider moves, which reads as the data changing.
    if (values && md_attribute_extract_f32(values, num_values, table, attr, md_attribute_slice_all(), md_unit_none()) == num_values) {
        float value_min =  FLT_MAX;
        float value_max = -FLT_MAX;
        size_t num_present = 0;
        for (size_t i = 0; i < num_values; ++i) {
            if (mask && !md_bitset_test(*mask, (uint32_t)(i % num_atoms))) continue;
            if (atom_attribute_value_absent(values[i])) continue;
            value_min = MD_MIN(value_min, values[i]);
            value_max = MD_MAX(value_max, values[i]);
            num_present += 1;
        }
        if (num_present > 0) {
            span = color_scale_span(value_min, value_max);
        }
    }
    md_temp_end(temp);

    return span;
}

void atom_attribute_legend_label(char* buf, size_t cap, const md_attribute_t* attr, int variant_idx) {
    MD_ASSERT(buf && cap > 0);
    const md_str_t label = atom_attribute_label(attr);
    const int num_variants = atom_attribute_variant_count(attr);
    if (num_variants > 1) {
        snprintf(buf, cap, "%.*s [%d/%d]", (int)label.len, label.ptr, MD_CLAMP(variant_idx, 0, num_variants - 1) + 1, num_variants);
    } else {
        snprintf(buf, cap, "%.*s", (int)label.len, label.ptr);
    }
}

void atom_attribute_select(AtomAttributeColoring* coloring, md_attribute_id_t key, const md_system_t& sys) {
    MD_ASSERT(coloring);

    coloring->key = key;
    coloring->variant_idx = 0;
    coloring->span_hash = 0;    // measured again, over the representation's atoms, with its colours

    // Started out from the values of every atom: the representation's own are not at hand here
    const ColorScaleSpan span = atom_attribute_span(&sys.attributes, md_attributes_get(&sys.attributes, key), nullptr);
    const bool signed_values = color_scale_span_signed(span);

    ColorScale& scale = coloring->scale;
    scale.symmetric  = signed_values;
    scale.auto_range = true;
    if (signed_values && scale.colormap == DEFAULT_COLORMAP) {
        scale.colormap = COLOR_SCALE_COLORMAP_RDBU;
    } else if (!signed_values && scale.colormap == COLOR_SCALE_COLORMAP_RDBU) {
        scale.colormap = DEFAULT_COLORMAP;
    }
    coloring->span = span;
    color_scale_update_range(&scale, span);
}

// ---------------------------------------------------------------------------
// Per dataset GPU data
// ---------------------------------------------------------------------------
// Built from the system's own basis/ attributes, so it needs no loader and no residue - which is
// the point of publishing the basis in the first place. Lives beside gl_mol because it is the same
// kind of thing: derived FROM the system so that something can draw or evaluate it, and never read
// back by mdlib.
//
// The device scratch is grown here rather than allocated per dataset. A wider basis loaded later
// grows it; a narrower one leaves it alone, because it is scratch and only the maximum matters.

#if MD_ENABLE_GPU
static bool dmat_build_take(ApplicationState* state, bool wait);
static void localization_gpu_free(ApplicationState* state);
static bool dpack_build_take(ApplicationState* state, bool wait);
static bool gpu_density_jobs_cancel(ApplicationState* state, uint32_t vol_tex);
#endif

void system_gpu_data_free(SystemContext& sys) {
    ApplicationState* state = sys.app;
    MD_ASSERT(state);
#if MD_ENABLE_GPU
    // the evaluations not yet all issued are of what goes: dropped, so that the drain waits for what is in flight only
    if (state->gpu_stream) gpu_density_jobs_cancel(state, 0);
#endif
    gpu_volume_jobs_drain(state);   // what the evaluations under way read goes below
#if MD_ENABLE_GPU
    if (sys.data.gpu_basis) {
        md_gto_gpu_basis_destroy(state->gpu_stream, sys.data.gpu_basis);
        sys.data.gpu_basis = nullptr;
    }
    if (sys.data.gpu_atoms) {
        md_gpu_free(state->gpu_stream, sys.data.gpu_atoms);
        sys.data.gpu_atoms = 0;
    }
    sys.data.gpu_atoms_hash = 0;
    // What workers have under way is the first system's (these read state->mold): taken, then freed below
    const bool first = &sys == &state->system_context;
    if (first) dmat_build_take(state, true);
    if (sys.data.gpu_dmat) {
        md_gto_gpu_dmat_destroy(state->gpu_stream, sys.data.gpu_dmat);
        sys.data.gpu_dmat = nullptr;
    }
    sys.data.gpu_dmat_key = 0;
    sys.data.gpu_dmat_failed_key = 0;
    if (first) {
        localization_gpu_free(state);
        dpack_build_take(state, true);
    }
    if (sys.data.gpu_dpack) {
        md_gpu_free(state->gpu_stream, sys.data.gpu_dpack);
        sys.data.gpu_dpack = 0;
    }
    sys.data.gpu_dpack_key = 0;
    sys.data.gpu_dpack_failed_key = 0;
#else
    (void)state;
#endif
}

bool system_gpu_data_update(ApplicationState* state, double cutoff) {
    MD_ASSERT(state);
#if MD_ENABLE_GPU
    system_gpu_data_free(state->system_context);

    if (!state->gpu_device) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &state->mold.sys.attributes, md_temp_allocator(temp))) {
        // A system with no basis published is the normal case, not a failure.
        return false;
    }

    md_gto_gpu_basis_desc_t desc = { .basis = &basis, .cutoff = cutoff };
    state->mold.gpu_basis = md_gto_gpu_basis_create(state->gpu_stream, &desc);
    if (!state->mold.gpu_basis) {
        MD_LOG_ERROR("Failed to upload the GTO basis to the device");
        return false;
    }

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    const size_t num_atoms = md_gto_gpu_basis_num_atoms(state->mold.gpu_basis);

    state->mold.gpu_atoms = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(num_atoms)).gpu;
    state->mold.gpu_atoms_hash = 0;

    // Density coefficients are the larger of the two packings, so one size covers both the density
    // and the MO evaluation paths.
    const size_t coeff_size = md_gto_gpu_coeff_size_density(num_cgtos);
    if (coeff_size > state->gpu_coeff_capacity) {
        if (state->gpu_coeff) {
            md_gpu_free(state->gpu_stream, state->gpu_coeff);
        }
        state->gpu_coeff = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_DEVICE, coeff_size).gpu;
        state->gpu_coeff_capacity = state->gpu_coeff ? coeff_size : 0;
    }

    return state->mold.gpu_atoms != 0 && state->gpu_coeff != 0;
#else
    (void)state; (void)cutoff;
    return false;
#endif
}

// ---------------------------------------------------------------------------
// Orbital evaluation
// ---------------------------------------------------------------------------
// Nothing here holds a reader. The basis and the AO coefficients are attributes on the system, the
// atom positions are the system's own state, and the grid, the destination texture and the
// evaluation parameters belong to the caller. There is no vlx pointer in scope, no event, and no
// residue to ask - which is the whole point of publishing the coefficients.
//
// The basis is rebuilt per call for now. That is one interleave over a few hundred shells, but it
// is the obvious thing for a representation to cache, keyed on the ids of the basis/ attributes it
// was built from.
double* orbital_coefficients_extract(size_t* out_num_ao, md_temp_scope_t temp, const md_system_t& sys, md_str_t coefficient_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, coefficient_path);
    if (!attr) {
        MD_LOG_DEBUG("No orbital coefficients published at '" MD_STR_FMT "'", MD_STR_ARG(coefficient_path));
        return nullptr;
    }

    // Whatever the attribute's own rank, the slice has to leave exactly ONE row of AO coefficients:
    // {M,A} sliced by the orbital, {S,L,A} sliced by state and lambda. Asking the slice for its
    // format rather than the attribute is what makes those the same call.
    md_attribute_format_t format = {};
    if (!md_attribute_slice_format(&format, attr, slice)) {
        MD_LOG_ERROR("The slice does not address '" MD_STR_FMT "'", MD_STR_ARG(coefficient_path));
        return nullptr;
    }
    if (format.rank != 1 || format.components != 1) {
        MD_LOG_ERROR("'" MD_STR_FMT "' does not slice down to one row of orbital coefficients", MD_STR_ARG(coefficient_path));
        return nullptr;
    }

    // Ask the slice how big it is, then allocate for exactly that. The size comes from the format
    // alone, so this same shape works whether the coefficients are stored or worked out on demand.
    const size_t num_ao = md_attribute_slice_count(attr, slice);
    if (num_ao == 0) {
        return nullptr;
    }

    double* dst = (double*)md_temp_alloc(temp, sizeof(double) * num_ao);
    if (!dst) {
        return nullptr;
    }

    // f64, because that is what md_gto takes: the coefficients are double at this boundary to keep
    // the QM code's precision, and extracting them through floats would spend it here.
    if (md_attribute_extract_f64(dst, num_ao, &sys.attributes, attr, slice, md_unit_none()) != num_ao) {
        return nullptr;
    }

    if (out_num_ao) *out_num_ao = num_ao;
    return dst;
}

// ---------------------------------------------------------------------------
// Which system atom each basis atom is
// ---------------------------------------------------------------------------
// A QM calculation may cover only PART of a loaded system - a chromophore inside a protein - and
// then the basis' atom indices are its own and not the system's. That map was the last thing an
// evaluation still needed a loader for, so it is published beside the basis and read from there.
//
// Absent means the identity AND that the QM atoms are the whole system, which is what a plain
// standalone load is, so nothing publishes it in the common case and every consumer needs the same
// one line to handle both. A standalone load that appends an embedding's sites after the QM atoms
// publishes the identity explicitly: the QM atoms are then only part of the system.
//
// It is VIAMD that publishes it and not the reader, because the file alone cannot decide: the same
// h5 carries a local-to-global map whether it is opened on its own - where the map must NOT be
// applied, since the system IS the QM atoms - or against a larger system, where it must. Only the
// side holding both knows which, and that is here.
static const md_str_t QM_ATOM_MAP_PATH = MD_STR_INIT("qm/atom/system_index");

// The map itself is written by the READER, on the load path: only the entry point that was called
// knows whether this file stands alone or supplements a larger system, and that is the whole of
// the decision. See vlx_publish_atom_system_index in md_vlx.c.

bool es_orbital_extent(const md_system_t& sys, size_t* out_num_mo, size_t* out_num_ao) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::alpha_coefficient);
    if (!attr || attr->format.rank != 2 || attr->format.components != 1) {
        return false;
    }
    if (out_num_mo) *out_num_mo = attr->format.shape[0];
    if (out_num_ao) *out_num_ao = attr->format.shape[1];
    return true;
}

bool es_orbital_frontier(OrbitalFrontier* out, const md_system_t& sys, md_str_t occupation_path) {
    MD_ASSERT(out);
    *out = {};

    size_t count = 0;
    const double* occ = attribute_view<double>(sys, occupation_path).get(&count);
    if (!occ) {
        return false;
    }

    const int num = (int)count;
    out->num_orbitals = num;

    // Scanning DOWN from the top rather than up to the first zero: an unoccupied orbital below an
    // occupied one is a non-aufbau ordering, not the frontier, and the last occupied orbital is
    // what the word means either way.
    int homo = -1;
    for (int i = num - 1; i >= 0; --i) {
        if (occ[i] > 0.0) {
            homo = i;
            break;
        }
    }
    out->homo_idx = homo;
    out->lumo_idx = homo + 1;   // == num when everything is occupied: out of range, and correct
    return true;
}

const double* es_orbital_energies(size_t* out_count, const md_system_t& sys, md_str_t energy_path) {
    return attribute_view<double>(sys, energy_path).get(out_count);
}

// NULL when there is none, which the callers read as the identity. Same explicitness as
// gto_attr_column in md_gto.c and for the same reason: the table is open, so a path with the wrong
// shape is a mistake to refuse rather than to reinterpret.
static const uint32_t* qm_atom_map_find(size_t* out_count, const md_system_t& sys) {
    size_t count = 0;
    const uint32_t* map = attribute_view<uint32_t>(sys, QM_ATOM_MAP_PATH).get(&count);
    if (!map && attribute_exists(sys, QM_ATOM_MAP_PATH)) {
        MD_LOG_ERROR("'" MD_STR_FMT "' is not a plain column of atom indices", MD_STR_ARG(QM_ATOM_MAP_PATH));
    }
    if (map && out_count) *out_count = count;
    return map;
}

bool es_has_distinct_beta_orbitals(const md_system_t& sys) {
    const md_attribute_t* a = md_attributes_find(&sys.attributes, es_path::alpha_coefficient);
    const md_attribute_t* b = md_attributes_find(&sys.attributes, es_path::beta_coefficient);
    if (!a || !b) {
        return false;
    }
    // md_attribute_same_data, never an id comparison: an alias has its own path and therefore its
    // own id, and this is the API's own answer to "are these two the same datum".
    return !md_attribute_same_data(a, b);
}

double es_electron_count(const md_system_t& sys, md_str_t occupation_path) {
    size_t count = 0;
    const double* occ = attribute_view<double>(sys, occupation_path).get(&count);
    if (!occ) {
        return 0.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < count; ++i) {
        sum += occ[i];
    }
    return sum;
}

bool es_has_distinct_beta_occupations(const md_system_t& sys) {
    const md_attribute_t* a = md_attributes_find(&sys.attributes, es_path::alpha_occupation);
    const md_attribute_t* b = md_attributes_find(&sys.attributes, es_path::beta_occupation);
    if (!a || !b) {
        return false;
    }
    return !md_attribute_same_data(a, b);
}

size_t es_excited_state_count(const md_system_t& sys) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::nto_lambda);
    if (!attr || attr->format.rank != 2 || attr->format.components != 1) {
        return 0;
    }
    return attr->format.shape[0];
}

// The states' NTOs read once on a worker thread: what es_nto_ready keeps. A state's status is written by its task
// and read by the frame. A task reads the system, so it ends before the system goes: free_system_data waits for
// every pool task (interrupt_async_tasks) before this is freed. The task is detached, so nothing cancels it before
// it has run (which would leave its state under way for good); it is listed in the overlay while it runs.
struct NtoWarm {
    enum : uint8_t { NotBegun = 0, UnderWay = 1, Done = 2 };
    size_t count = 0;
    std::atomic<uint8_t>* status = nullptr;
};

bool es_nto_ready(ApplicationState* state, size_t state_idx) {
    MD_ASSERT(state);
    const md_system_t& sys = state->mold.sys;
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::nto_lambda);
    const size_t count = es_excited_state_count(sys);
    if (!attr || attr->data || state_idx >= count) {
        return true;
    }
    NtoWarm* w = (NtoWarm*)state->mold.nto_warm;
    if (!w) {
        w = new NtoWarm();
        w->count = count;
        w->status = new std::atomic<uint8_t>[count];
        for (size_t i = 0; i < count; ++i) w->status[i].store(NtoWarm::NotBegun);
        state->mold.nto_warm = w;
    }
    if (state_idx >= w->count) {
        return true;
    }
    const uint8_t status = w->status[state_idx].load(std::memory_order_acquire);
    if (status == NtoWarm::Done) {
        return true;
    }
    if (status == NtoWarm::NotBegun) {
        w->status[state_idx].store(NtoWarm::UnderWay, std::memory_order_relaxed);
        const md_system_t* sys_ptr = &state->mold.sys;
        const uint32_t idx = (uint32_t)state_idx;
        task_system::detach(MD_STR_LIT("Natural transition orbitals: reading"), [w, sys_ptr, idx]() {
            // the weights: the reader computes the state's orbitals with them, and keeps all three
            const md_attribute_t* a = md_attributes_find(&sys_ptr->attributes, es_path::nto_lambda);
            if (a) {
                const md_attribute_slice_t slice = md_attribute_slice_1(idx);
                const size_t n = md_attribute_slice_count(a, slice);
                md_temp_scope_t temp = md_temp_begin();
                double* values = (double*)md_temp_alloc(temp, sizeof(double) * (n ? n : 1));
                if (values) md_attribute_extract_f64(values, n, &sys_ptr->attributes, a, slice, md_unit_none());
                md_temp_end(temp);
            }
            // done whether or not it could be read: where it is read next, it fails as loudly as it would have
            w->status[idx].store(NtoWarm::Done, std::memory_order_release);
            md_temp_trim();
        });
    }
    return false;
}

DensityCache* density_cache_of(ApplicationState* state) {
    MD_ASSERT(state);
    if (!state->mold.density_cache) state->mold.density_cache = density_cache_create();
    return state->mold.density_cache;
}

void es_nto_warm_free(SystemContext& sys) {
    NtoWarm* w = (NtoWarm*)sys.data.nto_warm;
    if (!w) return;
    delete[] w->status;
    delete w;
    sys.data.nto_warm = nullptr;
}

size_t es_nto_lambdas(double* out_values, size_t cap, const md_system_t& sys, size_t state_idx, double cutoff) {
    if (!out_values || cap == 0) {
        return 0;
    }
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::nto_lambda);
    if (!attr || state_idx >= es_excited_state_count(sys)) {
        return 0;
    }

    const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)state_idx);
    const size_t row = md_attribute_slice_count(attr, slice);
    if (row == 0) {
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    double* values = (double*)md_temp_alloc(temp, sizeof(double) * row);
    if (!values || md_attribute_extract_f64(values, row, &sys.attributes, attr, slice, md_unit_none()) != row) {
        return 0;
    }

    size_t count = 0;
    while (count < row && count < cap && values[count] >= cutoff) {
        out_values[count] = values[count];
        count += 1;
    }
    return count;
}

bool es_qm_atoms(QmAtoms* out, const md_system_t& sys) {
    MD_ASSERT(out);
    *out = {};

    // The element column defines the domain: it is published for every QM system, and its length
    // is the atom count. Without it there is no QM atom space to speak of.
    const md_attribute_t* z = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/atomic_number"));
    if (!z || !z->data || z->format.type != MD_ATTRIBUTE_TYPE_U8 ||
        z->format.rank != 1 || z->format.components != 1) {
        return false;
    }
    out->count         = z->format.shape[0];
    out->atomic_number = (const uint8_t*)z->data;

    if (const md_attribute_t* c = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"))) {
        if (c->data && c->format.type == MD_ATTRIBUTE_TYPE_F64 && c->format.rank == 1 &&
            c->format.components == 3 && c->format.shape[0] == out->count) {
            // md_double3_t is three doubles with no padding, and the {N,3} layout is exactly an
            // array of them, so this is a view and not a reinterpretation.
            out->coordinate = (const md_double3_t*)c->data;
        } else {
            MD_LOG_ERROR("'qm/atom/coordinate' is not %zu positions", out->count);
        }
    }

    size_t map_count = 0;
    if (const uint32_t* map = qm_atom_map_find(&map_count, sys)) {
        if (map_count == out->count) {
            out->system_index = map;
        } else {
            // Left null rather than half applied. The identity is wrong here, but sending an
            // evaluation to arbitrary atoms is worse than sending it to the first N.
            MD_LOG_ERROR("The QM domain holds %zu atoms and the system index maps %zu", out->count, map_count);
        }
    }
    return true;
}

// Gathers the BASIS atoms' positions out of the state, in the Bohr and the interleaved layout
// md_gto works in, in basis atom order. Returns the number written, 0 on failure.
//
// This conversion is the one input to an evaluation which is neither an attribute nor a caller
// parameter, and that is the right shape: the basis deliberately stores no coordinates, so that it
// survives a geometry change and the positions come from wherever the current ones live.
static size_t basis_atom_positions_gather(md_float3_t* dst, size_t cap, const md_system_t& sys, const md_system_state_t& state, size_t num_basis_atoms) {
    if (!dst || num_basis_atoms == 0 || num_basis_atoms > cap) {
        return 0;
    }
    // Checked here rather than at each call site: this is the only place the state is read, so this
    // is the only place that can be wrong about it.
    if (state.num_atoms == 0 || !state.xyz) {
        return 0;
    }

    size_t map_count = 0;
    const uint32_t* map = qm_atom_map_find(&map_count, sys);
    if (map && map_count < num_basis_atoms) {
        MD_LOG_ERROR("The basis spans %zu atoms and '" MD_STR_FMT "' maps %zu", num_basis_atoms, MD_STR_ARG(QM_ATOM_MAP_PATH), map_count);
        return 0;
    }

    for (size_t i = 0; i < num_basis_atoms; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= state.num_atoms) {
            MD_LOG_ERROR("Basis atom %zu is system atom %zu and the state holds %zu", i, idx, state.num_atoms);
            return 0;
        }
        dst[i] = md_float3_set(state.xyz[idx].x, state.xyz[idx].y, state.xyz[idx].z) * (float)ANGSTROM_TO_BOHR;
    }
    return num_basis_atoms;
}

// ---------------------------------------------------------------------------
// Volume readback
// ---------------------------------------------------------------------------
// Getting an evaluated volume out of the device scratch texture and into the GL texture a
// representation draws. Nothing about it is specific to what was evaluated or to who asked, and
// every handle it touches is on the ApplicationState above - which is why it lives here and not in
// whichever component happened to need it first.

#if MD_ENABLE_GPU

// Prefers a pixel unpack buffer: one write into memory the GPU already sees, and the transfer
// overlaps instead of blocking. Falls back to the plain client pointer upload when no buffer is
// available.
static void gpu_volume_upload_to_gl(uint32_t vol_tex, const void* src, size_t size) {
    if (!src) return;
    // The isosurface renderer's empty space grid follows the texels
    volume::notify_data_changed(vol_tex);
    if (void* dst = gl::pbo_upload_begin(size)) {
        md_memcpy(dst, src, size);
        if (gl::pbo_upload_end_texture_3D(vol_tex, 0, GL_R32F)) return;
    }
    gl::set_texture_3D_data(vol_tex, 0, src, GL_R32F);
}

// The scratch the evaluations write now (ApplicationState::gpu_eval_target)
static inline md_gpu_texture_t gpu_eval_target(const ApplicationState* state) {
    return state->gpu_eval_target ? state->gpu_eval_target : state->gpu_volume;
}

// Which normals a density volume's texels came with: set when both are uploaded (gpu_volume_job_complete), at the
// volume's data version then. Any later upload into the volume (notify_data_changed) leaves the link behind.
struct NormalLink {
    uint32_t vol_tex = 0, normal_tex = 0;
    uint64_t version = 0;
};
static NormalLink g_normal_links[32];

static void normal_link_set(uint32_t vol_tex, uint32_t normal_tex) {
    NormalLink* slot = &g_normal_links[0];
    for (NormalLink& l : g_normal_links) {
        if (l.vol_tex == vol_tex) { slot = &l; break; }
        if (l.version < slot->version) slot = &l;     // else the stalest
    }
    *slot = {vol_tex, normal_tex, volume::data_version(vol_tex)};
}

// The gradient scratch for grid (ApplicationState::gpu_grad): made, or grown, to hold it; false past GPU_GRAD_MAX_VOXELS
static bool gpu_grad_ensure(ApplicationState* state, const md_grid_t& grid) {
    if ((size_t)grid.dim[0] * grid.dim[1] * grid.dim[2] > GPU_GRAD_MAX_VOXELS) return false;
    if (state->gpu_grad && grid.dim[0] <= state->gpu_grad_dim[0] && grid.dim[1] <= state->gpu_grad_dim[1] && grid.dim[2] <= state->gpu_grad_dim[2]) return true;
    int dim[3];
    for (int i = 0; i < 3; ++i) dim[i] = MD_MAX(grid.dim[i], state->gpu_grad_dim[i]);
    if ((size_t)dim[0] * dim[1] * dim[2] > GPU_GRAD_MAX_VOXELS) {
        for (int i = 0; i < 3; ++i) dim[i] = grid.dim[i];
    }
    if (state->gpu_grad) md_gpu_texture_destroy(state->gpu_grad);
    if (state->gpu_normals) md_gpu_texture_destroy(state->gpu_normals);
    md_gpu_texture_desc_t desc = {
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_RGBA32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)dim[0], .height = (uint32_t)dim[1], .depth_or_layers = (uint32_t)dim[2], .label = "Evaluation gradient",
    };
    state->gpu_grad = md_gpu_texture_create(state->gpu_stream, &desc);
    desc.format = MD_GPU_FORMAT_R32_UINT;
    desc.label  = "Evaluation normals";
    state->gpu_normals = md_gpu_texture_create(state->gpu_stream, &desc);
    if (!state->gpu_grad || !state->gpu_normals) {
        MD_LOG_ERROR("Failed to create the GPU gradient scratch (%dx%dx%d): %s", dim[0], dim[1], dim[2], md_gpu_last_error());
        if (state->gpu_grad) md_gpu_texture_destroy(state->gpu_grad);
        if (state->gpu_normals) md_gpu_texture_destroy(state->gpu_normals);
        state->gpu_grad = state->gpu_normals = nullptr;
        for (int i = 0; i < 3; ++i) state->gpu_grad_dim[i] = 0;
        return false;
    }
    for (int i = 0; i < 3; ++i) state->gpu_grad_dim[i] = dim[i];
    return true;
}

// Where an evaluation launched now writes its gradient: gpu_grad when normals are asked for (gpu_eval_normals) of an
// evaluation into gpu_volume, noted for the readback (gpu_grad_written); else none
static md_gpu_texture_t gpu_grad_target(ApplicationState* state) {
    if (!state->gpu_eval_normals || !state->gpu_grad || gpu_eval_target(state) != state->gpu_volume) return nullptr;
    state->gpu_grad_written = true;
    return state->gpu_grad;
}

// The normals texture a readback queued now fills: the one asked for, when the evaluation wrote its gradient
static uint32_t gpu_normals_take(ApplicationState* state) {
    const uint32_t tex = state->gpu_grad_written ? state->gpu_eval_normals : 0;
    state->gpu_grad_written = false;
    return tex;
}

static ApplicationState::GpuVolumeJob* gpu_volume_job_acquire(ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (!state->gpu_volume_jobs[i].in_flight) {
            state->gpu_volume_jobs[i] = ApplicationState::GpuVolumeJob{};
            state->gpu_volume_jobs[i].in_flight = true;
            state->gpu_volume_jobs[i].owner = state;
            return &state->gpu_volume_jobs[i];
        }
    }
    return nullptr;
}

static bool gpu_volume_job_any_in_flight(const ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (state->gpu_volume_jobs[i].in_flight) return true;
    }
    return false;
}

static bool gpu_volume_job_in_flight_for(const ApplicationState* state, uint32_t tex_id) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (state->gpu_volume_jobs[i].in_flight && state->gpu_volume_jobs[i].tex_id == tex_id) return true;
    }
    return false;
}

// Runs on the GL thread, from md_gpu_device_poll() in the frame loop.
static void gpu_volume_job_complete(void* user) {
    ApplicationState::GpuVolumeJob* job = (ApplicationState::GpuVolumeJob*)user;
    ApplicationState* self = job->owner;
    if (self && job->tex_id) {
        gpu_volume_upload_to_gl(job->tex_id, job->rb.cpu, job->size);
        if (job->normal_tex && job->rb_normals.cpu) {
            gl::set_texture_3D_data(job->normal_tex, 0, job->rb_normals.cpu, GL_RGB10_A2);
            normal_link_set(job->tex_id, job->normal_tex);
        }
    }
    if (self && job->rb.gpu) md_gpu_free(self->gpu_stream, job->rb.gpu);
    if (self && job->rb_normals.gpu) md_gpu_free(self->gpu_stream, job->rb_normals.gpu);
    job->rb_normals = {};
    if (job->issued_job) {
        MD_LOG_DEBUG("Density volume %dx%dx%d: %.1f ms from the request to the texture, on the GPU in pieces over the frames (" MD_GTO_GPU_GEMM_PIECES_REVISION ")",
                     (int)job->grid.dim[0], (int)job->grid.dim[1], (int)job->grid.dim[2], md_tick_to_milliseconds(md_tick_now() - (md_tick_t)job->t0));
        md_gto_gpu_job_end((md_gto_gpu_job_t*)job->issued_job);
        job->issued_job = nullptr;
    }
    job->rb        = {};
    job->in_flight = false;
}

// Used when no job slot is free, and when a caller genuinely needs the data before it returns.
static bool gpu_volume_readback_blocking(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, size_t size, uint32_t normal_tex) {
    md_gpu_mem_t rb = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    const md_gpu_tex_region_t region = {
        .offset = {0, 0, 0},
        .extent = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] },
    };
    bool ok = md_gpu_copy_from_texture(state->gpu_stream, rb.gpu, gpu_eval_target(state), &region);
    md_gpu_mem_t rbn = {};
    if (ok && normal_tex) {
        rbn = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_HOST_READ, size);    // 4 bytes a voxel, as the values
        md_gto_gpu_gradient_pack(state->gpu_stream, state->gpu_grad, state->gpu_normals, &grid);
        if (!rbn.cpu || !md_gpu_copy_from_texture(state->gpu_stream, rbn.gpu, state->gpu_normals, &region)) normal_tex = 0;
    }
    md_gpu_stream_sync(state->gpu_stream);
    if (ok) {
        gpu_volume_upload_to_gl(vol_tex, rb.cpu, size);
        if (normal_tex) {
            gl::set_texture_3D_data(normal_tex, 0, rbn.cpu, GL_RGB10_A2);
            normal_link_set(vol_tex, normal_tex);
        }
    }
    md_gpu_free(state->gpu_stream, rb.gpu);
    if (rbn.gpu) md_gpu_free(state->gpu_stream, rbn.gpu);
    return ok;
}

static bool gpu_volume_readback_issue(ApplicationState* state, ApplicationState::GpuVolumeJob* job, uint32_t vol_tex, const md_grid_t& grid);

// The density jobs not yet all issued: each, once its counts are there (wait: waited for), issued and its
// readback queued. True when none is left.
static bool gpu_density_jobs_issue(ApplicationState* state, bool wait) {
    bool none = true;
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        if (!j.in_flight || !j.density_job) continue;
        md_gto_gpu_job_t* dj = (md_gto_gpu_job_t*)j.density_job;
        if (!md_gto_gpu_job_poll(dj, wait)) {
            none = false;
            continue;
        }
        j.density_job = nullptr;
        j.issued_job = dj;      // kept for its progress until the readback lands (gpu_volume_job_complete)
        gpu_volume_readback_issue(state, &j, j.tex_id, j.grid);
        if (!j.in_flight && j.issued_job) {     // failed, or read back on the spot
            md_gto_gpu_job_end(dj);
            j.issued_job = nullptr;
        }
    }
    return none;
}

static bool gpu_density_job_pending(const ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (state->gpu_volume_jobs[i].in_flight && state->gpu_volume_jobs[i].density_job) return true;
    }
    return false;
}

// What the GPU evaluations share - the atom positions, the coefficients, the dmat, the scratch volume - is
// only changed once every density job is issued: one not yet issued reads them when it is. Waits for those.
static void gpu_density_jobs_settle(ApplicationState* state) {
    if (gpu_density_job_pending(state)) gpu_density_jobs_issue(state, true);
}

// The density jobs not yet all issued, of vol_tex (0: every one), dropped: what they have in flight completes,
// nothing more is issued, and no readback is queued (the texture keeps what it held). True when any was.
static bool gpu_density_jobs_cancel(ApplicationState* state, uint32_t vol_tex) {
    bool any = false;
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        if (!j.in_flight || !j.density_job || (vol_tex && j.tex_id != vol_tex)) continue;
        md_gto_gpu_job_t* dj = (md_gto_gpu_job_t*)j.density_job;
        MD_LOG_DEBUG("Density volume: the evaluation under way no longer wanted, dropped at %.0f%% after %.1f ms", md_gto_gpu_job_progress(dj) * 100.0f,
                     md_tick_to_milliseconds(md_tick_now() - (md_tick_t)j.t0));
        md_gto_gpu_job_end(dj);
        j = ApplicationState::GpuVolumeJob{};
        any = true;
    }
    return any;
}

// Queues a readback of the evaluated region of gpu_volume. Returns once the copy is recorded; the
// GL texture is filled later, from gpu_volume_job_complete(), which md_gpu_device_poll() calls on
// the GL thread.
//
// Returns true when the work was QUEUED -- not that the texture holds data.
static bool gpu_volume_readback(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid) {
    if (!state->gpu_stream || !gpu_eval_target(state)) {
        return false;
    }
    const size_t size = sizeof(float) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    const uint32_t normal_tex = gpu_normals_take(state);

    // Never allow two outstanding readbacks for the same texture. Two reasons, both of which
    // produce a wrong image rather than a slow one:
    //
    //  - a staging block can be most of a gigabyte, so N in flight is N times that;
    //  - md_gpu_stream_sync does not fire user callbacks, only md_gpu_device_poll does. So a
    //    blocking fallback taken while an older job is still pending would upload new data now and
    //    let the older callback overwrite it with stale data next frame.
    //
    // Draining costs the stall we are trying to avoid, but only when the user outruns the GPU, and
    // it keeps uploads strictly ordered.
    if (gpu_volume_job_in_flight_for(state, vol_tex)) {
        gpu_volume_jobs_drain(state);
    }

    ApplicationState::GpuVolumeJob* job = gpu_volume_job_acquire(state);
    if (!job) {
        gpu_volume_jobs_drain(state);
        job = gpu_volume_job_acquire(state);
    }
    if (!job) {
        return gpu_volume_readback_blocking(state, vol_tex, grid, size, normal_tex);
    }
    job->src = gpu_eval_target(state);
    job->normal_tex = normal_tex;
    return gpu_volume_readback_issue(state, job, vol_tex, grid);
}

// The copy of gpu_volume into job's staging and its completion, queued; job is acquired
static bool gpu_volume_readback_issue(ApplicationState* state, ApplicationState::GpuVolumeJob* job, uint32_t vol_tex, const md_grid_t& grid) {
    const size_t size = sizeof(float) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    md_gpu_mem_t rb = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) {
        MD_LOG_ERROR("Failed to allocate volume readback staging (%zu bytes)", size);
        job->in_flight = false;
        return false;
    }

    const md_gpu_tex_region_t region = {
        .offset = {0, 0, 0},
        .extent = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] },
    };
    if (!md_gpu_copy_from_texture(state->gpu_stream, rb.gpu, job->src ? job->src : state->gpu_volume, &region)) {
        MD_LOG_ERROR("Failed to record the volume readback");
        md_gpu_free(state->gpu_stream, rb.gpu);
        job->in_flight = false;
        return false;
    }
    // the normals of its gradient beside it, packed on the device: 4 bytes a voxel
    if (job->normal_tex) {
        md_gto_gpu_gradient_pack(state->gpu_stream, state->gpu_grad, state->gpu_normals, &grid);
        job->rb_normals = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_HOST_READ, size);
        if (!job->rb_normals.cpu || !md_gpu_copy_from_texture(state->gpu_stream, job->rb_normals.gpu, state->gpu_normals, &region)) {
            if (job->rb_normals.gpu) md_gpu_free(state->gpu_stream, job->rb_normals.gpu);
            job->rb_normals = {};
            job->normal_tex = 0;
        }
    }

    job->tex_id = vol_tex;
    job->rb     = rb;
    job->size   = size;

    if (!md_gpu_sync_on_complete(state->gpu_device, md_gpu_stream_record(state->gpu_stream), gpu_volume_job_complete, job)) {
        // No callback means nothing would ever free rb or fill the texture, so finish this one
        // synchronously instead of leaking it.
        MD_LOG_ERROR("Failed to queue the volume completion; falling back to a blocking readback");
        md_gpu_stream_sync(state->gpu_stream);
        gpu_volume_upload_to_gl(vol_tex, rb.cpu, size);
        if (job->normal_tex && job->rb_normals.cpu) {
            gl::set_texture_3D_data(job->normal_tex, 0, job->rb_normals.cpu, GL_RGB10_A2);
            normal_link_set(vol_tex, job->normal_tex);
        }
        md_gpu_free(state->gpu_stream, rb.gpu);
        if (job->rb_normals.gpu) md_gpu_free(state->gpu_stream, job->rb_normals.gpu);
        job->rb_normals = {};
        job->in_flight = false;
        return true;
    }
    return true;
}

// The readback of a density job's evaluation (a density, ELF or LOL: md_gto_gpu_job_t): its slot now, the copy once the
// job is all issued (gpu_volume_jobs_poll). Takes the job. label: what it is, for the overlay (a static string).
static bool gpu_volume_readback_after(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, md_gto_gpu_job_t* dj, const char* label) {
    if (md_gto_gpu_job_poll(dj, false)) {
        md_gto_gpu_job_end(dj);
        return gpu_volume_readback(state, vol_tex, grid);
    }
    if (gpu_volume_job_in_flight_for(state, vol_tex)) {
        // not to be: the representations wait for the last one (gpu_volume_job_in_flight_for). Ordered as the
        // plain readback is
        md_gto_gpu_job_poll(dj, true);
        md_gto_gpu_job_end(dj);
        return gpu_volume_readback(state, vol_tex, grid);
    }
    ApplicationState::GpuVolumeJob* job = gpu_volume_job_acquire(state);
    if (!job) {
        md_gto_gpu_job_poll(dj, true);
        md_gto_gpu_job_end(dj);
        return gpu_volume_readback(state, vol_tex, grid);
    }
    job->tex_id = vol_tex;
    job->src = gpu_eval_target(state);
    job->normal_tex = gpu_normals_take(state);
    job->grid = grid;
    job->density_job = dj;
    job->t0 = md_tick_now();
    job->label = label;
    return true;
}

// Queues the atom upload on gpu_stream if it is dirty. Any evaluation launched afterwards on the
// same stream observes it. md_gpu_upload_begin writes straight into the destination when that is
// safe and stages through a transient arena otherwise, so there is one path regardless of whether
// the device is discrete.
static bool gpu_atoms_ensure_uploaded(ApplicationState* state, const md_system_t& sys, const md_system_state_t& sys_state) {
    if (!state->gpu_stream || !state->mold.gpu_basis || !state->mold.gpu_atoms) {
        return false;
    }
    const size_t num_atoms = md_gto_gpu_basis_num_atoms(state->mold.gpu_basis);
    const size_t sz = md_gto_gpu_atom_buffer_size(num_atoms);

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    // Same gather as the GL path, through the same map, so the two backends cannot disagree about
    // which atom a basis index means. Packed to vec4 because that is what the device wants.
    md_float3_t* xyz = (md_float3_t*)md_temp_alloc(temp, sizeof(md_float3_t) * MD_MAX(num_atoms, (size_t)1));
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, sys_state, num_atoms) != num_atoms) {
        return false;
    }

    // Gather first, then decide: the gather is what the comparison is ABOUT, and it is cheap next
    // to the transfer it may save.
    const uint64_t hash = md_hash64(xyz, sizeof(md_float3_t) * num_atoms, 0);
    if (hash != 0 && hash == state->mold.gpu_atoms_hash) {
        return true;
    }

    md_float4_t* xyzw = (md_float4_t*)md_temp_alloc(temp, sizeof(md_float4_t) * MD_MAX(num_atoms, (size_t)1));
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = md_float4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
    }

    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->mold.gpu_atoms, sz);
    if (!dst) {
        MD_LOG_ERROR("Failed to upload the atom positions to the device");
        return false;
    }
    md_gto_gpu_atom_pack(dst, (const float*)xyzw, sizeof(md_float4_t), num_atoms);
    if (!md_gpu_upload_end(state->gpu_stream)) {
        return false;
    }

    state->mold.gpu_atoms_hash = hash;
    return true;
}

// The volumes' evaluations in the overlay: a slot's activity while its job runs (issued or not) and until its readback
// lands, under the label it was given; not cancellable there (a representation's volume would keep the hash of what it
// never got)
static task_system::Activity g_volume_activity[ApplicationState::GPU_VOLUME_JOB_SLOTS];
static const char*           g_volume_activity_label[ApplicationState::GPU_VOLUME_JOB_SLOTS];

static void gpu_volume_activities_update(const ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        const ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        void* dj = j.density_job ? j.density_job : j.issued_job;
        task_system::Activity& a = g_volume_activity[i];
        const char* label = j.in_flight && dj ? (j.label ? j.label : "Volume") : nullptr;
        if (label != g_volume_activity_label[i]) {
            task_system::activity_end(&a);
            g_volume_activity_label[i] = label;
            if (label) {
                a.cancellable = false;
                a.progress.store(0.0f, std::memory_order_relaxed);
                task_system::activity_begin(&a, md_str_from_cstr(label));
            }
        }
        if (label) a.progress.store(md_gto_gpu_job_progress((md_gto_gpu_job_t*)dj), std::memory_order_relaxed);
    }
}

#endif // MD_ENABLE_GPU

uint32_t es_normals_of(uint32_t vol_tex) {
#if MD_ENABLE_GPU
    if (!vol_tex) return 0;
    for (const NormalLink& l : g_normal_links) {
        if (l.vol_tex == vol_tex) return l.version == volume::data_version(vol_tex) ? l.normal_tex : 0;
    }
#else
    (void)vol_tex;
#endif
    return 0;
}

void gpu_volume_jobs_poll(ApplicationState* state) {
    MD_ASSERT(state);
#if MD_ENABLE_GPU
    if (state->gpu_stream) gpu_density_jobs_issue(state, false);
    gpu_volume_activities_update(state);
#else
    (void)state;
#endif
}

bool gpu_volume_job_cancel(ApplicationState* state, uint32_t vol_tex) {
    MD_ASSERT(state);
#if MD_ENABLE_GPU
    return vol_tex != 0 && state->gpu_stream && gpu_density_jobs_cancel(state, vol_tex);
#else
    (void)vol_tex;
    return false;
#endif
}

void gpu_volume_jobs_drain(ApplicationState* state) {
    MD_ASSERT(state);
#if MD_ENABLE_GPU
    if (!state->gpu_stream || !gpu_volume_job_any_in_flight(state)) return;
    gpu_density_jobs_issue(state, true);
    md_gpu_stream_sync(state->gpu_stream);
    md_gpu_device_poll(state->gpu_device);   // this is what actually runs the callbacks

    // Backstop: if a callback somehow did not fire, release the staging block here rather than
    // leaking it, since the pool may be about to go.
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        if (j.in_flight) {
            if (j.density_job) md_gto_gpu_job_end((md_gto_gpu_job_t*)j.density_job);
            if (j.issued_job) md_gto_gpu_job_end((md_gto_gpu_job_t*)j.issued_job);
            j.density_job = nullptr;
            j.issued_job = nullptr;
            if (j.rb.gpu) md_gpu_free(state->gpu_stream, j.rb.gpu);
            j.rb = {};
            j.in_flight = false;
        }
    }
#else
    (void)state;
#endif
}

md_float3_t* basis_atom_positions_extract(size_t* out_count, md_temp_scope_t temp, const md_system_t& sys, const md_system_state_t& state) {
    if (out_count) *out_count = 0;

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return nullptr;
    }

    const size_t num_basis_atoms = md_gto_basis_num_atoms(&basis);
    md_float3_t* atom_xyz = (md_float3_t*)md_temp_alloc(temp, sizeof(md_float3_t) * MD_MAX(num_basis_atoms, (size_t)1));
    if (!atom_xyz || basis_atom_positions_gather(atom_xyz, num_basis_atoms, sys, state, num_basis_atoms) != num_basis_atoms) {
        return nullptr;
    }

    if (out_count) *out_count = num_basis_atoms;
    return atom_xyz;
}

// ---------------------------------------------------------------------------
// Density evaluation
// ---------------------------------------------------------------------------
// A density in this table is an AO x AO matrix. It is either that on its own - the SCF densities,
// the density properties - or the innermost two axes of something indexed by state, as the
// transition densities are at {S,A,A}. Both are addressed the same way here: a path plus a slice
// which narrows to one matrix. The caller names what it wants and this cares neither which of the
// two shapes it came from, nor whether the values were stored or worked out on demand, because
// md_attribute_slice_format answers the first and the extract answers the second.
// The context a density is extracted through: a density built from its orbitals on demand (md_qm's providers)
// spreads its tiles over the pool, at the priority of whoever asks. One per extraction, since a context is not
// to be used by two threads at once.
struct DensityIo {
    md_attribute_io_t* io = md_attribute_io_create(task_system::executor(), md_get_heap_allocator());
    DensityIo() = default;
    DensityIo(const DensityIo&) = delete;
    DensityIo& operator=(const DensityIo&) = delete;
    ~DensityIo() { md_attribute_io_destroy(io); }
};

double* density_matrix_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, md_str_t density_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" MD_STR_FMT "'", MD_STR_ARG(density_path));
        return nullptr;
    }

    // Square or packed, stored or computed on demand: md_qm reads either, and says why when the
    // slice does not narrow it to one symmetric matrix
    const size_t dim = md_qm_extract_symmetric_f64(nullptr, 0, &sys.attributes, attr, slice, nullptr);
    if (dim == 0) {
        return nullptr;
    }
    double* dst = (double*)md_temp_alloc(temp, sizeof(double) * dim * dim);
    DensityIo io;
    if (!dst || md_qm_extract_symmetric_f64(dst, dim * dim, &sys.attributes, attr, slice, io.io) != dim) {
        return nullptr;
    }

    if (out_dim) *out_dim = dim;
    return dst;
}

// The density a representation's volume is evaluated from, for the sources that are one: its path and slice
// (all of it: num_idx 0). False for the orbitals, and for a density property the system no longer has.
static bool es_volume_density(md_str_t* out_path, md_attribute_slice_t* out_slice, const ElectronicStructureRepresentation& es, const md_system_t& sys) {
    *out_slice = md_attribute_slice_all();
    switch (es.shown.source) {
    case ElectronicStructureSource::TransitionDensity: {
        md_str_t path = es_path::attachment_density;
        switch (es.shown.transition_density_component) {
        case ElectronicStructureTransitionDensityComponent::Detachment: path = es_path::detachment_density; break;
        case ElectronicStructureTransitionDensityComponent::Difference: path = es_path::transition_diff;    break;
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default: break;
        }
        // These are VIRTUAL: the slice is what tells the provider to reconstruct one state rather
        // than all of them, and it is the only reason asking for one is affordable.
        *out_path  = path;
        *out_slice = md_attribute_slice_1((uint32_t)es.shown.excited_state_idx);
        return true;
    }
    case ElectronicStructureSource::ElectronDensity:
        *out_path = es_electron_density_path(es.shown.spin);
        return true;
    case ElectronicStructureSource::DensityProperty: {
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, es.shown.density_property_key);
        if (!attr) return false;
        *out_path = attr->path;
        return true;
    }
    default:
        return false;
    }
}

// Evaluates whatever the representation is currently pointed at into its own volume texture.
// Everything it needs is an attribute on the system or a field of the representation; there is no
// reader, no component and no event in the path, which is the whole point of the port.
// ---------------------------------------------------------------------------
// ELF and LOL over a grid
// ---------------------------------------------------------------------------
// On the GPU by their own grid kernel (md_gto_gpu_localization_job_begin): the spins' density matrices in tiles,
// built on the device from the orbitals over a GPU basis at the AO cutoff 1e-7 (ELF is made of small differences
// in the tails), made by a pool task (localization_gpu_prepare) and kept until the orbitals change or the atoms
// move out of the tiles' margin; the kernel issued into the GPU stream as a job in pieces, as the density's is, and its
// volume read back as the density's (gpu_volume_readback_after): nothing waits. Without the kernel, as a surface field
// on the CPU, in double (surface_field.h, es_source_volume_poll); slow.
//
// No floor of rho. ELF is a ratio of two vanishing quantities far from the atoms (with a gaussian basis it goes to 1,
// or in float to noise, there: nothing of the molecule); it is taken as Multiwfn has it, D + MD_QM_ELF_D_MIN (1e-5),
// which follows the exact one to the 1e-3 isodensity and goes to 0 smoothly by rho 1e-4 (md_qm_elf). LOL is a plain
// ratio and goes to 0 there itself. The CPU path, in double from the same float data, alike.
#define ES_LOCALIZATION_RHO_MIN 0.0
#define ES_LOCALIZATION_CUTOFF 1.0e-7
// How far the atoms may move before the matrices' pattern is made again (bohr), as the density's (GPU_DMAT_MARGIN)
#define ES_LOCALIZATION_MARGIN 1.0

#if MD_ENABLE_GPU
static bool gpu_evaluation_ready(const ApplicationState* state);
static uint64_t density_dmat_key(const md_system_t& sys, md_str_t density_path);
#endif

#if MD_ENABLE_GPU
// What the kernel reads, made by a pool task on a stream of its own (destroyed by the task once its work is done:
// what it made is then ready for any stream), kept in state->mold.gpu_localization
struct LocalizationGpu {
    task_system::Job   task;            // a cancelled one (a load) made nothing to be taken
    uint64_t           key = 0;         // the orbitals and occupations it is made of (density_dmat_key)
    const md_system_t* sys = nullptr;
    md_float3_t*            xyz = nullptr;   // bohr, its copy: the positions the tiles' pattern is made for
    size_t             num_atoms = 0;
    md_gpu_stream_t    stream = nullptr;
    md_gto_basis_t     basis = {};      // heap
    md_gto_gpu_basis_t gpu_basis = nullptr;
    md_gto_gpu_dmat_t  dmat[2] = {};
    uint32_t           num_spins = 0;   // 1: dmat[0] the total density of a closed shell; 2: alpha and beta
    bool               ok = false;
    md_tick_t          t0 = 0;
};

// On a worker thread
static void localization_gpu_build(LocalizationGpu* L) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool ok = md_gto_basis_extract_attributes(&L->basis, &L->sys->attributes, heap);
    const md_gto_gpu_basis_desc_t bd = { .basis = &L->basis, .cutoff = ES_LOCALIZATION_CUTOFF };
    L->gpu_basis = ok ? md_gto_gpu_basis_create(L->stream, &bd) : nullptr;
    ok = ok && L->gpu_basis;
    // a closed shell: the same orbitals and occupations in both spins, which are then one density, the total
    double *ar = nullptr, *aw = nullptr, *br = nullptr, *bw = nullptr;
    size_t ana = 0, bna = 0;
    const size_t an = ok ? md_qm_density_orbitals(&ar, &aw, &ana, L->sys, es_path::alpha_density, alloc) : 0;
    const size_t bn = ok ? md_qm_density_orbitals(&br, &bw, &bna, L->sys, es_path::beta_density, alloc) : 0;
    const bool closed = an == bn && ana == bna && memcmp(ar, br, sizeof(double) * an * ana) == 0 && memcmp(aw, bw, sizeof(double) * an) == 0;
    ok = ok && an > 0 && bn > 0;
    L->num_spins = closed ? 1 : 2;
    for (uint32_t k = 0; ok && k < L->num_spins; ++k) {
        const md_gto_gpu_dmat_desc_t dd = {
            .basis = &L->basis, .atom_xyz = (const float*)L->xyz, .atom_xyz_stride = sizeof(md_float3_t), .cutoff = ES_LOCALIZATION_CUTOFF, .margin = ES_LOCALIZATION_MARGIN };
        L->dmat[k] = md_gto_gpu_dmat_create(L->stream, &dd);
        double* rows = k == 0 ? ar : br;
        double* w = k == 0 ? aw : bw;
        const size_t n = k == 0 ? an : bn;
        if (closed) for (size_t i = 0; i < n; ++i) w[i] *= 2.0;     // the total
        ok = L->dmat[k] && md_gto_gpu_dmat_build(L->stream, L->dmat[k], rows, w, n);
    }
    L->ok = ok;
    md_temp_end(temp);
    md_gpu_stream_destroy(L->stream);   // once its work is done
    L->stream = nullptr;
    md_temp_trim();
}

static void localization_gpu_destroy(LocalizationGpu* L, md_gpu_stream_t gl_stream) {
    // only once its task has ended
    L->task.cancel_and_wait();
    for (int k = 0; k < 2; ++k) if (L->dmat[k]) md_gto_gpu_dmat_destroy(gl_stream, L->dmat[k]);
    if (L->gpu_basis) md_gto_gpu_basis_destroy(gl_stream, L->gpu_basis);
    if (L->stream) md_gpu_stream_destroy(L->stream);    // of a task cancelled before it ran
    md_gto_basis_free(&L->basis, md_get_heap_allocator());
    free(L->xyz);
    delete L;
}

static void localization_gpu_free(ApplicationState* state) {
    if (state->mold.gpu_localization) {
        localization_gpu_destroy((LocalizationGpu*)state->mold.gpu_localization, state->gpu_stream);
        state->mold.gpu_localization = nullptr;
    }
    state->mold.gpu_localization_failed_key = 0;
}

// Whether ELF and LOL can be evaluated on the GPU now: 1 when what the kernel reads is there for the orbitals and
// the atoms as they are, 0 while it is being made (begun here when it is not), -1 when there is no GPU or it
// could not be made (then on the CPU)
static int localization_gpu_prepare(ApplicationState* state) {
    if (!gpu_evaluation_ready(state) || !state->gpu_device) return -1;
    const md_system_t& sys = state->mold.sys;
    const uint64_t key = density_dmat_key(sys, MD_STR_LIT("localization"));
    if (key == state->mold.gpu_localization_failed_key) return -1;
    LocalizationGpu* L = (LocalizationGpu*)state->mold.gpu_localization;
    if (L && L->task.running()) return 0;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    size_t num_atoms = 0;
    const md_float3_t* xyz = basis_atom_positions_extract(&num_atoms, temp, sys, system_instance(state->system_context).state);
    if (!xyz) return -1;
    if (L && L->key == key) {
        if (L->task.status() != task_system::Status::Completed || !L->ok) {
            MD_LOG_ERROR("ELF and LOL: what their GPU kernel reads could not be made; evaluating them on the CPU");
            localization_gpu_free(state);
            state->mold.gpu_localization_failed_key = key;
            return -1;
        }
        bool valid = num_atoms == L->num_atoms;
        for (uint32_t k = 0; valid && k < L->num_spins; ++k) valid = md_gto_gpu_dmat_pattern_valid(L->dmat[k], (const float*)xyz, sizeof(md_float3_t));
        if (valid) return 1;
    }
    // made (again): for other orbitals, or atoms beyond the margin
    if (L) localization_gpu_free(state);
    L = new LocalizationGpu();
    L->key = key;
    L->sys = &state->mold.sys;
    L->t0 = md_tick_now();
    L->xyz = (md_float3_t*)malloc(sizeof(md_float3_t) * num_atoms);
    L->num_atoms = num_atoms;
    L->stream = L->xyz ? md_gpu_stream_create(state->gpu_device, MD_GPU_STREAM_COMPUTE, "ELF/LOL matrices") : nullptr;
    if (!L->stream) {
        free(L->xyz);
        delete L;
        return -1;
    }
    md_memcpy(L->xyz, xyz, sizeof(md_float3_t) * num_atoms);
    state->mold.gpu_localization = L;
    L->task.start(MD_STR_LIT("ELF/LOL density matrices"), [L]() { localization_gpu_build(L); });
    return 0;
}

// ELF or LOL into vol_tex by the kernel, as a job read back as the density's is: its counting pass issued and nothing
// waited for (the launch waited on the GL thread for the counts and whatever the stream held before them, then issued
// every piece at once: minutes of GPU work ahead of the display for a protein), the pieces issued a few ahead of the
// GPU over the frames (gpu_volume_jobs_poll). What it reads is there (localization_gpu_prepare).
static bool localization_gpu_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, bool lol) {
    LocalizationGpu* L = (LocalizationGpu*)state->mold.gpu_localization;
    if (!L || !L->ok || L->task.running()) return false;
    gpu_density_jobs_settle(state);     // the scratch volume is written below
    if (!gpu_atoms_ensure_uploaded(state, state->mold.sys, system_instance(state->system_context).state)) return false;
    md_gto_gpu_localization_desc_t desc = {
        .basis = L->gpu_basis, .atom_xyz = state->mold.gpu_atoms, .num_spins = L->num_spins,
        .dmat = { L->dmat[0], L->num_spins > 1 ? L->dmat[1] : nullptr },
        .kind = lol ? MD_GTO_GPU_LOL : MD_GTO_GPU_ELF, .rho_min = ES_LOCALIZATION_RHO_MIN,
        .d_min = MD_QM_ELF_D_MIN,
        .out_tex = gpu_eval_target(state), .grid = &grid, .sample_offset = {0.5f, 0.5f, 0.5f},
    };
    md_gto_gpu_job_t* job = md_gto_gpu_localization_job_begin(state->gpu_stream, &desc);
    if (!job) return false;
    return gpu_volume_readback_after(state, vol_tex, grid, job, lol ? "Localized orbital locator (LOL)" : "Electron localization function (ELF)");
}
#endif

static bool es_evaluate_into(ApplicationState* state, const ElectronicStructureRepresentation& es, uint32_t tex_id, const md_grid_t& grid);

// es: what the volume is of
static bool electronic_structure_evaluate(ApplicationState* state, Representation* rep, const ElectronicStructureRepresentation& es) {

    const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)electronic_structure_volume_resolution(es)];

    md_grid_t grid = {};
    if (!electronic_structure_grid_init(&grid, state->mold.sys, system_instance(state->system_context).state, samples_per_angstrom * BOHR_TO_ANGSTROM)) {
        return false;
    }
    init_volume(&rep->electronic_structure.density_vol, grid, GL_R32F);
    rep->electronic_structure.grid = grid;
    rep->electronic_structure.box_to_world = rep->electronic_structure.density_vol.texture_to_world;
    rep->electronic_structure.has_box = true;

    // Raycast isosurfaces are shaded with the analytic normals of what the GPU evaluates with its gradient (the
    // orbitals and densities: gpu_grad_target), where it fits the scratch; the others with the volume's differences
    bool normals = false;
#if MD_ENABLE_GPU
    normals = es.display == ElectronicStructureDisplay::Isosurface && !es_mesh_mode(state, es) && gpu_evaluation_ready(state) &&
              gpu_grad_ensure(state, grid);
    if (normals) gl::init_texture_3D(&rep->electronic_structure.normal_tex, grid.dim[0], grid.dim[1], grid.dim[2], GL_RGB10_A2);
    state->gpu_eval_normals = normals ? rep->electronic_structure.normal_tex : 0;
    state->gpu_grad_written = false;
#endif
    const bool ok = es_evaluate_into(state, es, rep->electronic_structure.density_vol.tex_id, grid);
#if MD_ENABLE_GPU
    state->gpu_eval_normals = 0;
    state->gpu_grad_written = false;
#endif
    (void)normals;
    return ok;
}

// What the representation's source is, evaluated over grid into tex_id (a GL 3D texture of the grid's dims): the
// volume, or the image of a plane (one layer). On the GPU, written to the scratch gpu_eval_target and read back.
static bool es_evaluate_into(ApplicationState* state, const ElectronicStructureRepresentation& es, uint32_t tex_id, const md_grid_t& grid) {
    switch (es.shown.source) {
    case ElectronicStructureSource::MolecularOrbital: {
        const md_str_t path = (es.shown.spin == ElectronicStructureSpin::Beta) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.shown.orbital_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.shown.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::NaturalTransitionOrbital: {
        const md_str_t path = (es.shown.nto_component == ElectronicStructureNtoComponent::Particle) ? es_path::nto_particle : es_path::nto_hole;
        // {S,L,A}: the excited state and then the lambda pair within it. Two indices instead of the
        // one an ordinary orbital takes, which is exactly what a slice is for.
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)es.shown.excited_state_idx, (uint32_t)es.shown.nto_lambda_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.shown.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::TransitionDensity:
    case ElectronicStructureSource::ElectronDensity:
    case ElectronicStructureSource::DensityProperty: {
        md_str_t path = {};
        md_attribute_slice_t slice = {};
        if (!es_volume_density(&path, &slice, es, state->mold.sys)) {
            return false;
        }
        // A spin difference is signed, so it is the one density the magnitude toggle applies to.
        const bool magnitude = es.shown.source == ElectronicStructureSource::ElectronDensity && es.shown.spin == ElectronicStructureSpin::Difference && es.shown.use_magnitude;
        return density_evaluate(state, tex_id, grid, path, slice.num_idx > 0 ? &slice : nullptr, es_gto_op(magnitude), true);
    }
    case ElectronicStructureSource::ElectronLocalization:
    case ElectronicStructureSource::LocalizedOrbitalLocator:
#if MD_ENABLE_GPU
        // by their kernel when what it reads is there (the update waits for it); else as a surface field on the CPU,
        // which fills the volume itself (es_source_volume_poll)
        return localization_gpu_evaluate(state, tex_id, grid, es.shown.source == ElectronicStructureSource::LocalizedOrbitalLocator);
#else
        return false;
#endif
    default:
        // the fields of the system as a whole: as surface fields, into the volume (es_source_volume_poll)
        return false;
    }
}

// The colour volume beside the density one: the atoms' own colours splatted into a downsampled 3D
// texture, so that a surface can be shaded by whatever the representation is colouring atoms with.
static void electronic_structure_color_volume_update(ApplicationState* state, Representation* rep, const uint32_t* atom_colors) {
    if (!atom_colors) {
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = system_instance(state->system_context).state;

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return;
    }
    const size_t num_points = md_gto_basis_num_atoms(&basis);
    if (num_points == 0) {
        return;
    }

    // Half the density's resolution: the colours vary on the scale of atoms, and at full resolution the
    // RGBA8 texture is as large as the density itself (512 MB at the 512^3 limit)
    const int downsample_factor = 2;
    int dim[3] = {
        MD_MAX(1, MD_DIV_UP((int)rep->electronic_structure.density_vol.dim[0], downsample_factor)),
        MD_MAX(1, MD_DIV_UP((int)rep->electronic_structure.density_vol.dim[1], downsample_factor)),
        MD_MAX(1, MD_DIV_UP((int)rep->electronic_structure.density_vol.dim[2], downsample_factor)),
    };
    md_memcpy(rep->electronic_structure.color_vol.dim, dim, sizeof(dim));
    rep->electronic_structure.color_vol.world_to_model   = rep->electronic_structure.density_vol.world_to_model;
    rep->electronic_structure.color_vol.texture_to_world = rep->electronic_structure.density_vol.texture_to_world;
    rep->electronic_structure.color_vol.voxel_size       = md_float3_set(
        rep->electronic_structure.density_vol.voxel_size.x * rep->electronic_structure.density_vol.dim[0] / dim[0],
        rep->electronic_structure.density_vol.voxel_size.y * rep->electronic_structure.density_vol.dim[1] / dim[1],
        rep->electronic_structure.density_vol.voxel_size.z * rep->electronic_structure.density_vol.dim[2] / dim[2]);
    gl::init_texture_3D(&rep->electronic_structure.color_vol.tex_id, dim[0], dim[1], dim[2], GL_RGBA8);

    const md_float3_t& voxel_size     = rep->electronic_structure.color_vol.voxel_size;
    const md_float4x4_t& world_to_model = rep->electronic_structure.color_vol.world_to_model;
    md_float4x4_t index_to_world = rep->electronic_structure.color_vol.texture_to_world
        * md_float4x4_scale(1.0f / dim[0], 1.0f / dim[1], 1.0f / dim[2])
        * md_float4x4_translate(0.5f, 0.5f, 0.5f); // Center of the corner voxel should be at the origin

    // The colours are per SYSTEM atom and the splats are per BASIS atom, so the map is consulted
    // here too - and once more it is the same lookup in the same direction the evaluation used.
    size_t map_count = 0;
    const uint32_t* map = qm_atom_map_find(&map_count, sys);
    if (map && map_count < num_points) {
        return;
    }

    md_float4_t*   point_xyzw   = (md_float4_t*)md_temp_alloc(temp, sizeof(md_float4_t) * num_points);
    uint32_t* point_colors = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_points);
    if (!point_xyzw || !point_colors) {
        return;
    }
    if (sys_state.num_atoms == 0 || !sys_state.xyz) {
        return;
    }
    for (size_t i = 0; i < num_points; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= sys_state.num_atoms) {
            return;
        }
        // Angstrom and not Bohr: these are world coordinates for the splatting pass, which shares
        // the Volume's transforms, unlike the grid the density was evaluated on.
        const float radius = md_atom_radius(&sys.atom, idx);
        point_xyzw[i]   = md_float4_set(sys_state.xyz[idx].x, sys_state.xyz[idx].y, sys_state.xyz[idx].z, radius);
        point_colors[i] = atom_colors[idx];
    }

    volume::compute_point_color_volume(rep->electronic_structure.color_vol.tex_id, dim, voxel_size.elem, world_to_model.elem,
                                       index_to_world.elem, point_xyzw, point_colors, num_points,
                                       rep->electronic_structure.gaussian_splatting_power);
}

// The values of a volume texture, read back. A volume the md_gpu path is still filling is waited
// for first (the caller waits for it without blocking: gpu_volume_job_in_flight_for), and image stores
// of the GL compute path are made visible before the read.
static bool volume_read_values(float* dst, ApplicationState* state, const Volume& vol) {
    if (!vol.tex_id) return false;
#if MD_ENABLE_GPU
    // its own readback only: draining every one would wait for whatever else the stream holds
    if (gpu_volume_job_in_flight_for(state, vol.tex_id)) gpu_volume_jobs_drain(state);
#else
    (void)state;
#endif
    if (glMemoryBarrier) {
        glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
    }
    GLint prev = 0;
    glGetIntegerv(GL_TEXTURE_BINDING_3D, &prev);
    glBindTexture(GL_TEXTURE_3D, vol.tex_id);
    glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, dst);
    glBindTexture(GL_TEXTURE_3D, (GLuint)prev);
    return true;
}

// The field the isosurfaces are coloured by, over the box of the density volume and around the
// surfaces drawn from it (surface_field.h). The VALUES depend on the field and the geometry - the
// frame - and on nothing else the representation chooses, so they are kept while the surfaces
// change: a new isovalue, orbital, density or resolution evaluates only what its band adds, and on
// the GPU, which evaluates the whole field once, nothing at all.
// False while the density it reads is still being evaluated (to be called again)
static bool electronic_structure_field_update(ApplicationState* state, Representation* rep, const ElectronicStructureRepresentation& g) {
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    const md_system_t& sys = state->mold.sys;
    if (!surface_field_available(g.field_kind, sys)) {
        surface_field_free(&es.field_vol);
        return true;
    }

    IsoDesc iso;
    electronic_structure_iso_desc_init(&iso, g);

    // The values are of the field at this frame; the surfaces are what the density volume holds
    // (vol_hash names it) cut at the isovalues
    uint64_t source_hash = md_hash64(&g.field_kind, sizeof(g.field_kind), 0);
    source_hash = md_hash64(&state->animation.time, sizeof(state->animation.time), source_hash);
    uint64_t band_hash = md_hash64(iso.values, sizeof(float) * iso.count, (uint64_t)iso.count + 1);
    band_hash = md_hash64_combine(band_hash, es.vol_hash);
    if (source_hash == es.field_vol.source_hash && band_hash == es.field_vol.band_hash) {
        return true;    // what is under way is taken by electronic_structure_field_poll
    }

    const size_t num_voxels = md_grid_num_points(&es.grid);
    if (num_voxels == 0) return true;
#if MD_ENABLE_GPU
    // The density still being evaluated on the GPU: read once it is there, a later frame, rather than waited for
    if (es.density_vol.tex_id && gpu_volume_job_in_flight_for(state, es.density_vol.tex_id)) return false;
#endif

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    float* density = (float*)md_temp_alloc(temp, sizeof(float) * num_voxels);
    if (!density || !volume_read_values(density, state, es.density_vol)) {
        return true;
    }

    SurfaceFieldDesc desc = {};
    desc.kind        = g.field_kind;
    desc.sys         = &sys;
    desc.state       = &system_instance(state->system_context).state;
    desc.grid        = &es.grid;
    desc.density     = density;
    desc.iso_values  = iso.values;
    desc.num_iso     = iso.count;
    desc.source_hash = source_hash;
    desc.band_hash   = band_hash;
    desc.density_cache = density_cache_of(state);
    if (surface_field_uses_basis(g.field_kind)) {
        // The QM nuclei and basis where the density was drawn: the current frame, as for the density
        desc.basis_atom_xyz = basis_atom_positions_extract(&desc.num_basis_atoms, temp, sys, system_instance(state->system_context).state);
    }
#if MD_ENABLE_GPU
    desc.gpu_stream = state->gpu_stream;
#endif

    const md_tick_t t0 = md_tick_now();
    if (surface_field_update(&es.field_vol, desc)) {
        color_scale_update_range(&es.field_map, surface_field_span(es.field_vol));
        const SurfaceFieldVolume& fv = es.field_vol;
        MD_LOG_DEBUG("Surface field: %zu of %zu voxels (%dx%dx%d) evaluated%s%s, the update %.1f ms", fv.num_evaluated, fv.num_voxels,
                     fv.grid.dim[0], fv.grid.dim[1], fv.grid.dim[2], fv.on_gpu ? " on the GPU" : "", surface_field_busy(fv) ? ", more under way" : "",
                     md_tick_to_milliseconds(md_tick_now() - t0));
    }
    return true;
}

void electronic_structure_field_poll(Representation* rep) {
    MD_ASSERT(rep);
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    // the source evaluated as a surface field: taken by the update once complete (es_source_volume_poll), which this
    // has look again when it is
    if (surface_field_busy(es.source_vol) && surface_field_poll(&es.source_vol) && es.source_vol.complete) {
        rep->needs_update = true;
    }
    if (!surface_field_busy(es.field_vol)) return;
    if (surface_field_poll(&es.field_vol)) {
        color_scale_update_range(&es.field_map, surface_field_span(es.field_vol));
        const SurfaceFieldVolume& fv = es.field_vol;
        MD_LOG_DEBUG("Surface field: %zu of %zu voxels (%dx%dx%d) evaluated%s; on %zu surface samples min %g, 1%% %g, 99%% %g, max %g",
                     fv.num_evaluated, fv.num_voxels, fv.grid.dim[0], fv.grid.dim[1], fv.grid.dim[2], fv.on_gpu ? " on the GPU" : "",
                     fv.num_surface_samples, fv.surface_min, fv.surface_lo, fv.surface_hi, fv.surface_max);
    }
}

// ---------------------------------------------------------------------------
// Where to evaluate
// ---------------------------------------------------------------------------
// An object aligned box around the BASIS atoms, padded, at the requested sample density. Derived
// from the system on every call rather than cached at load, for two reasons: it is a PCA over the
// tens to hundreds of atoms a QM calculation covers, which is nothing next to the evaluation it
// precedes, and computing it here is what makes the box follow the geometry instead of pinning it
// to whatever the coordinates were when the file was opened.
bool electronic_structure_grid_init(md_grid_t* grid, const md_system_t& sys, const md_system_state_t& state, double samples_per_unit_length) {
    MD_ASSERT(grid);

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return false;
    }

    const size_t num_atoms = md_gto_basis_num_atoms(&basis);
    if (num_atoms == 0) {
        return false;
    }

    md_float3_t* xyz = (md_float3_t*)md_temp_alloc(temp, sizeof(md_float3_t) * num_atoms);
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, state, num_atoms) != num_atoms) {
        return false;
    }

    // mat3_PCA and calculate_bounds both want a homogeneous point, and w == 1 is what makes the
    // rotation in calculate_bounds a rotation of a POINT rather than of a direction.
    md_float4_t* xyzw = (md_float4_t*)md_temp_alloc(temp, sizeof(md_float4_t) * num_atoms);
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = md_float4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
    }

    OABB oabb = {};
    oabb.orientation = mat3_PCA(xyzw, num_atoms);
    calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, xyzw, num_atoms, oabb.orientation);

    init_grid(grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
    return true;
}

// ---------------------------------------------------------------------------
// Evaluation: the entry points a representation calls
// ---------------------------------------------------------------------------
// One function per KIND of thing being evaluated - an orbital, a density - and the choice of
// backend inside it. A caller names the attribute it wants drawn and gets a filled texture; which
// device did the work is not its business, and pushing that choice down here is what lets the two
// backends share the basis, the atom gather and the map they all agree on.
//
// The GPU path writes the device scratch volume and queues a readback into vol_tex; the texture is
// filled from the frame loop's md_gpu_device_poll rather than before this returns. Call
// gpu_volume_jobs_drain() where the contents are needed immediately.

#if MD_ENABLE_GPU
// True when the device scratch and this dataset's uploaded basis are both present, which is what
// every GPU path below needs before it can start.
static bool gpu_evaluation_ready(const ApplicationState* state) {
    return state->gpu_stream && state->mold.gpu_basis && state->mold.gpu_atoms && state->gpu_coeff && state->gpu_volume;
}
#endif

// Orbitals and densities are evaluated on a grid by md_gpu. Without a device for it, or in a build without
// MD_ENABLE_GPU, there is nothing to evaluate them with: said once, rather than every frame something asks.
static void gpu_evaluation_unavailable(const char* what) {
    static bool said = false;
    if (!said) {
        MD_LOG_ERROR("Cannot evaluate %s: evaluation on a grid needs a GPU device (md_gpu), and none is available", what);
        said = true;
    }
}

bool orbital_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, md_str_t coefficient_path,
                      const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    MD_ASSERT(state);
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    size_t num_ao = 0;
    const double* ao_coeffs = orbital_coefficients_extract(&num_ao, temp, state->mold.sys, coefficient_path, slice);
    if (!ao_coeffs) {
        return false;
    }
    return orbital_coefficients_evaluate(state, vol_tex, grid, ao_coeffs, num_ao, mode, op, cutoff);
}

bool orbital_coefficients_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, const double* ao_coeffs, size_t num_ao,
                                   md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    MD_ASSERT(state);
    if (!ao_coeffs || num_ao == 0) return false;

#if MD_ENABLE_GPU
    if (gpu_evaluation_ready(state)) {
        const md_system_t&       sys       = state->mold.sys;
        const md_system_state_t& sys_state = system_instance(state->system_context).state;

        const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
        if (num_cgtos != num_ao) {
            MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and the orbital %zu", num_cgtos, num_ao);
            return false;
        }
        gpu_density_jobs_settle(state);     // the inputs are changed below
        if (!gpu_atoms_ensure_uploaded(state, sys, sys_state)) {
            return false;
        }

        const double* coeff_ptrs[1] = { ao_coeffs };
        float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_mo(1, num_cgtos));
        if (!dst) {
            return false;
        }
        md_gto_gpu_coeff_pack_mo(dst, coeff_ptrs, nullptr, 1, num_cgtos);
        md_gpu_upload_end(state->gpu_stream);

        md_gto_gpu_orbital_desc_t desc = {
            .basis        = state->mold.gpu_basis,
            .atom_xyz     = state->mold.gpu_atoms,
            .coeff        = state->gpu_coeff,
            .out_tex      = gpu_eval_target(state),
            .out_grad     = gpu_grad_target(state),
            .grid         = &grid,
            .sample_offset = {0.5f, 0.5f, 0.5f},
            .num_orbitals = 1,
            .eval_mode    = mode,
            .op           = op,
        };
        md_gto_gpu_orbital_launch(state->gpu_stream, &desc);
        return gpu_volume_readback(state, vol_tex, grid);
    }
#endif

    (void)vol_tex; (void)grid; (void)mode; (void)op; (void)cutoff;
    gpu_evaluation_unavailable("an orbital");
    return false;
}

#if MD_ENABLE_GPU
static void density_gpu_launch(ApplicationState* state, const md_grid_t& grid, md_gto_op_t op, md_gto_gpu_dmat_t dmat = nullptr, md_gto_gpu_job_t** out_job = nullptr,
                               md_gpu_addr_t coeff = 0);

bool density_matrix_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid,
                                           const double* density_matrix, size_t dim, md_gto_op_t op) {
    MD_ASSERT(state);
    if (!density_matrix || dim == 0 || !gpu_evaluation_ready(state)) {
        return false;
    }

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = system_instance(state->system_context).state;

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    if (num_cgtos != dim) {
        MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and the density matrix %zu", num_cgtos, dim);
        return false;
    }
    gpu_density_jobs_settle(state);     // the inputs are changed below
    if (!gpu_atoms_ensure_uploaded(state, sys, sys_state)) {
        return false;
    }

    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_density(num_cgtos));
    if (!dst) {
        return false;
    }
    md_gto_gpu_coeff_pack_density(dst, density_matrix, num_cgtos);
    md_gpu_upload_end(state->gpu_stream);

    density_gpu_launch(state, grid, op);
    return true;
}

// The density kernel over whatever was last uploaded to the coefficient buffer, or over a dmat, or over coeff (a
// packed triangle of its own). out_job: as a job (md_gto_gpu_density_job_begin), which the caller completes, instead
// of a launch that waits for its counts
static void density_gpu_launch(ApplicationState* state, const md_grid_t& grid, md_gto_op_t op, md_gto_gpu_dmat_t dmat, md_gto_gpu_job_t** out_job,
                               md_gpu_addr_t coeff) {
    md_gto_gpu_density_desc_t desc = {
        .basis         = state->mold.gpu_basis,
        .atom_xyz      = state->mold.gpu_atoms,
        .coeff         = dmat ? 0 : (coeff ? coeff : state->gpu_coeff),
        .dmat          = dmat,
        .out_tex       = gpu_eval_target(state),
        .out_grad      = gpu_grad_target(state),
        .grid          = &grid,
        .sample_offset = {0.5f, 0.5f, 0.5f},
        .op            = op,
    };
    if (out_job) {
        *out_job = md_gto_gpu_density_job_begin(state->gpu_stream, &desc);
    } else {
        md_gto_gpu_density_launch(state->gpu_stream, &desc);
    }
}

// How far the atoms may move before the dmat's pattern is made again (bohr): the pattern keeps the AO pairs
// up to this much further apart, about 13% more tiles for a protein, and holds while no atom has moved
// more than half of it.
static constexpr double GPU_DMAT_MARGIN = 1.0;

static bool density_is_scf(md_str_t density_path) {
    return md_str_eq(density_path, es_path::alpha_density) || md_str_eq(density_path, es_path::beta_density) ||
           md_str_eq(density_path, es_path::total_density) || md_str_eq(density_path, es_path::difference_density);
}

// What a dmat holds: the density, and the versions of the orbitals and occupations it is made of
static uint64_t density_dmat_key(const md_system_t& sys, md_str_t density_path) {
    uint64_t key = md_hash64(density_path.ptr, density_path.len, 0x5cf0d3a7);
    const md_str_t sources[4] = {es_path::alpha_coefficient, es_path::alpha_occupation, es_path::beta_coefficient, es_path::beta_occupation};
    for (const md_str_t& p : sources) {
        const md_attribute_t* a = md_attributes_find(&sys.attributes, p);
        const uint64_t v = a ? md_attributes_version(&sys.attributes, a->id) : 0;
        key = md_hash64(&v, sizeof(v), key);
    }
    return key ? key : 1;
}

// The density at density_path into dm, from its orbitals (md_qm_density_orbitals), on stream s
static bool density_dmat_fill(md_gpu_stream_t s, md_gto_gpu_dmat_t dm, const md_system_t& sys, md_str_t density_path, size_t num_cgtos) {
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    double* rows = nullptr;
    double* weights = nullptr;
    size_t num_ao = 0;
    const size_t count = md_qm_density_orbitals(&rows, &weights, &num_ao, &sys, density_path, md_temp_allocator(temp));
    if (count == 0 || num_ao != num_cgtos) {
        return false;
    }
    return md_gto_gpu_dmat_build(s, dm, rows, weights, count);
}

// A dmat whose pattern is made where the atoms are (xyz, bohr), on stream s
static md_gto_gpu_dmat_t density_dmat_pattern(md_gpu_stream_t s, const md_system_t& sys, const md_float3_t* xyz) {
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return nullptr;
    }
    // The cutoff of gpu_basis (system_gpu_data_update), whose radii the pattern is made with
    const md_gto_gpu_dmat_desc_t desc = {
        .basis = &basis, .atom_xyz = (const float*)xyz, .atom_xyz_stride = sizeof(md_float3_t),
        .cutoff = DEFAULT_GTO_CUTOFF_VALUE, .margin = GPU_DMAT_MARGIN,
    };
    return md_gto_gpu_dmat_create(s, &desc);
}

// A dmat made on a worker thread, on a GPU stream of its own: the pattern, the orbitals and the build are for a
// protein's total density a third of a second of the CPU, more than a second the first time the orbitals are
// read, which the representations' frames do not wait for (density_dmat_prepare). Made anew rather than
// refilled, because the one in use may be read by an evaluation still on the GPU.
struct DmatBuild {
    task_system::Job task;              // a cancelled one (a load) made nothing to be taken
    uint64_t key = 0;
    md_str_t path = {};                    // an es_path constant
    const md_system_t* sys = nullptr;
    md_float3_t* xyz = nullptr;              // bohr, its copy
    size_t num_cgtos = 0;
    md_gpu_stream_t stream = nullptr;   // the task's, destroyed by it once its work is done (its staging given back)
    md_gto_gpu_dmat_t dmat = nullptr;
    md_tick_t t0 = 0;
};

static void dmat_build_destroy(DmatBuild* b, md_gpu_stream_t gl_stream) {
    // only once its task has ended (dmat_build_take)
    if (b->dmat) md_gto_gpu_dmat_destroy(gl_stream, b->dmat);
    if (b->stream) md_gpu_stream_destroy(b->stream);    // of a task cancelled before it ran: nothing on it
    free(b->xyz);
    delete b;
}

// The build under way, taken once it is done: its dmat becomes the one in use. wait: waited for if it is not
// done. Returns false while it is under way.
static bool dmat_build_take(ApplicationState* state, bool wait) {
    DmatBuild* b = (DmatBuild*)state->mold.gpu_dmat_build;
    if (!b) return true;
    if (wait) b->task.wait();
    if (b->task.running()) return false;
    state->mold.gpu_dmat_build = nullptr;
    const bool made = b->task.status() == task_system::Status::Completed;
    if (made && b->dmat) {
        gpu_density_jobs_settle(state);     // one not yet issued reads the dmat in use
        // its build is done on the GPU: the task waited for it before it ended (destroying its stream)
        if (state->mold.gpu_dmat) md_gto_gpu_dmat_destroy(state->gpu_stream, state->mold.gpu_dmat);
        state->mold.gpu_dmat = b->dmat;
        state->mold.gpu_dmat_key = b->key;
        b->dmat = nullptr;
        MD_LOG_DEBUG("The density matrix of '" MD_STR_FMT "' made in %.1f ms (on a worker thread)", MD_STR_ARG(b->path), md_tick_to_milliseconds(md_tick_now() - b->t0));
    } else if (made) {
        state->mold.gpu_dmat_failed_key = b->key;     // made where it is asked for (density_dmat_ensure), which may fail as loudly
    }
    dmat_build_destroy(b, state->gpu_stream);
    return true;
}

// Any other density a representation evaluates, as the packed float triangle the kernel reads: a transition density
// is made by its provider on every read (seconds for a protein: md_attributes.h, VIRTUAL ATTRIBUTES), a density
// property read from its file, either on the GL thread a stall of the frame as long. Extracted on a worker thread
// straight into a device buffer of its own (the upload's staging), on a GPU stream of its own, and kept for the
// evaluations that follow (state->mold.gpu_dpack): a new resolution or frame of the same density reads it again
// from the device.
struct DensityPackBuild {
    task_system::Job task;              // a cancelled one (a load) made nothing to be taken
    uint64_t key = 0;
    md_str_t path = {};                    // an es_path constant, or a path owned by the system's attribute table
    md_attribute_slice_t slice = {};
    const md_system_t* sys = nullptr;
    size_t num_cgtos = 0;
    md_gpu_stream_t stream = nullptr;   // the task's, destroyed by it once its uploads are done
    md_gpu_addr_t buf = 0;
    bool ok = false;
    md_tick_t t0 = 0;
};

// What a packed density holds: the attribute, the slice of it and its version
static uint64_t density_pack_key(const md_system_t& sys, const md_attribute_t& attr, const md_attribute_slice_t& slice) {
    uint64_t key = md_hash64(attr.path.ptr, attr.path.len, 0x3a9d07c1);
    key = md_hash64(&slice.num_idx, sizeof(slice.num_idx), key);
    if (slice.num_idx > 0) key = md_hash64(slice.idx, sizeof(uint32_t) * slice.num_idx, key);
    const uint64_t v = md_attributes_version(&sys.attributes, attr.id);
    key = md_hash64(&v, sizeof(v), key);
    return key ? key : 1;
}

// The build under way, taken once it is done: its buffer becomes the one in use. wait: waited for if it is not
// done. Returns false while it is under way.
static bool dpack_build_take(ApplicationState* state, bool wait) {
    DensityPackBuild* b = (DensityPackBuild*)state->mold.gpu_dpack_build;
    if (!b) return true;
    if (wait) b->task.wait();
    if (b->task.running()) return false;
    state->mold.gpu_dpack_build = nullptr;
    const bool made = b->task.status() == task_system::Status::Completed;
    if (made && b->ok) {
        gpu_density_jobs_settle(state);     // one not yet all issued reads the buffer in use
        if (state->mold.gpu_dpack) md_gpu_free(state->gpu_stream, state->mold.gpu_dpack);
        state->mold.gpu_dpack = b->buf;
        state->mold.gpu_dpack_key = b->key;
        b->buf = 0;
        MD_LOG_DEBUG("The density '" MD_STR_FMT "' extracted to the device in %.1f ms (on a worker thread)", MD_STR_ARG(b->path), md_tick_to_milliseconds(md_tick_now() - b->t0));
    } else if (made) {
        state->mold.gpu_dpack_failed_key = b->key;   // extracted where it is asked for, which fails as loudly
    }
    if (b->buf) md_gpu_free(state->gpu_stream, b->buf);
    if (b->stream) md_gpu_stream_destroy(b->stream);    // of a task cancelled before it ran: nothing on it
    delete b;
    return true;
}


// For the representations' evaluations, which wait for it without waiting on it: true when the dmat holds the
// SCF density at density_path for the atoms where they are (the evaluation can be issued), or when the density
// does not go through one; false while it is being made on a worker thread, which this begins when it is not.
static bool density_dmat_prepare(ApplicationState* state, md_str_t density_path) {
    if (!density_is_scf(density_path) || !gpu_evaluation_ready(state) || !state->gpu_device) {
        return true;
    }
    if (!dmat_build_take(state, false)) return false;
    const md_system_t& sys = state->mold.sys;
    const uint64_t key = density_dmat_key(sys, density_path);
    if (key == state->mold.gpu_dmat_failed_key) return true;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    size_t num_atoms = 0;
    const md_float3_t* xyz = basis_atom_positions_extract(&num_atoms, temp, sys, system_instance(state->system_context).state);
    if (!xyz || num_atoms != md_gto_gpu_basis_num_atoms(state->mold.gpu_basis)) {
        return true;    // the evaluation fails as it would have
    }
    if (state->mold.gpu_dmat && state->mold.gpu_dmat_key == key && md_gto_gpu_dmat_pattern_valid(state->mold.gpu_dmat, (const float*)xyz, sizeof(md_float3_t))) {
        return true;
    }

    DmatBuild* b = new DmatBuild();
    b->key = key;
    b->path = density_path;
    b->sys = &state->mold.sys;
    b->num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    b->t0 = md_tick_now();
    b->xyz = (md_float3_t*)malloc(sizeof(md_float3_t) * num_atoms);
    b->stream = b->xyz ? md_gpu_stream_create(state->gpu_device, MD_GPU_STREAM_COMPUTE, "density matrix") : nullptr;
    if (!b->stream) {
        free(b->xyz);
        delete b;
        return true;    // made where it is asked for
    }
    md_memcpy(b->xyz, xyz, sizeof(md_float3_t) * num_atoms);
    b->task.start(MD_STR_LIT("Density matrix: preparing"), [b]() {
        b->dmat = density_dmat_pattern(b->stream, *b->sys, b->xyz);
        if (b->dmat && !density_dmat_fill(b->stream, b->dmat, *b->sys, b->path, b->num_cgtos)) {
            md_gto_gpu_dmat_destroy(b->stream, b->dmat);
            b->dmat = nullptr;
        }
        md_gpu_stream_destroy(b->stream);   // once its work is done: the dmat is then ready for any stream
        b->stream = nullptr;
        md_temp_trim();     // the orbitals were read into this worker's temp arenas
    });
    state->mold.gpu_dmat_build = b;
    return false;
}

// For the representations' evaluations, as density_dmat_prepare and for any density: true when the density at
// density_path (sliced) can be issued now - an SCF density's dmat made, any other's packed triangle on the device
// (state->mold.gpu_dpack) - or when it is not to be had this way (evaluated as it would have been); false while it
// is being made on a worker thread, which this begins when it is not.
static bool density_pack_prepare(ApplicationState* state, md_str_t density_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t s = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    if (s.num_idx == 0 && density_is_scf(density_path)) {
        if (!density_dmat_prepare(state, density_path)) return false;
        // through the dmat, unless it could not be made of the orbitals (not published): through the packed triangle
        if (state->mold.gpu_dmat_failed_key != density_dmat_key(state->mold.sys, density_path)) return true;
    }
    if (!gpu_evaluation_ready(state) || !state->gpu_device) {
        return true;
    }
    if (!dpack_build_take(state, false)) return false;
    const md_system_t& sys = state->mold.sys;
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) return true;
    const uint64_t key = density_pack_key(sys, *attr, s);
    if ((state->mold.gpu_dpack && key == state->mold.gpu_dpack_key) || key == state->mold.gpu_dpack_failed_key) return true;
    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    if (md_qm_extract_packed_symmetric_f32(nullptr, 0, &sys.attributes, attr, s, nullptr) != num_cgtos) {
        return true;    // the evaluation fails as it would have, and says why
    }

    DensityPackBuild* b = new DensityPackBuild();
    b->key = key;
    b->path = attr->path;
    b->slice = s;
    b->sys = &state->mold.sys;
    b->num_cgtos = num_cgtos;
    b->t0 = md_tick_now();
    b->stream = md_gpu_stream_create(state->gpu_device, MD_GPU_STREAM_COMPUTE, "packed density");
    if (!b->stream) {
        delete b;
        return true;    // extracted where it is asked for
    }
    b->task.start(MD_STR_LIT("Density matrix: reading"), [b]() {
        const md_attribute_t* a = md_attributes_find(&b->sys->attributes, b->path);
        const size_t len = b->num_cgtos * (b->num_cgtos + 1) / 2;
        if (a) {
            b->buf = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, sizeof(float) * len).gpu;
            float* dst = b->buf ? (float*)md_gpu_upload_begin(b->stream, b->buf, sizeof(float) * len) : nullptr;
            if (dst) {
                DensityIo io;
                const bool extracted = md_qm_extract_packed_symmetric_f32(dst, len, &b->sys->attributes, a, b->slice, io.io) == b->num_cgtos;
                b->ok = md_gpu_upload_end(b->stream) && extracted;
            }
        }
        md_gpu_stream_destroy(b->stream);   // once its uploads are done: the buffer is then ready for any stream
        b->stream = nullptr;
        md_temp_trim();     // what the provider made in this worker's temp arenas
    });
    state->mold.gpu_dpack_build = b;
    return false;
}

// The SCF density at density_path (orbital/{alpha,beta,total,difference}/density) in state->mold.gpu_dmat,
// built on the device from its orbitals (md_qm_density_orbitals): on this thread unless a worker has made it
// (density_dmat_prepare), and one under way is waited for. False for any other density, or one whose orbitals
// are not published or that is zero: those go through the packed matrix.
static bool density_dmat_ensure(ApplicationState* state, md_str_t density_path) {
    if (!density_is_scf(density_path) || !state->gpu_stream || !state->mold.gpu_basis) {
        return false;
    }
    dmat_build_take(state, true);
    const md_system_t& sys = state->mold.sys;
    const uint64_t key = density_dmat_key(sys, density_path);

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    size_t num_atoms = 0;
    const md_float3_t* xyz = basis_atom_positions_extract(&num_atoms, temp, sys, system_instance(state->system_context).state);
    if (!xyz || num_atoms != md_gto_gpu_basis_num_atoms(state->mold.gpu_basis)) {
        return false;
    }
    // A pattern made where the atoms were: kept while they stay within the margin, else made again
    if (state->mold.gpu_dmat && !md_gto_gpu_dmat_pattern_valid(state->mold.gpu_dmat, (const float*)xyz, sizeof(md_float3_t))) {
        md_gto_gpu_dmat_destroy(state->gpu_stream, state->mold.gpu_dmat);
        state->mold.gpu_dmat = nullptr;
    }
    if (!state->mold.gpu_dmat) {
        state->mold.gpu_dmat = density_dmat_pattern(state->gpu_stream, sys, xyz);
        state->mold.gpu_dmat_key = 0;
        if (!state->mold.gpu_dmat) {
            return false;
        }
    }
    if (key != state->mold.gpu_dmat_key) {
        if (!density_dmat_fill(state->gpu_stream, state->mold.gpu_dmat, sys, density_path, md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis))) {
            return false;
        }
        state->mold.gpu_dmat_key = key;
    }
    return true;
}

static bool density_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, md_str_t density_path,
                                  const md_attribute_slice_t* slice, md_gto_op_t op, md_gto_gpu_job_t** out_job);

bool density_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, md_str_t density_path,
                                    const md_attribute_slice_t* slice, md_gto_op_t op) {
    return density_to_gpu_volume(state, grid, density_path, slice, op, nullptr);
}

// out_job: the evaluation as a job (density_gpu_launch), when it is one; nullptr: launched
static bool density_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, md_str_t density_path,
                                  const md_attribute_slice_t* slice, md_gto_op_t op, md_gto_gpu_job_t** out_job) {
    MD_ASSERT(state);
    if (out_job) *out_job = nullptr;
    if (!gpu_evaluation_ready(state)) {
        return false;
    }
    gpu_density_jobs_settle(state);     // the inputs are changed below

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = system_instance(state->system_context).state;
    const md_attribute_slice_t s = slice ? *slice : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" MD_STR_FMT "'", MD_STR_ARG(density_path));
        return false;
    }

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    const size_t dim = md_qm_extract_packed_symmetric_f32(nullptr, 0, &sys.attributes, attr, s, nullptr);
    if (dim == 0) {
        return false;
    }
    if (num_cgtos != dim) {
        MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and the density matrix %zu", num_cgtos, dim);
        return false;
    }
    if (!gpu_atoms_ensure_uploaded(state, sys, sys_state)) {
        return false;
    }

    // An SCF density through the dmat, built on the device from its orbitals: no matrix on the host, and
    // only the AO pairs the evaluation reads
    if (s.num_idx == 0 && density_dmat_ensure(state, density_path)) {
        density_gpu_launch(state, grid, op, state->mold.gpu_dmat, out_job);
        return true;
    }

    // Extracted on a worker thread for a representation (density_pack_prepare): on the device already
    dpack_build_take(state, false);
    if (state->mold.gpu_dpack && state->mold.gpu_dpack_key == density_pack_key(sys, *attr, s)) {
        density_gpu_launch(state, grid, op, nullptr, out_job, state->mold.gpu_dpack);
        return true;
    }

    // Straight into the upload buffer: the packed float triangle is exactly what the kernel reads, so
    // a packed attribute - a density property read from its file - never exists in any other form
    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_density(num_cgtos));
    if (!dst) {
        return false;
    }
    DensityIo io;
    const bool extracted = md_qm_extract_packed_symmetric_f32(dst, dim * (dim + 1) / 2, &sys.attributes, attr, s, io.io) == dim;
    md_gpu_upload_end(state->gpu_stream);
    if (!extracted) {
        return false;
    }

    density_gpu_launch(state, grid, op, nullptr, out_job);
    return true;
}
#endif

// The backend choice, over a density matrix the caller already holds. density_evaluate below is
// this plus the extract, and a caller which has cached the matrix - the QM residue does, because
// the transition densities are rebuilt on every read - calls this one and skips the rebuild.
bool gpu_volume_evaluation_available(const ApplicationState* state) {
#if MD_ENABLE_GPU
    return state && gpu_evaluation_ready(state);
#else
    (void)state;
    return false;
#endif
}

bool density_matrix_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid,
                             const double* density_matrix, size_t dim, md_gto_op_t op) {
    MD_ASSERT(state);
    if (!density_matrix || dim == 0) {
        return false;
    }

#if MD_ENABLE_GPU
    if (gpu_evaluation_ready(state)) {
        return density_matrix_evaluate_to_gpu_volume(state, grid, density_matrix, dim, op) && gpu_volume_readback(state, vol_tex, grid);
    }
#endif
    (void)vol_tex; (void)grid; (void)op;
    gpu_evaluation_unavailable("a density");
    return false;
}

bool density_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, md_str_t density_path,
                      const md_attribute_slice_t* slice, md_gto_op_t op, bool async) {
    MD_ASSERT(state);

#if MD_ENABLE_GPU
    // Evaluate, then read back. The split exists because a consumer that runs another kernel over
    // the result - the critical point extraction does - wants the first half and not the second.
    md_gto_gpu_job_t* job = nullptr;
    if (gpu_evaluation_ready(state)) {
        if (!density_to_gpu_volume(state, grid, density_path, slice, op, async ? &job : nullptr)) return false;
        return job ? gpu_volume_readback_after(state, vol_tex, grid, job, "Density volume") : gpu_volume_readback(state, vol_tex, grid);
    }
#endif
    (void)vol_tex; (void)grid; (void)density_path; (void)slice; (void)op; (void)async;
    gpu_evaluation_unavailable("a density");
    return false;
}

// The transfer function of a volume rendering, set up for the field when it is first drawn (or when the kind of its
// values changes), on a range set by hand: one to the largest values would put most of the field at the bottom of
// the curve. A signed field (an orbital, a difference) from -3 to 3 times the isovalue, the value it is significant
// at, through a diverging map, the curve on the magnitude from a tenth of the range; a density from 0 to 10 times it
// (its largest values are at the nuclei) through a sequential one; ELF and LOL over [0, 1], what they span, the curve
// from 1/2, the electron gas's value, so what shows is where the electrons are more localized.
// The colour map's range from the spread of the source's values (source_span), for what takes its range from them:
// the 1st to 99th percentile, about zero for a signed field
static void es_volume_map_range(ElectronicStructureRepresentation& es) {
    const ColorScaleSpan& sp = es.source_span;
    if (!sp.valid || !(sp.hi > sp.lo || electronic_structure_is_signed(es))) return;
    TransferFunction& tf = es.volume.tf;
    if (electronic_structure_is_signed(es)) {
        const float r = MD_MAX(fabsf(sp.lo), fabsf(sp.hi));
        if (!(r > 0.0f)) return;
        tf.range_beg = -r;
        tf.range_end =  r;
    } else {
        tf.range_beg = sp.lo;
        tf.range_end = sp.hi;
    }
    es.volume.map_auto = false;
    es.volume.view.fit = 1;
}

static void es_volume_map_update(ElectronicStructureRepresentation& es) {
    const ElectronicStructureValueKind kind = electronic_structure_value_kind(es);
    // a field of a range of its own is set up again for another one of the same kind (one potential for another)
    const int src = electronic_structure_source_is_surface_field(es.shown.source) ? (int)es.shown.source : -2;
    if (es.volume.map_init && es.volume.map_kind == (int)kind) {
        if (es.volume.map_source == -1) es.volume.map_source = src;
        if (es.volume.map_source == src) return;
    }
    TransferFunction& tf = es.volume.tf;
    const float iso = (float)fabs(es.iso_value);
    switch (kind) {
    case ElectronicStructureValueKind::Signed:
        tf.symmetric = true;
        tf.colormap = COLOR_SCALE_COLORMAP_RDBU;
        tf.range_beg = -3.0f * iso;
        tf.range_end =  3.0f * iso;
        transfer_function_ramp_magnitude(&tf, 0.1f);
        break;
    case ElectronicStructureValueKind::Localization:
        tf.symmetric = false;
        tf.colormap = COLOR_SCALE_COLORMAP_VIRIDIS;
        tf.range_beg = 0.0f;
        tf.range_end = 1.0f;
        transfer_function_ramp(&tf, 0.5f);
        break;
    case ElectronicStructureValueKind::Ranged:
        // a range of its own, from the values once they are there (es_source_volume_poll)
        tf.symmetric = false;
        tf.colormap = COLOR_SCALE_COLORMAP_VIRIDIS;
        tf.range_beg = 0.0f;
        tf.range_end = 1.0f;
        transfer_function_ramp(&tf, 0.0f);
        es.volume.map_auto = true;
        break;
    case ElectronicStructureValueKind::SignedRanged:
        tf.symmetric = true;
        tf.colormap = COLOR_SCALE_COLORMAP_RDBU;
        tf.range_beg = -1.0f;
        tf.range_end =  1.0f;
        transfer_function_ramp_magnitude(&tf, 0.1f);
        es.volume.map_auto = true;
        break;
    case ElectronicStructureValueKind::Density:
    default:
        tf.symmetric = false;
        tf.colormap = COLOR_SCALE_COLORMAP_VIRIDIS;
        tf.range_beg = 0.0f;
        tf.range_end = 10.0f * iso;
        transfer_function_ramp(&tf, 0.0f);
        break;
    }
    // a field of a colour map of its own (surface_field.h) on it, rather than the default of its value kind
    const SurfaceFieldKind fk = electronic_structure_source_field_kind(es.shown.source);
    if (fk != SurfaceFieldKind::Count) tf.colormap = surface_field_kind_info(fk).colormap;
    es.volume.view.fit = 1;
    es.volume.map_init = true;
    es.volume.map_kind = (int)kind;
    es.volume.map_source = src;
    // the range from the values that are there, when they are (else once they come)
    if (es.volume.map_auto) es_volume_map_range(es);
}

// What the representation's source is, of the frame: what an evaluation of it, on any grid, depends on
static uint64_t es_source_hash(const ApplicationState* state, const ElectronicStructureRepresentation& es) {
    uint64_t h = md_hash64(&state->animation.time, sizeof(state->animation.time), 0);
    h = md_hash64(&es.shown.source,                       sizeof(es.shown.source), h);
    h = md_hash64(&es.shown.use_magnitude,                sizeof(es.shown.use_magnitude), h);
    h = md_hash64(&es.shown.spin,                         sizeof(es.shown.spin), h);
    h = md_hash64(&es.shown.nto_component,                sizeof(es.shown.nto_component), h);
    h = md_hash64(&es.shown.transition_density_component, sizeof(es.shown.transition_density_component), h);
    h = md_hash64(&es.shown.orbital_idx,                  sizeof(es.shown.orbital_idx), h);
    h = md_hash64(&es.shown.excited_state_idx,            sizeof(es.shown.excited_state_idx), h);
    h = md_hash64(&es.shown.nto_lambda_idx,               sizeof(es.shown.nto_lambda_idx), h);
    h = md_hash64(&es.shown.density_property_key,         sizeof(es.shown.density_property_key), h);
    return h;
}

// Whether what an evaluation of the source reads is there: 1 yes, 0 not yet (being made: ask again), -1 not on the
// GPU (ELF and LOL on the CPU, which fills the volume itself)
static int es_source_ready(ApplicationState* state, const ElectronicStructureRepresentation& es) {
#if MD_ENABLE_GPU
    // a density job under way, of any representation: an evaluation now would wait for its counts
    if (gpu_density_job_pending(state)) return 0;
    // a density's matrix, made on a worker thread (an SCF density's dmat, another's packed triangle)
    md_str_t path = {};
    md_attribute_slice_t slice = {};
    if (es_volume_density(&path, &slice, es, state->mold.sys) && !density_pack_prepare(state, path, &slice)) return 0;
#endif
    if (es.shown.source == ElectronicStructureSource::NaturalTransitionOrbital && !es_nto_ready(state, (size_t)es.shown.excited_state_idx)) return 0;
    if (electronic_structure_source_is_surface_field(es.shown.source)) return -1;   // as a surface field, into the volume
    if (electronic_structure_is_localization(es)) {
#if MD_ENABLE_GPU
        const int g = localization_gpu_prepare(state);
        return g == 0 ? 0 : (g < 0 ? -1 : 1);
#else
        return -1;
#endif
    }
    return 1;
}

bool es_plane_images_supported(ApplicationState* state, const ElectronicStructureRepresentation& es) {
#if MD_ENABLE_GPU
    if (es.display != ElectronicStructureDisplay::Planes) return false;
    if (!gpu_evaluation_ready(state) || !state->gpu_plane) return false;
    // the fields evaluated as surface fields, and ELF and LOL without their kernel: into the volume, whole
    if (electronic_structure_source_is_surface_field(es.shown.source)) return false;
    return !electronic_structure_is_localization(es) || localization_gpu_prepare(state) >= 0;
#else
    (void)state; (void)es;
    return false;
#endif
}

// The box of the volume without the volume: its grid, at the resolution, and its unit cube in the world
static void es_box_update(ApplicationState* state, Representation* rep) {
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)electronic_structure_volume_resolution(es)];
    md_grid_t grid = {};
    if (!electronic_structure_grid_init(&grid, state->mold.sys, system_instance(state->system_context).state, samples_per_angstrom * BOHR_TO_ANGSTROM)) {
        es.has_box = false;
        return;
    }
    const float scl = (float)BOHR_TO_ANGSTROM;
    es.box_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, md_grid_extent(&grid) * scl);
    es.has_box = true;
}

bool electronic_structure_plane_images_update(ApplicationState* state, Representation* rep) {
#if MD_ENABLE_GPU
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    es_slices_conform(&es);
    if (!es.has_box || !es_plane_images_supported(state, es)) return true;
    const uint64_t src = es_source_hash(state, es);
    const float detail_spacing = plane_detail_spacing_angstrom[MD_CLAMP((int)es.plane_detail, 0, (int)PlaneDetail::Count - 1)];
    const double now = ImGui::GetTime();
    bool done = true;
    int ready = -2;     // asked when first needed, and again after each evaluation issued
    for (int i = 0; i < es.slices.count && i < ES_MAX_PLANES; ++i) {
        ElectronicStructurePlaneImage& im = es.plane_images[i];
        // the one being evaluated, there: drawn from now on
        if (im.pending) {
            if (gpu_volume_job_in_flight_for(state, im.img[1 - im.front].tex_id)) {
                done = false;
                continue;
            }
            im.front = 1 - im.front;
            im.pending = false;
            im.serial += 1;
            im.ms = (float)((now - im.issued) * 1000.0);
        }
        im.current = false;
        const ElectronicStructurePlane& p = es.slices.items[i];
        if (!p.enabled) continue;
        EsSliceFrame frame;
        if (!es_slice_frame(&frame, p)) continue;
        // what an image is of: the source and the rectangle (not the box: the volume's values past it are no part of
        // an image, evaluated wherever it is)
        uint64_t key = md_hash64(&frame, sizeof(frame), src);
        key = key ? key : 1;
        const int f = im.front;
        // all of the rectangle at the detail's spacing, at most GPU_PLANE_DIM a side
        const float ext = 2.0f * MD_MAX(frame.hu, frame.hv);
        const float spacing = MD_MAX(detail_spacing, ext / (float)GPU_PLANE_DIM);
        im.current = key == im.key[f];
        // while dragged: drawn from the volume, evaluated once let go
        if (es_planes_dragging(rep->node.id, i)) continue;
        if (im.current && fabsf(im.spacing[f] - spacing) <= 1.0e-3f * spacing) continue;
        if (ready == -2) ready = es_source_ready(state, es);
        if (ready <= 0) {
            done &= ready < 0;
            continue;
        }
        md_grid_t grid = {};
        es_slice_image_grid(&grid, frame, spacing);
        const int b = 1 - f;
        Volume& v = im.img[b];
        init_volume(&v, grid, GL_R32F);
        state->gpu_eval_target = state->gpu_plane;
        const bool ok = es_evaluate_into(state, es, v.tex_id, grid);
        state->gpu_eval_target = nullptr;
        if (!ok) continue;
        im.key[b] = key;
        im.spacing[b] = spacing;
        im.pending = true;
        im.issued = now;
        done = false;
        ready = -2;
    }
    return done;
#else
    (void)state; (void)rep;
    return true;
#endif
}

// The representation's source evaluated as a surface field (surface_field.h), all of the grid, for the volume hash
// 'key': 1 when it is in the volume (put there now), 0 while it is under way (started here if it was not), -1 when
// it cannot be made. g: what the volume is of.
static int es_source_volume_poll(ApplicationState* state, Representation* rep, const ElectronicStructureRepresentation& g) {
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    const md_system_t& sys = state->mold.sys;
    const SurfaceFieldKind kind = electronic_structure_source_field_kind(g.shown.source);
    if (kind == SurfaceFieldKind::Count || !surface_field_available(kind, sys)) return -1;

    md_grid_t grid = {};
    const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)electronic_structure_volume_resolution(g)];
    if (!electronic_structure_grid_init(&grid, sys, system_instance(state->system_context).state, samples_per_angstrom * BOHR_TO_ANGSTROM)) return -1;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    SurfaceFieldDesc desc = {};
    desc.kind        = kind;
    desc.sys         = &sys;
    desc.state       = &system_instance(state->system_context).state;
    desc.grid        = &grid;               // no density: all of it
    // what the values are of (the field, the frame); the grid is the volume's own concern: a new one keeps what the
    // field is made of (an electrostatic potential's distribution) and evaluates again
    desc.source_hash = es_source_hash(state, g);
    desc.band_hash   = 1;
    if (surface_field_uses_basis(kind)) {
        desc.basis_atom_xyz = basis_atom_positions_extract(&desc.num_basis_atoms, temp, sys, system_instance(state->system_context).state);
    }
#if MD_ENABLE_GPU
    desc.gpu_stream = state->gpu_stream;
#endif
    SurfaceFieldVolume& fv = es.source_vol;
    if (!surface_field_update(&fv, desc)) return -1;
    if (!fv.complete) return surface_field_busy(fv) ? 0 : -1;

    init_volume(&es.density_vol, fv.grid, GL_R32F);
    gl::set_texture_3D_data(es.density_vol.tex_id, 0, fv.values, GL_R32F);
    volume::notify_data_changed(es.density_vol.tex_id);
    es.grid = fv.grid;
    es.box_to_world = es.density_vol.texture_to_world;
    es.has_box = true;
    es.source_span = surface_field_span(fv);
    if (es.volume.map_auto) es_volume_map_range(es);
    MD_LOG_INFO("%s: %zu points (%dx%dx%d) on the %s; 1%% %g, 99%% %g", surface_field_kind_info(kind).label, fv.num_voxels,
                fv.grid.dim[0], fv.grid.dim[1], fv.grid.dim[2], fv.on_gpu ? "GPU" : "CPU", es.source_span.lo, es.source_span.hi);
    return 1;
}

bool es_mesh_mode(const ApplicationState* state, const ElectronicStructureRepresentation& es) {
    (void)state;
    return es.rendering == IsosurfaceRendering_Mesh &&
           electronic_structure_source_meshable(es.shown.source);
}

// Re-evaluating is expensive and everything it depends on is right here, so it is gated on a hash of
// exactly those inputs - the frame included, since the geometry moves the grid.
bool electronic_structure_representation_update(ApplicationState* state, Representation* rep, const uint32_t* atom_colors) {
    MD_ASSERT(state);
    MD_ASSERT(rep);
    ElectronicStructureRepresentation& es = rep->electronic_structure;

    uint64_t vol_hash = es_source_hash(state, es);
    const VolumeResolution resolution = electronic_structure_volume_resolution(es);
    vol_hash = md_hash64(&resolution, sizeof(resolution), vol_hash);

    // The mesh rendering mode draws the surfaces from the field itself (iso_meshes.h), and evaluates a surface
    // field at their vertices: the volume is only evaluated for what still reads it, the colours by atom. Of the
    // orbitals and densities only: the other fields' surfaces are raycast through their volume.
    const bool mesh_mode = es_mesh_mode(state, es);
    const bool iso_display = es.display == ElectronicStructureDisplay::Isosurface;
    // a volume rendering and a plane draw from (or over the box of) the volume, whatever renders the isosurfaces
    // (the planes display evaluated on the planes too: the volume is drawn where an image is not there, or not of the
    // plane as it is, as while it is dragged)
    const bool needs_volume = !mesh_mode || !iso_display || es.coloring == SurfaceColoring::AtomColors;
    vol_hash = md_hash64(&needs_volume, sizeof(needs_volume), vol_hash);

    bool done = true;
    if (vol_hash != es.vol_hash) {
        bool ready = true;
#if MD_ENABLE_GPU
        // The volume's evaluation under way, of what is no longer wanted: dropped where it is not yet all issued
        // (its pieces in flight complete), rather than waited for - for a large system seconds. Not while
        // playing, where each frame would drop the one before and none would be shown. The volume then holds
        // what it held before that one: of no hash (asked for again, it is evaluated again).
        if (es.density_vol.tex_id && state->animation.mode != PlaybackMode::Playing && gpu_volume_job_cancel(state, es.density_vol.tex_id)) {
            es.vol_hash = 0;
        }
        // What the evaluation needs, which the frame does not wait for: the volume's last evaluation off the GPU
        // (its readback in flight: a new one now would wait for it, and the volume's texture be remade under
        // it), and a density's matrix, made on a worker thread (an SCF density's dmat, another's packed
        // triangle). Until then the volume shown stays, and this is asked again next frame.
        if (needs_volume) {
            // (a density job under way of any representation, not yet all issued: an evaluation now would issue
            // all it has left at once)
            md_str_t path = {};
            md_attribute_slice_t slice = {};
            if ((es.density_vol.tex_id && gpu_volume_job_in_flight_for(state, es.density_vol.tex_id)) || gpu_density_job_pending(state)) {
                ready = false;
            } else if (es_volume_density(&path, &slice, es, state->mold.sys)) {
                ready = density_pack_prepare(state, path, &slice);
            }
        }
#endif
        // a natural transition orbital whose state is being read on a worker thread
        if (ready && needs_volume && es.shown.source == ElectronicStructureSource::NaturalTransitionOrbital) {
            ready = es_nto_ready(state, (size_t)es.shown.excited_state_idx);
        }
        // The fields of the system as a whole, and ELF and LOL without their kernel: as a surface field, all of the
        // grid, put in the volume once it is there. ELF and LOL on the GPU once what their kernel reads is made (by a
        // pool task). Until then the volume shown stays (a failure is taken as done, not tried again every frame).
        bool as_field = needs_volume && electronic_structure_source_is_surface_field(es.shown.source);
        if (ready && needs_volume && electronic_structure_is_localization(es)) {
#if MD_ENABLE_GPU
            const int k = localization_gpu_prepare(state);
            if (k == 0) ready = false;
            as_field = k < 0;
#else
            as_field = true;
#endif
        }
        if (ready && as_field) {
            ready = es_source_volume_poll(state, rep, es) != 0;
        }
        if (ready) {
            es.vol_hash = vol_hash;
            if (needs_volume && !as_field) {
                electronic_structure_evaluate(state, rep, es);
            } else if (!needs_volume) {
                es_box_update(state, rep);
            }
        } else {
            done = false;
        }
    }

    if (es.display != ElectronicStructureDisplay::Isosurface) {
        es_volume_map_update(es);
    }

    if (iso_display && es.coloring == SurfaceColoring::AtomColors) {
        uint64_t col_hash = md_hash64(&es.gaussian_splatting_power, sizeof(es.gaussian_splatting_power), (int)rep->color_mapping);
        col_hash = md_hash64_combine(col_hash, es.vol_hash);
        if (col_hash != es.col_hash) {
            es.col_hash = col_hash;
            electronic_structure_color_volume_update(state, rep, atom_colors);
        }
    } else if (iso_display && es.coloring == SurfaceColoring::Field && !mesh_mode) {
        done &= electronic_structure_field_update(state, rep, es);
    }
    return done;
}

// Every leaf under the density property group is one property, in the table's own path order. The
// label is the attribute's own, falling back to the last segment of its path when the file gave it
// none - both borrowed from the table, so nothing is copied and nothing has to be freed.
//
// Counts past 'cap' on purpose: the return value is how many the system HAS, so a caller with a
// fixed buffer can tell that it saw all of them.
size_t density_properties_gather(DensityProperty out[], size_t cap, const md_system_t& sys) {
    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, es_path::density_property); md_attributes_next(&it);) {
        const md_attribute_t* attr = it.attr;
        if (out && count < cap) {
            const md_str_t label = md_str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
            out[count] = { .key = attr->id, .label = label };
        }
        count += 1;
    }
    return count;
}

ElectronicStructureSourceFlags es_source_mask(const md_system_t& sys) {
    ElectronicStructureSourceFlags mask = 0;

    if (attribute_exists(sys, es_path::alpha_coefficient))  mask |= ElectronicStructureSourceFlag_MolecularOrbital;
    if (attribute_exists(sys, es_path::alpha_density))      mask |= ElectronicStructureSourceFlag_ElectronDensity;
    if (attribute_exists(sys, es_path::nto_particle))       mask |= ElectronicStructureSourceFlag_NaturalTransitionOrbital;
    if (attribute_exists(sys, es_path::attachment_density)) mask |= ElectronicStructureSourceFlag_TransitionDensity;
    if (density_properties_gather(nullptr, 0, sys) > 0)        mask |= ElectronicStructureSourceFlag_DensityProperty;
    // ELF and LOL take the occupied orbitals' gradients: the orbitals with their occupations, and a basis whose
    // derivatives are within what is evaluated (no g functions)
    if (attribute_exists(sys, es_path::alpha_coefficient) && md_qm_kinetic_available(&sys)) {
        mask |= ElectronicStructureSourceFlag_ElectronLocalization | ElectronicStructureSourceFlag_LocalizedOrbitalLocator;
    }
    // The fields of the system as a whole, where it has what they are made of (ELF and LOL above: the same)
    for (int i = 0; i < (int)ElectronicStructureSource::Count; ++i) {
        const ElectronicStructureSource src = (ElectronicStructureSource)i;
        if (electronic_structure_source_is_surface_field(src) && surface_field_available(electronic_structure_source_field_kind(src), sys)) {
            mask |= electronic_structure_source_flag(src);
        }
    }

    return mask;
}

bool electronic_structure_select_available_source(ElectronicStructureRepresentation* es, const md_system_t& sys) {
    MD_ASSERT(es);
    const ElectronicStructureSourceFlags mask = es_source_mask(sys);
    if (mask == 0) {
        return false;
    }

    if (!electronic_structure_source_supported(mask, es->shown.source)) {
        for (int n = 0; n < (int)ElectronicStructureSource::Count; ++n) {
            const ElectronicStructureSource source = (ElectronicStructureSource)n;
            if (electronic_structure_source_supported(mask, source)) {
                es->shown.source = source;
                electronic_structure_set_source_defaults(es);
                break;
            }
        }
    }

    // The same rule the representation window applies when it draws: a key naming none of the
    // properties this system has is replaced by the first one.
    if (es->shown.source == ElectronicStructureSource::DensityProperty) {
        DensityProperty props[DensityProperty::MaxCount];
        const size_t num_props = MD_MIN(density_properties_gather(props, MD_COUNTOF(props), sys), MD_COUNTOF(props));
        bool found = false;
        for (size_t i = 0; i < num_props; ++i) {
            if (props[i].key == es->shown.density_property_key) {
                found = true;
                break;
            }
        }
        if (!found && num_props > 0) {
            es->shown.density_property_key = props[0].key;
        }
    }

    return true;
}


// Per atom scalar fields are deliberately not gathered anywhere. They live in the system's attribute
// table under atom/, whoever loaded the data put them there, and the UI reads that table directly
// through atom_attribute_query.

int playback_stride(const ApplicationState* app) {
    MD_ASSERT(app);
    return (app->animation.mode == PlaybackMode::Playing) ? MD_MAX(1, app->animation.stride) : 1;
}

void interpolate_system_state(SystemInstance& inst) {
    SystemContext& sc = *inst.sys;
    ApplicationState* app = sc.app;
    MD_ASSERT(app);
	const auto& sys  = sc.data.sys;

    // Without a run, the system's own coordinates, as they are, with the backbone they carry
    SystemRun* run = instance_run(inst);
    if (!run) {
        if (sys.atom.count == 0 || !md_system_state_has_coords(&sc.data.state)) return;
        md_system_state_copy(&inst.state, &sc.data.state);
        if (md_system_state_backbone_compute(&inst.state, &sc.data.sys)) {
            upload_secondary_structure_weights(inst);
        }
        inst.state_rotation = md_float4x4_ident();
        md_float3x3_t A;
        md_unitcell_A_extract_float(A.elem, &inst.state.unitcell);
        const md_float3_t c = md_float3x3_mul_float3(A, md_float3_set(0.5f, 0.5f, 0.5f));
        inst.unitcell_transform = md_float4x4_translate(-c.x, -c.y, -c.z);
        inst.dirty_gpu_buffers |= MolBit_DirtyPosition;
        return;
    }

    // At a frame of its own, or at the time shown. Recorded before anything can return early, so that a change of
    // stride is acted on once (see main)
    const bool fixed = inst.frame >= 0;
    const int64_t stride = fixed ? 1 : playback_stride(app);
    inst.last_interpolated_stride = (int)stride;

    // Room for the coordinates of every atom: an instance made while the system had none (a topology, before its first
    // run) has none of its own yet
    if (inst.state.num_atoms != sys.atom.count || !inst.state.xyz) {
        md_system_state_init(&inst.state, sys.atom.count);
    }
    const size_t num_atoms = inst.state.num_atoms;
    const size_t num_frames = run_num_frames(sc, run);
    if (num_atoms == 0 || num_frames == 0) return;

    const int64_t last_frame = MD_MAX(0LL, (int64_t)num_frames - 1);
    // Not a time: the fractional frame of the run, its own or the one of the clock's time (clock.h)
    const double time = MD_CLAMP(instance_current_frame(inst), 0.0, double(last_frame));

    // The frames interpolated between, the keys: every stride:th frame counted from the first, and the last, which
    // closes a shorter final interval when the stride does not divide the trajectory. With a stride of 1 every frame
    // is a key. Counting from frame 0 rather than from where playback started keeps the path the same every time
    // it is played.
    const int64_t last_key = (last_frame + stride - 1) / stride;
    auto key_frame = [stride, last_key, last_frame](int64_t k) {
        return MD_MIN(MD_CLAMP(k, (int64_t)0, last_key) * stride, last_frame);
    };
    const int64_t key = MD_MIN((int64_t)(time / (double)stride), last_key);

    // This represents the frames that we would like to load into memory for interpolation (worst case).
    const int64_t frames[4] = {
        key_frame(key - 1),
        key_frame(key),
        key_frame(key + 1),
        key_frame(key + 2),
    };

    // Where in the interval between the two middle keys the time is, and which of them is nearest
    const double t = (frames[2] > frames[1]) ? (time - (double)frames[1]) / (double)(frames[2] - frames[1]) : 0.0;
    const int64_t nearest_frame = (t < 0.5) ? frames[1] : frames[2];

    if (app->animation.interpolation == InterpolationMode::Nearest) {
        if (inst.last_interpolated_nearest_frame == nearest_frame) {
            return;
        }
        inst.dirty_gpu_buffers |= MolBit_ClearVelocity;
    }
    inst.last_interpolated_nearest_frame = nearest_frame;

    // The backbone keeps the orientation of its cross sections coherent from one displayed state to the next. That
    // continuity only means something while the structure moves continuously, so a jump (seeking, skipping frames,
    // a new run) starts it over. Playback keeps it whatever the speed.
    {
        const double last = inst.last_interpolated_frame;
        const bool playing = app->animation.mode == PlaybackMode::Playing;
        if (last < 0.0 || (!playing && fabs(time - last) > 1.5)) {
            inst.dirty_gpu_buffers |= MolBit_ResetBackboneHistory;
        }
        inst.last_interpolated_frame = time;
    }

    // The number of atoms to be processed per thread when divided into chunks
    const uint32_t grain_size = 1024;

    md_allocator_i* temp_arena = app->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(temp_arena);
    MD_DEFER { md_temp_end(temp); };

    struct Payload {
        SystemContext*  sys;
        SystemInstance* inst;
        SystemRun*      run;
        float t;
        float w[4];     // The weights of the four frames, for the cubic modes
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];

        md_system_state_t* src_states[4];
		md_system_state_t* dst_state;

        // The backbone of the destination state, and the cartoon's weights for it (temporary, uploaded at the end)
        md_backbone_angles_t*        dst_angle;
        md_secondary_structure_t*    dst_ss;
        md_gl_secondary_structure_t* ss_weights;

        md_float4x4_t recenter_transform;
    };

    const InterpolationMode mode = (frames[1] == frames[2]) ? InterpolationMode::Nearest : app->animation.interpolation;

    Payload payload = {
        .sys = &sc,
        .inst = &inst,
        .run = run,
        .t = (float)t,
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .dst_state = &inst.state,
    };

    // The cubic modes are a weighted sum of the four frames, the same for everything interpolated: coordinates, cell and
    // backbone angles.
    if (mode == InterpolationMode::CubicSpline) {
        // Tangent scale of the cardinal spline: (p2 - p0) * s. Tension 0 gives s = 0.5, i.e. Catmull-Rom, which plays
        // uniform motion back uniformly; tension 1 eases in and out of every frame.
        const float s = 0.5f * (1.0f - MD_CLAMP(app->animation.tension, 0.0f, 1.0f));
        md_geom_cubic_spline_weights(payload.w, payload.t, s);
    } else if (mode == InterpolationMode::BSpline) {
        md_geom_cubic_bspline_weights(payload.w, payload.t);
        // At the ends of the trajectory the missing neighbour is the reflection of the next frame through the end one
        // (p[-1] = 2 p[0] - p[1]) rather than the end frame repeated: that is what makes the B-spline begin and end on
        // the first and last frames, which repeating them would leave it short of. Within the trajectory it lands on none.
        float* w = payload.w;
        if (frames[0] == frames[1]) { w[1] += 2.0f * w[0]; w[2] -= w[0]; w[0] = 0.0f; }
        if (frames[3] == frames[2]) { w[2] += 2.0f * w[3]; w[1] -= w[3]; w[3] = 0.0f; }
    }

    // The backbone of the displayed state is the run's at the frames around it: the angles and the secondary structure
    // go into the state's attributes, the cartoon's weights into temporary memory, which is uploaded at the end.
    const size_t num_segments = sc.data.sys.protein_backbone.segment.count;
    if (num_segments > 0 && run->trajectory_data.backbone_angles.data) {
        payload.dst_angle = md_system_state_backbone_angles_write(&inst.state, &sc.data.sys);
    }
    if (num_segments > 0 && run->trajectory_data.secondary_structure.data) {
        payload.dst_ss = md_system_state_secondary_structure_write(&inst.state, &sc.data.sys);
        payload.ss_weights = payload.dst_ss ? md_temp_alloc_array(temp, md_gl_secondary_structure_t, num_segments) : nullptr;
    }

    // Stamp the destination with the frame it is about to represent. The interpolated state is not
    // written by a run extraction, so nothing else would fill this in, and a
    // stale value is worse than an absent one. Nearest snaps to a whole frame; the other modes land
    // between two, which is exactly what the fractional part is for.
    inst.state.frame = (mode == InterpolationMode::Nearest) ? (double)nearest_frame : time;

    // Fresh coordinates, in the lattice frame of the cell: the turn the previous ones carried goes
    // with them. The System State Changed handler puts one back if the orientation is kept.
    inst.state_rotation = md_float4x4_ident();

    int requested_frames[4] = { 0 };
    int num_requested_frames = 0;

    switch (mode) {
        case InterpolationMode::Nearest:
            requested_frames[num_requested_frames++] = (int)nearest_frame;
            break;
        case InterpolationMode::Linear:
            requested_frames[num_requested_frames++] = (int)frames[1];
            requested_frames[num_requested_frames++] = (int)frames[2];
            break;
        case InterpolationMode::CubicSpline:
        case InterpolationMode::BSpline:
            requested_frames[num_requested_frames++] = (int)frames[0];
            requested_frames[num_requested_frames++] = (int)frames[1];
            requested_frames[num_requested_frames++] = (int)frames[2];
            requested_frames[num_requested_frames++] = (int)frames[3];
            break;
        default:
            MD_ASSERT(false);
            break;
    }

    // Represents the frame cache slot indices for the requested frames, -1 if not present in cache
    int frame_cache_slot_idx[4] = { -1, -1, -1, -1 };

    // Array of frame_cache slots which requires a load
    int frame_cache_load_slot[4] = {0};
    int num_frames_to_load = 0;

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = -1;
        if (!find_frame_in_cache(&slot_idx, requested_frames[i], &inst.frame_cache)) {
            slot_idx = find_lru_cache_slot(&inst.frame_cache);
            frame_cache_load_slot[num_frames_to_load++] = slot_idx;
        }
        set_mru_cache_slot(&inst.frame_cache, slot_idx);
        inst.frame_cache.frame_idx[slot_idx] = requested_frames[i];
        frame_cache_slot_idx[i] = slot_idx;
    }

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = frame_cache_slot_idx[i];
        payload.src_states[i] = &inst.frame_cache.states[slot_idx];
    }

    // The steps in order, each spread over the pool where it is a range; the frame waits for them all (frame priority:
    // the main thread never takes up background work while it does)
    Payload* data = &payload;

    task_system::parallel_for(num_frames_to_load, 1,
        [data = &payload, frame_cache_load_slot](uint32_t range_beg, uint32_t range_end) {
            const uint32_t thread_num = task_system::thread_index();
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int slot_idx = frame_cache_load_slot[i];
                int frame_idx = data->inst->frame_cache.frame_idx[slot_idx];
                md_system_state_t* state = &data->inst->frame_cache.states[slot_idx];
                // This thread's own context, kept for as long as the trajectory, so playback opens
                // the run's files once rather than once per frame.
                md_system_extract_t** ex = (thread_num < MD_COUNTOF(data->run->frame_extract)) ? &data->run->frame_extract[thread_num] : nullptr;
                if (ex && !*ex) {
                    *ex = md_system_extract_begin(&data->sys->data.sys, md_str_from_cstr(data->run->path),
                        frame_extract_paths, MD_COUNTOF(frame_extract_paths), md_get_heap_allocator());
                }
                if (ex && *ex) {
                    md_system_extract_frame(*ex, frame_idx, state);
                } else {
                    extract_frame(*data->sys, data->run, frame_idx, state);
                }
            }
        }
    );


    switch (mode) {
        case InterpolationMode::Nearest: {
            {
                data->dst_state->unitcell = data->src_states[0]->unitcell;
                md_memcpy(data->dst_state->xyz, data->src_states[0]->xyz, sizeof(md_float3_t) * data->dst_state->num_atoms);
            }
            break;
        }
        case InterpolationMode::Linear: {
            {
                double x  = lerp(data->src_states[0]->unitcell.x,  data->src_states[1]->unitcell.x,  data->t);
                double y  = lerp(data->src_states[0]->unitcell.y,  data->src_states[1]->unitcell.y,  data->t);
                double z  = lerp(data->src_states[0]->unitcell.z,  data->src_states[1]->unitcell.z,  data->t);
                double xy = lerp(data->src_states[0]->unitcell.xy, data->src_states[1]->unitcell.xy, data->t);
                double xz = lerp(data->src_states[0]->unitcell.xz, data->src_states[1]->unitcell.xz, data->t);
                double yz = lerp(data->src_states[0]->unitcell.yz, data->src_states[1]->unitcell.yz, data->t);
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
			}

            task_system::parallel_for((uint32_t)num_atoms, grain_size, [data = &payload](uint32_t range_beg, uint32_t range_end) {
                size_t count = range_end - range_beg;
                md_float3_t* dst = data->dst_state->xyz + range_beg;
                const md_float3_t* src[2] = { data->src_states[0]->xyz + range_beg, data->src_states[1]->xyz + range_beg };

                md_geom_interpolate_linear(dst, src, count, &data->dst_state->unitcell, data->t);
            });


            break;
        }
        case InterpolationMode::CubicSpline:
        case InterpolationMode::BSpline: {
            {
                const md_unitcell_t* c[4] = { &data->src_states[0]->unitcell, &data->src_states[1]->unitcell, &data->src_states[2]->unitcell, &data->src_states[3]->unitcell };
                const float* w = data->w;
                double x  = w[0] * c[0]->x  + w[1] * c[1]->x  + w[2] * c[2]->x  + w[3] * c[3]->x;
                double y  = w[0] * c[0]->y  + w[1] * c[1]->y  + w[2] * c[2]->y  + w[3] * c[3]->y;
                double z  = w[0] * c[0]->z  + w[1] * c[1]->z  + w[2] * c[2]->z  + w[3] * c[3]->z;
                double xy = w[0] * c[0]->xy + w[1] * c[1]->xy + w[2] * c[2]->xy + w[3] * c[3]->xy;
                double xz = w[0] * c[0]->xz + w[1] * c[1]->xz + w[2] * c[2]->xz + w[3] * c[3]->xz;
                double yz = w[0] * c[0]->yz + w[1] * c[1]->yz + w[2] * c[2]->yz + w[3] * c[3]->yz;
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
            }

            task_system::parallel_for((uint32_t)num_atoms, grain_size, [data = &payload](uint32_t range_beg, uint32_t range_end) {
                size_t count = range_end - range_beg;
                md_float3_t* dst = data->dst_state->xyz + range_beg;
                const md_float3_t* src[4] = { data->src_states[0]->xyz + range_beg, data->src_states[1]->xyz + range_beg, data->src_states[2]->xyz + range_beg, data->src_states[3]->xyz + range_beg };

                md_geom_interpolate_cubic_weighted(dst, src, count, &data->dst_state->unitcell, data->w);
            });

            
            break;
        }
        default:
            MD_ASSERT(false);
            break;
    }

    if (payload.dst_angle) {
        switch (mode) {
            case InterpolationMode::Nearest: {
                {
                    const md_backbone_angles_t* src_angles[2] = {
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    md_memcpy(data->dst_angle, src_angle, data->sys->data.sys.protein_backbone.segment.count * sizeof(md_backbone_angles_t));
                }

                break;
            }
            case InterpolationMode::Linear: {
                task_system::parallel_for((uint32_t)sys.protein_backbone.segment.count, 1, [data = &payload](uint32_t range_beg, uint32_t range_end) {
                    const md_backbone_angles_t* src_angles[2] = {
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[2] = {src_angles[0][i].phi, src_angles[1][i].phi};
                        float psi[2] = {src_angles[0][i].psi, src_angles[1][i].psi};

                        phi[1] = md_f32_deperiodize_ortho(phi[1], phi[0], (float)MD_TWO_PI);
                        psi[1] = md_f32_deperiodize_ortho(psi[1], psi[0], (float)MD_TWO_PI);

                        float final_phi = lerp(phi[0], phi[1], data->t);
                        float final_psi = lerp(psi[0], psi[1], data->t);
                        data->dst_angle[i] = {md_f32_deperiodize_ortho(final_phi, 0, (float)MD_TWO_PI), md_f32_deperiodize_ortho(final_psi, 0, (float)MD_TWO_PI)};
                    }
                });

                break;
            }
            case InterpolationMode::CubicSpline:
            case InterpolationMode::BSpline: {
                task_system::parallel_for((uint32_t)sys.protein_backbone.segment.count, 1, [data = &payload](uint32_t range_beg, uint32_t range_end) {
                    const md_backbone_angles_t* src_angles[4] = {
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[0],
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[2],
                        data->run->trajectory_data.backbone_angles.data + data->run->trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[4] = {src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi};
                        float psi[4] = {src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi};

                        phi[0] = md_f32_deperiodize_ortho(phi[0], phi[1], (float)MD_TWO_PI);
                        phi[2] = md_f32_deperiodize_ortho(phi[2], phi[1], (float)MD_TWO_PI);
                        phi[3] = md_f32_deperiodize_ortho(phi[3], phi[2], (float)MD_TWO_PI);

                        psi[0] = md_f32_deperiodize_ortho(psi[0], psi[1], (float)MD_TWO_PI);
                        psi[2] = md_f32_deperiodize_ortho(psi[2], psi[1], (float)MD_TWO_PI);
                        psi[3] = md_f32_deperiodize_ortho(psi[3], psi[2], (float)MD_TWO_PI);

                        const float* w = data->w;
                        float final_phi = w[0] * phi[0] + w[1] * phi[1] + w[2] * phi[2] + w[3] * phi[3];
                        float final_psi = w[0] * psi[0] + w[1] * psi[1] + w[2] * psi[2] + w[3] * psi[3];
                        data->dst_angle[i] = {md_f32_deperiodize_ortho(final_phi, 0, (float)MD_TWO_PI), md_f32_deperiodize_ortho(final_psi, 0, (float)MD_TWO_PI)};
                    }
                });

                break;
            }
            default:
                MD_ASSERT(false);
                break;
        }
    }

    if (payload.dst_ss) {
        task_system::parallel_for((uint32_t)num_segments, 1, [data = &payload, mode](uint32_t range_beg, uint32_t range_end) {
            // The state carries the secondary structure of the nearest frame as it was assigned. The cartoon blends the
            // denoised copy (a presentation smoothing, see secondary_structure_render) where it exists.
            const md_secondary_structure_t* label_nearest = data->run->trajectory_data.secondary_structure.data +
                data->run->trajectory_data.secondary_structure.stride * (data->t < 0.5f ? data->frames[1] : data->frames[2]);
            const md_secondary_structure_t* ss_data = data->run->trajectory_data.secondary_structure_render.data ?
                data->run->trajectory_data.secondary_structure_render.data :
                data->run->trajectory_data.secondary_structure.data;
            const size_t ss_stride = data->run->trajectory_data.secondary_structure_render.data ?
                data->run->trajectory_data.secondary_structure_render.stride :
                data->run->trajectory_data.secondary_structure.stride;
            const md_secondary_structure_t* src_ss[4] = {
                ss_data + ss_stride * data->frames[0],
                ss_data + ss_stride * data->frames[1],
                ss_data + ss_stride * data->frames[2],
                ss_data + ss_stride * data->frames[3],
            };
            const md_secondary_structure_t* src_ss_nearest = data->t < 0.5f ? src_ss[1] : src_ss[2];
            auto normalize_ss = [](md_gl_secondary_structure_t ss) {
                ss.helix = MD_CLAMP(ss.helix, 0.0f, 1.0f);
                ss.sheet = MD_CLAMP(ss.sheet, 0.0f, 1.0f);
                float sum = ss.helix + ss.sheet;
                if (sum > 1.0f) {
                    ss.helix /= sum;
                    ss.sheet /= sum;
                }
                return ss;
            };
            auto blend_ss = [normalize_ss](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b, float t) {
                md_gl_secondary_structure_t ss = {
                    .helix = lerp(a.helix, b.helix, t),
                    .sheet = lerp(a.sheet, b.sheet, t),
                };
                return normalize_ss(ss);
            };
            auto smoothstep5 = [](float t) {
                t = MD_CLAMP(t, 0.0f, 1.0f);
                return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
            };
            auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
                return a.helix == b.helix && a.sheet == b.sheet;
            };
            switch (mode) {
            default:
                MD_LOG_DEBUG("Unsupported interpolation mode for secondary structure interpolation");
                [[fallthrough]];
            case InterpolationMode::Nearest: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    data->dst_ss[i] = label_nearest[i];
                    data->ss_weights[i] = md_gl_secondary_structure_convert(src_ss_nearest[i]);
                }
                break;
            }
            case InterpolationMode::Linear: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[2] = { src_ss[1][i], src_ss[2][i] };
                    md_gl_secondary_structure_t ss_gl[2] = { md_gl_secondary_structure_convert(ss[0]), md_gl_secondary_structure_convert(ss[1]) };
                    data->dst_ss[i] = label_nearest[i];
                    data->ss_weights[i] = blend_ss(ss_gl[0], ss_gl[1], data->t);
                }
                break;
            }
            case InterpolationMode::CubicSpline:
            case InterpolationMode::BSpline: {
                
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[4] = { src_ss[0][i], src_ss[1][i], src_ss[2][i], src_ss[3][i] };
                    
                    md_gl_secondary_structure_t ss_gl[4] = {
                        md_gl_secondary_structure_convert(ss[0]),
                        md_gl_secondary_structure_convert(ss[1]),
                        md_gl_secondary_structure_convert(ss[2]),
                        md_gl_secondary_structure_convert(ss[3]),
                    };

                    // Cleanup isolated temporal assignments to reduce noise during transitions.
                    if (is_eq(ss_gl[0], ss_gl[2]) && !is_eq(ss_gl[1], ss_gl[0])) {
                        ss_gl[1] = ss_gl[0];
                    }
                    if (is_eq(ss_gl[1], ss_gl[3]) && !is_eq(ss_gl[2], ss_gl[1])) {
                        ss_gl[2] = ss_gl[1];
                    }

                    data->dst_ss[i] = label_nearest[i];
                    data->ss_weights[i] = blend_ss(ss_gl[1], ss_gl[2], smoothstep5(data->t));
                }
                break;
            }
            }
        });

        // Isolated coils between matching structured segments are filled in, to reduce the noise during transitions.
        // A non temporal filtering step, after the blend above.
        {
            secondary_structure_weights_fill_isolated_coils(data->ss_weights, &data->sys->data.sys.protein_backbone);
        }
    }

    // The weights are the renderer's input and nothing else's: straight to the GPU, from the temporary memory above
    if (payload.ss_weights) {
        md_gl_mol_set_backbone_secondary_structure(inst.gl_mol, 0, (uint32_t)num_segments, payload.ss_weights, 0);
    }

    // unitcell transform is essentially just a translation to place the center of the unitcell at the origin
    md_float3x3_t A;
    md_unitcell_A_extract_float(A.elem, &inst.state.unitcell);
    md_float3_t c = md_float3x3_mul_float3(A, md_float3_set(0.5f, 0.5f, 0.5f));
    inst.unitcell_transform = md_float4x4_translate(-c.x, -c.y, -c.z);

#if 0
    if (sys.unitcell.flags) {
        md_float3_t c = sys.unitcell.basis * md_float3_set1(0.5f);
        sc.data.model_mat = md_float4x4_translate_float3(-c);
    }
#endif

    inst.dirty_gpu_buffers |= MolBit_DirtyPosition;
}

// Exact, element by element. The identities compared against here are only ever assigned, never computed.
static bool is_identity(const md_float4x4_t& M) {
    const md_float4x4_t I = md_float4x4_ident();
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            if (M.elem[c][r] != I.elem[c][r]) return false;
        }
    }
    return true;
}

void recenter_mark_query_dirty(SystemContext& sys) {
    sys.operations.recenter_query.version += 1;
    if (sys.operations.recenter_query.version == 0) {
        sys.operations.recenter_query.version = 1;
    }
}

const md_bitset_t& recenter_get_active_target_mask(const SystemContext& sys) {
    return sys.operations.recenter_query.enabled ? sys.operations.recenter_query.mask : sys.operations.selection_mask;
}

bool recenter_update_query_mask(SystemContext& sys) {
    ApplicationState* state = sys.app;
    MD_ASSERT(state);

    auto& query = sys.operations.recenter_query;
    const uint64_t ir_fingerprint = state->script.ir ? state->script.ir_fingerprint : 0;
    if (query.ir_fingerprint != ir_fingerprint) {
        query.ir_fingerprint = ir_fingerprint;
        recenter_mark_query_dirty(sys);
    }

    if (!query.enabled) {
        return false;
    }

    if (query.filter.result.is_dynamic) {
        recenter_mark_query_dirty(sys);
    }

    if (query.evaluated_version == query.version) {
        return false;
    }

    filter_evaluate(&query.filter, &query.mask, state->script.language, filter_context(sys), state->allocator.persistent);
    query.evaluated_version = query.version;
    return true;
}

void recenter_update(SystemContext& sys) {
    recenter_update_query_mask(sys);
    recenter_update_target_data(sys);
}

// Exactly the same atoms, whatever form each set is in
static bool same_atoms(const md_bitset_t* a, const md_bitset_t* b) {
    return md_bitset_equal(*a, *b);
}

void recenter_update_target_data(SystemContext& sys) {
    ApplicationState* state = sys.app;
    // The reference is of the system's own coordinates, which every instance of it is fitted to: what superposes them
    if (!md_system_state_has_coords(&sys.data.state)) return;

    // The reference is only for keeping the orientation. Centering alone needs nothing from the first frame,
    // so without this a target picked by a query which depends on the frame read it from disk in every frame.
    if (!sys.operations.fixate_orientation) return;

    auto& ref = sys.operations.initial_frame;
    const md_bitset_t& target_mask = recenter_get_active_target_mask(sys);

    // Keyed on the atoms themselves. The versions this used to compare change every time a query is
    // evaluated, which for one that depends on the frame is every frame, whether or not its atoms changed.
    if (!ref.valid || !same_atoms(&ref.target_mask, &target_mask)) {
        md_bitset_assign(&ref.target_mask, target_mask, state->allocator.persistent);
        ref.valid = true;
        size_t count = md_bitset_count(target_mask);

        md_array_resize(sys.operations.initial_frame.rel_xyzw, count, state->allocator.persistent);

        sys.operations.initial_frame.com = md_float3_zero();

        if (sys.operations.initial_frame.rel_xyzw && count > 0) {
            // The system's own coordinates: those its files gave, or the first frame of its first run
            const md_system_state_t& temp_state = sys.data.state;
            const md_float3_t* temp_xyz = temp_state.xyz;

            md_bitset_iter_t it = md_bitset_iter(target_mask);
            int dst_idx = 0;
            while (md_bitset_iter_next(&it)) {
                uint64_t src_idx = it.idx;
                float mass = md_atom_mass(&sys.data.sys.atom, src_idx);
                sys.operations.initial_frame.rel_xyzw[dst_idx++] = md_float4_from_float3(temp_xyz[src_idx], mass);
            }

            // Mutually consistent images and the plain weighted mean of them, then relative to that mean.
            // The circular mean (md_geom_com with a cell) is not the mean of the placed
            // points, and a fit against coordinates which are not centred is biased.
            md_float3_t com = md_float3_zero();
            md_geom_deperiodize_self(sys.operations.initial_frame.rel_xyzw, count, &temp_state.unitcell, &com);
            const md_float4_t com4 = md_float4_from_float3(com, 0);
            for (size_t i = 0; i < count; ++i) {
                sys.operations.initial_frame.rel_xyzw[i] = md_float4_sub(sys.operations.initial_frame.rel_xyzw[i], com4);
            }
            sys.operations.initial_frame.com = com;
        }
    }
}

bool recenter_calculate_transform(md_float4x4_t* translation, md_float4x4_t* rotation, const SystemInstance& inst) {
    const SystemContext& sys = *inst.sys;
    const ApplicationState* app = sys.app;
    MD_ASSERT(translation);
    MD_ASSERT(rotation);
    MD_ASSERT(app);

    *translation = md_float4x4_ident();
    *rotation    = md_float4x4_ident();
    bool turn = false;

    const md_bitset_t& target_mask = recenter_get_active_target_mask(sys);
    size_t count = md_bitset_count(target_mask);

    if (count > 0) {
        md_temp_scope_t temp = md_temp_begin_in(app->allocator.frame);
        MD_DEFER { md_temp_end(temp); };

        // Extract xyzw subset of target
        md_float4_t* target_xyzw = md_temp_alloc_array(temp, md_float4_t, count);

		extract_xyzw_from_mask(target_xyzw, &target_mask, &sys.data.sys, &inst.state);

        // Calculate target
        md_float3_t target = {0};
        if (md_unitcell_flags(&inst.state.unitcell) != 0) {
            md_float3x3_t A = {0};
            md_unitcell_A_extract_float(A.elem, &inst.state.unitcell);
            target = md_float3x3_mul_float3(A, md_float3_set1(0.5f));
        } 

        // Place the target in mutually consistent images and take its centre. When there is a
        // reference to hold the orientation against, the rotation, the centre and each point's
        // periodic image are solved for together, so the images are chosen to minimise the alignment
        // residual rather than being committed to beforehand by a criterion unrelated to it.
        //
        // Both calls report the centre in the image the STATE coordinates actually occupy, which is
        // what makes the translation below valid - it is applied to those same coordinates, untouched.
        // A centre folded into the reference cell would place a target living outside that cell one
        // lattice vector off, and with a rotation in play, off by R times a lattice vector.
        md_float3x3_t R = md_float3x3_ident();
        md_float3_t target_com = md_float3_zero();

        // The reference has to have been built from the SAME target that is being fitted now.
        // A size match is not sufficient: the selection can change to a different set of equal
        // size between recenter_update() and this call, which would silently pair up unrelated
        // atoms and yield a garbage rotation. The atoms themselves are the identity of the target.
        const bool reference_valid =
            sys.operations.fixate_orientation &&
            sys.operations.initial_frame.valid &&
            sys.operations.initial_frame.rel_xyzw &&
            md_array_size(sys.operations.initial_frame.rel_xyzw) == count &&
            same_atoms(&sys.operations.initial_frame.target_mask, &target_mask);

        // R maps the CURRENT target onto the reference: R * (q - target_com) ~= p. The relative fit this
        // replaced had its operands the other way around, which yields the rotation carrying the reference
        // onto the current frame - applied to the current frame, it doubled the rotation it was meant to
        // cancel. It also centred on the circular mean, which lands in the reference cell whatever image
        // the target occupies and is not the mean of the placed points.
        if (sys.operations.fixate_orientation && reference_valid) {
            // The reference is stored relative to its own centre, so its centre here is the origin
            md_geom_optimal_rotation_pbc(&R, &target_com, target_xyzw, sys.operations.initial_frame.rel_xyzw, md_float3_zero(),
                                                   target_xyzw, count, &inst.state.unitcell, 8, 1.0e-6f);
            R = md_float3x3_orthonormalize(R);
        } else {
            md_geom_deperiodize_self(target_xyzw, count, &inst.state.unitcell, &target_com);
        }

        // Split at the cell centre. The translation preserves the lattice, so it is valid for coordinates
        // in whatever image they arrived in. The turn does not: R times a lattice vector is not a lattice
        // vector, so it may only be applied once every image is settled relative to the target. The
        // product is the transform this used to return in one piece,
        //     translate(target) * A * R * translate(-target_com)
        const md_float4x4_t A = sys.operations.alignment_mat;
        turn = (sys.operations.fixate_orientation && reference_valid) || !is_identity(A);

        *translation = md_float4x4_translate_float3(md_float3_sub(target, target_com));
        if (turn) {
            *rotation = md_float4x4_translate_float3(target) * A * md_float4x4_from_float3x3(R) * md_float4x4_translate_float3(-target);
        }
    }
    return turn;
}

// Every atom of the instance's coordinates through M
static void state_transform(SystemInstance& inst, const md_float4x4_t& M) {
    md_system_state_t& s = inst.state;
    task_system::parallel_for((uint32_t)s.num_atoms, 1024, [&s, M](uint32_t range_beg, uint32_t range_end) {
        md_float4x4_batch_transform_inplace(s.xyz + range_beg, 1.0f, range_end - range_beg, M);
    });
}

bool apply_state_operations(SystemInstance& inst, bool recenter, bool pbc, bool unwrap) {
    SystemContext& sc = *inst.sys;
    md_system_t& sys = sc.data.sys;
    md_system_state_t& s = inst.state;

    if (s.num_atoms == 0 || !s.xyz) return false;
    if (!recenter && !pbc && !unwrap) return false;

    const bool periodic = md_unitcell_flags(&s.unitcell) != 0;

    // Into the lattice frame first. Wrapping or making whole turned coordinates against the unturned
    // cell moves atoms by vectors that are not lattice vectors, to places where they have no image.
    const md_float4x4_t prev_rotation = inst.state_rotation;
    if (!is_identity(prev_rotation)) {
        state_transform(inst, md_float4x4_inverse(prev_rotation));
    }

    // A wrap or a make whole on its own keeps the turn the coordinates had
    md_float4x4_t rotation = prev_rotation;
    bool fresh_turn = false;
    if (recenter && !md_bitset_empty(recenter_get_active_target_mask(sc))) {
        // Current for this target before it is fitted against. The main loop does this every frame too, but
        // ticking 'keep orientation' applies at once, before the loop has had a chance to build it.
        recenter_update_target_data(sc);

        // Measured on the lattice frame coordinates: the fit is against the reference, not against a
        // previous fit's output
        md_float4x4_t translation = md_float4x4_ident();
        fresh_turn = recenter_calculate_transform(&translation, &rotation, inst);
        state_transform(inst, translation);
    }
    const bool turn = !is_identity(rotation);

    // Settle the periodic images before turning: every atom into the image nearest the cell centre,
    // which is where the translation just put the target. Without it the turn carries the image each
    // atom happened to be written in by the trajectory into the result, as a displacement of R times
    // a lattice vector. That differs from frame to frame for every atom crossing the trajectory's own
    // box, and for the whole system at once whenever the target's centre lands in another image.
    if (periodic && (pbc || fresh_turn)) {
        task_system::parallel_for((uint32_t)s.num_atoms, 1, [&s](uint32_t range_beg, uint32_t range_end) {
            md_geom_pbc((float*)(s.xyz + range_beg), 0, NULL, range_end - range_beg, &s.unitcell);
        });
    }

    if (unwrap) {
        const size_t num_molecules = md_system_molecule_count(&sys);
        if (num_molecules > 0) {
            task_system::parallel_for((uint32_t)num_molecules, 1, [&s, &sys](uint32_t range_beg, uint32_t range_end) {
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    const md_molecule_t mol = md_system_molecule(&sys, i);
                    md_geom_unwrap((float*)s.xyz, sizeof(md_float3_t), mol.atom_idx, mol.parent_idx, mol.count, &s.unitcell);
                }
            });
        }
    }

    if (turn) {
        state_transform(inst, rotation);
    }
    inst.state_rotation = turn ? rotation : md_float4x4_ident();

    return true;
}

// The coordinates bonds are found in: the whole frame of the run nearest frame, or the coordinates shown when there is
// no run, taken back out of any turn from keeping the orientation. The bond search takes its minimum images in the cell,
// so they have to be in the lattice frame of the cell, which a frame of the run always is.
static bool bond_coordinates(md_system_state_t* state, const SystemContext& sys, int64_t frame) {
    const size_t num_atoms = sys.data.sys.atom.count;
    if (!md_system_state_init(state, num_atoms)) return false;

    const size_t num_frames = run_num_frames(sys);
    if (num_frames > 0) {
        if (!extract_frame(sys, MD_CLAMP(frame, (int64_t)0, (int64_t)num_frames - 1), state)) {
            MD_LOG_ERROR("Failed to extract frame data");
            return false;
        }
        return true;
    }
    const md_system_state_t& shown = system_instance(sys).state;
    if (!shown.xyz || shown.num_atoms != num_atoms) return false;
    md_memcpy(state->xyz, shown.xyz, num_atoms * sizeof(md_float3_t));
    state->unitcell = shown.unitcell;
    if (!is_identity(system_instance(sys).state_rotation)) {
        md_float4x4_batch_transform_inplace(state->xyz, 1.0f, num_atoms, md_float4x4_inverse(system_instance(sys).state_rotation));
    }
    return true;
}

bool recompute_bonds(SystemContext& sc, int64_t frame) {
    md_system_t& sys = sc.data.sys;
    if (sys.atom.count == 0) return false;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    md_system_state_t state = { .alloc = temp.arena };
    if (!bond_coordinates(&state, sc, frame)) return false;

    // All of it inferred again from these coordinates (copied by the work), which become the reference: in the
    // background, as after an edit of the bonds
    system_inference_start(sc, &state);
    return true;
}

void system_core_edit_begin(SystemContext& sys) {
    ApplicationState* app = sys.app;
    MD_ASSERT(app);
    SystemInference& inf = sys.data.inference;
    script_evals_interrupt(app);
    // The work of an inference under way shares the core: stopped, and kept for the coordinates it was to be inferred
    // from (system_core_changed starts it again)
    inf.job.cancel_and_wait();
    task_system::wait_all();
    // What the system has derived is of the core before the edit
    inf.pending = true;
    inf.restart = false;
}

void system_core_changed(SystemContext& sys) {
    system_inference_start(sys, nullptr);
}

// Of the bonds the system has now: not of those of an inference before (BondTracking::topology_version)
static bool bond_tracking_valid(const SystemContext& sys) {
    const auto& bt = sys.operations.bond_tracking;
    return bt.enabled && bt.key != 0 && bt.topology_version == sys.data.sys.topology_version;
}

const md_bond_delta_t* bond_tracking_delta(const SystemContext& sys) {
    return bond_tracking_valid(sys) ? &sys.operations.bond_tracking.delta : nullptr;
}

const md_chem_delta_t* bond_tracking_chem(const SystemContext& sys) {
    return bond_tracking_valid(sys) ? &sys.operations.bond_tracking.chem : nullptr;
}

md_atom_flags_t frame_atom_flags(const SystemContext& sc, size_t atom_idx, int* formal_charge, int* hydrogen_count) {
    const md_system_t& sys = sc.data.sys;
    const md_chem_delta_t* chem = bond_tracking_chem(sc);
    if (chem && chem->num_atoms) {
        // The atoms of the region are ascending
        size_t lo = 0, hi = chem->num_atoms;
        while (lo < hi) {
            const size_t mid = (lo + hi) / 2;
            if ((size_t)chem->atom_idx[mid] < atom_idx) lo = mid + 1; else hi = mid;
        }
        if (lo < chem->num_atoms && (size_t)chem->atom_idx[lo] == atom_idx) {
            if (formal_charge)  *formal_charge  = chem->formal_charge[lo];
            if (hydrogen_count) *hydrogen_count = chem->hydrogen_count[lo];
            return chem->atom_flags[lo];
        }
    }
    if (formal_charge)  *formal_charge  = md_atom_formal_charge(&sys.atom, atom_idx);
    if (hydrogen_count) *hydrogen_count = md_atom_hydrogen_count(&sys.atom, atom_idx);
    return md_atom_flags(&sys.atom, atom_idx);
}

bool picked_ring(const SystemContext& sc, size_t idx, const int32_t** atoms, size_t* count, md_ring_flags_t* flags) {
    const md_system_t& sys = sc.data.sys;
    const size_t num = md_ring_data_count(&sys.ring);
    if (idx < num) {
        *atoms = md_ring_atoms(&sys.ring, idx);
        *count = md_ring_size(&sys.ring, idx);
        *flags = md_system_ring_flags(&sys, idx);
        return true;
    }
    const md_chem_delta_t* chem = bond_tracking_chem(sc);
    idx -= num;
    if (chem && idx < md_ring_data_count(&chem->ring)) {
        *atoms = md_ring_atoms(&chem->ring, idx);
        *count = md_ring_size(&chem->ring, idx);
        *flags = md_ring_flags(&chem->ring, idx);
        return true;
    }
    return false;
}

static void bond_tracking_clear(SystemContext& sys) {
    ApplicationState* app = sys.app;
    auto& bt = sys.operations.bond_tracking;
    if (bt.key != 0) {
        md_bond_delta_free(&bt.delta, app->allocator.persistent);
        md_chem_delta_free(&bt.chem, app->allocator.persistent);
        bt.delta = {};
        bt.key = 0;
        bt.formed_capacity = 0;
        bt.ring_capacity = 0;
        system_instance(sys).dirty_gpu_buffers |= MolBit_DirtyBonds;
    }
}

void bond_tracking_update(SystemContext& sc) {
    ApplicationState* app = sc.app;
    MD_ASSERT(app);
    auto& bt = sc.operations.bond_tracking;
    const md_system_t& sys = sc.data.sys;
    if (!bt.enabled || sys.atom.count == 0) {
        bond_tracking_clear(sc);
        return;
    }

    // The scope, evaluated again when the script whose names it may use changes, and in every frame when what it picks
    // depends on the frame. Empty for every atom.
    const bool has_scope = !md_str_empty(md_str_trim(md_str_from_cstr(bt.scope.text)));
    const uint64_t ir_fingerprint = app->script.ir ? app->script.ir_fingerprint : 0;
    if (bt.scope_ir_fingerprint != ir_fingerprint) {
        bt.scope_ir_fingerprint = ir_fingerprint;
        bt.scope_dirty = true;
    }
    if (bt.scope.result.is_dynamic) {
        bt.scope_dirty = true;
    }
    if (bt.scope_dirty) {
        if (has_scope) {
            filter_evaluate(&bt.scope, &bt.scope_mask, app->script.language, filter_context(sc), app->allocator.persistent);
        } else {
            md_bitset_free(&bt.scope_mask, app->allocator.persistent);
            bt.scope.result = {};
            bt.scope.result.success = true;
        }
        bt.scope_dirty = false;
    }
    // A scope that does not evaluate tracks nothing, rather than everything
    if (!bt.scope.result.success) {
        bond_tracking_clear(sc);
        return;
    }

    // What the delta is of. The frame and nothing finer: the bonds of a state between two frames are those of the
    // nearest one, as the script has them.
    const int64_t frame = run_num_frames(sc) > 0 ? (int64_t)(instance_current_frame(system_instance(sc)) + 0.5) : -1;
    uint64_t key = md_hash64(&frame, sizeof(frame), sys.topology_version);
    key = md_hash64(&bt.tolerance, sizeof(bt.tolerance), key);
    key = md_hash64(&bt.breakable, sizeof(bt.breakable), key);
    key = has_scope ? md_bitset_hash64(bt.scope_mask, key) : md_hash64(&has_scope, sizeof(has_scope), key);
    // Without a run the coordinates shown are those of the system, which the operations move only as a whole or into
    // other images: the minimum images the delta is found in are the same
    key = key ? key : 1;
    if (key == bt.key) return;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_system_state_t state = { .alloc = temp.arena };
    md_bond_delta_t delta = {};
    const md_tick_t t0 = md_tick_now();
    if (bond_coordinates(&state, sc, frame)) {
        const md_bond_delta_desc_t desc = {
            .scope = has_scope ? &bt.scope_mask : nullptr,
            .tolerance = bt.tolerance,
            .breakable = bt.breakable,
        };
        md_bond_delta_compute(&delta, &sys, &state, &desc, app->allocator.persistent);
    }
    const double time_ms = md_tick_to_milliseconds(md_tick_now() - t0);

    // Uploaded again only when it changed: most frames of most runs have the bonds of the one before
    const bool same = bt.key != 0 && delta.num_formed == bt.delta.num_formed && md_bitset_equal(delta.broken, bt.delta.broken) &&
        (delta.num_formed == 0 || (md_memcmp(delta.formed_pair, bt.delta.formed_pair, delta.num_formed * sizeof(md_atom_pair_t)) == 0 &&
                                   md_memcmp(delta.formed_flags, bt.delta.formed_flags, delta.num_formed * sizeof(md_bond_flags_t)) == 0));
    md_bond_delta_free(&bt.delta, app->allocator.persistent);
    bt.delta = delta;
    bt.key = key;
    bt.topology_version = sys.topology_version;
    bt.time_ms = time_ms;
    bt.formed_capacity = MD_MAX(bt.formed_capacity, delta.num_formed);
    if (!same) {
        // The chemistry where the bonds changed: what is drawn there, the orders and rings of the frame. A region too
        // large to perceive leaves it empty, and the bonds formed plain.
        if (!md_chem_perceive_delta(&bt.chem, &sys, &state, &bt.delta, MD_CHEM_FLAG_PROTONATE_PH7, app->allocator.persistent)) {
            md_chem_delta_free(&bt.chem, app->allocator.persistent);
        }
        bt.ring_capacity = MD_MAX(bt.ring_capacity, md_ring_data_count(&bt.chem.ring));
        bt.time_ms = md_tick_to_milliseconds(md_tick_now() - t0);
        system_instance(sc).dirty_gpu_buffers |= MolBit_DirtyBonds;
    }
}

// The flags of a bond in the frame's chemistry (bond_tracking_chem), whose bonds are ascending by index: those of the
// region where it is in it, flags otherwise
static md_bond_flags_t frame_bond_flags(const SystemContext& sys, size_t idx, md_bond_flags_t flags) {
    const md_chem_delta_t* chem = bond_tracking_chem(sys);
    if (!chem || !chem->num_bonds) return flags;
    size_t lo = 0, hi = chem->num_bonds;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (chem->bond_idx[mid] < idx) lo = mid + 1; else hi = mid;
    }
    return (lo < chem->num_bonds && chem->bond_idx[lo] == idx) ? chem->bond_flags[lo] : flags;
}

// The bond between the two atoms of a pair, as of the frame shown: of the system (the graph, the bonds outside of it in
// its range, the interactions), or formed in the frame
static PickedBondSource bond_between(const SystemContext& sc, md_atom_pair_t pair, md_bond_flags_t* flags) {
    const md_system_t& sys = sc.data.sys;
    *flags = MD_BOND_FLAG_NONE;
    const md_atom_idx_t a = pair.idx[0], b = pair.idx[1];
    if (a < 0 || b < 0 || (size_t)a >= sys.atom.count || (size_t)b >= sys.atom.count) return PickedBondSource::Unbonded;
    md_bond_idx_t idx = md_system_bond_find(&sys, a, b);
    if (idx == -1) idx = md_bond_find_interaction(&sys.bond, a, b);
    if (idx == -1 && sys.bond.pairs && sys.bond.flags) {
        // The elastic and metallic bonds have no connectivity: all of the range is looked through, as for the tooltip
        // of an atom of a metal
        for (size_t i = 0; i < sys.bond.count && idx == -1; ++i) {
            if (md_bond_is_graph(sys.bond.flags[i])) continue;
            const md_atom_pair_t p = sys.bond.pairs[i];
            if ((p.idx[0] == a && p.idx[1] == b) || (p.idx[0] == b && p.idx[1] == a)) idx = (md_bond_idx_t)i;
        }
    }
    if (idx != -1) {
        *flags = (size_t)idx < sys.bond.count ? frame_bond_flags(sc, (size_t)idx, sys.bond.flags[idx]) : md_system_bond_flags(&sys, idx);
        return PickedBondSource::System;
    }
    const md_bond_delta_t* delta = bond_tracking_delta(sc);
    for (size_t k = 0; delta && k < delta->num_formed; ++k) {
        const md_atom_pair_t p = delta->formed_pair[k];
        if ((p.idx[0] == a && p.idx[1] == b) || (p.idx[0] == b && p.idx[1] == a)) {
            *flags = frame_bond_flags(sc, sys.bond.count + k, delta->formed_flags[k]);
            return PickedBondSource::Formed;
        }
    }
    return PickedBondSource::Unbonded;
}

bool picked_bond(const ApplicationState& app, const PickingHit& hit, md_atom_pair_t* pair, md_bond_flags_t* flags, PickedBondSource* source) {
    if (hit.domain != PickingDomain_Bond) return false;
    // Of the system the hit names
    const SystemContext* sc = picking_hit_system(app, hit);
    if (!sc) return false;
    if (hit.key != 0) {
        // A weak bond of a representation: its atoms are what it has of its own, and a bond its filter found for the
        // frame where the system has none (a hydrogen, ionic or halogen bond bond() finds)
        if (!representation_weak_bond_pair(pair, *sc, (uint32_t)hit.key, hit.local_idx)) return false;
        *source = bond_between(*sc, *pair, flags);
        md_bond_flags_t found = MD_BOND_FLAG_NONE;
        if (*source == PickedBondSource::Unbonded && representation_weak_bond_flags(&found, *sc, (uint32_t)hit.key, hit.local_idx) && !(found & MD_BOND_FLAG_INACTIVE)) {
            *flags = found;
            *source = PickedBondSource::Found;
        }
        return true;
    }
    const md_system_t& sys = sc->data.sys;
    const size_t idx = hit.local_idx;
    if (idx < sys.bond.count) {
        *pair = sys.bond.pairs[idx];
        *flags = frame_bond_flags(*sc, idx, sys.bond.flags[idx]);
        *source = PickedBondSource::System;
        return true;
    }
    const md_bond_delta_t* delta = bond_tracking_delta(*sc);
    const size_t k = idx - sys.bond.count;
    if (delta && k < delta->num_formed) {
        *pair = delta->formed_pair[k];
        *flags = frame_bond_flags(*sc, idx, delta->formed_flags[k]);
        *source = PickedBondSource::Formed;
        return true;
    }
    return false;
}

bool picking_range_reserve(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count, uint64_t key, uint32_t system) {
    MD_ASSERT(space);

    if (count > 0 && space->num_ranges < MD_COUNTOF(space->ranges)) {
        PickingRange* curr_range = &space->ranges[space->num_ranges++];
        PickingRange* prev_range = space->num_ranges > 1 ? &space->ranges[space->num_ranges - 2] : NULL;
        curr_range->domain = domain;
        curr_range->beg = prev_range ? prev_range->end : 0;
        curr_range->end = curr_range->beg + (uint32_t)count;
        curr_range->key = key;
        curr_range->system = system;
        if (out_range) {
            md_memcpy(out_range, curr_range, sizeof(PickingRange));
        }
        return true;
    }

    return false;
}

void picking_handler_new_frame(PickingHandler* handler) {
    MD_ASSERT(handler);
    handler->frame_idx += 1;

    const uint32_t slot_idx = handler->frame_idx % MD_COUNTOF(handler->history);
    handler->history[slot_idx].submitted_frame_idx = handler->frame_idx;
    handler->history[slot_idx].space = PickingSpace{};
}

PickingSpace* picking_handler_current_space(PickingHandler* handler) {
    MD_ASSERT(handler);
    return &handler->history[handler->frame_idx % MD_COUNTOF(handler->history)].space;
}

const PickingSpace* picking_handler_find_space(const PickingHandler& handler, uint32_t submitted_frame_idx) {
    for (size_t i = 0; i < MD_COUNTOF(handler.history); ++i) {
        const auto& hist = handler.history[i];
        if (hist.submitted_frame_idx == submitted_frame_idx) {
            return &hist.space;
        }
    }

    return nullptr;
}

void picking_surface_init(PickingSurface* surface, PickingSourceID source) {
    MD_ASSERT(surface);

    *surface = PickingSurface{};
    surface->source = source;

    for (size_t i = 0; i < MD_COUNTOF(surface->slots); ++i) {
        auto& slot = surface->slots[i];
        glGenBuffers(1, &slot.color_pbo);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, nullptr, GL_DYNAMIC_READ);

        glGenBuffers(1, &slot.depth_pbo);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, nullptr, GL_DYNAMIC_READ);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
}

void picking_surface_free(PickingSurface* surface) {
    MD_ASSERT(surface);

    for (size_t i = 0; i < MD_COUNTOF(surface->slots); ++i) {
        auto& slot = surface->slots[i];
        if (slot.color_pbo) glDeleteBuffers(1, &slot.color_pbo);
        if (slot.depth_pbo) glDeleteBuffers(1, &slot.depth_pbo);
    }

    *surface = PickingSurface{};
}

bool picking_surface_submit_readback(
    PickingSurface* surface,
    uint32_t fbo,
    uint32_t width,
    uint32_t height,
    uint32_t submitted_frame_idx,
    md_float2_t surface_coord,
    md_float2_t screen_coord,
    const md_float4x4_t& clip_to_world
) {
    MD_ASSERT(surface);

    if (!fbo || width == 0 || height == 0) {
        return false;
    }

    const int x = (int)surface_coord.x;
    const int y = (int)surface_coord.y;
    if (x < 0 || y < 0 || x >= (int)width || y >= (int)height) {
        return false;
    }

    const uint32_t queue_idx = surface->slot_cursor % MD_COUNTOF(surface->slots);
    surface->slot_cursor += 1;

    auto& slot = surface->slots[queue_idx];
    if (slot.color_pbo == 0 || slot.depth_pbo == 0) {
        MD_LOG_ERROR("Invalid PBOs in picking surface slot");
        return false;
    }

    slot.submitted_frame_idx = submitted_frame_idx;
    slot.pending = true;
    slot.viewport_width = width;
    slot.viewport_height = height;
    slot.surface_coord = surface_coord;
    slot.screen_coord = screen_coord;
    slot.clip_to_world = clip_to_world;

    PUSH_GPU_SECTION("QUEUE PICKING READBACK")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
    glReadPixels(x, y, 1, 1, GL_BGRA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    POP_GPU_SECTION()

    return true;
}

const PickingRange* picking_space_find_range(const PickingSpace& space, PickingDomainID domain, uint64_t key, uint32_t system) {
    for (size_t i = 0; i < space.num_ranges; ++i) {
        if (space.ranges[i].domain == domain && space.ranges[i].key == key && space.ranges[i].system == system) {
            return &space.ranges[i];
        }
    }

    return nullptr;
}

static const PickingRange* find_picking_range(const PickingSpace* space, uint32_t raw_idx) {
    MD_ASSERT(space);

    for (size_t i = 0; i < space->num_ranges; ++i) {
        const PickingRange* range = &space->ranges[i];
        if (range->beg <= raw_idx && raw_idx < range->end) {
            return range;
        }
    }

    return nullptr;
}

bool picking_surface_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler
) {
    MD_ASSERT(out_hit);
    MD_ASSERT(surface);

    *out_hit = {};

    if (surface->slot_cursor < MD_COUNTOF(surface->slots)) {
        return false;
    }

    const uint32_t read_idx = surface->slot_cursor % MD_COUNTOF(surface->slots);
    auto& slot = surface->slots[read_idx];
    if (!slot.pending) {
        return false;
    }

    uint8_t color[4] = {};
    float depth = 1.0f;

    PUSH_GPU_SECTION("POLL PICKING READBACK")
    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
    glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(color), color);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
    glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(depth), &depth);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    POP_GPU_SECTION()

    slot.pending = false;

    const uint32_t raw_idx = (color[0] << 16) | (color[1] << 8) | (color[2] << 0) | (color[3] << 24);

    const PickingSpace* space = picking_handler_find_space(handler, slot.submitted_frame_idx);
    if (!space) {
        return false;
    }

    PickingRange range = {0};
    if (const PickingRange* space_range = find_picking_range(space, raw_idx)) {
        range = *space_range;
    }

    out_hit->source = surface->source;
    out_hit->domain = range.domain;
    out_hit->key = range.key;
    out_hit->system = range.system;
    out_hit->frame_idx = slot.submitted_frame_idx;
    out_hit->raw_idx = raw_idx;
    out_hit->local_idx = raw_idx - range.beg;
    out_hit->surface_coord = slot.surface_coord;
    out_hit->screen_coord = slot.screen_coord;
    out_hit->depth = depth;

    const md_float4_t viewport = {0, 0, (float)slot.viewport_width, (float)slot.viewport_height};
    out_hit->world_pos = md_float4x4_unproject({slot.surface_coord.x, slot.surface_coord.y, depth}, slot.clip_to_world, viewport);

    return true;
}

bool picking_surface_submit_readback_and_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler,
    const PickingReadbackRequest& request
) {
    MD_ASSERT(out_hit);
    MD_ASSERT(surface);

    picking_surface_submit_readback(
        surface,
        request.fbo,
        request.width,
        request.height,
        handler.frame_idx,
        request.surface_coord,
        request.screen_coord,
        request.clip_to_world
    );

    return picking_surface_poll_hit(out_hit, surface, handler);
}

InteractionSurfaceState interaction_surface(InteractionSurfaceID id, const md_float2_t& size, InteractionSurfaceFlags flags) {
    InteractionSurfaceState state = {};

    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (!window) {
        MD_LOG_ERROR("No current ImGui window for interaction surface");
        return state;
    }

    state.surface_id = id;

    static const ImGuiButtonFlags btn_flags = ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap;
    ImGui::InvisibleButton("interaction surface button", vec_cast(size), btn_flags);
    state.item_id = ImGui::GetItemID();

    state.hovered = ImGui::IsItemHovered();
    state.active = ImGui::IsItemActive();
    state.activated = ImGui::IsItemActivated();
    state.deactivated = ImGui::IsItemDeactivated();

    ImDrawList* draw_list = window->DrawList;
    MD_ASSERT(draw_list);

    const ImVec2 canvas_min = ImGui::GetItemRectMin();
    const ImVec2 canvas_max = ImGui::GetItemRectMax();
    const ImVec2 canvas_size = ImGui::GetItemRectSize();

    state.surface_size = vec_cast(canvas_size);

    const ImVec2 mouse_pos = ImGui::GetMousePos();
    const ImVec2 local_mouse = mouse_pos - canvas_min;
    state.mouse_local = vec_cast(local_mouse);

    if (state.activated || state.active || state.deactivated) {
        if (ImGui::IsKeyPressed(ImGuiMod_Shift, false)) {
            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
        }

        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
            ImGuiMouseButton button = ImGuiMouseButton_Left;
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                state.selection_mode = InteractionSelectionMode::Append;
                button = ImGuiMouseButton_Left;
            } else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                state.selection_mode = InteractionSelectionMode::Remove;
                button = ImGuiMouseButton_Right;
            }

            if (state.selection_mode != InteractionSelectionMode::None) {
                const ImVec2 ext = ImGui::GetMouseDragDelta(button);
                const ImVec2 pos = ImGui::GetMousePos() - ext;

                ImVec2 sel_min = ImClamp(ImMin(pos, pos + ext), canvas_min, canvas_max);
                ImVec2 sel_max = ImClamp(ImMax(pos, pos + ext), canvas_min, canvas_max);

                draw_list->AddRectFilled(sel_min, sel_max, 0x22222222);
                draw_list->AddRect(sel_min, sel_max, 0x88888888);

                sel_min -= canvas_min;
                sel_max -= canvas_min;

                state.region_min = vec_cast(sel_min);
                state.region_max = vec_cast(sel_max);
            }
        }
    }

    if (flags & InteractionSurfaceFlags_NoRegionSelect) {
        state.region_min = { 0, 0 };
        state.region_max = { 0, 0 };
    }

    return state;
}

bool interaction_surface_hit_extract(PickingHit* out_hit, const InteractionSurfaceState& state, const InteractionSurfaceHitArgs& args) {
    MD_ASSERT(out_hit);
    if (state.hovered && args.fbo && args.width && args.height) {
        const ImVec2 local_coord = ImVec2(state.mouse_local.x, state.surface_size.y - state.mouse_local.y) * ImGui::GetIO().DisplayFramebufferScale * args.pixel_scale;

        PickingReadbackRequest request = {
            .fbo = args.fbo,
            .width = args.width,
            .height = args.height,
            .surface_coord = {local_coord.x, local_coord.y},
            .screen_coord = {state.mouse_local.x, state.mouse_local.y},
            .clip_to_world = args.clip_to_world,
        };

        return picking_surface_submit_readback_and_poll_hit(out_hit, args.picking_surface, args.picking_handler, request);
    }
    return false;
}

InteractionSurfaceViewTransformResult interaction_surface_view_transform_apply(ViewTransform* target, const InteractionSurfaceState& state, const InteractionSurfaceViewTransformArgs& args) {
    MD_ASSERT(target);
    InteractionSurfaceViewTransformResult result = {};
    if (state.active || state.hovered) {
        if (state.selection_mode == InteractionSelectionMode::None) {
            const md_float2_t delta = vec_cast(ImGui::GetIO().MouseDelta);
            const md_float2_t coord = state.mouse_local;
            const float scroll_delta = ImGui::GetIO().MouseWheel;

            TrackballControllerInput input = {};
            input.rotate_button = ImGui::IsMouseDown(ImGuiMouseButton_Left);
            input.pan_button = ImGui::IsMouseDown(ImGuiMouseButton_Right);
            input.dolly_button = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
            input.mouse_coord_curr = coord;
            input.mouse_coord_prev = coord - delta;
            input.screen_size = state.surface_size;
            input.dolly_delta = scroll_delta;
            input.fov_y = args.camera.fov_y;
            input.ndc_center = args.view_shift;

            TrackballFlags flags = TrackballFlags_None;
            if (state.active) {
                flags |= TrackballFlags_EnableAllInteractions;
            } else {
                flags |= TrackballFlags_DollyEnabled;
            }

            camera_controller_trackball(target, input, args.trackball_param, flags);

            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                result.reset_requested = true;
            }
        }
    }
    return result;
}

void interaction_surface_event_extract(InteractionSurfaceEvent* event, const InteractionSurfaceState& state, const PickingHit& hit) {
    MD_ASSERT(event);
    *event = InteractionSurfaceEvent{};

    event->surface_id   = state.surface_id;
    event->item_id      = state.item_id;
    event->mouse_local  = state.mouse_local;
    event->surface_size = state.surface_size;
    event->region_min   = state.region_min;
    event->region_max   = state.region_max;
    event->hit = hit;

    if (ImGui::IsKeyDown(ImGuiMod_Shift) && state.region_max != state.region_min) {
        event->selection_mode = state.selection_mode;
        event->kind = InteractionSurfaceEventKind::RegionSelect;
        if (state.active) {
            event->region_phase = InteractionSurfaceEventPhase::Update;
        } else if (state.deactivated) {
            event->region_phase = InteractionSurfaceEventPhase::Commit;
        }
    } else if (state.hovered) {
		bool left_click  = ImGui::IsMouseReleased(ImGuiMouseButton_Left)  && ImGui::GetMouseDragDelta(ImGuiMouseButton_Left)  == ImVec2(0, 0);
        bool right_click = ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::GetMouseDragDelta(ImGuiMouseButton_Right) == ImVec2(0, 0);
        if (right_click && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
            event->kind = InteractionSurfaceEventKind::ContextMenu;
		} else if (left_click || right_click) {
			event->kind = InteractionSurfaceEventKind::Click;
            event->selection_mode = state.selection_mode;
        } else if (!ImGui::IsMouseDown(ImGuiMouseButton_Left) && !ImGui::IsMouseDown(ImGuiMouseButton_Right) && !ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
            // No hover while a button is held: the view is being rotated, panned or zoomed
            event->kind = InteractionSurfaceEventKind::Hover;
        }
    }
}

void point_set_region_mask_compute(md_bitmask_t* mask,
    const md_float3_t xyz[],
    size_t count,
    const md_bitset_t* candidate_mask,
    const md_float4x4_t& world_to_clip,
    const md_float2_t& region_min,
    const md_float2_t& region_max,
    const md_float2_t& surface_size)
{
    MD_ASSERT(mask);
    MD_ASSERT(xyz);
    // The candidates are read while the mask is written: they must not be a view of it
    MD_ASSERT(!candidate_mask || !candidate_mask->words || !mask->words || candidate_mask->words < mask->words || candidate_mask->words >= mask->words + mask->cap);

    md_bitmask_clear(mask);

    // Transform visible atoms using supplied world_to_clip and set if within region

    if (candidate_mask) {
        // Do for candidate set
        md_bitset_iter_t it = md_bitset_iter(*candidate_mask);

        while (md_bitset_iter_next(&it)) {
            size_t idx = it.idx;
            md_float4_t xyz1 = md_float4_from_float3(xyz[idx], 1.0f);
            md_float4_t coord = md_float4x4_mul_float4(world_to_clip, xyz1);
            md_float2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitmask_set_bit(mask, (uint32_t)idx);
            }
        }
    } else {
        // Do for full set
        for (size_t i = 0; i < count; ++i) {
            md_float4_t xyz1 = md_float4_from_float3(xyz[i], 1.0f);
            md_float4_t coord = md_float4x4_mul_float4(world_to_clip, xyz1);
            md_float2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitmask_set_bit(mask, (uint32_t)i);
            }
        }
    }
}

bool file_queue_empty(const FileQueue* queue) {
    return queue->head == queue->tail;
}

bool file_queue_full(const FileQueue* queue) {
    return (queue->head + 1) % MD_COUNTOF(queue->arr) == queue->tail;
}

void file_queue_push(FileQueue* queue, md_str_t path, FileFlags flags, uint32_t system) {
    MD_ASSERT(queue);
    MD_ASSERT(!file_queue_full(queue));
    int prio = 5;


    md_str_t ext = md_path_ext(path);
    loader::LoaderState detected = {};
    if (!md_str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
        // By the contents where the extension does not tell (Amber's mdout, prmtop and inpcrd have none)
        loader::init(&detected, path);
    }
    if (md_str_eq(ext, WORKSPACE_FILE_EXTENSION) || detected.type != LoaderType_Undefined) {
        LoaderFlags loader_flags = detected.flags;
        if (md_str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
            prio = 1;
        } else if (loader_flags & LoaderFlag_System) {
            prio = 2;
        } else if (loader_flags & LoaderFlag_Temporal) {
            // Joins the trajectory's run, so it goes after the trajectory when both are dropped
            prio = 4;
        } else if (loader_flags & (LoaderFlag_Trajectory | LoaderFlag_Supplemental)) {
            // A supplemental file needs the system in place, same as a trajectory
            prio = 3;
        } else {
            flags |= FileFlags_ShowDialogue;
        }
    } else {
        // Unknown extension
        flags |= FileFlags_ShowDialogue;
    }

    uint32_t i = queue->head;
    queue->arr[queue->head] = {md_alloc_str(queue->ring, path), flags, prio, system};
    queue->head = (queue->head + 1) % MD_COUNTOF(queue->arr);

    // Sort queue based on prio
     while (i != queue->tail && queue->arr[i].prio < queue->arr[(i - 1) % MD_COUNTOF(queue->arr)].prio) {
        FileQueue::Entry tmp = queue->arr[i];
        queue->arr[i] = queue->arr[(i - 1) % MD_COUNTOF(queue->arr)];
        queue->arr[(i - 1) % MD_COUNTOF(queue->arr)] = tmp;
        i = (i - 1) % MD_COUNTOF(queue->arr);
     }
}

FileQueue::Entry file_queue_front(const FileQueue* queue) {
    MD_ASSERT(!file_queue_empty(queue));
    return queue->arr[queue->tail];
}

FileQueue::Entry file_queue_pop(FileQueue* queue) {
    MD_ASSERT(queue);
    MD_ASSERT(!file_queue_empty(queue));
    FileQueue::Entry front = file_queue_front(queue);
    queue->tail = (queue->tail + 1) % MD_COUNTOF(queue->arr);
    return front;
}

// What follows a file loaded into sys: a system built shown with the default representations (a system loaded again
// keeps its own if asked to), and every system framed when a system or a run is new. A system new for it which nothing
// went into is taken out again.
static void load_into_done(ApplicationState* app, SystemContext* sys, bool fresh, LoadResult result) {
    if (!system_loaded(*sys)) {
        if (fresh) system_release(app, sys);
        return;
    }
    if (fresh) {
        // A name of its own, which the systems are told apart by (a copy of one is named after the same file)
        const size_t len = strnlen(sys->name, sizeof(sys->name));
        for (int n = 2; n < 100; ++n) {
            bool taken = false;
            for (size_t i = 0; i < app->num_systems; ++i) {
                taken |= app->systems[i] != sys && system_loaded(*app->systems[i]) && strcmp(app->systems[i]->name, sys->name) == 0;
            }
            if (!taken) break;
            snprintf(sys->name + MD_MIN(len, sizeof(sys->name) - 8), 8, " (%d)", n);
        }
    }
    if (result == LoadResult_System) {
        if (fresh || !app->settings.keep_representations) {
            remove_all_representations(*sys);
            if (system_inference_pending(*sys)) {
                // The defaults are decided from what inference tells of the system (what its residues are, its
                // strands, its lipids): its atoms meanwhile
                Representation* rep = create_representation(*sys, RepresentationType::SpaceFill, ColorMapping::Type);
                snprintf(rep->name, sizeof(rep->name), "default");
                md_str_copy_to_char_buf(sys->data.inference.placeholder_name, sizeof(sys->data.inference.placeholder_name), md_str_from_cstr(rep->name));
                sys->data.inference.defaults_pending = true;
            } else {
                create_default_representations(*sys);
            }
        }
        recompute_atom_visibility_mask(*sys);
        system_instances_update(*sys, MolBit_ClearVelocity);
    }
    if (result == LoadResult_System || result == LoadResult_Trajectory) {
        clock_update(app);
        clock_set_frame(app, 0.0);
        app->animation.mode = PlaybackMode::Stopped;
        reset_view_scene(&app->view.target, *app, app->view.camera.fov_y, app->view.region);
    }
}

LoadResult load_file_into(ApplicationState* app, SystemContext* target, md_str_t path, const loader::LoaderState& load_state) {
    MD_ASSERT(app);
    SystemContext* sys = target ? target : system_add(app);
    if (!sys) {
        VIAMD_LOG_ERROR("'" MD_STR_FMT "' is not loaded: at most %zu systems are loaded at once", MD_STR_ARG(path), SYSTEM_MAX);
        return LoadResult_Failed;
    }
    const bool fresh = !system_loaded(*sys);
    const bool own_run = (load_state.flags & LoaderFlag_System) && (load_state.flags & LoaderFlag_Trajectory);
    const LoadResult result = load_data(app, *sys, path, load_state, FileRole_None, own_run);
    load_into_done(app, sys, fresh, result);
    return result;
}

// The files of a system but those of its runs, loaded again into another: a copy of it (PlacementAction_Copy). What it
// made of them: LoadResult_System, or LoadResult_Failed when nothing loaded.
static LoadResult load_system_copy(ApplicationState* app, const SystemContext& of, SystemContext& into) {
    bool any = false;
    for (size_t i = 0; i < file_list_count(&of.files); ++i) {
        const FileEntry e = *file_list_get(&of.files, i);
        if (e.run >= 0 || e.missing) continue;
        const loader::LoaderState ls = file_entry_loader_state(e);
        any |= load_data(app, into, md_str_from_cstr(e.path), ls, file_list_role(&of.files, i)) != LoadResult_Failed;
    }
    return any ? LoadResult_System : LoadResult_Failed;
}

// The systems loaded, as a drop sees them: their runs' paths in run_paths
static size_t drop_systems(DropSystem* out, md_str_t (*run_paths)[RUN_MAX], const ApplicationState& app) {
    size_t n = 0;
    for (size_t i = 0; i < app.num_systems; ++i) {
        const SystemContext& sc = *app.systems[i];
        if (!system_loaded(sc)) continue;
        for (size_t r = 0; r < sc.num_runs; ++r) run_paths[n][r] = md_str_from_cstr(sc.runs[r]->path);
        out[n] = { sc.id, &sc.data.sys, &sc.data.state, run_paths[n], sc.num_runs };
        n += 1;
    }
    return n;
}

// The files a system of the plan loads, into sys (fresh: added for them, and shown what was loaded into it before them,
// a copy), in its order
static void load_placed(ApplicationState* app, DropFiles* drop, const PlacementSystem& ps, SystemContext* sys, bool fresh, LoadResult shown = LoadResult_Failed) {
    for (size_t i = 0; i < ps.num_entries; ++i) {
        DropFile& f = drop->files[ps.entries[i].file];
        const md_str_t path = md_str_from_cstr(f.entry.path);
        if (file_list_has(&sys->files, path)) {
            VIAMD_LOG_INFO("'%s' is loaded into '%s' already", f.entry.path, sys->name);
            continue;
        }
        // The system it is, loaded to be matched, is taken rather than read again
        const DropFile* loaded = f.loaded > 0 ? &f : nullptr;
        const FileRole role = f.kind == PlacementKind_Series ? FileRole_Series : FileRole_None;
        const bool own_run = (f.ls.flags & LoaderFlag_System) && (f.ls.flags & LoaderFlag_Trajectory);
        // A series along the run it was matched along: one loaded, or the one a trajectory of the drop opened
        const PlacementRef along = ps.entries[i].run;
        int run = -1;
        if (along.file >= 0 && (size_t)along.file < drop->num_files) {
            run = file_list_run_of(&sys->files, md_str_from_cstr(drop->files[along.file].entry.path));
        } else if (along.run >= 0 && along.system == sys->id) {
            run = along.run;
        }
        const LoadResult r = load_data(app, *sys, path, f.ls, role, own_run, loaded ? &loaded->sys : nullptr, loaded ? &loaded->state : nullptr, run);
        if (r != LoadResult_Failed) f.into = sys->id;
        if (r == LoadResult_System || (r == LoadResult_Trajectory && shown != LoadResult_System)) shown = r;
    }
    load_into_done(app, sys, fresh, shown);
}

// The files of one drop where the placement puts them (drop_files.h, file_placement.h). True when some are of the atoms
// of more than one system: they are left for the user to pick where (the plan's systems of PlacementAction_Ask).
static bool drop_files_load(ApplicationState* app, DropFiles* drop, const PlacementPlan& plan) {
    for (size_t i = 0; i < plan.num_refused; ++i) {
        const PlacementRefusal& r = plan.refused[i];
        const char* path = drop->files[r.file].entry.path;
        switch (r.reason) {
        case PlacementReason_NeedsSystem: VIAMD_LOG_ERROR("'%s' adds to a system's atoms, and is of none of those loaded", path); break;
        case PlacementReason_NeedsRun:    VIAMD_LOG_ERROR("'%s' is along a trajectory, and along none of those loaded", path); break;
        case PlacementReason_NotOfTarget: VIAMD_LOG_ERROR("'%s' is not of the atoms of the system it was put on", path); break;
        case PlacementReason_Cap:         VIAMD_LOG_ERROR("'%s' is not loaded: at most %zu systems are loaded at once", path, SYSTEM_MAX); break;
        case PlacementReason_TooMany:     VIAMD_LOG_ERROR("'%s' is not loaded: at most %zu files at once", path, PLACEMENT_MAX_FILES); break;
        }
    }

    bool asks = false;
    for (size_t k = 0; k < plan.num_systems; ++k) {
        const PlacementSystem& ps = plan.systems[k];
        if (ps.action == PlacementAction_Ask) {
            asks = true;
            continue;
        }

        SystemContext* sys = nullptr;
        bool fresh = false;
        LoadResult shown = LoadResult_Failed;
        if (ps.action == PlacementAction_Join) {
            sys = system_find(app, ps.system);
        } else {
            sys = system_add(app);
            fresh = sys != nullptr;
            if (sys && ps.action == PlacementAction_Copy) {
                // A copy of the system the files are of, which they would replace something of
                if (const SystemContext* of = system_find(app, ps.system)) {
                    VIAMD_LOG_INFO("The files dropped are of the atoms of '%s' and would replace some of it: they go into a copy of it", of->name);
                    shown = load_system_copy(app, *of, *sys);
                }
            }
        }
        if (!sys) {
            VIAMD_LOG_ERROR("At most %zu systems are loaded at once: %zu files are not loaded", SYSTEM_MAX, ps.num_entries);
            continue;
        }
        load_placed(app, drop, ps, sys, fresh, shown);
    }
    return asks;
}

// ### WHERE FILES GO WHICH ARE OF MORE THAN ONE SYSTEM ###
//
// A question per system of the plan the user is asked about: its files, and where they may go - each of its candidates,
// and a new system where one is offered. Answered (or not, with Escape), file_queue_process loads them as told, in the
// order of the plan: a candidate which is a system the drop makes is where its files went.

// The system a candidate is: a loaded one, or the one the files of the drop it names went into. NULL when it is neither.
static SystemContext* placement_ref_system(ApplicationState* app, const DropFiles* drop, PlacementRef ref) {
    if (ref.file >= 0) {
        return (size_t)ref.file < drop->num_files ? system_find(app, drop->files[ref.file].into) : nullptr;
    }
    return system_find(app, ref.system);
}

static void placement_prompt_close(ApplicationState* app) {
    auto& p = app->placement_prompt;
    if (p.drop) {
        drop_files_free(p.drop);
        delete p.drop;
    }
    delete p.plan;
    p.drop = nullptr;
    p.plan = nullptr;
    p.answered = false;
}

// The drop's files of several systems, asked about: the first candidate of each picked until the user says otherwise
static void placement_prompt_open(ApplicationState* app, DropFiles* drop, PlacementPlan* plan) {
    auto& p = app->placement_prompt;
    MD_ASSERT(!p.drop);
    p.drop = drop;
    p.plan = plan;
    p.answered = false;
    for (size_t k = 0; k < MD_COUNTOF(p.choice); ++k) {
        p.choice[k] = -1;
        if (k >= plan->num_systems) continue;
        const PlacementSystem& ps = plan->systems[k];
        for (size_t j = 0; j < ps.num_candidates && p.choice[k] < 0; ++j) {
            if (placement_ref_system(app, drop, ps.candidates[j])) p.choice[k] = (int)j;
        }
        if (p.choice[k] < 0 && ps.new_allowed) p.choice[k] = (int)ps.num_candidates;
    }
}

// The answers loaded, in the order of the plan
static void placement_prompt_load(ApplicationState* app) {
    auto& p = app->placement_prompt;
    for (size_t k = 0; k < p.plan->num_systems; ++k) {
        const PlacementSystem& ps = p.plan->systems[k];
        if (ps.action != PlacementAction_Ask) continue;
        const int c = p.choice[k];
        if (c < 0) {
            for (size_t i = 0; i < ps.num_entries; ++i) VIAMD_LOG_INFO("'%s' is not loaded", p.drop->files[ps.entries[i].file].entry.path);
            continue;
        }
        SystemContext* sys = nullptr;
        bool fresh = false;
        if ((size_t)c >= ps.num_candidates) {
            sys = system_add(app);
            fresh = sys != nullptr;
        } else {
            sys = placement_ref_system(app, p.drop, ps.candidates[c]);
        }
        if (!sys) {
            VIAMD_LOG_ERROR("The system picked for '%s' is not there any more: %zu files are not loaded", p.drop->files[ps.entries[0].file].entry.path, ps.num_entries);
            continue;
        }
        load_placed(app, p.drop, ps, sys, fresh);
    }
    placement_prompt_close(app);
}

// What a candidate is called: the system, with the run of it a series goes along
static void placement_ref_label(char* buf, size_t cap, ApplicationState* app, const DropFiles* drop, PlacementRef ref) {
    const SystemContext* sys = placement_ref_system(app, drop, ref);
    if (!sys) {
        const md_str_t name = ref.file >= 0 && (size_t)ref.file < drop->num_files ? md_path_filename(md_str_from_cstr(drop->files[ref.file].entry.path)) : MD_STR_LIT("?");
        snprintf(buf, cap, MD_STR_FMT " (not loaded)", MD_STR_ARG(name));
        return;
    }
    char label[128];
    system_label(label, sizeof(label), *sys);
    // A run: the file that opened it
    const char* run_file = nullptr;
    if (ref.run >= 0) {
        const FileEntry* e = file_list_run_file_of(&sys->files, (size_t)ref.run);
        run_file = e ? e->path : nullptr;
    } else if (ref.file >= 0 && drop->files[ref.file].kind == PlacementKind_Trajectory) {
        run_file = drop->files[ref.file].entry.path;
    }
    if (run_file) {
        snprintf(buf, cap, "%s, along " MD_STR_FMT, label, MD_STR_ARG(md_path_filename(md_str_from_cstr(run_file))));
    } else {
        snprintf(buf, cap, "%s", label);
    }
}

void draw_placement_prompt(ApplicationState* state) {
    auto& p = state->placement_prompt;
    if (!p.drop || p.answered) return;

    const char* title = "Where do these files go?##placement";
    if (!ImGui::IsPopupOpen(title)) ImGui::OpenPopup(title);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal(title, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;

    ImGui::TextDisabled("They are of the atoms of more than one system.");
    ImGui::Spacing();
    if (ImGui::BeginTable("##questions", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg | ImGuiTableFlags_PadOuterX)) {
        for (size_t k = 0; k < p.plan->num_systems; ++k) {
            const PlacementSystem& ps = p.plan->systems[k];
            if (ps.action != PlacementAction_Ask) continue;
            ImGui::PushID((int)k);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::AlignTextToFramePadding();
            for (size_t i = 0; i < ps.num_entries; ++i) {
                const DropFile& f = p.drop->files[ps.entries[i].file];
                const md_str_t name = md_path_filename(md_str_from_cstr(f.entry.path));
                if (i > 0) ImGui::SameLine(0, 0);
                ImGui::Text("%s" MD_STR_FMT, i > 0 ? ", " : "", MD_STR_ARG(name));
                ImGui::SetItemTooltip("%s", f.entry.path);
            }

            ImGui::TableNextColumn();
            char preview[256] = "Not loaded";
            int& c = p.choice[k];
            if (c >= 0 && (size_t)c < ps.num_candidates) placement_ref_label(preview, sizeof(preview), state, p.drop, ps.candidates[c]);
            else if (c >= 0) snprintf(preview, sizeof(preview), "A new system");
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 18.0f);
            if (ImGui::BeginCombo("##into", preview)) {
                for (size_t j = 0; j < ps.num_candidates; ++j) {
                    char label[256];
                    placement_ref_label(label, sizeof(label), state, p.drop, ps.candidates[j]);
                    ImGui::BeginDisabled(!placement_ref_system(state, p.drop, ps.candidates[j]));
                    ImGui::PushID((int)j);
                    if (ImGui::Selectable(label, c == (int)j)) c = (int)j;
                    ImGui::PopID();
                    ImGui::EndDisabled();
                }
                if (ps.new_allowed && ImGui::Selectable("A new system", c == (int)ps.num_candidates)) c = (int)ps.num_candidates;
                if (ImGui::Selectable("Not loaded", c < 0)) c = -1;
                ImGui::EndCombo();
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    ImGui::Spacing();
    if (ImGui::Button("Load")) {
        p.answered = true;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Load none") || ImGui::IsKeyPressed(ImGuiKey_Escape)) {
        for (size_t k = 0; k < MD_COUNTOF(p.choice); ++k) p.choice[k] = -1;
        p.answered = true;
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}

// ### A FILE MOVED BETWEEN SYSTEMS ###

// The system built again of the files of list, in their places: what is left of it as a file goes. What is the user's of
// it stays: its name, and when it is of as many atoms again, its representations and the edits of its bonds.
static void system_rebuild(ApplicationState* app, SystemContext& sys, const FileList& list) {
    char name[sizeof(sys.name)];
    md_memcpy(name, sys.name, sizeof(name));
    const size_t num_atoms = sys.data.sys.atom.count;
    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
    md_bond_list_t edits = {};
    for (size_t i = 0; i < md_bond_list_count(&sys.data.sys.bond_edit); ++i) {
        const md_atom_pair_t pair = sys.data.sys.bond_edit.pairs[i];
        md_bond_list_push(&edits, pair.idx[0], pair.idx[1], sys.data.sys.bond_edit.flags[i], temp);
    }

    // The run each instance shows, by its path: what it is called among the runs opened again, which are those of list
    char inst_runs[INSTANCE_MAX][sizeof(SystemRun::path)] = {};
    int64_t inst_frames[INSTANCE_MAX] = {};
    for (size_t k = 0; k < sys.num_instances; ++k) {
        const SystemRun* r = instance_run(*sys.instances[k]);
        if (r) md_memcpy(inst_runs[k], r->path, sizeof(r->path));
        inst_frames[k] = sys.instances[k]->frame;
    }

    interrupt_async_tasks(app);
    free_trajectory_data(app, sys);
    free_system_data(app, sys);

    sys.instances_held = true;
    LoadResult shown = LoadResult_Failed;
    for (size_t i = 0; i < file_list_count(&list); ++i) {
        FileEntry e = *file_list_get(&list, i);
        const FileRole role = file_list_role(&list, i);
        LoadResult r = LoadResult_Failed;
        if (!e.missing) {
            r = load_data(app, sys, md_str_from_cstr(e.path), file_entry_loader_state(e), role);
        }
        if (r == LoadResult_Failed) {
            // Kept in its place, as a workspace keeps a file it could not load
            e.missing = true;
            file_list_record(&sys.files, e, role, app->allocator.persistent);
        }
        if (r == LoadResult_System || (r == LoadResult_Trajectory && shown != LoadResult_System)) shown = r;
    }
    md_memcpy(sys.name, name, sizeof(name));
    // The runs' offsets along the clock: of list, whose runs are recorded again in its order, a file which did not load
    // in its place
    for (size_t r = 0; r < md_array_size(list.runs) && r < md_array_size(sys.files.runs); ++r) {
        sys.files.runs[r].time_offset = list.runs[r].time_offset;
    }
    sys.instances_held = false;
    for (size_t k = sys.num_instances; k-- > 0;) {
        int run = -1;
        for (size_t i = 0; i < sys.num_runs && inst_runs[k][0]; ++i) {
            if (strcmp(sys.runs[i]->path, inst_runs[k]) == 0) run = (int)i;
        }
        // One whose run went goes with it, as when the run is closed (close_run); the first shows the first run left
        if (run < 0 && inst_runs[k][0]) {
            if (k > 0) {
                instance_remove(sys, k);
                continue;
            }
            if (sys.num_runs > 0) run = 0;
        }
        instance_show(*sys.instances[k], run, run >= 0 ? inst_frames[k] : -1);
    }

    if (num_atoms > 0 && sys.data.sys.atom.count == num_atoms) {
        apply_bond_edits(sys, edits);
        recompute_atom_visibility_mask(sys);
        system_instances_update(sys, MolBit_ClearVelocity);
    } else {
        load_into_done(app, &sys, false, shown);
    }
    md_arena_allocator_destroy(temp);
}

bool system_move_file(ApplicationState* app, uint32_t from_id, size_t i, uint32_t to_id) {
    MD_ASSERT(app);
    SystemContext* from = system_find(app, from_id);
    SystemContext* to = system_find(app, to_id);
    if (!from || i >= file_list_count(&from->files) || to == from || (to_id != SYSTEM_ID_NONE && !to)) return false;
    const FileEntry file = *file_list_get(&from->files, i);
    if (file.missing) {
        VIAMD_LOG_ERROR("'%s' is not loaded: there is nothing of it to move", file.path);
        return false;
    }

    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
    FileList rest = {}, moved = {};
    file_list_split(&rest, &moved, &from->files, i, temp);

    // The file where it fits, as a file dropped onto the system (or below the systems, for a new one) is
    DropSystem loaded[SYSTEM_MAX];
    md_str_t run_paths[SYSTEM_MAX][RUN_MAX];
    const size_t num_loaded = drop_systems(loaded, run_paths, *app);
    DropFiles* drop = new DropFiles();
    drop_files_init(drop, loaded, num_loaded);
    drop_files_add(drop, md_str_from_cstr(file.path), file_entry_loader_state(file), to ? to->id : PLACEMENT_NEW);
    PlacementPlan* plan = new PlacementPlan();
    drop_files_plan(plan, drop, SYSTEM_MAX);
    drop->loaded = nullptr;
    drop->num_loaded = 0;

    bool ok = plan->num_refused == 0;
    if (ok) {
        // One file onto a system named, or into a new one: nothing to ask about
        drop_files_load(app, drop, *plan);
        SystemContext* into = system_find(app, drop->files[0].into);
        ok = into != nullptr;
        if (into) {
            // What goes along with it: the series of a run it opened (the run opened again with it)
            for (size_t j = 0; j < file_list_count(&moved); ++j) {
                const FileEntry& e = *file_list_get(&moved, j);
                if (strcmp(e.path, file.path) == 0 || e.missing) continue;
                if (load_data(app, *into, md_str_from_cstr(e.path), file_entry_loader_state(e), file_list_role(&moved, j)) == LoadResult_Failed) {
                    VIAMD_LOG_ERROR("'%s' went with '%s' and was not loaded into '%s'", e.path, file.path, into->name);
                }
            }
            VIAMD_LOG_SUCCESS("Moved '%s' from '%s' to '%s'", file.path, from->name, into->name);
        }
    } else {
        drop_files_load(app, drop, *plan);   // Says why
    }
    drop_files_free(drop);
    delete drop;
    delete plan;

    // Runs of it, and nothing else of it (a trajectory, the series along it): closed, rather than the system built again
    size_t runs_moved[RUN_MAX];
    size_t num_runs_moved = 0;
    bool only_runs = true;
    for (size_t j = 0; j < file_list_count(&from->files); ++j) {
        const FileEntry& e = *file_list_get(&from->files, j);
        if (strcmp(e.path, file.path) != 0) continue;
        if (e.run < 0 || file_list_role(&from->files, j) != FileRole_Trajectory || num_runs_moved == RUN_MAX) only_runs = false;
        else runs_moved[num_runs_moved++] = (size_t)e.run;
    }
    if (ok && only_runs) {
        for (size_t k = num_runs_moved; k-- > 0;) close_run(app, *from, runs_moved[k]);
        // A trajectory which was the system (points): nothing of it is left
        if (file_list_count(&from->files) == 0) {
            system_remove(*from);
            system_release(app, from);
        }
    } else if (ok) {
        // What is left is a system while a file of it stands alone; what only adds to one goes with it
        bool stands = false;
        for (size_t j = 0; j < file_list_count(&rest); ++j) {
            const PlacementKind kind = drop_file_kind(loader::type_flags(file_list_get(&rest, j)->type));
            stands |= kind == PlacementKind_System || kind == PlacementKind_Topology || kind == PlacementKind_Trajectory;
        }
        if (stands) {
            system_rebuild(app, *from, rest);
        } else {
            for (size_t j = 0; j < file_list_count(&rest); ++j) {
                VIAMD_LOG_INFO("'%s' is unloaded: it adds to the atoms of '%s', of which nothing is left", file_list_get(&rest, j)->path, from->name);
            }
            system_remove(*from);
            system_release(app, from);
        }
    }
    md_arena_allocator_destroy(temp);
    return ok;
}

void file_queue_process(ApplicationState* state) {
    MD_ASSERT(state);
    if (state->load_dataset.show_window) return;

    // Files of several systems are asked about: what is queued waits for the answer
    if (state->placement_prompt.drop) {
        if (!state->placement_prompt.answered) return;
        placement_prompt_load(state);
        return;
    }
    if (state->run_close.system != SYSTEM_ID_NONE) {
        const auto close = state->run_close;
        state->run_close = {};
        if (SystemContext* sys = system_find(state, close.system)) close_run(state, *sys, close.run);
        return;
    }
    if (state->file_move.from != SYSTEM_ID_NONE) {
        const auto move = state->file_move;
        state->file_move = {};
        system_move_file(state, move.from, move.file, move.to);
        return;
    }
    if (file_queue_empty(&state->file_queue)) return;

    // Files dropped on the window: where, once the windows of the frame they were dropped in are drawn
    FileQueue& queue = state->file_queue;
    if ((file_queue_front(&queue).flags & FileFlags_Dropped) && ImGui::GetFrameCount() <= queue.drop_frame) return;

    // A workspace, and a file to ask about (the load dialog), on their own
    {
        const FileQueue::Entry e = file_queue_front(&state->file_queue);
        if (md_str_eq_ignore_case(md_path_ext(e.path), WORKSPACE_FILE_EXTENSION)) {
            file_queue_pop(&state->file_queue);
            load_workspace(state, e.path);
            return;
        }
        loader::LoaderState ls = {};
        loader::init(&ls, e.path);
        if ((e.flags & FileFlags_ShowDialogue) || (ls.flags & LoaderFlag_RequiresDialogue)) {
            file_queue_pop(&state->file_queue);
            state->load_dataset = LoadDatasetWindowState();
            md_str_copy_to_char_buf(state->load_dataset.path_buf, sizeof(state->load_dataset.path_buf), e.path);
            state->load_dataset.path_changed = true;
            state->load_dataset.show_window = true;
            const uint32_t system = (e.flags & FileFlags_Dropped) && queue.drop_system != SYSTEM_ID_NONE ? queue.drop_system : e.system;
            if (system != SYSTEM_ID_NONE) {
                state->load_dataset.system = system;
            } else if ((e.flags & FileFlags_NewSystem) || ((e.flags & FileFlags_Dropped) && queue.drop_new) || system_count(*state) == 0) {
                state->load_dataset.new_system = true;
            } else {
                state->load_dataset.system = scene_frame_system(*state).id;
            }
            return;
        }
    }

    // What is queued together up to the next of those: one drop, placed as a whole (drop_files.h)
    DropSystem loaded[SYSTEM_MAX];
    md_str_t run_paths[SYSTEM_MAX][RUN_MAX];
    const size_t num_loaded = drop_systems(loaded, run_paths, *state);
    DropFiles* drop = new DropFiles();
    drop_files_init(drop, loaded, num_loaded);
    while (!file_queue_empty(&state->file_queue)) {
        const FileQueue::Entry e = file_queue_front(&state->file_queue);
        if (md_str_eq_ignore_case(md_path_ext(e.path), WORKSPACE_FILE_EXTENSION)) break;
        loader::LoaderState ls = {};
        loader::init(&ls, e.path);
        if ((e.flags & FileFlags_ShowDialogue) || (ls.flags & LoaderFlag_RequiresDialogue)) break;
        file_queue_pop(&state->file_queue);
        ls.flags |= (e.flags & FileFlags_DisableCacheWrite) ? LoaderFlag_DisableCacheWrite : 0;
        uint32_t target = (e.flags & FileFlags_NewSystem) ? PLACEMENT_NEW : e.system != SYSTEM_ID_NONE ? e.system : PLACEMENT_AUTO;
        if (e.flags & FileFlags_Dropped) {
            if (queue.drop_new) target = PLACEMENT_NEW;
            else if (queue.drop_system != SYSTEM_ID_NONE) target = queue.drop_system;
        }
        drop_files_add(drop, e.path, ls, target);
    }
    PlacementPlan* plan = new PlacementPlan();
    drop_files_plan(plan, drop, SYSTEM_MAX);
    // The systems are loaded into from here: what the drop was matched against is not looked at again
    drop->loaded = nullptr;
    drop->num_loaded = 0;
    if (drop_files_load(state, drop, *plan)) {
        placement_prompt_open(state, drop, plan);
        return;
    }
    delete plan;
    drop_files_free(drop);
    delete drop;
}

const SystemContext& scene_frame_system(const ApplicationState& app) {
    for (size_t i = 0; i < app.num_systems; ++i) {
        if (system_loaded(*app.systems[i])) return *app.systems[i];
    }
    return app.system_context;
}

void reset_view_scene(ViewTransform* transform, const ApplicationState& app, float fov_y, const ViewRegion& region) {
    MD_ASSERT(transform);
    // Every instance shown, of every system shown
    const SystemInstance* shown[SYSTEM_MAX * INSTANCE_MAX];
    size_t num_shown = 0;
    for (size_t i = 0; i < app.num_systems; ++i) {
        const SystemContext& sc = *app.systems[i];
        if (!sc.visible) continue;
        for (size_t k = 0; k < sc.num_instances; ++k) {
            const SystemInstance& inst = *sc.instances[k];
            if (inst.scene.visible && inst.state.num_atoms > 0) shown[num_shown++] = &inst;
        }
    }
    if (num_shown == 0) return;
    if (num_shown == 1 && md_float4x4_equal(shown[0]->scene.placement, md_float4x4_ident())) {
        const SystemInstance& inst = *shown[0];
        const SystemContext& sc = *inst.sys;
        const md_bitset_t visible = md_bitmask_view(&sc.representation.visibility_mask);
        reset_view(transform, inst.state, &visible, fov_y, region, &sc.data.sys, inst.scene.assembly);
        if (&sc != &scene_frame_system(app)) {
            // Its cell is not the frame of the world: reset_view placed it as if it were
            const md_float3_t own   = md_float4x4_mul_float3(inst.unitcell_transform, md_float3_set1(0.0f), 1.0f);
            const md_float3_t world = md_float4x4_mul_float3(system_instance(scene_frame_system(app)).unitcell_transform, md_float3_set1(0.0f), 1.0f);
            transform->position = transform->position - own + world;
        }
        return;
    }

    // The visible atoms of every instance shown, as they are drawn: in their own coordinates, in the frame of the
    // world. Of large systems every k-th, at most about a million in all.
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    size_t total = 0;
    for (size_t i = 0; i < num_shown; ++i) {
        const md_bitset_t visible = md_bitmask_view(&shown[i]->sys->representation.visibility_mask);
        const size_t n = md_bitset_count(visible);
        total += n ? n : shown[i]->state.num_atoms;
    }
    const size_t max_points = 1000000;
    const size_t stride = MD_MAX((size_t)1, MD_DIV_UP(total, max_points));
    md_array(md_float3_t) xyz = 0;
    md_array_ensure(xyz, MD_MIN(total, max_points) + num_shown, md_temp_allocator(temp));
    for (size_t i = 0; i < num_shown; ++i) {
        const SystemInstance& inst = *shown[i];
        const md_system_state_t& st = inst.state;
        const md_bitset_t visible = md_bitmask_view(&inst.sys->representation.visibility_mask);
        const bool all = md_bitset_count(visible) == 0;
        size_t k = 0;
        for (size_t a = 0; a < st.num_atoms; ++a) {
            if (!all && !md_bitset_test(visible, a)) continue;
            if (k++ % stride) continue;
            // Where the instance is placed
            md_array_push(xyz, md_float4x4_mul_float3(inst.scene.placement, st.xyz[a], 1.0f), md_temp_allocator(temp));
        }
    }
    const size_t n = md_array_size(xyz);
    if (n == 0) return;
    *transform = camera_compute_default_view(xyz, n, nullptr, n, nullptr, fov_y, region);
    transform->position = md_float4x4_mul_float3(system_instance(scene_frame_system(app)).unitcell_transform, transform->position, 1.0f);
}

void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitset_t* mask, float fov_y, const ViewRegion& region, const md_system_t* sys, int assembly) {
    MD_ASSERT(transform);
    if (!state.num_atoms) return;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    if (sys && assembly >= 0 && (size_t)assembly < sys->assembly.count && sys->atom.count == state.num_atoms) {
        // The atoms of the instances of each group where each of its operators puts them, framed as one: the cell has
        // nothing to say about an assembly. Of a large one (a virus) every k-th atom, at most about a million in all.
        const md_assembly_data_t& a = sys->assembly;
        const md_bitset_t* m = (mask && md_bitset_count(*mask) > 0) ? mask : nullptr;
        // The groups, their instances and operators from the offsets as they are (assembly_group)
        const size_t part_beg = a.group_offset[assembly];
        const size_t part_end = a.group_offset[assembly + 1];
        size_t total = 0;
        for (size_t p = part_beg; p < part_end; ++p) {
            const AssemblyGroup group = assembly_group(a, p);
            size_t atoms = 0;
            for (const md_instance_idx_t* k = group.inst; k < group.inst_end; ++k) {
                const md_urange_t r = md_system_instance_atom_range(sys, (size_t)*k);
                atoms += r.end - r.beg;
            }
            total += atoms * (size_t)(group.oper_end - group.oper);
        }
        const size_t max_points = 1000000;
        const uint32_t stride = (uint32_t)MD_MAX((size_t)1, MD_DIV_UP(total, max_points));
        md_array(md_float3_t) xyz = 0;
        md_array_ensure(xyz, MD_MIN(total, max_points) + 1, md_temp_allocator(temp));
        for (size_t p = part_beg; p < part_end; ++p) {
            const AssemblyGroup group = assembly_group(a, p);
            for (const uint32_t* o = group.oper; o < group.oper_end; ++o) {
                const md_float4x4_t M = a.oper_transform[*o];
                for (const md_instance_idx_t* k = group.inst; k < group.inst_end; ++k) {
                    const md_urange_t r = md_system_instance_atom_range(sys, (size_t)*k);
                    for (uint32_t i = r.beg; i < r.end; i += stride) {
                        if (m && !md_bitset_test(*m, i)) continue;
                        md_array_push(xyz, md_float4x4_mul_float3(M, state.xyz[i], 1.0f), md_temp_allocator(temp));
                    }
                }
            }
        }
        const size_t n = md_array_size(xyz);
        if (n > 0) {
            md_float3x3_t A = {};
            md_unitcell_A_extract(A.elem, &state.unitcell);
            const md_float3_t cell_offset = -md_float3x3_mul_float3(A, md_float3_set1(0.5f));
            *transform = camera_compute_default_view(xyz, n, nullptr, n, nullptr, fov_y, region);
            transform->position = transform->position + cell_offset;
            return;
        }
        // None of the atoms shown is in the assembly: those shown
    }

    // A mask that selects nothing or everything is the whole system
    const int32_t* indices = nullptr;
    size_t count = state.num_atoms;
    if (mask) {
        const size_t popcount = md_bitset_count(*mask);
        if (0 < popcount && popcount < state.num_atoms) {
            int32_t* idx = md_temp_alloc_array(temp, int32_t, popcount);
            size_t len = md_bitset_iter_extract_indices(idx, popcount, md_bitset_iter(*mask));
            if (len > popcount || len > state.num_atoms) {
                MD_LOG_DEBUG("Error: Invalid number of indices");
                len = MD_MIN(popcount, state.num_atoms);
            }
            indices = idx;
            count   = len;
        }
    }

    // The world is drawn translated so that the center of the unit cell is at the origin
    md_float3x3_t A = {};
    md_unitcell_A_extract(A.elem, &state.unitcell);
    const bool   has_cell    = (md_unitcell_flags(&state.unitcell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) != 0;
    const md_float3_t cell_offset = -md_float3x3_mul_float3(A, md_float3_set1(0.5f));

    if (indices && count <= 4) {
        // A handful of atoms has no meaningful shape: center on them, keep the current orientation, fit the distance
        md_float3_t aabb_min = {}, aabb_max = {};
        md_geom_aabb(aabb_min.elem, aabb_max.elem, (const float*)state.xyz, 0, indices, count);
        const md_float3_t center = (aabb_min + aabb_max) * 0.5f;
        transform->distance = camera_fit_distance(state.xyz, indices, count, center, transform->orientation, fov_y, region);
        transform->position = center + cell_offset + transform->orientation * md_float3_set(0, 0, transform->distance);
        return;
    }

    // See camera_compute_default_view for what it considers a good view of what
    *transform = camera_compute_default_view(state.xyz, state.num_atoms, indices, count, has_cell ? &A : nullptr, fov_y, region);
    transform->position = transform->position + cell_offset;
}

// Sets the bits of what a hit stands for, grown by the granularity: what hovering it highlights and clicking it selects
void picking_hit_mask_set(md_bitmask_t* mask, const SystemContext& sys, const PickingHit& hit, SelectionGranularity granularity) {
    MD_ASSERT(mask);
    if (hit.domain == PickingDomain_Atom) {
        int32_t atom_idx = hit.local_idx;
        if (atom_idx >= 0 && (size_t)atom_idx < sys.data.sys.atom.count) {
            // Grow the hovered atom by the current granularity so hover feedback matches what a click would select
            mask_set_atom_by_selection_granularity(mask, (size_t)atom_idx, granularity, sys.data.sys);
        }
    } else if (hit.domain == PickingDomain_Bond) {
        md_atom_pair_t pair;
        md_bond_flags_t flags;
        PickedBondSource source;
        if (picked_bond(*sys.app, hit, &pair, &flags, &source)) {
            // Grow both bond atoms by the current granularity, same as for atom hits
            mask_set_atom_by_selection_granularity(mask, (size_t)pair.idx[0], granularity, sys.data.sys);
            mask_set_atom_by_selection_granularity(mask, (size_t)pair.idx[1], granularity, sys.data.sys);
        }
    } else if (hit.domain == PickingDomain_BackboneSegment) {
        const md_protein_backbone_data_t& bb = sys.data.sys.protein_backbone;
        const size_t seg_idx = hit.local_idx;
        if (seg_idx < bb.segment.count && bb.segment.res_idx) {
            const md_residue_idx_t res_idx = bb.segment.res_idx[seg_idx];
            if (res_idx >= 0) {
                // A cartoon segment stands for the residue it stems from: at least that is hovered and selected
                mask_set_residue_by_selection_granularity(mask, (size_t)res_idx, granularity, sys.data.sys);
            } else {
                // Without one, for the atoms of its amino acid by their roles (md_system_protein_segment_atoms)
                md_atom_idx_t atoms[512];
                const size_t n = MD_MIN(md_system_protein_segment_atoms(atoms, MD_COUNTOF(atoms), &sys.data.sys, seg_idx), MD_COUNTOF(atoms));
                for (size_t k = 0; k < n; ++k) {
                    mask_set_atom_by_selection_granularity(mask, (size_t)atoms[k], granularity, sys.data.sys);
                }
            }
        }
    } else if (hit.domain == PickingDomain_Ring) {
        // The circle of an aromatic ring stands for the atoms of the ring
        const int32_t* ring_atoms = nullptr;
        size_t ring_count = 0;
        md_ring_flags_t ring_flags = MD_RING_FLAG_NONE;
        if (picked_ring(sys, hit.local_idx, &ring_atoms, &ring_count, &ring_flags)) {
            mask_set_ring_by_selection_granularity(mask, ring_atoms, ring_count, granularity, sys.data.sys);
        }
    } else if (hit.domain == PickingDomain_ChainGap) {
        // The gap of a chain stands for the residues on either side of its break
        const size_t res_idx = hit.local_idx;
        if (res_idx + 1 < sys.data.sys.residue.count) {
            mask_set_residue_by_selection_granularity(mask, res_idx,        granularity, sys.data.sys);
            mask_set_residue_by_selection_granularity(mask, res_idx + 1, granularity, sys.data.sys);
        }
    } else if (hit.domain == PickingDomain_CartoonGlyph) {
        // A point of the cartoon of a lipid or a nucleic acid stands for its residue, as a cartoon segment does
        const uint32_t point_idx = cartoon_glyphs_picked_point(system_instance(sys).cartoon_glyphs, hit.local_idx);
        const CartoonGlyphs& glyphs = system_instance(sys).cartoon_glyphs;
        const int32_t res_idx = cartoon_glyphs_point_residue(glyphs, point_idx);
        if (res_idx >= 0 && (size_t)res_idx < sys.data.sys.residue.count) {
            mask_set_residue_by_selection_granularity(mask, (size_t)res_idx, granularity, sys.data.sys);
        } else if (point_idx < glyphs.num_points) {
            // Without one (a nucleic backbone found from the bonds alone) for the atoms of the point
            for (uint32_t k = glyphs.point_offset[point_idx]; k < glyphs.point_offset[point_idx + 1]; ++k) {
                mask_set_atom_by_selection_granularity(mask, (size_t)glyphs.point_atom[k], granularity, sys.data.sys);
            }
        }
    }
}

void ViamdEventHandler::process_events(const viamd::Event* events, size_t num_events) {
    for (size_t i = 0; i < num_events; ++i) {
        const viamd::Event& event = events[i];
        switch (event.type) {
        case viamd::EventType_ViamdFrameTick:
            break;
        case viamd::EventType_ViamdPickingRangeReserve: {
			MD_ASSERT(event.payload_type == viamd::EventPayloadType_PickingSpace);
			PickingSpace* space = (PickingSpace*)event.payload;
            // Every system its own ranges, in the order listed: a hit names the system by the range it is in
            for (size_t s = 0; s < state->num_systems; ++s) {
                SystemContext& sc = *state->systems[s];
                for (size_t k = 0; k < sc.num_instances; ++k) sc.instances[k]->picking = {};
                if (!system_loaded(sc)) continue;
                const md_system_t& sys = sc.data.sys;
                const uint32_t id = sc.id;
                const auto& bt = sc.operations.bond_tracking;
                // The bonds of the graph, then those formed in the frame shown, as they are drawn (picked_bond). Which
                // frame that is is not settled yet: as many as have formed in any frame are reserved for.
                const md_bond_delta_t* delta = bond_tracking_delta(sc);
                const size_t num_bonds = sys.bond.count + (delta ? MD_MAX(delta->num_formed, bt.formed_capacity) : 0);
                const size_t num_frame_rings = delta ? bt.ring_capacity : 0;
                // Each instance its own: the same atom in another instance is another object to pick. Keyed by the
                // instance's place in the system (picking_instance_key), which a hit names it by.
                for (size_t k = 0; k < sc.num_instances; ++k) {
                    SystemInstance& inst = *sc.instances[k];
                    const uint64_t key = picking_instance_key(k);
                    picking_range_reserve(&inst.picking.atom, space, PickingDomain_Atom, sys.atom.count, key, id);
                    picking_range_reserve(&inst.picking.bond, space, PickingDomain_Bond, num_bonds, key, id);
                    // One index per segment, in the order md_gl draws them: the cartoon of segment i writes beg + i
                    picking_range_reserve(&inst.picking.backbone, space, PickingDomain_BackboneSegment, sys.protein_backbone.segment.count, key, id);
                    // One index per ring of the system, which the circle of an aromatic ring writes (only those are
                    // drawn), then the rings of the frame's chemistry where its bonds changed, reserved for as the bonds
                    // formed are
                    picking_range_reserve(&inst.picking.ring, space, PickingDomain_Ring, md_ring_data_count(&sys.ring) + num_frame_rings, key, id);
                    // One index per point of the cartoon of the lipids and nucleic acids, which its spheres and capsules
                    // write, and one per control point of its tubes after them (cartoon_glyphs_picked_point)
                    picking_range_reserve(&inst.picking.glyph, space, PickingDomain_CartoonGlyph, cartoon_glyphs_picking_count(inst.cartoon_glyphs), key, id);
                    // One index per residue, which the gap over a break after it writes: only those are drawn
                    picking_range_reserve(&inst.picking.gap, space, PickingDomain_ChainGap, sys.residue.count, key, id);
                }

                // One range per dipole group, keyed by the group's vector attribute, sized by that attribute's own
                // shape. An index within a range is then the element index inside the group, so a hit carries the
                // whole identity and nobody has to reproduce a flat ordering over every (group, element) pair to
                // speak about one dipole.
                DipoleGroup groups[16];
                const size_t num_groups = MD_MIN(dipole_groups_gather(groups, MD_COUNTOF(groups), sys), MD_COUNTOF(groups));
                for (size_t j = 0; j < num_groups; ++j) {
                    picking_range_reserve(NULL, space, PickingDomain_Dipole, groups[j].count, groups[j].key, id);
                }

                // Each glyph representation shown: its points, keyed by the representation. Read and placed on the frame
                // here, so that what it reserves is what it draws this frame.
                // One range per shown representation which draws weak bonds (those of a Bonds representation, the
                // interactions over the sticks of a licorice), keyed by its id: weak bond i picks as beg + i
                // (representation_weak_bond_draw_ops, picked_bond)
                for (size_t j = 0; j < md_array_size(sc.representation.reps); ++j) {
                    Representation& rep = sc.representation.reps[j];
                    if (!rep.type_is_valid || !representation_shown(sc, rep)) continue;
                    if (rep.type == RepresentationType::Glyphs && glyphs_update(&rep.glyphs, sc.data.sys, representation_instance(sc, rep).state)) {
                        picking_range_reserve(NULL, space, PickingDomain_Glyph, glyphs_count(rep.glyphs), rep.node.id, id);
                    }
                    const size_t count = representation_weak_bond_count(sc, rep);
                    if (count > 0) picking_range_reserve(NULL, space, PickingDomain_Bond, count, rep.node.id, id);
                }
            }
            break;
        }
        case viamd::EventType_ViamdInteractionSurface:
            if (event.payload) {
                MD_ASSERT(event.payload_type == viamd::EventPayloadType_InteractionSurfaceEvent);
                InteractionSurfaceEvent* surf = (InteractionSurfaceEvent*)event.payload;
                switch (surf->kind) {
                case InteractionSurfaceEventKind::Hover:
                    draw_picking_tooltip_window(surf->hit, *state);
					[[fallthrough]];
                case InteractionSurfaceEventKind::Click: {
                    // The glyph under the cursor, marked; fading in where there was none
                    {
                        const uint32_t prev = state->glyph_hover.rep_id;
                        if (surf->hit.domain == PickingDomain_Glyph) {
                            state->glyph_hover.rep_id = (uint32_t)surf->hit.key;
                            state->glyph_hover.system = surf->hit.system;
                            state->glyph_hover.point  = surf->hit.local_idx;
                            if (!prev) state->selection.highlight_time = ImGui::GetTime();
                        } else {
                            state->glyph_hover = {};
                        }
                    }
                    // The highlight is of what is under the cursor: of the system hit, none of the others'. It is the
                    // intermediate mask of what a click selects, too.
                    SystemContext* hit_sys = picking_hit_system(state, surf->hit);
                    for (size_t k = 0; k < state->num_systems; ++k) md_bitmask_clear(&state->systems[k]->selection.highlight_mask);
                    if (hit_sys) picking_hit_mask_set(&hit_sys->selection.highlight_mask, *hit_sys, surf->hit, state->selection.granularity);

                    // Commit to selection mask upon shift+click release. A plain click leaves the selection alone (double-click recenters the view).
                    if (surf->kind == InteractionSurfaceEventKind::Click && surf->selection_mode != InteractionSelectionMode::None) {
                        // The single selection sequence records the order in which individual atoms were picked
                        // and is what the context menu turns into script suggestions. Only a real click advances
                        // it: a Hover event always carries selection_mode None, so doing this above would reset
                        // the sequence to the atom under the cursor on every frame the mouse crosses the molecule.
                        // Its atoms are of one system: a pick in another starts it over.
                        SingleSelectionSequence& seq = state->selection.single_selection_sequence;
                        if (hit_sys && surf->hit.domain == PickingDomain_Atom) {
                            int32_t atom_idx = surf->hit.local_idx;
                            if (atom_idx >= 0 && (size_t)atom_idx < hit_sys->data.sys.atom.count) {
                                if (seq.system != hit_sys->id) {
                                    single_selection_sequence_clear(&seq);
                                    seq.system = hit_sys->id;
                                }
                                if (surf->selection_mode == InteractionSelectionMode::Append) {
                                    single_selection_sequence_push_idx(&seq, atom_idx);
                                }
                                else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                                    single_selection_sequence_pop_idx(&seq, atom_idx);
                                }
                            }
                        }

                        if (hit_sys && picking_domain_resolves_to_atoms(surf->hit.domain)) {
                            md_bitmask_t& sel = hit_sys->selection.selection_mask;
                            const md_bitset_t hl = md_bitmask_view(&hit_sys->selection.highlight_mask);
                            if (surf->selection_mode == InteractionSelectionMode::Append) {
                                md_bitmask_or(&sel, hl);
                            }
                            else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                                md_bitmask_andnot(&sel, hl);
                            }
                            else if (surf->selection_mode == InteractionSelectionMode::None) {
                                // What is selected is what was clicked, in whichever system
                                for (size_t k = 0; k < state->num_systems; ++k) md_bitmask_clear(&state->systems[k]->selection.selection_mask);
                                md_bitmask_or(&sel, hl);
                            }
                        }
                        else if (surf->hit.domain == 0) {
                            // A plain click (no modifier, no drag) or a remove-click on empty space clears the selection
                            if (surf->selection_mode == InteractionSelectionMode::None || surf->selection_mode == InteractionSelectionMode::Remove) {
                                for (size_t k = 0; k < state->num_systems; ++k) md_bitmask_clear(&state->systems[k]->selection.selection_mask);
                                single_selection_sequence_clear(&seq);
                            }
                        }
                    }
                    break;
                }
                case InteractionSurfaceEventKind::RegionSelect:
                    break;
                case InteractionSurfaceEventKind::ContextMenu:
                    break;
                case InteractionSurfaceEventKind::None: [[fallthrough]];
                default:
                    break;
                }
            }
            break;
        case viamd::EventType_ViamdPickingTooltipTextRequest: {
            MD_ASSERT(event.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
            PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)event.payload;
            if (const SystemContext* hit_sys = picking_hit_system(*state, req->hit)) {
                // Which system it is of, and which instance, when there are several to tell apart
                const bool systems = system_count(*state) > 1;
                const bool instances = hit_sys->num_instances > 1;
                if (systems || instances) {
                    char label[64] = "";
                    char inst_label[96] = "";
                    if (systems) system_label(label, sizeof(label), *hit_sys);
                    if (instances) instance_label(inst_label, sizeof(inst_label), picking_hit_instance(*hit_sys, req->hit));
                    tooltip_title(req, "%s%s%s", label, systems && instances ? TOOLTIP_SEP : "", inst_label);
                }
                fill_picking_tooltip(req, *hit_sys, req->hit);
            }
            break;
        }
        case viamd::EventType_ViamdViewFit: {
            MD_ASSERT(event.payload_type == viamd::EventPayloadType_ViewFitRequest);
            ViewFitRequest* req = (ViewFitRequest*)event.payload;
            if (req) {
                // Of every system shown, in its own coordinates: they are drawn in one frame (scene_frame_system)
                for (size_t k = 0; k < state->num_systems; ++k) {
                    const SystemContext& sc = *state->systems[k];
                    if (!sc.visible) continue;
                    const md_bitmask_t* m = nullptr;
                    switch (req->round) {
                    case ViewFitRound_Highlight:
                        m = &sc.selection.highlight_mask; break;
                    case ViewFitRound_Selection:
                        m = &sc.selection.selection_mask; break;
                    case ViewFitRound_Visible:
                        m = &sc.representation.visibility_mask; break;
                    default:
                        break;
                    }

                    const md_bitset_t bf = m ? md_bitmask_view(m) : md_bitset_empty_set();
                    size_t popcount = md_bitset_count(bf);
                    // In every instance shown, where it is placed
                    for (size_t j = 0; popcount > 0 && j < sc.num_instances; ++j) {
                        const SystemInstance& inst = *sc.instances[j];
                        if (!inst.scene.visible || inst.state.num_atoms != sc.data.sys.atom.count) continue;
                        md_float4_t* dst_xyzw = md_array_extend(req->xyzw, popcount, req->alloc);
                        if (!dst_xyzw) break;
                        extract_xyzw_from_mask(dst_xyzw, &bf, &sc.data.sys, &inst.state);
                        for (size_t i = 0; i < popcount; ++i) {
                            const md_float3_t p = md_float4x4_mul_float3(inst.scene.placement, md_float3_from_float4(dst_xyzw[i]), 1.0f);
                            dst_xyzw[i] = md_float4_from_float3(p, dst_xyzw[i].w);
                        }
                    }
                }
            }
            break;
        }
        case viamd::EventType_ViamdSystemStateChanged: {
			// Apply operators if toggled on system change

            MD_ASSERT(event.payload_type == viamd::EventPayloadType_System);
            SystemContext& sys = *(SystemContext*)event.payload;

            // The event means the first instance's state was just written from its source, in the lattice frame of
            // its cell: whatever turn the previous coordinates carried went with them.
            SystemInstance& inst = system_instance(sys);
            inst.state_rotation = md_float4x4_ident();

            // The bond changes of the frame are not here: bond_tracking_update keeps them up to date with the frame
            // shown once a frame, and they leave the bonds of the system as they are

            // Recenter, wrap, make whole and turn, in that order (see apply_state_operations). They rewrite
            // the coordinates that update_md_buffers uploads, and nothing else flags them: the synchronous
            // broadcasts of this event do not pass through the interpolation step that would otherwise
            // have set the bit.
            // Not while an assembly is shown, whose operators place the coordinates as the file gives them
            const auto& op = sys.operations;
            if (inst.scene.assembly < 0 && apply_state_operations(inst, op.recenter, op.apply_pbc, op.unwrap_structures)) {
                inst.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
            break;
        }
        default:
            break;
        }
    }
}

script::Program* script_ir_create(ApplicationState* state, md_str_t src, const script::CompileDesc& desc) {
    MD_ASSERT(state);
    script::Program* ir = script::compile(state->script.language, src, desc, state->allocator.persistent);
    md_array_push(state->script.all_irs, ir, state->allocator.persistent);
    return ir;
}

void script_set_language(ApplicationState* state, script::Language lang) {
    MD_ASSERT(state);
    if (state->script.language == lang) return;
    state->script.language = lang;
    state->editor.SetLanguage(script_editor::language(lang));
    script_editor::enable_autocomplete(state->editor, lang, &state->mold.sys);
    state->script.compile_ir = true;
    state->script.time_since_last_change = 1.0e9;     // now, not after the pause that follows an edit
    flag_all_representations_as_dirty(state);
}

void script_switch_language(ApplicationState* state, script::Language lang) {
    MD_ASSERT(state);
    if (state->script.language == lang) return;
    std::string text = state->editor.GetText();
    state->editor.SetText(state->script.other_text);
    state->script.other_text = std::move(text);
    script_set_language(state, lang);
}

void script_eval_free(ApplicationState* state) {
    // Nothing runs on them now: the caller stopped the evaluation task
    state->script.evaluating = nullptr;
    if (state->script.eval || md_array_size(state->script.run_evals)) {
        viamd::event_system_broadcast_event(viamd::EventType_ViamdScriptEvalFree, viamd::EventPayloadType_ApplicationState, state);
    }
    for (size_t i = 0; i < md_array_size(state->script.run_evals); ++i) {
        script::eval_free(state->script.run_evals[i].eval);
    }
    md_array_shrink(state->script.run_evals, 0);
    if (!state->script.eval) return;
    script::eval_free(state->script.eval);
    state->script.eval = nullptr;
    state->script.eval_run[0] = '\0';
    state->script.eval_queued = false;
}

void script_evals_interrupt(ApplicationState* state) {
    // What was still to be evaluated after the one running is not either: whatever stops them starts them again
    state->script.eval_queued = false;
    if (state->script.eval) script::eval_interrupt(state->script.eval);
    for (size_t i = 0; i < md_array_size(state->script.run_evals); ++i) {
        state->script.run_evals[i].queued = false;
        script::eval_interrupt(state->script.run_evals[i].eval);
    }
}

script::Eval* script_eval_of_run(const ApplicationState* state, md_str_t run) {
    if (md_str_empty(run) || md_str_eq_cstr(run, state->script.eval_run)) return state->script.eval;
    for (size_t i = 0; i < md_array_size(state->script.run_evals); ++i) {
        if (md_str_eq_cstr(run, state->script.run_evals[i].run)) return state->script.run_evals[i].eval;
    }
    return nullptr;
}

// Whether an instance of the system shows the run
static bool run_shown(const SystemContext& sys, md_str_t run) {
    for (size_t k = 0; k < sys.num_instances; ++k) {
        const SystemRun* r = instance_run(*sys.instances[k]);
        if (r && md_str_eq_cstr(run, r->path)) return true;
    }
    return false;
}

void script_run_evals_sync(ApplicationState* state) {
    auto& sc = state->script;
    const SystemContext& sys = state->system_context;
    // Of the script evaluated, as it was made: an evaluation of a script being made anew is made with it (main.cpp)
    if (!sc.eval || !script::valid(sc.eval_ir) || script::eval_fingerprint(sc.eval) != script::fingerprint(sc.eval_ir)) return;

    for (size_t i = 0; i < md_array_size(sc.run_evals);) {
        ScriptRunEval& e = sc.run_evals[i];
        const md_str_t run = md_str_from_cstr(e.run);
        if (run_shown(sys, run) && !md_str_eq_cstr(run, sc.eval_run)) {
            ++i;
            continue;
        }
        if (sc.evaluating == e.eval) {
            // Freed once it has stopped
            script::eval_interrupt(e.eval);
            ++i;
            continue;
        }
        viamd::event_system_broadcast_event(viamd::EventType_ViamdScriptEvalFree, viamd::EventPayloadType_ApplicationState, state);
        script::eval_free(e.eval);
        md_array_swap_back_and_pop(sc.run_evals, i);
        // The caches are keyed on versions within its table, and the next may be allocated where it was
        series_cache_free(state);
    }

    for (size_t k = 0; k < sys.num_instances; ++k) {
        const SystemRun* r = instance_run(*sys.instances[k]);
        if (!r || r->path[0] == '\0') continue;
        const md_str_t run = md_str_from_cstr(r->path);
        if (script_eval_of_run(state, run)) continue;
        const script::EvalDesc desc = { .num_frames = run_num_frames(sys, r), .sys = &state->mold.sys, .run = run,
                                        .executor = task_system::executor() };
        if (desc.num_frames == 0) continue;
        ScriptRunEval e;
        md_str_copy_to_char_buf(e.run, sizeof(e.run), run);
        e.eval = script::eval_create(nullptr, sc.eval_ir, desc, state->allocator.persistent);
        if (!e.eval) continue;
        md_array_push(sc.run_evals, e, state->allocator.persistent);
        series_cache_free(state);
    }
}

void script_ir_collect(ApplicationState* state) {
    MD_ASSERT(state);

    // Freed once neither field points at it and no task holds it: a long evaluation keeps the IR it runs on, not every
    // IR compiled while it runs
    for (size_t i = 0; i < md_array_size(state->script.all_irs);) {
        script::Program* ir = state->script.all_irs[i];
        if (ir != state->script.ir && ir != state->script.eval_ir && script::program_users(ir) == 0) {
            script::program_free(ir);
            md_array_swap_back_and_pop(state->script.all_irs, i);
        } else {
            ++i;
        }
    }
}

const script::Program* script_ir_of(const ApplicationState* state, script::Ref ref) {
    MD_ASSERT(state);
    if (script::ref_valid(state->script.ir, ref))      return state->script.ir;
    if (script::ref_valid(state->script.eval_ir, ref)) return state->script.eval_ir;
    return nullptr;
}

void script_visualize_ref(ApplicationState* state, script::Ref ref, int subidx, uint32_t flags) {
    MD_ASSERT(state);
    // Resolved in the IR which made it: an editor token and a plotted property may come from different ones
    if (!script_ir_of(state, ref)) {
        return;
    }
    ScriptVisTarget& t = state->script.vis_target;
    t = {};
    t.kind   = ScriptVisTarget::Ref;
    t.ref    = ref;
    t.subidx = subidx;
    t.flags  = flags;
}

void script_visualize_str(ApplicationState* state, md_str_t str, uint32_t flags) {
    MD_ASSERT(state);
    if (md_str_empty(str)) {
        return;
    }
    ScriptVisTarget& t = state->script.vis_target;
    t = {};
    t.kind  = ScriptVisTarget::Str;
    t.str   = md_alloc_str(state->allocator.frame, str);
    t.flags = flags;
}

// A visualization computed off the main thread. Everything it reads is its own or outlives it: a copy of the
// atoms' state, the IR (retained until the job is freed, see script_ir_collect), the system (which stays until the
// pool has finished, see interrupt_async_tasks). Everything it makes is in its arena.
struct ScriptVisJob {
    uint64_t key = 0;               // what it is of: the target, and the IR the target is resolved in
    uint64_t state_key = 0;         // where the atoms were
    md_allocator_i* arena = nullptr;
    ScriptVisTarget target = {};    // str in the arena
    const script::Program* ir = nullptr;
    const md_system_t* sys = nullptr;
    md_str_t run = {};                 // the run v2 reads frames of (at(), an sdf's reference), in the arena
    md_system_state_t state = {};
    script::Drawing vis = {};
    md_str_t value_text = {};          // Ref: what the value is where the atoms were (script::value_text), in the arena
    script::ValuePlot plot = {};    // Ref: its numbers plotted beside that (script::value_plot), in the arena
    double value_scale = 1.0;       // Ref: into the display unit of its unit, as the mark's text shows it
    bool ok = false;
};

static void vis_job_free(ScriptVisJob* job) {
    if (job) {
        script::program_release(job->ir);     // held from dispatch: its result is keyed by it
        md_arena_allocator_destroy(job->arena);
        delete job;
    }
}

// The IR a target is resolved in: the one which made the reference, or the script being edited for an expression
static const script::Program* vis_target_ir(const ApplicationState* state, const ScriptVisTarget& t) {
    return t.kind == ScriptVisTarget::Ref ? script_ir_of(state, t.ref) : state->script.ir;
}

static uint64_t vis_target_key(const ApplicationState* state, const ScriptVisTarget& t) {
    const script::Program* ir = vis_target_ir(state, t);
    const uintptr_t ir_ptr = (uintptr_t)ir;
    const uint64_t fingerprint = script::fingerprint(ir);
    uint64_t h = md_hash64(&t.kind, sizeof(t.kind), 0);
    h = md_hash64(&t.subidx, sizeof(t.subidx), h);
    h = md_hash64(&t.flags, sizeof(t.flags), h);
    h = md_hash64(&t.ref.program, sizeof(t.ref.program), h);
    h = md_hash64(&t.ref.node, sizeof(t.ref.node), h);
    if (t.kind == ScriptVisTarget::Str) {
        h = md_hash64(t.str.ptr, t.str.len, h);
    }
    h = md_hash64(&ir_ptr, sizeof(ir_ptr), h);
    h = md_hash64(&fingerprint, sizeof(fingerprint), h);
    // The value's text is in the display units
    const uint64_t units = display_units::version();
    h = md_hash64(&units, sizeof(units), h);
    return h;
}

// Where the atoms are. Hashed rather than tracked: nothing that moves them can be missed. Only while something
// is visualized.
static uint64_t vis_state_key(const md_system_state_t& st) {
    // The cell field by field: the struct has padding
    const md_unitcell_t& uc = st.unitcell;
    const double cell[6] = { uc.x, uc.xy, uc.xz, uc.y, uc.yz, uc.z };
    const uint32_t cell_flags = (uint32_t)uc.flags;
    uint64_t h = md_hash64(cell, sizeof(cell), 0);
    h = md_hash64(&cell_flags, sizeof(cell_flags), h);
    if (st.xyz && st.num_atoms) {
        h = md_hash64(st.xyz, st.num_atoms * sizeof(md_float3_t), h);
    }
    return h;
}

static void vis_job_dispatch(ApplicationState* state, uint64_t key, uint64_t state_key) {
    const ScriptVisTarget& t = state->script.vis_target;
    const script::Program* ir = vis_target_ir(state, t);
    if (t.kind == ScriptVisTarget::Ref && !ir) {
        return;
    }

    ScriptVisJob* job = new ScriptVisJob();
    job->key = key;
    job->state_key = state_key;
    job->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    job->target = t;
    job->target.str = md_alloc_str(job->arena, t.str);
    job->ir = ir;
    script::program_retain(ir);
    job->sys = &state->mold.sys;
    job->run = md_alloc_str(job->arena, run_path(state->system_context));
    job->state.alloc = job->arena;
    md_system_state_copy(&job->state, &system_instance(state->system_context).state);
    script::drawing_init(&job->vis, job->arena);
    if (t.kind == ScriptVisTarget::Ref) {
        // The value in the unit its mark shows it in
        const size_t nm = script::num_marks(ir);
        const script::Mark* mk = script::marks(ir);
        for (size_t i = 0; i < nm; ++i) {
            if (mk[i].ref.program == t.ref.program && mk[i].ref.node == t.ref.node && !md_unit_is_none(mk[i].unit)) {
                job->value_scale = display_units::factor(nullptr, mk[i].unit);
                break;
            }
        }
    }

    state->script.vis_task.start({}, [job]() {
        script::DrawCtx ctx = {
            .prog  = job->ir,
            .sys   = job->sys,
            .state = &job->state,
            .run   = job->run,
        };
        ctx.value_scale = job->value_scale;
        const ScriptVisTarget& tgt = job->target;
        job->ok = (tgt.kind == ScriptVisTarget::Ref)
            ? script::draw_ref(&job->vis, tgt.ref, tgt.subidx, ctx, tgt.flags)
            : script::draw_text(&job->vis, tgt.str, ctx, tgt.flags);
        if (job->ok && tgt.kind == ScriptVisTarget::Ref) {
            char buf[2048];
            const size_t len = script::value_text(buf, sizeof(buf), tgt.ref, tgt.subidx, ctx);
            if (len) {
                job->value_text = md_alloc_str(job->arena, {buf, len});
            }
        }
        if (tgt.kind == ScriptVisTarget::Ref) {
            // Without the evaluation, which may be freed or run meanwhile: a value over the trajectory is plotted on
            // the main thread (script_vis_value_plot)
            script::value_plot(&job->plot, tgt.ref, tgt.subidx, ctx, job->arena);
        }
    });

    state->script.vis_running = job;
}

bool energy_term_groups(md_edr_term_groups_t* out, const ApplicationState* state, const md_attribute_t* attr) {
    MD_ASSERT(out);
    MD_ASSERT(state);
    if (!attr || state->mold.sys.group.count == 0) return false;
    // The label is the name as the file has it; an unlabelled series has its leaf as its name
    const md_str_t name = md_str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
    return md_edr_term_groups(out, &state->mold.sys, name);
}

bool energy_term_names_groups(const md_attribute_t* attr) {
    if (!attr) return false;
    return md_edr_term_is_of_groups(md_str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label);
}

size_t energy_term_groups_text(char* buf, size_t cap, const ApplicationState* state, const md_edr_term_groups_t& tg) {
    MD_ASSERT(state);
    const md_group_data_t& g = state->mold.sys.group;
    if (!cap) return 0;
    buf[0] = '\0';
    if (tg.group[0] >= g.count || tg.group[1] >= g.count) return 0;
    const md_str_t a = g.name[tg.group[0]], b = g.name[tg.group[1]];
    int n = (tg.count == 2 && tg.group[0] != tg.group[1]) ? snprintf(buf, cap, MD_STR_FMT " and " MD_STR_FMT, MD_STR_ARG(a), MD_STR_ARG(b))
          : tg.count == 2 ? snprintf(buf, cap, MD_STR_FMT " with itself", MD_STR_ARG(a))
          : snprintf(buf, cap, MD_STR_FMT, MD_STR_ARG(a));
    return n > 0 ? MD_MIN((size_t)n, cap - 1) : 0;
}

bool highlight_energy_term_groups(ApplicationState* state, const md_attribute_t* attr) {
    MD_ASSERT(state);
    md_edr_term_groups_t tg;
    if (!energy_term_groups(&tg, state, attr)) return false;
    const md_group_data_t& g = state->mold.sys.group;
    const uint32_t count = (tg.count == 2 && tg.group[0] != tg.group[1]) ? 2 : 1;
    // Group 0 in the highlight's own colour, group 1 in the first colour of the groups after it
    md_bitset_t* sets = (md_bitset_t*)md_alloc(state->allocator.frame, 2 * sizeof(md_bitset_t));
    sets[0] = g.atoms[tg.group[0]];
    sets[1] = g.atoms[tg.group[1]];
    md_bitmask_assign(&state->system_context.selection.highlight_mask, sets[0]);
    if (count == 2) md_bitmask_or(&state->system_context.selection.highlight_mask, sets[1]);
    if (count == 2) {
        // Indices of groups: of the topology and of the groups as they are now
        uint64_t key = md_hash64(tg.group, sizeof(tg.group), state->mold.sys.topology_version ^ (state->mold.sys.group.version * 0x9E3779B97F4A7C15ull));
        state->selection.highlight_groups = {sets, count, key ? key : 1};
    }
    state->selection.term_groups_highlighted = true;
    return true;
}

void script_vis_begin_frame(ApplicationState* state) {
    MD_ASSERT(state);
    state->script.vis_target = {};
    state->script.vis_shown = nullptr;
    state->selection.highlight_groups = {};
}

void script_vis_update(ApplicationState* state) {
    MD_ASSERT(state);
    auto& sc = state->script;

    // The groups of an energy term are hovered no more
    auto& sel = state->selection;
    if (sel.term_groups_highlighted_prev && !sel.term_groups_highlighted) {
        md_bitmask_clear(&state->system_context.selection.highlight_mask);
    }
    sel.term_groups_highlighted_prev = sel.term_groups_highlighted;
    sel.term_groups_highlighted = false;

    // A finished computation becomes the result
    if (sc.vis_running && !sc.vis_task.running()) {
        vis_job_free(sc.vis_result);
        sc.vis_result = sc.vis_running;
        sc.vis_running = nullptr;
    }

    sc.vis_shown = nullptr;
    if (sc.vis_target.kind == ScriptVisTarget::None) {
        return;
    }

    const uint64_t key = vis_target_key(state, sc.vis_target);
    const uint64_t state_key = vis_state_key(system_instance(state->system_context).state);

    // At most one at a time. While one runs, whatever is set when it has finished is what comes next: the
    // targets in between are never computed.
    const bool computed = sc.vis_result && sc.vis_result->key == key && sc.vis_result->state_key == state_key;
    if (!sc.vis_running && !computed) {
        vis_job_dispatch(state, key, state_key);
    }

    // Shown only if it is of what is set. For where the atoms were, while it is computed for where they are.
    if (sc.vis_result && sc.vis_result->ok && sc.vis_result->key == key) {
        sc.vis_shown = &sc.vis_result->vis;
        if (!md_bitset_empty(sc.vis_shown->atom_mask)) {
            md_bitmask_assign(&state->system_context.selection.highlight_mask, sc.vis_shown->atom_mask);
            // Its elements, a colour each, when there is more than one. Known by what is visualized, not where the
            // atoms are: playing, the atoms' flags are not made again every frame for them
            const size_t num_groups = md_array_size(sc.vis_shown->groups);
            if (num_groups > 1) {
                const uint64_t key = md_hash64(&num_groups, sizeof(num_groups), sc.vis_result->key);
                state->selection.highlight_groups = {sc.vis_shown->groups, num_groups, key ? key : 1};
            }
        }
    }
}

uint64_t script_vis_shown_key(const ApplicationState* state) {
    MD_ASSERT(state);
    const auto& sc = state->script;
    if (!sc.vis_shown || !sc.vis_result) return 0;
    const uint64_t key = md_hash64(&sc.vis_result->state_key, sizeof(sc.vis_result->state_key), sc.vis_result->key);
    return key ? key : 1;
}

md_str_t script_vis_value_text(ApplicationState* state) {
    MD_ASSERT(state);
    const auto& sc = state->script;
    if (sc.vis_target.kind != ScriptVisTarget::Ref || !sc.vis_result || !sc.vis_result->ok) {
        return {};
    }
    // Of what is set: while the next one is computed, the text is that of the atoms where they were
    if (sc.vis_result->key != vis_target_key(state, sc.vis_target)) {
        return {};
    }
    return sc.vis_result->value_text;
}

const script::ValuePlot* script_vis_value_plot(ApplicationState* state) {
    MD_ASSERT(state);
    auto& sc = state->script;
    if (sc.vis_target.kind != ScriptVisTarget::Ref) return nullptr;
    const uint64_t key = vis_target_key(state, sc.vis_target);
    if (sc.vis_result && sc.vis_result->key == key && sc.vis_result->plot.kind != script::PlotKind::None) {
        return &sc.vis_result->plot;
    }
    // Over the trajectory: the evaluation's, read here, where it is freed and its runs are started, while none is
    // queued or running, and only if it is of the program the reference was made in
    const script::Program* ir = script_ir_of(state, sc.vis_target.ref);
    script::Eval* eval = sc.eval;
    if (!ir || ir != sc.eval_ir || !eval || state->tasks.evaluate.running() || script::eval_fingerprint(eval) != script::fingerprint(ir)) {
        return nullptr;
    }
    // A run's errors change version with every run that ends: what the plot is of
    const uint64_t run = script::eval_errors(eval).version;
    if (!run) return nullptr;
    const int64_t frame = (int64_t)floor(system_instance(state->system_context).state.frame);
    uint64_t h = md_hash64(&run, sizeof(run), key);
    h = md_hash64(&eval, sizeof(eval), h);
    h = md_hash64(&frame, sizeof(frame), h);
    if (h != sc.trajectory_plot_key) {
        sc.trajectory_plot_key = h;
        if (!sc.trajectory_plot_arena) sc.trajectory_plot_arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
        md_arena_allocator_reset(sc.trajectory_plot_arena);
        script::DrawCtx ctx = {
            .prog  = ir,
            .sys   = &state->mold.sys,
            .state = &system_instance(state->system_context).state,
            .eval  = eval,
            .run   = run_path(state->system_context),
        };
        script::value_plot(&sc.trajectory_plot, sc.vis_target.ref, sc.vis_target.subidx, ctx, sc.trajectory_plot_arena);
    }
    return sc.trajectory_plot.kind != script::PlotKind::None ? &sc.trajectory_plot : nullptr;
}

void script_vis_reset(ApplicationState* state) {
    MD_ASSERT(state);
    auto& sc = state->script;
    if (sc.vis_running) {
        sc.vis_task.wait();
        vis_job_free(sc.vis_running);
        sc.vis_running = nullptr;
    }
    vis_job_free(sc.vis_result);
    sc.vis_result = nullptr;
    sc.vis_shown = nullptr;
    sc.vis_target = {};
    if (sc.trajectory_plot_arena) md_arena_allocator_destroy(sc.trajectory_plot_arena);
    sc.trajectory_plot_arena = nullptr;
    sc.trajectory_plot = {};
    sc.trajectory_plot_key = 0;
}

void script_set_hovered_property(ApplicationState* state, md_str_t label, int population_idx) {
    // The label, a script property's identifier, is copied. One too long to fit is taken as none: the views match
    // it against identifiers, and a cut one would match the wrong property or nothing.
    char* dst = state->hovered_property_label;
    const size_t cap = sizeof(state->hovered_property_label);
    if (label.ptr && label.len < cap) {
        md_memcpy(dst, label.ptr, label.len);
        dst[label.len] = '\0';
    } else {
        dst[0] = '\0';
    }
    state->hovered_property_pop_idx = population_idx;
}
