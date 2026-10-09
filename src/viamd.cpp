
#include <md_util.h>
#include <md_qm.h>

#include <stdarg.h>

#include <md_filter.h>

#include <core/md_log.h>
#include <core/md_str_builder.h>
#include <core/md_arena_allocator.h>

#include <viamd.h>
#include <viamd_event.h>
#include <event.h>
#include <loader.h>
#include <color_utils.h>
#include <display_units.h>
#include <serialization_utils.h>

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>

#include <imgui.h>
#include <imgui_notify.h>
#include <implot.h>
#include <implot_internal.h>

mat3_t mat3_PCA(const vec4_t* xyzw, size_t count) {
    vec3_t acc = vec3_zero();
    for (size_t i = 0; i < count; ++i) {
        acc = vec3_add(acc, vec3_from_vec4(xyzw[i]));
    }
    vec3_t mean = acc / (float)count;

    mat3_t cov = mat3_covariance_matrix_vec4(xyzw, nullptr, count, mean);
    mat3_eigen_t eigen = mat3_eigen(cov);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    return PCA;
}

void calculate_bounds(float out_min[3], float out_max[3], const vec4_t* xyzw, size_t count, const mat3_t& orientation) {
    vec4_t min_v = vec4_set1(FLT_MAX);
    vec4_t max_v = vec4_set1(-FLT_MAX);

    mat4_t rot = mat4_from_mat3(mat3_transpose(orientation));

    for (size_t i = 0; i < count; ++i) {
        vec4_t v = mat4_mul_vec4(rot, xyzw[i]);
        min_v = vec4_min(min_v, v);
        max_v = vec4_max(max_v, v);
    }

    // Padding
    const float pad = 6.0f;
    min_v -= pad;
    max_v += pad;

    MEMCPY(out_min, &min_v, sizeof(float) * 3);
    MEMCPY(out_max, &max_v, sizeof(float) * 3);
}

// Construct texture to world transformation matrix for Volume
// extent is the extent of the volume (dim * voxel_size)
mat4_t compute_texture_to_world_mat(const mat3_t& orientation, const vec3_t& origin, const vec3_t& extent) {
    mat4_t T = mat4_translate_vec3(origin);
    mat4_t R = mat4_from_mat3(orientation);
    mat4_t S = mat4_scale_vec3(extent);
    return T * R * S;
}

mat4_t compute_world_to_model_mat(const mat3_t& orientation, const vec3_t& origin) {
    mat4_t world_to_model = mat4_from_mat3(mat3_transpose(orientation)) * mat4_translate_vec3(-origin);
    return world_to_model;
}

mat4_t compute_index_to_world_mat(const mat3_t& orientation, const vec3_t& in_origin, const vec3_t& stepsize) {
    vec3_t step_x = orientation.col[0] * stepsize.x;
    vec3_t step_y = orientation.col[1] * stepsize.y;
    vec3_t step_z = orientation.col[2] * stepsize.z;
    // Shift origin by half voxel
    vec3_t origin = in_origin + orientation * (stepsize * 0.5f);

    mat4_t index_to_world = {
        step_x.x, step_x.y, step_x.z, 0.0f,
        step_y.x, step_y.y, step_y.z, 0.0f,
        step_z.x, step_z.y, step_z.z, 0.0f,
        origin.x, origin.y, origin.z, 1.0f,
    };

    return index_to_world;
}

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
void compute_dim(int out_dim[3], const vec3_t& in_ext, double samples_per_unit_length) {
    out_dim[0] = CLAMP(ALIGN_TO((int)(in_ext.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = CLAMP(ALIGN_TO((int)(in_ext.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = CLAMP(ALIGN_TO((int)(in_ext.z * samples_per_unit_length), 8), 8, 512);
}

// Grid units are BOHR, matching what md_gto evaluates in; a Volume's transforms are Angstrom,
// matching the world the camera lives in. init_volume is where the two meet, and it is the only
// place that conversion belongs.
void init_grid(md_grid_t* grid, const mat3_t& orientation, const vec3_t& min_ext, const vec3_t& max_ext, double samples_per_unit_length) {
    ASSERT(grid);
    vec3_t extent = max_ext - min_ext;
    compute_dim(grid->dim, extent, samples_per_unit_length);
    vec3_t voxel_size = vec3_div(extent, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
    grid->orientation = orientation;
    grid->origin = orientation * min_ext;
    grid->spacing = voxel_size;
}

void init_volume(Volume* vol, const md_grid_t& grid, GLenum format) {
    ASSERT(vol);
    MEMCPY(vol->dim, grid.dim, sizeof(vol->dim));

    // A grid is in BOHR and a Volume's transforms are in Angstrom; this is where the two meet.
    const float scl = (float)BOHR_TO_ANGSTROM;

    vec3_t extent = md_grid_extent(&grid);
    vol->world_to_model   = compute_world_to_model_mat(grid.orientation, grid.origin * scl);
    vol->texture_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
    vol->voxel_size       = grid.spacing * scl;
    gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
    volume::notify_data_changed(vol->tex_id);
}

static void init_all_representations(ApplicationState* state);

// PICKING TOOLTIP
//
// The tooltip describes the object under the cursor: an atom, a bond, or a backbone segment of the cartoon or ribbons,
// which stands for its residue. What a click selects is the highlight's to show, as that follows the selection
// granularity; the tooltip does not change with it, so the same atom always reads the same.
//
// The title names the object and where it sits (residue, chain). Each row gives one property, grouped by what it means
// rather than by the bitfield it is stored in: kind, role, chemistry, region, source. Keys, and values that are our own
// words, are lower case; names from the file (atoms, residues, chains, entities) and element symbols are printed as
// given. A row with nothing to say is left out. The index row spells the object in the script's terms, from 1.
//
// The font carries Latin-1, Greek and U+2010-2027 (application.cpp): the middle dot, the en dash, Å, °, ², φ and ψ
// print; an arrow, the minus sign U+2212 and the prime U+2032 do not, so '-', '/' and "'" stand in for them.

static str_t tooltip_vprintf(md_allocator_i* alloc, const char* format, va_list args) {
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
    ASSERT(req && req->alloc);
    va_list args;
    va_start(args, format);
    const TooltipLine line = {.kind = TooltipLineKind::Title, .text = tooltip_vprintf(req->alloc, format, args)};
    va_end(args);
    md_array_push(req->lines, line, req->alloc);
}

void tooltip_row(PickingTooltipTextRequest* req, str_t key, const char* format, ...) {
    ASSERT(req && req->alloc);
    va_list args;
    va_start(args, format);
    const TooltipLine line = {.kind = TooltipLineKind::Row, .key = key, .text = tooltip_vprintf(req->alloc, format, args)};
    va_end(args);
    md_array_push(req->lines, line, req->alloc);
}

// A row of a list built in sb, left out when the list is empty
static void tooltip_row_list(PickingTooltipTextRequest* req, str_t key, const md_strb_t& sb) {
    if (!md_strb_empty(sb)) {
        tooltip_row(req, key, "%s", md_strb_to_cstr(sb));
    }
}

static void tooltip_list_add(md_strb_t* sb, const char* item) {
    if (!md_strb_empty(*sb)) *sb += ", ";
    *sb += item;
}

// A lower case copy, for names the system spells capitalised ("Carbon") where the tooltip uses them as words
static str_t tooltip_lower(md_allocator_i* alloc, str_t str) {
    if (str_empty(str)) return {};
    char* buf = (char*)md_alloc(alloc, str.len + 1);
    MEMCPY(buf, str.ptr, str.len);
    buf[str.len] = '\0';
    convert_to_lower(buf, str.len);
    return {buf, str.len};
}

// " e", " Å", or nothing for a dimensionless value
static void tooltip_unit_suffix(char* buf, size_t cap, md_unit_t unit) {
    ASSERT(buf && cap > 1);
    buf[0] = '\0';
    if (md_unit_print(buf + 1, cap - 1, unit) > 0) {
        buf[0] = ' ';
    }
}

// Where an atom sits: its residue and its chain, -1 for either it is not in
struct TooltipPlace {
    int comp_idx = -1;
    int inst_idx = -1;
};

static TooltipPlace tooltip_place_of_atom(const md_system_t& sys, size_t atom_idx) {
    return {
        .comp_idx = md_component_find_by_atom_idx(&sys.component, atom_idx),
        .inst_idx = md_system_instance_find_by_atom_idx(&sys, atom_idx),
    };
}

// "ALA 42", or the component's index when it has no name
static void tooltip_append_residue(md_strb_t* sb, const md_system_t& sys, int comp_idx) {
    const str_t name = md_component_name(&sys.component, comp_idx);
    if (!str_empty(name)) {
        md_strb_fmt(sb, STR_FMT " %d", STR_ARG(name), md_component_seq_id(&sys.component, comp_idx));
    } else {
        md_strb_fmt(sb, "component %d", comp_idx + 1);
    }
}

// "chain A", with the author's id where it differs ("chain C (auth A)"), or the instance's index when it has no id
static void tooltip_append_chain(md_strb_t* sb, const md_system_t& sys, int inst_idx) {
    const str_t id   = md_instance_id(&sys.instance, inst_idx);
    const str_t auth = md_instance_auth_id(&sys.instance, inst_idx);
    if (!str_empty(id)) {
        md_strb_fmt(sb, "chain " STR_FMT, STR_ARG(id));
    } else {
        md_strb_fmt(sb, "instance %d", inst_idx + 1);
    }
    if (!str_empty(auth) && !str_eq(auth, id)) {
        md_strb_fmt(sb, " (auth " STR_FMT ")", STR_ARG(auth));
    }
}

// The rest of a title: " · ALA 42 · chain A" for one place, and for the two ends of a bond that crosses over
// " · ALA 42 / GLY 43 · chain A" or " · CYS 10 · chain A / CYS 50 · chain B"
static void tooltip_append_place(md_strb_t* sb, const md_system_t& sys, const TooltipPlace& a, const TooltipPlace* b = nullptr) {
    auto sep = [sb]() { if (!md_strb_empty(*sb)) *sb += TOOLTIP_SEP; };
    const bool two_insts = b && b->inst_idx != a.inst_idx;
    const bool two_comps = b && b->comp_idx != a.comp_idx;

    if (a.comp_idx != -1) {
        sep();
        tooltip_append_residue(sb, sys, a.comp_idx);
        if (two_comps && !two_insts && b->comp_idx != -1) {
            *sb += " / ";
            tooltip_append_residue(sb, sys, b->comp_idx);
        }
    }
    if (a.inst_idx != -1) {
        sep();
        tooltip_append_chain(sb, sys, a.inst_idx);
    }
    if (two_insts) {
        *sb += " / ";
        if (b->comp_idx != -1) {
            tooltip_append_residue(sb, sys, b->comp_idx);
            if (b->inst_idx != -1) *sb += TOOLTIP_SEP;
        }
        if (b->inst_idx != -1) {
            tooltip_append_chain(sb, sys, b->inst_idx);
        }
    }
}

// "component(56) instance(1)", after whatever sb holds
static void tooltip_append_index(md_strb_t* sb, const TooltipPlace& place) {
    if (place.comp_idx != -1) md_strb_fmt(sb, "%scomponent(%d)", md_strb_empty(*sb) ? "" : " ", place.comp_idx + 1);
    if (place.inst_idx != -1) md_strb_fmt(sb, "%sinstance(%d)",  md_strb_empty(*sb) ? "" : " ", place.inst_idx + 1);
}

// The atom's name, else its element's symbol, else its index
static void tooltip_append_atom_name(md_strb_t* sb, const md_system_t& sys, size_t atom_idx) {
    str_t name = md_atom_name(&sys.atom, atom_idx);
    if (str_empty(name)) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom_idx);
        if (z) name = md_util_element_symbol(z);
    }
    if (!str_empty(name)) {
        *sb += name;
    } else {
        md_strb_fmt(sb, "atom %zu", atom_idx + 1);
    }
}

// Values of the attribute table for one atom: the field a visible representation is coloured by, at the variant it
// shows, which is the value behind the colour; then every field with one value per atom (a partial charge, a B-factor,
// an occupancy). Fields the atom's rows already give from the system itself are skipped, as are absent values.
static void tooltip_atom_attribute_rows(PickingTooltipTextRequest* req, const ApplicationState& state, uint32_t atom_idx) {
    const md_system_t& sys = state.mold.sys;

    struct Field {
        const md_attribute_t* attr;
        int variant;    // -1 for a field without a variant axis
    };
    Field fields[6];
    size_t num_fields = 0;

    auto add = [&](const md_attribute_t* attr, int variant) {
        if (!attr || num_fields == ARRAY_SIZE(fields)) return;
        for (size_t i = 0; i < num_fields; ++i) {
            if (fields[i].attr == attr && fields[i].variant == variant) return;
        }
        fields[num_fields++] = {attr, variant};
    };

    for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
        const Representation& rep = state.representation.reps[i];
        if (!rep.enabled || rep.color_mapping != ColorMapping::Attribute) continue;
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, rep.atom_attribute.key);
        if (!attr) continue;
        const int variant = attr->format.rank > 1 ? CLAMP(rep.atom_attribute.variant_idx, 0, atom_attribute_variant_count(attr) - 1) : -1;
        add(attr, variant);
    }

    md_attribute_id_t ids[32];
    const size_t num_ids = MIN(atom_attribute_query(ids, ARRAY_SIZE(ids), sys), ARRAY_SIZE(ids));
    for (size_t i = 0; i < num_ids; ++i) {
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, ids[i]);
        if (!attr || attr->format.rank != 1) continue;
        if (str_eq(attr->path, STR_LIT("atom/formal_charge")) || str_eq(attr->path, STR_LIT("atom/mass"))) continue;
        add(attr, -1);
    }

    for (size_t i = 0; i < num_fields; ++i) {
        const md_attribute_t* attr = fields[i].attr;
        const int variant = fields[i].variant;
        const md_attribute_slice_t slice = variant >= 0 ? md_attribute_slice_2((uint32_t)variant, atom_idx) : md_attribute_slice_1(atom_idx);

        float value = 0.0f;
        if (md_attribute_extract_f32(&value, 1, attr, slice, md_unit_none()) != 1 || atom_attribute_value_absent(value)) continue;

        str_t key = atom_attribute_label(attr);
        if (variant >= 0) {
            key = str_printf(req->alloc, STR_FMT " (%d)", STR_ARG(key), variant + 1);
        }
        char unit[32];
        tooltip_unit_suffix(unit, sizeof(unit), attr->unit);
        // A charge reads with its sign and to a thousandth of e; anything else by its significant digits
        const bool charge = md_unit_base_equal(attr->unit, md_unit_elementary_charge());
        tooltip_row(req, key, charge ? "%+.3f%s" : "%.4g%s", value, unit);
    }
}

static void tooltip_atom(PickingTooltipTextRequest* req, const ApplicationState& state, uint32_t atom_idx) {
    const md_system_t& sys = state.mold.sys;
    const TooltipPlace place = tooltip_place_of_atom(sys, atom_idx);

    md_strb_t title = md_strb_create(req->alloc);
    tooltip_append_atom_name(&title, sys, atom_idx);
    tooltip_append_place(&title, sys, place);
    tooltip_title(req, "%s", md_strb_to_cstr(title));

    const md_atomic_number_t z        = md_atom_atomic_number(&sys.atom, atom_idx);
    const md_particle_kind_t particle = md_atom_particle_kind(&sys.atom, atom_idx);
    const md_atom_flags_t    flags    = md_system_atom_flags(&sys, atom_idx);
    const bool nucleotide = place.comp_idx != -1 && md_component_kind(&sys.component, place.comp_idx) == MD_COMPONENT_KIND_NUCLEOTIDE;

    // What it is
    if (particle != MD_PARTICLE_ATOM) {
        tooltip_row(req, STR_LIT("particle"), "%s", md_particle_kind_name(particle));
    }
    if (z) {
        const str_t symbol = md_util_element_symbol(z);
        const str_t name   = tooltip_lower(req->alloc, md_util_element_name(z));
        tooltip_row(req, STR_LIT("element"), STR_FMT " (" STR_FMT ")", STR_ARG(symbol), STR_ARG(name));
    }
    const str_t ff_type = md_atom_type_ff_type(&sys.atom.type, (size_t)md_atom_type_idx(&sys.atom, atom_idx));
    if (!str_empty(ff_type)) {
        tooltip_row(req, STR_LIT("type"), STR_FMT, STR_ARG(ff_type));
    }
    if (particle != MD_PARTICLE_ATOM) {
        // An atom weighs what its element does, a bead what its force field says
        const float mass = md_atom_mass(&sys.atom, atom_idx);
        if (mass > 0.0f) tooltip_row(req, STR_LIT("mass"), "%.3f Da", mass);
    }

    // Its part in the residue. A nucleoside is the sugar and the base: what of it is not the base is the sugar.
    md_strb_t role = md_strb_create(req->alloc);
    if (flags & MD_ATOM_FLAG_BACKBONE)   tooltip_list_add(&role, "backbone");
    if (flags & MD_ATOM_FLAG_SIDE_CHAIN) tooltip_list_add(&role, "side chain");
    if (flags & MD_ATOM_FLAG_NUCLEOBASE) {
        tooltip_list_add(&role, "base");
    } else if (flags & MD_ATOM_FLAG_NUCLEOSIDE) {
        tooltip_list_add(&role, "sugar");
    }
    if (flags & MD_ATOM_FLAG_TERMINAL_BEG) tooltip_list_add(&role, nucleotide ? "5' end" : "N-terminus");
    if (flags & MD_ATOM_FLAG_TERMINAL_END) tooltip_list_add(&role, nucleotide ? "3' end" : "C-terminus");
    tooltip_row_list(req, STR_LIT("role"), role);

    // Its chemistry
    md_strb_t chem = md_strb_create(req->alloc);
    const md_hybridization_t hyb = md_atom_flags_hybridization(flags);
    if (hyb != MD_HYBRIDIZATION_UNKNOWN) tooltip_list_add(&chem, md_hybridization_name(hyb));
    if (flags & MD_ATOM_FLAG_AROMATIC)   tooltip_list_add(&chem, "aromatic");
    const int num_h = md_atom_hydrogen_count(&sys.atom, atom_idx);
    if (num_h > 0) md_strb_fmt(&chem, "%s%d H", md_strb_empty(chem) ? "" : ", ", num_h);
    tooltip_row_list(req, STR_LIT("chemistry"), chem);

    const int formal_charge = md_atom_formal_charge(&sys.atom, atom_idx);
    if (formal_charge) {
        tooltip_row(req, STR_LIT("formal charge"), "%+d", formal_charge);
    }
    if (flags & MD_ATOM_FLAG_QM) {
        tooltip_row(req, STR_LIT("region"), "QM");
    }

    tooltip_atom_attribute_rows(req, state, atom_idx);

    const vec3_t pos = md_state_coord(&state.mold.state, atom_idx);
    tooltip_row(req, STR_LIT("position"), "%.3f, %.3f, %.3f " TOOLTIP_ANGSTROM, pos.x, pos.y, pos.z);

    md_strb_t index = md_strb_create(req->alloc);
    md_strb_fmt(&index, "atom(%u)", atom_idx + 1);
    tooltip_append_index(&index, place);
    tooltip_row_list(req, STR_LIT("index"), index);
}

static void tooltip_bond(PickingTooltipTextRequest* req, const ApplicationState& state, uint32_t bond_idx) {
    const md_system_t& sys = state.mold.sys;
    const md_system_state_t& sys_state = state.mold.state;
    const md_atom_pair_t  pair  = sys.bond.pairs[bond_idx];
    const md_bond_flags_t flags = sys.bond.flags[bond_idx];
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

    md_strb_t title = md_strb_create(req->alloc);
    tooltip_append_atom_name(&title, sys, pair.idx[0]);
    title += symbol;
    tooltip_append_atom_name(&title, sys, pair.idx[1]);
    tooltip_append_place(&title, sys, place[0], &place[1]);
    tooltip_title(req, "%s", md_strb_to_cstr(title));

    // What it is: an aromatic or delocalized bond is that before it is of an order
    static const char* order_names[] = {"", "single", "double", "triple", "quadruple"};
    md_strb_t kind = md_strb_create(req->alloc);
    if (flags & MD_BOND_FLAG_AROMATIC) {
        tooltip_list_add(&kind, "aromatic");
    } else if (flags & MD_BOND_FLAG_DELOCALIZED) {
        tooltip_list_add(&kind, "delocalized");
    } else if (0 < order && order < (int)ARRAY_SIZE(order_names)) {
        tooltip_list_add(&kind, order_names[order]);
    }
    if (flags & MD_BOND_FLAG_COORDINATE) tooltip_list_add(&kind, "coordinate");
    tooltip_row_list(req, STR_LIT("kind"), kind);

    // Across a periodic boundary the bond is the short way round
    vec3_t d = vec3_sub(md_state_coord(&sys_state, pair.idx[1]), md_state_coord(&sys_state, pair.idx[0]));
    md_util_min_image_vec3(&d, 1, &sys_state.unitcell);
    tooltip_row(req, STR_LIT("length"), "%.3f " TOOLTIP_ANGSTROM, vec3_length(d));

    md_strb_t source = md_strb_create(req->alloc);
    tooltip_list_add(&source, md_bond_origin_name(md_bond_origin(flags)));
    if (flags & MD_BOND_FLAG_ORDER_PERCEIVED) tooltip_list_add(&source, "order perceived");
    tooltip_row_list(req, STR_LIT("source"), source);

    tooltip_row(req, STR_LIT("index"), "atom(%d) atom(%d)", pair.idx[0] + 1, pair.idx[1] + 1);
}

// A backbone segment stands for its residue: the cartoon and ribbons have no atoms to speak of. What the segment adds is
// the secondary structure shown for it (what shapes the cartoon and colours it by secondary structure) and its
// dihedrals. phi and psi need a neighbour on either side within the chain, so the ends of a chain, and chains too
// short to assign, have none (md_util_backbone_angles_compute).
static void tooltip_backbone_segment(PickingTooltipTextRequest* req, const ApplicationState& state, uint32_t seg_idx) {
    const md_system_t& sys = state.mold.sys;
    const md_protein_backbone_data_t& bb = sys.protein_backbone;
    ASSERT(seg_idx < bb.segment.count && bb.segment.comp_idx);

    const int comp_idx = bb.segment.comp_idx[seg_idx];
    if (comp_idx < 0 || (size_t)comp_idx >= sys.component.count) return;
    const md_urange_t range = md_component_atom_range(&sys.component, comp_idx);
    const TooltipPlace place = {
        .comp_idx = comp_idx,
        .inst_idx = range.beg < range.end ? md_system_instance_find_by_atom_idx(&sys, range.beg) : -1,
    };

    md_strb_t title = md_strb_create(req->alloc);
    tooltip_append_place(&title, sys, place);
    tooltip_title(req, "%s", md_strb_to_cstr(title));

    const md_component_kind_t  kind  = md_component_kind(&sys.component, comp_idx);
    const md_component_flags_t flags = md_system_component_flags(&sys, comp_idx);
    const bool nucleotide = kind == MD_COMPONENT_KIND_NUCLEOTIDE;

    md_strb_t kind_list = md_strb_create(req->alloc);
    if (kind != MD_COMPONENT_KIND_OTHER) tooltip_list_add(&kind_list, md_component_kind_name(kind));
    if ((kind == MD_COMPONENT_KIND_AMINO_ACID || nucleotide) && !(flags & MD_COMPONENT_FLAG_RESOLVED)) tooltip_list_add(&kind_list, "unresolved");
    if (flags & MD_COMPONENT_FLAG_TERMINAL_BEG) tooltip_list_add(&kind_list, nucleotide ? "5' end" : "N-terminus");
    if (flags & MD_COMPONENT_FLAG_TERMINAL_END) tooltip_list_add(&kind_list, nucleotide ? "3' end" : "C-terminus");
    tooltip_row_list(req, STR_LIT("kind"), kind_list);

    if (const md_secondary_structure_t* ss = displayed_secondary_structure(&state)) {
        tooltip_row(req, STR_LIT("sec. structure"), "%s", md_secondary_structure_name(ss[seg_idx]));
    }

    const md_backbone_angles_t* angles = md_util_state_backbone_angles(&state.mold.state, &sys);
    if (angles && bb.range.offset) {
        for (size_t i = 0; i < bb.range.count; ++i) {
            const uint32_t beg = bb.range.offset[i];
            const uint32_t end = bb.range.offset[i + 1];
            if (seg_idx < beg || end <= seg_idx) continue;
            if (end - beg >= 4 && beg < seg_idx && seg_idx + 1 < end) {
                // phi, psi and the degree sign in UTF-8 escapes
                tooltip_row(req, STR_LIT("\xCF\x86, \xCF\x88"), "%.1f\xC2\xB0, %.1f\xC2\xB0", RAD_TO_DEG(angles[seg_idx].phi), RAD_TO_DEG(angles[seg_idx].psi));
            }
            break;
        }
    }

    tooltip_row(req, STR_LIT("atoms"), "%u", range.end - range.beg);

    // The molecule it is part of, which the title names only by its chain
    if (place.inst_idx != -1) {
        const md_entity_idx_t ent_idx = md_instance_entity_idx(&sys.instance, place.inst_idx);
        if (ent_idx != -1) {
            const str_t desc = md_entity_description(&sys.entity, ent_idx);
            const md_entity_flags_t ent_flags = md_entity_flags(&sys.entity, ent_idx);
            const char* ent_kind = md_entity_kind_name(md_entity_kind(&sys.entity, ent_idx));
            const char* inferred = (ent_flags & MD_ENTITY_FLAG_INFERRED) ? ", inferred" : "";
            if (!str_empty(desc)) {
                tooltip_row(req, STR_LIT("entity"), STR_FMT " (%s%s)", STR_ARG(desc), ent_kind, inferred);
            } else {
                tooltip_row(req, STR_LIT("entity"), "%s%s", ent_kind, inferred);
            }
        }
    }

    md_strb_t index = md_strb_create(req->alloc);
    tooltip_append_index(&index, place);
    tooltip_row_list(req, STR_LIT("index"), index);
}

// The hit names the group and the element outright - its range was reserved against that group's attribute - so there
// is no list to rebuild and nothing to look the index up in.
static void tooltip_dipole(PickingTooltipTextRequest* req, const ApplicationState& state, const PickingHit& hit) {
    const md_system_t& sys = state.mold.sys;
    DipoleGroup group = {};
    if (!dipole_group_from_key(&group, sys, hit.key) || hit.local_idx >= group.count) return;

    char label[64];
    const int label_len = dipole_entry_label(label, sizeof(label), group, hit.local_idx);
    tooltip_title(req, "dipole" TOOLTIP_SEP "%.*s", label_len, label);

    vec3_t vec = {0, 0, 0};
    if (dipole_moment_read(&vec, nullptr, sys, group.key, hit.local_idx)) {
        char unit[32];
        if (md_unit_is_atomic(group.unit)) {
            snprintf(unit, sizeof(unit), " a.u.");
        } else {
            tooltip_unit_suffix(unit, sizeof(unit), group.unit);
        }
        tooltip_row(req, STR_LIT("magnitude"), "%.3f%s", vec3_length(vec), unit);
        tooltip_row(req, STR_LIT("vector"), "%.3f, %.3f, %.3f%s", vec.x, vec.y, vec.z, unit);
    }
}

static void fill_picking_tooltip(PickingTooltipTextRequest* req, const ApplicationState& state, const PickingHit& hit) {
    ASSERT(req);
    const md_system_t& sys = state.mold.sys;

    if (hit.domain == PickingDomain_Atom) {
        if (hit.local_idx < sys.atom.count) tooltip_atom(req, state, hit.local_idx);
    } else if (hit.domain == PickingDomain_Bond) {
        if (hit.local_idx < sys.bond.count) tooltip_bond(req, state, hit.local_idx);
    } else if (hit.domain == PickingDomain_BackboneSegment) {
        if (hit.local_idx < sys.protein_backbone.segment.count && sys.protein_backbone.segment.comp_idx) tooltip_backbone_segment(req, state, hit.local_idx);
    } else if (hit.domain == PickingDomain_Dipole) {
        tooltip_dipole(req, state, hit);
    }
}

static void tooltip_text(str_t str) {
    if (str_empty(str)) {
        ImGui::TextUnformatted("");
    } else {
        ImGui::TextUnformatted(str.ptr, str.ptr + str.len);
    }
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state) {
    if (hit.raw_idx == INVALID_PICKING_IDX) return;
    
	md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
    defer { md_temp_end(temp); };

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
    task_system::pool_interrupt_running_tasks();

    if (state->script.eval) md_script_eval_interrupt(state->script.eval);

    task_system::pool_wait_for_completion();
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
    size_t capacity = ALIGN_TO(num_atoms, 16);
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
    ASSERT(out_slot_idx);
#if FRAME_CACHE_SIZE == 4
    md_128i frame_indices = md_mm_loadu_epi32(cache->frame_idx);
    md_128i cmp_mask = md_mm_cmpeq_epi32(frame_indices, md_mm_set1_epi32((int32_t)frame_idx));
    int mask = md_mm_movemask_epi8(cmp_mask);
    *out_slot_idx = ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
    return mask != 0;
#elif FRAME_CACHE_SIZE == 8
    md_256i frame_indices = md_mm256_loadu_epi32(cache->frame_idx);
    md_256i cmp_mask = md_mm256_cmpeq_epi32(frame_indices, md_mm256_set1_epi32((int32_t)frame_idx));
    int mask = md_mm256_movemask_epi8(cmp_mask);
    *out_slot_idx = ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
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

void clear_system_frame_cache(ApplicationState* state) {
    ASSERT(state);
    clear_frame_cache(&state->mold.frame_cache);
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
static void upload_secondary_structure_weights(ApplicationState* data) {
    const md_system_t& sys = data->mold.sys;
    const size_t num_segments = sys.protein_backbone.segment.count;
    const md_secondary_structure_t* ss = md_util_state_secondary_structure(&data->mold.state, &sys);
    if (!ss || num_segments == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_gl_secondary_structure_t* weights = md_temp_alloc_array(temp, md_gl_secondary_structure_t, num_segments);
    for (size_t i = 0; i < num_segments; ++i) {
        weights[i] = md_gl_secondary_structure_convert(ss[i]);
    }
    secondary_structure_weights_fill_isolated_coils(weights, &sys.protein_backbone);
    md_gl_mol_set_backbone_secondary_structure(data->mold.gl_mol, 0, (uint32_t)num_segments, weights, 0);
}

static void secondary_structure_render_denoise(md_secondary_structure_t* dst, const md_secondary_structure_t* src, size_t num_frames, size_t stride) {
    ASSERT(dst);
    ASSERT(src);

    if (num_frames == 0 || stride == 0) {
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    defer { md_temp_end(temp_scope); };

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
            size_t end = MIN(frame_idx + window_radius + 1, num_frames);
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

str_t run_attribute_path(char* buf, size_t cap, const ApplicationState* app, str_t leaf) {
    ASSERT(buf && cap > 0);
    buf[0] = '\0';
    if (app->mold.run[0] == '\0') {
        return {};
    }
    const int len = snprintf(buf, cap, "%s/" STR_FMT, app->mold.run, STR_ARG(leaf));
    if (len <= 0 || (size_t)len >= cap) {
        buf[0] = '\0';
        return {};
    }
    return {buf, (size_t)len};
}

const md_attribute_t* run_attribute(const ApplicationState* app, str_t leaf) {
    ASSERT(app);
    if (app->mold.run[0] == '\0') return nullptr;
    return md_attributes_find_in(&app->mold.sys.attributes, str_from_cstr(app->mold.run), leaf);
}

const md_attribute_t* run_time_axis(const ApplicationState* app) {
    const md_attribute_t* axis = run_attribute(app, STR_LIT("time"));
    return md_attribute_view(axis, MD_ATTRIBUTE_TYPE_F64, 1, 1) ? axis : nullptr;
}

size_t run_num_frames(const ApplicationState* app) {
    const md_attribute_t* axis = run_time_axis(app);
    return axis ? axis->format.shape[0] : 0;
}

const double* run_frame_times(const ApplicationState* app) {
    return (const double*)md_attribute_view(run_time_axis(app), MD_ATTRIBUTE_TYPE_F64, 1, 1);
}

md_unit_t run_time_unit(const ApplicationState* app) {
    const md_attribute_t* axis = run_time_axis(app);
    return axis ? axis->unit : md_unit_none();
}

const md_secondary_structure_t* displayed_secondary_structure(const ApplicationState* app) {
    ASSERT(app);
    const md_system_state_t* state = &app->mold.state;
    const auto& render = app->trajectory_data.secondary_structure_render;
    if (render.data && render.stride == app->mold.sys.protein_backbone.segment.count && md_state_has_frame(state)) {
        const size_t frame = (size_t)md_state_frame_nearest(state);
        if ((frame + 1) * render.stride <= render.count) {
            return render.data + frame * render.stride;
        }
    }
    return md_util_state_secondary_structure(state, &app->mold.sys);
}

static const str_t frame_extract_paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };

bool extract_frame(const ApplicationState* app, int64_t frame, md_system_state_t* out) {
    ASSERT(app && out);
    md_system_extract_t* ex = md_system_extract_begin(&app->mold.sys, str_from_cstr(app->mold.run),
        frame_extract_paths, ARRAY_SIZE(frame_extract_paths), md_get_heap_allocator());
    if (!ex) {
        return false;
    }
    const bool ok = md_system_extract_frame(ex, frame, out);
    md_system_extract_end(ex);
    return ok;
}

static void end_frame_extracts(ApplicationState* app) {
    for (size_t i = 0; i < ARRAY_SIZE(app->mold.frame_extract); ++i) {
        md_system_extract_end(app->mold.frame_extract[i]);
        app->mold.frame_extract[i] = nullptr;
    }
}

// "run/<stem>" for a trajectory file. The stem is folded to letters, digits, '_' and '-' so that it
// is one path segment whatever the file is called, and it is taken from the file name alone so the
// same file gives the same run - and the same attribute ids - every time it is loaded.
static void run_path_from_file(char* buf, size_t cap, str_t path) {
    str_t file = path;
    extract_file(&file, path);
    size_t dot;
    if (str_rfind_char(&dot, file, '.') && dot > 0) {
        file = str_substr(file, 0, dot);
    }
    size_t len = (size_t)snprintf(buf, cap, "run/");
    for (size_t i = 0; i < file.len && len + 1 < cap; ++i) {
        const char c = file.ptr[i];
        const bool keep = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
        buf[len++] = keep ? c : '_';
    }
    if (len == 4) {
        len += (size_t)snprintf(buf + len, cap - len, "trajectory");
    }
    buf[MIN(len, cap - 1)] = '\0';
}

void free_trajectory_data(ApplicationState* state) {
    ASSERT(state);

    // The orientation reference was read from this run's first frame. A run loaded in its place, on the
    // same topology, has a first frame of its own.
    state->operations.initial_frame.valid = false;

    // Before anything they read goes: the contexts hold the run's files open.
    end_frame_extracts(state);

    state->files.trajectory[0] = '\0';

    md_array_free(state->timeline.x_values,  state->allocator.persistent);
    // Converted axes and histograms of what is about to go. The attribute versions they are keyed
    // on only mean something within one table, and the next one may be allocated where this was.
    series_cache_free(state);

    // backbone_angles.data and secondary_structure.data are VIEWS into sys.attributes and
    // are not md_arrays - freeing them here would hand the wrong pointer to the array allocator. The
    // table owns them and releases them with the system; dropping the views is all that is owed.
    state->trajectory_data.backbone_angles.data       = nullptr;
    state->trajectory_data.backbone_angles.stride     = 0;
    state->trajectory_data.backbone_angles.count      = 0;
    state->trajectory_data.secondary_structure.data   = nullptr;
    state->trajectory_data.secondary_structure.stride = 0;
    state->trajectory_data.secondary_structure.count  = 0;

    // The denoised render copy is still an ordinary array owned here.
    md_array_free(state->trajectory_data.secondary_structure_render.data,    state->allocator.persistent);

    // A stale match here would make interpolate_system_state silently skip pushing the new run's
    // first frame to the renderer whenever it happens to resolve to the same nearest frame index
    // as whatever was last displayed (frame 0 is the common case).
    state->mold.last_interpolated_nearest_frame = -1;
    state->mold.last_interpolated_frame = -1.0;

    // The views above are dropped first; now the storage they pointed at goes with the run - its
    // frame axis, the quantities derived per frame, and whatever was loaded against it.
    if (state->mold.run[0]) {
        md_attributes_remove_prefix(&state->mold.sys.attributes, str_from_cstr(state->mold.run));
        state->mold.run[0] = '\0';
    }

    free_frame_cache(&state->mold.frame_cache, state->allocator.persistent);
}

void init_trajectory_data(ApplicationState* data, uint32_t traj_flags) {
    // The trajectory is a RUN in the attribute table: "<run>/time" is its frame axis, and everything
    // sampled along it lives below the same prefix - the positions streamed from the file, the per
    // frame backbone data below, an energy file loaded later - so that freeing the trajectory is
    // removing the prefix. A trajectory that came inside the structure file (a multi model PDB) is
    // named after that; a structure of one frame publishes nothing and has no run.
    {
        const char* traj_file = data->files.trajectory[0] ? data->files.trajectory : data->files.molecule;
        run_path_from_file(data->mold.run, sizeof(data->mold.run), str_from_cstr(traj_file));
        const str_t run = str_from_cstr(data->mold.run);
        if (!loader::publish_run(&data->mold.sys, str_from_cstr(traj_file), run, traj_flags)) {
            md_attributes_remove_prefix(&data->mold.sys.attributes, run);
            data->mold.run[0] = '\0';
        }
    }

    size_t num_frames = run_num_frames(data);
    if (num_frames > 0) {
        size_t min_frame = 0;
        size_t max_frame = num_frames - 1;
        const double* frame_times = run_frame_times(data);
        char path_buf[256];

        init_frame_cache(&data->mold.frame_cache, data->mold.sys.atom.count, data->allocator.persistent);

        ASSERT(frame_times);

        // The timeline carries time in the unit the user asked to see it in. Everything downstream
        // reads x_values and view_range, so this is the one place the conversion happens; a later
        // change to the preference is picked up by update_timeline_time_unit in main.cpp.
        const double time_scl = display_units::factor(&data->timeline.time_unit, run_time_unit(data));
        data->timeline.time_scale    = time_scl;
        data->timeline.units_version = display_units::version();

        double min_time = frame_times[0] * time_scl;
        double max_time = frame_times[num_frames - 1] * time_scl;

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = (double)min_frame;
        data->timeline.filter.end_frame = (double)max_frame;

        md_array_resize(data->timeline.x_values, num_frames, data->allocator.persistent);
        for (size_t i = 0; i < num_frames; ++i) {
            data->timeline.x_values[i] = (float)(frame_times[i] * time_scl);
        }

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, (int64_t)max_frame);

        extract_frame(data, frame_idx, &data->mold.state);

        if (data->mold.sys.protein_backbone.segment.count > 0) {
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
            const size_t num_segments = data->mold.sys.protein_backbone.segment.count;
            md_attributes_t* attributes = &data->mold.sys.attributes;

            // Both are checked against "<run>/time", published above, so a declaration whose
            // outermost extent disagrees is refused at the point the mistake is made rather than
            // surviving until something reads past the end of it.

            // md_secondary_structure_t is a 4 byte enum, so I32 is the storage and the PATH is what
            // tells a consumer these are labels rather than numbers to average - which is the rule
            // md_system.h already states for integral attributes.
            const md_attribute_desc_t ss_desc = {
                .path   = run_attribute_path(path_buf, sizeof(path_buf), data, STR_LIT("backbone/secondary_structure")),
                .format = {
                    .type = MD_ATTRIBUTE_TYPE_I32, .components = 1,
                    .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
                },
                .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit   = md_unit_none(),
                .label  = STR_INIT("Secondary Structure"),
            };
            md_attribute_id_t ss_id = md_attributes_replace(attributes, &ss_desc);

            data->trajectory_data.secondary_structure.data   = (md_secondary_structure_t*)md_attributes_data(attributes, ss_id, MD_ATTRIBUTE_TYPE_I32);
            data->trajectory_data.secondary_structure.stride = num_segments;
            data->trajectory_data.secondary_structure.count  = num_segments * num_frames;

            const md_attribute_format_t angle_format = {
                .type = MD_ATTRIBUTE_TYPE_F32, .components = 2,   // phi, psi - one value, two parts
                .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
            };
            const md_attribute_desc_t angle_desc = {
                .path   = run_attribute_path(path_buf, sizeof(path_buf), data, STR_LIT("backbone/angle")),
                .format = angle_format,
                .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit   = md_unit_radian(),
                .label  = STR_INIT("Backbone Angles"),
            };
            md_attribute_id_t angle_id = md_attributes_replace(attributes, &angle_desc);

            // md_backbone_angles_t is exactly {float phi; float psi;}, so the table's storage IS a
            // md_backbone_angles_t array - no repacking, and the existing consumers keep their type.
            data->trajectory_data.backbone_angles.data   = (md_backbone_angles_t*)md_attributes_data(attributes, angle_id, MD_ATTRIBUTE_TYPE_F32);
            data->trajectory_data.backbone_angles.stride = num_segments;
            data->trajectory_data.backbone_angles.count  = num_segments * num_frames;

            // The denoised copy stays an ordinary array: it is a PRESENTATION smoothing of the one
            // above, read only by the renderer, and publishing it would put two answers to "what is
            // the secondary structure at frame f" in the same table.
            data->trajectory_data.secondary_structure_render.stride = data->mold.sys.protein_backbone.segment.count;
            data->trajectory_data.secondary_structure_render.count = data->mold.sys.protein_backbone.segment.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure_render.data, data->mold.sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(data->trajectory_data.secondary_structure_render.data, 0, md_array_bytes(data->trajectory_data.secondary_structure_render.data));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::create_pool_task(STR_LIT("Backbone Operations"), (uint32_t)num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                const md_system_t* sys = &data->mold.sys;

                md_temp_scope_t temp = md_temp_begin();
                defer { md_temp_end(temp); };

                md_system_state_t frame_state = { .alloc = temp.arena };
				md_system_state_init(&frame_state, sys->atom.count);

                // One context for the range, so the run's files stay open across its frames.
                md_system_extract_t* ex = md_system_extract_begin(sys, str_from_cstr(data->mold.run),
                    frame_extract_paths, ARRAY_SIZE(frame_extract_paths), md_get_heap_allocator());
                if (!ex) {
                    return;
                }
                for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                    md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                    md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;

                    if (!md_system_extract_frame(ex, frame_idx, &frame_state)) {
                        continue;
                    }
                    md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, frame_state.xyz, &frame_state.unitcell, &sys->protein_backbone);
                    md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, frame_state.xyz, &frame_state.unitcell, &sys->protein_backbone);
                }
                md_system_extract_end(ex);
            });

            uint64_t time = (uint64_t)md_tick_now();
            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data, t0 = time, num_frames]() {
                secondary_structure_render_denoise(
                    data->trajectory_data.secondary_structure_render.data,
                    data->trajectory_data.secondary_structure.data,
                    num_frames,
                    data->trajectory_data.secondary_structure.stride);

                uint64_t t1 = (uint64_t)md_tick_now();
                double elapsed = md_tick_to_seconds(t1 - t0);
                MD_LOG_INFO("Finished computing trajectory data (%.2fs)", elapsed);
                // The fill is finished, so the two temporal attributes have changed - which is the
                // one moment a consumer caching something derived from them needs to hear about.
                // The table carries that now; the hand rolled fingerprints beside the data are gone,
                // because a stamp that lives next to storage the table owns is a second answer
                // waiting to disagree with it.
                md_attributes_t* attributes = &data->mold.sys.attributes;
                char buf[256];
                md_attributes_touch(attributes, md_attributes_id_from_path(run_attribute_path(buf, sizeof(buf), data, STR_LIT("backbone/angle"))));
                md_attributes_touch(attributes, md_attributes_id_from_path(run_attribute_path(buf, sizeof(buf), data, STR_LIT("backbone/secondary_structure"))));
                data->trajectory_data.secondary_structure_render.fingerprint = generate_fingerprint();

				// The data just went from placeholder to real values for every frame - force the
				// next interpolate_system_state to actually push it, even if playback is paused on
				// the same nearest frame it was on before this task started.
				data->mold.last_interpolated_nearest_frame = -1;
				data->mold.interpolate_system_state = true;
                data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                flag_all_representations_as_dirty(data);
            });

            task_system::set_task_dependency(main_task, data->tasks.backbone_computations);
            task_system::enqueue_task(data->tasks.backbone_computations);
        }

        data->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
        data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;

        // Prefetch frames
        //launch_prefetch_job(data);
    }
}

void init_system_data(ApplicationState* data) {
    if (data->mold.sys.atom.count) {
        md_bitfield_clear(&data->operations.recenter_query.mask);
        md_bitfield_clear(&data->operations.selection_mask);
        data->operations.recenter_query.valid = false;
        data->operations.recenter_query.dynamic = false;
        data->operations.recenter_query.evaluated_version = 0;
        data->operations.recenter_query.ir_fingerprint = 0;
        data->operations.initial_frame.valid = false;
        data->operations.state_rotation = mat4_ident();
        recenter_mark_query_dirty(data);

        data->mold.gl_mol = md_gl_mol_create(&data->mold.sys);

        // The backbone of the loaded coordinates, carried by the state, which stands until a run replaces it frame by
        // frame (interpolate_system_state)
        if (md_util_state_backbone_compute(&data->mold.state, &data->mold.sys)) {
            upload_secondary_structure_weights(data);
        }

        mat3_t A;
		md_unitcell_A_extract_float(A.elem, &data->mold.state.unitcell);
		vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
        data->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

        vec3_t aabb_min = {};
        vec3_t aabb_max = {};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, data->mold.state.xyz, nullptr, nullptr, data->mold.state.num_atoms);

        const vec3_t cell_ext = mat3_mul_vec3(A, vec3_set1(1.0f));
        const float max_cell_ext = vec3_reduce_max(cell_ext);
        const float max_aabb_ext = vec3_reduce_max(vec3_sub(aabb_max, aabb_min));

        // Calculate a default view transform to use later as a reset target
        ViewTransform default_view = {};
        reset_view(&default_view, data->mold.state);

        data->view.camera.near_plane = 1.0f;
        data->view.camera.far_plane = 100000.0f;
        data->view.trackball_param.max_distance = MAX(max_cell_ext, max_aabb_ext) * 10.0f;

        data->view.target = default_view;
        data->view.camera = default_view;
        

#if EXPERIMENTAL_GFX_API
        const md_system_t& mol = data->mold.sys;
        vec3_t& aabb_min = data->mold.sys_aabb_min;
        vec3_t& aabb_max = data->mold.sys_aabb_max;
        md_util_compute_aabb_soa(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

        data->mold.gfx_structure = md_gfx_structure_create(mol.atom.count, mol.covalent.count, mol.backbone.count, mol.backbone.range_count, mol.residue.count, mol.instance.count);
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.radius, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.sys_aabb_min, &data->mold.sys_aabb_max);
        if (mol.instance.count > 0) {
            md_gfx_structure_set_instance_atom_ranges(data->mold.gfx_structure, 0, mol.instance.count, (md_gfx_range_t*)mol.instance.atom_range, 0);
            md_gfx_structure_set_instance_transforms(data->mold.gfx_structure, 0, mol.instance.count, mol.instance.transform, 0);
        }
#endif
    }
    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemInit, viamd::EventPayloadType_ApplicationState, data);

    init_all_representations(data);
    data->script.compile_ir = true;
}

void free_system_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    // The arena OUTLIVES the system it backs, and the zeroing below is why there used to be a
    // second handle to it: MEMSET over md_system_t clears sys.alloc along with everything else, so
    // the only pointer to the arena has to be carried across by hand. Rewound, not destroyed - the
    // next load reuses it.
    md_allocator_i* sys_alloc = data->mold.sys.alloc;
    md_arena_allocator_reset(sys_alloc);
    MEMSET(&data->mold.sys, 0, sizeof(data->mold.sys));
    MEMSET(&data->mold.state, 0, sizeof(data->mold.state));
    data->mold.sys.alloc   = sys_alloc;
    data->mold.state.alloc = sys_alloc;

    md_array_free(data->operations.initial_frame.rel_xyzw, data->allocator.persistent);
    data->operations.initial_frame.rel_xyzw = nullptr;
    data->operations.initial_frame.valid = false;

    md_gl_mol_destroy(data->mold.gl_mol);

    // The dataset's GPU data goes with the dataset, for the same reason gl_mol does.
    system_gpu_data_free(data);

    MEMSET(data->files.molecule, 0, sizeof(data->files.molecule));


    MEMSET(data->mold.frame_cache.states, 0, sizeof(data->mold.frame_cache.states));
    clear_frame_cache(&data->mold.frame_cache);

    md_bitfield_clear(&data->selection.selection_mask);
    md_bitfield_clear(&data->selection.highlight_mask);
    // Computed for this system
    script_vis_reset(data);

    // The tasks were waited for above: with the fields cleared, nothing uses any IR
    data->script.ir = nullptr;
    data->script.eval_ir = nullptr;
    script_ir_collect(data);
    if (data->script.eval) {
        md_script_eval_free(data->script.eval);
        data->script.eval = nullptr;
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemFree, viamd::EventPayloadType_ApplicationState, data);
}

bool load_data_from_file(ApplicationState* state, str_t filepath, const loader::LoaderState& load_state) {
    ASSERT(state);

    bool success = false;
    str_t path_to_file = md_path_make_canonical(filepath, state->allocator.frame);
    if (path_to_file) {
        if (load_state.flags & LoaderFlag_Supplemental) {
            // Neither the system nor the trajectory is touched - this file only adds to what is
            // already loaded, and the loader merges it into the system's attribute table. A
            // component picks it up from there, off the ViamdLoadData event below.
            //
            // 'success' stays false on purpose: it is what tells the caller a system was loaded,
            // and it resets the camera and the animation when it is true.
            const str_t run = str_from_cstr(state->mold.run);
            if ((load_state.flags & LoaderFlag_Temporal) && str_empty(run)) {
                VIAMD_LOG_ERROR("'" STR_FMT "' holds data along a trajectory; load the trajectory first", STR_ARG(path_to_file));
                return false;
            }

            // The table is about to be written to, and worker threads read it: script evaluation
            // (attr() and its frames), the backbone task and playback's frame loads, all through
            // extraction contexts. The evaluations are interrupted, since the script is recompiled
            // below and they would start over anyway; the rest are left to FINISH rather than
            // interrupted - an interrupted backbone task leaves a just loaded trajectory without its
            // backbone data.
            if (state->script.eval) md_script_eval_interrupt(state->script.eval);
            task_system::pool_wait_for_completion();

            if (loader::load_supplemental(&state->mold.sys, path_to_file, load_state, run)) {
                // A topology replaces bonds (and with them structures), which the GPU holds a copy of
                state->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
                // attr() in the script resolves against the table when it is compiled.
                state->script.compile_ir = true;
                if (load_state.flags & LoaderFlag_Temporal) {
                    const char* hint = (load_state.type == LoaderType_EDR) ? "edr/<term>" :
                                       (load_state.type == LoaderType_XVG) ? "xvg/<file>/<legend>" : "csv/<file>/<column>";
                    VIAMD_LOG_SUCCESS("Loaded '" STR_FMT "' into '" STR_FMT "'; plot it from System > Series, or read it in the script with attr(\"%s\")", STR_ARG(path_to_file), STR_ARG(run), hint);
                } else {
                    VIAMD_LOG_SUCCESS("Successfully loaded supplemental data from file '" STR_FMT "'", STR_ARG(path_to_file));
                }
            } else {
                VIAMD_LOG_ERROR("Failed to load supplemental data from file '" STR_FMT "'", STR_ARG(path_to_file));
            }
        } else if (load_state.flags & LoaderFlag_System) {
            interrupt_async_tasks(state);
            free_trajectory_data(state);
            free_system_data(state);

            // sys.alloc and state.alloc are already the dataset arena - set at startup and put
            // back by free_system_data above. The loaded coordinates are the initial current state;
            // inference copies them into sys.reference, so the two stay independent.
            ASSERT(state->mold.sys.alloc && state->mold.state.alloc == state->mold.sys.alloc);
            if (!loader::load(&state->mold.sys, &state->mold.state, path_to_file, load_state)) {
                VIAMD_LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            success = true;
            VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            str_copy_to_char_buf(state->files.molecule, sizeof(state->files.molecule), path_to_file);
            state->files.coarse_grained = load_state.flags & LoaderFlag_CoarseGrained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_infer_flags_t flags = state->files.coarse_grained ? MD_UTIL_INFER_NONE : MD_UTIL_INFER_ALL;
            if ((load_state.flags & LoaderFlag_Topology) && state->mold.sys.bond.count > 0) {
                // The file's bonds are the force field's; inferring would replace them with a guess.
                // Structures and rings are still derived, from those bonds. A topology format that
                // carried no bonds at all (an H5MD file without connectivity) has none to protect, and
                // gets them inferred like any other structure.
                flags &= ~MD_UTIL_INFER_BOND_BIT;
                flags |= MD_UTIL_INFER_STRUCTURE_BIT;
            }
            md_util_system_infer(&state->mold.sys, &state->mold.state, flags);
            init_system_data(state);

            // A structure file of several frames (a multi model PDB, an XYZ trajectory) is its own run.
            init_trajectory_data(state, (load_state.flags & LoaderFlag_DisableCacheWrite) ? MD_RUN_FLAG_DISABLE_CACHE_WRITE : 0);
        } else if (load_state.flags & LoaderFlag_Trajectory) {
            if (!state->mold.sys.atom.count) {
                VIAMD_LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }
            interrupt_async_tasks(state);
            free_trajectory_data(state);
            state->animation.frame = 0;

            // Publishing the run IS opening the trajectory: the run is named after the file, and its
            // positions are read from it frame by frame from here on.
            str_copy_to_char_buf(state->files.trajectory, sizeof(state->files.trajectory), path_to_file);
            init_trajectory_data(state, (load_state.flags & LoaderFlag_DisableCacheWrite) ? MD_RUN_FLAG_DISABLE_CACHE_WRITE : 0);
            success = run_num_frames(state) > 0;
            if (success) {
                VIAMD_LOG_SUCCESS("Successfully opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            } else {
                state->files.trajectory[0] = '\0';
                VIAMD_LOG_ERROR("Failed to open trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            }
        }
#if MD_ENABLE_GPU
        // The dataset's uploaded GTO basis, built here rather than by whichever component happens to
        // want it first: it is derived from the system's own basis/ attributes, every consumer of
        // that system wants the same one, and a system's GPU basis must not depend on which UI
        // component is compiled in. Returns false when the system publishes no basis, which is the
        // normal case rather than a failure.
        system_gpu_data_update(state, DEFAULT_GTO_CUTOFF_VALUE);
#endif

        LoadDataPayload data = {
            .app_state = state,
            .loader_state = load_state,
            .path_to_file = path_to_file,
        };
        viamd::event_system_broadcast_event(viamd::EventType_ViamdLoadData, viamd::EventPayloadType_LoadData, &data);
    }

    return success;
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

struct WorkspaceWindow {
    char  name[32];
    bool* show;
};
static WorkspaceWindow workspace_windows[64];
static size_t num_workspace_windows = 0;

void workspace_register_window(const char* name, bool* show) {
    ASSERT(name && show);
    for (size_t i = 0; i < num_workspace_windows; ++i) {
        if (strcmp(workspace_windows[i].name, name) == 0) {
            workspace_windows[i].show = show;
            return;
        }
    }
    if (num_workspace_windows < ARRAY_SIZE(workspace_windows)) {
        WorkspaceWindow& w = workspace_windows[num_workspace_windows++];
        snprintf(w.name, sizeof(w.name), "%s", name);
        w.show = show;
    } else {
        MD_LOG_ERROR("Workspace: too many windows registered, '%s' is not stored", name);
    }
}

// Frame <-> time along the timeline, for the parts of it stored as frames (see TimelineView)
static double workspace_frame_to_time(const ApplicationState* app, double frame) {
    const int64_t n = (int64_t)md_array_size(app->timeline.x_values);
    if (n == 0) return frame;
    const int64_t f0 = CLAMP((int64_t)frame, (int64_t)0, n - 1);
    const int64_t f1 = CLAMP(f0 + 1, (int64_t)0, n - 1);
    const double t = CLAMP(frame - (double)f0, 0.0, 1.0);
    return lerp((double)app->timeline.x_values[f0], (double)app->timeline.x_values[f1], t);
}

static double workspace_time_to_frame(const ApplicationState* app, double time) {
    const int64_t n = (int64_t)md_array_size(app->timeline.x_values);
    if (n == 0) return time;
    const float* x = app->timeline.x_values;
    if (time <= x[0]) return 0.0;
    if (time >= x[n - 1]) return (double)(n - 1);
    int64_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        const int64_t mid = (lo + hi) / 2;
        if (x[mid] <= time) lo = mid; else hi = mid;
    }
    const double dx = (double)x[hi] - (double)x[lo];
    return (double)lo + (dx > 0.0 ? (time - x[lo]) / dx : 0.0);
}

// What is read before the data it refers to is loaded, applied once it is
struct WorkspacePending {
    str_t molecule_file;
    str_t trajectory_file;
    str_t energy_file;
    md_array(str_t) series_files;
    bool  coarse_grained;

    bool   has_frame;
    double frame;

    // Kept aside rather than written to the camera: loading the molecule puts the camera at the
    // default view of what was loaded (init_system_data), which would overwrite it
    bool          has_camera;
    ViewTransform camera;

    md_array(md_atom_pair_t) user_bonds;

    bool has_selection_mask;
    md_bitfield_t selection_mask;

    bool has_recenter_target;
    md_bitfield_t recenter_target;

    // Timeline filter and zoom, in frames: a time would depend on the unit it is shown in
    bool   has_filter_range;
    double filter_beg_frame;
    double filter_end_frame;
    bool   has_view_range;
    double view_beg_frame;
    double view_end_frame;
};

// A file named in the workspace: relative to it, unless no relative path existed when it was saved
static str_t workspace_file_path(str_t folder, str_t arg, md_allocator_i* alloc) {
    str_t file;
    viamd::extract_str(file, arg);
    if (str_empty(file)) return {};
    md_strb_t path = md_strb_create(alloc);
    if (!md_path_is_absolute(file)) {
        path += folder;
    }
    path += file;
    return md_path_make_canonical(path, alloc);
}

static void workspace_reset(ApplicationState* data) {
    remove_all_selections(data);
    remove_all_representations(data);
    data->editor.SetText("");
    data->files.workspace[0] = '\0';

    // The fields rather than the struct: whether the window is open is decided by [Windows]
    data->animation.frame = 0.0;
    data->animation.fps = 10.0f;
    data->animation.tension = 0.0f;
    data->animation.interpolation = InterpolationMode::CubicSpline;
    data->animation.mode = PlaybackMode::Stopped;

    plot_clear(data->timeline.subplots, PLOT_MAX_SUBPLOTS);
    plot_clear(data->distributions.subplots, PLOT_MAX_SUBPLOTS);
    data->timeline.num_subplots = 1;
    data->distributions.num_subplots = 1;
    data->timeline.filter.enabled = false;
    data->timeline.filter.temporal_window.enabled = false;
    data->timeline.filter.temporal_window.extent_in_frames = 10;

    data->visuals = {};
    data->simulation_box = {};
    data->view.mode = CameraMode::Perspective;
    data->view.camera.fov_y = Camera{}.fov_y;

    data->selection.granularity = SelectionGranularity::Atom;
    md_bitfield_clear(&data->selection.selection_mask);
    single_selection_sequence_clear(&data->selection.single_selection_sequence);

    data->operations.recenter = false;
    data->operations.fixate_orientation = false;
    data->operations.apply_pbc = false;
    data->operations.unwrap_structures = false;
    data->operations.recalc_bonds = false;
    data->operations.recenter_query.enabled = false;
    data->operations.recenter_query.query[0] = '\0';
    recenter_mark_query_dirty(data);
    md_bitfield_clear(&data->operations.selection_mask);
}

static void deserialize_files(viamd::deserialization_state_t& state, WorkspacePending& pending, str_t folder, md_allocator_i* alloc) {
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        if (str_eq(ident, STR_LIT("MoleculeFile"))) {
            pending.molecule_file = workspace_file_path(folder, arg, alloc);
        } else if (str_eq(ident, STR_LIT("TrajectoryFile"))) {
            pending.trajectory_file = workspace_file_path(folder, arg, alloc);
        } else if (str_eq(ident, STR_LIT("EnergyFile"))) {
            pending.energy_file = workspace_file_path(folder, arg, alloc);
        } else if (str_eq(ident, STR_LIT("SeriesFile"))) {
            const str_t path = workspace_file_path(folder, arg, alloc);
            if (!str_empty(path)) md_array_push(pending.series_files, path, alloc);
        } else if (str_eq(ident, STR_LIT("CoarseGrained"))) {
            viamd::extract_bool(pending.coarse_grained, arg);
        }
    }
}

static void deserialize_representation(ApplicationState* data, viamd::deserialization_state_t& state) {
    Representation* rep = create_representation(data);
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        if (str_eq(ident, STR_LIT("Name"))) {
            viamd::extract_to_char_buf(rep->name, sizeof(rep->name), arg);
        } else if (str_eq(ident, STR_LIT("Filter"))) {
            viamd::extract_to_char_buf(rep->filt, sizeof(rep->filt), arg);
        } else if (str_eq(ident, STR_LIT("Enabled"))) {
            viamd::extract_bool(rep->enabled, arg);
        } else if (str_eq(ident, STR_LIT("Type"))) {
            viamd::extract_enum(rep->type, arg, (int)RepresentationType::Count);
        } else if (str_eq(ident, STR_LIT("ColorMapping"))) {
            viamd::extract_enum(rep->color_mapping, arg, (int)ColorMapping::Count);
        } else if (str_eq(ident, STR_LIT("StaticColor")) || str_eq(ident, STR_LIT("BaseColor"))) {
            viamd::extract_vec4(rep->base_color, arg);
        } else if (str_eq(ident, STR_LIT("Saturation"))) {
            viamd::extract_flt(rep->saturation, arg);
        } else if (str_eq(ident, STR_LIT("TintColor"))) {
            viamd::extract_vec4(rep->tint_color, arg);
        } else if (str_eq(ident, STR_LIT("TintScale"))) {
            viamd::extract_flt(rep->tint_scale, arg);
        } else if (str_eq(ident, STR_LIT("SecondaryStructureColorUnknown"))) {
            viamd::extract_vec4(rep->secondary_structure.color_unknown, arg);
        } else if (str_eq(ident, STR_LIT("SecondaryStructureColorCoil"))) {
            viamd::extract_vec4(rep->secondary_structure.color_coil, arg);
        } else if (str_eq(ident, STR_LIT("SecondaryStructureColorHelix"))) {
            viamd::extract_vec4(rep->secondary_structure.color_helix, arg);
        } else if (str_eq(ident, STR_LIT("SecondaryStructureColorSheet"))) {
            viamd::extract_vec4(rep->secondary_structure.color_sheet, arg);
        } else if (str_eq(ident, STR_LIT("BondColor"))) {
            viamd::extract_enum(rep->bond_color, arg, (int)BondColorMode::Count);
        } else if (str_eq(ident, STR_LIT("BondSharpness"))) {
            viamd::extract_flt(rep->bond_sharpness, arg);
        } else if (str_eq(ident, STR_LIT("BondBaseColor"))) {
            viamd::extract_vec4(rep->bond_base_color, arg);
        } else if (str_eq(ident, STR_LIT("Radius"))) {
            // DEPRECATED
            viamd::extract_flt(rep->scale.x, arg);
        } else if (str_eq(ident, STR_LIT("Tension"))) {
            // DEPRECATED
        } else if (str_eq(ident, STR_LIT("Width"))) {
            viamd::extract_flt(rep->scale.x, arg);
        } else if (str_eq(ident, STR_LIT("Thickness"))) {
            viamd::extract_flt(rep->scale.y, arg);
        } else if (str_eq(ident, STR_LIT("Param"))) {
            viamd::extract_vec4(rep->scale, arg);
        } else if (str_eq(ident, STR_LIT("DynamicEval"))) {
            viamd::extract_bool(rep->dynamic_evaluation, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyPath"))) {
            // Ids are a function of the path: this resolves before the dataset is loaded
            rep->atom_attribute.key = md_attributes_id_from_path(arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyVariant"))) {
            viamd::extract_int(rep->atom_attribute.variant_idx, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyColormap"))) {
            viamd::extract_int(rep->atom_attribute.scale.colormap, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyRange"))) {
            // A range that was written is one somebody set: a workspace from before the range could
            // follow the values has no AtomicPropertyAutoRange, and keeps it. One that has, has it
            // after this line and says for itself.
            float r[2];
            if (viamd::extract_flt_vec(r, 2, arg)) {
                rep->atom_attribute.scale.range_beg = r[0];
                rep->atom_attribute.scale.range_end = r[1];
                rep->atom_attribute.scale.auto_range = false;
            }
        } else if (str_eq(ident, STR_LIT("AtomicPropertyAutoRange"))) {
            viamd::extract_bool(rep->atom_attribute.scale.auto_range, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertySymmetric"))) {
            viamd::extract_bool(rep->atom_attribute.scale.symmetric, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyLegend"))) {
            viamd::extract_bool(rep->atom_attribute.scale.show_legend, arg);
        } else if (str_eq(ident, STR_LIT("AtomicPropertyDataRange")) || str_eq(ident, STR_LIT("AtomicPropertySymmetricZero"))) {
            // Written by earlier versions: the span is measured from the data now, and the old
            // symmetric flag never reached the colours
        } else if (str_eq(ident, STR_LIT("DipolePath"))) {
            rep->dipole.dipole_key = md_attributes_id_from_path(arg);
        } else if (str_eq(ident, STR_LIT("DipoleIndex"))) {
            int i;
            if (viamd::extract_int(i, arg) && i >= 0) rep->dipole.dipole_index = (uint32_t)i;
        } else if (str_eq(ident, STR_LIT("DipoleColor"))) {
            viamd::extract_vec4(rep->dipole.color, arg);
        } else if (str_eq(ident, STR_LIT("DipoleOffset"))) {
            viamd::extract_vec3(rep->dipole.offset, arg);
        } else if (str_eq(ident, STR_LIT("DipoleScale"))) {
            viamd::extract_dbl(rep->dipole.scale, arg);
        } else if (str_eq(ident, STR_LIT("DipoleRadius"))) {
            viamd::extract_flt(rep->dipole.radius, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureMoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureOrbitalIdx"))) {
            viamd::extract_int(rep->electronic_structure.orbital_idx, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureExcitedStateIdx"))) {
            viamd::extract_int(rep->electronic_structure.excited_state_idx, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoLambdaIdx"))) {
            viamd::extract_int(rep->electronic_structure.nto_lambda_idx, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyPath"))) {
            // The id is a function of the path, so this resolves without the dataset being
            // loaded yet - and names the same property after a reload, which the index it
            // replaces did not.
            rep->electronic_structure.density_property_key = md_attributes_id_from_path(arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIdx"))) {
            // DEPRECATED. A position in a list that only existed while a particular file was
            // open; there is nothing to resolve it against. Workspaces written before the
            // path was stored fall back to the first property available.
        } else if (str_eq(ident, STR_LIT("ElectronicStructureRes"))) {
            int res;
            if (viamd::extract_int(res, arg)) rep->electronic_structure.resolution = (VolumeResolution)res;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureType"))) {
            int type;
            if (viamd::extract_int(type, arg)) electronic_structure_set_legacy_type(&rep->electronic_structure, type);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureSource"))) {
            int source;
            if (viamd::extract_int(source, arg)) {
                rep->electronic_structure.source = (ElectronicStructureSource)source;
                electronic_structure_set_source_defaults(&rep->electronic_structure);
            }
        } else if (str_eq(ident, STR_LIT("ElectronicStructureField"))) {
            int field;
            if (viamd::extract_int(field, arg)) electronic_structure_set_legacy_field(&rep->electronic_structure, field);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureUseMagnitude"))) {
            viamd::extract_bool(rep->electronic_structure.use_magnitude, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureSpin"))) {
            int spin;
            if (viamd::extract_int(spin, arg)) rep->electronic_structure.spin = (ElectronicStructureSpin)spin;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoComponent"))) {
            int component;
            if (viamd::extract_int(component, arg)) rep->electronic_structure.nto_component = (ElectronicStructureNtoComponent)component;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureTransitionDensityComponent"))) {
            int component;
            if (viamd::extract_int(component, arg)) rep->electronic_structure.transition_density_component = (ElectronicStructureTransitionDensityComponent)component;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureIso"))) {
            viamd::extract_dbl(rep->electronic_structure.iso_value, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColoring"))) {
            int c;
            if (viamd::extract_int(c, arg) && c >= 0 && c < (int)SurfaceColoring::Count) rep->electronic_structure.coloring = (SurfaceColoring)c;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldKind"))) {
            int k;
            if (viamd::extract_int(k, arg) && k >= 0 && k < (int)SurfaceFieldKind::Count) rep->electronic_structure.field_kind = (SurfaceFieldKind)k;
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldColormap"))) {
            viamd::extract_int(rep->electronic_structure.field_map.colormap, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldRange"))) {
            float r[2];
            if (viamd::extract_flt_vec(r, 2, arg)) { rep->electronic_structure.field_map.range_beg = r[0]; rep->electronic_structure.field_map.range_end = r[1]; }
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldSymmetric"))) {
            viamd::extract_bool(rep->electronic_structure.field_map.symmetric, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldAutoRange"))) {
            viamd::extract_bool(rep->electronic_structure.field_map.auto_range, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldLegend"))) {
            viamd::extract_bool(rep->electronic_structure.field_map.show_legend, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureTintPos"))) {
            viamd::extract_vec4(rep->electronic_structure.tint_psi_pos, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureTintNeg"))) {
            viamd::extract_vec4(rep->electronic_structure.tint_psi_neg, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureTintDen"))) {
            viamd::extract_vec4(rep->electronic_structure.tint_den, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColPos"))) {
            viamd::extract_vec4(rep->electronic_structure.col_psi_pos, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColNeg"))) {
            viamd::extract_vec4(rep->electronic_structure.col_psi_neg, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColDen"))) {
            viamd::extract_vec4(rep->electronic_structure.col_den, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColAtt"))) {
            viamd::extract_vec4(rep->electronic_structure.col_att, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureColDet"))) {
            viamd::extract_vec4(rep->electronic_structure.col_det, arg);
        } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIsoCount"))) {
            int n;
            if (viamd::extract_int(n, arg)) {
                rep->electronic_structure.density_property.num_isos = CLAMP(n, 0, (int)ARRAY_SIZE(rep->electronic_structure.density_property.values));
            }
        } else {
            for (int i = 0; i < (int)ARRAY_SIZE(rep->electronic_structure.density_property.values); ++i) {
                char key[64];
                snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", i);
                if (str_eq(ident, str_from_cstr(key))) {
                    viamd::extract_dbl(rep->electronic_structure.density_property.values[i], arg);
                    break;
                }
                snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", i);
                if (str_eq(ident, str_from_cstr(key))) {
                    viamd::extract_vec4(rep->electronic_structure.density_property.colors[i], arg);
                    break;
                }
            }
        }
    }
}

static void serialize_representation(viamd::serialization_state_t& state, const ApplicationState* app_state, const Representation& rep) {
    const md_attributes_t* attributes = &app_state->mold.sys.attributes;

    viamd::write_section_header(state, STR_LIT("Representation"));
    viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
    viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filt));
    viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
    viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
    viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
    viamd::write_vec4(state, STR_LIT("BaseColor"), rep.base_color);
    viamd::write_flt(state,  STR_LIT("Saturation"), rep.saturation);
    viamd::write_vec4(state, STR_LIT("TintColor"), rep.tint_color);
    viamd::write_flt(state,  STR_LIT("TintScale"), rep.tint_scale);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorUnknown"), rep.secondary_structure.color_unknown);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorCoil"),    rep.secondary_structure.color_coil);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorHelix"),   rep.secondary_structure.color_helix);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorSheet"),   rep.secondary_structure.color_sheet);
    viamd::write_int(state,  STR_LIT("BondColor"), (int)rep.bond_color);
    viamd::write_flt(state,  STR_LIT("BondSharpness"), rep.bond_sharpness);
    viamd::write_vec4(state, STR_LIT("BondBaseColor"), rep.bond_base_color);

    viamd::write_vec4(state, STR_LIT("Param"), rep.scale);
    viamd::write_bool(state, STR_LIT("DynamicEval"), rep.dynamic_evaluation);

    // Attributes by path: an id is a hash of it, and the path is what a reader can see
    if (const md_attribute_t* prop = md_attributes_get(attributes, rep.atom_attribute.key)) {
        viamd::write_str(state, STR_LIT("AtomicPropertyPath"), prop->path);
        viamd::write_int(state, STR_LIT("AtomicPropertyVariant"), rep.atom_attribute.variant_idx);
        const ColorScale& scale = rep.atom_attribute.scale;
        const float range[2] = { scale.range_beg, scale.range_end };
        viamd::write_int(state, STR_LIT("AtomicPropertyColormap"), scale.colormap);
        viamd::write_flt_vec(state, STR_LIT("AtomicPropertyRange"), range, 2);
        // After the range, which on its own reads as a range somebody set
        viamd::write_bool(state, STR_LIT("AtomicPropertyAutoRange"), scale.auto_range);
        viamd::write_bool(state, STR_LIT("AtomicPropertySymmetric"), scale.symmetric);
        viamd::write_bool(state, STR_LIT("AtomicPropertyLegend"), scale.show_legend);
    }

    if (rep.type == RepresentationType::DipoleMoment) {
        if (const md_attribute_t* dip = md_attributes_get(attributes, rep.dipole.dipole_key)) {
            viamd::write_str(state, STR_LIT("DipolePath"), dip->path);
        }
        viamd::write_int(state,  STR_LIT("DipoleIndex"), (int64_t)rep.dipole.dipole_index);
        viamd::write_vec4(state, STR_LIT("DipoleColor"), rep.dipole.color);
        viamd::write_vec3(state, STR_LIT("DipoleOffset"), rep.dipole.offset);
        viamd::write_dbl(state,  STR_LIT("DipoleScale"), rep.dipole.scale);
        viamd::write_flt(state,  STR_LIT("DipoleRadius"), rep.dipole.radius);
    }

    if (rep.type == RepresentationType::ElectronicStructure) {
        viamd::write_int(state,  STR_LIT("ElectronicStructureMoIdx"),    rep.electronic_structure.orbital_idx);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoIdx"),   rep.electronic_structure.excited_state_idx);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoLambdaIdx"), rep.electronic_structure.nto_lambda_idx);
        if (const md_attribute_t* prop = md_attributes_get(attributes, rep.electronic_structure.density_property_key)) {
            viamd::write_str(state, STR_LIT("ElectronicStructureDensityPropertyPath"), prop->path);
        }
        viamd::write_int(state,  STR_LIT("ElectronicStructureType"),     electronic_structure_legacy_type(rep.electronic_structure));
        viamd::write_int(state,  STR_LIT("ElectronicStructureSource"),   (int)rep.electronic_structure.source);
        viamd::write_int(state,  STR_LIT("ElectronicStructureField"),    (int)electronic_structure_legacy_field(rep.electronic_structure));
        viamd::write_bool(state, STR_LIT("ElectronicStructureUseMagnitude"), rep.electronic_structure.use_magnitude);
        viamd::write_int(state,  STR_LIT("ElectronicStructureSpin"),     (int)rep.electronic_structure.spin);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoComponent"), (int)rep.electronic_structure.nto_component);
        viamd::write_int(state,  STR_LIT("ElectronicStructureTransitionDensityComponent"), (int)rep.electronic_structure.transition_density_component);
        viamd::write_int(state,  STR_LIT("ElectronicStructureRes"), (int)rep.electronic_structure.resolution);
        viamd::write_dbl(state,  STR_LIT("ElectronicStructureIso"),      rep.electronic_structure.iso_value);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColPos"),   rep.electronic_structure.col_psi_pos);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColNeg"),   rep.electronic_structure.col_psi_neg);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColDen"),   rep.electronic_structure.col_den);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColAtt"),   rep.electronic_structure.col_att);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColDet"),   rep.electronic_structure.col_det);
        viamd::write_int(state,  STR_LIT("ElectronicStructureColoring"), (int)rep.electronic_structure.coloring);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintPos"),  rep.electronic_structure.tint_psi_pos);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintNeg"),  rep.electronic_structure.tint_psi_neg);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintDen"),  rep.electronic_structure.tint_den);
        if (rep.electronic_structure.coloring == SurfaceColoring::Field) {
            const float range[2] = { rep.electronic_structure.field_map.range_beg, rep.electronic_structure.field_map.range_end };
            viamd::write_int(state,      STR_LIT("ElectronicStructureFieldKind"),      (int)rep.electronic_structure.field_kind);
            viamd::write_int(state,      STR_LIT("ElectronicStructureFieldColormap"),  rep.electronic_structure.field_map.colormap);
            viamd::write_flt_vec(state,  STR_LIT("ElectronicStructureFieldRange"),     range, 2);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldSymmetric"), rep.electronic_structure.field_map.symmetric);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldAutoRange"), rep.electronic_structure.field_map.auto_range);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldLegend"),    rep.electronic_structure.field_map.show_legend);
        }
        viamd::write_int(state,  STR_LIT("ElectronicStructureDensityPropertyIsoCount"), rep.electronic_structure.density_property.num_isos);
        for (int j = 0; j < rep.electronic_structure.density_property.num_isos; ++j) {
            char key[64];
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", j);
            viamd::write_dbl(state, str_from_cstr(key), rep.electronic_structure.density_property.values[j]);
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", j);
            viamd::write_vec4(state, str_from_cstr(key), rep.electronic_structure.density_property.colors[j]);
        }
    }
}

// A mask by the atoms it holds. Base64 now; a workspace written before that holds the raw bytes,
// which is read when it survived the text format (a mask whose bytes held no line break).
static bool deserialize_mask(md_bitfield_t* mask, str_t arg) {
    if (str_begins_with(arg, STR_LIT("###"))) {
        return viamd::extract_bitfield(mask, arg);
    }
    str_t raw;
    viamd::extract_str(raw, arg);
    return !str_empty(raw) && md_bitfield_deserialize(mask, raw.ptr, raw.len);
}

void load_workspace(ApplicationState* data, str_t filename) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    str_t txt = load_textfile(filename, temp_alloc);

    if (str_empty(txt)) {
        VIAMD_LOG_ERROR("Could not open workspace file: '" STR_FMT "'", STR_ARG(filename));
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
    md_bitfield_init(&pending.selection_mask, temp_alloc);
    md_bitfield_init(&pending.recenter_target, temp_alloc);

    str_t folder = {};
    extract_folder_path(&folder, filename);

    str_t section;
    while (viamd::next_section_header(section, state)) {
        str_t ident, arg;
        if (str_eq(section, STR_LIT("Files")) || str_eq(section, STR_LIT("File"))) {
            deserialize_files(state, pending, folder, temp_alloc);
        } else if (str_eq(section, STR_LIT("Animation"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Frame"))) {
                    pending.has_frame = viamd::extract_dbl(pending.frame, arg);
                } else if (str_eq(ident, STR_LIT("Fps"))) {
                    viamd::extract_flt(data->animation.fps, arg);
                } else if (str_eq(ident, STR_LIT("Tension"))) {
                    viamd::extract_flt(data->animation.tension, arg);
                } else if (str_eq(ident, STR_LIT("Interpolation"))) {
                    viamd::extract_enum(data->animation.interpolation, arg, (int)InterpolationMode::Count);
                }
            }
        } else if (str_eq(section, STR_LIT("TimelineView"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("FilterEnabled"))) {
                    viamd::extract_bool(data->timeline.filter.enabled, arg);
                } else if (str_eq(ident, STR_LIT("FilterFrames"))) {
                    float r[2];
                    if (viamd::extract_flt_vec(r, 2, arg)) {
                        pending.has_filter_range = true;
                        pending.filter_beg_frame = r[0];
                        pending.filter_end_frame = r[1];
                    }
                } else if (str_eq(ident, STR_LIT("TemporalWindowEnabled"))) {
                    viamd::extract_bool(data->timeline.filter.temporal_window.enabled, arg);
                } else if (str_eq(ident, STR_LIT("TemporalWindowExtent"))) {
                    viamd::extract_dbl(data->timeline.filter.temporal_window.extent_in_frames, arg);
                } else if (str_eq(ident, STR_LIT("ViewFrames"))) {
                    float r[2];
                    if (viamd::extract_flt_vec(r, 2, arg)) {
                        pending.has_view_range = true;
                        pending.view_beg_frame = r[0];
                        pending.view_end_frame = r[1];
                    }
                }
            }
        } else if (str_eq(section, STR_LIT("RenderSettings"))) {
            auto& v = data->visuals;
            while (viamd::next_entry(ident, arg, state)) {
                if      (str_eq(ident, STR_LIT("BackgroundColor")))      viamd::extract_vec3(v.background.color, arg);
                else if (str_eq(ident, STR_LIT("BackgroundIntensity")))  viamd::extract_flt(v.background.intensity, arg);
                else if (str_eq(ident, STR_LIT("SsaoEnabled")))          viamd::extract_bool(v.ssao.enabled, arg);
                // SsaoIntensity, SsaoRadius, SsaoBias and DofFocusScale belonged to the previous scale-dependent SSAO / DOF
                // and are deliberately ignored: their values do not translate to the new parameters.
                else if (str_eq(ident, STR_LIT("SsaoStrength")))         viamd::extract_flt(v.ssao.intensity, arg);
                else if (str_eq(ident, STR_LIT("TonemapEnabled")))       viamd::extract_bool(v.tonemapping.enabled, arg);
                else if (str_eq(ident, STR_LIT("Tonemapper")))           viamd::extract_enum(v.tonemapping.tonemapper, arg, (int)postprocess_pipeline::Tonemapper_ACES + 1);
                else if (str_eq(ident, STR_LIT("TonemapExposure")))      viamd::extract_flt(v.tonemapping.exposure, arg);
                else if (str_eq(ident, STR_LIT("TonemapGamma")))         viamd::extract_flt(v.tonemapping.gamma, arg);
                else if (str_eq(ident, STR_LIT("DofEnabled")))           viamd::extract_bool(v.dof.enabled, arg);
                else if (str_eq(ident, STR_LIT("DofAperture")))          viamd::extract_flt(v.dof.aperture, arg);
                else if (str_eq(ident, STR_LIT("FxaaEnabled")))          viamd::extract_bool(v.fxaa.enabled, arg);
                else if (str_eq(ident, STR_LIT("TaaEnabled")))           viamd::extract_bool(v.temporal_aa.enabled, arg);
                else if (str_eq(ident, STR_LIT("TaaJitter")))            viamd::extract_bool(v.temporal_aa.jitter, arg);
                else if (str_eq(ident, STR_LIT("TaaFeedbackMin")))       viamd::extract_flt(v.temporal_aa.feedback_min, arg);
                else if (str_eq(ident, STR_LIT("TaaFeedbackMax")))       viamd::extract_flt(v.temporal_aa.feedback_max, arg);
                else if (str_eq(ident, STR_LIT("MotionBlurEnabled")))    viamd::extract_bool(v.temporal_aa.motion_blur.enabled, arg);
                else if (str_eq(ident, STR_LIT("MotionBlurScale")))      viamd::extract_flt(v.temporal_aa.motion_blur.motion_scale, arg);
                else if (str_eq(ident, STR_LIT("SharpenEnabled")))       viamd::extract_bool(v.sharpen.enabled, arg);
                else if (str_eq(ident, STR_LIT("SharpenWeight")))        viamd::extract_flt(v.sharpen.weight, arg);
                else if (str_eq(ident, STR_LIT("SimulationBoxEnabled"))) viamd::extract_bool(data->simulation_box.enabled, arg);
                else if (str_eq(ident, STR_LIT("SimulationBoxColor")))   viamd::extract_vec4(data->simulation_box.color, arg);
            }
        } else if (str_eq(section, STR_LIT("Camera"))) {
            // A transform is only restored whole: half of one mixed with the default view is neither
            bool has_position = false, has_orientation = false, has_distance = false;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Position"))) {
                    has_position = viamd::extract_vec3(pending.camera.position, arg);
                } else if (str_eq(ident, STR_LIT("Orientation")) || str_eq(ident, STR_LIT("Rotation"))) {
                    // Rotation: DEPRECATED name
                    has_orientation = viamd::extract_quat(pending.camera.orientation, arg);
                } else if (str_eq(ident, STR_LIT("Distance"))) {
                    has_distance = viamd::extract_flt(pending.camera.distance, arg);
                } else if (str_eq(ident, STR_LIT("Mode"))) {
                    viamd::extract_enum(data->view.mode, arg, (int)CameraMode::Count);
                } else if (str_eq(ident, STR_LIT("FovY"))) {
                    viamd::extract_flt(data->view.camera.fov_y, arg);
                }
            }
            pending.has_camera = has_position && has_orientation && has_distance;
        } else if (str_eq(section, STR_LIT("Operations"))) {
            auto& op = data->operations;
            while (viamd::next_entry(ident, arg, state)) {
                if      (str_eq(ident, STR_LIT("Recenter")))            viamd::extract_bool(op.recenter, arg);
                else if (str_eq(ident, STR_LIT("FixateOrientation")))   viamd::extract_bool(op.fixate_orientation, arg);
                else if (str_eq(ident, STR_LIT("ApplyPbc")))            viamd::extract_bool(op.apply_pbc, arg);
                else if (str_eq(ident, STR_LIT("UnwrapStructures")))    viamd::extract_bool(op.unwrap_structures, arg);
                else if (str_eq(ident, STR_LIT("RecalcBonds")))         viamd::extract_bool(op.recalc_bonds, arg);
                else if (str_eq(ident, STR_LIT("RecenterQueryEnabled"))) viamd::extract_bool(op.recenter_query.enabled, arg);
                else if (str_eq(ident, STR_LIT("RecenterQuery"))) {
                    viamd::extract_to_char_buf(op.recenter_query.query, sizeof(op.recenter_query.query), arg);
                    recenter_mark_query_dirty(data);
                }
                else if (str_eq(ident, STR_LIT("RecenterTarget")))      pending.has_recenter_target = deserialize_mask(&pending.recenter_target, arg);
            }
        } else if (str_eq(section, STR_LIT("Representation"))) {
            deserialize_representation(data, state);
        } else if (str_eq(section, STR_LIT("UserBonds"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("atoms"))) {
                    int atom_indices[2] = {-1, -1};
                    if (viamd::extract_int_vec(atom_indices, 2, arg) && atom_indices[0] >= 0 && atom_indices[1] >= 0) {
                        md_atom_pair_t pair = { atom_indices[0], atom_indices[1] };
                        md_array_push(pending.user_bonds, pair, temp_alloc);
                    }
                }
            }
        } else if (str_eq(section, STR_LIT("Script"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Text"))) {
                    str_t str;
                    viamd::extract_str(str, arg);
                    data->editor.SetText(std::string(str.ptr, str.len));
                }
            }
        } else if (str_eq(section, STR_LIT("Selection"))) {
            str_t label = {};
            md_bitfield_t mask = {};
            md_bitfield_init(&mask, temp_alloc);
            bool has_mask = false;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Label"))) {
                    viamd::extract_str(label, arg);
                } else if (str_eq(ident, STR_LIT("Mask"))) {
                    has_mask = deserialize_mask(&mask, arg);
                }
            }
            if (!str_empty(label) && has_mask) {
                Selection* sel = create_selection(data, label);
                md_bitfield_copy(&sel->atom_mask, &mask);
            }
        } else if (str_eq(section, STR_LIT("ActiveSelection"))) {
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Granularity"))) {
                    viamd::extract_enum(data->selection.granularity, arg, (int)SelectionGranularity::Count);
                } else if (str_eq(ident, STR_LIT("Mask"))) {
                    pending.has_selection_mask = deserialize_mask(&pending.selection_mask, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Windows"))) {
            // Windows not named keep their state: a workspace older than this section says nothing about them
            while (viamd::next_entry(ident, arg, state)) {
                for (size_t i = 0; i < num_workspace_windows; ++i) {
                    if (str_eq_cstr(ident, workspace_windows[i].name)) {
                        viamd::extract_bool(*workspace_windows[i].show, arg);
                        break;
                    }
                }
            }
        } else if (plot_layout_deserialize(state, STR_LIT("Timeline"), STR_LIT("TimelineSeries"), data, data->timeline.subplots, &data->timeline.num_subplots)) {
        } else if (plot_layout_deserialize(state, STR_LIT("Distributions"), STR_LIT("DistributionSeries"), data, data->distributions.subplots, &data->distributions.num_subplots)) {
        } else {
            const char* before = state.text.ptr;
            viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserialize, viamd::EventPayloadType_DeserializationState, &state);
            if (state.text.ptr == before && viamd::next_entry(ident, arg, state)) {
                // Written by a component this build does not have, or by a newer version
                MD_LOG_DEBUG("Workspace: nothing reads the section [" STR_FMT "], it is skipped", STR_ARG(section));
            }
        }
    }

    // ## 3. Load the files, then apply what refers to them
    str_copy_to_char_buf(data->files.workspace, sizeof(data->files.workspace), filename);
    data->files.coarse_grained = pending.coarse_grained;

    loader::LoaderState loader_state = {};
    if (!str_empty(pending.molecule_file)) {
        loader::init(&loader_state, pending.molecule_file);
        if (pending.coarse_grained) {
            loader_state.flags |= LoaderFlag_CoarseGrained;
        }
        if (load_data_from_file(data, pending.molecule_file, loader_state)) {
            str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), pending.molecule_file);
        } else {
            data->files.molecule[0] = '\0';
        }
    } else {
        data->files.molecule[0] = '\0';
    }

    if (!str_empty(pending.trajectory_file)) {
        loader::init(&loader_state, pending.trajectory_file);
        if (load_data_from_file(data, pending.trajectory_file, loader_state)) {
            str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), pending.trajectory_file);
        }
    } else {
        data->files.trajectory[0] = '\0';
    }

    // Joins the trajectory's run, so it can only go in once the trajectory is there
    if (!str_empty(pending.energy_file)) {
        loader::init(&loader_state, pending.energy_file, &data->mold.sys);
        load_data_from_file(data, pending.energy_file, loader_state);
    }
    for (size_t i = 0; i < md_array_size(pending.series_files); ++i) {
        loader::init(&loader_state, pending.series_files[i], &data->mold.sys);
        load_data_from_file(data, pending.series_files[i], loader_state);
    }

    const size_t num_atoms  = data->mold.sys.atom.count;
    const size_t num_frames = run_num_frames(data);

    if (pending.has_frame && num_frames > 0) {
        data->animation.frame = CLAMP(pending.frame, 0.0, (double)(num_frames - 1));
    }

    // Loading a trajectory sets the filter and the view to all of it; the workspace's come after
    if (num_frames > 0) {
        const double last = (double)(num_frames - 1);
        if (pending.has_filter_range) {
            data->timeline.filter.beg_frame = CLAMP(pending.filter_beg_frame, 0.0, last);
            data->timeline.filter.end_frame = CLAMP(pending.filter_end_frame, data->timeline.filter.beg_frame, last);
        }
        if (pending.has_view_range && pending.view_end_frame > pending.view_beg_frame) {
            data->timeline.view_range.beg_x = workspace_frame_to_time(data, CLAMP(pending.view_beg_frame, 0.0, last));
            data->timeline.view_range.end_x = workspace_frame_to_time(data, CLAMP(pending.view_end_frame, 0.0, last));
        }
    }

    // Masks of atoms that are not there are dropped rather than applied to whatever holds those indices now
    auto mask_fits = [num_atoms](const md_bitfield_t* mask) {
        uint64_t first = 0, last = 0;
        return num_atoms > 0 && (md_bitfield_empty(mask) || (md_bitfield_get_range(&first, &last, mask) && last < num_atoms));
    };
    if (pending.has_selection_mask && mask_fits(&pending.selection_mask)) {
        md_bitfield_copy(&data->selection.selection_mask, &pending.selection_mask);
    }
    if (pending.has_recenter_target && mask_fits(&pending.recenter_target)) {
        md_bitfield_copy(&data->operations.selection_mask, &pending.recenter_target);
        recenter_update_target_data(data);
    }

    const size_t num_user_bonds = md_array_size(pending.user_bonds);
    if (num_user_bonds > 0) {
        for (size_t i = 0; i < num_user_bonds; ++i) {
            const md_atom_pair_t& pair = pending.user_bonds[i];
            if ((size_t)pair.idx[0] < num_atoms && (size_t)pair.idx[1] < num_atoms) {
                md_system_bond_insert(&data->mold.sys, pair.idx[0], pair.idx[1], md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER));
            }
        }
        md_util_system_infer_coordination(&data->mold.sys);
        data->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
    }

    // The camera flies in: from the whole system to where the workspace looked from. Without a
    // stored camera, both are the whole system.
    reset_view(&data->view.camera, data->mold.state, &data->representation.visibility_mask);
    data->view.target = pending.has_camera ? pending.camera : (ViewTransform)data->view.camera;

    viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserializeEnd, viamd::EventPayloadType_DeserializationState, &state);
}

bool save_workspace(ApplicationState* app_state, str_t filename) {
    md_allocator_i* temp_alloc = app_state->allocator.frame;

    viamd::serialization_state_t state {
        .filename = filename,
        .sb = md_strb_create(temp_alloc),
    };

    constexpr str_t header_snippet = STR_LIT(
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
    state.sb += header_snippet;
    state.sb += '\n';

    // Files are stored relative to the workspace so that a dataset can be moved as a whole.
    // If no relative path exists (a separate volume on windows) we fall back to the absolute
    // path, and a slot which holds no file is written as an empty string.
    auto workspace_relative_path = [filename, temp_alloc](str_t path) -> str_t {
        if (str_empty(path)) {
            return {};
        }
        str_t rel = md_path_make_relative(filename, path, temp_alloc);
        return str_empty(rel) ? md_path_make_canonical(path, temp_alloc) : rel;
    };

    const md_attributes_t* attributes = &app_state->mold.sys.attributes;

    viamd::write_section_header(state, STR_LIT("Files"));
    viamd::write_str(state, STR_LIT("MoleculeFile"),   workspace_relative_path(str_from_cstr(app_state->files.molecule)));
    viamd::write_str(state, STR_LIT("TrajectoryFile"), workspace_relative_path(str_from_cstr(app_state->files.trajectory)));

    // An energy file is only in the session while its data is: the source path is published beside
    // the energies and removed with them, so what is written here cannot name a file that was
    // dropped with an earlier trajectory.
    {
        if (const md_attribute_t* src = run_attribute(app_state, STR_LIT("edr/source"))) {
            viamd::write_str(state, STR_LIT("EnergyFile"), workspace_relative_path(md_attribute_str(attributes, src, 0)));
        }

        // Series loaded along the run (.xvg, .csv), the same way: "<run>/<kind>/<name>/source"
        const str_t kinds[] = { STR_INIT("xvg"), STR_INIT("csv") };
        for (size_t k = 0; k < ARRAY_SIZE(kinds); ++k) {
            char group_buf[256];
            const str_t group = run_attribute_path(group_buf, sizeof(group_buf), app_state, kinds[k]);
            if (str_empty(group)) continue;
            for (md_attribute_iter_t it = md_attributes_iter_children(attributes, group); md_attributes_next(&it);) {
                if (const md_attribute_t* series_src = md_attributes_find_in(attributes, it.child_path, STR_LIT("source"))) {
                    viamd::write_str(state, STR_LIT("SeriesFile"), workspace_relative_path(md_attribute_str(attributes, series_src, 0)));
                }
            }
        }
    }
    viamd::write_bool(state, STR_LIT("CoarseGrained"), app_state->files.coarse_grained);

    viamd::write_section_header(state, STR_LIT("Animation"));
    viamd::write_dbl(state, STR_LIT("Frame"), app_state->animation.frame);
    viamd::write_flt(state, STR_LIT("Fps"), app_state->animation.fps);
    viamd::write_flt(state, STR_LIT("Tension"), app_state->animation.tension);
    viamd::write_int(state, STR_LIT("Interpolation"), (int)app_state->animation.interpolation);

    plot_layout_serialize(state, STR_LIT("Timeline"), STR_LIT("TimelineSeries"), app_state->timeline.subplots, app_state->timeline.num_subplots);
    plot_layout_serialize(state, STR_LIT("Distributions"), STR_LIT("DistributionSeries"), app_state->distributions.subplots, app_state->distributions.num_subplots);

    {
        // In frames: a time is in whatever unit it happens to be shown in
        const auto& tl = app_state->timeline;
        viamd::write_section_header(state, STR_LIT("TimelineView"));
        viamd::write_bool(state, STR_LIT("FilterEnabled"), tl.filter.enabled);
        const float filter[2] = { (float)tl.filter.beg_frame, (float)tl.filter.end_frame };
        viamd::write_flt_vec(state, STR_LIT("FilterFrames"), filter, 2);
        viamd::write_bool(state, STR_LIT("TemporalWindowEnabled"), tl.filter.temporal_window.enabled);
        viamd::write_dbl(state, STR_LIT("TemporalWindowExtent"), tl.filter.temporal_window.extent_in_frames);
        if (md_array_size(tl.x_values) > 0) {
            const float view[2] = {
                (float)workspace_time_to_frame(app_state, tl.view_range.beg_x),
                (float)workspace_time_to_frame(app_state, tl.view_range.end_x),
            };
            viamd::write_flt_vec(state, STR_LIT("ViewFrames"), view, 2);
        }
    }

    {
        const auto& v = app_state->visuals;
        viamd::write_section_header(state, STR_LIT("RenderSettings"));
        viamd::write_vec3(state, STR_LIT("BackgroundColor"), v.background.color);
        viamd::write_flt(state,  STR_LIT("BackgroundIntensity"), v.background.intensity);
        viamd::write_bool(state, STR_LIT("SsaoEnabled"), v.ssao.enabled);
        viamd::write_flt(state,  STR_LIT("SsaoStrength"), v.ssao.intensity);
        viamd::write_bool(state, STR_LIT("TonemapEnabled"), v.tonemapping.enabled);
        viamd::write_int(state,  STR_LIT("Tonemapper"), (int)v.tonemapping.tonemapper);
        viamd::write_flt(state,  STR_LIT("TonemapExposure"), v.tonemapping.exposure);
        viamd::write_flt(state,  STR_LIT("TonemapGamma"), v.tonemapping.gamma);
        viamd::write_bool(state, STR_LIT("DofEnabled"), v.dof.enabled);
        viamd::write_flt(state,  STR_LIT("DofAperture"), v.dof.aperture);
        viamd::write_bool(state, STR_LIT("FxaaEnabled"), v.fxaa.enabled);
        viamd::write_bool(state, STR_LIT("TaaEnabled"), v.temporal_aa.enabled);
        viamd::write_bool(state, STR_LIT("TaaJitter"), v.temporal_aa.jitter);
        viamd::write_flt(state,  STR_LIT("TaaFeedbackMin"), v.temporal_aa.feedback_min);
        viamd::write_flt(state,  STR_LIT("TaaFeedbackMax"), v.temporal_aa.feedback_max);
        viamd::write_bool(state, STR_LIT("MotionBlurEnabled"), v.temporal_aa.motion_blur.enabled);
        viamd::write_flt(state,  STR_LIT("MotionBlurScale"), v.temporal_aa.motion_blur.motion_scale);
        viamd::write_bool(state, STR_LIT("SharpenEnabled"), v.sharpen.enabled);
        viamd::write_flt(state,  STR_LIT("SharpenWeight"), v.sharpen.weight);
        viamd::write_bool(state, STR_LIT("SimulationBoxEnabled"), app_state->simulation_box.enabled);
        viamd::write_vec4(state, STR_LIT("SimulationBoxColor"), app_state->simulation_box.color);
    }

    // The target, not the camera: the camera eases towards it over several frames (camera_animate),
    // so saving during a fly-in or right after a drag would store a point along the way
    viamd::write_section_header(state, STR_LIT("Camera"));
    viamd::write_vec3(state, STR_LIT("Position"), app_state->view.target.position);
    viamd::write_quat(state, STR_LIT("Orientation"), app_state->view.target.orientation);
    viamd::write_flt(state,  STR_LIT("Distance"), app_state->view.target.distance);
    viamd::write_int(state,  STR_LIT("Mode"), (int)app_state->view.mode);
    viamd::write_flt(state,  STR_LIT("FovY"), app_state->view.camera.fov_y);

    {
        const auto& op = app_state->operations;
        viamd::write_section_header(state, STR_LIT("Operations"));
        viamd::write_bool(state, STR_LIT("Recenter"), op.recenter);
        viamd::write_bool(state, STR_LIT("FixateOrientation"), op.fixate_orientation);
        viamd::write_bool(state, STR_LIT("ApplyPbc"), op.apply_pbc);
        viamd::write_bool(state, STR_LIT("UnwrapStructures"), op.unwrap_structures);
        viamd::write_bool(state, STR_LIT("RecalcBonds"), op.recalc_bonds);
        viamd::write_bool(state, STR_LIT("RecenterQueryEnabled"), op.recenter_query.enabled);
        viamd::write_str(state,  STR_LIT("RecenterQuery"), str_from_cstr(op.recenter_query.query));
        if (!md_bitfield_empty(&op.selection_mask)) {
            viamd::write_bitfield(state, STR_LIT("RecenterTarget"), &op.selection_mask);
        }
    }

    for (size_t i = 0; i < md_array_size(app_state->representation.reps); ++i) {
        serialize_representation(state, app_state, app_state->representation.reps[i]);
    }

    {
        std::string text = app_state->editor.GetText();
        viamd::write_section_header(state, STR_LIT("Script"));
        viamd::write_str(state, STR_LIT("Text"), str_t{text.c_str(), text.size()});
    }

    for (size_t i = 0; i < md_array_size(app_state->selection.stored_selections); ++i) {
        const Selection& sel = app_state->selection.stored_selections[i];
        viamd::write_section_header(state, STR_LIT("Selection"));
        viamd::write_str(state, STR_LIT("Label"), str_from_cstr(sel.name));
        viamd::write_bitfield(state, STR_LIT("Mask"), &sel.atom_mask);
    }

    viamd::write_section_header(state, STR_LIT("ActiveSelection"));
    viamd::write_int(state, STR_LIT("Granularity"), (int)app_state->selection.granularity);
    if (!md_bitfield_empty(&app_state->selection.selection_mask)) {
        viamd::write_bitfield(state, STR_LIT("Mask"), &app_state->selection.selection_mask);
    }

    // Save user defined bonds
    bool has_user_bonds = false;
    for (size_t i = 0; i < app_state->mold.sys.bond.count; ++i) {
        if (md_bond_origin(app_state->mold.sys.bond.flags[i]) == MD_BOND_ORIGIN_USER) {
            if (!has_user_bonds) {
                viamd::write_section_header(state, STR_LIT("UserBonds"));
                has_user_bonds = true;
            }
            viamd::write_int_vec(state, STR_LIT("atoms"), app_state->mold.sys.bond.pairs[i].idx, 2);
        }
    }

    viamd::write_section_header(state, STR_LIT("Windows"));
    for (size_t i = 0; i < num_workspace_windows; ++i) {
        viamd::write_bool(state, str_from_cstr(workspace_windows[i].name), *workspace_windows[i].show);
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSerialize, viamd::EventPayloadType_SerializationState, &state);

    // Only now, with all of it in hand, is the file touched
    const str_t text = md_strb_to_str(state.sb);
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Could not open workspace file for writing: '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    const size_t written = md_file_write(file, str_ptr(text), str_len(text));
    md_file_close(&file);
    if (written != str_len(text)) {
        VIAMD_LOG_ERROR("Could not write all of the workspace file: '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    // Saved there, it is the workspace from here on: the next plain save goes to the same file
    str_copy_to_char_buf(app_state->files.workspace, sizeof(app_state->files.workspace), filename);
    VIAMD_LOG_SUCCESS("Saved workspace '" STR_FMT "'", STR_ARG(filename));
    return true;
}

// --- SELECTION ---
Selection* create_selection(ApplicationState* state, str_t name, md_bitfield_t* atom_mask) {
    ASSERT(state);
    Selection sel;
    str_copy_to_char_buf(sel.name, sizeof(sel.name), name);
    md_bitfield_init(&sel.atom_mask, state->allocator.persistent);
    if (atom_mask) {
        md_bitfield_copy(&sel.atom_mask, atom_mask);
    }
    md_array_push(state->selection.stored_selections, sel, state->allocator.persistent);
    return md_array_last(state->selection.stored_selections);
}

void remove_all_selections(ApplicationState* state) {
    md_array_shrink(state->selection.stored_selections, 0);
}

void remove_selection(ApplicationState* state, size_t idx) {
    ASSERT(state);
    if (md_array_size(state->selection.stored_selections) <= idx) {
        VIAMD_LOG_ERROR("Index [%zu] out of range when trying to remove selection", idx);
    }
    auto item = &state->selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);

    state->selection.stored_selections[idx] = *md_array_last(state->selection.stored_selections);
    md_array_pop(state->selection.stored_selections);
}

// --- Representation ---

static void init_representation(ApplicationState* state, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(state->mold.sys.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(state->mold.gl_mol);
    md_bitfield_init(&rep->atom_mask, state->allocator.persistent);

    // Default to the first per atom field the system offers, if it offers any - unless the
    // representation already names one the system has. This runs again for every representation
    // when a system is loaded, which is AFTER a workspace's representations were read, and for a
    // clone: neither may lose the attribute it was coloured by, or how. There is no list to consult:
    // the attribute table is the list.
    if (!md_attributes_get(&state->mold.sys.attributes, rep->atom_attribute.key)) {
        md_attribute_id_t first_attribute = MD_ATTRIBUTE_INVALID;
        if (atom_attribute_query(&first_attribute, 1, state->mold.sys) > 0) {
            atom_attribute_select(&rep->atom_attribute, first_attribute, state->mold.sys);
        }
    }

    // A system loaded under a representation that colours by a field: the field's values were of the
    // old one, whatever the frame says
    surface_field_invalidate(&rep->electronic_structure.field_vol);

    flag_representation_as_dirty(rep);
}

Representation* create_representation(ApplicationState* state, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(state);
    md_array_push(state->representation.reps, Representation(), state->allocator.persistent);
    Representation* rep = md_array_last(state->representation.reps);
    rep->type = type;
    rep->color_mapping = color_mapping;
    if (!str_empty(filter)) {
        str_copy_to_char_buf(rep->filt, sizeof(rep->filt), filter);
    }
    // Opens on the HOMO, derived from the occupations rather than handed over by whoever loaded
    // the file. -1 (nothing occupied) clamps to 0, which is the only orbital there is to show.
    OrbitalFrontier frontier = {};
    es_orbital_frontier(&frontier, state->mold.sys, es_path::alpha_occupation);
    rep->electronic_structure.orbital_idx = MAX(frontier.homo_idx, 0);
    init_representation(state, rep);
    return rep;
}

Representation* clone_representation(ApplicationState* state, const Representation& rep) {
    ASSERT(state);
    md_array_push(state->representation.reps, rep, state->allocator.persistent);
    Representation* clone = md_array_last(state->representation.reps);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    // The volumes' textures and buffers are the original's: the clone evaluates its own, and must
    // neither write into nor free those (each representation frees its own on removal)
    clone->electronic_structure.density_vol.tex_id = 0;
    clone->electronic_structure.color_vol.tex_id   = 0;
    clone->electronic_structure.field_vol          = SurfaceFieldVolume{};
    init_representation(state, clone);
    return clone;
}

void remove_representation(ApplicationState* state, size_t idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(state->representation.reps));
    auto& rep = state->representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_rep_destroy(rep.md_rep);
    // A readback queued for this representation's volume would land in a texture that no longer
    // exists, so let the queue run out first.
    gpu_volume_jobs_drain(state);
    if (rep.electronic_structure.density_vol.tex_id) gl::free_texture(&rep.electronic_structure.density_vol.tex_id);
    if (rep.electronic_structure.color_vol.tex_id)   gl::free_texture(&rep.electronic_structure.color_vol.tex_id);
    surface_field_free(&rep.electronic_structure.field_vol);
    md_array_swap_back_and_pop(state->representation.reps, idx);
    recompute_atom_visibility_mask(state);
}

void recompute_atom_visibility_mask(ApplicationState* state) {
    ASSERT(state);
    auto& mask = state->representation.visibility_mask;

    md_bitfield_clear(&mask);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }
    state->representation.visibility_mask_hash = md_bitfield_hash64(&mask, 0);
}

// "ground_state" reads as "Ground State" in a menu. The path segment is the identity; this is
// presentation only, which is why it lives here and not in mdlib. Writes into a caller buffer so
// gathering stays allocation free.
int dipole_label_pretty(char* buf, size_t cap, str_t group) {
    if (!buf || cap == 0) return 0;
    size_t n = MIN(group.len, cap - 1);
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
    return md_attributes_sibling(&sys.attributes, vec, STR_LIT("origin"));
}

// A vector and an origin form a group only if both are there, both are 3 component, and the origin
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
    for (md_attribute_iter_t it = md_attributes_iter_children(&sys.attributes, STR_LIT("dipole")); md_attributes_next(&it);) {
        const md_attribute_t* vec = md_attributes_find_in(&sys.attributes, it.child_path, STR_LIT("vector"));
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
        str_t label = vec->path;
        size_t sep = 0;
        if (str_rfind_char(&sep, label, '/')) label = str_substr(label, 0, sep);
        if (str_rfind_char(&sep, label, '/')) label = str_substr(label, sep + 1, SIZE_MAX);

        *out = {
            .key   = key,
            .label = label,
            .count = (uint32_t)md_attribute_value_count(&vec->format),
            .unit  = vec->unit,
        };
    }
    return true;
}

bool dipole_moment_read(vec3_t* out_vec, vec3_t* out_origin, const md_system_t& sys, md_attribute_id_t key, uint32_t index) {
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
    md_attribute_extract_f32(v, ARRAY_SIZE(v), vec, vec_slice, md_unit_none());
    md_attribute_extract_f32(o, ARRAY_SIZE(o), org, org_slice, md_unit_angstrom());

    if (out_vec)    *out_vec    = vec3_set(v[0], v[1], v[2]);
    if (out_origin) *out_origin = vec3_set(o[0], o[1], o[2]);
    return true;
}


// A per atom scalar field: values one component wide, the atom axis LAST, and at most one axis of
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
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, STR_LIT("atom")); md_attributes_next(&it);) {
        if (!atom_attribute_qualifies(it.attr, sys)) continue;

        if (out_ids && count < cap) {
            out_ids[count] = it.attr->id;
        }
        count += 1;
    }

    return count;
}

str_t atom_attribute_label(const md_attribute_t* attr) {
    if (!attr) return str_t{};
    // An empty label is a valid state, and the leaf is what the path spells for itself.
    return str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
}

int atom_attribute_variant_count(const md_attribute_t* attr) {
    if (!attr) return 0;
    return attr->format.rank > 1 ? (int)attr->format.shape[0] : 1;
}

ColorScaleSpan atom_attribute_span(const md_attribute_t* attr, const md_bitfield_t* mask) {
    ColorScaleSpan span;
    if (!attr) return span;

    const size_t num_values = md_attribute_element_count(&attr->format);
    const size_t num_atoms  = attr->format.rank > 0 ? attr->format.shape[attr->format.rank - 1] : 0;
    if (num_values == 0 || num_atoms == 0) return span;

    md_temp_scope_t temp = md_temp_begin();
    float* values = (float*)md_temp_alloc(temp, sizeof(float) * num_values);

    // Deliberately the whole attribute and not one variant: a span measured per variant would make
    // the colours shift as the index slider moves, which reads as the data changing.
    if (values && md_attribute_extract_f32(values, num_values, attr, md_attribute_slice_all(), md_unit_none()) == num_values) {
        float value_min =  FLT_MAX;
        float value_max = -FLT_MAX;
        size_t num_present = 0;
        for (size_t i = 0; i < num_values; ++i) {
            if (mask && !md_bitfield_test_bit(mask, i % num_atoms)) continue;
            if (atom_attribute_value_absent(values[i])) continue;
            value_min = MIN(value_min, values[i]);
            value_max = MAX(value_max, values[i]);
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
    ASSERT(buf && cap > 0);
    const str_t label = atom_attribute_label(attr);
    const int num_variants = atom_attribute_variant_count(attr);
    if (num_variants > 1) {
        snprintf(buf, cap, "%.*s [%d/%d]", (int)label.len, label.ptr, CLAMP(variant_idx, 0, num_variants - 1) + 1, num_variants);
    } else {
        snprintf(buf, cap, "%.*s", (int)label.len, label.ptr);
    }
}

void atom_attribute_select(AtomAttributeColoring* coloring, md_attribute_id_t key, const md_system_t& sys) {
    ASSERT(coloring);

    coloring->key = key;
    coloring->variant_idx = 0;
    coloring->span_hash = 0;    // measured again, over the representation's atoms, with its colours

    // Started out from the values of every atom: the representation's own are not at hand here
    const ColorScaleSpan span = atom_attribute_span(md_attributes_get(&sys.attributes, key), nullptr);
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
// Built from the system's own basis/ attributes, so it needs no loader and no component - which is
// the point of publishing the basis in the first place. Lives beside gl_mol because it is the same
// kind of thing: derived FROM the system so that something can draw or evaluate it, and never read
// back by mdlib.
//
// The device scratch is grown here rather than allocated per dataset. A wider basis loaded later
// grows it; a narrower one leaves it alone, because it is scratch and only the maximum matters.

void system_gpu_data_free(ApplicationState* state) {
    ASSERT(state);
#if MD_ENABLE_GPU
    if (state->mold.gpu_basis) {
        md_gto_gpu_basis_destroy(state->gpu_stream, state->mold.gpu_basis);
        state->mold.gpu_basis = nullptr;
    }
    if (state->mold.gpu_atoms) {
        md_gpu_free(state->gpu_stream, state->mold.gpu_atoms);
        state->mold.gpu_atoms = 0;
    }
    state->mold.gpu_atoms_hash = 0;
#else
    (void)state;
#endif
}

bool system_gpu_data_update(ApplicationState* state, double cutoff) {
    ASSERT(state);
#if MD_ENABLE_GPU
    system_gpu_data_free(state);

    if (!state->gpu_device) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

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
// component to ask - which is the whole point of publishing the coefficients.
//
// The basis is rebuilt per call for now. That is one interleave over a few hundred shells, but it
// is the obvious thing for a representation to cache, keyed on the ids of the basis/ attributes it
// was built from.
double* orbital_coefficients_extract(size_t* out_num_ao, md_temp_scope_t temp, const md_system_t& sys, str_t coefficient_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, coefficient_path);
    if (!attr) {
        MD_LOG_DEBUG("No orbital coefficients published at '" STR_FMT "'", STR_ARG(coefficient_path));
        return nullptr;
    }

    // Whatever the attribute's own rank, the slice has to leave exactly ONE row of AO coefficients:
    // {M,A} sliced by the orbital, {S,L,A} sliced by state and lambda. Asking the slice for its
    // format rather than the attribute is what makes those the same call.
    md_attribute_format_t format = {};
    if (!md_attribute_slice_format(&format, attr, slice)) {
        MD_LOG_ERROR("The slice does not address '" STR_FMT "'", STR_ARG(coefficient_path));
        return nullptr;
    }
    if (format.rank != 1 || format.components != 1) {
        MD_LOG_ERROR("'" STR_FMT "' does not slice down to one row of orbital coefficients", STR_ARG(coefficient_path));
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
    if (md_attribute_extract_f64(dst, num_ao, attr, slice, md_unit_none()) != num_ao) {
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
static const str_t QM_ATOM_MAP_PATH = STR_INIT("qm/atom/system_index");

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

bool es_orbital_frontier(OrbitalFrontier* out, const md_system_t& sys, str_t occupation_path) {
    ASSERT(out);
    *out = {};

    const md_attribute_t* attr = md_attributes_find(&sys.attributes, occupation_path);
    const double* occ = (const double*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F64, 1, 1);
    if (!occ) {
        return false;
    }

    const int num = (int)attr->format.shape[0];
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

const double* es_orbital_energies(size_t* out_count, const md_system_t& sys, str_t energy_path) {
    if (out_count) *out_count = 0;

    const md_attribute_t* attr = md_attributes_find(&sys.attributes, energy_path);
    const double* energies = (const double*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F64, 1, 1);
    if (energies && out_count) *out_count = attr->format.shape[0];
    return energies;
}

// NULL when there is none, which the callers read as the identity. Same explicitness as
// gto_attr_column in md_gto.c and for the same reason: the table is open, so a path with the wrong
// shape is a mistake to refuse rather than to reinterpret.
static const uint32_t* qm_atom_map_find(size_t* out_count, const md_system_t& sys) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, QM_ATOM_MAP_PATH);
    const uint32_t* map = (const uint32_t*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_U32, 1, 1);
    if (attr && !map) {
        MD_LOG_ERROR("'" STR_FMT "' is not a plain column of atom indices", STR_ARG(QM_ATOM_MAP_PATH));
    }
    if (map && out_count) *out_count = attr->format.shape[0];
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

double es_electron_count(const md_system_t& sys, str_t occupation_path) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, occupation_path);
    const double* occ = (const double*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F64, 1, 1);
    if (!occ) {
        return 0.0;
    }
    double sum = 0.0;
    for (uint32_t i = 0; i < attr->format.shape[0]; ++i) {
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
    defer { md_temp_end(temp); };
    double* values = (double*)md_temp_alloc(temp, sizeof(double) * row);
    if (!values || md_attribute_extract_f64(values, row, attr, slice, md_unit_none()) != row) {
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
    ASSERT(out);
    *out = {};

    // The element column defines the domain: it is published for every QM system, and its length
    // is the atom count. Without it there is no QM atom space to speak of.
    const md_attribute_t* z = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/atomic_number"));
    if (!z || !z->data || z->format.type != MD_ATTRIBUTE_TYPE_U8 ||
        z->format.rank != 1 || z->format.components != 1) {
        return false;
    }
    out->count         = z->format.shape[0];
    out->atomic_number = (const uint8_t*)z->data;

    if (const md_attribute_t* c = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/coordinate"))) {
        if (c->data && c->format.type == MD_ATTRIBUTE_TYPE_F64 && c->format.rank == 1 &&
            c->format.components == 3 && c->format.shape[0] == out->count) {
            // dvec3_t is three doubles with no padding, and the {N,3} layout is exactly an
            // array of them, so this is a view and not a reinterpretation.
            out->coordinate = (const dvec3_t*)c->data;
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
static size_t basis_atom_positions_gather(vec3_t* dst, size_t cap, const md_system_t& sys, const md_system_state_t& state, size_t num_basis_atoms) {
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
        MD_LOG_ERROR("The basis spans %zu atoms and '" STR_FMT "' maps %zu", num_basis_atoms, STR_ARG(QM_ATOM_MAP_PATH), map_count);
        return 0;
    }

    for (size_t i = 0; i < num_basis_atoms; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= state.num_atoms) {
            MD_LOG_ERROR("Basis atom %zu is system atom %zu and the state holds %zu", i, idx, state.num_atoms);
            return 0;
        }
        dst[i] = vec3_set(state.xyz[idx].x, state.xyz[idx].y, state.xyz[idx].z) * (float)ANGSTROM_TO_BOHR;
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
        MEMCPY(dst, src, size);
        if (gl::pbo_upload_end_texture_3D(vol_tex, 0, GL_R32F)) return;
    }
    gl::set_texture_3D_data(vol_tex, 0, src, GL_R32F);
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
    }
    if (self && job->rb.gpu) md_gpu_free(self->gpu_stream, job->rb.gpu);
    job->rb        = {};
    job->in_flight = false;
}

// Used when no job slot is free, and when a caller genuinely needs the data before it returns.
static bool gpu_volume_readback_blocking(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, size_t size) {
    md_gpu_mem_t rb = md_gpu_malloc(state->gpu_stream, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    const md_gpu_tex_region_t region = {
        .offset = {0, 0, 0},
        .extent = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] },
    };
    bool ok = md_gpu_copy_from_texture(state->gpu_stream, rb.gpu, state->gpu_volume, &region);
    md_gpu_stream_sync(state->gpu_stream);
    if (ok) gpu_volume_upload_to_gl(vol_tex, rb.cpu, size);
    md_gpu_free(state->gpu_stream, rb.gpu);
    return ok;
}

// Queues a readback of the evaluated region of gpu_volume. Returns once the copy is recorded; the
// GL texture is filled later, from gpu_volume_job_complete(), which md_gpu_device_poll() calls on
// the GL thread.
//
// Returns true when the work was QUEUED -- not that the texture holds data.
static bool gpu_volume_readback(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid) {
    if (!state->gpu_stream || !state->gpu_volume) {
        return false;
    }
    const size_t size = sizeof(float) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];

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
        return gpu_volume_readback_blocking(state, vol_tex, grid, size);
    }

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
    if (!md_gpu_copy_from_texture(state->gpu_stream, rb.gpu, state->gpu_volume, &region)) {
        MD_LOG_ERROR("Failed to record the volume readback");
        md_gpu_free(state->gpu_stream, rb.gpu);
        job->in_flight = false;
        return false;
    }

    job->tex_id = vol_tex;
    job->rb     = rb;
    job->size   = size;

    if (!md_gpu_launch_host_fn(state->gpu_stream, gpu_volume_job_complete, job)) {
        // No callback means nothing would ever free rb or fill the texture, so finish this one
        // synchronously instead of leaking it.
        MD_LOG_ERROR("Failed to queue the volume completion; falling back to a blocking readback");
        md_gpu_stream_sync(state->gpu_stream);
        gpu_volume_upload_to_gl(vol_tex, rb.cpu, size);
        md_gpu_free(state->gpu_stream, rb.gpu);
        job->in_flight = false;
        return true;
    }
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
    defer { md_temp_end(temp); };

    // Same gather as the GL path, through the same map, so the two backends cannot disagree about
    // which atom a basis index means. Packed to vec4 because that is what the device wants.
    vec3_t* xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * MAX(num_atoms, (size_t)1));
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, sys_state, num_atoms) != num_atoms) {
        return false;
    }

    // Gather first, then decide: the gather is what the comparison is ABOUT, and it is cheap next
    // to the transfer it may save.
    const uint64_t hash = md_hash64(xyz, sizeof(vec3_t) * num_atoms, 0);
    if (hash != 0 && hash == state->mold.gpu_atoms_hash) {
        return true;
    }

    vec4_t* xyzw = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * MAX(num_atoms, (size_t)1));
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = vec4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
    }

    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->mold.gpu_atoms, sz);
    if (!dst) {
        MD_LOG_ERROR("Failed to upload the atom positions to the device");
        return false;
    }
    md_gto_gpu_atom_pack(dst, (const float*)xyzw, sizeof(vec4_t), num_atoms);
    if (!md_gpu_upload_end(state->gpu_stream)) {
        return false;
    }

    state->mold.gpu_atoms_hash = hash;
    return true;
}

#endif // MD_ENABLE_GPU

void gpu_volume_jobs_drain(ApplicationState* state) {
    ASSERT(state);
#if MD_ENABLE_GPU
    if (!state->gpu_stream || !gpu_volume_job_any_in_flight(state)) return;
    md_gpu_stream_sync(state->gpu_stream);
    md_gpu_device_poll(state->gpu_device);   // this is what actually runs the callbacks

    // Backstop: if a callback somehow did not fire, release the staging block here rather than
    // leaking it, since the pool may be about to go.
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        if (j.in_flight) {
            if (j.rb.gpu) md_gpu_free(state->gpu_stream, j.rb.gpu);
            j.rb = {};
            j.in_flight = false;
        }
    }
#else
    (void)state;
#endif
}

vec3_t* basis_atom_positions_extract(size_t* out_count, md_temp_scope_t temp, const md_system_t& sys, const md_system_state_t& state) {
    if (out_count) *out_count = 0;

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return nullptr;
    }

    const size_t num_basis_atoms = md_gto_basis_num_atoms(&basis);
    vec3_t* atom_xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * MAX(num_basis_atoms, (size_t)1));
    if (!atom_xyz || basis_atom_positions_gather(atom_xyz, num_basis_atoms, sys, state, num_basis_atoms) != num_basis_atoms) {
        return nullptr;
    }

    if (out_count) *out_count = num_basis_atoms;
    return atom_xyz;
}

// The basis an evaluation runs against, rebuilt from the system's basis/ attributes. Positions are
// NOT part of it: they come in as a parameter, because which geometry an orbital is drawn at is the
// caller's decision and nothing down here is in a position to make it.
static bool gto_basis_context(md_gto_basis_t* out_basis, md_temp_scope_t temp, const md_system_t& sys) {
    ASSERT(out_basis);
    MEMSET(out_basis, 0, sizeof(*out_basis));
    return md_gto_basis_extract_attributes(out_basis, &sys.attributes, md_temp_allocator(temp));
}

// Every evaluation checks its positions the same way, and the count is the only thing that can be
// wrong about them: an array of the wrong length would put shells on the wrong nuclei.
static bool gto_positions_match_basis(const md_gto_basis_t* basis, const vec3_t* atom_pos, size_t num_atom_pos) {
    const size_t expected = md_gto_basis_num_atoms(basis);
    if (!atom_pos || num_atom_pos != expected) {
        MD_LOG_ERROR("The basis spans %zu atoms and %zu positions were supplied", expected, num_atom_pos);
        return false;
    }
    return true;
}

bool orbital_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                         const vec3_t* atom_pos, size_t num_atom_pos,
                         str_t coefficient_path, const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t  num_ao    = 0;
    double* ao_coeffs = orbital_coefficients_extract(&num_ao, temp, sys, coefficient_path, slice);
    if (!ao_coeffs) {
        return false;
    }

    md_gto_basis_t basis = {};
    if (!gto_basis_context(&basis, temp, sys) || !gto_positions_match_basis(&basis, atom_pos, num_atom_pos)) {
        return false;
    }

    // The coefficients and the basis are separate attributes, so nothing guarantees they agree
    // until it is checked here.
    if (md_gto_basis_num_ao(&basis) != num_ao) {
        MD_LOG_ERROR("The basis spans %zu atomic orbitals and '" STR_FMT "' %zu", md_gto_basis_num_ao(&basis), STR_ARG(coefficient_path), num_ao);
        return false;
    }

    md_gto_grid_evaluate_mo_GL(vol_tex, &grid, &basis, (const float*)atom_pos, sizeof(vec3_t), ao_coeffs, cutoff, mode, op);
    volume::notify_data_changed(vol_tex);
    return true;
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
double* density_matrix_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, str_t density_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" STR_FMT "'", STR_ARG(density_path));
        return nullptr;
    }

    // Square or packed, stored or computed on demand: md_qm reads either, and says why when the
    // slice does not narrow it to one symmetric matrix
    const size_t dim = md_qm_extract_symmetric_f64(nullptr, 0, attr, slice);
    if (dim == 0) {
        return nullptr;
    }
    double* dst = (double*)md_temp_alloc(temp, sizeof(double) * dim * dim);
    if (!dst || md_qm_extract_symmetric_f64(dst, dim * dim, attr, slice) != dim) {
        return nullptr;
    }

    if (out_dim) *out_dim = dim;
    return dst;
}

// The same matrix as its packed upper triangle in float - what both density paths hand the shader.
// A packed attribute is extracted as it is and nothing larger is ever made; a square one is
// extracted whole into scratch first.
float* density_packed_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, str_t density_path, const md_attribute_slice_t* slice_ptr) {
    const md_attribute_slice_t slice = slice_ptr ? *slice_ptr : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" STR_FMT "'", STR_ARG(density_path));
        return nullptr;
    }

    const size_t dim = md_qm_extract_packed_symmetric_f32(nullptr, 0, attr, slice);
    if (dim == 0) {
        return nullptr;
    }
    const size_t len = dim * (dim + 1) / 2;
    float* dst = (float*)md_temp_alloc(temp, sizeof(float) * len);
    if (!dst || md_qm_extract_packed_symmetric_f32(dst, len, attr, slice) != dim) {
        return nullptr;
    }

    if (out_dim) *out_dim = dim;
    return dst;
}

static bool density_packed_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                                       const vec3_t* atom_pos, size_t num_atom_pos,
                                       const float* packed, size_t dim, md_gto_op_t op) {
    if (!packed || dim == 0) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!gto_basis_context(&basis, temp, sys) || !gto_positions_match_basis(&basis, atom_pos, num_atom_pos)) {
        return false;
    }

    if (md_gto_basis_num_ao(&basis) != dim) {
        MD_LOG_ERROR("The basis spans %zu atomic orbitals and the density matrix %zu", md_gto_basis_num_ao(&basis), dim);
        return false;
    }

    md_gto_grid_evaluate_density_packed_GL(vol_tex, &grid, &basis, (const float*)atom_pos, sizeof(vec3_t), packed, false, op);
    volume::notify_data_changed(vol_tex);
    return true;
}

bool density_matrix_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                                const vec3_t* atom_pos, size_t num_atom_pos,
                                const double* density_matrix, size_t dim, md_gto_op_t op) {
    if (!density_matrix || dim == 0) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!gto_basis_context(&basis, temp, sys) || !gto_positions_match_basis(&basis, atom_pos, num_atom_pos)) {
        return false;
    }

    if (md_gto_basis_num_ao(&basis) != dim) {
        MD_LOG_ERROR("The basis spans %zu atomic orbitals and the density matrix %zu", md_gto_basis_num_ao(&basis), dim);
        return false;
    }

    md_gto_grid_evaluate_density_GL(vol_tex, &grid, &basis, (const float*)atom_pos, sizeof(vec3_t), density_matrix, false, op);
    volume::notify_data_changed(vol_tex);
    return true;
}

bool density_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                         const vec3_t* atom_pos, size_t num_atom_pos,
                         str_t density_path, const md_attribute_slice_t* slice, md_gto_op_t op) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    // Packed from the start: the shader reads the upper triangle as float, so that is all that is made
    size_t dim = 0;
    const float* packed = density_packed_extract(&dim, temp, sys, density_path, slice);
    if (!packed) {
        return false;
    }

    return density_packed_evaluate_gl(vol_tex, grid, sys, atom_pos, num_atom_pos, packed, dim, op);
}

static bool es_attribute_exists(const md_system_t& sys, str_t path) {
    return md_attributes_find(&sys.attributes, path) != nullptr;
}

// Evaluates whatever the representation is currently pointed at into its own volume texture.
// Everything it needs is an attribute on the system or a field of the representation; there is no
// reader, no component and no event in the path, which is the whole point of the port.
static bool electronic_structure_evaluate(ApplicationState* state, Representation* rep) {
    const ElectronicStructureRepresentation& es = rep->electronic_structure;

    const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)es.resolution];

    md_grid_t grid = {};
    if (!electronic_structure_grid_init(&grid, state->mold.sys, state->mold.state, samples_per_angstrom * BOHR_TO_ANGSTROM)) {
        return false;
    }
    init_volume(&rep->electronic_structure.density_vol, grid, GL_R32F);
    rep->electronic_structure.grid = grid;
    const uint32_t tex_id = rep->electronic_structure.density_vol.tex_id;

    switch (es.source) {
    case ElectronicStructureSource::MolecularOrbital: {
        const str_t path = (es.spin == ElectronicStructureSpin::Beta) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.orbital_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::NaturalTransitionOrbital: {
        const str_t path = (es.nto_component == ElectronicStructureNtoComponent::Particle) ? es_path::nto_particle : es_path::nto_hole;
        // {S,L,A}: the excited state and then the lambda pair within it. Two indices instead of the
        // one an ordinary orbital takes, which is exactly what a slice is for.
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)es.excited_state_idx, (uint32_t)es.nto_lambda_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::TransitionDensity: {
        str_t path = es_path::attachment_density;
        switch (es.transition_density_component) {
        case ElectronicStructureTransitionDensityComponent::Detachment: path = es_path::detachment_density; break;
        case ElectronicStructureTransitionDensityComponent::Difference: path = es_path::transition_diff;    break;
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default: break;
        }
        // These are VIRTUAL: the slice is what tells the provider to reconstruct one state rather
        // than all of them, and it is the only reason asking for one is affordable.
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.excited_state_idx);
        return density_evaluate(state, tex_id, grid, path, &slice, MD_GTO_OP_SET);
    }
    case ElectronicStructureSource::ElectronDensity: {
        // A spin difference is signed, so it is the one density the magnitude toggle applies to.
        const bool magnitude = (es.spin == ElectronicStructureSpin::Difference) && es.use_magnitude;
        return density_evaluate(state, tex_id, grid, es_electron_density_path(es.spin), nullptr, es_gto_op(magnitude));
    }
    case ElectronicStructureSource::DensityProperty: {
        const md_attribute_t* attr = md_attributes_get(&state->mold.sys.attributes, es.density_property_key);
        if (!attr) {
            return false;
        }
        return density_evaluate(state, tex_id, grid, attr->path, nullptr, MD_GTO_OP_SET);
    }
    default:
        MD_LOG_ERROR("Unknown electronic structure source");
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
    defer { md_temp_end(temp); };

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = state->mold.state;

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
        MAX(1, DIV_UP((int)rep->electronic_structure.density_vol.dim[0], downsample_factor)),
        MAX(1, DIV_UP((int)rep->electronic_structure.density_vol.dim[1], downsample_factor)),
        MAX(1, DIV_UP((int)rep->electronic_structure.density_vol.dim[2], downsample_factor)),
    };
    MEMCPY(rep->electronic_structure.color_vol.dim, dim, sizeof(dim));
    rep->electronic_structure.color_vol.world_to_model   = rep->electronic_structure.density_vol.world_to_model;
    rep->electronic_structure.color_vol.texture_to_world = rep->electronic_structure.density_vol.texture_to_world;
    rep->electronic_structure.color_vol.voxel_size       = vec3_set(
        rep->electronic_structure.density_vol.voxel_size.x * rep->electronic_structure.density_vol.dim[0] / dim[0],
        rep->electronic_structure.density_vol.voxel_size.y * rep->electronic_structure.density_vol.dim[1] / dim[1],
        rep->electronic_structure.density_vol.voxel_size.z * rep->electronic_structure.density_vol.dim[2] / dim[2]);
    gl::init_texture_3D(&rep->electronic_structure.color_vol.tex_id, dim[0], dim[1], dim[2], GL_RGBA8);

    const vec3_t& voxel_size     = rep->electronic_structure.color_vol.voxel_size;
    const mat4_t& world_to_model = rep->electronic_structure.color_vol.world_to_model;
    mat4_t index_to_world = rep->electronic_structure.color_vol.texture_to_world
        * mat4_scale(1.0f / dim[0], 1.0f / dim[1], 1.0f / dim[2])
        * mat4_translate(0.5f, 0.5f, 0.5f); // Center of the corner voxel should be at the origin

    // The colours are per SYSTEM atom and the splats are per BASIS atom, so the map is consulted
    // here too - and once more it is the same lookup in the same direction the evaluation used.
    size_t map_count = 0;
    const uint32_t* map = qm_atom_map_find(&map_count, sys);
    if (map && map_count < num_points) {
        return;
    }

    vec4_t*   point_xyzw   = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * num_points);
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
        point_xyzw[i]   = vec4_set(sys_state.xyz[idx].x, sys_state.xyz[idx].y, sys_state.xyz[idx].z, radius);
        point_colors[i] = atom_colors[idx];
    }

    volume::compute_point_color_volume(rep->electronic_structure.color_vol.tex_id, dim, voxel_size.elem, world_to_model.elem,
                                       index_to_world.elem, point_xyzw, point_colors, num_points,
                                       rep->electronic_structure.gaussian_splatting_power);
}

// The values of a volume texture, read back. A volume the md_gpu path is still filling is waited
// for first, and image stores of the GL compute path are made visible before the read.
static bool volume_read_values(float* dst, ApplicationState* state, const Volume& vol) {
    if (!vol.tex_id) return false;
    gpu_volume_jobs_drain(state);
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
static void electronic_structure_field_update(ApplicationState* state, Representation* rep) {
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    const md_system_t& sys = state->mold.sys;
    if (!surface_field_available(es.field_kind, sys)) {
        surface_field_free(&es.field_vol);
        return;
    }

    IsoDesc iso;
    electronic_structure_iso_desc_init(&iso, es);

    // The values are of the field at this frame; the surfaces are what the density volume holds
    // (vol_hash names it) cut at the isovalues
    uint64_t source_hash = md_hash64(&es.field_kind, sizeof(es.field_kind), 0);
    source_hash = md_hash64(&state->animation.frame, sizeof(state->animation.frame), source_hash);
    uint64_t band_hash = md_hash64(iso.values, sizeof(float) * iso.count, (uint64_t)iso.count + 1);
    band_hash = md_hash64_combine(band_hash, es.vol_hash);
    if (es.field_vol.tex_id && source_hash == es.field_vol.source_hash && band_hash == es.field_vol.band_hash) {
        return;
    }

    const size_t num_voxels = md_grid_num_points(&es.grid);
    if (num_voxels == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    float* density = (float*)md_temp_alloc(temp, sizeof(float) * num_voxels);
    if (!density || !volume_read_values(density, state, es.density_vol)) {
        return;
    }

    SurfaceFieldDesc desc = {};
    desc.kind        = es.field_kind;
    desc.sys         = &sys;
    desc.state       = &state->mold.state;
    desc.grid        = &es.grid;
    desc.density     = density;
    desc.iso_values  = iso.values;
    desc.num_iso     = iso.count;
    desc.source_hash = source_hash;
    desc.band_hash   = band_hash;
    if (es.field_kind == SurfaceFieldKind::ElectrostaticPotential) {
        // The QM nuclei and basis where the density was drawn: the current frame, as for the density
        desc.basis_atom_xyz = basis_atom_positions_extract(&desc.num_basis_atoms, temp, sys, state->mold.state);
    }
#if MD_ENABLE_GPU
    if (state->gpu_stream && state->gpu_volume) {
        desc.gpu_stream  = state->gpu_stream;
        desc.gpu_scratch = state->gpu_volume;
        desc.gpu_scratch_dim[0] = desc.gpu_scratch_dim[1] = desc.gpu_scratch_dim[2] = GPU_VOLUME_DIM;
    }
#endif

    const md_tick_t t0 = md_tick_now();
    if (surface_field_update(&es.field_vol, desc)) {
        color_scale_update_range(&es.field_map, surface_field_span(es.field_vol));
        const SurfaceFieldVolume& fv = es.field_vol;
        MD_LOG_DEBUG("Surface field: %zu of %zu voxels (%dx%dx%d) evaluated%s, %.1f ms; on %zu surface samples min %g, 1%% %g, 99%% %g, max %g",
                     fv.num_evaluated, fv.num_voxels, fv.grid.dim[0], fv.grid.dim[1], fv.grid.dim[2], fv.on_gpu ? " on the GPU" : "",
                     md_tick_to_milliseconds(md_tick_now() - t0), fv.num_surface_samples, fv.surface_min, fv.surface_lo, fv.surface_hi, fv.surface_max);
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
    ASSERT(grid);

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return false;
    }

    const size_t num_atoms = md_gto_basis_num_atoms(&basis);
    if (num_atoms == 0) {
        return false;
    }

    vec3_t* xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * num_atoms);
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, state, num_atoms) != num_atoms) {
        return false;
    }

    // mat3_PCA and calculate_bounds both want a homogeneous point, and w == 1 is what makes the
    // rotation in calculate_bounds a rotation of a POINT rather than of a direction.
    vec4_t* xyzw = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * num_atoms);
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = vec4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
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

bool orbital_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t coefficient_path,
                      const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    ASSERT(state);

#if MD_ENABLE_GPU
    if (gpu_evaluation_ready(state)) {
        const md_system_t&       sys       = state->mold.sys;
        const md_system_state_t& sys_state = state->mold.state;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        size_t  num_ao    = 0;
        const double* ao_coeffs = orbital_coefficients_extract(&num_ao, temp, sys, coefficient_path, slice);
        if (!ao_coeffs) {
            return false;
        }

        const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
        if (num_cgtos != num_ao) {
            MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and '" STR_FMT "' %zu", num_cgtos, STR_ARG(coefficient_path), num_ao);
            return false;
        }
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
            .out_tex      = state->gpu_volume,
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

    // The call site's choice of geometry, made explicit: this wrapper draws at the system's CURRENT
    // state. A caller wanting another one gathers its own positions and calls the _gl form.
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t num_atom_pos = 0;
    const vec3_t* atom_pos = basis_atom_positions_extract(&num_atom_pos, temp, state->mold.sys, state->mold.state);
    return orbital_evaluate_gl(vol_tex, grid, state->mold.sys, atom_pos, num_atom_pos, coefficient_path, slice, mode, op, cutoff);
}

#if MD_ENABLE_GPU
static void density_gpu_launch(ApplicationState* state, const md_grid_t& grid, md_gto_op_t op);

bool density_matrix_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid,
                                           const double* density_matrix, size_t dim, md_gto_op_t op) {
    ASSERT(state);
    if (!density_matrix || dim == 0 || !gpu_evaluation_ready(state)) {
        return false;
    }

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = state->mold.state;

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    if (num_cgtos != dim) {
        MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and the density matrix %zu", num_cgtos, dim);
        return false;
    }
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

// The density kernel over whatever was last uploaded to the coefficient buffer
static void density_gpu_launch(ApplicationState* state, const md_grid_t& grid, md_gto_op_t op) {
    md_gto_gpu_density_desc_t desc = {
        .basis         = state->mold.gpu_basis,
        .atom_xyz      = state->mold.gpu_atoms,
        .coeff         = state->gpu_coeff,
        .out_tex       = state->gpu_volume,
        .grid          = &grid,
        .sample_offset = {0.5f, 0.5f, 0.5f},
        .op            = op,
    };
    md_gto_gpu_density_launch(state->gpu_stream, &desc);
}

bool density_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, str_t density_path,
                                    const md_attribute_slice_t* slice, md_gto_op_t op) {
    ASSERT(state);
    if (!gpu_evaluation_ready(state)) {
        return false;
    }

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = state->mold.state;
    const md_attribute_slice_t s = slice ? *slice : md_attribute_slice_all();
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" STR_FMT "'", STR_ARG(density_path));
        return false;
    }

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    const size_t dim = md_qm_extract_packed_symmetric_f32(nullptr, 0, attr, s);
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

    // Straight into the upload buffer: the packed float triangle is exactly what the kernel reads, so
    // a packed attribute - a density property read from its file - never exists in any other form
    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_density(num_cgtos));
    if (!dst) {
        return false;
    }
    const bool extracted = md_qm_extract_packed_symmetric_f32(dst, dim * (dim + 1) / 2, attr, s) == dim;
    md_gpu_upload_end(state->gpu_stream);
    if (!extracted) {
        return false;
    }

    density_gpu_launch(state, grid, op);
    return true;
}
#endif

// The backend choice, over a density matrix the caller already holds. density_evaluate below is
// this plus the extract, and a caller which has cached the matrix - the QM component does, because
// the transition densities are rebuilt on every read - calls this one and skips the rebuild.
bool density_matrix_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid,
                             const double* density_matrix, size_t dim, md_gto_op_t op) {
    ASSERT(state);
    if (!density_matrix || dim == 0) {
        return false;
    }

#if MD_ENABLE_GPU
    if (density_matrix_evaluate_to_gpu_volume(state, grid, density_matrix, dim, op)) {
        return gpu_volume_readback(state, vol_tex, grid);
    }
    if (gpu_evaluation_ready(state)) {
        return false;   // the GPU path was available and failed; the GL one would fail the same way
    }
#endif

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t num_atom_pos = 0;
    const vec3_t* atom_pos = basis_atom_positions_extract(&num_atom_pos, temp, state->mold.sys, state->mold.state);
    return density_matrix_evaluate_gl(vol_tex, grid, state->mold.sys, atom_pos, num_atom_pos, density_matrix, dim, op);
}

bool density_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t density_path,
                      const md_attribute_slice_t* slice, md_gto_op_t op) {
    ASSERT(state);

#if MD_ENABLE_GPU
    // Evaluate, then read back. The split exists because a consumer that runs another kernel over
    // the result - the critical point extraction does - wants the first half and not the second.
    if (density_evaluate_to_gpu_volume(state, grid, density_path, slice, op)) {
        return gpu_volume_readback(state, vol_tex, grid);
    }
    if (gpu_evaluation_ready(state)) {
        return false;   // the GPU path was available and failed; the GL one would fail the same way
    }
#endif

    // As above: the current state, chosen here rather than assumed further down.
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t num_atom_pos = 0;
    const vec3_t* atom_pos = basis_atom_positions_extract(&num_atom_pos, temp, state->mold.sys, state->mold.state);
    return density_evaluate_gl(vol_tex, grid, state->mold.sys, atom_pos, num_atom_pos, density_path, slice, op);
}

void update_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        update_representation(state, &state->representation.reps[i]);
    }
}

bool representation_uses_atom_colors(const Representation& rep) {
    switch (rep.type) {
        case RepresentationType::ElectronicStructure:
            return rep.electronic_structure.coloring == SurfaceColoring::AtomColors;
        case RepresentationType::DipoleMoment:
            return false;
        default:
            return true;
    }
}

void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    if (!rep->enabled) return;
    if (!rep->needs_update) return;

    const auto& sys = state->mold.sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(temp); };

    const size_t bytes = num_atoms * sizeof(uint32_t);

    //md_script_property_t prop = {0};
    //if (rep->color_mapping == ColorMapping::Attribute) {
    //rep->prop_is_valid = md_script_compile_and_eval_property(&prop, rep->prop, &data->mold.sys, frame_allocator, &data->script.ir, rep->prop_error.beg(), rep->prop_error.capacity());
    //}

    uint32_t* colors = 0;
    if (representation_uses_atom_colors(*rep)) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * num_atoms);

        // The atoms first: the colours are applied through them, and a colour scale's range can
        // follow the values of exactly the atoms shown
        if (rep->dynamic_evaluation) {
            rep->filt_is_dirty = true;
        }
        if (rep->filt_is_dirty) {
            rep->filt_is_valid = md_filter(&rep->atom_mask, str_from_cstr(rep->filt), &state->mold.sys, &state->mold.state, state->script.ir, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
            rep->filt_is_dirty = false;
        }

        switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, num_atoms, convert_color(rep->base_color));
            break;
        case ColorMapping::Type:
            color_atoms_type(colors, num_atoms, sys);
            break;
        case ColorMapping::Serial:
            color_atoms_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::CompName:
            color_atoms_comp_name(colors, num_atoms, sys);
            break;
        case ColorMapping::CompSeqId:
            color_atoms_comp_seq_id(colors, num_atoms, sys);
            break;
        case ColorMapping::CompIndex:
            color_atoms_comp_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::InstId:
            color_atoms_inst_id(colors, num_atoms, sys);
            break;
        case ColorMapping::InstIndex:
            color_atoms_inst_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::SecondaryStructure: {
            SecondaryStructurePalette palette = {
                .coil  = convert_color(rep->secondary_structure.color_coil),
                .helix = convert_color(rep->secondary_structure.color_helix),
                .sheet = convert_color(rep->secondary_structure.color_sheet),
            };
            color_atoms_secondary_structure(colors, num_atoms, sys, displayed_secondary_structure(state), palette);
            break;
        }
        case ColorMapping::Attribute: {
            AtomAttributeColoring& coloring = rep->atom_attribute;
            const md_attribute_t* attr = md_attributes_get(&sys.attributes, coloring.key);
            size_t num_extracted = 0;
            float* values = nullptr;

            if (attr) {
                values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);

                // Fixing the variant axis hands back exactly the atom axis, so there is no
                // offset arithmetic here to get wrong. A field with no variant axis is rank 1
                // and takes no indices at all.
                const uint32_t variant = (uint32_t)CLAMP(coloring.variant_idx, 0, MAX(atom_attribute_variant_count(attr) - 1, 0));
                const md_attribute_slice_t slice = attr->format.rank > 1 ? md_attribute_slice_1(variant) : md_attribute_slice_all();
                num_extracted = md_attribute_extract_f32(values, num_atoms, attr, slice, md_unit_none());
            }

            if (num_extracted == num_atoms) {
                // The span of the atoms shown, measured again only when the values or the atoms
                // changed: it scans every variant of the attribute
                const md_bitfield_t* shown = rep->filt_is_valid ? &rep->atom_mask : nullptr;
                uint64_t span_hash = md_hash64(&coloring.key, sizeof(coloring.key), 1);
                const uint64_t version = md_attributes_version(&sys.attributes, coloring.key);
                span_hash = md_hash64(&version, sizeof(version), span_hash);
                span_hash = shown ? md_bitfield_hash64(shown, span_hash) : md_hash64(&num_atoms, sizeof(num_atoms), span_hash);
                if (span_hash != coloring.span_hash) {
                    coloring.span = atom_attribute_span(attr, shown);
                    coloring.span_hash = span_hash;
                }
                color_scale_update_range(&coloring.scale, coloring.span);

                for (size_t i = 0; i < num_atoms; ++i) {
                    // An atom without a value is not a point on the ramp: NAN would reach the
                    // colormap lookup as an index, so it gets a neutral grey of its own instead
                    colors[i] = atom_attribute_value_absent(values[i]) ? IM_COL32(128, 128, 128, 255) : color_scale_color_u32(coloring.scale, values[i]);
                }
            } else {
                if (attr) {
                    MD_LOG_DEBUG("Failed to extract values for the selected atom attribute");
                }
                MEMSET(colors, 0xFFFFFFFFu, bytes);
            }
        }
#if 0
            if (rep->prop) {
                MEMSET(colors, 0xFFFFFFFF, bytes);
                md_script_pro
                    const float* values = rep->prop->data.values;
                if (rep->prop->data.aggregate) {
                    const int dim = rep->prop->data.dim[0];
                    md_script_vis_t vis = {0};
                    bool result = false;

                    //if (md_semaphore_aquire(&data->script.ir_semaphore)) {
                    //    defer { md_semaphore_release(&data->script.ir_semaphore); };

                    if (md_script_ir_valid(state->script.eval_ir)) {
                        md_script_vis_init(&vis, frame_alloc);
                        md_script_vis_ctx_t ctx = {
                            .ir = state->script.eval_ir,
                            .mol = &state->mold.sys,
                            .traj = state->mold.sys.trajectory,
                        };
                        result = md_script_vis_eval_payload(&vis, rep->prop->vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS);
                    }
                    //}
                    if (result) {
                        if (dim == (int)md_array_size(vis.structure)) {
                            int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)state->animation.frame);

                            md_bitfield_t mask = {0};
                            md_bitfield_init(&mask, frame_alloc);
                            for (int i = 0; i < dim; ++i) {
                                md_bitfield_and(&mask, &rep->atom_mask, &vis.structure[i]);
                                float value = lerpf(values[i0 * dim + i], values[i1 * dim + i], frame_fract);
                                float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                                ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                                color_atoms_uniform(colors, mol.atom.count, vec_cast(color), &mask);
                            }
                        }
                    }
                } else {
                    int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values - 1);
                    int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values - 1);
                    float value = lerpf(values[i0], values[i1], fractf((float)state->animation.frame));
                    float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                    ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                    color_atoms_uniform(colors, mol.atom.count, vec_cast(color));
                }
            } else {
                color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            }
#endif
            break;
        default:
            ASSERT(false);
            break;
        }
    }

    if (colors && (rep->tint_scale > 0.0f || rep->saturation < 1.0f)) {
        uint32_t tint_color = convert_color(rep->tint_color);
        tint_colors(colors, num_atoms, tint_color, rep->tint_scale, rep->saturation);
    }

    switch (rep->type) {
    case RepresentationType::SpaceFill:
        rep->type_is_valid = sys.atom.count > 0;
        break;
    case RepresentationType::Licorice:
        rep->type_is_valid = sys.bond.count > 0;
        break;
    case RepresentationType::BallAndStick:
        rep->type_is_valid = sys.atom.count > 0;
        break;
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = sys.protein_backbone.range.count > 0;
        break;
    case RepresentationType::ElectronicStructure: {
        rep->type_is_valid = electronic_structure_source_supported(es_source_mask(sys), rep->electronic_structure.source);
        if (rep->type_is_valid && rep->enabled) {
            // Re-evaluating is expensive and everything it depends on is right here, so it is
            // gated on a hash of exactly those inputs - the frame included, since the geometry
            // moves the grid.
            uint64_t vol_hash = md_hash64(&state->animation.frame, sizeof(state->animation.frame), 0);
            vol_hash = md_hash64(&rep->electronic_structure.source,                       sizeof(rep->electronic_structure.source), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.use_magnitude,                sizeof(rep->electronic_structure.use_magnitude), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.spin,                         sizeof(rep->electronic_structure.spin), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.nto_component,                sizeof(rep->electronic_structure.nto_component), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.transition_density_component, sizeof(rep->electronic_structure.transition_density_component), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.resolution,                   sizeof(rep->electronic_structure.resolution), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.orbital_idx,                  sizeof(rep->electronic_structure.orbital_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.excited_state_idx,            sizeof(rep->electronic_structure.excited_state_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.nto_lambda_idx,               sizeof(rep->electronic_structure.nto_lambda_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.density_property_key,         sizeof(rep->electronic_structure.density_property_key), vol_hash);

            if (vol_hash != rep->electronic_structure.vol_hash) {
                rep->electronic_structure.vol_hash = vol_hash;
                electronic_structure_evaluate(state, rep);
            }

            if (rep->electronic_structure.coloring == SurfaceColoring::AtomColors) {
                uint64_t col_hash = md_hash64(&rep->electronic_structure.gaussian_splatting_power, sizeof(rep->electronic_structure.gaussian_splatting_power), (int)rep->color_mapping);
                col_hash = md_hash64_combine(col_hash, rep->electronic_structure.vol_hash);
                if (col_hash != rep->electronic_structure.col_hash) {
                    rep->electronic_structure.col_hash = col_hash;
                    electronic_structure_color_volume_update(state, rep, colors);
                }
            } else if (rep->electronic_structure.coloring == SurfaceColoring::Field) {
                electronic_structure_field_update(state, rep);
            }
        }
        break;
    }
	case RepresentationType::DipoleMoment:
		rep->type_is_valid = md_attributes_get(&state->mold.sys.attributes, rep->dipole.dipole_key) != NULL;
		break;
    default:
        ASSERT(false);
        break;
    }

    if (colors) {
        if (rep->filt_is_valid) {
            filter_colors(colors, num_atoms, &rep->atom_mask);
            state->representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_atom_colors(rep->md_rep, 0, (uint32_t)num_atoms, colors, 0);

#if EXPERIMENTAL_GFX_API
            md_gfx_rep_attr_t attributes = {};
            attributes.spacefill.radius_scale = 1.0f;
            md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
            md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
        }
    }

    rep->needs_update = false;
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
            const str_t label = str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
            out[count] = { .key = attr->id, .label = label };
        }
        count += 1;
    }
    return count;
}

ElectronicStructureSourceFlags es_source_mask(const md_system_t& sys) {
    ElectronicStructureSourceFlags mask = 0;

    if (es_attribute_exists(sys, es_path::alpha_coefficient))  mask |= ElectronicStructureSourceFlag_MolecularOrbital;
    if (es_attribute_exists(sys, es_path::alpha_density))      mask |= ElectronicStructureSourceFlag_ElectronDensity;
    if (es_attribute_exists(sys, es_path::nto_particle))       mask |= ElectronicStructureSourceFlag_NaturalTransitionOrbital;
    if (es_attribute_exists(sys, es_path::attachment_density)) mask |= ElectronicStructureSourceFlag_TransitionDensity;
    if (density_properties_gather(nullptr, 0, sys) > 0)        mask |= ElectronicStructureSourceFlag_DensityProperty;

    return mask;
}

bool electronic_structure_select_available_source(ElectronicStructureRepresentation* es, const md_system_t& sys) {
    ASSERT(es);
    const ElectronicStructureSourceFlags mask = es_source_mask(sys);
    if (mask == 0) {
        return false;
    }

    if (!electronic_structure_source_supported(mask, es->source)) {
        for (int n = 0; n < (int)ElectronicStructureSource::Count; ++n) {
            const ElectronicStructureSource source = (ElectronicStructureSource)n;
            if (electronic_structure_source_supported(mask, source)) {
                es->source = source;
                electronic_structure_set_source_defaults(es);
                break;
            }
        }
    }

    // The same rule the representation window applies when it draws: a key naming none of the
    // properties this system has is replaced by the first one.
    if (es->source == ElectronicStructureSource::DensityProperty) {
        DensityProperty props[64];
        const size_t num_props = MIN(density_properties_gather(props, ARRAY_SIZE(props), sys), ARRAY_SIZE(props));
        bool found = false;
        for (size_t i = 0; i < num_props; ++i) {
            if (props[i].key == es->density_property_key) {
                found = true;
                break;
            }
        }
        if (!found && num_props > 0) {
            es->density_property_key = props[0].key;
        }
    }

    return true;
}


// Per atom scalar fields are deliberately not gathered anywhere. They live in the system's attribute
// table under atom/, whoever loaded the data put them there, and the UI reads that table directly
// through atom_attribute_query.

static void init_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        init_representation(state, &rep);
    }
}

void flag_representation_as_dirty(Representation* rep) {
    ASSERT(rep);
    rep->filt_is_dirty = true;
    rep->needs_update  = true;
    rep->electronic_structure.vol_hash = 0;
}

void flag_all_representations_as_dirty(ApplicationState* state) {
    ASSERT(state);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        flag_representation_as_dirty(&state->representation.reps[i]);
    }
}

void remove_all_representations(ApplicationState* state) {
    while (md_array_size(state->representation.reps) > 0) {
        remove_representation(state, (int32_t)md_array_size(state->representation.reps) - 1);
    }
}


void create_default_representations(ApplicationState* state) {
    bool amino_acid_present = false;
    bool nucleic_present = false;
    bool ion_present = false;
    bool water_present = false;
    bool ligand_present = false;
    size_t num_coarse_grained = 0;
    // Any electronic structure, not orbitals in particular: a file can carry density properties and
    // no SCF block, and its density properties are as much something to show as orbitals would be.
    const bool electronic_structure_present = es_source_mask(state->mold.sys) != 0;


    if (state->mold.sys.atom.count > 3'000'000) {
        VIAMD_LOG_INFO("Large system detected, creating default representation for all atoms");
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (state->mold.sys.component.count == 0) {
        // No residues present
        Representation* rep = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    // What is present, by what the components are. The selections of the representations below are by the same.
    for (size_t i = 0; i < state->mold.sys.component.count; ++i) {
        switch (md_component_kind(&state->mold.sys.component, i)) {
        case MD_COMPONENT_KIND_AMINO_ACID: amino_acid_present = true; break;
        case MD_COMPONENT_KIND_NUCLEOTIDE: nucleic_present = true;    break;
        case MD_COMPONENT_KIND_ION:        ion_present = true;        break;
        case MD_COMPONENT_KIND_WATER:      water_present = true;      break;
        default:                           ligand_present = true;     break;
        }
    }
    for (size_t i = 0; i < state->mold.sys.atom.count; ++i) {
        num_coarse_grained += md_atom_particle_kind(&state->mold.sys.atom, i) == MD_PARTICLE_BEAD;
    }

    // Coarse grained when most of it is. A handful of beads, or atoms a loader could not assign an
    // element to, should not take the protein, ligand and water representations from the rest.
    if (2 * num_coarse_grained > state->mold.sys.atom.count) {
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (amino_acid_present) {
        RepresentationType type = RepresentationType::Cartoon;
        ColorMapping color = ColorMapping::SecondaryStructure;

        // Several chains are told apart by color; a short single chain (or none: free amino acids) is shown atom by atom
        size_t num_chains = 0;
        size_t chain_res_count = 0;
        for (size_t i = 0; i < state->mold.sys.instance.count; ++i) {
            if (md_system_instance_entity_kind(&state->mold.sys, i) == MD_ENTITY_KIND_PEPTIDE) {
                if (num_chains++ == 0) chain_res_count = md_instance_comp_count(&state->mold.sys.instance, i);
            }
        }
        if (num_chains > 1) {
            color = ColorMapping::InstId;
        } else if (chain_res_count < 20) {
            type = RepresentationType::BallAndStick;
            color = ColorMapping::Type;
        }

        Representation* prot = create_representation(state, type, color, STR_LIT("protein"));
        snprintf(prot->name, sizeof(prot->name), "protein");
    }
    if (nucleic_present) {
        Representation* nucl = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("nucleic"));
        snprintf(nucl->name, sizeof(nucl->name), "nucleic");
    }
    if (ion_present) {
        Representation* ion = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("ion"));
        snprintf(ion->name, sizeof(ion->name), "ion");
    }
    if (ligand_present) {
        Representation* ligand = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("not (protein or nucleic or water or ion)"));
        snprintf(ligand->name, sizeof(ligand->name), "ligand");
    }
    if (water_present) {
        Representation* water = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("water"));
        water->scale.x = 0.5f;
        snprintf(water->name, sizeof(water->name), "water");
        water->enabled = false;
        if (!amino_acid_present && !nucleic_present && !ligand_present) {
            water->enabled = true;
        }
    }

done:
    if (electronic_structure_present) {
        Representation* rep = create_representation(state, RepresentationType::ElectronicStructure);
        snprintf(rep->name, sizeof(rep->name), "electronic structure");
        rep->enabled = true;
        electronic_structure_select_available_source(&rep->electronic_structure, state->mold.sys);


		// ONE representation, for the ground state moment. A group each meant a VeloxChem file with
		// transition dipoles opened with a representation per group, and with the electric,
		// magnetic and velocity sets that is three nobody asked for on top of the one they wanted.
		// The ground state is what a dipole means to someone who has not said otherwise; the rest
		// are a representation away, and the index slider covers the excited states within each.
		DipoleGroup groups[16];
		size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), state->mold.sys), ARRAY_SIZE(groups));
        for (size_t i = 0; i < num_groups; ++i) {
            // The group name is the identity here, the same string the path spells.
            if (!str_eq(groups[i].label, STR_LIT("ground_state"))) continue;

            vec3_t vec = {0, 0, 0};
            if (!dipole_moment_read(&vec, nullptr, state->mold.sys, groups[i].key, 0)) continue;
            if (vec3_length(vec) <= 1e-3f) continue;

            Representation* dipole_rep = create_representation(state, RepresentationType::DipoleMoment);
            dipole_rep->dipole.dipole_key   = groups[i].key;
            dipole_rep->dipole.dipole_index = 0;

            // Lower case, like every other auto created representation - "protein", "water",
            // "electronic structure". dipole_label_pretty title cases for menus and tooltips,
            // which is a different job and stays as it is.
            snprintf(dipole_rep->name, sizeof(dipole_rep->name), "dipole moment");
            dipole_rep->enabled = true;
            break;
        }
    }

    recompute_atom_visibility_mask(state);
}

void interpolate_system_state(ApplicationState* app) {
    ASSERT(app);
	const auto& sys  = app->mold.sys;

	size_t num_atoms = app->mold.state.num_atoms;
    const size_t num_frames = run_num_frames(app);
    if (num_atoms == 0 || num_frames == 0) return;

    const int64_t last_frame = MAX(0LL, (int64_t)num_frames - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(app->animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    if (app->animation.interpolation == InterpolationMode::Nearest) {
        if (app->mold.last_interpolated_nearest_frame == nearest_frame) {
            return;
        }
        app->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
    }
    app->mold.last_interpolated_nearest_frame = nearest_frame;

    // The backbone keeps the orientation of its cross sections coherent from one displayed state to the next. That
    // continuity only means something while the structure moves continuously, so a jump (seeking, skipping frames,
    // a new run) starts it over. Playback keeps it whatever the speed.
    {
        const double last = app->mold.last_interpolated_frame;
        const bool playing = app->animation.mode == PlaybackMode::Playing;
        if (last < 0.0 || (!playing && fabs(time - last) > 1.5)) {
            app->mold.dirty_gpu_buffers |= MolBit_ResetBackboneHistory;
        }
        app->mold.last_interpolated_frame = time;
    }

    // This represents the frames that we would like to load into memory for interpolation (worst case).
    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame),
    };

    const size_t num_threads = task_system::pool_num_threads();

    // The number of atoms to be processed per thread when divided into chunks
    const uint32_t grain_size = 1024;

    md_allocator_i* temp_arena = app->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(temp_arena);
    defer { md_temp_end(temp); };

    struct Payload {
        ApplicationState* app;
        float s;
        float t;
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];

        md_system_state_t* src_states[4];
		md_system_state_t* dst_state;

        // The backbone of the destination state, and the cartoon's weights for it (temporary, uploaded at the end)
        md_backbone_angles_t*        dst_angle;
        md_secondary_structure_t*    dst_ss;
        md_gl_secondary_structure_t* ss_weights;

        vec3_t* aabb_min;
        vec3_t* aabb_max;

        mat4_t recenter_transform;
    };

    const InterpolationMode mode = (frames[1] == frames[2]) ? InterpolationMode::Nearest : app->animation.interpolation;

    Payload payload = {
        .app = app,
        // Tangent scale of the cardinal spline: (p2 - p0) * s. Tension 0 gives s = 0.5, i.e. Catmull-Rom, which plays
        // uniform motion back uniformly; tension 1 eases in and out of every frame.
        .s = 0.5f * (1.0f - CLAMP(app->animation.tension, 0.0f, 1.0f)),
        .t = (float)fract(time),
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .dst_state = &app->mold.state,
        .aabb_min = md_temp_alloc_array(temp, vec3_t, num_threads),
        .aabb_max = md_temp_alloc_array(temp, vec3_t, num_threads),
    };

    // The backbone of the displayed state is the run's at the frames around it: the angles and the secondary structure
    // go into the state's attributes, the cartoon's weights into temporary memory, which is uploaded at the end.
    const size_t num_segments = app->mold.sys.protein_backbone.segment.count;
    if (num_segments > 0 && app->trajectory_data.backbone_angles.data) {
        payload.dst_angle = md_util_state_backbone_angles_write(&app->mold.state, &app->mold.sys);
    }
    if (num_segments > 0 && app->trajectory_data.secondary_structure.data) {
        payload.dst_ss = md_util_state_secondary_structure_write(&app->mold.state, &app->mold.sys);
        payload.ss_weights = payload.dst_ss ? md_temp_alloc_array(temp, md_gl_secondary_structure_t, num_segments) : nullptr;
    }

    // Stamp the destination with the frame it is about to represent. The interpolated state is not
    // written by a run extraction, so nothing else would fill this in, and a
    // stale value is worse than an absent one. Nearest snaps to a whole frame; the other modes land
    // between two, which is exactly what the fractional part is for.
    app->mold.state.frame = (mode == InterpolationMode::Nearest) ? (double)nearest_frame : time;

    // Fresh coordinates, in the lattice frame of the cell: the turn the previous ones carried goes
    // with them. The System State Changed handler puts one back if the orientation is kept.
    app->operations.state_rotation = mat4_ident();

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
            requested_frames[num_requested_frames++] = (int)frames[0];
            requested_frames[num_requested_frames++] = (int)frames[1];
            requested_frames[num_requested_frames++] = (int)frames[2];
            requested_frames[num_requested_frames++] = (int)frames[3];
            break;
        default:
            ASSERT(false);
            break;
    }

    // Represents the frame cache slot indices for the requested frames, -1 if not present in cache
    int frame_cache_slot_idx[4] = { -1, -1, -1, -1 };

    // Array of frame_cache slots which requires a load
    int frame_cache_load_slot[4] = {0};
    int num_frames_to_load = 0;

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = -1;
        if (!find_frame_in_cache(&slot_idx, requested_frames[i], &app->mold.frame_cache)) {
            slot_idx = find_lru_cache_slot(&app->mold.frame_cache);
            frame_cache_load_slot[num_frames_to_load++] = slot_idx;
        }
        set_mru_cache_slot(&app->mold.frame_cache, slot_idx);
        app->mold.frame_cache.frame_idx[slot_idx] = requested_frames[i];
        frame_cache_slot_idx[i] = slot_idx;
    }

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = frame_cache_slot_idx[i];
        payload.src_states[i] = &app->mold.frame_cache.states[slot_idx];
    }

    // This holds the chain of tasks we are about to submit
    task_system::ID tasks[16] = {0};
    int num_tasks = 0;

    task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frames"), num_frames_to_load,
        [data = &payload, frame_cache_load_slot](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int slot_idx = frame_cache_load_slot[i];
                int frame_idx = data->app->mold.frame_cache.frame_idx[slot_idx];
                md_system_state_t* state = &data->app->mold.frame_cache.states[slot_idx];
                // This thread's own context, kept for as long as the trajectory, so playback opens
                // the run's files once rather than once per frame.
                md_system_extract_t** ex = (thread_num < ARRAY_SIZE(data->app->mold.frame_extract)) ? &data->app->mold.frame_extract[thread_num] : nullptr;
                if (ex && !*ex) {
                    *ex = md_system_extract_begin(&data->app->mold.sys, str_from_cstr(data->app->mold.run),
                        frame_extract_paths, ARRAY_SIZE(frame_extract_paths), md_get_heap_allocator());
                }
                if (ex && *ex) {
                    md_system_extract_frame(*ex, frame_idx, state);
                } else {
                    extract_frame(data->app, frame_idx, state);
                }
            }
        }
    );

    tasks[num_tasks++] = load_task;

    switch (mode) {
        case InterpolationMode::Nearest: {
            task_system::ID interp_task = task_system::create_pool_task(STR_LIT("## Interpolate"), [data = &payload]() {
                data->dst_state->unitcell = data->src_states[0]->unitcell;
                MEMCPY(data->dst_state->xyz, data->src_states[0]->xyz, sizeof(vec3_t) * data->dst_state->num_atoms);
            });
            tasks[num_tasks++] = interp_task;
            break;
        }
        case InterpolationMode::Linear: {
            task_system::ID iterp_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unitcell"), [data = &payload]() {
                double x  = lerp(data->src_states[0]->unitcell.x,  data->src_states[1]->unitcell.x,  data->t);
                double y  = lerp(data->src_states[0]->unitcell.y,  data->src_states[1]->unitcell.y,  data->t);
                double z  = lerp(data->src_states[0]->unitcell.z,  data->src_states[1]->unitcell.z,  data->t);
                double xy = lerp(data->src_states[0]->unitcell.xy, data->src_states[1]->unitcell.xy, data->t);
                double xz = lerp(data->src_states[0]->unitcell.xz, data->src_states[1]->unitcell.xz, data->t);
                double yz = lerp(data->src_states[0]->unitcell.yz, data->src_states[1]->unitcell.yz, data->t);
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
			});

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                vec3_t* dst = data->dst_state->xyz + range_beg;
                const vec3_t* src[2] = { data->src_states[0]->xyz + range_beg, data->src_states[1]->xyz + range_beg };

                md_util_interpolate_linear(dst, src, count, &data->dst_state->unitcell, data->t);
            }, grain_size);

			tasks[num_tasks++] = iterp_cell_task;
            tasks[num_tasks++] = interp_coord_task;

            break;
        }
        case InterpolationMode::CubicSpline: {
            task_system::ID iterp_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unitcell"), [data = &payload]() {
                double x  = cubic_spline(data->src_states[0]->unitcell.x,  data->src_states[1]->unitcell.x,  data->src_states[2]->unitcell.x,  data->src_states[3]->unitcell.x,  data->t, data->s);
                double y  = cubic_spline(data->src_states[0]->unitcell.y,  data->src_states[1]->unitcell.y,  data->src_states[2]->unitcell.y,  data->src_states[3]->unitcell.y,  data->t, data->s);
                double z  = cubic_spline(data->src_states[0]->unitcell.z,  data->src_states[1]->unitcell.z,  data->src_states[2]->unitcell.z,  data->src_states[3]->unitcell.z,  data->t, data->s);
                double xy = cubic_spline(data->src_states[0]->unitcell.xy, data->src_states[1]->unitcell.xy, data->src_states[2]->unitcell.xy, data->src_states[3]->unitcell.xy, data->t, data->s);
                double xz = cubic_spline(data->src_states[0]->unitcell.xz, data->src_states[1]->unitcell.xz, data->src_states[2]->unitcell.xz, data->src_states[3]->unitcell.xz, data->t, data->s);
                double yz = cubic_spline(data->src_states[0]->unitcell.yz, data->src_states[1]->unitcell.yz, data->src_states[2]->unitcell.yz, data->src_states[3]->unitcell.yz, data->t, data->s);
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                vec3_t* dst = data->dst_state->xyz + range_beg;
                const vec3_t* src[4] = { data->src_states[0]->xyz + range_beg, data->src_states[1]->xyz + range_beg, data->src_states[2]->xyz + range_beg, data->src_states[3]->xyz + range_beg };

                md_util_interpolate_cubic_spline(dst, src, count, &data->dst_state->unitcell, data->t, data->s);
            }, grain_size);

            tasks[num_tasks++] = iterp_cell_task;
            tasks[num_tasks++] = interp_coord_task;
            
            break;
        }
        default:
            ASSERT(false);
            break;
    }

    {
        // Calculate a global AABB for the molecule
        task_system::ID aabb_task = task_system::create_pool_task(STR_LIT("## Compute AABB"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            size_t range_len = range_end - range_beg;
            const vec3_t* xyz = data->dst_state->xyz + range_beg;

            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };
            float* r = md_temp_alloc_array(temp, float, range_len);
            md_atom_extract_radii(r, range_beg, range_len, &data->app->mold.sys.atom);

            vec3_t aabb_min = vec3_set1(FLT_MAX);
            vec3_t aabb_max = vec3_set1(-FLT_MAX);
            md_util_aabb_compute(aabb_min.elem, aabb_max.elem, xyz, r, 0, range_len);

            data->aabb_min[thread_num] = aabb_min;
            data->aabb_max[thread_num] = aabb_max;
        });
        tasks[num_tasks++] = aabb_task;
    }

    if (payload.dst_angle) {
        switch (mode) {
            case InterpolationMode::Nearest: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), [data = &payload]() {
                    const md_backbone_angles_t* src_angles[2] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    MEMCPY(data->dst_angle, src_angle, data->app->mold.sys.protein_backbone.segment.count * sizeof(md_backbone_angles_t));
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[2] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[2] = {src_angles[0][i].phi, src_angles[1][i].phi};
                        float psi[2] = {src_angles[0][i].psi, src_angles[1][i].psi};

                        phi[1] = deperiodize_orthof(phi[1], phi[0], (float)TWO_PI);
                        psi[1] = deperiodize_orthof(psi[1], psi[0], (float)TWO_PI);

                        float final_phi = lerp(phi[0], phi[1], data->t);
                        float final_psi = lerp(psi[0], psi[1], data->t);
                        data->dst_angle[i] = {deperiodize_orthof(final_phi, 0, (float)TWO_PI), deperiodize_orthof(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::CubicSpline: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Interpolate Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[4] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[0],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[4] = {src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi};
                        float psi[4] = {src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi};

                        phi[0] = deperiodize_orthof(phi[0], phi[1], (float)TWO_PI);
                        phi[2] = deperiodize_orthof(phi[2], phi[1], (float)TWO_PI);
                        phi[3] = deperiodize_orthof(phi[3], phi[2], (float)TWO_PI);

                        psi[0] = deperiodize_orthof(psi[0], psi[1], (float)TWO_PI);
                        psi[2] = deperiodize_orthof(psi[2], psi[1], (float)TWO_PI);
                        psi[3] = deperiodize_orthof(psi[3], psi[2], (float)TWO_PI);

                        float final_phi = cubic_spline(phi[0], phi[1], phi[2], phi[3], data->t, data->s);
                        float final_psi = cubic_spline(psi[0], psi[1], psi[2], psi[3], data->t, data->s);
                        data->dst_angle[i] = {deperiodize_orthof(final_phi, 0, (float)TWO_PI), deperiodize_orthof(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            default:
                ASSERT(false);
                break;
        }
    }

    if (payload.dst_ss) {
        task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)num_segments, [data = &payload, mode](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            // The state carries the secondary structure of the nearest frame as it was assigned. The cartoon blends the
            // denoised copy (a presentation smoothing, see secondary_structure_render) where it exists.
            const md_secondary_structure_t* label_nearest = data->app->trajectory_data.secondary_structure.data +
                data->app->trajectory_data.secondary_structure.stride * (data->t < 0.5f ? data->frames[1] : data->frames[2]);
            const md_secondary_structure_t* ss_data = data->app->trajectory_data.secondary_structure_render.data ?
                data->app->trajectory_data.secondary_structure_render.data :
                data->app->trajectory_data.secondary_structure.data;
            const size_t ss_stride = data->app->trajectory_data.secondary_structure_render.data ?
                data->app->trajectory_data.secondary_structure_render.stride :
                data->app->trajectory_data.secondary_structure.stride;
            const md_secondary_structure_t* src_ss[4] = {
                ss_data + ss_stride * data->frames[0],
                ss_data + ss_stride * data->frames[1],
                ss_data + ss_stride * data->frames[2],
                ss_data + ss_stride * data->frames[3],
            };
            const md_secondary_structure_t* src_ss_nearest = data->t < 0.5f ? src_ss[1] : src_ss[2];
            auto normalize_ss = [](md_gl_secondary_structure_t ss) {
                ss.helix = CLAMP(ss.helix, 0.0f, 1.0f);
                ss.sheet = CLAMP(ss.sheet, 0.0f, 1.0f);
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
                t = CLAMP(t, 0.0f, 1.0f);
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
            case InterpolationMode::CubicSpline: {
                
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
        tasks[num_tasks++] = ss_task;

        // Isolated coils between matching structured segments are filled in, to reduce the noise during transitions.
        // A non temporal filtering step, after the blend above.
        task_system::ID ss_cleanup_task = task_system::create_pool_task(STR_LIT("## Cleanup Secondary Structures"), [data = &payload]() {
            secondary_structure_weights_fill_isolated_coils(data->ss_weights, &data->app->mold.sys.protein_backbone);
        });
        tasks[num_tasks++] = ss_cleanup_task;
    }

    if (num_tasks > 0) {
        for (int i = 1; i < num_tasks; ++i) {
            task_system::set_task_dependency(tasks[i], tasks[i-1]);
        }
        task_system::enqueue_task(tasks[0]);
        task_system::task_wait_for(tasks[num_tasks - 1]);
    }

    // The weights are the renderer's input and nothing else's: straight to the GPU, from the temporary memory above
    if (payload.ss_weights) {
        md_gl_mol_set_backbone_secondary_structure(app->mold.gl_mol, 0, (uint32_t)num_segments, payload.ss_weights, 0);
    }

    vec3_t aabb_min = payload.aabb_min[0];
    vec3_t aabb_max = payload.aabb_max[0];
    for (size_t i = 1; i < task_system::pool_num_threads(); ++i) {
        aabb_min = vec3_min(aabb_min, payload.aabb_min[i]);
        aabb_max = vec3_max(aabb_max, payload.aabb_max[i]);
    }
    app->mold.sys_aabb_min = aabb_min;
    app->mold.sys_aabb_max = aabb_max;

    // unitcell transform is essentially just a translation to place the center of the unitcell at the origin
    mat3_t A;
    md_unitcell_A_extract_float(A.elem, &app->mold.state.unitcell);
    vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
    app->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

#if 0
    if (sys.unitcell.flags) {
        vec3_t c = sys.unitcell.basis * vec3_set1(0.5f);
        app->mold.model_mat = mat4_translate_vec3(-c);
    }
#endif

    app->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
}

// Exact, element by element. The identities compared against here are only ever assigned, never computed.
static bool is_identity(const mat4_t& M) {
    const mat4_t I = mat4_ident();
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            if (M.elem[c][r] != I.elem[c][r]) return false;
        }
    }
    return true;
}

void recenter_mark_query_dirty(ApplicationState* state) {
    ASSERT(state);
    state->operations.recenter_query.version += 1;
    if (state->operations.recenter_query.version == 0) {
        state->operations.recenter_query.version = 1;
    }
}

const md_bitfield_t& recenter_get_active_target_mask(const ApplicationState* state) {
    ASSERT(state);
    return state->operations.recenter_query.enabled ? state->operations.recenter_query.mask : state->operations.selection_mask;
}

bool recenter_update_query_mask(ApplicationState* state) {
    ASSERT(state);

    auto& query = state->operations.recenter_query;
    const uint64_t ir_fingerprint = state->script.ir ? state->script.ir_fingerprint : 0;
    if (query.ir_fingerprint != ir_fingerprint) {
        query.ir_fingerprint = ir_fingerprint;
        recenter_mark_query_dirty(state);
    }

    if (!query.enabled) {
        return false;
    }

    if (query.dynamic) {
        recenter_mark_query_dirty(state);
    }

    if (query.evaluated_version == query.version) {
        return false;
    }

    md_bitfield_clear(&query.mask);
    query.dynamic = false;
    query.valid = md_filter(&query.mask, str_from_cstr(query.query), &state->mold.sys, &state->mold.state, state->script.ir, &query.dynamic, query.error, sizeof(query.error));
    query.evaluated_version = query.version;
    return true;
}

void recenter_update(ApplicationState* state) {
    ASSERT(state);
    recenter_update_query_mask(state);
    recenter_update_target_data(state);
}

// Exactly the same atoms, whatever range of bits each bitfield happens to be stored over. A hash of the
// stored blocks is neither: it misses their offset, and the same atoms can be stored over different ranges.
static bool same_atoms(const md_bitfield_t* a, const md_bitfield_t* b) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_bitfield_t diff = md_bitfield_create(md_temp_allocator(temp));
    md_bitfield_xor(&diff, a, b);
    return md_bitfield_empty(&diff);
}

void recenter_update_target_data(ApplicationState* state) {
    if (run_num_frames(state) == 0) return;

    // The reference is only for keeping the orientation. Centering alone needs nothing from the first frame,
    // so without this a target picked by a query which depends on the frame read it from disk in every frame.
    if (!state->operations.fixate_orientation) return;

    auto& ref = state->operations.initial_frame;
    const md_bitfield_t& target_mask = recenter_get_active_target_mask(state);

    // Keyed on the atoms themselves. The versions this used to compare change every time a query is
    // evaluated, which for one that depends on the frame is every frame, whether or not its atoms changed.
    if (!ref.valid || !same_atoms(&ref.target_mask, &target_mask)) {
        if (!ref.target_mask.alloc) {
            md_bitfield_init(&ref.target_mask, state->allocator.persistent);
        }
        md_bitfield_copy(&ref.target_mask, &target_mask);
        ref.valid = true;
        size_t count = md_bitfield_popcount(&target_mask);

        md_array_resize(state->operations.initial_frame.rel_xyzw, count, state->allocator.persistent);

        state->operations.initial_frame.com = vec3_zero();

        if (state->operations.initial_frame.rel_xyzw && count > 0) {
            // Fetch initial frame data required for orienting the structure
            size_t num_atoms = state->mold.sys.atom.count;
            vec3_t* temp_xyz = (vec3_t*)md_vm_arena_push(state->allocator.frame, sizeof(vec3_t) * ALIGN_TO(num_atoms, 16));

            md_system_state_t temp_state = { 0 };
            temp_state.num_atoms = num_atoms;
            temp_state.xyz = temp_xyz;
            extract_frame(state, 0, &temp_state);

            md_bitfield_iter_t it = md_bitfield_iter_create(&target_mask);
            int dst_idx = 0;
            while (md_bitfield_iter_next(&it)) {
                uint64_t src_idx = md_bitfield_iter_idx(&it);
                float mass = md_atom_mass(&state->mold.sys.atom, src_idx);
                state->operations.initial_frame.rel_xyzw[dst_idx++] = vec4_from_vec3(temp_xyz[src_idx], mass);
            }

            // Mutually consistent images and the plain weighted mean of them, then relative to that mean.
            // The circular mean (md_util_com_compute_vec4 with a cell) is not the mean of the placed
            // points, and a fit against coordinates which are not centred is biased.
            vec3_t com = vec3_zero();
            md_util_deperiodize_self_vec4(state->operations.initial_frame.rel_xyzw, count, &temp_state.unitcell, &com);
            const vec4_t com4 = vec4_from_vec3(com, 0);
            for (size_t i = 0; i < count; ++i) {
                state->operations.initial_frame.rel_xyzw[i] = vec4_sub(state->operations.initial_frame.rel_xyzw[i], com4);
            }
            state->operations.initial_frame.com = com;
        }
    }
}

bool recenter_calculate_transform(mat4_t* translation, mat4_t* rotation, const ApplicationState* app) {
    ASSERT(translation);
    ASSERT(rotation);
    ASSERT(app);

    *translation = mat4_ident();
    *rotation    = mat4_ident();
    bool turn = false;

    const md_bitfield_t& target_mask = recenter_get_active_target_mask(app);
    size_t count = md_bitfield_popcount(&target_mask);

    if (count > 0) {
        md_temp_scope_t temp = md_temp_begin_in(app->allocator.frame);
        defer { md_temp_end(temp); };

        // Extract xyzw subset of target
        vec4_t* target_xyzw = md_temp_alloc_array(temp, vec4_t, count);

		md_util_system_extract_xyzw_from_mask(target_xyzw, &target_mask, &app->mold.sys, &app->mold.state);

        // Calculate target
        vec3_t target = {0};
        if (md_unitcell_flags(&app->mold.state.unitcell) != 0) {
            mat3_t A = {0};
            md_unitcell_A_extract_float(A.elem, &app->mold.state.unitcell);
            target = mat3_mul_vec3(A, vec3_set1(0.5f));
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
        mat3_t R = mat3_ident();
        vec3_t target_com = vec3_zero();

        // The reference has to have been built from the SAME target that is being fitted now.
        // A size match is not sufficient: the selection can change to a different set of equal
        // size between recenter_update() and this call, which would silently pair up unrelated
        // atoms and yield a garbage rotation. The atoms themselves are the identity of the target.
        const bool reference_valid =
            app->operations.fixate_orientation &&
            app->operations.initial_frame.valid &&
            app->operations.initial_frame.rel_xyzw &&
            md_array_size(app->operations.initial_frame.rel_xyzw) == count &&
            same_atoms(&app->operations.initial_frame.target_mask, &target_mask);

        // R maps the CURRENT target onto the reference: R * (q - target_com) ~= p. The relative fit this
        // replaced had its operands the other way around, which yields the rotation carrying the reference
        // onto the current frame - applied to the current frame, it doubled the rotation it was meant to
        // cancel. It also centred on the circular mean, which lands in the reference cell whatever image
        // the target occupies and is not the mean of the placed points.
        if (app->operations.fixate_orientation && reference_valid) {
            // The reference is stored relative to its own centre, so its centre here is the origin
            md_util_optimal_rotation_pbc_vec4_iter(&R, &target_com, target_xyzw, app->operations.initial_frame.rel_xyzw, vec3_zero(),
                                                   target_xyzw, count, &app->mold.state.unitcell, 8, 1.0e-6f);
            R = mat3_orthonormalize(R);
        } else {
            md_util_deperiodize_self_vec4(target_xyzw, count, &app->mold.state.unitcell, &target_com);
        }

        // Split at the cell centre. The translation preserves the lattice, so it is valid for coordinates
        // in whatever image they arrived in. The turn does not: R times a lattice vector is not a lattice
        // vector, so it may only be applied once every image is settled relative to the target. The
        // product is the transform this used to return in one piece,
        //     translate(target) * A * R * translate(-target_com)
        const mat4_t A = app->operations.alignment_mat;
        turn = (app->operations.fixate_orientation && reference_valid) || !is_identity(A);

        *translation = mat4_translate_vec3(vec3_sub(target, target_com));
        if (turn) {
            *rotation = mat4_translate_vec3(target) * A * mat4_from_mat3(R) * mat4_translate_vec3(-target);
        }
    }
    return turn;
}

// Every atom of mold.state through M
static void state_transform(ApplicationState* app, const mat4_t& M) {
    md_system_state_t& s = app->mold.state;
    task_system::ID task = task_system::create_pool_task(STR_LIT("## Transform"), (uint32_t)s.num_atoms, [&s, M](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
        (void)thread_num;
        mat4_batch_transform_inplace(s.xyz + range_beg, 1.0f, range_end - range_beg, M);
    }, 1024);
    task_system::enqueue_task(task);
    task_system::task_wait_for(task);
}

bool apply_state_operations(ApplicationState* app, bool recenter, bool pbc, bool unwrap) {
    ASSERT(app);
    md_system_t& sys = app->mold.sys;
    md_system_state_t& s = app->mold.state;

    if (s.num_atoms == 0 || !s.xyz) return false;
    if (!recenter && !pbc && !unwrap) return false;

    const bool periodic = md_unitcell_flags(&s.unitcell) != 0;

    // Into the lattice frame first. Wrapping or making whole turned coordinates against the unturned
    // cell moves atoms by vectors that are not lattice vectors, to places where they have no image.
    const mat4_t prev_rotation = app->operations.state_rotation;
    if (!is_identity(prev_rotation)) {
        state_transform(app, mat4_inverse(prev_rotation));
    }

    // A wrap or a make whole on its own keeps the turn the coordinates had
    mat4_t rotation = prev_rotation;
    bool fresh_turn = false;
    if (recenter && !md_bitfield_empty(&recenter_get_active_target_mask(app))) {
        // Current for this target before it is fitted against. The main loop does this every frame too, but
        // ticking 'keep orientation' applies at once, before the loop has had a chance to build it.
        recenter_update_target_data(app);

        // Measured on the lattice frame coordinates: the fit is against the reference, not against a
        // previous fit's output
        mat4_t translation = mat4_ident();
        fresh_turn = recenter_calculate_transform(&translation, &rotation, app);
        state_transform(app, translation);
    }
    const bool turn = !is_identity(rotation);

    // Settle the periodic images before turning: every atom into the image nearest the cell centre,
    // which is where the translation just put the target. Without it the turn carries the image each
    // atom happened to be written in by the trajectory into the result, as a displacement of R times
    // a lattice vector. That differs from frame to frame for every atom crossing the trajectory's own
    // box, and for the whole system at once whenever the target's centre lands in another image.
    if (periodic && (pbc || fresh_turn)) {
        task_system::ID task = task_system::create_pool_task(STR_LIT("## Apply PBC"), (uint32_t)s.num_atoms, [&s](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            md_util_pbc(s.xyz + range_beg, NULL, range_end - range_beg, &s.unitcell);
        });
        task_system::enqueue_task(task);
        task_system::task_wait_for(task);
    }

    if (unwrap) {
        const size_t num_structures = md_structure_count(&sys.structure);
        if (num_structures > 0) {
            task_system::ID task = task_system::create_pool_task(STR_LIT("## Unwrap Structures"), (uint32_t)num_structures, [&s, &sys](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    md_structure_t structure = {};
                    md_structure_extract(&structure, &sys.structure, i);
                    md_util_unwrap_structure(&s, &structure);
                }
            });
            task_system::enqueue_task(task);
            task_system::task_wait_for(task);
        }
    }

    if (turn) {
        state_transform(app, rotation);
    }
    app->operations.state_rotation = turn ? rotation : mat4_ident();

    return true;
}

bool recompute_covalent_bonds(ApplicationState* app, int64_t frame) {
    ASSERT(app);
    md_system_t& sys = app->mold.sys;
    const size_t num_atoms = sys.atom.count;
    if (num_atoms == 0) return false;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_system_state_t state = { .alloc = temp.arena };
    if (!md_system_state_init(&state, num_atoms)) return false;

    // The bond search takes its minimum images in the cell, so its coordinates have to be in the cell's
    // lattice frame. A frame of the run always is.
    const size_t num_frames = run_num_frames(app);
    if (num_frames > 0) {
        if (!extract_frame(app, CLAMP(frame, (int64_t)0, (int64_t)num_frames - 1), &state)) {
            MD_LOG_ERROR("Failed to extract frame data");
            return false;
        }
    } else {
        // The coordinates shown, taken back out of any turn from keeping the orientation
        const md_system_state_t& shown = app->mold.state;
        if (!shown.xyz || shown.num_atoms != num_atoms) return false;
        MEMCPY(state.xyz, shown.xyz, num_atoms * sizeof(vec3_t));
        state.unitcell = shown.unitcell;
        if (!is_identity(app->operations.state_rotation)) {
            mat4_batch_transform_inplace(state.xyz, 1.0f, num_atoms, mat4_inverse(app->operations.state_rotation));
        }
    }

    MD_LOG_DEBUG("RECALCULATING BONDS");
    md_util_infer_covalent_bonds(&sys.bond, &state, &sys, sys.alloc);
    md_bond_build_connectivity(&sys.bond, num_atoms, sys.alloc);
    app->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
    return true;
}

bool picking_range_reserve(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count, uint64_t key) {
    ASSERT(space);

    if (count > 0 && space->num_ranges < ARRAY_SIZE(space->ranges)) {
        PickingRange* curr_range = &space->ranges[space->num_ranges++];
        PickingRange* prev_range = space->num_ranges > 1 ? &space->ranges[space->num_ranges - 2] : NULL;
        curr_range->domain = domain;
        curr_range->beg = prev_range ? prev_range->end : 0;
        curr_range->end = curr_range->beg + (uint32_t)count;
        curr_range->key = key;
        if (out_range) {
            MEMCPY(out_range, curr_range, sizeof(PickingRange));
        }
        return true;
    }

    return false;
}

void picking_handler_new_frame(PickingHandler* handler) {
    ASSERT(handler);
    handler->frame_idx += 1;

    const uint32_t slot_idx = handler->frame_idx % ARRAY_SIZE(handler->history);
    handler->history[slot_idx].submitted_frame_idx = handler->frame_idx;
    handler->history[slot_idx].space = PickingSpace{};
}

PickingSpace* picking_handler_current_space(PickingHandler* handler) {
    ASSERT(handler);
    return &handler->history[handler->frame_idx % ARRAY_SIZE(handler->history)].space;
}

const PickingSpace* picking_handler_find_space(const PickingHandler& handler, uint32_t submitted_frame_idx) {
    for (size_t i = 0; i < ARRAY_SIZE(handler.history); ++i) {
        const auto& hist = handler.history[i];
        if (hist.submitted_frame_idx == submitted_frame_idx) {
            return &hist.space;
        }
    }

    return nullptr;
}

void picking_surface_init(PickingSurface* surface, PickingSourceID source) {
    ASSERT(surface);

    *surface = PickingSurface{};
    surface->source = source;

    for (size_t i = 0; i < ARRAY_SIZE(surface->slots); ++i) {
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
    ASSERT(surface);

    for (size_t i = 0; i < ARRAY_SIZE(surface->slots); ++i) {
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
    vec2_t surface_coord,
    vec2_t screen_coord,
    const mat4_t& clip_to_world
) {
    ASSERT(surface);

    if (!fbo || width == 0 || height == 0) {
        return false;
    }

    const int x = (int)surface_coord.x;
    const int y = (int)surface_coord.y;
    if (x < 0 || y < 0 || x >= (int)width || y >= (int)height) {
        return false;
    }

    const uint32_t queue_idx = surface->slot_cursor % ARRAY_SIZE(surface->slots);
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

const PickingRange* picking_space_find_range(const PickingSpace& space, PickingDomainID domain, uint64_t key) {
    for (size_t i = 0; i < space.num_ranges; ++i) {
        if (space.ranges[i].domain == domain && space.ranges[i].key == key) {
            return &space.ranges[i];
        }
    }

    return nullptr;
}

static const PickingRange* find_picking_range(const PickingSpace* space, uint32_t raw_idx) {
    ASSERT(space);

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
    ASSERT(out_hit);
    ASSERT(surface);

    *out_hit = {};

    if (surface->slot_cursor < ARRAY_SIZE(surface->slots)) {
        return false;
    }

    const uint32_t read_idx = surface->slot_cursor % ARRAY_SIZE(surface->slots);
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
    out_hit->frame_idx = slot.submitted_frame_idx;
    out_hit->raw_idx = raw_idx;
    out_hit->local_idx = raw_idx - range.beg;
    out_hit->surface_coord = slot.surface_coord;
    out_hit->screen_coord = slot.screen_coord;
    out_hit->depth = depth;

    const vec4_t viewport = {0, 0, (float)slot.viewport_width, (float)slot.viewport_height};
    out_hit->world_pos = mat4_unproject({slot.surface_coord.x, slot.surface_coord.y, depth}, slot.clip_to_world, viewport);

    return true;
}

bool picking_surface_submit_readback_and_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler,
    const PickingReadbackRequest& request
) {
    ASSERT(out_hit);
    ASSERT(surface);

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

InteractionSurfaceState interaction_surface(InteractionSurfaceID id, const vec2_t& size, InteractionSurfaceFlags flags) {
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
    ASSERT(draw_list);

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
    ASSERT(out_hit);
    if (state.hovered && args.fbo && args.width && args.height) {
        const ImVec2 local_coord = ImVec2(state.mouse_local.x, state.surface_size.y - state.mouse_local.y) * ImGui::GetIO().DisplayFramebufferScale;

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
    ASSERT(target);
    InteractionSurfaceViewTransformResult result = {};
    if (state.active || state.hovered) {
        if (state.selection_mode == InteractionSelectionMode::None) {
            const vec2_t delta = vec_cast(ImGui::GetIO().MouseDelta);
            const vec2_t coord = state.mouse_local;
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
    ASSERT(event);
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
        } else {
            event->kind = InteractionSurfaceEventKind::Hover;
        }
    }
}

void point_set_region_mask_compute(md_bitfield_t* mask,
    const vec3_t xyz[],
    size_t count,
    const md_bitfield_t* candidate_mask,
    const mat4_t& world_to_clip,
    const vec2_t& region_min,
    const vec2_t& region_max,
    const vec2_t& surface_size)
{
    ASSERT(mask);
    ASSERT(xyz);

    md_bitfield_clear(mask);

    // Transform visible atoms using supplied world_to_clip and set if within region

    if (candidate_mask) {
        // Do for candidate set
        md_bitfield_iter_t it = md_bitfield_iter_create(candidate_mask);

        while (md_bitfield_iter_next(&it)) {
            size_t idx = md_bitfield_iter_idx(&it);
            vec4_t xyz1 = vec4_from_vec3(xyz[idx], 1.0f);
            vec4_t coord = mat4_mul_vec4(world_to_clip, xyz1);
            vec2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitfield_set_bit(mask, idx);
            }
        }
    } else {
        // Do for full set
        for (size_t i = 0; i < count; ++i) {
            vec4_t xyz1 = vec4_from_vec3(xyz[i], 1.0f);
            vec4_t coord = mat4_mul_vec4(world_to_clip, xyz1);
            vec2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitfield_set_bit(mask, i);
            }
        }
    }
}

bool file_queue_empty(const FileQueue* queue) {
    return queue->head == queue->tail;
}

bool file_queue_full(const FileQueue* queue) {
    return (queue->head + 1) % ARRAY_SIZE(queue->arr) == queue->tail;
}

void file_queue_push(FileQueue* queue, str_t path, FileFlags flags) {
    ASSERT(queue);
    ASSERT(!file_queue_full(queue));
    int prio = 5;


    str_t ext;
    if (extract_ext(&ext, path)) {
        LoaderType type = loader::type_from_ext(ext);
        LoaderFlags loader_flags = loader::type_flags(type);
        if (str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
            prio = 1;
        } else if (loader_flags & LoaderFlag_System) {
            prio = 2;
        } else if (loader_flags & LoaderFlag_Temporal) {
            // Joins the trajectory's run, so it goes after the trajectory when both are dropped
            prio = 4;
        } else if (loader_flags & (LoaderFlag_Trajectory | LoaderFlag_Supplemental)) {
            // A supplemental file (a topology) needs the system in place, same as a trajectory
            prio = 3;
        } else {
            flags |= FileFlags_ShowDialogue;
        }
    } else {
        // Unknown extension
        flags |= FileFlags_ShowDialogue;
    }

    uint32_t i = queue->head;
    queue->arr[queue->head] = {str_copy(path, queue->ring), flags, prio};
    queue->head = (queue->head + 1) % ARRAY_SIZE(queue->arr);

    // Sort queue based on prio
     while (i != queue->tail && queue->arr[i].prio < queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)].prio) {
        FileQueue::Entry tmp = queue->arr[i];
        queue->arr[i] = queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)];
        queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)] = tmp;
        i = (i - 1) % ARRAY_SIZE(queue->arr);
     }
}

FileQueue::Entry file_queue_front(const FileQueue* queue) {
    ASSERT(!file_queue_empty(queue));
    return queue->arr[queue->tail];
}

FileQueue::Entry file_queue_pop(FileQueue* queue) {
    ASSERT(queue);
    ASSERT(!file_queue_empty(queue));
    FileQueue::Entry front = file_queue_front(queue);
    queue->tail = (queue->tail + 1) % ARRAY_SIZE(queue->arr);
    return front;
}

void file_queue_process(ApplicationState* state) {
    ASSERT(state);
    if (!file_queue_empty(&state->file_queue) && !state->load_dataset.show_window) {
        FileQueue::Entry e = file_queue_pop(&state->file_queue);

        str_t ext;
        extract_ext(&ext, e.path);

        if (str_eq_ignore_case(ext, WORKSPACE_FILE_EXTENSION)) {
            load_workspace(state, e.path);
        } else {
            loader::LoaderState loader_state = {};
            loader::init(&loader_state, e.path, &state->mold.sys);
                
            if ((e.flags & FileFlags_ShowDialogue) || (loader_state.flags & LoaderFlag_RequiresDialogue)) {
                state->load_dataset = LoadDatasetWindowState();
                str_copy_to_char_buf(state->load_dataset.path_buf, sizeof(state->load_dataset.path_buf), e.path);
                state->load_dataset.path_changed = true;
                state->load_dataset.show_window = true;
                state->load_dataset.coarse_grained = e.flags & FileFlags_CoarseGrained;
            } else {
                loader_state.flags |= (e.flags & FileFlags_DisableCacheWrite) ? LoaderFlag_DisableCacheWrite : 0;
                loader_state.flags |= (e.flags & FileFlags_CoarseGrained) ? LoaderFlag_CoarseGrained : 0;
                if (load_data_from_file(state, e.path, loader_state)) {
                    state->animation = {};
                    // @TODO @FIX: This is hacky, just because the loader CAN set system state does not mean it always will.
                    // This should be instead captured and performed by the Event that signals when a new system is loaded.
                    if (loader_state.flags & LoaderFlag_System) {
                        md_bitfield_reset(&state->representation.visibility_mask);

                        if (!state->settings.keep_representations) {
                            remove_all_representations(state);
                            create_default_representations(state);
                        }
                        recompute_atom_visibility_mask(state);
                        state->mold.interpolate_system_state = true;
                        state->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                        reset_view(&state->view.camera, state->mold.state, &state->representation.visibility_mask);
                    }
                }
            }
        }
    }
}

void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitfield_t* mask) {
    ASSERT(transform);
    if (!state.num_atoms) return;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    // A mask that selects nothing or everything is the whole system
    const int32_t* indices = nullptr;
    size_t count = state.num_atoms;
    if (mask) {
        const size_t popcount = md_bitfield_popcount(mask);
        if (0 < popcount && popcount < state.num_atoms) {
            int32_t* idx = md_temp_alloc_array(temp, int32_t, popcount);
            size_t len = md_bitfield_iter_extract_indices(idx, popcount, md_bitfield_iter_create(mask));
            if (len > popcount || len > state.num_atoms) {
                MD_LOG_DEBUG("Error: Invalid number of indices");
                len = MIN(popcount, state.num_atoms);
            }
            indices = idx;
            count   = len;
        }
    }

    // The world is drawn translated so that the center of the unit cell is at the origin
    mat3_t A = {};
    md_unitcell_A_extract(A.elem, &state.unitcell);
    const bool   has_cell    = (md_unitcell_flags(&state.unitcell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) != 0;
    const vec3_t cell_offset = -mat3_mul_vec3(A, vec3_set1(0.5f));
    const float  fov_y       = Camera().fov_y;

    if (indices && count <= 4) {
        // A handful of atoms has no meaningful shape: center on them, keep the current orientation, fit the distance
        vec3_t aabb_min = {}, aabb_max = {};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, state.xyz, nullptr, indices, count);
        const vec3_t center = (aabb_min + aabb_max) * 0.5f;
        transform->distance = camera_fit_distance(state.xyz, indices, count, center, transform->orientation, fov_y);
        transform->position = center + cell_offset + transform->orientation * vec3_set(0, 0, transform->distance);
        return;
    }

    // See camera_compute_default_view for what it considers a good view of what
    *transform = camera_compute_default_view(state.xyz, state.num_atoms, indices, count, has_cell ? &A : nullptr, fov_y);
    transform->position = transform->position + cell_offset;
}

void ViamdEventHandler::process_events(const viamd::Event* events, size_t num_events) {
    for (size_t i = 0; i < num_events; ++i) {
        const viamd::Event& event = events[i];
        switch (event.type) {
        case viamd::EventType_ViamdFrameTick:
            break;
        case viamd::EventType_ViamdPickingRangeReserve: {
			ASSERT(event.payload_type == viamd::EventPayloadType_PickingSpace);
			PickingSpace* space = (PickingSpace*)event.payload;
            size_t num_atoms = state->mold.sys.atom.count;
            size_t num_bonds = state->mold.sys.bond.count;
            picking_range_reserve(&state->picking_range_atom, space, PickingDomain_Atom, num_atoms);
            picking_range_reserve(&state->picking_range_bond, space, PickingDomain_Bond, num_bonds);
            // One index per segment, in the order md_gl draws them: the cartoon of segment i writes beg + i
            picking_range_reserve(&state->picking_range_backbone, space, PickingDomain_BackboneSegment, state->mold.sys.protein_backbone.segment.count);

            // One range per dipole group, keyed by the group's vector attribute, sized by that
            // attribute's own shape. An index within a range is then the element index inside the
            // group, so a hit carries the whole identity and nobody has to reproduce a flat
            // ordering over every (group, element) pair to speak about one dipole.
            DipoleGroup groups[16];
            const size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), state->mold.sys), ARRAY_SIZE(groups));
            for (size_t j = 0; j < num_groups; ++j) {
                picking_range_reserve(NULL, space, PickingDomain_Dipole, groups[j].count, groups[j].key);
            }
            break;
        }
        case viamd::EventType_ViamdInteractionSurface:
            if (event.payload) {
                ASSERT(event.payload_type == viamd::EventPayloadType_InteractionSurfaceEvent);
                InteractionSurfaceEvent* surf = (InteractionSurfaceEvent*)event.payload;
                switch (surf->kind) {
                case InteractionSurfaceEventKind::Hover:
                    draw_picking_tooltip_window(surf->hit, *state);
					[[fallthrough]];
                case InteractionSurfaceEventKind::Click:
                    // Use highlight mask as intermediate mask for selection operations and to provide hover feedback
                    md_bitfield_clear(&state->selection.highlight_mask);
                    if (surf->hit.domain == PickingDomain_Atom) {
                        int32_t atom_idx = surf->hit.local_idx;
                        if (atom_idx >= 0 && (size_t)atom_idx < state->mold.sys.atom.count) {
                            // Grow the hovered atom by the current granularity so hover feedback matches what a click would select
                            mask_set_atom_by_selection_granularity(&state->selection.highlight_mask, (size_t)atom_idx, state->selection.granularity, state->mold.sys);
                        }
                    } else if (surf->hit.domain == PickingDomain_Bond) {
                        size_t bond_idx = surf->hit.local_idx;
                        if (bond_idx < state->mold.sys.bond.count) {
                            // Grow both bond atoms by the current granularity, same as for atom hits
                            const md_atom_pair_t& pair = state->mold.sys.bond.pairs[bond_idx];
                            mask_set_atom_by_selection_granularity(&state->selection.highlight_mask, (size_t)pair.idx[0], state->selection.granularity, state->mold.sys);
                            mask_set_atom_by_selection_granularity(&state->selection.highlight_mask, (size_t)pair.idx[1], state->selection.granularity, state->mold.sys);
                        }
                    } else if (surf->hit.domain == PickingDomain_BackboneSegment) {
                        const md_protein_backbone_data_t& bb = state->mold.sys.protein_backbone;
                        const size_t seg_idx = surf->hit.local_idx;
                        if (seg_idx < bb.segment.count && bb.segment.comp_idx) {
                            // A cartoon segment stands for the component it stems from: at least that is hovered and selected
                            mask_set_component_by_selection_granularity(&state->selection.highlight_mask, (size_t)bb.segment.comp_idx[seg_idx], state->selection.granularity, state->mold.sys);
                        }
                    }
                    
                    // Commit to selection mask upon click release, for hover we only update the highlight mask
                    if (surf->kind == InteractionSurfaceEventKind::Click) {
                        // The single selection sequence records the order in which individual atoms were picked
                        // and is what the context menu turns into script suggestions. Only a real click advances
                        // it: a Hover event always carries selection_mode None, so doing this above would reset
                        // the sequence to the atom under the cursor on every frame the mouse crosses the molecule.
                        if (surf->hit.domain == PickingDomain_Atom) {
                            int32_t atom_idx = surf->hit.local_idx;
                            if (atom_idx >= 0 && (size_t)atom_idx < state->mold.sys.atom.count) {
                                if (surf->selection_mode == InteractionSelectionMode::Append) {
                                    single_selection_sequence_push_idx(&state->selection.single_selection_sequence, atom_idx);
                                }
                                else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                                    single_selection_sequence_pop_idx(&state->selection.single_selection_sequence, atom_idx);
                                }
                                else if (surf->selection_mode == InteractionSelectionMode::None) {
                                    single_selection_sequence_clear(&state->selection.single_selection_sequence);
                                    single_selection_sequence_push_idx(&state->selection.single_selection_sequence, atom_idx);
                                }
                            }
                        }

                        if (surf->hit.domain == PickingDomain_Atom || surf->hit.domain == PickingDomain_Bond || surf->hit.domain == PickingDomain_BackboneSegment) {
                            if (surf->selection_mode == InteractionSelectionMode::Append) {
                                md_bitfield_or_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                            }
                            else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                                md_bitfield_andnot_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                            }
                            else if (surf->selection_mode == InteractionSelectionMode::None) {
                                md_bitfield_clear(&state->selection.selection_mask);
                                md_bitfield_or_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                            }
                        }
                        else if (surf->hit.domain == 0) {
                            // A plain click (no modifier, no drag) or a remove-click on empty space clears the selection
                            if (surf->selection_mode == InteractionSelectionMode::None || surf->selection_mode == InteractionSelectionMode::Remove) {
                                md_bitfield_clear(&state->selection.selection_mask);
                                single_selection_sequence_clear(&state->selection.single_selection_sequence);
                            }
                        }
                    }
                    break;
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
            ASSERT(event.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
            PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)event.payload;
            fill_picking_tooltip(req, *state, req->hit);
            break;
        }
        case viamd::EventType_ViamdViewFit: {
            ASSERT(event.payload_type == viamd::EventPayloadType_ViewFitRequest);
            ViewFitRequest* req = (ViewFitRequest*)event.payload;
            if (req) {
                md_bitfield_t* bf = nullptr;
                switch (req->round) {
                case ViewFitRound_Highlight:
                    bf = &state->selection.highlight_mask; break;
                case ViewFitRound_Selection:
                    bf = &state->selection.selection_mask; break;
                case ViewFitRound_Visible:
                    bf = &state->representation.visibility_mask; break;
                default:
                    break;
                }

                size_t popcount = bf ? md_bitfield_popcount(bf) : 0;
                if (popcount > 0) {
                    vec4_t* dst_xyzw = md_array_extend(req->xyzw, popcount, req->alloc);
                    if (dst_xyzw) {
                        md_util_system_extract_xyzw_from_mask(dst_xyzw, bf, &state->mold.sys, &state->mold.state);
                    }
                }
            }
            break;
        }
        case viamd::EventType_ViamdSystemStateChanged: {
			// Apply operators if toggled on system change

            // EXTRACT STATE HERE FROM PAYLOAD
            ASSERT(event.payload_type == viamd::EventPayloadType_ApplicationState);
            ApplicationState* app = (ApplicationState*)event.payload;

            int num_tasks = 0;
            task_system::ID tasks[16];
            
            // The event means mold.state was just written from its source, in the lattice frame of its
            // cell: whatever turn the previous coordinates carried went with them.
            app->operations.state_rotation = mat4_ident();

            if (app->operations.recalc_bonds) {
                static int64_t cur_nearest_frame = -1;

                // We cannot recalculate bonds while the evaluation is running
                // because it would overwrite the bond data while we are reading it
                int64_t nearest_frame = (int64_t)(app->animation.frame + 0.5);
                if (!task_system::task_is_running(app->tasks.evaluate) && !task_system::task_is_running(app->script.vis_task)) {
                    const bool has_frames = run_num_frames(app) > 0;
                    if (!has_frames || (cur_nearest_frame != nearest_frame)) {
                        cur_nearest_frame = nearest_frame;
                        // From the whole frame nearest the animation time, not the interpolated state shown
                        task_system::ID recalc_bond_task = task_system::create_pool_task(STR_LIT("## Recalc bond task"), [app, nearest_frame]() {
                            recompute_covalent_bonds(app, nearest_frame);
                        });
                        tasks[num_tasks++] = recalc_bond_task;
                    }
                }
            }

            if (num_tasks > 0) {
                for (int j = 1; j < num_tasks; ++j) {
                    task_system::set_task_dependency(tasks[j], tasks[j-1]);
                }
                task_system::enqueue_task(tasks[0]);
                task_system::task_wait_for(tasks[num_tasks - 1]);
            }

            // Recenter, wrap, make whole and turn, in that order (see apply_state_operations). They rewrite
            // the coordinates that update_md_buffers uploads, and nothing else flags them: the synchronous
            // broadcasts of this event do not pass through the interpolation step that would otherwise
            // have set the bit.
            if (apply_state_operations(app, app->operations.recenter, app->operations.apply_pbc, app->operations.unwrap_structures)) {
                app->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
            break;
        }
        default:
            break;
        }
    }
}

md_script_ir_t* script_ir_create(ApplicationState* state) {
    ASSERT(state);
    md_script_ir_t* ir = md_script_ir_create(state->allocator.persistent);
    md_array_push(state->script.all_irs, ir, state->allocator.persistent);
    return ir;
}

void script_ir_collect(ApplicationState* state) {
    ASSERT(state);

    // A running task may be on an IR the fields have since moved on from
    if (task_system::task_is_running(state->tasks.evaluate) || task_system::task_is_running(state->script.vis_task)) {
        return;
    }

    for (size_t i = 0; i < md_array_size(state->script.all_irs);) {
        md_script_ir_t* ir = state->script.all_irs[i];
        if (ir != state->script.ir && ir != state->script.eval_ir) {
            md_script_ir_free(ir);
            md_array_swap_back_and_pop(state->script.all_irs, i);
        } else {
            ++i;
        }
    }
}

const md_script_ir_t* script_ir_of(const ApplicationState* state, md_script_vis_ref_t ref) {
    ASSERT(state);
    if (md_script_vis_ref_valid(state->script.ir, ref))      return state->script.ir;
    if (md_script_vis_ref_valid(state->script.eval_ir, ref)) return state->script.eval_ir;
    return nullptr;
}

void script_visualize_ref(ApplicationState* state, md_script_vis_ref_t ref, int subidx, md_script_vis_flags_t flags) {
    ASSERT(state);
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

void script_visualize_str(ApplicationState* state, str_t str, md_script_vis_flags_t flags) {
    ASSERT(state);
    if (str_empty(str)) {
        return;
    }
    ScriptVisTarget& t = state->script.vis_target;
    t = {};
    t.kind  = ScriptVisTarget::Str;
    t.str   = str_copy(str, state->allocator.frame);
    t.flags = flags;
}

// A visualization computed off the main thread. Everything it reads is its own or outlives it: a copy of the
// atoms' state, the IR (not freed while the task runs, see script_ir_collect), the system (which stays until the
// pool has finished, see interrupt_async_tasks). Everything it makes is in its arena.
struct ScriptVisJob {
    uint64_t key = 0;               // what it is of: the target, and the IR the target is resolved in
    uint64_t state_key = 0;         // where the atoms were
    md_allocator_i* arena = nullptr;
    ScriptVisTarget target = {};    // str in the arena
    const md_script_ir_t* ir = nullptr;
    const md_system_t* sys = nullptr;
    md_system_state_t state = {};
    md_script_vis_t vis = {};
    bool ok = false;
};

static void vis_job_free(ScriptVisJob* job) {
    if (job) {
        md_arena_allocator_destroy(job->arena);
        delete job;
    }
}

// The IR a target is resolved in: the one which made the reference, or the script being edited for an expression
static const md_script_ir_t* vis_target_ir(const ApplicationState* state, const ScriptVisTarget& t) {
    return t.kind == ScriptVisTarget::Ref ? script_ir_of(state, t.ref) : state->script.ir;
}

static uint64_t vis_target_key(const ApplicationState* state, const ScriptVisTarget& t) {
    const md_script_ir_t* ir = vis_target_ir(state, t);
    const uintptr_t ir_ptr = (uintptr_t)ir;
    const uint64_t fingerprint = md_script_ir_fingerprint(ir);
    uint64_t h = md_hash64(&t.kind, sizeof(t.kind), 0);
    h = md_hash64(&t.subidx, sizeof(t.subidx), h);
    h = md_hash64(&t.flags, sizeof(t.flags), h);
    h = md_hash64(&t.ref.ir_id, sizeof(t.ref.ir_id), h);
    h = md_hash64(&t.ref.node_idx, sizeof(t.ref.node_idx), h);
    if (t.kind == ScriptVisTarget::Str) {
        h = md_hash64(t.str.ptr, t.str.len, h);
    }
    h = md_hash64(&ir_ptr, sizeof(ir_ptr), h);
    h = md_hash64(&fingerprint, sizeof(fingerprint), h);
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
        h = md_hash64(st.xyz, st.num_atoms * sizeof(vec3_t), h);
    }
    return h;
}

static void vis_job_dispatch(ApplicationState* state, uint64_t key, uint64_t state_key) {
    const ScriptVisTarget& t = state->script.vis_target;
    const md_script_ir_t* ir = vis_target_ir(state, t);
    if (t.kind == ScriptVisTarget::Ref && !ir) {
        return;
    }

    ScriptVisJob* job = new ScriptVisJob();
    job->key = key;
    job->state_key = state_key;
    job->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    job->target = t;
    job->target.str = str_copy(t.str, job->arena);
    job->ir = ir;
    job->sys = &state->mold.sys;
    job->state.alloc = job->arena;
    md_system_state_copy(&job->state, &state->mold.state);
    md_script_vis_init(&job->vis, job->arena);

    const task_system::ID task = task_system::create_pool_task(STR_LIT("##Visualize script"), [job]() {
        md_script_vis_ctx_t ctx = {
            .ir    = job->ir,
            .sys   = job->sys,
            .state = &job->state,
        };
        const ScriptVisTarget& tgt = job->target;
        job->ok = (tgt.kind == ScriptVisTarget::Ref)
            ? md_script_vis_eval_ref(&job->vis, tgt.ref, tgt.subidx, &ctx, tgt.flags)
            : md_script_vis_eval_string(&job->vis, tgt.str, &ctx, tgt.flags);
    });
    task_system::enqueue_task(task);

    state->script.vis_running = job;
    state->script.vis_task = task;
}

void script_vis_begin_frame(ApplicationState* state) {
    ASSERT(state);
    state->script.vis_target = {};
    state->script.vis_shown = nullptr;
}

void script_vis_update(ApplicationState* state) {
    ASSERT(state);
    auto& sc = state->script;

    // A finished computation becomes the result
    if (sc.vis_running && !task_system::task_is_running(sc.vis_task)) {
        vis_job_free(sc.vis_result);
        sc.vis_result = sc.vis_running;
        sc.vis_running = nullptr;
    }

    sc.vis_shown = nullptr;
    if (sc.vis_target.kind == ScriptVisTarget::None) {
        return;
    }

    const uint64_t key = vis_target_key(state, sc.vis_target);
    const uint64_t state_key = vis_state_key(state->mold.state);

    // At most one at a time. While one runs, whatever is set when it has finished is what comes next: the
    // targets in between are never computed.
    const bool computed = sc.vis_result && sc.vis_result->key == key && sc.vis_result->state_key == state_key;
    if (!sc.vis_running && !computed) {
        vis_job_dispatch(state, key, state_key);
    }

    // Shown only if it is of what is set. For where the atoms were, while it is computed for where they are.
    if (sc.vis_result && sc.vis_result->ok && sc.vis_result->key == key) {
        sc.vis_shown = &sc.vis_result->vis;
        if (!md_bitfield_empty(&sc.vis_shown->atom_mask)) {
            md_bitfield_copy(&state->selection.highlight_mask, &sc.vis_shown->atom_mask);
        }
    }
}

void script_vis_reset(ApplicationState* state) {
    ASSERT(state);
    auto& sc = state->script;
    if (sc.vis_running) {
        task_system::task_wait_for(sc.vis_task);
        vis_job_free(sc.vis_running);
        sc.vis_running = nullptr;
    }
    vis_job_free(sc.vis_result);
    sc.vis_result = nullptr;
    sc.vis_shown = nullptr;
    sc.vis_target = {};
}

void script_set_hovered_property(ApplicationState* state, str_t label, int population_idx) {
    // The label, a script property's identifier, is copied. One too long to fit is taken as none: the views match
    // it against identifiers, and a cut one would match the wrong property or nothing.
    char* dst = state->hovered_property_label;
    const size_t cap = sizeof(state->hovered_property_label);
    if (label.ptr && label.len < cap) {
        MEMCPY(dst, label.ptr, label.len);
        dst[label.len] = '\0';
    } else {
        dst[0] = '\0';
    }
    state->hovered_property_pop_idx = population_idx;
}
