#include <viamd.h>

#include <core/md_str_builder.h>
#include <md_util.h>

#include <imgui.h>

void draw_info_window(const ApplicationState& state, uint32_t picking_idx) {
    const auto& sys = state.mold.sys;

    if (picking_idx == INVALID_PICKING_IDX) return;

    md_strb_t sb = md_strb_create(state.allocator.frame);
    defer {
        md_strb_free(&sb);
    };

    if (picking_idx < sys.atom.count) {
        int atom_idx = picking_idx;
        int local_idx = atom_idx;
        const vec3_t pos = md_atom_coord(&sys.atom, atom_idx);
        str_t type = md_atom_name(&sys.atom, atom_idx);
        md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom_idx);
        str_t elem = z ? md_util_element_name(z)   : str_t{};
        str_t symb = z ? md_util_element_symbol(z) : str_t{};

        int comp_idx = md_comp_find_by_atom_idx(&sys.comp, atom_idx);
        str_t comp_name = {};
        int comp_seq_id = 0;
        if (comp_idx != -1) {
            comp_name   = md_comp_name(&sys.comp, comp_idx);
            comp_seq_id = md_comp_seq_id(&sys.comp, comp_idx);
            md_urange_t range = md_comp_atom_range(&sys.comp, comp_idx);
            local_idx = atom_idx - range.beg;
        }

        int inst_idx = md_system_inst_find_by_atom_idx(&sys, atom_idx);
        str_t chain_id = {};
        if (inst_idx != -1) {
            chain_id = md_inst_id(&sys.inst, inst_idx);
        }

        // External indices begin with 1 not 0
        md_strb_fmt(&sb, "atom[%i][%i]: %.*s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx + 1, local_idx + 1, STR_ARG(type), STR_ARG(elem), STR_ARG(symb), pos.x, pos.y, pos.z);
        if (comp_idx != -1) {
            md_strb_fmt(&sb, "res[%i]: %.*s %i\n", comp_idx + 1, STR_ARG(comp_name), comp_seq_id);
        }
        if (inst_idx != -1) {
            md_strb_fmt(&sb, "chain[%i]: %.*s\n", inst_idx + 1, STR_ARG(chain_id));
        }

        uint32_t flags = sys.atom.flags[atom_idx];
        if (flags) {
            sb += "flags: ";
            if (flags & MD_FLAG_HETERO)             { sb += "HETERO "; }
            if (flags & MD_FLAG_AMINO_ACID)         { sb += "AMINO "; }
            //if (flags & MD_FLAG_SIDE_CHAIN)         { sb += "SIDE_CHAIN "; }
            if (flags & MD_FLAG_NUCLEOTIDE)         { sb += "NUCLEOTIDE "; }
            if (flags & MD_FLAG_NUCLEOBASE)         { sb += "NUCLEOBASE "; }
            //if (flags & MD_FLAG_NUCLEOSIDE)         { sb += "NUCLEOSIDE "; }
            if (flags & MD_FLAG_WATER)              { sb += "WATER "; }
            if (flags & MD_FLAG_ION)                { sb += "ION "; }
            if (flags & MD_FLAG_BACKBONE)           { sb += "BACKBONE "; }
            if (flags & MD_FLAG_SP)                 { sb += "SP "; }
            if (flags & MD_FLAG_SP2)                { sb += "SP2 "; }
            if (flags & MD_FLAG_SP3)                { sb += "SP3 "; }
            if (flags & MD_FLAG_AROMATIC)           { sb += "AROMATIC "; }
            sb += "\n";
        }
        /*
        // @TODO: REIMPLEMENT THIS
        if (res_idx < sys.backbone.segment.angleangles.size() && res_idx < sys.backbone.segments.size() && valid_backbone_atoms(sys.backbone.segments[res_idx])) {
        const auto angles = RAD_TO_DEG((vec2)sys.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
        }
        */
    }
    else if (picking_idx >= 0x80000000) {
        int bond_idx = picking_idx & 0x7FFFFFFF;
        if (0 <= bond_idx && bond_idx < (int)sys.bond.count) {
            md_atom_pair_t pair = sys.bond.pairs[bond_idx];
            md_flags_t    flags = sys.bond.flags[bond_idx];
            char bond_flags_buf[256] = {};
            int  len = 0;

            typedef struct {
                md_flags_t flag;
                const char* label;
            } bond_flag_label_t;

            bond_flag_label_t bond_flag_map[] = {
                {MD_BOND_FLAG_COVALENT,     "COVALENT"},
                {MD_BOND_FLAG_DOUBLE,       "DOUBLE"},
                {MD_BOND_FLAG_TRIPLE,       "TRIPLE"},
                {MD_BOND_FLAG_QUADRUPLE,    "QUADRUPLE"},
                {MD_BOND_FLAG_AROMATIC,     "AROMATIC"},
                {MD_BOND_FLAG_INTER,        "INTER"},
                {MD_BOND_FLAG_COORDINATE,   "COORD"},
				{MD_BOND_FLAG_METAL,        "METAL"},
				{MD_BOND_FLAG_USER_DEFINED, "USER"},
            };

            for (size_t i = 0; i < ARRAY_SIZE(bond_flag_map); ++i) {
                if (flags & bond_flag_map[i].flag) {
                    len += snprintf(bond_flags_buf + len, 256 - len, "%s ", bond_flag_map[i].label);
                }
            }
            
            char bond_type = '-';
            if (flags & MD_BOND_FLAG_DOUBLE) bond_type = '=';
            if (flags & MD_BOND_FLAG_TRIPLE) bond_type = '#';
            if (flags & MD_BOND_FLAG_QUADRUPLE) bond_type = '$';
            if (flags & MD_BOND_FLAG_AROMATIC) bond_type = ':';

            vec3_t p0 = md_atom_coord(&sys.atom, pair.idx[0]);
            vec3_t p1 = md_atom_coord(&sys.atom, pair.idx[1]);
            float d = vec3_distance(p0, p1);

            str_t type0 = md_atom_name(&sys.atom, pair.idx[0]);
            str_t type1 = md_atom_name(&sys.atom, pair.idx[1]);

            md_strb_fmt(&sb, "bond: " STR_FMT "%c" STR_FMT "\n", STR_ARG(type0), bond_type, STR_ARG(type1));
            md_strb_fmt(&sb, "flags: %.*s\n", len, bond_flags_buf);
            md_strb_fmt(&sb, "length: %.3f\n", d);
        }
    }

    const ImVec2 offset = { 10.f, 18.f };
    const ImVec2 new_pos = {ImGui::GetMousePos().x + offset.x, ImGui::GetMousePos().y + offset.y};
    ImGui::SetNextWindowPos(new_pos);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0,
        ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
    ImGui::Text("%s", md_strb_to_cstr(sb));
    ImGui::End();
    ImGui::PopStyleColor();
}

void extract_picking_data(PickingData& out_picking, GBuffer& gbuffer, const vec2_t& coord, const mat4_t& inv_MVP) {
    out_picking = {};

#if MD_PLATFORM_OSX
    coord = coord * vec_cast(ImGui::GetIO().DisplayFramebufferScale);
#endif
    if (0.f < coord.x && coord.x < (float)gbuffer.width && 0.f < coord.y && coord.y < (float)gbuffer.height) {
        extract_gbuffer_picking_idx_and_depth(&out_picking.idx, &out_picking.depth, &gbuffer, (int)coord.x, (int)coord.y);
        const vec4_t viewport = {0, 0, (float)gbuffer.width, (float)gbuffer.height};
        out_picking.world_coord = mat4_unproject({coord.x, coord.y, out_picking.depth}, inv_MVP, viewport);
        out_picking.screen_coord = {coord.x, coord.y};
    }
}

md_atom_idx_t atom_idx_from_picking_idx(uint32_t picking_idx) {
    if (picking_idx < 0x80000000) {
        return (md_atom_idx_t)picking_idx;
    }
    return -1;
}

md_bond_idx_t bond_idx_from_picking_idx(uint32_t picking_idx) {
    if (picking_idx >= 0x80000000) {
        return (md_bond_idx_t)(picking_idx & 0x7FFFFFFF);
    }
    return -1;
}
