#include <viamd.h>

#include <core/md_str_builder.h>
#include <md_util.h>

#include <imgui.h>

void draw_info_window(const ApplicationState& state, uint32_t picking_idx) {
    const auto& mol = state.mold.mol;

    if (picking_idx == INVALID_PICKING_IDX) return;

    md_strb_t sb = md_strb_create(state.allocator.frame);
    defer {
        md_strb_free(&sb);
    };

    if (picking_idx < mol.atom.count) {
        int atom_idx = picking_idx;
        int local_idx = atom_idx;
        const vec3_t pos = md_atom_coord(&mol.atom, atom_idx);
        str_t type = md_atom_name(&mol.atom, atom_idx);
        md_atomic_number_t z = md_atom_atomic_number(&mol.atom, atom_idx);
        str_t elem = z ? md_util_element_name(z)   : str_t{};
        str_t symb = z ? md_util_element_symbol(z) : str_t{};

        int res_idx = md_residue_find_by_atom_idx(&mol.residue, atom_idx);
        str_t res_name = {};
        int res_id = 0;
        if (res_idx != -1) {
            res_name = md_residue_name(&mol.residue, res_idx);
            res_id   = md_residue_id(&mol.residue, res_idx);
            md_range_t range = md_residue_atom_range(&mol.residue, res_idx);
            local_idx = atom_idx - range.beg;
        }

        int chain_idx = md_chain_find_by_atom_idx(&mol.chain, atom_idx);
        str_t chain_id = {};
        if (chain_idx != -1) {
            chain_id = md_chain_id(&mol.chain, chain_idx);
        }

        // External indices begin with 1 not 0
        md_strb_fmt(&sb, "atom[%i][%i]: %.*s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx + 1, local_idx + 1, STR_ARG(type), STR_ARG(elem), STR_ARG(symb), pos.x, pos.y, pos.z);
        if (res_idx != -1) {
            md_strb_fmt(&sb, "res[%i]: %.*s %i\n", res_idx + 1, STR_ARG(res_name), res_id);
        }
        if (chain_idx != -1) {
            md_strb_fmt(&sb, "chain[%i]: %.*s\n", chain_idx + 1, STR_ARG(chain_id));
        }

        uint32_t flags = mol.atom.flags[atom_idx];
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
            //if (flags & MD_FLAG_BACKBONE)           { sb += "BACKBONE "; }
            if (flags & MD_FLAG_SP)                 { sb += "SP "; }
            if (flags & MD_FLAG_SP2)                { sb += "SP2 "; }
            if (flags & MD_FLAG_SP3)                { sb += "SP3 "; }
            if (flags & MD_FLAG_AROMATIC)           { sb += "AROMATIC "; }
            sb += "\n";
        }
        /*
        // @TODO: REIMPLEMENT THIS
        if (res_idx < mol.backbone.segment.angleangles.size() && res_idx < mol.backbone.segments.size() && valid_backbone_atoms(mol.backbone.segments[res_idx])) {
        const auto angles = RAD_TO_DEG((vec2)mol.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
        }
        */
    }
    else if (picking_idx >= 0x80000000) {
        int bond_idx = picking_idx & 0x7FFFFFFF;
        if (0 <= bond_idx && bond_idx < (int)mol.bond.count) {
            md_bond_pair_t b = mol.bond.pairs[bond_idx];
            char bond_type = ' ';
            switch (mol.bond.order[bond_idx] & MD_BOND_ORDER_MASK) {
                case 1: bond_type = '-'; break;
                case 2: bond_type = '='; break;
                case 3: bond_type = '#'; break;
                case 4: bond_type = '$'; break;
                default: bond_type = '?'; break;
            }
            vec3_t p0 = md_atom_coord(&mol.atom, b.idx[0]);
            vec3_t p1 = md_atom_coord(&mol.atom, b.idx[1]);
            float d = vec3_distance(p0, p1);

            str_t type0 = md_atom_name(&mol.atom, b.idx[0]);
            str_t type1 = md_atom_name(&mol.atom, b.idx[1]);

            md_strb_fmt(&sb, "bond: " STR_FMT "%c" STR_FMT "\n", STR_ARG(type0), bond_type, STR_ARG(type1));
            md_strb_fmt(&sb, "order: %i\n", mol.bond.order[bond_idx]);
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