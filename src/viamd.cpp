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
        const vec3_t pos = { mol.atom.x[atom_idx], mol.atom.y[atom_idx], mol.atom.z[atom_idx] };
        str_t type = mol.atom.type ? mol.atom.type[atom_idx] : str_t{};
        str_t elem = mol.atom.element ? md_util_element_name(mol.atom.element[atom_idx]) : str_t{};
        str_t symbol = mol.atom.element ? md_util_element_symbol(mol.atom.element[atom_idx]) : str_t{};

        int res_idx = -1;
        str_t res_name = {};
        int res_id = 0;
        if (mol.residue.count && mol.atom.res_idx) {
            res_idx = mol.atom.res_idx[atom_idx];
            res_name = mol.residue.name[res_idx];
            res_id = mol.residue.id[res_idx];
            md_range_t range = md_residue_atom_range(mol.residue, res_idx);
            local_idx = atom_idx - range.beg;
        }

        int chain_idx = -1;
        str_t chain_id = {};
        if (mol.chain.count && mol.atom.chain_idx) {
            chain_idx = mol.atom.chain_idx[atom_idx];
            if (0 <= chain_idx && chain_idx < (int)mol.chain.count) {
                chain_id = mol.chain.id[chain_idx];
            }
        }

        // External indices begin with 1 not 0
        md_strb_fmt(&sb, "atom[%i][%i]: %.*s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx + 1, local_idx + 1, STR_ARG(type), STR_ARG(elem), STR_ARG(symbol), pos.x, pos.y, pos.z);
        if (res_idx != -1) {
            md_strb_fmt(&sb, "res[%i]: %.*s %i\n", res_idx + 1, STR_ARG(res_name), res_id);
        }
        if (chain_idx != -1) {
            md_strb_fmt(&sb, "chain[%i]: %.*s\n", chain_idx + 1, STR_ARG(chain_id));
        }

        uint32_t flags = mol.atom.flags[atom_idx];
        if (flags) {
            sb += "flags: ";
            if (flags & MD_FLAG_RES_BEG)            { sb += "RES_BEG "; }
            if (flags & MD_FLAG_RES)                { sb += "RES "; }
            if (flags & MD_FLAG_RES_END)            { sb += "RES_END "; }
            if (flags & MD_FLAG_CHAIN_BEG)          { sb += "CHAIN_BEG "; }
            if (flags & MD_FLAG_CHAIN)              { sb += "CHAIN "; }
            if (flags & MD_FLAG_CHAIN_END)          { sb += "CHAIN_END "; }
            if (flags & MD_FLAG_HETATM)             { sb += "HETATM "; }
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
            char bond_type;
            switch (mol.bond.order[bond_idx]) {
            case 1: bond_type = '-'; break;
            case 2: bond_type = '='; break;
            case 3: bond_type = '#'; break;
            case 4: bond_type = '$'; break;
            default: bond_type = '?'; break;
            }
            const float dx = mol.atom.x[b.idx[0]] - mol.atom.x[b.idx[1]];
            const float dy = mol.atom.y[b.idx[0]] - mol.atom.y[b.idx[1]];
            const float dz = mol.atom.z[b.idx[0]] - mol.atom.z[b.idx[1]];
            const float d = sqrtf(dx*dx + dy*dy + dz*dz);
            md_strb_fmt(&sb, "bond: %s%c%s\n", mol.atom.type[b.idx[0]].buf, bond_type, mol.atom.type[b.idx[1]].buf);
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