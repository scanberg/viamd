#include <event.h>
#include <viamd.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <md_system.h>
#include <md_util.h>

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#include <string>

namespace dataset {

// Helper function to convert amino acids and nucleotides from three-letter to single-letter codes
static str_t convert_to_short(str_t in_str) {
    
    // Standard amino acids
    if (str_eq_cstr(in_str, "ALA")) return STR_LIT("A");
    if (str_eq_cstr(in_str, "ARG")) return STR_LIT("R");
    if (str_eq_cstr(in_str, "ASN")) return STR_LIT("N");
    if (str_eq_cstr(in_str, "ASP")) return STR_LIT("D");
    if (str_eq_cstr(in_str, "CYS")) return STR_LIT("C");
    if (str_eq_cstr(in_str, "GLU")) return STR_LIT("E");
    if (str_eq_cstr(in_str, "GLN")) return STR_LIT("Q");
    if (str_eq_cstr(in_str, "GLY")) return STR_LIT("G");
    if (str_eq_cstr(in_str, "HIS")) return STR_LIT("H");
    if (str_eq_cstr(in_str, "ILE")) return STR_LIT("I");
    if (str_eq_cstr(in_str, "LEU")) return STR_LIT("L");
    if (str_eq_cstr(in_str, "LYS")) return STR_LIT("K");
    if (str_eq_cstr(in_str, "MET")) return STR_LIT("M");
    if (str_eq_cstr(in_str, "PHE")) return STR_LIT("F");
    if (str_eq_cstr(in_str, "PRO")) return STR_LIT("P");
    if (str_eq_cstr(in_str, "SER")) return STR_LIT("S");
    if (str_eq_cstr(in_str, "THR")) return STR_LIT("T");
    if (str_eq_cstr(in_str, "TRP")) return STR_LIT("W");
    if (str_eq_cstr(in_str, "TYR")) return STR_LIT("Y");
    if (str_eq_cstr(in_str, "VAL")) return STR_LIT("V");

    // Standard nucleotides
    if (str_eq_cstr(in_str, "DA")) return STR_LIT("A");
    if (str_eq_cstr(in_str, "DC")) return STR_LIT("C");
    if (str_eq_cstr(in_str, "DG")) return STR_LIT("G");
    if (str_eq_cstr(in_str, "DT")) return STR_LIT("T");
    if (str_eq_cstr(in_str, "DU")) return STR_LIT("U");
    
    // Failed to map it, return itself
    return in_str;
}

// Essentially a bswap to go from RGBA to ABGR
#define RGBA_HEX(x) BSWAP32(x)

// Helper function to assign colors to components based on their type
static uint32_t component_color(str_t str) {
    if (str_eq_cstr(str, "DG")) return RGBA_HEX(0xD5B3EFAA); // Purple
    if (str_eq_cstr(str, "DT")) return RGBA_HEX(0x9FF3A0AA); // Green
    if (str_eq_cstr(str, "DC")) return RGBA_HEX(0xF8EE5CAA); // Yellow
    if (str_eq_cstr(str, "DU")) return RGBA_HEX(0xFE9D2DAA); // Orange
    if (str_eq_cstr(str, "DA")) return RGBA_HEX(0xFC697AAA); // Red

    return 0;
}

struct AtomElementMapping {
    char lbl[31] = "";
    md_element_t elem = 0;
};

// We use this to represent a single entity within the loaded system, e.g. a residue type
// This is used to represent multiple types, so all fields are not used in all cases
struct DatasetItem {
    char label[32] = "";
    uint32_t count = 0;
    float fraction = 0;
    
    // Extended metadata for popups
    uint64_t key = 0;            // Unique key of the type
    md_array(int) indices = 0;   // Indices into the corresponding structures which are represented by this item: i.e. chain or residue indices (for highlighting)
    md_array(int) sub_items = 0; // Indices into the items of the subcatagories: i.e. for chain -> unique residues types within that chain

    // Atom type only
    bool use_defaults = true; // Flag if this particular atom type should be linked to the default values stemming for the element (only applicable to atom types with an element, i.e. not coarse grained)
};

struct ElementDefault {
    vec4_t color;
    float radius;
    float mass;
};

struct Dataset : viamd::EventHandler {
    bool show_window = false;
    
    // Dataset data (moved from ApplicationState.dataset)
    md_array(AtomElementMapping) atom_element_remappings = 0;
    md_array(DatasetItem) inst_types = 0;
    md_array(DatasetItem) comp_types = 0;
    md_array(DatasetItem) atom_types = 0;
    md_allocator_i* arena = 0;

    ElementDefault element_defaults[MD_Z_Count] = {};

    Dataset() { 
        viamd::event_system_register_handler(*this); 
    }
    
    ~Dataset() {
        // The arena will be cleaned up when the persistent allocator is destroyed
        // No need for explicit cleanup here
    }

    void init_element_defaults() {
        for (md_atomic_number_t z = 0; z < MD_Z_Count; ++z) {
            element_defaults[z].color   = vec4_from_u32(md_atomic_number_cpk_color(z));
            element_defaults[z].mass    = md_atomic_number_mass(z);
            element_defaults[z].radius  = md_atomic_number_vdw_radius(z);
        }
    }

    void clear_dataset_items() {
        inst_types = 0;
        comp_types = 0;
        atom_types = 0;
        if (arena) {
            md_arena_allocator_reset(arena);
        }
    }

    void init_dataset_items(ApplicationState& data) {
        if (!arena) {
            arena = md_arena_allocator_create(data.allocator.persistent, MEGABYTES(1));
        }
        clear_dataset_items();

        md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));

        const md_system_t& sys = data.mold.sys;

        size_t atom_type_count = md_atom_type_count(&sys.atom.type);
        size_t atom_count = md_atom_count(&sys.atom);
        size_t comp_count = md_comp_count(&sys.comp);
        size_t inst_count = md_inst_count(&sys.inst);
		size_t entity_count = md_entity_count(&sys.entity);

        if (atom_count == 0) return;

        // Map atom types into dataset items
        for (size_t i = 0; i < atom_type_count; ++i) {
            str_t atom_type_name = md_atom_type_name(&sys.atom.type, i);
            DatasetItem item = { .key = i };
            snprintf(item.label, sizeof(item.label), STR_FMT, STR_ARG(atom_type_name));
            md_array_push(atom_types, item, arena);
        }

        // Count and set indices for each atom type
        for (size_t i = 0; i < atom_count; ++i) {
            md_atom_type_idx_t type_idx = sys.atom.type_idx[i]; 
            atom_types[type_idx].count += 1;
            md_array_push(atom_types[type_idx].indices, (int)i, arena);
        }
        
        // Calculate fractions
        for (size_t i = 0; i < atom_type_count; ++i) {
            atom_types[i].fraction = atom_types[i].count / (float)atom_count;
        }

        md_array(int) sequence = 0;
        md_array(int) comp_idx_type = 0;

        if (comp_count > 0) {
			md_array_resize(comp_idx_type, comp_count, temp_arena);
			MEMSET(comp_idx_type, -1, md_array_bytes(comp_idx_type));
        }

        // Process components - group by name + atom type sequence
        for (size_t i = 0; i < comp_count; ++i) {
            str_t comp_name = md_comp_name(&sys.comp, i);

            md_array_shrink(sequence, 0);
            md_urange_t range = md_comp_atom_range(&sys.comp, i);
            for (int j = range.beg; j < range.end; ++j) {
                int ai = sys.atom.type_idx[j];
                md_array_push(sequence, ai, temp_arena);
            }

            // Create combined string for hash (label + sequence of atom types)
            uint64_t hash = md_hash64_str(comp_name, md_hash64(sequence, md_array_bytes(sequence), 0));

            // Check if we already have this chain type
            DatasetItem* item = nullptr;
            for (size_t j = 0; j < md_array_size(comp_types); ++j) {
                if (comp_types[j].key == hash) {
                    item = &comp_types[j];
                    md_array_push(comp_idx_type, (int)j, temp_arena);
                    break;
                }
            }
            if (!item) {
                DatasetItem it = { .key = hash };
                snprintf(it.label, sizeof(it.label), STR_FMT, STR_ARG(comp_name));
                md_array_push(comp_idx_type, (int)md_array_size(comp_types), temp_arena);
                md_array_push(comp_types, it, arena);
                item = md_array_last(comp_types);
                md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
            }

            size_t comp_atom_count = md_comp_atom_count(&sys.comp, i);

            item->count += 1;
            item->fraction += (float)(comp_atom_count / (double)atom_count);
            md_array_push(item->indices, (int)i, arena);
        }

        // Process chains - key = residue type sequence
        if (comp_count > 0 && inst_count > 0) {
            for (size_t i = 0; i < inst_count; ++i) {
                str_t inst_id = md_inst_id(&sys.inst, i);

                md_array_shrink(sequence, 0);
                md_urange_t range = md_inst_comp_range(&sys.inst, i);
                for (int j = range.beg; j < range.end; ++j) {
                    int res_type_idx = comp_idx_type[j];
                    md_array_push(sequence, res_type_idx, temp_arena);
                }

                // Create combined string for hash (label + sequence of residue types)
                uint64_t hash = md_hash64(sequence, md_array_bytes(sequence), 0);

                // Check if we already have this chain type
                DatasetItem* item = nullptr;
                for (size_t j = 0; j < md_array_size(inst_types); ++j) {
                    if (inst_types[j].key == hash) {
                        item = &inst_types[j];
                        break;
                    }
                }
                if (!item) {
					int c_idx = (int)md_array_size(inst_types);
                    DatasetItem it = { .key = hash };
                    snprintf(it.label, sizeof(it.label), "Type %i", c_idx + 1);
                    md_array_push(inst_types, it, arena);
                    item = md_array_last(inst_types);
                    md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
                }

                size_t chain_atom_count = md_system_inst_atom_count(&sys, i);

                item->count += 1;
                item->fraction += (float)(chain_atom_count / (double)atom_count);
                md_array_push(item->indices, (int)i, arena);
            }
        }

		md_vm_arena_destroy(temp_arena);
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                // Initialize component
                init_element_defaults();
                break;
            }
            case viamd::EventType_ViamdShutdown:
                // Cleanup
                clear_dataset_items();
                break;
            case viamd::EventType_ViamdTopologyInit: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                init_dataset_items(state);
                break;
            }
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Dataset", &show_window);
                break;
            default:
                break;
            }
        }
    }

    static inline void handle_item_click(ApplicationState& state) {
        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
            // Shift + click on element header to select all atom types that use this element
            // Since the hover should already contain the correct selection we can just copy it to the filter mask to achieve this
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
            } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
            }
        }
    }

    // Return the row in the periodic table for a given atomic number
    static inline int element_row(int z) {
        if (z == 0) return 9; // Placeholder for unknown elements

        if (z <= 2) return 0;
        if (z <= 10) return 1;
        if (z <= 18) return 2;
        if (z <= 36) return 3;
        if (z <= 54) return 4;

        if (57 <= z && z <= 71) return 8; // Lanthanides
        if (89 <= z && z <= 103) return 9; // Actinides

        if (z <= 86) return 5;

        return 6; // For actinides and beyond
    }

    // Return the column in the periodic table for a given atomic number
    static inline int element_col(int z) {
        if (z == 0) return 0; // Placeholder for unknown elements
        if (z == 1) return 0;
        if (z == 2) return 17;
        if (z <= 4)  return z - 3;
        if (z <= 10) return z + 7;

        if (z <= 12) return z - 11;
        if (z <= 18) return z - 1;
        if (z <= 36) return (z - 1) % 18;
        if (z <= 54) return (z - 1) % 18;
        if (z <= 56) return (z - 1) % 18;
        if (z <= 71) return (z - 57) % 18 + 2; // Lanthanides
        if (z <= 86) return z - 69;
        if (z <= 88) return z - 87;
        if (z <= 103) return (z - 89) % 18 + 2; // Actinides
        if (z <= 118) return (z - 101) % 18;

        return 0;
    }

    struct PeriodicTableResult {
        bool hovered = false;
        bool clicked = false;
        int  z = -1;
    };

    static bool element_button(const char* lbl, vec4_t color) {
        const float button_size = ImGui::GetFontSize() * 1.5f;
        const vec4_t color_hover  = vec4_clamp(color + vec4_set1(0.2f), vec4_set1(0.0f), vec4_set1(1.0f));
        const vec4_t color_active = color * vec4_set(0.8f, 0.8f, 0.8f, 1.0f);

        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(color.x, color.y, color.z, color.w));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0, 0, 0, 1));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(color_hover.x,  color_hover.y,  color_hover.z,  1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(color_active.x, color_active.y, color_active.z, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0, 0, 0, 1));

        bool clicked = ImGui::Button(lbl, ImVec2(button_size, button_size));

        ImGui::PopStyleColor(5);
		ImGui::PopStyleVar();

		return clicked;
    }

    // If `enabled_mask` is non-null: bit=1 means enabled.
    // Returns per-frame interaction info.
    static PeriodicTableResult periodic_table_widget(const ElementDefault* elem_defs, const uint64_t enabled_mask[2] = nullptr) {
        PeriodicTableResult res = {};

        const int width  = 18;
        const int height = 10;
        const float button_size = ImGui::GetFontSize() * 1.5f;
        const float spacing = ImGui::GetStyle().ItemSpacing.x * 0.25f;

        ImVec2 offset = ImGui::GetCursorPos();
        ImGui::Dummy(ImVec2(width * (button_size + spacing), height * (button_size + spacing)));

        for (int z = 0; z < MD_Z_Count; ++z) {
            int row = element_row(z);
            int col = element_col(z);

            ImGui::PushID(z);
            ImGui::SetCursorPos(offset + ImVec2(col * (button_size + spacing), row * (button_size + spacing)));

            bool disabled = enabled_mask ? !(enabled_mask[z / 64] & (1ULL << (z % 64))) : false;

            const ElementDefault& def = elem_defs[z];
            const vec4_t color = def.color;
            const str_t sym = md_atomic_number_symbol((md_atomic_number_t)z);

			ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.25f);

            if (disabled) ImGui::BeginDisabled();
			if (element_button(str_ptr(sym), color)) {
                res.clicked = true;
                res.z = z;
            }
            if (disabled) ImGui::EndDisabled();

            ImGui::PopStyleVar();

            // Interaction is queried after the item
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                res.hovered = true;
                res.z = z;
            }

            str_t name = md_atomic_number_name((md_atomic_number_t)z);
            ImGui::SetItemTooltip("%d: " STR_FMT, z, STR_ARG(name));

            ImGui::PopID();
        }

        return res;
    }

    void draw(ApplicationState& data) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(500, 600), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Dataset", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::Text("System: %s", data.files.molecule);
            ImGui::Text("Num entities:   %9zu", data.mold.sys.entity.count);
            ImGui::Text("Num instances:  %9zu", data.mold.sys.inst.count);
            ImGui::Text("Num components: %9zu", data.mold.sys.comp.count);
            ImGui::Text("Num atoms:      %9zu", data.mold.sys.atom.count);

            ImGui::Separator();

            if (data.mold.sys.unitcell.flags) {
                md_flags_t flags = data.mold.sys.unitcell.flags;
                bool ortho = flags & MD_UNITCELL_ORTHO;
                bool tricl = flags & MD_UNITCELL_TRICLINIC;
                ImGui::Text("Unitcell %s", ortho ? "Orthorhombic" : tricl ? "Triclinic" : "");
                if (flags & MD_UNITCELL_ORTHO) {
                    ImGui::Indent();
                    ImGui::Text("X: %f %s", data.mold.sys.unitcell.x, flags & MD_UNITCELL_PBC_X ? "(pbc)" : "");
                    ImGui::Text("Y: %f %s", data.mold.sys.unitcell.y, flags & MD_UNITCELL_PBC_Y ? "(pbc)" : "");
                    ImGui::Text("Z: %f %s", data.mold.sys.unitcell.z, flags & MD_UNITCELL_PBC_Z ? "(pbc)" : "");
                    ImGui::Unindent();
                } else if (flags & MD_UNITCELL_TRICLINIC) {
                    ImGui::Indent();
                    ImGui::Text("X:  %f", data.mold.sys.unitcell.x);
                    ImGui::Text("XY: %f", data.mold.sys.unitcell.xy);
                    ImGui::Text("XZ: %f", data.mold.sys.unitcell.xz);
                    ImGui::Text("Y:  %f", data.mold.sys.unitcell.y);
                    ImGui::Text("YZ: %f", data.mold.sys.unitcell.yz);
                    ImGui::Text("Z:  %f", data.mold.sys.unitcell.z);
                    ImGui::Unindent();
                } 
                ImGui::Separator();
            }

            if (data.files.trajectory[0] != '\0') {
                ImGui::Text("Trajectory data: %s", data.files.trajectory);
                ImGui::Text("Num frames:    %9zu", md_trajectory_num_frames(data.mold.traj));
                ImGui::Text("Num atoms:     %9zu", md_trajectory_num_atoms(data.mold.traj));
                ImGui::Separator();
            }

            if (ImGui::IsWindowHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
            }

            static bool use_short_labels = true;
            static bool use_auth_labels  = true;

            ImGui::Checkbox("Use shorthand component labels", &use_short_labels);
            ImGui::Checkbox("Use author instance labels", &use_auth_labels);

            size_t num_entities  = md_system_entity_count(&data.mold.sys);
            size_t num_instances = md_system_inst_count(&data.mold.sys);
			size_t num_atom_types = md_system_atom_type_count(&data.mold.sys);

            if (num_entities && ImGui::CollapsingHeader("Entities", ImGuiTreeNodeFlags_DefaultOpen)) {
                const ImVec2 item_size = ImVec2(ImGui::GetFontSize() * 1.4f, ImGui::GetFontSize() * 1.1f);

                for (size_t ent_idx = 0; ent_idx < num_entities; ++ent_idx) {
                    char buf[256];
                    str_t entity_id   = md_entity_id(&data.mold.sys.entity, ent_idx);
                    str_t entity_desc = md_entity_description(&data.mold.sys.entity, ent_idx);
                    md_flags_t entity_flags = md_entity_flags(&data.mold.sys.entity, ent_idx);

                    snprintf(buf, sizeof(buf), STR_FMT ": " STR_FMT, STR_ARG(entity_id), STR_ARG(entity_desc));
                    ImGui::Indent();

                    bool expand_entity = ImGui::CollapsingHeader(buf);
                    if (ImGui::IsItemHovered()) {
                        md_bitfield_clear(&data.selection.highlight_mask);
                        for (size_t inst_idx = 0; inst_idx < num_instances; ++inst_idx) {
                            if (md_inst_entity_idx(&data.mold.sys.inst, inst_idx) == (int)ent_idx) {
                                md_urange_t range = md_system_inst_atom_range(&data.mold.sys, inst_idx);
                                md_bitfield_set_range(&data.selection.highlight_mask, range.beg, range.end);
                            }
                        }
                        handle_item_click(data);
                    }
                    if (expand_entity) {
                        ImGui::Indent();
                        for (size_t inst_idx = 0; inst_idx < num_instances; ++inst_idx) {
                            if (md_inst_entity_idx(&data.mold.sys.inst, inst_idx) != (int)ent_idx) continue;
                            ImGui::PushID((int)inst_idx);  // avoid ID collisions if names repeat
                            defer { ImGui::PopID(); };

                            str_t inst_id = STR_LIT("");
                            if (use_auth_labels) {
                                inst_id = md_inst_auth_id(&data.mold.sys.inst, inst_idx);
                            } else {
                                inst_id = md_inst_id(&data.mold.sys.inst, inst_idx);
                            }

                            snprintf(buf, sizeof(buf), STR_FMT, STR_ARG(inst_id));
                            bool expand_inst = ImGui::CollapsingHeader(buf);
                            if (ImGui::IsItemHovered()) {
                                md_bitfield_clear(&data.selection.highlight_mask);
                                md_urange_t range = md_system_inst_atom_range(&data.mold.sys, inst_idx);
                                md_bitfield_set_range(&data.selection.highlight_mask, range.beg, range.end);
                                handle_item_click(data);
                            }
                            if (expand_inst) {
                                const ImGuiStyle& style = ImGui::GetStyle();
                                md_urange_t range = md_system_inst_comp_range(&data.mold.sys, inst_idx);

                                for (size_t comp_idx = range.beg; comp_idx < range.end; ++comp_idx) {
                                    str_t comp_name = md_comp_name(&data.mold.sys.comp, comp_idx);
                                    if (use_short_labels) {
                                        uint32_t color = component_color(comp_name);
                                        comp_name = convert_to_short(comp_name);
                                        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
                                        ImGui::PushStyleColor(ImGuiCol_Header, color);
                                        ImGui::PushStyleColor(ImGuiCol_HeaderActive, color);
                                    }

                                    // Measure label + padding (match what we pass to Selectable)
                                    ImVec2 text_sz = ImGui::CalcTextSize(str_beg(comp_name), str_end(comp_name));
                                    ImVec2 item_sz(text_sz.x + style.ItemSpacing.x * 2.0f,
                                        text_sz.y + style.ItemSpacing.x * 2.0f);

                                    // Flow: stay on same line only if it fits to the right of the previous item
                                    if (comp_idx != range.beg) {
                                        float last_x = ImGui::GetItemRectMax().x; // right edge of previous item (screen coords)
                                        float max_x  = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x; // window content right edge (screen coords)
                                        if (last_x + style.ItemSpacing.x + item_sz.x <= max_x) {
                                            ImGui::SameLine();
                                        }
                                    }

                                    ImGui::PushID((int)comp_idx); // avoid ID collisions if names repeat
                                    ImGui::Selectable(comp_name.ptr, true, 0, text_sz);
                                    ImGui::PopID();

                                    if (ImGui::IsItemHovered()) {
                                        ImGui::SetTooltip("%d", md_comp_seq_id(&data.mold.sys.comp, comp_idx));
                                    }

                                    if (use_short_labels) {
                                        ImGui::PopStyleVar();
                                        ImGui::PopStyleColor(2);
                                    }

                                    if (ImGui::IsItemHovered()) {
                                        md_bitfield_clear(&data.selection.highlight_mask);
                                        md_urange_t atom_range = md_system_comp_atom_range(&data.mold.sys, comp_idx);
                                        md_bitfield_set_range(&data.selection.highlight_mask, atom_range.beg, atom_range.end);
                                        handle_item_click(data);
                                    }
                                }
                            }
                        }
                        ImGui::Unindent();
                    }
                    ImGui::Unindent();
                }
                ImGui::Separator();
            }

            if (num_atom_types) {
                const float min_mass = 1.0f;
                const float max_mass = 500.0f;

                const float min_radius = 0.1f;
                const float max_radius = 20.0f;

                bool radius_changed = false;
                bool color_changed  = false;
                bool mass_changed   = false;

                if (ImGui::CollapsingHeader("Atom Types")) {
                    ImGui::Indent();
                    for (size_t i = 0; i < num_atom_types; ++i) {
                        DatasetItem& item = atom_types[i];
                        if (i == 0 && item.count == 0) {
                            // Skip sentinel "unknown" atom type if unused
                            continue;
                        }

                        ImGui::PushID((int)i);
                        defer{ ImGui::PopID(); };
                        char buf_num[32];
						snprintf(buf_num, sizeof(buf_num), "%d (%.2f%%)", item.count, item.fraction * 100.0f);
                        char buf_tot[256];
                        snprintf(buf_tot, sizeof(buf_tot), "%-4s %10s", item.label, buf_num);
                        bool expand = ImGui::CollapsingHeader(buf_tot);
                        if (ImGui::IsItemHovered()) {
                            md_bitfield_clear(&data.selection.highlight_mask);
                            md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, (uint32_t)md_array_size(item.indices));
                            if (ImGui::IsItemClicked()) {
                                handle_item_click(data);
                            }
                        }
                        if (expand) {
                            bool coarse_grained = data.mold.sys.atom.type.flags[i] & MD_FLAG_COARSE_GRAINED;
                            if (ImGui::Checkbox("Coarse Grained", &coarse_grained)) {
                                if (coarse_grained) {
                                    data.mold.sys.atom.type.flags[i] |=  MD_FLAG_COARSE_GRAINED;
                                } else {
                                    data.mold.sys.atom.type.flags[i] &= ~MD_FLAG_COARSE_GRAINED;
                                }
                            }
                            if (!(data.mold.sys.atom.type.flags[i] & MD_FLAG_COARSE_GRAINED)) {
								str_t symbol = md_atomic_number_symbol((md_atomic_number_t)data.mold.sys.atom.type.z[i]);
                                if (ImGui::Checkbox("Use element defaults", &item.use_defaults)) {
                                    if (item.use_defaults) {
										// If the user enables the use_defaults flag then we should set the values back to the element defaults
                                        md_element_t elem = data.mold.sys.atom.type.z[i];
										data.mold.sys.atom.type.radius[i] = md_util_element_vdw_radius(elem);
										data.mold.sys.atom.type.mass[i]   = md_util_element_atomic_mass(elem);
										data.mold.sys.atom.type.color[i]  = md_util_element_cpk_color(elem);
										radius_changed = true;
                                        color_changed = true;
										mass_changed = true;
                                    }
                                }
                                if (item.use_defaults) {
                                    ImGui::SameLine();
                                    int z = data.mold.sys.atom.type.z[i];
                                    vec4_t color = element_defaults[z].color;
                                    if (element_button(str_ptr(symbol), color)) {
                                        ImGui::OpenPopup("Element Popup");
                                    }
                                }
                                if (ImGui::BeginPopup("Element Popup")) {
                                    PeriodicTableResult table_res = periodic_table_widget(element_defaults);
                                    if (table_res.clicked) {
                                        int z = table_res.z;
                                        data.mold.sys.atom.type.z[i] = (md_atomic_number_t)table_res.z;
                                        data.mold.sys.atom.type.color[i] = u32_from_vec4(element_defaults[z].color);
                                        data.mold.sys.atom.type.radius[i] = element_defaults[z].radius;
                                        data.mold.sys.atom.type.mass[i] = element_defaults[z].mass;

										radius_changed = true;
										color_changed = true;
										mass_changed = true;
                                        ImGui::CloseCurrentPopup();
                                    }
									ImGui::EndPopup();
                                }
                            } else {
								item.use_defaults = false; // Coarse grained types always have custom properties
                            }

                            if (item.use_defaults) {
                                ImGui::PushDisabled();
                            }
                            float* radius = &data.mold.sys.atom.type.radius[i];
                            radius_changed |= ImGui::SliderFloat("Radius", radius, min_radius, max_radius);

                            float* mass = &data.mold.sys.atom.type.mass[i];
                            mass_changed |= ImGui::SliderFloat("Mass", mass, min_mass, max_mass);

                            ImVec4 color = ImColor(data.mold.sys.atom.type.color[i]);
                            if (ImGui::ColorEdit4("Color", &color.x)) {
                                data.mold.sys.atom.type.color[i] = ImColor(color);
                                color_changed = true;
                            }

                            if (item.use_defaults) {
                                ImGui::PopDisabled();
							}
                        }
                    }
                    ImGui::Unindent();
                }

                // Keep track of what elements are used in the system
				uint64_t elem_mask[2] = { 0 };

                for (size_t i = 0; i < num_atom_types; ++i) {
                    int z = data.mold.sys.atom.type.z[i];
                    if (atom_types[i].use_defaults && atom_types[i].count > 0) {
						elem_mask[z / 64] |= (1ULL << (z % 64));
                    }
                }

                if ((elem_mask[0] || elem_mask[1]) && ImGui::CollapsingHeader("Element Defaults")) {
                    ImGui::Indent();

                    static int z = -1;
					PeriodicTableResult table_res = periodic_table_widget(element_defaults, elem_mask);
					if (table_res.hovered) {
                        md_bitfield_clear(&data.selection.highlight_mask);
                        for (size_t i = 0; i < num_atom_types; ++i) {
							const DatasetItem& item = atom_types[i];
                            if (data.mold.sys.atom.type.z[i] == table_res.z) {
                                md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, (uint32_t)md_array_size(item.indices));
                            }
                        }
                        handle_item_click(data);

                        if (table_res.clicked && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                            ImGui::OpenPopup("Element Popup");
							z = table_res.z;
                        }
                    }

                    if (ImGui::BeginPopup("Element Popup")) {
                        str_t sym = md_atomic_number_symbol((md_atomic_number_t)z);
                        str_t name = md_atomic_number_name((md_atomic_number_t)z);
                        char buf[64];
                        snprintf(buf, sizeof(buf), "%d: %s (%s)", z, str_ptr(name), str_ptr(sym));
                        ImGui::Text("%s", buf);
                        ImGui::Separator();

                        ElementDefault& elem_def = element_defaults[z];

                        if (ImGui::ColorEdit3("Color", elem_def.color.elem)) {
                            // Iterate and set color for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.color[i] = u32_from_vec4(elem_def.color);
                                }
                            }
                            color_changed = true;
                        }
                        if (ImGui::InputFloat("Van der Waals Radius", &elem_def.radius)) {
                            // Iterate and set radius for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.radius[i] = elem_def.radius;
                                }
                            }
                            radius_changed = true;
                        }
                        if (ImGui::InputFloat("Atomic Mass", &elem_def.mass)) {
                            // Iterate and set mass for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.mass[i] = elem_def.mass;
                                }
                            }
                            mass_changed = true;
                        }

                        ImGui::EndPopup();
                    }
                    ImGui::Unindent();
                }

                if (radius_changed) {
                    data.mold.dirty_gpu_buffers |= MolBit_DirtyRadius;
                }

                if (color_changed) {
                    // @NOTE: Only the color within representations needs to be updated, not the filter.
                    flag_all_representations_as_dirty(&data);
                }
            }


            // Draw the three sections
            //draw_dataset_section("Chain Types",     inst_types,      md_array_size(inst_types),   0);
            //draw_dataset_section("Residue Types",   comp_types,    md_array_size(comp_types), 1);  
            //draw_dataset_section("Atom Types",      atom_types,       md_array_size(atom_types),    2);

            // Atom Element Mappings section (keep existing functionality)
            const size_t num_mappings = md_array_size(atom_element_remappings);
            if (num_mappings) {
                if (ImGui::CollapsingHeader("Atom Element Mappings")) {
                    for (size_t i = 0; i < num_mappings; ++i) {
                        const auto& mapping = atom_element_remappings[i];
                        ImGui::Text("%s -> %s (%s)", mapping.lbl, md_util_element_name(mapping.elem).ptr, md_util_element_symbol(mapping.elem).ptr);
                    }
                }
            }
        }
        ImGui::End();
    }
};

static Dataset dataset_instance;

}  // namespace dataset