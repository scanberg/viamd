#include <event.h>
#include <viamd.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <md_molecule.h>
#include <md_util.h>

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

struct Dataset : viamd::EventHandler {
    bool show_window = false;
    
    // Dataset data (moved from ApplicationState.dataset)
    md_array(AtomElementMapping) atom_element_remappings = 0;
    md_array(DatasetItem) inst_types = 0;
    md_array(DatasetItem) comp_types = 0;
    md_array(DatasetItem) atom_types = 0;
    md_allocator_i* arena = 0;

    Dataset() { 
        viamd::event_system_register_handler(*this); 
    }
    
    ~Dataset() {
        // The arena will be cleaned up when the persistent allocator is destroyed
        // No need for explicit cleanup here
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

        size_t temp_pos = md_temp_get_pos();
        md_allocator_i* temp_alloc = md_get_temp_allocator();
        md_array(int) sequence = 0;
        md_array(int) comp_idx_type = 0;

        if (comp_count > 0) {
			md_array_resize(comp_idx_type, comp_count, temp_alloc);
			MEMSET(comp_idx_type, -1, md_array_bytes(comp_idx_type));
        }

        // Process components - group by name + atom type sequence
        for (size_t i = 0; i < comp_count; ++i) {
            str_t comp_name = md_comp_name(&sys.comp, i);

            md_array_shrink(sequence, 0);
            md_urange_t range = md_comp_atom_range(&sys.comp, i);
            for (int j = range.beg; j < range.end; ++j) {
                int ai = sys.atom.type_idx[j];
                md_array_push(sequence, ai, temp_alloc);
            }

            // Create combined string for hash (label + sequence of atom types)
            uint64_t hash = md_hash64_str(comp_name, md_hash64(sequence, md_array_bytes(sequence), 0));

            // Check if we already have this chain type
            DatasetItem* item = nullptr;
            for (size_t j = 0; j < md_array_size(comp_types); ++j) {
                if (comp_types[j].key == hash) {
                    item = &comp_types[j];
                    md_array_push(comp_idx_type, (int)j, temp_alloc);
                    break;
                }
            }
            if (!item) {
                DatasetItem it = { .key = hash };
                snprintf(it.label, sizeof(it.label), STR_FMT, STR_ARG(comp_name));
                md_array_push(comp_idx_type, (int)md_array_size(comp_types), temp_alloc);
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
                    md_array_push(sequence, res_type_idx, temp_alloc);
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

        md_temp_set_pos_back(temp_pos);
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                // Initialize component
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

    void draw(ApplicationState& data) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(500, 600), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Dataset", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::Text("System: %s", data.files.molecule);
            ImGui::Text("Num entities:   %9zu", data.mold.sys.entity.count);
            ImGui::Text("Num instances:  %9zu", data.mold.sys.inst.count);
            ImGui::Text("Num components: %9zu", data.mold.sys.comp.count);
            ImGui::Text("Num atoms:      %9zu", data.mold.sys.atom.count);

            if (data.files.trajectory[0] != '\0') {
                ImGui::Separator();
                ImGui::Text("Trajectory data: %s", data.files.trajectory);
                ImGui::Text("Num frames:    %9zu", md_trajectory_num_frames(data.mold.traj));
                ImGui::Text("Num atoms:     %9zu", md_trajectory_num_atoms(data.mold.traj));
            }

            ImGui::Separator();

            if (ImGui::IsWindowHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
            }

            static bool use_short_labels = true;
            ImGui::Checkbox("Use short labels", &use_short_labels);

            size_t num_entities = md_system_entity_count(&data.mold.sys);
            size_t num_instances = md_system_inst_count(&data.mold.sys);

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
                    }
                    if (expand_entity) {
                        ImGui::Indent();
                        for (size_t inst_idx = 0; inst_idx < num_instances; ++inst_idx) {
                            if (md_inst_entity_idx(&data.mold.sys.inst, inst_idx) != (int)ent_idx) continue;
                            str_t inst_id = md_inst_id(&data.mold.sys.inst, inst_idx);
                            snprintf(buf, sizeof(buf), STR_FMT, STR_ARG(inst_id));
                            bool expand_inst = ImGui::CollapsingHeader(buf);
                            if (ImGui::IsItemHovered()) {
                                md_bitfield_clear(&data.selection.highlight_mask);
                                md_urange_t range = md_system_inst_atom_range(&data.mold.sys, inst_idx);
                                md_bitfield_set_range(&data.selection.highlight_mask, range.beg, range.end);
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
                                        //ImGui::PushStyleColor(ImGuiCol_HeaderHovered, color);

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