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

// Helper function to convert 3-letter amino acid codes to single-letter codes
static char amino_acid_to_single_letter(const char* aa3) {
    if (!aa3) return 'X';
    
    // Standard amino acids
    if (strcmp(aa3, "ALA") == 0) return 'A';
    if (strcmp(aa3, "ARG") == 0) return 'R';
    if (strcmp(aa3, "ASN") == 0) return 'N';
    if (strcmp(aa3, "ASP") == 0) return 'D';
    if (strcmp(aa3, "CYS") == 0) return 'C';
    if (strcmp(aa3, "GLU") == 0) return 'E';
    if (strcmp(aa3, "GLN") == 0) return 'Q';
    if (strcmp(aa3, "GLY") == 0) return 'G';
    if (strcmp(aa3, "HIS") == 0) return 'H';
    if (strcmp(aa3, "ILE") == 0) return 'I';
    if (strcmp(aa3, "LEU") == 0) return 'L';
    if (strcmp(aa3, "LYS") == 0) return 'K';
    if (strcmp(aa3, "MET") == 0) return 'M';
    if (strcmp(aa3, "PHE") == 0) return 'F';
    if (strcmp(aa3, "PRO") == 0) return 'P';
    if (strcmp(aa3, "SER") == 0) return 'S';
    if (strcmp(aa3, "THR") == 0) return 'T';
    if (strcmp(aa3, "TRP") == 0) return 'W';
    if (strcmp(aa3, "TYR") == 0) return 'Y';
    if (strcmp(aa3, "VAL") == 0) return 'V';
    
    // Non-standard/modified amino acids - return 'X' for unknown
    return 'X';
}

struct Dataset : viamd::EventHandler {
    bool show_window = false;
    
    // Dataset data (moved from ApplicationState.dataset)
    md_array(AtomElementMapping) atom_element_remappings = 0;
    md_array(DatasetItem) chain_types = 0;
    md_array(DatasetItem) residue_types = 0;
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
        chain_types = 0;
        residue_types = 0;
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

        const md_molecule_t& mol = data.mold.mol;
        size_t atom_count = md_atom_count(&mol.atom);
        if (atom_count == 0) return;

        size_t atom_type_count = md_atom_type_count(&mol.atom.type);

        // Map atom types into dataset items
        for (size_t i = 0; i < atom_type_count; ++i) {
            str_t atom_type_name = md_atom_type_name(&mol.atom.type, i);
            DatasetItem item = { .key = i };
            snprintf(item.label, sizeof(item.label), STR_FMT, STR_ARG(atom_type_name));
            md_array_push(atom_types, item, arena);
        }

        // Count and set indices for each atom type
        for (size_t i = 0; i < atom_count; ++i) {
            md_atom_type_idx_t type_idx = mol.atom.type_idx[i]; 
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
        md_array(int) residue_idx_to_type = 0;

        // Process residues - group by name + atom type sequence
        size_t res_count = md_residue_count(&mol.residue);
        for (size_t i = 0; i < res_count; ++i) {
            str_t res_name = md_residue_name(&mol.residue, i);

            md_array_shrink(sequence, 0);
            md_range_t range = md_residue_atom_range(&mol.residue, i);
            for (int j = range.beg; j < range.end; ++j) {
                int ai = mol.atom.type_idx[j];
                md_array_push(sequence, ai, temp_alloc);
            }

            // Create combined string for hash (label + sequence of atom types)
            uint64_t hash = md_hash64_str(res_name, md_hash64(sequence, md_array_bytes(sequence), 0));

            // Check if we already have this chain type
            DatasetItem* item = nullptr;
            for (size_t j = 0; j < md_array_size(residue_types); ++j) {
                if (residue_types[j].key == hash) {
                    item = &residue_types[j];
                    md_array_push(residue_idx_to_type, (int)j, temp_alloc);
                    break;
                }
            }
            if (!item) {
                DatasetItem it = { .key = hash };
                snprintf(it.label, sizeof(it.label), STR_FMT, STR_ARG(res_name));
                md_array_push(residue_idx_to_type, (int)md_array_size(residue_types), temp_alloc);
                md_array_push(residue_types, it, arena);
                item = md_array_last(residue_types);
                md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
            }

            size_t res_atom_count = md_residue_atom_count(&mol.residue, i);

            item->count += 1;
            item->fraction += (float)(res_atom_count / (double)atom_count);
            md_array_push(item->indices, (int)i, arena);
        }

        // Process chains - group by label + residue sequence
        size_t chain_count = md_chain_count(&mol.chain);
        for (size_t i = 0; i < chain_count; ++i) {
            str_t chain_id = md_chain_id(&mol.chain, i);

            md_array_shrink(sequence, 0);
            md_range_t range = md_chain_residue_range(&mol.chain, i);
            for (int j = range.beg; j < range.end; ++j) {
                int res_type_idx = residue_idx_to_type[j];
                md_array_push(sequence, res_type_idx, temp_alloc);
            }

            // Create combined string for hash (label + sequence of residue types)
            uint64_t hash = md_hash64_str(chain_id,  md_hash64(sequence, md_array_bytes(sequence), 0));
            
            // Check if we already have this chain type
            DatasetItem* item = nullptr;
            for (size_t j = 0; j < md_array_size(chain_types); ++j) {
                if (chain_types[j].key == hash) {
                    item = &chain_types[j];
                    break;
                }
            }
            if (!item) {
                DatasetItem it = { .key = hash };
                snprintf(it.label, sizeof(it.label), STR_FMT, STR_ARG(chain_id));
                md_array_push(chain_types, it, arena);
                item = md_array_last(chain_types);
                md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
            }
            
            size_t chain_atom_count = md_chain_atom_count(&mol.chain, i);

            item->count += 1;
            item->fraction += (float)(chain_atom_count / (double)atom_count);
            md_array_push(item->indices, (int)i, arena);
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
            ImGui::Text("Molecular data: %s", data.files.molecule);
            ImGui::Text("Num atoms:    %9i", (int)data.mold.mol.atom.count);
            ImGui::Text("Num residues: %9i", (int)data.mold.mol.residue.count);
            ImGui::Text("Num chains:   %9i", (int)data.mold.mol.chain.count);

            if (data.files.trajectory[0] != '\0') {
                ImGui::Separator();
                ImGui::Text("Trajectory data: %s", data.files.trajectory);
                ImGui::Text("Num frames:    %9i", (int)md_trajectory_num_frames(data.mold.traj));
                ImGui::Text("Num atoms:     %9i", (int)md_trajectory_num_atoms(data.mold.traj));
            }

            ImGui::Separator();

            if (ImGui::IsWindowHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
            }
            
            // Helper lambda for displaying dataset sections
            auto draw_dataset_section = [this, &data](const char* title, DatasetItem* items, size_t item_count, int section_type) {
                const size_t count = item_count;
                if (count && ImGui::CollapsingHeader(title, ImGuiTreeNodeFlags_DefaultOpen)) {
                    const ImVec2 item_size = ImVec2(ImGui::GetFontSize() * 1.8f, ImGui::GetFontSize() * 1.1f);
                    const float window_x_max = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;

                    auto handle_item_hover = [&data](const DatasetItem& item, int section_type) {
                        if (ImGui::IsItemHovered()) {
                            ImGui::SetTooltip("%s: count %d (%.2f%%)", item.label, item.count, item.fraction * 100.f);
                            switch (section_type) {
                            case 0: {
                                // Chains
                                for (int* it = md_array_beg(item.indices); it != md_array_end(item.indices); ++it) {
                                    int chain_idx = *it;
                                    md_range_t range = md_chain_atom_range(&data.mold.mol.chain, chain_idx);
                                    md_bitfield_set_range(&data.selection.highlight_mask, range.beg, range.end);
                                }
                                break;
                            }
                            case 1:
                                // Residues
                                for (int* it = md_array_beg(item.indices); it != md_array_end(item.indices); ++it) {
                                    int res_idx = *it;
                                    md_range_t range = md_residue_atom_range(&data.mold.mol.residue, res_idx);
                                    md_bitfield_set_range(&data.selection.highlight_mask, range.beg, range.end);
                                }
                                break;
                            case 2: // Atom Types
                                md_bitfield_set_indices_u32(&data.selection.highlight_mask, (const uint32_t*)item.indices, md_array_size(item.indices));
                                break;
                            }

                            //Select
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_or_inplace(&data.selection.selection_mask, &data.selection.highlight_mask);
                            }
                            //Deselect
                            else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_andnot_inplace(&data.selection.selection_mask, &data.selection.highlight_mask);
                            }
                        }
                    };
                    
                    for (size_t j = 0; j < count; ++j) {
                        const DatasetItem& item = items[j];
                        const float t = powf(item.fraction, 0.2f) * 0.5f;

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImVec4(1, 1, 0.5, 0.3));
                        ImGui::PushStyleColor(ImGuiCol_Header, ImPlot::SampleColormap(t, ImPlotColormap_Plasma));
                        
                        ImGui::Selectable(item.label, false, 0, item_size);
                        ImGui::PopStyleColor(2);

                        handle_item_hover(item, section_type);
                        
                        // Right-click popup
                        if (!ImGui::IsKeyDown(ImGuiKey_LeftShift) && ImGui::IsItemClicked(ImGuiMouseButton_Right)) {
                            ImGui::OpenPopup(("##popup_" + std::to_string(section_type) + "_" + std::to_string(j)).c_str());
                        }
                        
                        std::string popup_id = "##popup_" + std::to_string(section_type) + "_" + std::to_string(j);
                        if (ImGui::BeginPopup(popup_id.c_str())) {

                            if (ImGui::IsWindowHovered()) {
                                md_bitfield_clear(&data.selection.highlight_mask);
                            }

                            if (section_type == 0) { // Chains
                                ImGui::Text("Chain: %s", item.label);
                                ImGui::Separator();
                                
                                // First, get the actual chain indices to access residue sequence
                                if (md_array_size(item.indices) > 0) {
                                    int chain_idx = item.indices[0]; // Get first chain of this type
                                    md_range_t chain_res_range = md_chain_residue_range(&data.mold.mol.chain, chain_idx);
                                    
                                    // 1. Show unique residue types (not sequence)
                                    ImGui::Text("Unique Residue Types:");
                                    size_t seq_len = md_array_size(item.sub_items);
                                    if (seq_len > 0) {
                                        int button_count = 0;
                                        for (size_t i = 0; i < seq_len; ++i) {
                                            const DatasetItem& sub = residue_types[item.sub_items[i]];
                                            const float button_t = 0.3f;
                                            ImGui::PushStyleColor(ImGuiCol_Button, ImPlot::SampleColormap(button_t, ImPlotColormap_Plasma));
                                            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImPlot::SampleColormap(button_t + 0.1f, ImPlotColormap_Plasma));
                                            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImPlot::SampleColormap(button_t + 0.2f, ImPlotColormap_Plasma));
                                            
                                            ImGui::SmallButton(sub.label);
                                            handle_item_hover(sub, section_type + 1);
                                            ImGui::PopStyleColor(3);
                                            
                                            if ((button_count + 1) % 6 != 0) { // 6 buttons per row
                                                ImGui::SameLine();
                                            }
                                            button_count++;
                                        }
                                    }
                                    
                                    ImGui::Separator();
                                    
                                    // 2. Show sequence using single-letter codes with hover highlighting
                                    ImGui::Text("Sequence:");
                                    
                                    // Build the sequence string and create hoverable buttons
                                    if (chain_res_range.end > chain_res_range.beg) {
                                        ImGui::BeginGroup();
                                        for (int res_idx = chain_res_range.beg; res_idx < chain_res_range.end; ++res_idx) {
                                            str_t res_name = LBL_TO_STR(data.mold.mol.residue.name[res_idx]);
                                            
                                            // Convert to single-letter code
                                            char aa_code[2] = {amino_acid_to_single_letter(res_name.ptr), '\0'};
                                            
                                            // Create a small hoverable button for each residue
                                            ImGui::PushID(res_idx);
                                            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
                                            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(1, 1));
                                            
                                            // Use different color for each residue type
                                            float res_t = (float)(res_idx - chain_res_range.beg) / (float)(chain_res_range.end - chain_res_range.beg) * 0.8f + 0.1f;
                                            ImGui::PushStyleColor(ImGuiCol_Button, ImPlot::SampleColormap(res_t, ImPlotColormap_Viridis));
                                            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImPlot::SampleColormap(res_t + 0.1f, ImPlotColormap_Viridis));
                                            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImPlot::SampleColormap(res_t + 0.2f, ImPlotColormap_Viridis));
                                            
                                            if (ImGui::SmallButton(aa_code)) {
                                                // Could implement additional functionality here
                                            }
                                            
                                            // Highlight the corresponding residue when hovering
                                            if (ImGui::IsItemHovered()) {
                                                ImGui::SetTooltip("%.*s (residue %d)", (int)res_name.len, res_name.ptr, res_idx + 1);
                                                // Highlight this specific residue
                                                md_range_t res_atom_range = md_residue_atom_range(&data.mold.mol.residue, res_idx);
                                                md_bitfield_set_range(&data.selection.highlight_mask, res_atom_range.beg, res_atom_range.end);
                                            }
                                            
                                            ImGui::PopStyleColor(3);
                                            ImGui::PopStyleVar(2);
                                            ImGui::PopID();
                                            
                                            // Arrange in rows of ~20 amino acids
                                            if ((res_idx - chain_res_range.beg + 1) % 20 != 0 && res_idx + 1 < chain_res_range.end) {
                                                ImGui::SameLine();
                                            }
                                        }
                                        ImGui::EndGroup();
                                    }
                                } // Close the if (md_array_size(item.indices) > 0) block
                                
                                ImGui::Separator();
                            } else if (section_type == 1) { // Residues
                                ImGui::Text("Residue: %s", item.label);
                                ImGui::Separator();                            
                                
                                // Show subcomponents as buttons (atoms in this residue)
                                size_t seq_len = md_array_size(item.sub_items);
                                if (seq_len > 0) {
                                    int button_count = 0;
                                    for (size_t i = 0; i < seq_len; ++i) {
                                        const DatasetItem& sub = atom_types[item.sub_items[i]];

                                        const float button_t = 0.5f; // Use consistent color for atom buttons
                                        ImGui::PushStyleColor(ImGuiCol_Button, ImPlot::SampleColormap(button_t, ImPlotColormap_Plasma));
                                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImPlot::SampleColormap(button_t + 0.1f, ImPlotColormap_Plasma));
                                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImPlot::SampleColormap(button_t + 0.2f, ImPlotColormap_Plasma));
                                        
                                        ImGui::SmallButton(sub.label);
                                        handle_item_hover(sub, section_type + 1);
                                        ImGui::PopStyleColor(3);
                                        
                                        if ((button_count + 1) % 8 != 0) { // 8 buttons per row
                                            ImGui::SameLine();
                                        }
                                        button_count++;
                                    }
                                }
                            } else if (section_type == 2) { // Atom Types
                                ImGui::Text("Atom Type: %s", item.label);
                                ImGui::Separator();
                                
                                // Create a copy for editing
                                static float  editing_item_radius = 0;
                                static float  editing_item_mass = 0;
                                static int    editing_item_z = 0;
                                static int    editing_item_atom_type_idx = -1;
                                static vec4_t editing_item_color = vec4_set(1, 1, 1, 1);

                                if (editing_item_atom_type_idx != item.key) {
                                    editing_item_atom_type_idx  = item.key;
                                    editing_item_radius = md_atom_type_radius(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                    editing_item_mass = md_atom_type_mass(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                    editing_item_z = md_atom_type_atomic_number(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                }
                                
                                ImGui::Text("Properties:");
                                ImGui::InputFloat("Radius", &editing_item_radius, 0.01f, 0.1f, "%.3f");
                                ImGui::InputFloat("Mass",   &editing_item_mass,   0.01f, 1.0f, "%.3f");
                                
                                // Element dropdown
                                char lbl[64];
                                snprintf(lbl, sizeof(lbl), STR_FMT " (" STR_FMT ")", STR_ARG(md_util_element_name(editing_item_z)), STR_ARG(md_util_element_symbol(editing_item_z)));
                                if (ImGui::BeginCombo("Element", lbl)) {
                                    for (int i = 0; i <= 118; ++i) {
                                        str_t name = md_util_element_name((int)i);
                                        str_t symbol = md_util_element_symbol((int)i);
                                        snprintf(lbl, sizeof(lbl), STR_FMT " (" STR_FMT ")", STR_ARG(name), STR_ARG(symbol));
                                        if (ImGui::Selectable(lbl, editing_item_z == (int)i)) {
                                            editing_item_z = (int)i;
                                        }
                                    }
                                    ImGui::EndCombo();
                                }
                                
                                if (ImGui::Button("Apply Changes")) {
                                    data.mold.mol.atom.type.mass[editing_item_atom_type_idx] = editing_item_mass;
                                    data.mold.mol.atom.type.radius[editing_item_atom_type_idx] = editing_item_radius;
                                    data.mold.mol.atom.type.z[editing_item_atom_type_idx] = (md_atomic_number_t)editing_item_z;
                                    data.mold.dirty_buffers |= MolBit_DirtyRadius;
                                    ImGui::CloseCurrentPopup();
                                }
                                ImGui::SameLine();
                                if (ImGui::Button("Reset")) {
                                    editing_item_radius = md_atom_type_radius(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                    editing_item_mass = md_atom_type_mass(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                    editing_item_z = md_atom_type_atomic_number(&data.mold.mol.atom.type, editing_item_atom_type_idx);
                                }
                            }
                            ImGui::EndPopup();
                        }

                        // Arrange buttons horizontally when possible (like original)
                        float last_item_x = ImGui::GetItemRectMax().x;
                        float next_button_x = last_item_x + item_size.x;
                        if (j + 1 < count && next_button_x < window_x_max) {
                            ImGui::SameLine();
                        }
                    }
                }
            };

            // Draw the three sections
            draw_dataset_section("Chain Types",     chain_types,      md_array_size(chain_types),   0);
            draw_dataset_section("Residue Types",   residue_types,    md_array_size(residue_types), 1);  
            draw_dataset_section("Atom Types",      atom_types,       md_array_size(atom_types),    2);

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