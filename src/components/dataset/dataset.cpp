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

struct Dataset : viamd::EventHandler {
    bool show_window = false;

    Dataset() { 
        viamd::event_system_register_handler(*this); 
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
                break;
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
            auto draw_dataset_section = [&data](const char* title, DatasetItem* items, size_t item_count, int section_type) {
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
                                
                                // Show subcomponents as buttons (residues in this chain)
                                size_t seq_len = md_array_size(item.sub_items);
                                if (seq_len > 0) {
                                    int button_count = 0;
                                    for (size_t i = 0; i < seq_len; ++i) {
                                        const DatasetItem& sub = data.dataset.residue_types[item.sub_items[i]];
                                        // Possibly do not use the same color for all subitems
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
                            } else if (section_type == 1) { // Residues
                                ImGui::Text("Residue: %s", item.label);
                                ImGui::Separator();                            
                                
                                // Show subcomponents as buttons (atoms in this residue)
                                size_t seq_len = md_array_size(item.sub_items);
                                if (seq_len > 0) {
                                    int button_count = 0;
                                    for (size_t i = 0; i < seq_len; ++i) {
                                        const DatasetItem& sub = data.dataset.atom_types[item.sub_items[i]];

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
            draw_dataset_section("Chain Types",     data.dataset.chain_types,      md_array_size(data.dataset.chain_types),   0);
            draw_dataset_section("Residue Types",   data.dataset.residue_types,    md_array_size(data.dataset.residue_types), 1);  
            draw_dataset_section("Atom Types",      data.dataset.atom_types,       md_array_size(data.dataset.atom_types),    2);

            // Atom Element Mappings section (keep existing functionality)
            const size_t num_mappings = md_array_size(data.dataset.atom_element_remappings);
            if (num_mappings) {
                if (ImGui::CollapsingHeader("Atom Element Mappings")) {
                    for (size_t i = 0; i < num_mappings; ++i) {
                        const auto& mapping = data.dataset.atom_element_remappings[i];
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