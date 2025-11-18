// TREXIO UI Panel for VIAMD
// Provides file picker, summary window, and orbital grid evaluation interface

#ifdef MD_TREXIO

#define IMGUI_DEFINE_MATH_OPERATORS

#include "trexio_panel.h"
#include <trexio/trexio_orbital_grid.h>

#include <event.h>
#include <viamd.h>
#include <task_system.h>

#include <md_trexio.h>
#include <md_util.h>
#include <core/md_log.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>

#include <app/IconsFontAwesome6.h>

#include <nfd.h>

#include <ctime>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <atomic>

#define HARTREE_TO_EV 27.211386245988

namespace trexio_ui {

struct State {
    md_allocator_i* allocator = nullptr;
    md_trexio_t* trexio_data = nullptr;
    
    // UI state
    struct {
        bool show_main_window = false;
        bool show_summary = false;
        bool show_orbital_grid = false;
        bool verbose_logging = false;
    } ui;
    
    // File info
    char current_file[512] = "";
    std::vector<std::string> recent_files;
    
    // Summary data (cached for UI display)
    struct {
        int num_atoms = 0;
        int num_electrons = 0;
        int num_ao = 0;
        int num_mo = 0;
        int num_basis_shells = 0;
        int num_basis_prims = 0;
        std::vector<std::string> atom_labels;
        std::vector<double> atom_coords; // x,y,z for each atom
    } summary;
    
    // Orbital grid state
    struct {
        int selected_mo_index = 0;
        
        // Grid parameters
        struct {
            float min_x = -5.0f, max_x = 5.0f;
            float min_y = -5.0f, max_y = 5.0f;
            float min_z = -5.0f, max_z = 5.0f;
            int res_x = 32, res_y = 32, res_z = 32;
        } params;
        
        // Grid computation state
        bool computing = false;
        std::atomic<float> progress{0.0f};
        task_system::ID compute_task = task_system::INVALID_ID;
        std::atomic<bool> cancel_requested{false};
        
        // Grid data
        float* grid_data = nullptr;
        size_t grid_data_size = 0;
        
        // Grid results
        struct {
            bool available = false;
            float min_value = 0.0f;
            float max_value = 0.0f;
            float mean_value = 0.0f;
            char output_file[512] = "";
        } result;
        
        // Log buffer
        char log_buffer[4096] = "";
    } orbital_grid;
    
    // Debug logging
    FILE* debug_log_file = nullptr;
};

static State state = {};

void log_debug(const char* fmt, ...) {
    if (!state.ui.verbose_logging || !state.debug_log_file) {
        return;
    }
    
    va_list args;
    va_start(args, fmt);
    
    // Write timestamp
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char time_str[64];
    strftime(time_str, sizeof(time_str), "[%Y-%m-%d %H:%M:%S] ", timeinfo);
    fprintf(state.debug_log_file, "%s", time_str);
    
    // Write message
    vfprintf(state.debug_log_file, fmt, args);
    fprintf(state.debug_log_file, "\n");
    fflush(state.debug_log_file);
    
    va_end(args);
}

void initialize() {
    state.allocator = md_get_heap_allocator();
    log_debug("TREXIO UI panel initialized");
}

void shutdown() {
    if (state.trexio_data) {
        md_trexio_destroy(state.trexio_data);
        state.trexio_data = nullptr;
    }
    
    if (state.orbital_grid.grid_data) {
        md_free(state.allocator, state.orbital_grid.grid_data, state.orbital_grid.grid_data_size);
        state.orbital_grid.grid_data = nullptr;
        state.orbital_grid.grid_data_size = 0;
    }
    
    if (state.debug_log_file) {
        fclose(state.debug_log_file);
        state.debug_log_file = nullptr;
    }
    
    log_debug("TREXIO UI panel shutdown");
}

bool load_file(const char* filepath) {
    log_debug("Loading TREXIO file: %s", filepath);
    
    // Clean up existing data
    if (state.trexio_data) {
        md_trexio_destroy(state.trexio_data);
        state.trexio_data = nullptr;
    }
    
    // Create new TREXIO object
    state.trexio_data = md_trexio_create(state.allocator);
    if (!state.trexio_data) {
        MD_LOG_ERROR("Failed to create TREXIO object");
        log_debug("ERROR: Failed to create TREXIO object");
        return false;
    }
    
    // Parse file
    str_t filename_str = {filepath, (int64_t)strlen(filepath)};
    if (!md_trexio_parse_file(state.trexio_data, filename_str)) {
        MD_LOG_ERROR("Failed to load TREXIO file: %s", filepath);
        log_debug("ERROR: Failed to parse TREXIO file");
        md_trexio_destroy(state.trexio_data);
        state.trexio_data = nullptr;
        return false;
    }
    
    // Update file info
    strncpy(state.current_file, filepath, sizeof(state.current_file) - 1);
    state.current_file[sizeof(state.current_file) - 1] = '\0';
    
    // Cache summary data
    state.summary.num_atoms = (int)md_trexio_number_of_atoms(state.trexio_data);
    state.summary.num_electrons = (int)(md_trexio_num_up_electrons(state.trexio_data) + 
                                        md_trexio_num_down_electrons(state.trexio_data));
    state.summary.num_ao = (int)md_trexio_number_of_aos(state.trexio_data);
    state.summary.num_mo = (int)md_trexio_mo_num(state.trexio_data);
    state.summary.num_basis_shells = (int)md_trexio_basis_shell_num(state.trexio_data);
    state.summary.num_basis_prims = (int)md_trexio_basis_prim_num(state.trexio_data);
    
    // Cache atom labels and coordinates
    const char** labels = md_trexio_atom_labels(state.trexio_data);
    const double* coords = md_trexio_atom_coordinates(state.trexio_data);
    
    state.summary.atom_labels.clear();
    state.summary.atom_coords.clear();
    
    if (labels && coords) {
        for (int i = 0; i < state.summary.num_atoms; ++i) {
            state.summary.atom_labels.push_back(labels[i] ? labels[i] : "X");
            state.summary.atom_coords.push_back(coords[i * 3 + 0]);
            state.summary.atom_coords.push_back(coords[i * 3 + 1]);
            state.summary.atom_coords.push_back(coords[i * 3 + 2]);
        }
    }
    
    // Auto-suggest grid bounds from molecule
    if (coords && state.summary.num_atoms > 0) {
        float min_x = coords[0], max_x = coords[0];
        float min_y = coords[1], max_y = coords[1];
        float min_z = coords[2], max_z = coords[2];
        
        for (int i = 0; i < state.summary.num_atoms; ++i) {
            float x = coords[i * 3 + 0];
            float y = coords[i * 3 + 1];
            float z = coords[i * 3 + 2];
            
            if (x < min_x) min_x = x;
            if (x > max_x) max_x = x;
            if (y < min_y) min_y = y;
            if (y > max_y) max_y = y;
            if (z < min_z) min_z = z;
            if (z > max_z) max_z = z;
        }
        
        // Add padding (5 Angstrom)
        const float padding = 5.0f;
        state.orbital_grid.params.min_x = min_x - padding;
        state.orbital_grid.params.max_x = max_x + padding;
        state.orbital_grid.params.min_y = min_y - padding;
        state.orbital_grid.params.max_y = max_y + padding;
        state.orbital_grid.params.min_z = min_z - padding;
        state.orbital_grid.params.max_z = max_z + padding;
    }
    
    log_debug("Successfully loaded: %d atoms, %d MOs", state.summary.num_atoms, state.summary.num_mo);
    MD_LOG_INFO("Loaded TREXIO file: %s (%d atoms, %d MOs)", filepath, 
                state.summary.num_atoms, state.summary.num_mo);
    
    // Show summary window after successful load
    state.ui.show_summary = true;
    
    return true;
}

void draw_file_picker() {
    if (ImGui::Button(ICON_FA_FOLDER_OPEN " Open TREXIO File...")) {
        nfdchar_t* outPath = nullptr;
        nfdresult_t result = NFD_OpenDialog("h5,trexio", nullptr, &outPath);
        
        if (result == NFD_OKAY) {
            load_file(outPath);
            free(outPath);
        } else if (result != NFD_CANCEL) {
            MD_LOG_ERROR("File dialog error");
        }
    }
    
    if (state.current_file[0] != '\0') {
        ImGui::Text("Current file: %s", state.current_file);
    }
}

void draw_summary_window() {
    if (!state.ui.show_summary) return;
    
    if (ImGui::Begin(ICON_FA_INFO " TREXIO Summary", &state.ui.show_summary)) {
        if (!state.trexio_data) {
            ImGui::Text("No TREXIO file loaded");
        } else {
            ImGui::SeparatorText("System Information");
            ImGui::Text("Atoms: %d", state.summary.num_atoms);
            ImGui::Text("Electrons: %d", state.summary.num_electrons);
            
            ImGui::Spacing();
            ImGui::SeparatorText("Basis Set");
            ImGui::Text("Shells: %d", state.summary.num_basis_shells);
            ImGui::Text("Primitives: %d", state.summary.num_basis_prims);
            ImGui::Text("AO dims: %d", state.summary.num_ao);
            
            ImGui::Spacing();
            ImGui::SeparatorText("Molecular Orbitals");
            ImGui::Text("MO dims: %d", state.summary.num_mo);
            
            if (state.summary.num_mo > 0) {
                if (ImGui::Button(ICON_FA_ATOM " Open Orbital Grid Tool")) {
                    state.ui.show_orbital_grid = true;
                }
            }
            
            ImGui::Spacing();
            ImGui::SeparatorText("Atom List (first 50)");
            
            if (ImGui::BeginChild("AtomList", ImVec2(0, 200), true)) {
                int max_atoms = state.summary.num_atoms < 50 ? state.summary.num_atoms : 50;
                for (int i = 0; i < max_atoms; ++i) {
                    const char* label = state.summary.atom_labels[i].c_str();
                    double x = state.summary.atom_coords[i * 3 + 0];
                    double y = state.summary.atom_coords[i * 3 + 1];
                    double z = state.summary.atom_coords[i * 3 + 2];
                    
                    ImGui::Text("%3d: %-4s  %8.3f %8.3f %8.3f", i + 1, label, x, y, z);
                }
                if (state.summary.num_atoms > 50) {
                    ImGui::TextDisabled("... and %d more atoms", state.summary.num_atoms - 50);
                }
            }
            ImGui::EndChild();
        }
    }
    ImGui::End();
}

void draw_orbital_grid_window() {
    if (!state.ui.show_orbital_grid) return;
    
    if (ImGui::Begin(ICON_FA_ATOM " TREXIO Orbital Grid", &state.ui.show_orbital_grid)) {
        if (!state.trexio_data || state.summary.num_mo == 0) {
            ImGui::Text("No molecular orbital data available");
        } else {
            ImGui::SeparatorText("Orbital Selection");
            
            // MO index dropdown
            ImGui::SliderInt("MO Index", &state.orbital_grid.selected_mo_index, 
                           0, state.summary.num_mo - 1);
            
            ImGui::Spacing();
            ImGui::SeparatorText("Grid Parameters");
            
            ImGui::Text("Bounding Box:");
            ImGui::DragFloatRange2("X Range", &state.orbital_grid.params.min_x, 
                                  &state.orbital_grid.params.max_x, 0.1f, -50.0f, 50.0f);
            ImGui::DragFloatRange2("Y Range", &state.orbital_grid.params.min_y, 
                                  &state.orbital_grid.params.max_y, 0.1f, -50.0f, 50.0f);
            ImGui::DragFloatRange2("Z Range", &state.orbital_grid.params.min_z, 
                                  &state.orbital_grid.params.max_z, 0.1f, -50.0f, 50.0f);
            
            ImGui::Spacing();
            ImGui::Text("Resolution:");
            ImGui::SliderInt("Res X", &state.orbital_grid.params.res_x, 8, 128);
            ImGui::SliderInt("Res Y", &state.orbital_grid.params.res_y, 8, 128);
            ImGui::SliderInt("Res Z", &state.orbital_grid.params.res_z, 8, 128);
            
            ImGui::Spacing();
            ImGui::SeparatorText("Computation");
            
            if (state.orbital_grid.computing) {
                float prog = state.orbital_grid.progress.load();
                ImGui::ProgressBar(prog);
                ImGui::SameLine();
                ImGui::Text("%.0f%%", prog * 100.0f);
                
                if (ImGui::Button("Cancel")) {
                    state.orbital_grid.cancel_requested = true;
                    task_system::task_interrupt(state.orbital_grid.compute_task);
                }
            } else {
                if (ImGui::Button(ICON_FA_PLAY " Compute Grid")) {
                    // Clear previous results
                    state.orbital_grid.result.available = false;
                    state.orbital_grid.log_buffer[0] = '\0';
                    state.orbital_grid.cancel_requested = false;
                    
                    // Allocate grid memory
                    size_t grid_size = state.orbital_grid.params.res_x * 
                                      state.orbital_grid.params.res_y * 
                                      state.orbital_grid.params.res_z;
                    size_t grid_bytes = grid_size * sizeof(float);
                    
                    if (state.orbital_grid.grid_data) {
                        md_free(state.allocator, state.orbital_grid.grid_data, state.orbital_grid.grid_data_size);
                    }
                    
                    state.orbital_grid.grid_data = (float*)md_alloc(state.allocator, grid_bytes);
                    state.orbital_grid.grid_data_size = grid_bytes;
                    
                    if (!state.orbital_grid.grid_data) {
                        snprintf(state.orbital_grid.log_buffer, sizeof(state.orbital_grid.log_buffer),
                                "ERROR: Failed to allocate grid memory (%zu MB)\n",
                                grid_bytes / (1024 * 1024));
                        MD_LOG_ERROR("Failed to allocate grid memory");
                    } else {
                        state.orbital_grid.computing = true;
                        state.orbital_grid.progress = 0.0f;
                        
                        // Create background task
                        state.orbital_grid.compute_task = task_system::create_pool_task(
                            STR_LIT("Compute Orbital Grid"),
                            []() {
                                // Progress callback
                                auto progress_cb = [](float prog, void* user_data) -> bool {
                                    state.orbital_grid.progress = prog;
                                    return !state.orbital_grid.cancel_requested.load();
                                };
                                
                                // Log callback
                                auto log_cb = [](const char* msg, void* user_data) {
                                    size_t len = strlen(state.orbital_grid.log_buffer);
                                    size_t remaining = sizeof(state.orbital_grid.log_buffer) - len - 1;
                                    if (remaining > 0) {
                                        strncat(state.orbital_grid.log_buffer, msg, remaining);
                                        strncat(state.orbital_grid.log_buffer, "\n", remaining - strlen(msg));
                                    }
                                    log_debug("%s", msg);
                                };
                                
                                // Setup grid parameters
                                trexio_orbital_grid::BoundingBox bbox;
                                bbox.min_x = state.orbital_grid.params.min_x;
                                bbox.max_x = state.orbital_grid.params.max_x;
                                bbox.min_y = state.orbital_grid.params.min_y;
                                bbox.max_y = state.orbital_grid.params.max_y;
                                bbox.min_z = state.orbital_grid.params.min_z;
                                bbox.max_z = state.orbital_grid.params.max_z;
                                
                                trexio_orbital_grid::GridResolution resolution;
                                resolution.nx = state.orbital_grid.params.res_x;
                                resolution.ny = state.orbital_grid.params.res_y;
                                resolution.nz = state.orbital_grid.params.res_z;
                                
                                trexio_orbital_grid::GridStats stats;
                                
                                // Compute grid
                                bool success = trexio_orbital_grid::compute_orbital_grid(
                                    state.trexio_data,
                                    state.orbital_grid.selected_mo_index,
                                    bbox,
                                    resolution,
                                    state.orbital_grid.grid_data,
                                    &stats,
                                    progress_cb,
                                    log_cb,
                                    nullptr
                                );
                                
                                if (success && !state.orbital_grid.cancel_requested) {
                                    // Save results
                                    state.orbital_grid.result.min_value = stats.min_value;
                                    state.orbital_grid.result.max_value = stats.max_value;
                                    state.orbital_grid.result.mean_value = stats.mean_value;
                                    
                                    // Generate output filename
                                    time_t now = time(nullptr);
                                    struct tm* timeinfo = localtime(&now);
                                    char timestamp[64];
                                    strftime(timestamp, sizeof(timestamp), "%Y%m%d-%H%M%S", timeinfo);
                                    
                                    snprintf(state.orbital_grid.result.output_file,
                                            sizeof(state.orbital_grid.result.output_file),
                                            "orbital_mo%d_%s.grid",
                                            state.orbital_grid.selected_mo_index,
                                            timestamp);
                                    
                                    // Save grid to file
                                    trexio_orbital_grid::save_grid_to_file(
                                        state.orbital_grid.result.output_file,
                                        state.orbital_grid.grid_data,
                                        resolution,
                                        bbox
                                    );
                                    
                                    state.orbital_grid.result.available = true;
                                    
                                    char msg[256];
                                    snprintf(msg, sizeof(msg), "Grid saved to: %s", 
                                            state.orbital_grid.result.output_file);
                                    log_cb(msg, nullptr);
                                } else if (state.orbital_grid.cancel_requested) {
                                    log_cb("Computation cancelled by user", nullptr);
                                }
                                
                                state.orbital_grid.computing = false;
                                state.orbital_grid.progress = 1.0f;
                            }
                        );
                        
                        task_system::enqueue_task(state.orbital_grid.compute_task);
                    }
                }
            }
            
            // Results
            if (state.orbital_grid.result.available) {
                ImGui::Spacing();
                ImGui::SeparatorText("Grid Statistics");
                ImGui::Text("Min value: %.6e", state.orbital_grid.result.min_value);
                ImGui::Text("Max value: %.6e", state.orbital_grid.result.max_value);
                ImGui::Text("Mean value: %.6e", state.orbital_grid.result.mean_value);
                ImGui::Text("Output file: %s", state.orbital_grid.result.output_file);
            }
            
            // Log window
            ImGui::Spacing();
            ImGui::SeparatorText("Log");
            ImGui::TextWrapped("%s", state.orbital_grid.log_buffer);
        }
    }
    ImGui::End();
}

void draw_main_window() {
    if (!state.ui.show_main_window) return;
    
    if (ImGui::Begin(ICON_FA_FLASK " TREXIO Tools", &state.ui.show_main_window)) {
        draw_file_picker();
        
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        
        if (ImGui::Checkbox("Verbose Logging", &state.ui.verbose_logging)) {
            if (state.ui.verbose_logging && !state.debug_log_file) {
                // Open debug log file
                char filename[512];
                time_t now = time(nullptr);
                struct tm* timeinfo = localtime(&now);
                char timestamp[64];
                strftime(timestamp, sizeof(timestamp), "%Y%m%d-%H%M%S", timeinfo);
                snprintf(filename, sizeof(filename), "tools/trexio/debug-%s.log", timestamp);
                
                state.debug_log_file = fopen(filename, "w");
                if (state.debug_log_file) {
                    log_debug("=== TREXIO Debug Log Started ===");
                }
            } else if (!state.ui.verbose_logging && state.debug_log_file) {
                log_debug("=== TREXIO Debug Log Ended ===");
                fclose(state.debug_log_file);
                state.debug_log_file = nullptr;
            }
        }
        
        ImGui::Spacing();
        
        if (ImGui::Button(ICON_FA_INFO " Show Summary")) {
            state.ui.show_summary = true;
        }
        
        if (state.summary.num_mo > 0) {
            ImGui::SameLine();
            if (ImGui::Button(ICON_FA_ATOM " Show Orbital Grid")) {
                state.ui.show_orbital_grid = true;
            }
        }
    }
    ImGui::End();
}

void set_main_window_visible(bool visible) {
    state.ui.show_main_window = visible;
}

void update() {
    draw_main_window();
    draw_summary_window();
    draw_orbital_grid_window();
}

} // namespace trexio_ui

#endif // MD_TREXIO
