// TREXIO Component for VIAMD
// Provides visualization of quantum chemistry data from TREXIO files
// Reuses existing infrastructure for orbital and summary displays

#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>

#ifdef MD_TREXIO
#include <md_trexio.h>
#endif

#include <md_util.h>
#include <core/md_log.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>

#include <app/IconsFontAwesome6.h>

#define HARTREE_TO_EV 27.211386245988

struct Trexio {
    md_allocator_i* allocator = NULL;
    
#ifdef MD_TREXIO
    md_trexio_t* trexio_data = NULL;
#endif

    struct {
        int num_atoms = 0;
        int num_electrons = 0;
        int num_mo = 0;
        int homo_idx = -1;
        int lumo_idx = -1;
        
        double* mo_energies = NULL;
        double* mo_occupations = NULL;
    } qc_data;

    struct {
        bool show_summary = true;
        bool show_orbital = false;
    } ui;

    bool init(md_allocator_i* alloc) {
        allocator = alloc;
        
#ifdef MD_TREXIO
        MD_LOG_INFO("TREXIO Component initialized");
        return true;
#else
        MD_LOG_INFO("TREXIO Component: Compiled without TREXIO support");
        return false;
#endif
    }

    void shutdown() {
#ifdef MD_TREXIO
        if (trexio_data) {
            md_trexio_free(trexio_data);
            trexio_data = NULL;
        }
#endif
    }

    bool load_trexio_file(const char* filename) {
#ifdef MD_TREXIO
        if (trexio_data) {
            md_trexio_free(trexio_data);
            trexio_data = NULL;
        }

        trexio_data = md_trexio_parse_file(filename, allocator);
        if (!trexio_data) {
            MD_LOG_ERROR("Failed to load TREXIO file: %s", filename);
            return false;
        }

        qc_data.num_atoms = (int)trexio_data->num_atoms;
        qc_data.num_electrons = trexio_data->num_electrons_up + trexio_data->num_electrons_down;
        
        // TODO: Extract MO data when available
        qc_data.num_mo = 0;
        qc_data.homo_idx = -1;
        qc_data.lumo_idx = -1;

        MD_LOG_INFO("Loaded TREXIO file: %s (%d atoms, %d electrons)", 
                    filename, qc_data.num_atoms, qc_data.num_electrons);
        
        return true;
#else
        MD_LOG_ERROR("TREXIO support not compiled in");
        return false;
#endif
    }

    void draw_summary_window() {
        if (!ui.show_summary) return;

        if (ImGui::Begin(ICON_FA_INFO " TREXIO Summary", &ui.show_summary)) {
#ifdef MD_TREXIO
            if (trexio_data) {
                ImGui::SeparatorText("System Information");
                ImGui::Text("Atoms: %d", qc_data.num_atoms);
                ImGui::Text("Electrons: %d", qc_data.num_electrons);
                ImGui::Text("Molecular Orbitals: %d", qc_data.num_mo);
                
                if (qc_data.homo_idx >= 0) {
                    ImGui::Text("HOMO: MO %d", qc_data.homo_idx);
                }
                if (qc_data.lumo_idx >= 0) {
                    ImGui::Text("LUMO: MO %d", qc_data.lumo_idx);
                }
                
                ImGui::Spacing();
                ImGui::Checkbox("Show Orbital Viewer", &ui.show_orbital);
                
            } else {
                ImGui::Text("No TREXIO data loaded");
            }
#else
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), 
                "TREXIO support not compiled in");
            ImGui::Text("Rebuild with -DVIAMD_ENABLE_TREXIO=ON");
#endif
        }
        ImGui::End();
    }

    void draw_orbital_window() {
        if (!ui.show_orbital) return;

        if (ImGui::Begin(ICON_FA_ATOM " TREXIO Orbitals", &ui.show_orbital)) {
#ifdef MD_TREXIO
            if (trexio_data && qc_data.num_mo > 0) {
                ImGui::Text("Orbital visualization - Coming soon");
                ImGui::Text("Will reuse VeloxChem grid orbital renderer");
            } else {
                ImGui::Text("No molecular orbital data available");
            }
#else
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), 
                "TREXIO support not compiled in");
#endif
        }
        ImGui::End();
    }

    void update() {
        draw_summary_window();
        draw_orbital_window();
    }
};

static Trexio instance = {};
