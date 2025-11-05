// TREXIO Loader Component for VIAMD
// This component provides support for loading TREXIO files into VIAMD
//
// TREXIO (https://github.com/TREX-CoE/trexio) is a file format and library
// for quantum chemistry data developed by the TREX European Centre of Excellence.
//
// Current Implementation Status:
// - Basic file recognition and component structure implemented
// - Full molecule/wavefunction loading requires mdlib integration (see TODO below)
//
// TODO: Complete Implementation (see docs/TREXIO_INTEGRATION.md for details)
// 1. Create md_trexio.c/.h in ext/mdlib/src/ following md_vlx pattern
// 2. Implement molecule loader interface: init_from_file()
// 3. Implement basis set and orbital data structures
// 4. Add visualization support similar to VeloxChem component

#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <core/md_vec_math.h>
#include <core/md_log.h>
#include <core/md_arena_allocator.h>
#include <core/md_util.h>
#include <md_molecule.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#ifdef MD_TREXIO
#include <trexio.h>
#endif

#include <algorithm>
#include <cstring>

// Conversion constants
#define BOHR_TO_ANGSTROM 0.529177210903

namespace trexio_loader {

struct TrexioData {
    bool initialized = false;
    str_t current_file = {0};
    
    // Molecule data read from TREXIO file
    int64_t num_atoms = 0;
    int64_t num_electrons = 0;
    int64_t num_alpha_electrons = 0;
    int64_t num_beta_electrons = 0;
    
    // Basis set info
    int64_t num_shells = 0;
    int64_t num_primitives = 0;
    
    // MO info
    int64_t num_mo = 0;
    
    // Status info
    char info_text[1024] = "No TREXIO file loaded";
    char error_text[512] = "";
};

static TrexioData data = {};

#ifdef MD_TREXIO

// Helper function to convert TREXIO coordinates from Bohr to Angstrom
static void convert_coords_bohr_to_angstrom(double* coords, int64_t num_atoms) {
    for (int64_t i = 0; i < num_atoms * 3; ++i) {
        coords[i] *= BOHR_TO_ANGSTROM;
    }
}

static bool load_trexio_file(str_t filename) {
    MD_LOG_INFO("Attempting to load TREXIO file: %.*s", (int)filename.len, filename.ptr);
    
    // Clear previous state
    data.num_atoms = 0;
    data.num_electrons = 0;
    data.num_alpha_electrons = 0;
    data.num_beta_electrons = 0;
    data.num_shells = 0;
    data.num_primitives = 0;
    data.num_mo = 0;
    data.error_text[0] = '\0';
    
    // Open TREXIO file in read mode
    trexio_t* file = nullptr;
    trexio_exit_code rc;
    
    // Convert str_t to null-terminated string
    char* fname = (char*)malloc(filename.len + 1);
    if (!fname) {
        MD_LOG_ERROR("Failed to allocate memory for filename");
        snprintf(data.error_text, sizeof(data.error_text), "Memory allocation failed");
        return false;
    }
    memcpy(fname, filename.ptr, filename.len);
    fname[filename.len] = '\0';
    
    file = trexio_open(fname, 'r', TREXIO_AUTO, &rc);
    free(fname);
    
    if (file == nullptr || rc != TREXIO_SUCCESS) {
        snprintf(data.error_text, sizeof(data.error_text), 
                 "Failed to open: %s", trexio_string_of_error(rc));
        MD_LOG_ERROR("Failed to open TREXIO file: %s", trexio_string_of_error(rc));
        return false;
    }
    
    // Read nucleus count
    rc = trexio_read_nucleus_num(file, &data.num_atoms);
    if (rc != TREXIO_SUCCESS) {
        snprintf(data.error_text, sizeof(data.error_text), 
                 "Failed to read nucleus count: %s", trexio_string_of_error(rc));
        MD_LOG_ERROR("Failed to read nucleus count: %s", trexio_string_of_error(rc));
        trexio_close(file);
        return false;
    }
    
    // Read electron counts
    rc = trexio_read_electron_up_num(file, &data.num_alpha_electrons);
    if (rc == TREXIO_SUCCESS) {
        rc = trexio_read_electron_dn_num(file, &data.num_beta_electrons);
        if (rc == TREXIO_SUCCESS) {
            data.num_electrons = data.num_alpha_electrons + data.num_beta_electrons;
        }
    }
    
    if (rc != TREXIO_SUCCESS) {
        MD_LOG_WARN("Failed to read electron count: %s", trexio_string_of_error(rc));
        // Try alternative method
        rc = trexio_read_electron_num(file, &data.num_electrons);
        if (rc != TREXIO_SUCCESS) {
            MD_LOG_WARN("Could not read electron count via alternative method");
            data.num_electrons = 0;
        }
    }
    
    // Try to read basis set info (optional)
    rc = trexio_read_basis_shell_num(file, &data.num_shells);
    if (rc == TREXIO_SUCCESS) {
        MD_LOG_INFO("Found %lld basis shells", (long long)data.num_shells);
    }
    
    rc = trexio_read_basis_prim_num(file, &data.num_primitives);
    if (rc == TREXIO_SUCCESS) {
        MD_LOG_INFO("Found %lld basis primitives", (long long)data.num_primitives);
    }
    
    // Try to read MO info (optional)
    rc = trexio_read_mo_num(file, &data.num_mo);
    if (rc == TREXIO_SUCCESS) {
        MD_LOG_INFO("Found %lld molecular orbitals", (long long)data.num_mo);
    }
    
    MD_LOG_INFO("TREXIO file loaded successfully: %lld atoms, %lld electrons (α=%lld, β=%lld)", 
                (long long)data.num_atoms, (long long)data.num_electrons,
                (long long)data.num_alpha_electrons, (long long)data.num_beta_electrons);
    
    snprintf(data.info_text, sizeof(data.info_text), 
             "Loaded: %lld atoms, %lld electrons\nBasis: %lld shells, %lld primitives\nMOs: %lld orbitals", 
             (long long)data.num_atoms, (long long)data.num_electrons,
             (long long)data.num_shells, (long long)data.num_primitives,
             (long long)data.num_mo);
    
    trexio_close(file);
    data.initialized = true;
    data.current_file = filename;
    
    return true;
}

#else

static bool load_trexio_file(str_t filename) {
    (void)filename;
    MD_LOG_ERROR("TREXIO support not enabled at compile time");
    snprintf(data.error_text, sizeof(data.error_text), 
             "TREXIO support not compiled. Rebuild with -DVIAMD_ENABLE_TREXIO=ON");
    snprintf(data.info_text, sizeof(data.info_text), "TREXIO support not compiled");
    return false;
}

#endif // MD_TREXIO

} // namespace trexio_loader

// Component interface implementation
// This follows the pattern used by other VIAMD components

static struct TrexioComponent {
    
    void init() {
        MD_LOG_INFO("TREXIO component initialized");
#ifdef MD_TREXIO
        MD_LOG_INFO("TREXIO library support enabled (version: %s)", TREXIO_PACKAGE_VERSION);
#else
        MD_LOG_WARN("TREXIO library support not compiled in");
#endif
    }
    
    void draw_menu() {
        if (ImGui::BeginMenu("TREXIO")) {
            if (ImGui::MenuItem("Info")) {
                ImGui::OpenPopup("TREXIO Info");
            }
            
#ifdef MD_TREXIO
            if (ImGui::MenuItem("About TREXIO")) {
                ImGui::OpenPopup("About TREXIO");
            }
#endif
            ImGui::EndMenu();
        }
    }
    
    void draw_window() {
        // Info popup
        if (ImGui::BeginPopupModal("TREXIO Info", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("TREXIO File Loader Component");
            ImGui::Separator();
            
#ifdef MD_TREXIO
            ImGui::Text("Status: TREXIO support enabled (v%s)", TREXIO_PACKAGE_VERSION);
            ImGui::Spacing();
            
            if (trexio_loader::data.initialized) {
                ImGui::TextColored(ImVec4(0, 1, 0, 1), "File loaded successfully");
                ImGui::Text("File: %.*s", (int)trexio_loader::data.current_file.len, 
                           trexio_loader::data.current_file.ptr);
                ImGui::Separator();
                ImGui::Text("%s", trexio_loader::data.info_text);
                
                ImGui::Spacing();
                ImGui::Separator();
                ImGui::TextWrapped("Note: Full integration requires mdlib loader implementation. "
                                   "See docs/TREXIO_INTEGRATION.md for details.");
            } else {
                ImGui::TextColored(ImVec4(1, 0.5, 0, 1), "No file loaded");
                if (trexio_loader::data.error_text[0]) {
                    ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error: %s", trexio_loader::data.error_text);
                }
            }
#else
            ImGui::Text("Status: TREXIO support not compiled");
            ImGui::Spacing();
            ImGui::TextWrapped("To enable TREXIO support:");
            ImGui::BulletText("Install TREXIO library (https://github.com/TREX-CoE/trexio)");
            ImGui::BulletText("Rebuild with: cmake -DVIAMD_ENABLE_TREXIO=ON");
            ImGui::BulletText("See README.md for detailed instructions");
#endif
            
            ImGui::Spacing();
            if (ImGui::Button("Close")) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
        
#ifdef MD_TREXIO
        // About TREXIO popup
        if (ImGui::BeginPopupModal("About TREXIO", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("TREXIO File Format");
            ImGui::Separator();
            ImGui::Spacing();
            
            ImGui::TextWrapped(
                "TREXIO is an open-source file format and library for quantum chemistry data "
                "developed by the TREX Centre of Excellence.");
            
            ImGui::Spacing();
            ImGui::Text("Features:");
            ImGui::BulletText("Standardized format for quantum chemistry");
            ImGui::BulletText("Stores geometry, basis sets, and wavefunctions");
            ImGui::BulletText("Multiple backend support (HDF5, text)");
            ImGui::BulletText("Language bindings (C, Fortran, Python)");
            
            ImGui::Spacing();
            ImGui::Text("More info: https://github.com/TREX-CoE/trexio");
            
            ImGui::Spacing();
            if (ImGui::Button("Close")) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
#endif
    }
    
    void shutdown() {
        MD_LOG_INFO("TREXIO component shutdown");
    }
    
    // Event handlers
    void on_file_drop(const char* filename) {
        // Check if this is a TREXIO file
        const char* ext = strrchr(filename, '.');
        if (ext && (strcmp(ext, ".trexio") == 0 || strcmp(ext, ".h5") == 0)) {
            str_t file_str = {filename, strlen(filename)};
            if (trexio_loader::load_trexio_file(file_str)) {
                MD_LOG_INFO("TREXIO file loaded via drag-and-drop");
            }
        }
    }
    
} instance;

// Static registration happens via static initialization
// The component is automatically available when the module is linked

