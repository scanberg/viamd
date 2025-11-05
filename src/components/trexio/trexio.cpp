// TREXIO Loader Component for VIAMD
// This component provides support for loading TREXIO files into VIAMD

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

// TREXIO component implementation
// This is a minimal stub that registers the component but does not yet implement full functionality

namespace trexio_loader {

struct TrexioData {
    bool initialized = false;
    str_t current_file = {0};
    
    // Molecule data
    int64_t num_atoms = 0;
    int64_t num_electrons = 0;
    
    // Basic info
    char info_text[1024] = "No TREXIO file loaded";
};

static TrexioData data = {};

#ifdef MD_TREXIO
static bool load_trexio_file(str_t filename) {
    MD_LOG_INFO("Attempting to load TREXIO file: %.*s", (int)filename.len, filename.ptr);
    
    // Open TREXIO file in read mode
    trexio_t* file = nullptr;
    trexio_exit_code rc;
    
    // Convert str_t to null-terminated string
    char* fname = (char*)malloc(filename.len + 1);
    if (!fname) {
        MD_LOG_ERROR("Failed to allocate memory for filename");
        return false;
    }
    memcpy(fname, filename.ptr, filename.len);
    fname[filename.len] = '\0';
    
    file = trexio_open(fname, 'r', TREXIO_AUTO, &rc);
    free(fname);
    
    if (file == nullptr || rc != TREXIO_SUCCESS) {
        MD_LOG_ERROR("Failed to open TREXIO file: %s", trexio_string_of_error(rc));
        return false;
    }
    
    // Read nucleus count
    rc = trexio_read_nucleus_num(file, &data.num_atoms);
    if (rc != TREXIO_SUCCESS) {
        MD_LOG_ERROR("Failed to read nucleus count: %s", trexio_string_of_error(rc));
        trexio_close(file);
        return false;
    }
    
    // Read electron count
    rc = trexio_read_electron_num(file, &data.num_electrons);
    if (rc != TREXIO_SUCCESS) {
        MD_LOG_WARN("Failed to read electron count: %s", trexio_string_of_error(rc));
        data.num_electrons = 0; // Not critical, continue
    }
    
    MD_LOG_INFO("TREXIO file loaded successfully: %lld atoms, %lld electrons", 
                (long long)data.num_atoms, (long long)data.num_electrons);
    
    snprintf(data.info_text, sizeof(data.info_text), 
             "Loaded: %lld atoms, %lld electrons", 
             (long long)data.num_atoms, (long long)data.num_electrons);
    
    trexio_close(file);
    data.initialized = true;
    return true;
}
#else
static bool load_trexio_file(str_t filename) {
    (void)filename;
    MD_LOG_ERROR("TREXIO support not enabled at compile time");
    snprintf(data.info_text, sizeof(data.info_text), "TREXIO support not compiled");
    return false;
}
#endif

} // namespace trexio_loader

// Component interface implementation
// This follows the pattern used by other VIAMD components

static struct TrexioComponent {
    
    void init() {
        MD_LOG_INFO("TREXIO component initialized");
#ifdef MD_TREXIO
        MD_LOG_INFO("TREXIO library support enabled");
#else
        MD_LOG_WARN("TREXIO library support not compiled in");
#endif
    }
    
    void draw_menu() {
        if (ImGui::BeginMenu("TREXIO")) {
            if (ImGui::MenuItem("Info")) {
                ImGui::OpenPopup("TREXIO Info");
            }
            ImGui::EndMenu();
        }
    }
    
    void draw_window() {
        if (ImGui::BeginPopupModal("TREXIO Info", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("TREXIO File Loader Component");
            ImGui::Separator();
            
#ifdef MD_TREXIO
            ImGui::Text("Status: TREXIO support enabled");
            ImGui::Text("%s", trexio_loader::data.info_text);
            
            if (trexio_loader::data.initialized) {
                ImGui::Separator();
                ImGui::Text("Atoms: %lld", (long long)trexio_loader::data.num_atoms);
                ImGui::Text("Electrons: %lld", (long long)trexio_loader::data.num_electrons);
            }
#else
            ImGui::Text("Status: TREXIO support not compiled");
            ImGui::TextWrapped("To enable TREXIO support, rebuild with -DVIAMD_ENABLE_TREXIO=ON");
#endif
            
            if (ImGui::Button("Close")) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
    }
    
    void shutdown() {
        MD_LOG_INFO("TREXIO component shutdown");
    }
    
    // Event handlers
    void on_file_drop(const char* filename) {
        // Check if this is a TREXIO file
        const char* ext = strrchr(filename, '.');
        if (ext && strcmp(ext, ".trexio") == 0) {
            str_t file_str = {filename, strlen(filename)};
            trexio_loader::load_trexio_file(file_str);
        }
    }
    
} instance;

// Static registration happens via static initialization
// The component is automatically available when the module is linked
