// TREXIO Component for VIAMD
// Provides visualization of quantum chemistry data from TREXIO files
// Integrates with the event system and TREXIO UI panel

#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>

#ifdef MD_TREXIO
#include <md_trexio.h>
#include <ui/trexio_panel.h>
#endif

#include <md_util.h>
#include <core/md_log.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>

#include <app/IconsFontAwesome6.h>

#define HARTREE_TO_EV 27.211386245988

#ifdef MD_TREXIO

struct TrexioComponent : viamd::EventHandler {
    bool initialized = false;
    bool show_tools_window = false;
    
    void init() {
        if (initialized) return;
        
        md_allocator_i* alloc = md_get_heap_allocator();
        trexio_ui::initialize();
        
        viamd::event_system_register_handler(*this);
        initialized = true;
        
        MD_LOG_INFO("TREXIO Component initialized with UI panel");
    }
    
    void shutdown() {
        if (!initialized) return;
        
        trexio_ui::shutdown();
        initialized = false;
    }
    
    void process_events(const viamd::Event* events, size_t num_events) override {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& event = events[i];
            
            switch (event.type) {
                case viamd::EventType_ViamdInitialize:
                    init();
                    break;
                    
                case viamd::EventType_ViamdShutdown:
                    shutdown();
                    break;
                    
                case viamd::EventType_ViamdFrameTick:
                    if (initialized) {
                        trexio_ui::update();
                    }
                    break;
                    
                case viamd::EventType_ViamdWindowDrawMenu:
                    if (initialized) {
                        draw_menu_item();
                    }
                    break;
            }
        }
    }
    
    void draw_menu_item() {
        if (ImGui::MenuItem(ICON_FA_FLASK " TREXIO Tools")) {
            show_tools_window = !show_tools_window;
            trexio_ui::set_main_window_visible(show_tools_window);
        }
    }
};

static TrexioComponent instance = {};

#else // MD_TREXIO

// Stub implementation when TREXIO is not enabled
struct TrexioComponent : viamd::EventHandler {
    void init() {
        MD_LOG_INFO("TREXIO Component: Not compiled in (rebuild with -DVIAMD_ENABLE_TREXIO=ON)");
    }
    
    void shutdown() {}
    
    void process_events(const viamd::Event* events, size_t num_events) override {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& event = events[i];
            
            if (event.type == viamd::EventType_ViamdInitialize) {
                init();
                viamd::event_system_register_handler(*this);
            }
        }
    }
};

static TrexioComponent instance = {};

#endif // MD_TREXIO

