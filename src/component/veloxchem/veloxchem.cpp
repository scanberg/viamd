#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>

#include <core/md_log.h>

#include <imgui_widgets.h>
#include <implot_widgets.h>

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }

    bool show_window = false;

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                break;
            case viamd::EventType_ViamdShutdown:
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("VeloxChem", &show_window);
                break;
            case viamd::EventType_ViamdDeserialize:
                break;
            case viamd::EventType_ViamdSerialize:
                break;
            case viamd::EventType_ViamdFileOpen:

                break;
            default:
                break;
            }
        }
    }

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("VeloxChem", &show_window)) {

        }
        ImGui::End();
    }
} instance;
