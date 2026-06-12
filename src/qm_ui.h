#pragma once

#include <imgui.h>

#include <cstdio>

namespace qm_ui {

static constexpr const char* summary_label = "Summary";
static constexpr const char* orbital_grid_label = "Orbital Grid";

static inline bool begin_summary_window(const char* source, const char* unique_id, bool* show_window, ImGuiWindowFlags flags = ImGuiWindowFlags_NoFocusOnAppearing) {
    char label[128];
    const char* id = unique_id ? unique_id : (source ? source : "QM");
    std::snprintf(label, sizeof(label), "%s Summary###%s", source ? source : "QM", id);
    return ImGui::Begin(label, show_window, flags);
}

static inline bool begin_orbital_grid_window(const char* source, const char* unique_id, bool* show_window, ImGuiWindowFlags flags = ImGuiWindowFlags_NoFocusOnAppearing) {
    char label[128];
    const char* id = unique_id ? unique_id : (source ? source : "QM");
    std::snprintf(label, sizeof(label), "%s Orbital Grid###%s", source ? source : "QM", id);
    return ImGui::Begin(label, show_window, flags);
}

static inline void draw_source_label(const char* source) {
    ImGui::Text("Source:                %s", source ? source : "<unknown>");
}

} // namespace qm_ui
