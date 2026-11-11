// The display units' side in the application: the .ini and the Settings menu. Kept apart from display_units.cpp so
// that the conversions build and are tested without ImGui.

#include <display_units.h>
#include <app_settings.h>

#include <imgui.h>

namespace display_units {

static void adopt(void*) { adopt_ini(); }

void register_settings() {
    // The defaults have to stand before the .ini is read, since a file which says nothing about a
    // quantity leaves its string untouched.
    reset();
    for (int i = 0; i < Quantity_Count; ++i) {
        size_t cap = 0;
        char* buf = ini_buffer((Quantity)i, &cap);
        app_settings::bind(md_str_from_cstr(ini_key((Quantity)i)), buf, cap);
    }
    app_settings::on_apply(adopt, nullptr);
}

void draw_settings_menu_items() {
    for (int i = 0; i < Quantity_Count; ++i) {
        const Quantity q = (Quantity)i;
        const size_t cur = current_option(q);
        size_t cap = 0;
        const char* shown = cur < num_options(q) ? option_label(q, cur) : ini_buffer(q, &cap);
        if (ImGui::BeginCombo(quantity_label(q), shown)) {
            for (size_t j = 0; j < num_options(q); ++j) {
                if (ImGui::Selectable(option_label(q, j), j == cur)) {
                    select(q, j);
                    app_settings::mark_dirty();
                }
            }
            ImGui::EndCombo();
        }
    }
}

}  // namespace display_units
