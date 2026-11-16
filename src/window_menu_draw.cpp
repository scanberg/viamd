#include <window_menu.h>

#include <imgui.h>

#include <stdio.h>

// The ImGui half of window_menu: the tree itself is built in window_menu.cpp, which has no ImGui in it and is
// tested on its own (test/test_window_menu.cpp)

namespace viamd {

static void node_label(char* buf, size_t cap, const WindowMenuNode& node) {
    snprintf(buf, cap, "%.*s", (int)node.name_len, node.name);
}

static const char* first_disabled_reason(const WindowMenuTree& tree, const WindowMenuRequest& req, int32_t idx) {
    const WindowMenuNode& node = tree.nodes[idx];
    if (node.item != -1) {
        const WindowMenuItem& item = req.items[node.item];
        return item.opts.enabled ? nullptr : item.opts.disabled_reason;
    }
    for (int32_t c = node.first_child; c != -1; c = tree.nodes[c].next_sibling) {
        if (const char* reason = first_disabled_reason(tree, req, c)) return reason;
    }
    return nullptr;
}

static void draw_children(const WindowMenuTree& tree, const WindowMenuRequest& req, int32_t parent, int depth) {
    bool separated = false;
    for (int32_t c = tree.nodes[parent].first_child; c != -1; c = tree.nodes[c].next_sibling) {
        const WindowMenuNode& node = tree.nodes[c];
        if (node.num_windows == 0) continue;

        char label[128];
        node_label(label, sizeof(label), node);
        ImGui::PushID(node.name, node.name + node.name_len);

        if (node.item != -1) {
            // A window without a section is listed last, apart
            if (depth == 0 && !separated) {
                ImGui::Separator();
                separated = true;
            }
            const WindowMenuItem& item = req.items[node.item];
            ImGui::MenuItem(label, item.opts.shortcut, item.open, item.opts.enabled);
            if (!item.opts.enabled && item.opts.disabled_reason) {
                ImGui::SetItemTooltip("%s", item.opts.disabled_reason);
            }
        } else if (depth == 0) {
            ImGui::SeparatorText(label);
            draw_children(tree, req, c, depth + 1);
        } else {
            // How many of its windows are open is shown on a submenu, so that none is out of sight. It is left out
            // of the ID, or the submenu would close as a window in it is opened or closed.
            char menu_label[160];
            if (node.num_open > 0) {
                snprintf(menu_label, sizeof(menu_label), "%s (%u)###%s", label, node.num_open, label);
            } else {
                snprintf(menu_label, sizeof(menu_label), "%s###%s", label, label);
            }
            const bool enabled = node.num_enabled > 0;
            if (ImGui::BeginMenu(menu_label, enabled)) {
                draw_children(tree, req, c, depth + 1);
                ImGui::EndMenu();
            }
            if (!enabled) {
                if (const char* reason = first_disabled_reason(tree, req, c)) {
                    ImGui::SetItemTooltip("%s", reason);
                }
            }
        }

        ImGui::PopID();
    }
}

void window_menu_draw(const WindowMenuRequest* req) {
    static WindowMenuTree tree;  // Not on the stack: it is several KB, and only the frame's menu uses it
    window_menu_build(&tree, req);

    // As with the checkboxes the menu had: a click opens or closes a window and leaves the menu open, so that
    // several can be opened in one go
    ImGui::PushItemFlag(ImGuiItemFlags_AutoClosePopups, false);
    draw_children(tree, *req, 0, 0);
    ImGui::PopItemFlag();
}

}
