#pragma once

#include <stdint.h>
#include <stddef.h>

#include <event.h>

// The Windows menu.
//
// A window says where it goes with a path, "Analysis/Conformation/Ramachandran": the last part is its label,
// the first is the section it is listed under (a header in the menu), and any parts in between are submenus.
// It is added while the menu is drawn, in answer to EventType_ViamdWindowMenuRequest:
//
//     case viamd::EventType_ViamdWindowMenuRequest:
//         viamd::window_menu_add(e, "Analysis/Contacts", &show_window);
//         break;
//
// Order among siblings is WindowMenuOpts::order, ascending. A window without one comes after those with one,
// by label, so a new window is listed without anyone deciding where; the order is polish. A submenu sits where
// its first window would: give the windows in it orders from the numbering of the level the submenu is on, not
// from 0. The sections come in the order of the table in window_menu.cpp. A section that is not in it is listed
// last and logged, as it is most likely a misspelling.
//
// A submenu whose windows are all disabled is shown disabled, with the first reason given as its tooltip, and
// one with no windows is not shown: which is what happens to the windows of a component that is not built.

namespace viamd {

struct WindowMenuOpts {
    int16_t     order           = INT16_MAX;   // among siblings, ascending. Unset: after those set, by label
    const char* shortcut        = nullptr;     // shown beside the label, e.g. ImGui::GetKeyChordName(chord) (imgui_internal.h)
    bool        enabled         = true;
    const char* disabled_reason = nullptr;     // the tooltip while disabled
};

// A window as window_menu_add stores it: the strings point into the request's text
struct WindowMenuItem {
    const char*    path = nullptr;
    bool*          open = nullptr;
    WindowMenuOpts opts;
};

// Filled while the menu is drawn: main.cpp adds its own windows, then broadcasts it to the components
struct WindowMenuRequest {
    WindowMenuItem items[64];
    size_t         count = 0;
    char           text[4096];
    size_t         text_len = 0;
};

// Add a window, in answer to EventType_ViamdWindowMenuRequest. Any other event is ignored. The strings are
// copied, so they only have to live for the call (ImGui::GetKeyChordName's does not live longer).
void window_menu_add(const Event& e, const char* path, bool* open, const WindowMenuOpts& opts = {});
void window_menu_add(WindowMenuRequest* req, const char* path, bool* open, const WindowMenuOpts& opts = {});

// The menu as a tree: what window_menu_draw draws, built apart from it so that it can be tested without ImGui.
// Node 0 is the root. Its children are the sections, in the order of the table, then any others in the order
// they were first named. The children of every other node are sorted by (order, label).
struct WindowMenuNode {
    const char* name         = nullptr;    // one part of a path, not NUL-terminated
    uint32_t    name_len     = 0;
    int32_t     order        = INT16_MAX;  // a window's own; a submenu's is that of its first window
    int32_t     item         = -1;         // a window: index into the request's items. -1: a section or submenu
    int32_t     first_child  = -1;
    int32_t     next_sibling = -1;
    uint32_t    num_windows  = 0;          // windows at or under the node
    uint32_t    num_enabled  = 0;          // of those, enabled
    uint32_t    num_open     = 0;          // of those, enabled and open
};

struct WindowMenuTree {
    WindowMenuNode nodes[192];
    size_t         count = 0;
};

void window_menu_build(WindowMenuTree* tree, const WindowMenuRequest* req);

// The content of the menu, between ImGui::BeginMenu("Windows") and ImGui::EndMenu()
void window_menu_draw(const WindowMenuRequest* req);

}
