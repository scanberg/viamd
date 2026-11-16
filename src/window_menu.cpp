#include <window_menu.h>
#include <viamd_event.h>

#include <core/md_common.h>
#include <core/md_log.h>

#include <string.h>
#include <stdio.h>
#include <algorithm>

namespace viamd {

// The sections of the menu, in order. A section is named by the first part of a window's path; one not listed
// here is listed after these and logged.
static const char* const sections[] = {
    "View",
    "Data",
    "Analysis",
    "Scripting",
    "Output",
    "Developer",
};

// The request is made every frame the menu is open, so whatever is wrong with it would be logged every frame.
// Each message is logged once per run instead.
static bool first_report(const char* key) {
    static char   seen[32][96];
    static size_t num_seen = 0;
    for (size_t i = 0; i < num_seen; ++i) {
        if (strncmp(seen[i], key, sizeof(seen[i]) - 1) == 0) return false;
    }
    if (num_seen < MD_COUNTOF(seen)) {
        snprintf(seen[num_seen++], sizeof(seen[0]), "%s", key);
    }
    return true;
}

static const char* request_copy(WindowMenuRequest* req, const char* str) {
    if (!str) return nullptr;
    const size_t len = strlen(str);
    if (req->text_len + len + 1 > sizeof(req->text)) return nullptr;
    char* dst = req->text + req->text_len;
    memcpy(dst, str, len + 1);
    req->text_len += len + 1;
    return dst;
}

void window_menu_add(WindowMenuRequest* req, const char* path, bool* open, const WindowMenuOpts& opts) {
    // Not ASSERTs: a release build compiles them out and would go on into the NULL
    if (!req) return;
    if (!path || !path[0] || !open) {
        if (first_report("window_menu_add: no path or no flag")) {
            MD_LOG_ERROR("Window menu: a window was added without a path or without its open flag, and is not listed");
        }
        return;
    }
    if (req->count == MD_COUNTOF(req->items)) {
        if (first_report("window_menu_add: full")) {
            MD_LOG_ERROR("Window menu: more than %zu windows, '%s' and those after it are not listed", MD_COUNTOF(req->items), path);
        }
        return;
    }

    const size_t text_len = req->text_len;
    WindowMenuItem item;
    item.open = open;
    item.opts = opts;
    item.path                 = request_copy(req, path);
    item.opts.shortcut        = request_copy(req, opts.shortcut);
    item.opts.disabled_reason = request_copy(req, opts.disabled_reason);
    if (!item.path || (opts.shortcut && !item.opts.shortcut) || (opts.disabled_reason && !item.opts.disabled_reason)) {
        req->text_len = text_len;
        if (first_report("window_menu_add: text full")) {
            MD_LOG_ERROR("Window menu: out of room for text, '%s' and those after it are not listed", path);
        }
        return;
    }
    req->items[req->count++] = item;
}

void window_menu_add(const Event& e, const char* path, bool* open, const WindowMenuOpts& opts) {
    if (e.type != EventType_ViamdWindowMenuRequest || e.payload_type != EventPayloadType_WindowMenuRequest) return;
    // The request is const in the event as any payload is; it is there to be filled
    window_menu_add((WindowMenuRequest*)e.payload, path, open, opts);
}

static int32_t node_new(WindowMenuTree* tree, const char* name, uint32_t name_len) {
    if (tree->count == MD_COUNTOF(tree->nodes)) return -1;
    const int32_t idx = (int32_t)tree->count++;
    tree->nodes[idx] = WindowMenuNode{};
    tree->nodes[idx].name     = name;
    tree->nodes[idx].name_len = name_len;
    return idx;
}

static void node_append_child(WindowMenuTree* tree, int32_t parent, int32_t child) {
    int32_t* link = &tree->nodes[parent].first_child;
    while (*link != -1) link = &tree->nodes[*link].next_sibling;
    *link = child;
}

static bool name_eq(const WindowMenuNode& node, const char* name, uint32_t len) {
    return node.name_len == len && memcmp(node.name, name, len) == 0;
}

// A section or submenu: a window of the same name is a different thing, and is not merged with
static int32_t node_find_group(const WindowMenuTree* tree, int32_t parent, const char* name, uint32_t len) {
    for (int32_t c = tree->nodes[parent].first_child; c != -1; c = tree->nodes[c].next_sibling) {
        if (tree->nodes[c].item == -1 && name_eq(tree->nodes[c], name, len)) return c;
    }
    return -1;
}

static void node_aggregate(WindowMenuTree* tree, int32_t idx) {
    WindowMenuNode& node = tree->nodes[idx];
    if (node.item != -1) return;
    for (int32_t c = node.first_child; c != -1; c = tree->nodes[c].next_sibling) {
        node_aggregate(tree, c);
        const WindowMenuNode& child = tree->nodes[c];
        node.order        = MD_MIN(node.order, child.order);
        node.num_windows += child.num_windows;
        node.num_enabled += child.num_enabled;
        node.num_open    += child.num_open;
    }
}

// By label, ignoring case; where the labels are equal, in the order they were added
static bool node_less(const WindowMenuTree* tree, int32_t a, int32_t b) {
    const WindowMenuNode& na = tree->nodes[a];
    const WindowMenuNode& nb = tree->nodes[b];
    if (na.order != nb.order) return na.order < nb.order;
    const uint32_t len = MD_MIN(na.name_len, nb.name_len);
    for (uint32_t i = 0; i < len; ++i) {
        char ca = na.name[i], cb = nb.name[i];
        if (ca >= 'A' && ca <= 'Z') ca += 'a' - 'A';
        if (cb >= 'A' && cb <= 'Z') cb += 'a' - 'A';
        if (ca != cb) return ca < cb;
    }
    if (na.name_len != nb.name_len) return na.name_len < nb.name_len;
    return a < b;
}

static void node_sort_children(WindowMenuTree* tree, int32_t idx) {
    int32_t children[sizeof(tree->nodes) / sizeof(tree->nodes[0])];
    size_t  count = 0;
    for (int32_t c = tree->nodes[idx].first_child; c != -1; c = tree->nodes[c].next_sibling) {
        children[count++] = c;
    }
    std::sort(children, children + count, [tree](int32_t a, int32_t b) { return node_less(tree, a, b); });
    int32_t* link = &tree->nodes[idx].first_child;
    for (size_t i = 0; i < count; ++i) {
        *link = children[i];
        link = &tree->nodes[children[i]].next_sibling;
    }
    *link = -1;
}

static void node_sort_subtree(WindowMenuTree* tree, int32_t idx) {
    node_sort_children(tree, idx);
    for (int32_t c = tree->nodes[idx].first_child; c != -1; c = tree->nodes[c].next_sibling) {
        if (tree->nodes[c].item == -1) node_sort_subtree(tree, c);
    }
}

void window_menu_build(WindowMenuTree* tree, const WindowMenuRequest* req) {
    tree->count = 0;
    node_new(tree, "", 0);  // the root

    for (const char* section : sections) {
        node_append_child(tree, 0, node_new(tree, section, (uint32_t)strlen(section)));
    }

    for (size_t i = 0; i < req->count; ++i) {
        const WindowMenuItem& item = req->items[i];

        // The parts of the path, with empty ones ("a//b", a leading or trailing '/') skipped
        const char* part[8];
        uint32_t    part_len[8];
        size_t      num_parts = 0;
        for (const char* c = item.path; *c; ) {
            const char* end = strchr(c, '/');
            if (!end) end = c + strlen(c);
            if (end > c) {
                if (num_parts == MD_COUNTOF(part)) break;
                part[num_parts] = c;
                part_len[num_parts] = (uint32_t)(end - c);
                num_parts += 1;
            }
            c = *end ? end + 1 : end;
        }
        if (num_parts == 0) continue;
        if (num_parts == 1) {
            if (first_report(item.path)) {
                MD_LOG_ERROR("Window menu: '%s' has no section and is listed last; give it a path like 'Analysis/%s'", item.path, item.path);
            }
        }

        int32_t parent = 0;
        for (size_t k = 0; k + 1 < num_parts && parent != -1; ++k) {
            int32_t group = node_find_group(tree, parent, part[k], part_len[k]);
            if (group == -1) {
                if (parent == 0 && first_report(item.path)) {
                    MD_LOG_ERROR("Window menu: '%.*s' (in '%s') is not one of the sections in window_menu.cpp and is listed last",
                        (int)part_len[k], part[k], item.path);
                }
                group = node_new(tree, part[k], part_len[k]);
                if (group != -1) node_append_child(tree, parent, group);
            }
            parent = group;
        }

        const uint32_t leaf_len = part_len[num_parts - 1];
        const char*    leaf_name = part[num_parts - 1];
        const int32_t  leaf = parent != -1 ? node_new(tree, leaf_name, leaf_len) : -1;
        if (leaf == -1) {
            if (first_report("window_menu_build: full")) {
                MD_LOG_ERROR("Window menu: out of nodes, '%s' is not listed", item.path);
            }
            continue;
        }
        for (int32_t c = tree->nodes[parent].first_child; c != -1; c = tree->nodes[c].next_sibling) {
            if (tree->nodes[c].item != -1 && name_eq(tree->nodes[c], leaf_name, leaf_len) && first_report(item.path)) {
                MD_LOG_ERROR("Window menu: '%s' is added more than once", item.path);
                break;
            }
        }

        WindowMenuNode& node = tree->nodes[leaf];
        node.item        = (int32_t)i;
        node.order       = item.opts.order;
        node.num_windows = 1;
        node.num_enabled = item.opts.enabled ? 1 : 0;
        node.num_open    = (item.opts.enabled && *item.open) ? 1 : 0;
        node_append_child(tree, parent, leaf);
    }

    node_aggregate(tree, 0);

    // The root keeps the order of the sections; everything under it is sorted
    for (int32_t c = tree->nodes[0].first_child; c != -1; c = tree->nodes[c].next_sibling) {
        if (tree->nodes[c].item == -1) node_sort_subtree(tree, c);
    }
}

}
