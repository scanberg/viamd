#include "utest.h"

#include <window_menu.h>
#include <viamd_event.h>

#include <string>
#include <string.h>

/* The Windows menu is the tree window_menu_build makes of the paths the windows give. These pin down where a
 * window ends up: sections in the order of the table, siblings by (order, label), a submenu where its first
 * window would be, and what a submenu shows of the windows in it. */

namespace {

// The tree as text, without what is not drawn (a section or submenu with no windows):
// "View{Representations,Animation},Analysis{Timelines,Conformation{Ramachandran}}"
std::string dump(const viamd::WindowMenuTree& tree, int32_t idx = 0) {
    std::string out;
    for (int32_t c = tree.nodes[idx].first_child; c != -1; c = tree.nodes[c].next_sibling) {
        const viamd::WindowMenuNode& node = tree.nodes[c];
        if (node.num_windows == 0) continue;
        if (!out.empty()) out += ",";
        out.append(node.name, node.name_len);
        if (node.item == -1) out += "{" + dump(tree, c) + "}";
    }
    return out;
}

int32_t find(const viamd::WindowMenuTree& tree, const char* name) {
    for (size_t i = 0; i < tree.count; ++i) {
        const viamd::WindowMenuNode& node = tree.nodes[i];
        if (node.name_len == strlen(name) && memcmp(node.name, name, node.name_len) == 0) return (int32_t)i;
    }
    return -1;
}

struct Fixture {
    viamd::WindowMenuRequest req;
    viamd::WindowMenuTree    tree;
    bool flags[64] = {};
    size_t next = 0;

    void add(const char* path, viamd::WindowMenuOpts opts = {}) {
        viamd::window_menu_add(&req, path, &flags[next++], opts);
    }
    std::string text;  // kept here: utest keeps the pointer it is given past the statement that made it

    const char* build() {
        viamd::window_menu_build(&tree, &req);
        text = dump(tree);
        return text.c_str();
    }
};

}

UTEST(window_menu, sections_in_table_order) {
    Fixture f;
    f.add("Developer/Profiler");
    f.add("Analysis/Timelines");
    f.add("View/Animation");
    f.add("Data/System");
    EXPECT_STREQ("View{Animation},Data{System},Analysis{Timelines},Developer{Profiler}", f.build());
}

UTEST(window_menu, unknown_section_last) {
    Fixture f;
    f.add("Analysys/Contacts");  // a misspelling: listed, but after the known ones
    f.add("Output/Structure Export");
    f.add("View/Animation");
    EXPECT_STREQ("View{Animation},Output{Structure Export},Analysys{Contacts}", f.build());
}

UTEST(window_menu, siblings_by_order_then_label) {
    Fixture f;
    f.add("Analysis/Contacts",      {.order = 30});
    f.add("Analysis/zeta");                         // unordered: after the ordered ones, by label ignoring case
    f.add("Analysis/Timelines",     {.order = 10});
    f.add("Analysis/Alpha");
    f.add("Analysis/Distributions", {.order = 20});
    EXPECT_STREQ("Analysis{Timelines,Distributions,Contacts,Alpha,zeta}", f.build());
}

UTEST(window_menu, submenu_where_its_first_window_is) {
    Fixture f;
    f.add("Analysis/Timelines",                  {.order = 10});
    f.add("Analysis/Shapespace",                 {.order = 50});
    f.add("Analysis/Conformation/Ring Pucker",   {.order = 70});
    f.add("Analysis/Conformation/Ramachandran",  {.order = 60});
    f.add("Analysis/Quantum Chemistry/Summary",  {.order = 100});
    f.add("Analysis/Contacts",                   {.order = 30});
    EXPECT_STREQ("Analysis{Timelines,Contacts,Shapespace,Conformation{Ramachandran,Ring Pucker},Quantum Chemistry{Summary}}",
        f.build());
    EXPECT_EQ(60, f.tree.nodes[find(f.tree, "Conformation")].order);
}

UTEST(window_menu, empty_parts_are_skipped) {
    Fixture f;
    f.add("/Analysis//Contacts/");
    EXPECT_STREQ("Analysis{Contacts}", f.build());
}

UTEST(window_menu, window_and_submenu_of_the_same_name_are_apart) {
    Fixture f;
    f.add("Analysis/Export", {.order = 10});
    f.add("Analysis/Export/Settings", {.order = 20});
    EXPECT_STREQ("Analysis{Export,Export{Settings}}", f.build());
}

UTEST(window_menu, counts_of_a_submenu) {
    Fixture f;
    f.add("Analysis/Quantum Chemistry/Summary", {.order = 10, .enabled = false, .disabled_reason = "no data"});
    f.add("Analysis/Quantum Chemistry/Export",  {.order = 20, .enabled = false});
    f.flags[0] = true;  // open, but disabled: not drawn, so not counted as open
    f.build();

    const viamd::WindowMenuNode& qm = f.tree.nodes[find(f.tree, "Quantum Chemistry")];
    EXPECT_EQ(2u, qm.num_windows);
    EXPECT_EQ(0u, qm.num_enabled);  // drawn disabled
    EXPECT_EQ(0u, qm.num_open);

    Fixture g;
    g.add("Analysis/Conformation/Ramachandran");
    g.add("Analysis/Conformation/Ring Pucker");
    g.flags[1] = true;
    g.build();
    const viamd::WindowMenuNode& conf = g.tree.nodes[find(g.tree, "Conformation")];
    EXPECT_EQ(2u, conf.num_enabled);
    EXPECT_EQ(1u, conf.num_open);
}

UTEST(window_menu, empty_sections_have_no_windows) {
    Fixture f;
    f.add("View/Animation");
    f.build();
    EXPECT_EQ(0u, f.tree.nodes[find(f.tree, "Analysis")].num_windows);
}

UTEST(window_menu, strings_are_copied) {
    Fixture f;
    char path[64] = "Analysis/Contacts";
    char shortcut[16] = "Ctrl+K";
    f.add(path, {.shortcut = shortcut});
    strcpy(path, "Overwritten/Path");
    strcpy(shortcut, "xxxxxx");
    EXPECT_STREQ("Analysis{Contacts}", f.build());
    EXPECT_STREQ("Ctrl+K", f.req.items[0].opts.shortcut);
}

UTEST(window_menu, add_through_the_event) {
    viamd::WindowMenuRequest req;
    bool open = false;

    viamd::Event e = {};
    e.type         = viamd::EventType_ViamdWindowMenuRequest;
    e.payload_type = viamd::EventPayloadType_WindowMenuRequest;
    e.payload      = &req;
    viamd::window_menu_add(e, "Analysis/Contacts", &open);
    EXPECT_EQ(1u, req.count);

    // Any other event is ignored, whatever its payload
    e.type = viamd::EventType_ViamdFrameTick;
    viamd::window_menu_add(e, "Analysis/Contacts", &open);
    EXPECT_EQ(1u, req.count);
}

UTEST(window_menu, no_path_or_flag_is_not_listed) {
    viamd::WindowMenuRequest req;
    bool open = false;
    viamd::window_menu_add(&req, nullptr, &open);
    viamd::window_menu_add(&req, "", &open);
    viamd::window_menu_add(&req, "Analysis/Contacts", nullptr);
    EXPECT_EQ(0u, req.count);
}

UTEST(window_menu, full_request_drops_the_rest) {
    viamd::WindowMenuRequest req;
    bool open = false;
    char path[32];
    for (int i = 0; i < 80; ++i) {
        snprintf(path, sizeof(path), "Analysis/W%d", i);
        viamd::window_menu_add(&req, path, &open);
    }
    EXPECT_EQ(sizeof(req.items) / sizeof(req.items[0]), req.count);
}
