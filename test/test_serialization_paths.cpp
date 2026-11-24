#include "utest.h"

#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>

#include <string.h>

// The paths of the workspace format: a struct written by one function under a prefix, and read by another from the
// entries of its section, by path. A toy object stands in for a representation: fields of its own, a block, an array
// of blocks, lists, enums by name, a field that sets others up when it is read, and the keys of an older workspace.

namespace {

enum class Shade { Light, Dark, Count };
const char* const shade_names[] = { "Light", "Dark" };

struct Scale {
    int   colormap = 3;
    float range[2] = { 0.0f, 1.0f };
    bool  symmetric = false;
};

struct Item {
    md_float3_t pos = {};
    bool  on = true;
    float opacity = 1.0f;
};

struct Thing {
    char   name[32] = "thing";
    int    kind = 0;
    Shade  shade = Shade::Light;
    double iso = 0.05;
    int    set_up = 0;          // what kind sets up when it is read

    bool   has_scale = true;
    Scale  scale;

    float  values[8] = {};
    int    num_values = 0;

    float  curve[2 * 4] = { 0, 0, 1, 1 };
    int    num_curve = 2;       // not empty to begin with

    Item   items[4];
    int    num_items = 0;
};

// Colour maps by a name of their own, as ImPlot's are
const char* cmap_name(int i) {
    static const char* names[] = { "Deep", "Dark", "Pastel", "Paired", "Viridis" };
    return 0 <= i && i < (int)MD_COUNTOF(names) ? names[i] : nullptr;
}
int cmap_index(md_str_t name) {
    for (int i = 0; cmap_name(i); ++i) {
        if (md_str_eq_cstr(name, cmap_name(i))) return i;
    }
    return -1;
}

void write(viamd::serialization_state_t& s, viamd::PathPrefix p, const Scale& sc) {
    viamd::write_named(s, p("Colormap"), sc.colormap, cmap_name);
    viamd::write_flt_vec(s, p("Range"), sc.range, 2);
    viamd::write_bool(s, p("Symmetric"), sc.symmetric);
}

void read(viamd::Entries& e, viamd::PathPrefix p, Scale& sc) {
    viamd::get_named(e, p("Colormap"), sc.colormap, cmap_index);
    viamd::get(e, p("Range"), sc.range, 2);
    viamd::get(e, p("Symmetric"), sc.symmetric);
}

void write(viamd::serialization_state_t& s, const Thing& t) {
    viamd::PathPrefix p;
    viamd::write_str(s, p("Name"), md_str_from_cstr(t.name));
    viamd::write_int(s, p("Kind"), t.kind);
    viamd::write_enum(s, p("Shade"), (int)t.shade, viamd::enum_names(shade_names));
    viamd::write_dbl(s, p("Iso"), t.iso);
    viamd::write_int(s, p("SetUp"), t.set_up);
    if (t.has_scale) write(s, p.sub("Scale"), t.scale);
    viamd::write_list(s, p("Values"), t.values, t.num_values);
    viamd::write_list(s, p("Curve"), t.curve, t.num_curve, 2);
    viamd::write_int(s, p("Item"), t.num_items);
    for (int i = 0; i < t.num_items; ++i) {
        viamd::PathPrefix q = p.sub("Item", i);
        viamd::write_vec3(s, q("Pos"), t.items[i].pos);
        viamd::write_bool(s, q("On"), t.items[i].on);
        viamd::write_flt(s, q("Opacity"), t.items[i].opacity);
    }
}

void read(viamd::Entries& e, Thing& t) {
    viamd::PathPrefix p;
    viamd::get(e, p("Name"), t.name, sizeof(t.name));
    // Set up by the kind first: what the file has besides, it keeps, wherever it is in the file
    if (viamd::get(e, p("Kind"), t.kind)) {
        t.set_up = 100 + t.kind;
        t.shade = Shade::Dark;
    }
    viamd::get_enum(e, p("Shade"), t.shade, viamd::enum_names(shade_names));
    viamd::get(e, p("Iso"), t.iso);
    viamd::get(e, p("SetUp"), t.set_up);
    read(e, p.sub("Scale"), t.scale);
    viamd::get_list(e, p("Values"), t.values, t.num_values, 8);
    viamd::get_list(e, p("Curve"), t.curve, t.num_curve, 4, 2);
    viamd::get_count(e, p("Item"), t.num_items, 4);
    for (int i = 0; i < t.num_items; ++i) {
        viamd::PathPrefix q = p.sub("Item", i);
        viamd::get(e, q("Pos"), t.items[i].pos);
        viamd::get(e, q("On"), t.items[i].on);
        if (viamd::get(e, q("Opacity"), t.items[i].opacity)) t.items[i].opacity = MD_CLAMP(t.items[i].opacity, 0.0f, 1.0f);
    }
}

struct text_t {
    viamd::serialization_state_t state = {};
    text_t() {
        state.filename = MD_STR_LIT("test.via");
        md_str_writer_init(&state.sw, md_get_heap_allocator());
    }
    ~text_t() { md_str_writer_free(&state.sw); }
    md_str_t str() { return md_str_writer_to_str(&state.sw); }
    bool has_line(const char* line) {
        const md_str_t s = str();
        const size_t n = strlen(line);
        for (const char* p = s.ptr; p + n <= s.ptr + s.len; ++p) {
            if (memcmp(p, line, n) == 0 && (p == s.ptr || p[-1] == '\n') && (p + n == s.ptr + s.len || p[n] == '\n')) return true;
        }
        return false;
    }
};

viamd::deserialization_state_t reader_of(md_str_t text) {
    viamd::deserialization_state_t d = {};
    d.filename = MD_STR_LIT("test.via");
    d.text = text;
    return d;
}

// The entries of the one section of the text
viamd::Entries entries_of(md_str_t text, viamd::deserialization_state_t& d, md_allocator_i* alloc, const viamd::PathAlias* aliases = nullptr, size_t num_aliases = 0) {
    d = reader_of(text);
    md_str_t section;
    viamd::next_section_header(section, d);
    return viamd::read_entries(d, aliases, num_aliases, alloc);
}

Thing a_thing() {
    Thing t;
    snprintf(t.name, sizeof(t.name), "orbital 12");
    t.kind = 2;
    t.shade = Shade::Light;          // after the kind, which set it up as Dark
    t.iso = 0.125;
    t.set_up = 7;
    t.scale.colormap = 4;
    t.scale.range[0] = -0.5f;
    t.scale.range[1] = 0.5f;
    t.scale.symmetric = true;
    t.num_values = 3;
    t.values[0] = -0.25f; t.values[1] = 0.0f; t.values[2] = 0.75f;
    t.num_curve = 3;
    const float curve[] = { 0.0f, 0.0f, 0.5f, 0.25f, 1.0f, 1.0f };
    memcpy(t.curve, curve, sizeof(curve));
    t.num_items = 2;
    t.items[0].pos = md_float3_t{1, 2, 3};
    t.items[0].on = false;
    t.items[0].opacity = 0.5f;
    t.items[1].pos = md_float3_t{-1, 0, 4};
    t.items[1].opacity = 0.25f;
    return t;
}

struct arena_t {
    md_allocator_i* a = md_arena_allocator_create(md_get_heap_allocator(), MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    ~arena_t() { md_arena_allocator_destroy(a); }
};

}  // namespace

UTEST(viamd_serialization_paths, round_trip) {
    const Thing src = a_thing();
    text_t txt;
    viamd::write_section_header(txt.state, MD_STR_LIT("Thing"));
    write(txt.state, src);

    EXPECT_TRUE(txt.has_line("Scale.Colormap=Viridis"));
    EXPECT_TRUE(txt.has_line("Shade=Light"));
    EXPECT_TRUE(txt.has_line("Item=2"));
    EXPECT_TRUE(txt.has_line("Item.1.Opacity=0.25"));
    EXPECT_TRUE(txt.has_line("Curve=0,0,0.5,0.25,1,1"));

    arena_t arena;
    viamd::deserialization_state_t d;
    viamd::Entries e = entries_of(txt.str(), d, arena.a);
    Thing dst;
    read(e, dst);
    for (size_t i = 0; i < e.count; ++i) EXPECT_TRUE(e.items[i].read);

    EXPECT_STREQ(src.name, dst.name);
    EXPECT_EQ(src.kind, dst.kind);
    EXPECT_EQ((int)src.shade, (int)dst.shade);
    EXPECT_EQ(src.iso, dst.iso);
    EXPECT_EQ(src.set_up, dst.set_up);
    EXPECT_EQ(src.scale.colormap, dst.scale.colormap);
    EXPECT_EQ(src.scale.range[0], dst.scale.range[0]);
    EXPECT_EQ(src.scale.symmetric, dst.scale.symmetric);
    EXPECT_EQ(src.num_values, dst.num_values);
    EXPECT_EQ(src.values[2], dst.values[2]);
    EXPECT_EQ(src.num_curve, dst.num_curve);
    EXPECT_EQ(src.curve[3], dst.curve[3]);
    EXPECT_EQ(src.num_items, dst.num_items);
    EXPECT_EQ(src.items[1].pos.z, dst.items[1].pos.z);
    EXPECT_EQ(src.items[0].on, dst.items[0].on);
    EXPECT_EQ(src.items[0].opacity, dst.items[0].opacity);
}

UTEST(viamd_serialization_paths, read_in_the_order_of_the_code) {
    // The kind sets up the shade: it is read first whatever the file's order, and the shade the file has stands
    arena_t arena;
    viamd::deserialization_state_t d;
    viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nShade=Light\nSetUp=5\nKind=1\n"), d, arena.a);
    Thing t;
    read(e, t);
    EXPECT_EQ(1, t.kind);
    EXPECT_EQ((int)Shade::Light, (int)t.shade);
    EXPECT_EQ(5, t.set_up);
}

UTEST(viamd_serialization_paths, enums_by_name_and_by_index) {
    arena_t arena;
    viamd::deserialization_state_t d;
    {
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nShade=Dark\n"), d, arena.a);
        Thing t;
        read(e, t);
        EXPECT_EQ((int)Shade::Dark, (int)t.shade);
    }
    {
        // An older workspace's index
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nShade=1\n"), d, arena.a);
        Thing t;
        read(e, t);
        EXPECT_EQ((int)Shade::Dark, (int)t.shade);
    }
    {
        // Neither a name nor an index: left as it was
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nShade=Dusk\n"), d, arena.a);
        Thing t;
        read(e, t);
        EXPECT_EQ((int)Shade::Light, (int)t.shade);
        viamd::Entries f = entries_of(MD_STR_LIT("[Thing]\nShade=2\n"), d, arena.a);
        read(f, t);
        EXPECT_EQ((int)Shade::Light, (int)t.shade);
    }
    {
        // A colour map by its name, or its number
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nScale.Colormap=Pastel\n"), d, arena.a);
        Thing t;
        read(e, t);
        EXPECT_EQ(2, t.scale.colormap);
        viamd::Entries f = entries_of(MD_STR_LIT("[Thing]\nScale.Colormap=1\n"), d, arena.a);
        read(f, t);
        EXPECT_EQ(1, t.scale.colormap);
    }
}

UTEST(viamd_serialization_paths, paths_are_matched_whole) {
    arena_t arena;
    viamd::deserialization_state_t d;
    viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nIsoX=1\nScale.RangeX=1,2\nItem.0.PosX=1,2,3\nScaleRange=3,4\n"), d, arena.a);
    Thing t;
    read(e, t);
    EXPECT_EQ(0.05, t.iso);
    EXPECT_EQ(0.0f, t.scale.range[0]);
    // Item.0.PosX is an item's entry: the item is there, with nothing of it read
    EXPECT_EQ(1, t.num_items);
    EXPECT_EQ(0.0f, t.items[0].pos.x);
    int unread = 0;
    for (size_t i = 0; i < e.count; ++i) unread += !e.items[i].read;
    EXPECT_EQ(4, unread);
}

UTEST(viamd_serialization_paths, arrays) {
    arena_t arena;
    viamd::deserialization_state_t d;
    {
        // The count says what there is, also none
        Thing t;
        t.num_items = 3;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nItem=0\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(0, t.num_items);
    }
    {
        // Items past the count make it larger; past the capacity they are not read
        Thing t;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nItem=1\nItem.2.Opacity=0.5\nItem.7.Opacity=0.5\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(3, t.num_items);
        EXPECT_EQ(0.5f, t.items[2].opacity);
    }
    {
        // Clamped as it is read
        Thing t;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nItem.0.Opacity=4\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(1, t.num_items);
        EXPECT_EQ(1.0f, t.items[0].opacity);
    }
    {
        // No count and no items: as it was
        Thing t;
        t.num_items = 2;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nIso=1\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(2, t.num_items);
    }
}

UTEST(viamd_serialization_paths, lists) {
    arena_t arena;
    viamd::deserialization_state_t d;
    {
        // An empty list is no items
        Thing t;
        text_t txt;
        viamd::write_section_header(txt.state, MD_STR_LIT("Thing"));
        write(txt.state, t);
        EXPECT_TRUE(txt.has_line("Values="));
        t.num_values = 5;
        viamd::Entries e = entries_of(txt.str(), d, arena.a);
        read(e, t);
        EXPECT_EQ(0, t.num_values);
    }
    {
        // Pairs: an odd count does not read
        Thing t;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nCurve=0,0,1\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(2, t.num_curve);
    }
    {
        // More than fit does not read
        Thing t;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nValues=1,2,3,4,5,6,7,8,9\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(0, t.num_values);
    }
    {
        // An older workspace's items, a key each
        Thing t;
        viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nValues.0=0.5\nValues.2=1.5\n"), d, arena.a);
        read(e, t);
        EXPECT_EQ(3, t.num_values);
        EXPECT_EQ(1.5f, t.values[2]);
        viamd::Entries f = entries_of(MD_STR_LIT("[Thing]\nCurve.2=0.5,0.75\n"), d, arena.a);
        read(f, t);
        EXPECT_EQ(3, t.num_curve);
        EXPECT_EQ(0.75f, t.curve[5]);
    }
}

UTEST(viamd_serialization_paths, aliases_of_older_keys) {
    static const viamd::PathAlias aliases[] = {
        { "ThingColormap",       "Scale.Colormap" },
        { "ThingItem#",          "Item.#.Pos" },
        { "ThingItemOpacity#",   "Item.#.Opacity" },
        { "ThingValue#",         "Values.#" },
        { "ThingIso",            "Iso" },
    };
    char buf[64];
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingColormap"), aliases, MD_COUNTOF(aliases)), MD_STR_LIT("Scale.Colormap")));
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItem12"), aliases, MD_COUNTOF(aliases)), MD_STR_LIT("Item.12.Pos")));
    // The longer family is not taken for the shorter one, whichever comes first
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItemOpacity3"), aliases, MD_COUNTOF(aliases)), MD_STR_LIT("Item.3.Opacity")));
    EXPECT_TRUE(md_str_empty(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItem"), aliases, MD_COUNTOF(aliases))));
    EXPECT_TRUE(md_str_empty(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItemX3"), aliases, MD_COUNTOF(aliases))));
    // What does not fit is none
    char small[8];
    EXPECT_TRUE(md_str_empty(viamd::path_alias(small, sizeof(small), MD_STR_LIT("ThingItemOpacity3"), aliases, MD_COUNTOF(aliases))));

    // Read as the paths they now are; a newer key at the same path stands over an older one before it
    arena_t arena;
    viamd::deserialization_state_t d;
    viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nThingColormap=Paired\nThingItem1=1,2,3\nThingItemOpacity1=0.5\nThingValue1=0.5\nThingIso=0.25\nIso=0.5\nThingOld=1\n"),
                                  d, arena.a, aliases, MD_COUNTOF(aliases));
    Thing t;
    read(e, t);
    EXPECT_EQ(3, t.scale.colormap);
    EXPECT_EQ(2, t.num_items);
    EXPECT_EQ(3.0f, t.items[1].pos.z);
    EXPECT_EQ(0.5f, t.items[1].opacity);
    EXPECT_EQ(2, t.num_values);
    EXPECT_EQ(0.5f, t.values[1]);
    EXPECT_EQ(0.5, t.iso);
    // A key whose meaning changed is read by its old name, by whoever reads the struct; else it is unread
    EXPECT_TRUE(viamd::get_value(e, MD_STR_LIT("ThingOld")) != nullptr);
    for (size_t i = 0; i < e.count; ++i) EXPECT_TRUE(e.items[i].read);
}

UTEST(viamd_serialization_paths, values_that_do_not_read) {
    arena_t arena;
    viamd::deserialization_state_t d;
    viamd::Entries e = entries_of(MD_STR_LIT("[Thing]\nIso=abc\nScale.Range=1;2\nItem.0.Pos=1,2\n"), d, arena.a);
    Thing t;
    read(e, t);
    EXPECT_EQ(0.05, t.iso);
    EXPECT_EQ(1.0f, t.scale.range[1]);
    EXPECT_EQ(0.0f, t.items[0].pos.x);
    // Read, but left out: none of them is reported as read by nothing
    for (size_t i = 0; i < e.count; ++i) EXPECT_TRUE(e.items[i].read);
}

UTEST(viamd_serialization_paths, extract_flt_list) {
    float v[4];
    size_t n = 99;
    EXPECT_TRUE(viamd::extract_flt_list(v, &n, 4, MD_STR_LIT(" 1, -2.5 ,3e-2")));
    EXPECT_EQ(3u, (uint32_t)n);
    EXPECT_EQ(-2.5f, v[1]);
    EXPECT_TRUE(viamd::extract_flt_list(v, &n, 4, MD_STR_LIT("")));
    EXPECT_EQ(0u, (uint32_t)n);
    n = 7;
    EXPECT_FALSE(viamd::extract_flt_list(v, &n, 2, MD_STR_LIT("1,2,3")));
    EXPECT_FALSE(viamd::extract_flt_list(v, &n, 4, MD_STR_LIT("1,,2")));
    EXPECT_EQ(7u, (uint32_t)n);
}
