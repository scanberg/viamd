#include "utest.h"

#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_writer.h>

#include <string.h>

// The paths of the workspace format and the Visitor that writes and reads them: one list of fields per struct, written
// in its order and read back in the order of the file. A toy object stands in for a representation: fields of its
// own, blocks, a block in a section of its own, an array of blocks, lists, enums by name, and a field that sets
// others up when it is read.

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

    md_float4_t color = { 1, 1, 1, 1 };   // in a section of its own
    int    detail = 2;
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

void visit(viamd::Visitor& v, Scale& s) {
    v.field_named("Colormap", s.colormap, cmap_name, cmap_index);
    v.field("Range", s.range, 2);
    v.field("Symmetric", s.symmetric);
}

void visit(viamd::Visitor& v, Thing& t) {
    v.field("Name", t.name, sizeof(t.name));
    if (v.field("Kind", t.kind)) {
        // Set up by the kind: what the file has after it, it keeps
        t.set_up = 100 + t.kind;
        t.shade = Shade::Dark;
    }
    v.field("Shade", t.shade, viamd::enum_names(shade_names));
    v.field("Iso", t.iso);
    v.field("SetUp", t.set_up);
    if (v.when(t.has_scale) && v.block("Scale")) {
        visit(v, t.scale);
        v.end();
    }
    v.list("Values", t.values, t.num_values, 8);
    v.list("Curve", t.curve, t.num_curve, 4, 2);
    for (int i = 0; v.each("Item", i, t.num_items, 4); ) {
        Item& it = t.items[i];
        v.field("Pos", it.pos);
        v.field("On", it.on);
        if (v.field("Opacity", it.opacity)) it.opacity = MD_CLAMP(it.opacity, 0.0f, 1.0f);
    }
    if (v.section("Extra")) {
        v.field("Color", t.color);
        v.field("Detail", t.detail);
        v.end();
    }
}

// What the reader of [Thing] does: its entries and those of [Thing.Extra], each through the list
size_t read_thing(Thing& t, viamd::deserialization_state_t& deser) {
    size_t unknown = 0;
    viamd::read_section(deser, MD_STR_LIT("Thing"), [&](md_str_t path, md_str_t arg) {
        viamd::Visitor v = viamd::Visitor::reader(path, arg);
        visit(v, t);
        if (!v.consumed()) unknown += 1;
        return v.consumed();
    });
    return unknown;
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

void write_thing(text_t& txt, const Thing& t) {
    viamd::write_section_header(txt.state, MD_STR_LIT("Thing"));
    viamd::Visitor v = viamd::Visitor::writer(txt.state, MD_STR_LIT("Thing"));
    visit(v, const_cast<Thing&>(t));
}

viamd::deserialization_state_t reader_of(md_str_t text) {
    viamd::deserialization_state_t d = {};
    d.filename = MD_STR_LIT("test.via");
    d.text = text;
    return d;
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
    t.color = md_float4_t{0.1f, 0.2f, 0.3f, 0.4f};
    t.detail = 1;
    return t;
}

}  // namespace

UTEST(viamd_serialization_visitor, round_trip) {
    const Thing src = a_thing();
    text_t txt;
    write_thing(txt, src);

    // Paths: short keys within the blocks, the section header once
    EXPECT_TRUE(txt.has_line("[Thing]"));
    EXPECT_TRUE(txt.has_line("Name=orbital 12"));
    EXPECT_TRUE(txt.has_line("Shade=Light"));
    EXPECT_TRUE(txt.has_line("Scale.Colormap=Viridis"));
    EXPECT_TRUE(txt.has_line("Scale.Range=-0.5,0.5"));
    EXPECT_TRUE(txt.has_line("Values=-0.25,0,0.75"));
    EXPECT_TRUE(txt.has_line("Curve=0,0,0.5,0.25,1,1"));
    EXPECT_TRUE(txt.has_line("Item=2"));
    EXPECT_TRUE(txt.has_line("Item.1.Opacity=0.25"));
    EXPECT_TRUE(txt.has_line("[Thing.Extra]"));
    EXPECT_TRUE(txt.has_line("Color=0.100000001,0.200000003,0.300000012,0.400000006"));

    viamd::deserialization_state_t deser = reader_of(txt.str());
    md_str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, deser));
    ASSERT_TRUE(md_str_eq(section, MD_STR_LIT("Thing")));

    Thing dst;
    dst.has_scale = false;
    EXPECT_EQ(0u, (uint32_t)read_thing(dst, deser));
    EXPECT_FALSE(viamd::next_section_header(section, deser));

    EXPECT_STREQ(src.name, dst.name);
    EXPECT_EQ(src.kind, dst.kind);
    // The kind set the shade up; the shade read after it is the one kept
    EXPECT_EQ((int)Shade::Light, (int)dst.shade);
    EXPECT_EQ(src.set_up, dst.set_up);
    EXPECT_NEAR(src.iso, dst.iso, 1e-12);
    EXPECT_EQ(src.scale.colormap, dst.scale.colormap);
    EXPECT_EQ(src.scale.range[0], dst.scale.range[0]);
    EXPECT_EQ(src.scale.range[1], dst.scale.range[1]);
    EXPECT_TRUE(dst.scale.symmetric);
    ASSERT_EQ(3, dst.num_values);
    EXPECT_EQ(0.75f, dst.values[2]);
    ASSERT_EQ(3, dst.num_curve);
    EXPECT_EQ(0.25f, dst.curve[3]);
    ASSERT_EQ(2, dst.num_items);
    EXPECT_FALSE(dst.items[0].on);
    EXPECT_EQ(3.0f, dst.items[0].pos.z);
    EXPECT_EQ(0.25f, dst.items[1].opacity);
    EXPECT_EQ(0.4f, dst.color.w);
    EXPECT_EQ(1, dst.detail);
}

UTEST(viamd_serialization_visitor, written_only_when) {
    Thing src = a_thing();
    src.has_scale = false;
    src.num_items = 0;
    text_t txt;
    write_thing(txt, src);
    EXPECT_FALSE(txt.has_line("Scale.Range=-0.5,0.5"));
    // An empty array says so
    EXPECT_TRUE(txt.has_line("Item=0"));
}

UTEST(viamd_serialization_visitor, file_order_is_read_order) {
    // The kind sets the shade up: written before it, the shade of the file wins; after it, the kind's
    Thing t;
    viamd::Visitor a = viamd::Visitor::reader(MD_STR_LIT("Shade"), MD_STR_LIT("Light"));
    visit(a, t);
    viamd::Visitor b = viamd::Visitor::reader(MD_STR_LIT("Kind"), MD_STR_LIT("1"));
    visit(b, t);
    EXPECT_EQ((int)Shade::Dark, (int)t.shade);
    EXPECT_EQ(101, t.set_up);
}

UTEST(viamd_serialization_visitor, enums_by_name_and_by_index) {
    Thing t;
    viamd::Visitor a = viamd::Visitor::reader(MD_STR_LIT("Shade"), MD_STR_LIT("Dark"));
    visit(a, t);
    EXPECT_EQ((int)Shade::Dark, (int)t.shade);

    // The index, as an older workspace has it
    viamd::Visitor b = viamd::Visitor::reader(MD_STR_LIT("Shade"), MD_STR_LIT("0"));
    visit(b, t);
    EXPECT_EQ((int)Shade::Light, (int)t.shade);

    // Neither: the entry is the field's, and the value stays
    viamd::Visitor c = viamd::Visitor::reader(MD_STR_LIT("Shade"), MD_STR_LIT("Purple"));
    visit(c, t);
    EXPECT_TRUE(c.consumed());
    EXPECT_EQ((int)Shade::Light, (int)t.shade);
    viamd::Visitor d = viamd::Visitor::reader(MD_STR_LIT("Shade"), MD_STR_LIT("2"));
    visit(d, t);
    EXPECT_EQ((int)Shade::Light, (int)t.shade);

    // A colour map by its name, or by an index that has one
    viamd::Visitor e = viamd::Visitor::reader(MD_STR_LIT("Scale.Colormap"), MD_STR_LIT("Paired"));
    visit(e, t);
    EXPECT_EQ(3, t.scale.colormap);
    viamd::Visitor f = viamd::Visitor::reader(MD_STR_LIT("Scale.Colormap"), MD_STR_LIT("1"));
    visit(f, t);
    EXPECT_EQ(1, t.scale.colormap);
    viamd::Visitor g = viamd::Visitor::reader(MD_STR_LIT("Scale.Colormap"), MD_STR_LIT("17"));
    visit(g, t);
    EXPECT_EQ(1, t.scale.colormap);
}

UTEST(viamd_serialization_visitor, paths_are_matched_whole) {
    Thing t;
    // A field's name within a block is not the object's field of that name
    viamd::Visitor a = viamd::Visitor::reader(MD_STR_LIT("Scale.Name"), MD_STR_LIT("x"));
    visit(a, t);
    EXPECT_FALSE(a.consumed());
    EXPECT_STREQ("thing", t.name);

    // Nor the other way around
    viamd::Visitor b = viamd::Visitor::reader(MD_STR_LIT("Symmetric"), MD_STR_LIT("1"));
    visit(b, t);
    EXPECT_FALSE(b.consumed());
    EXPECT_FALSE(t.scale.symmetric);

    // A block is not a value, and a value has no steps after it
    viamd::Visitor c = viamd::Visitor::reader(MD_STR_LIT("Scale"), MD_STR_LIT("1"));
    visit(c, t);
    EXPECT_FALSE(c.consumed());
    viamd::Visitor d = viamd::Visitor::reader(MD_STR_LIT("Iso.Value"), MD_STR_LIT("1"));
    visit(d, t);
    EXPECT_FALSE(d.consumed());
    EXPECT_NEAR(0.05, t.iso, 1e-12);

    // Unknown
    viamd::Visitor e = viamd::Visitor::reader(MD_STR_LIT("Nonsense"), MD_STR_LIT("1"));
    visit(e, t);
    EXPECT_FALSE(e.consumed());
}

UTEST(viamd_serialization_visitor, arrays) {
    Thing t;
    // The count first, then items: fewer items than the count leaves the rest as they were
    viamd::Visitor a = viamd::Visitor::reader(MD_STR_LIT("Item"), MD_STR_LIT("3"));
    visit(a, t);
    EXPECT_TRUE(a.consumed());
    EXPECT_EQ(3, t.num_items);

    // An item beyond the count makes it at least that long
    viamd::Visitor b = viamd::Visitor::reader(MD_STR_LIT("Item.3.Opacity"), MD_STR_LIT("2"));
    visit(b, t);
    EXPECT_TRUE(b.consumed());
    EXPECT_EQ(4, t.num_items);
    EXPECT_EQ(1.0f, t.items[3].opacity);     // clamped as it was read

    // Beyond the capacity, or not an index: in the array but nothing of it
    viamd::Visitor c = viamd::Visitor::reader(MD_STR_LIT("Item.4.Opacity"), MD_STR_LIT("0.5"));
    visit(c, t);
    EXPECT_FALSE(c.consumed());
    viamd::Visitor d = viamd::Visitor::reader(MD_STR_LIT("Item.x.Opacity"), MD_STR_LIT("0.5"));
    visit(d, t);
    EXPECT_FALSE(d.consumed());
    viamd::Visitor e = viamd::Visitor::reader(MD_STR_LIT("Item.0"), MD_STR_LIT("0.5"));
    visit(e, t);
    EXPECT_FALSE(e.consumed());

    // A field an item does not have makes no item
    viamd::Visitor h = viamd::Visitor::reader(MD_STR_LIT("Item"), MD_STR_LIT("1"));
    visit(h, t);
    viamd::Visitor k = viamd::Visitor::reader(MD_STR_LIT("Item.2.Bogus"), MD_STR_LIT("1"));
    visit(k, t);
    EXPECT_FALSE(k.consumed());
    EXPECT_EQ(1, t.num_items);

    // A count is kept within the capacity
    viamd::Visitor f = viamd::Visitor::reader(MD_STR_LIT("Item"), MD_STR_LIT("99"));
    visit(f, t);
    EXPECT_EQ(4, t.num_items);
    viamd::Visitor g = viamd::Visitor::reader(MD_STR_LIT("Item"), MD_STR_LIT("0"));
    visit(g, t);
    EXPECT_EQ(0, t.num_items);
}

UTEST(viamd_serialization_visitor, lists) {
    Thing t;
    // A whole list replaces the one there was, whatever its length
    viamd::Visitor a = viamd::Visitor::reader(MD_STR_LIT("Curve"), MD_STR_LIT("0,0.5"));
    visit(a, t);
    EXPECT_EQ(1, t.num_curve);
    EXPECT_EQ(0.5f, t.curve[1]);

    // Of pairs, an odd count is not one
    viamd::Visitor b = viamd::Visitor::reader(MD_STR_LIT("Curve"), MD_STR_LIT("0,0.5,1"));
    visit(b, t);
    EXPECT_TRUE(b.consumed());
    EXPECT_EQ(1, t.num_curve);

    // One item, as the numbered keys of an older workspace have it
    viamd::Visitor c = viamd::Visitor::reader(MD_STR_LIT("Curve.2"), MD_STR_LIT("0.75,0.5"));
    visit(c, t);
    EXPECT_EQ(3, t.num_curve);
    EXPECT_EQ(0.75f, t.curve[4]);

    // Not an item of it: not the list's
    viamd::Visitor g = viamd::Visitor::reader(MD_STR_LIT("Curve.x"), MD_STR_LIT("1,1"));
    visit(g, t);
    EXPECT_FALSE(g.consumed());
    EXPECT_EQ(3, t.num_curve);

    // Too many, or not numbers: as it was
    viamd::Visitor d = viamd::Visitor::reader(MD_STR_LIT("Values"), MD_STR_LIT("1,2,3,4,5,6,7,8,9"));
    visit(d, t);
    EXPECT_EQ(0, t.num_values);
    viamd::Visitor e = viamd::Visitor::reader(MD_STR_LIT("Values"), MD_STR_LIT("1,two"));
    visit(e, t);
    EXPECT_EQ(0, t.num_values);

    // The empty list
    t.num_values = 2;
    viamd::Visitor f = viamd::Visitor::reader(MD_STR_LIT("Values"), MD_STR_LIT(""));
    visit(f, t);
    EXPECT_EQ(0, t.num_values);
}

UTEST(viamd_serialization_visitor, continuations_follow_what_they_continue) {
    // A [Thing.Extra] continues the [Thing] just before it; one after another section is nobody's
    const md_str_t text = MD_STR_LIT(
        "[Thing]\n"
        "Name=first\n"
        "[Thing.Extra]\n"
        "Detail=5\n"
        "[Other]\n"
        "A=1\n"
        "[Thing.Extra]\n"
        "Detail=9\n"
        "[Thing]\n"
        "Name=second\n");
    viamd::deserialization_state_t deser = reader_of(text);
    md_str_t section;

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    Thing first;
    read_thing(first, deser);
    EXPECT_EQ(5, first.detail);

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(md_str_eq(section, MD_STR_LIT("Other")));
    // Nothing continues [Other]
    md_str_t sub;
    EXPECT_FALSE(viamd::next_subsection(sub, MD_STR_LIT("Other"), deser));

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(md_str_eq(section, MD_STR_LIT("Thing.Extra")));

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(md_str_eq(section, MD_STR_LIT("Thing")));
    Thing second;
    read_thing(second, deser);
    EXPECT_STREQ("second", second.name);
    EXPECT_EQ(2, second.detail);

    // A section whose name only begins like it does not continue it
    const md_str_t text2 = MD_STR_LIT("[Thing]\n[ThingExtra]\nDetail=3\n");
    viamd::deserialization_state_t deser2 = reader_of(text2);
    ASSERT_TRUE(viamd::next_section_header(section, deser2));
    EXPECT_FALSE(viamd::next_subsection(sub, MD_STR_LIT("Thing"), deser2));
    ASSERT_TRUE(viamd::next_section_header(section, deser2));
    EXPECT_TRUE(md_str_eq(section, MD_STR_LIT("ThingExtra")));
}

UTEST(viamd_serialization_visitor, aliases_of_older_keys) {
    static const viamd::PathAlias aliases[] = {
        { "ThingColor",          "Extra.Color" },
        { "ThingItem#",          "Item.#.Pos" },
        { "ThingItemOpacity#",   "Item.#.Opacity" },
        { "ThingValue#",         "Values.#" },
    };
    char buf[64];
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingColor"), aliases), MD_STR_LIT("Extra.Color")));
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItem12"), aliases), MD_STR_LIT("Item.12.Pos")));
    // The longer family is not taken for the shorter one, whichever comes first
    EXPECT_TRUE(md_str_eq(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItemOpacity3"), aliases), MD_STR_LIT("Item.3.Opacity")));
    EXPECT_TRUE(md_str_empty(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItem"), aliases)));
    EXPECT_TRUE(md_str_empty(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingItemX3"), aliases)));
    EXPECT_TRUE(md_str_empty(viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingColour"), aliases)));
    // What does not fit is none
    char small[8];
    EXPECT_TRUE(md_str_empty(viamd::path_alias(small, sizeof(small), MD_STR_LIT("ThingItemOpacity3"), aliases)));

    // And read through the list, as the key it now is
    Thing t;
    const md_str_t path = viamd::path_alias(buf, sizeof(buf), MD_STR_LIT("ThingValue1"), aliases);
    viamd::Visitor v = viamd::Visitor::reader(path, MD_STR_LIT("0.5"));
    visit(v, t);
    EXPECT_EQ(2, t.num_values);
    EXPECT_EQ(0.5f, t.values[1]);
}

UTEST(viamd_serialization_visitor, path_helpers) {
    md_str_t path = MD_STR_LIT("Slice.0.Opacity");
    md_str_t step;
    ASSERT_TRUE(viamd::path_next(step, path));
    EXPECT_TRUE(md_str_eq(step, MD_STR_LIT("Slice")));
    ASSERT_TRUE(viamd::path_next(step, path));
    EXPECT_TRUE(md_str_eq(step, MD_STR_LIT("0")));
    ASSERT_TRUE(viamd::path_next(step, path));
    EXPECT_TRUE(md_str_eq(step, MD_STR_LIT("Opacity")));
    EXPECT_FALSE(viamd::path_next(step, path));

    char buf[16];
    EXPECT_TRUE(md_str_eq(viamd::path_join(buf, sizeof(buf), MD_STR_LIT("Field"), MD_STR_LIT("Range")), MD_STR_LIT("Field.Range")));
    EXPECT_TRUE(md_str_eq(viamd::path_join(buf, sizeof(buf), md_str_t{}, MD_STR_LIT("Range")), MD_STR_LIT("Range")));
    EXPECT_TRUE(md_str_empty(viamd::path_join(buf, sizeof(buf), MD_STR_LIT("ElectronicStructure"), MD_STR_LIT("Iso"))));

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

UTEST(viamd_serialization_visitor, read_entry_falls_back_to_aliases_then_legacy) {
    static const viamd::PathAlias aliases[] = {
        { "ThingIso",        "Iso" },
        { "ThingValue#",     "Values.#" },
    };
    Thing t;
    int legacy_calls = 0;
    auto read = [&](md_str_t path, md_str_t arg) {
        return viamd::read_entry(path, arg, aliases,
            [&](viamd::Visitor& v) { visit(v, t); },
            [&](md_str_t key, md_str_t) {
                legacy_calls += 1;
                // A key whose meaning changed: read as it was meant
                if (md_str_eq_cstr(key, "ThingDoubleIso")) { t.iso = 2 * 0.5; return true; }
                return false;
            });
    };
    // A path first, then an alias: neither reaches the legacy keys
    EXPECT_TRUE(read(MD_STR_LIT("Iso"), MD_STR_LIT("0.25")));
    EXPECT_NEAR(0.25, t.iso, 1e-12);
    EXPECT_TRUE(read(MD_STR_LIT("ThingIso"), MD_STR_LIT("0.5")));
    EXPECT_NEAR(0.5, t.iso, 1e-12);
    EXPECT_TRUE(read(MD_STR_LIT("ThingValue0"), MD_STR_LIT("3")));
    EXPECT_EQ(3.0f, t.values[0]);
    EXPECT_EQ(0, legacy_calls);

    EXPECT_TRUE(read(MD_STR_LIT("ThingDoubleIso"), MD_STR_LIT("x")));
    EXPECT_NEAR(1.0, t.iso, 1e-12);
    EXPECT_FALSE(read(MD_STR_LIT("Unknown"), MD_STR_LIT("1")));
    EXPECT_EQ(2, legacy_calls);

    // A path that is no field is not taken for a key from before the paths
    EXPECT_FALSE(read(MD_STR_LIT("Scale.Unknown"), MD_STR_LIT("1")));
    EXPECT_EQ(2, legacy_calls);
}
