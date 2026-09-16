#include "utest.h"

#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_str_builder.h>

// The workspace file format: sections in brackets, then ident=value lines under them. These tests
// go through the writer and back through the reader, because that is the pairing that has to hold.

namespace {

struct writer_t {
    viamd::serialization_state_t state = {};

    writer_t(md_allocator_i* alloc) {
        state.filename = STR_LIT("test.viamd");
        md_strb_init(&state.sb, alloc);
    }
    ~writer_t() { md_strb_free(&state.sb); }

    viamd::deserialization_state_t reader() {
        viamd::deserialization_state_t d = {};
        d.filename = STR_LIT("test.viamd");
        d.text = md_strb_to_str(state.sb);
        return d;
    }
};

// Collect every ident=arg pair of the next section, in file order
struct entry_t {
    char ident[64];
    char arg[512];
};

size_t read_entries(entry_t* out, size_t cap, viamd::deserialization_state_t& deser) {
    size_t count = 0;
    str_t ident, arg;
    while (count < cap && viamd::next_entry(ident, arg, deser)) {
        snprintf(out[count].ident, sizeof(out[count].ident), STR_FMT, STR_ARG(ident));
        snprintf(out[count].arg,   sizeof(out[count].arg),   STR_FMT, STR_ARG(arg));
        count += 1;
    }
    return count;
}

}  // namespace

UTEST(viamd_serialization, scalar_round_trip) {
    md_allocator_i* alloc = md_get_heap_allocator();
    writer_t w(alloc);

    viamd::write_section_header(w.state, STR_LIT("Camera"));
    viamd::write_int (w.state, STR_LIT("Frame"),   1234);
    viamd::write_bool(w.state, STR_LIT("Ortho"),   true);
    viamd::write_flt (w.state, STR_LIT("Near"),    0.25f);
    viamd::write_dbl (w.state, STR_LIT("Far"),     1000.5);
    viamd::write_vec3(w.state, STR_LIT("Pos"),     vec3_t{1.5f, -2.25f, 3.0f});
    viamd::write_quat(w.state, STR_LIT("Ori"),     quat_t{0.0f, 0.0f, 0.0f, 1.0f});
    viamd::write_str (w.state, STR_LIT("Label"),   STR_LIT("a simple label"));

    viamd::deserialization_state_t deser = w.reader();

    str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(str_eq(section, STR_LIT("Camera")));

    entry_t e[16];
    const size_t n = read_entries(e, 16, deser);
    ASSERT_EQ(7u, (uint32_t)n);

    int i = 0;
    EXPECT_TRUE(viamd::extract_int(i, str_from_cstr(e[0].arg)));
    EXPECT_EQ(1234, i);

    bool b = false;
    EXPECT_TRUE(viamd::extract_bool(b, str_from_cstr(e[1].arg)));
    EXPECT_TRUE(b);

    float f = 0;
    EXPECT_TRUE(viamd::extract_flt(f, str_from_cstr(e[2].arg)));
    EXPECT_NEAR(0.25, (double)f, 1.0e-6);

    double d = 0;
    EXPECT_TRUE(viamd::extract_dbl(d, str_from_cstr(e[3].arg)));
    EXPECT_NEAR(1000.5, d, 1.0e-6);

    vec3_t v = {0, 0, 0};
    viamd::extract_vec3(v, str_from_cstr(e[4].arg));
    EXPECT_NEAR( 1.50, (double)v.x, 1.0e-5);
    EXPECT_NEAR(-2.25, (double)v.y, 1.0e-5);
    EXPECT_NEAR( 3.00, (double)v.z, 1.0e-5);

    quat_t q = {1, 1, 1, 1};
    viamd::extract_quat(q, str_from_cstr(e[5].arg));
    EXPECT_NEAR(0.0, (double)q.x, 1.0e-5);
    EXPECT_NEAR(1.0, (double)q.w, 1.0e-5);

    str_t s = {};
    EXPECT_TRUE(viamd::extract_str(s, str_from_cstr(e[6].arg)));
    EXPECT_TRUE(str_eq(s, STR_LIT("a simple label")));
}

UTEST(viamd_serialization, sections_are_walked_in_order) {
    md_allocator_i* alloc = md_get_heap_allocator();
    writer_t w(alloc);

    viamd::write_section_header(w.state, STR_LIT("First"));
    viamd::write_int(w.state, STR_LIT("A"), 1);
    viamd::write_section_header(w.state, STR_LIT("Second"));
    viamd::write_int(w.state, STR_LIT("B"), 2);
    viamd::write_int(w.state, STR_LIT("C"), 3);

    viamd::deserialization_state_t deser = w.reader();
    entry_t e[8];
    str_t section;

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(str_eq(section, STR_LIT("First")));
    // next_entry must stop at the following section header rather than running past it
    EXPECT_EQ(1u, (uint32_t)read_entries(e, 8, deser));

    ASSERT_TRUE(viamd::next_section_header(section, deser));
    EXPECT_TRUE(str_eq(section, STR_LIT("Second")));
    EXPECT_EQ(2u, (uint32_t)read_entries(e, 8, deser));

    EXPECT_FALSE(viamd::next_section_header(section, deser));
}

UTEST(viamd_serialization, multiline_string_round_trip) {
    md_allocator_i* alloc = md_get_heap_allocator();
    writer_t w(alloc);

    // Script text is written this way, and it is the only value that may contain newlines
    const str_t script = STR_LIT("d1 = distance(1,2);\nsel = residue(1:10);\n");

    viamd::write_section_header(w.state, STR_LIT("Script"));
    viamd::write_str(w.state, STR_LIT("Text"), script);
    viamd::write_int(w.state, STR_LIT("After"), 7);

    viamd::deserialization_state_t deser = w.reader();
    str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, deser));

    entry_t e[8];
    const size_t n = read_entries(e, 8, deser);
    ASSERT_EQ(2u, (uint32_t)n);
    EXPECT_STREQ("Text", e[0].ident);
    EXPECT_STREQ("d1 = distance(1,2);\nsel = residue(1:10);\n", e[0].arg);

    // The entry after a multiline string has to survive it
    EXPECT_STREQ("After", e[1].ident);
    int i = 0;
    EXPECT_TRUE(viamd::extract_int(i, str_from_cstr(e[1].arg)));
    EXPECT_EQ(7, i);
}

UTEST(viamd_serialization, bitfield_round_trip) {
    md_allocator_i* alloc = md_get_heap_allocator();
    writer_t w(alloc);

    md_bitfield_t src = md_bitfield_create(alloc);
    md_bitfield_t dst = md_bitfield_create(alloc);
    const int bits[] = { 0, 1, 5, 63, 64, 65, 1000, 4095 };
    for (size_t i = 0; i < ARRAY_SIZE(bits); ++i) md_bitfield_set_bit(&src, bits[i]);

    viamd::write_section_header(w.state, STR_LIT("Selection"));
    viamd::write_bitfield(w.state, STR_LIT("Mask"), &src);

    viamd::deserialization_state_t deser = w.reader();
    str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, deser));

    entry_t e[4];
    ASSERT_EQ(1u, (uint32_t)read_entries(e, 4, deser));
    EXPECT_STREQ("Mask", e[0].ident);
    ASSERT_TRUE(viamd::extract_bitfield(&dst, str_from_cstr(e[0].arg)));

    EXPECT_EQ(md_bitfield_popcount(&src), md_bitfield_popcount(&dst));
    for (size_t i = 0; i < ARRAY_SIZE(bits); ++i) {
        EXPECT_TRUE(md_bitfield_test_bit(&dst, bits[i]));
    }

    md_bitfield_free(&src);
    md_bitfield_free(&dst);
}

UTEST(viamd_serialization, entry_after_a_bitfield_survives) {
    // A bitfield is written as one long base64 blob between ### markers. It used to be written
    // without a trailing newline, which glued the next entry onto the same line - and since an
    // entry is parsed by splitting at the first '=', that entry was silently swallowed. Every NTO
    // group label after the first was lost from a saved workspace this way.
    md_allocator_i* alloc = md_get_heap_allocator();
    writer_t w(alloc);

    md_bitfield_t mask = md_bitfield_create(alloc);
    md_bitfield_set_range(&mask, 0, 32);

    viamd::write_section_header(w.state, STR_LIT("Groups"));
    viamd::write_str     (w.state, STR_LIT("GroupLabel"), STR_LIT("first"));
    viamd::write_bitfield(w.state, STR_LIT("GroupAtoms"), &mask);
    viamd::write_str     (w.state, STR_LIT("GroupLabel"), STR_LIT("second"));
    viamd::write_bitfield(w.state, STR_LIT("GroupAtoms"), &mask);

    viamd::deserialization_state_t deser = w.reader();
    str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, deser));

    entry_t e[8];
    const size_t n = read_entries(e, 8, deser);
    ASSERT_EQ(4u, (uint32_t)n);
    EXPECT_STREQ("GroupLabel", e[0].ident);
    EXPECT_STREQ("first",      e[0].arg);
    EXPECT_STREQ("GroupAtoms", e[1].ident);
    EXPECT_STREQ("GroupLabel", e[2].ident);
    EXPECT_STREQ("second",     e[2].arg);
    EXPECT_STREQ("GroupAtoms", e[3].ident);

    md_bitfield_free(&mask);
}

UTEST(viamd_serialization, malformed_values_are_rejected) {
    // A value that does not parse must be reported and must leave the destination alone, so a
    // corrupt workspace degrades to defaults instead of to garbage.
    int i = 42;
    EXPECT_FALSE(viamd::extract_int(i, STR_LIT("not_a_number")));
    EXPECT_EQ(42, i);

    bool b = false;
    EXPECT_FALSE(viamd::extract_bool(b, STR_LIT("2")));
    EXPECT_FALSE(b);

    double d = 1.5;
    EXPECT_FALSE(viamd::extract_dbl(d, STR_LIT("")));
    EXPECT_NEAR(1.5, d, 1.0e-9);

    md_allocator_i* alloc = md_get_heap_allocator();
    md_bitfield_t bf = md_bitfield_create(alloc);
    EXPECT_FALSE(viamd::extract_bitfield(&bf, STR_LIT("no markers here")));
    md_bitfield_free(&bf);
}
