#include "utest.h"

#include <script_engine.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>

#include <string.h>

/* The script interface under v2 (VIAL), as viamd calls it: what a program says when compiled, its outputs as plots
 * see them, the marks of the editor and the references they draw by. On the alanine peptide, its models published as
 * the run "run". */

namespace {

struct Peptide {
    md_allocator_i*   alloc = nullptr;
    md_system_t       sys = {};
    md_system_state_t state = {};
    size_t            num_frames = 0;
};

const Peptide* peptide() {
    static Peptide p;
    static bool loaded = false;
    if (!loaded) {
        loaded = true;
        p.alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
        p.sys = { .alloc = p.alloc };
        p.state = { .alloc = p.alloc };
        const str_t path = STR_LIT(VIAMD_TEST_DATA_DIR "/1ALA-560ns.pdb");
        if (md_pdb_system_init_from_file(&p.sys, &p.state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) &&
            md_util_system_infer(&p.sys, &p.state, MD_UTIL_INFER_ALL) &&
            md_pdb_system_publish_run(&p.sys, path, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
            const md_attribute_t* t = md_attributes_find(&p.sys.attributes, STR_LIT("run/time"));
            p.num_frames = t ? (size_t)t->format.shape[0] : 0;
        }
    }
    return p.num_frames ? &p : nullptr;
}

script::Program* compile_v2(const char* src, md_allocator_i* alloc, const script::NamedSet* sets = nullptr, size_t num_sets = 0) {
    script::CompileDesc d;
    d.sys = &peptide()->sys;
    d.sets = sets;
    d.num_sets = num_sets;
    d.run = STR_LIT("run");
    return script::compile(script::Language::V2, str_from_cstr(src), d, alloc);
}

bool has_name(const str_t* names, size_t n, const char* name) {
    for (size_t i = 0; i < n; ++i) if (str_eq(names[i], str_from_cstr(name))) return true;
    return false;
}

}  // namespace

UTEST(script_v2, compile_and_outputs) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "m  = mean(d, axis='frame')\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(p != nullptr);
    EXPECT_EQ((int)script::Language::V2, (int)script::language(p));
    EXPECT_TRUE(script::valid(p));
    EXPECT_EQ(0u, script::num_errors(p));
    EXPECT_NE(0u, script::fingerprint(p));

    // Outputs as plots see them: values, not the selection
    const size_t no = script::num_outputs(p);
    const str_t* names = script::output_names(p);
    EXPECT_EQ(3u, no);
    EXPECT_TRUE(has_name(names, no, "d") && has_name(names, no, "h") && has_name(names, no, "m"));
    EXPECT_FALSE(has_name(names, no, "ca"));
    EXPECT_TRUE(script::output_flags(p, STR_LIT("d")) & script::OUTPUT_TEMPORAL);
    EXPECT_TRUE(script::output_flags(p, STR_LIT("h")) & script::OUTPUT_DISTRIBUTION);
    EXPECT_EQ(0u, script::output_flags(p, STR_LIT("ca")));
    const script::Range r = script::output_range(p, STR_LIT("d"));
    EXPECT_EQ(0, strncmp(src + r.beg, "d", (size_t)(r.end - r.beg)));
    EXPECT_NE(0u, script::output_fingerprint(p, STR_LIT("d")));
    const script::Ref ref = script::output_ref(p, STR_LIT("d"));
    EXPECT_TRUE(script::ref_valid(p, ref));
    EXPECT_TRUE(str_eq(STR_LIT("d"), script::ref_name(p, ref)));
    EXPECT_EQ(1, script::ref_dim(p, ref));
    EXPECT_TRUE(has_name(script::identifiers(p), script::num_identifiers(p), "ca"));

    // Marks: every expression, with its type, each drawable by its ref
    const size_t nm = script::num_marks(p);
    const script::Mark* mk = script::marks(p);
    ASSERT_GT(nm, 5u);
    bool call = false;
    for (size_t i = 0; i < nm; ++i) {
        EXPECT_TRUE(script::ref_valid(p, mk[i].ref));
        if (strncmp(src + mk[i].range.beg, "distance(ca[0], ca[1])", (size_t)(mk[i].range.end - mk[i].range.beg)) == 0 &&
            mk[i].range.end - mk[i].range.beg == 22) {
            call = true;
            EXPECT_TRUE(strstr(mk[i].text.ptr, "float") != nullptr);
        }
    }
    EXPECT_TRUE(call);

    // Refs belong to their compilation; renamed and spaced, the same fingerprints
    script::Program* q = compile_v2("ca = particles(name='CA')\ndist  = distance(ca[0],   ca[1])\nh = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\nm = mean(dist, axis='frame')\n", alloc);
    EXPECT_FALSE(script::ref_valid(q, ref));
    EXPECT_EQ(script::output_fingerprint(p, STR_LIT("d")), script::output_fingerprint(q, STR_LIT("dist")));
    EXPECT_NE(script::fingerprint(p), script::fingerprint(q));     // d is called dist
    script::program_free(q);
    script::program_free(p);
}

UTEST(script_v2, diagnostics_and_sets) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    script::Program* p = compile_v2("x = particles(nme='CA')\ny = 1 +\n", alloc);
    EXPECT_FALSE(script::valid(p));
    ASSERT_GT(script::num_errors(p), 0u);
    const script::Diagnostic* e = script::errors(p);
    EXPECT_EQ('E', e[0].text.ptr[0]);
    EXPECT_GE(e[0].range.beg, 0);
    script::program_free(p);

    // Stored selections: names the script reads; empty source compiles nothing and has their names alone
    const uint32_t idx[] = {0, 4};
    md_bitset_t atoms = md_bitset_from_indices(idx, 2, alloc);
    script::NamedSet sets[1] = {{STR_LIT("mysel"), atoms}};
    p = compile_v2("n = count(mysel)\nc = mysel.com\n", alloc, sets, 1);
    EXPECT_TRUE(script::valid(p));
    EXPECT_TRUE(has_name(script::identifiers(p), script::num_identifiers(p), "mysel"));
    script::program_free(p);
    p = compile_v2("", alloc, sets, 1);
    EXPECT_FALSE(script::valid(p));
    EXPECT_EQ(0u, script::num_errors(p));
    EXPECT_EQ(1u, script::num_identifiers(p));
    script::program_free(p);
    md_bitset_free(&atoms, alloc);

    // v1 is as it was, through the same calls
    script::CompileDesc d;
    d.sys = &peptide()->sys;
    p = script::compile(script::Language::V1, STR_LIT("x = distance(1, 2);"), d, alloc);
    EXPECT_EQ((int)script::Language::V1, (int)script::language(p));
    EXPECT_TRUE(script::valid(p));
    script::program_free(p);
}

UTEST(script_v2, editor) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    const script::Language V2 = script::Language::V2;
    EXPECT_TRUE(has_name(script::keywords(V2), script::num_keywords(V2), "def"));
    const size_t nb = script::builtin_names(V2, nullptr, 0);
    str_t* names = (str_t*)md_alloc(alloc, nb * sizeof(str_t));
    EXPECT_EQ(nb, script::builtin_names(V2, names, nb));
    EXPECT_TRUE(has_name(names, nb, "distance") && has_name(names, nb, "particles"));
    EXPECT_TRUE(script::name_valid(V2, STR_LIT("my_sel")));
    EXPECT_FALSE(script::name_valid(V2, STR_LIT("protein")));
    EXPECT_FALSE(script::name_valid(V2, STR_LIT("for")));

    // Completion: a selector's keyword, and the system's values in a string
    const str_t src = STR_LIT("x = particles(res");
    script::Completions c = script::complete(V2, src, (int)src.len, &peptide()->sys, alloc);
    bool found = false;
    for (size_t i = 0; i < c.count; ++i) {
        if (str_eq(c.items[i].label, STR_LIT("resname"))) {
            found = true;
            EXPECT_TRUE(str_eq(c.items[i].text, STR_LIT("resname=")));
            EXPECT_EQ((int)script::CompletionKind::Parameter, (int)c.items[i].kind);
            EXPECT_TRUE(c.items[i].fits);
        }
    }
    EXPECT_TRUE(found);
    EXPECT_EQ(14, c.range.beg);
    const str_t s2 = STR_LIT("x = particles(resname='AL')");
    c = script::complete(V2, s2, (int)s2.len - 2, &peptide()->sys, alloc);
    found = false;
    for (size_t i = 0; i < c.count; ++i) found |= str_eq(c.items[i].label, STR_LIT("ALA"));
    EXPECT_TRUE(found);

    // attr(): the path past the run's prefix, for an attribute of the run; nothing for another
    const md_attributes_t* table = &peptide()->sys.attributes;
    bool tried = false;
    for (size_t i = 0; i < md_array_size(table->attr); ++i) {
        const md_attribute_t* a = &table->attr[i];
        if (a->path.len > 4 && !strncmp(a->path.ptr, "run/", 4) && !str_eq(a->path, STR_LIT("run/time"))) {
            const str_t p = script::attr_path(V2, table, a);
            EXPECT_TRUE(str_eq(p, str_t{a->path.ptr + 4, a->path.len - 4}));
            tried = true;
        }
    }
    EXPECT_TRUE(tried);
    md_arena_allocator_destroy(alloc);
}
