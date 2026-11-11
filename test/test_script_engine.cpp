#include "utest.h"
#include "test_executor.h"

#include <script_engine.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_writer.h>

#include <core/md_os.h>

#include <string.h>
#include <vector>

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
        p.alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
        p.sys = { .alloc = p.alloc };
        p.state = { .alloc = p.alloc };
        const md_str_t path = MD_STR_LIT(VIAMD_TEST_DATA_DIR "/1ALA-560ns.pdb");
        if (md_pdb_system_init_from_file(&p.sys, &p.state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL) &&
            md_util_system_infer(&p.sys, &p.state, MD_UTIL_INFER_ALL, NULL) &&
            md_pdb_system_publish_run(&p.sys, path, MD_STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
            const md_attribute_t* t = md_attributes_find(&p.sys.attributes, MD_STR_LIT("run/time"));
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
    d.run = MD_STR_LIT("run");
    return script::compile(script::Language::V2, md_str_from_cstr(src), d, alloc);
}

bool has_name(const md_str_t* names, size_t n, const char* name) {
    for (size_t i = 0; i < n; ++i) if (md_str_eq(names[i], md_str_from_cstr(name))) return true;
    return false;
}

}  // namespace

UTEST(script_v2, compile_and_outputs) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    const char* src =
        "ca = particle(name='CA')\n"
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
    const md_str_t* names = script::output_names(p);
    EXPECT_EQ(3u, no);
    EXPECT_TRUE(has_name(names, no, "d") && has_name(names, no, "h") && has_name(names, no, "m"));
    EXPECT_FALSE(has_name(names, no, "ca"));
    EXPECT_TRUE(script::output_flags(p, MD_STR_LIT("d")) & script::OUTPUT_TEMPORAL);
    EXPECT_TRUE(script::output_flags(p, MD_STR_LIT("h")) & script::OUTPUT_DISTRIBUTION);
    EXPECT_EQ(0u, script::output_flags(p, MD_STR_LIT("ca")));
    const script::Range r = script::output_range(p, MD_STR_LIT("d"));
    EXPECT_EQ(0, strncmp(src + r.beg, "d", (size_t)(r.end - r.beg)));
    EXPECT_NE(0u, script::output_fingerprint(p, MD_STR_LIT("d")));
    const script::Ref ref = script::output_ref(p, MD_STR_LIT("d"));
    EXPECT_TRUE(script::ref_valid(p, ref));
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("d"), script::ref_name(p, ref)));
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
    script::Program* q = compile_v2("ca = particle(name='CA')\ndist  = distance(ca[0],   ca[1])\nh = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\nm = mean(dist, axis='frame')\n", alloc);
    EXPECT_FALSE(script::ref_valid(q, ref));
    EXPECT_EQ(script::output_fingerprint(p, MD_STR_LIT("d")), script::output_fingerprint(q, MD_STR_LIT("dist")));
    EXPECT_NE(script::fingerprint(p), script::fingerprint(q));     // d is called dist
    script::program_free(q);
    script::program_free(p);
}

UTEST(script_v2, diagnostics_and_sets) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    script::Program* p = compile_v2("x = particle(nme='CA')\ny = 1 +\n", alloc);
    EXPECT_FALSE(script::valid(p));
    ASSERT_GT(script::num_errors(p), 0u);
    const script::Diagnostic* e = script::errors(p);
    EXPECT_EQ('E', e[0].text.ptr[0]);
    EXPECT_GE(e[0].range.beg, 0);
    script::program_free(p);

    // Stored selections: names the script reads; empty source compiles nothing and has their names alone
    const uint32_t idx[] = {0, 4};
    md_bitset_t atoms = md_bitset_from_indices(idx, 2, alloc);
    script::NamedSet sets[1] = {{MD_STR_LIT("mysel"), atoms}};
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
    p = script::compile(script::Language::V1, MD_STR_LIT("x = distance(1, 2);"), d, alloc);
    EXPECT_EQ((int)script::Language::V1, (int)script::language(p));
    EXPECT_TRUE(script::valid(p));
    script::program_free(p);
}

UTEST(script_v2, editor) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    const script::Language V2 = script::Language::V2;
    EXPECT_TRUE(has_name(script::keywords(V2), script::num_keywords(V2), "def"));
    const size_t nb = script::builtin_names(V2, nullptr, 0);
    md_str_t* names = (md_str_t*)md_alloc(alloc, nb * sizeof(md_str_t));
    EXPECT_EQ(nb, script::builtin_names(V2, names, nb));
    EXPECT_TRUE(has_name(names, nb, "distance") && has_name(names, nb, "particle"));
    EXPECT_TRUE(script::name_valid(V2, MD_STR_LIT("my_sel")));
    EXPECT_FALSE(script::name_valid(V2, MD_STR_LIT("protein")));
    EXPECT_FALSE(script::name_valid(V2, MD_STR_LIT("for")));

    // Completion: a selector's keyword, and the system's values in a string
    const md_str_t src = MD_STR_LIT("x = particle(res");
    script::Completions c = script::complete(V2, src, (int)src.len, &peptide()->sys, alloc);
    bool found = false;
    for (size_t i = 0; i < c.count; ++i) {
        if (md_str_eq(c.items[i].label, MD_STR_LIT("resname"))) {
            found = true;
            EXPECT_TRUE(md_str_eq(c.items[i].text, MD_STR_LIT("resname=")));
            EXPECT_EQ((int)script::CompletionKind::Parameter, (int)c.items[i].kind);
            EXPECT_TRUE(c.items[i].fits);
        }
    }
    EXPECT_TRUE(found);
    EXPECT_EQ(13, c.range.beg);
    const md_str_t s2 = MD_STR_LIT("x = particle(resname='AL')");
    c = script::complete(V2, s2, (int)s2.len - 2, &peptide()->sys, alloc);
    found = false;
    for (size_t i = 0; i < c.count; ++i) found |= md_str_eq(c.items[i].label, MD_STR_LIT("ALA"));
    EXPECT_TRUE(found);

    // attr(): the path past the run's prefix, for an attribute of the run; nothing for another
    const md_attributes_t* table = &peptide()->sys.attributes;
    bool tried = false;
    for (size_t i = 0; i < md_array_size(table->attr); ++i) {
        const md_attribute_t* a = &table->attr[i];
        if (a->path.len > 4 && !strncmp(a->path.ptr, "run/", 4) && !md_str_eq(a->path, MD_STR_LIT("run/time"))) {
            const md_str_t p = script::attr_path(V2, table, a);
            EXPECT_TRUE(md_str_eq(p, md_str_t{a->path.ptr + 4, a->path.len - 4}));
            tried = true;
        }
    }
    EXPECT_TRUE(tried);
    md_arena_allocator_destroy(alloc);
}

namespace {
script::EvalDesc eval_desc() {
    script::EvalDesc d;
    d.num_frames = peptide()->num_frames;
    d.sys = &peptide()->sys;
    d.run = MD_STR_LIT("run");
    // The frames one thread at a time, out of order and over three thread_idx: what the evaluation spreads over
    // threads is all done, whichever order and worker the frames come in
    static test_executor_t shuffled = test_executor_shuffled(3);
    d.executor = test_executor(&shuffled);
    return d;
}
}  // namespace

UTEST(script_v2, evaluate) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t nf = peptide()->num_frames;
    const char* src =
        "ca = particle(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "dc = distance(ca[0], ca.xyz)\n";
    script::Program* p = compile_v2(src, alloc);
    script::Eval* e = script::eval_create(nullptr, p, eval_desc(), alloc);
    ASSERT_TRUE(e != nullptr);
    EXPECT_TRUE(script::eval_whole(e));
    EXPECT_EQ(script::fingerprint(p), script::eval_fingerprint(e));
    EXPECT_EQ(nf, script::eval_num_frames(e));
    md_bitset_t pending = script::eval_pending_frames(e, alloc);
    EXPECT_EQ(nf, (size_t)md_bitset_count(pending));
    md_bitset_free(&pending, alloc);

    // One call does every frame, whichever frames it is given
    EXPECT_TRUE(script::eval_frames(e, p, &peptide()->sys, MD_STR_LIT("run"), nullptr, 0));
    EXPECT_EQ(1.0f, script::eval_progress(e));
    pending = script::eval_pending_frames(e, alloc);
    EXPECT_EQ(0u, (size_t)md_bitset_count(pending));
    md_bitset_free(&pending, alloc);
    const md_attributes_t* attrs = script::attributes(e);
    ASSERT_TRUE(attrs != nullptr);
    const md_attribute_t* ad = md_attributes_find(attrs, MD_STR_LIT("script/d"));
    ASSERT_TRUE(ad != nullptr);
    EXPECT_EQ(nf, (size_t)ad->format.shape[0]);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/h/range")) != nullptr);
    // A population: its summary per frame, which the Timelines' bands read, and the names of its members
    const md_attribute_t* mean = md_attributes_find(attrs, MD_STR_LIT("script/dc/mean"));
    ASSERT_TRUE(mean != nullptr);
    EXPECT_EQ(nf, (size_t)mean->format.shape[0]);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/dc/variance")) != nullptr);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/dc/extent")) != nullptr);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/d/mean")) == nullptr);
    const md_attribute_t* labels = md_attributes_find(attrs, MD_STR_LIT("script/dc/labels"));
    ASSERT_TRUE(labels != nullptr);
    EXPECT_TRUE(md_str_eq(md_attribute_str(attrs, labels, 1), MD_STR_LIT("CA/ALA2")));
    ad = md_attributes_find(attrs, MD_STR_LIT("script/d"));       // the table grew: pointers into it move
    md_bitset_t held = script::eval_output_frames(e, MD_STR_LIT("d"), alloc);
    EXPECT_EQ(nf, (size_t)md_bitset_count(held));
    md_bitset_free(&held, alloc);

    // A period: the histogram over its frames, from the values kept per frame
    const md_bitset_t first_half = md_bitset_range(0, (uint32_t)(nf / 2));
    ASSERT_TRUE(script::aggregate_request(e, 7, MD_STR_LIT("h"), &first_half));
    EXPECT_FALSE(script::aggregate_request(e, 7, MD_STR_LIT("d"), &first_half));     // not a distribution
    script::AggregateInfo info;
    ASSERT_TRUE(script::aggregate_info(&info, e, 7, MD_STR_LIT("h")));
    EXPECT_EQ(10u, info.num_values);
    EXPECT_EQ((uint32_t)(nf / 2), info.num_frames);
    float values[10], weights[10];
    ASSERT_TRUE(script::aggregate_read(&info, values, weights, 10, e, 7, MD_STR_LIT("h")));
    double total = 0;
    size_t within = 0;
    for (int i = 0; i < 10; ++i) total += values[i];
    for (size_t f = 0; f < nf / 2; ++f) {
        const float v = ((const float*)ad->data)[f];
        within += v >= 0.0f && v < 20.0f;
    }
    EXPECT_EQ((double)within, total);
    const uint64_t v0 = info.version;
    const md_bitset_t all = md_bitset_range(0, (uint32_t)nf);
    ASSERT_TRUE(script::aggregate_request(e, 7, MD_STR_LIT("h"), &all));
    ASSERT_TRUE(script::aggregate_info(&info, e, 7, MD_STR_LIT("h")));
    EXPECT_NE(v0, info.version);
    script::aggregate_release(e, 7, MD_STR_LIT("h"));
    EXPECT_FALSE(script::aggregate_info(&info, e, 7, MD_STR_LIT("h")));

    // An edit: renamed, d is taken over and nothing is pending; its values are there at once
    script::Program* q = compile_v2("ca = particle(name='CA')\ndist = distance(ca[0], ca[1])\nh = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\n", alloc);
    script::Plan plans[4];
    md_bitset_t would = md_bitset_empty_set();
    ASSERT_TRUE(script::plan(plans, 4, &would, e, q, nf, alloc));
    EXPECT_EQ((int)script::Status::Held, (int)plans[0].status);
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("d"), plans[0].from));
    EXPECT_EQ(0u, (size_t)md_bitset_count(would));
    script::Eval* e2 = script::eval_create(e, q, eval_desc(), alloc);
    ASSERT_TRUE(e2 != nullptr);
    pending = script::eval_pending_frames(e2, alloc);
    EXPECT_EQ(0u, (size_t)md_bitset_count(pending));
    md_bitset_free(&pending, alloc);
    const md_attribute_t* a2 = md_attributes_find(script::attributes(e2), MD_STR_LIT("script/dist"));
    ASSERT_TRUE(a2 != nullptr);
    EXPECT_EQ(0, memcmp(ad->data, a2->data, nf * sizeof(float)));
    script::eval_free(e2);
    script::eval_free(e);
    script::program_free(q);
    script::program_free(p);
}

UTEST(script_v2, drawing) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca = particle(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "a  = angle(ca[0], ca[1], ca[2])\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(script::valid(p));
    script::DrawCtx ctx;
    ctx.prog = p;
    ctx.sys = &peptide()->sys;
    ctx.state = &peptide()->state;
    ctx.run = MD_STR_LIT("run");
    auto mark_of = [&](const char* text) {
        const size_t n = script::num_marks(p);
        const script::Mark* m = script::marks(p);
        const int beg = (int)(strstr(src, text) - src), end = beg + (int)strlen(text);
        for (size_t i = 0; i < n; ++i) if (m[i].range.beg == beg && m[i].range.end == end) return m[i].ref;
        return script::Ref{};
    };

    // A distance: a line between the two, labelled with the value
    script::Drawing dr;
    script::drawing_init(&dr, alloc);
    const script::Ref rd = mark_of("distance(ca[0], ca[1])");
    ASSERT_TRUE(script::ref_valid(p, rd));
    ASSERT_TRUE(script::draw_ref(&dr, rd, -1, ctx, script::DRAW_DEFAULT));
    EXPECT_EQ(2u, md_array_size(dr.lines));
    ASSERT_EQ(1u, md_array_size(dr.text));
    EXPECT_FALSE(md_unit_is_none(dr.text[0].unit));
    char buf[128];
    EXPECT_GT(script::value_text(buf, sizeof(buf), rd, -1, ctx), 0u);
    EXPECT_NEAR(dr.text[0].value, atof(buf), 1e-3);
    // Geometry alone: no text
    script::drawing_clear(&dr);
    ASSERT_TRUE(script::draw_ref(&dr, rd, -1, ctx, script::DRAW_GEOMETRY));
    EXPECT_EQ(0u, md_array_size(dr.text));
    EXPECT_EQ(2u, md_array_size(dr.lines));

    // An angle: a filled arc at the vertex
    script::drawing_clear(&dr);
    ASSERT_TRUE(script::draw_ref(&dr, mark_of("angle(ca[0], ca[1], ca[2])"), -1, ctx, script::DRAW_DEFAULT));
    EXPECT_GT(md_array_size(dr.triangles), 0u);
    EXPECT_GT(md_array_size(dr.lines), 2u);

    // A selection: its atoms; an expression of the language: what it picks
    script::drawing_clear(&dr);
    ASSERT_TRUE(script::draw_ref(&dr, mark_of("particle(name='CA')"), -1, ctx, script::DRAW_ATOMS));
    ASSERT_EQ(1u, md_array_size(dr.groups));
    const uint64_t nca = md_bitset_count(dr.atom_mask);
    EXPECT_GT(nca, 2u);
    script::drawing_clear(&dr);
    ASSERT_TRUE(script::draw_text(&dr, MD_STR_LIT("residue()"), ctx, script::DRAW_DEFAULT));
    EXPECT_EQ(md_system_residue_count(&peptide()->sys), md_array_size(dr.groups));
    script::drawing_clear(&dr);
    EXPECT_FALSE(script::draw_ref(&dr, script::Ref{1, 1}, -1, ctx, script::DRAW_DEFAULT));   // another compilation's
    script::program_free(p);
    md_arena_allocator_destroy(alloc);
}

UTEST(script_v2, filters) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const script::Language V2 = script::Language::V2;
    script::Program* p = compile_v2("ca = particle(name='CA')\nres = residue()\n", alloc);
    script::FilterCtx ctx;
    ctx.sys = &peptide()->sys;
    ctx.state = &peptide()->state;
    ctx.prog = p;
    ctx.run = MD_STR_LIT("run");

    md_bitset_t atoms = md_bitset_empty_set();
    script::FilterResult r = script::filter(V2, &atoms, MD_STR_LIT("protein"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_FALSE(r.is_dynamic);
    EXPECT_EQ((uint64_t)peptide()->sys.atom.count, md_bitset_count(atoms));
    // The program's selection, and the language's set algebra
    r = script::filter(V2, &atoms, MD_STR_LIT("ca"), ctx, alloc);
    ASSERT_TRUE(r.success);
    const uint64_t nca = md_bitset_count(atoms);
    EXPECT_GT(nca, 2u);
    r = script::filter(V2, &atoms, MD_STR_LIT("protein and not ca"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ((uint64_t)peptide()->sys.atom.count - nca, md_bitset_count(atoms));

    md_array(md_bitset_t) sets = nullptr;
    r = script::filter_sets(V2, &sets, MD_STR_LIT("residue()"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ(md_system_residue_count(&peptide()->sys), md_array_size(sets));
    // An array the program publishes, by name: all its sets, or one
    sets = nullptr;
    r = script::filter_sets(V2, &sets, MD_STR_LIT("res"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ(md_system_residue_count(&peptide()->sys), md_array_size(sets));
    r = script::filter(V2, &atoms, MD_STR_LIT("res[1]"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ((uint64_t)md_system_residue_atom_count(&peptide()->sys, 1), md_bitset_count(atoms));

    // Bonds: the pairs bond() gives, or the bonds among atoms; the same for one residue
    md_array(md_atom_pair_t) b1 = nullptr;
    md_array(md_atom_pair_t) b2 = nullptr;
    ASSERT_TRUE(script::filter_bonds(V2, &b1, MD_STR_LIT("bond(residue()[0])"), ctx, alloc).success);
    ASSERT_TRUE(script::filter_bonds(V2, &b2, MD_STR_LIT("residue()[0]"), ctx, alloc).success);
    EXPECT_GT(md_array_size(b1), 0u);
    EXPECT_EQ(md_array_size(b1), md_array_size(b2));

    // Dynamic, within a scope, and an error with its range in the expression
    r = script::filter(V2, &atoms, MD_STR_LIT("within(3.0, of=ca)"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_TRUE(r.is_dynamic);
    md_bitset_t scope = md_bitset_empty_set();
    ASSERT_TRUE(script::filter(V2, &scope, MD_STR_LIT("residue()[0]"), ctx, alloc).success);
    script::FilterCtx sc = ctx;
    sc.atoms = &scope;
    ASSERT_TRUE(script::filter(V2, &atoms, MD_STR_LIT("protein"), sc, alloc).success);
    EXPECT_EQ(md_bitset_count(scope), md_bitset_count(atoms));
    r = script::filter(V2, &atoms, MD_STR_LIT("particle(nme='CA')"), ctx, alloc);
    EXPECT_FALSE(r.success);
    EXPECT_EQ('E', r.error[0]);
    EXPECT_GE(r.error_range.beg, 0);
    EXPECT_LE(r.error_range.end, 19);

    script::Completions c = script::filter_complete(V2, MD_STR_LIT("particle(res"), 12, ctx, alloc);
    bool found = false;
    for (size_t i = 0; i < c.count; ++i) found |= md_str_eq(c.items[i].label, MD_STR_LIT("resname"));
    EXPECT_TRUE(found);
    script::program_free(p);
    md_arena_allocator_destroy(alloc);
}

UTEST(script_v2, runtime_errors) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    // A position past the CA particles from some frame on: compiles, and is missing there (R0101), a warning
    const char* src = "ca = particle(name='CA')\nx = distance(ca[0], ca[frame()])\nd = distance(ca[0], ca[1])\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(script::valid(p));
    script::Eval* e = script::eval_create(nullptr, p, eval_desc(), alloc);
    ASSERT_TRUE(e != nullptr);
    EXPECT_EQ(0u, script::eval_errors(e).num_warnings);
    EXPECT_TRUE(script::eval_frames(e, p, &peptide()->sys, MD_STR_LIT("run"), nullptr, 0));
    const script::EvalErrors er = script::eval_errors(e);
    EXPECT_NE(0u, er.version);
    EXPECT_EQ(0u, er.count);
    ASSERT_EQ(1u, er.num_warnings);
    EXPECT_EQ(0, strncmp(er.warnings[0].text.ptr, "R0101", 5));
    // About the position, within the call it is in
    const int at = (int)(strstr(src, "ca[frame()]") - src);
    EXPECT_GE(er.warnings[0].range.beg, at);
    EXPECT_LE(er.warnings[0].range.end, at + 11);
    EXPECT_TRUE(md_attributes_find(script::attributes(e), MD_STR_LIT("script/d")) != nullptr);
    script::eval_free(e);
    script::program_free(p);
}

// Who uses a program off its owner's thread: counted, from any thread, and nothing freed by the count
UTEST(script_v2, program_users) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    script::Program* p = compile_v2("ca = particle(name='CA')\n", alloc);
    EXPECT_EQ(0u, script::program_users(p));
    script::program_retain(p);
    script::program_retain(p);
    EXPECT_EQ(2u, script::program_users(p));
    script::program_release(p);
    EXPECT_EQ(1u, script::program_users(p));
    script::program_release(p);
    EXPECT_EQ(0u, script::program_users(p));
    EXPECT_TRUE(script::valid(p));
    script::program_retain(nullptr);
    EXPECT_EQ(0u, script::program_users(nullptr));
    script::program_free(p);
    md_arena_allocator_destroy(alloc);
}

// Volumes: sdf and density, published with what places them; over a period made again on a thread of their own, the
// UI asking until it is done; an sdf's structures drawn with the poses that place the volume on them
UTEST(script_v2, volumes) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t nf = peptide()->num_frames;
    const char* src =
        "s = sdf(residue[1:5], particle(name='O'), 6.0)\n"
        "d = density(particle(), box=((-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)), samples=0.25)\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(script::valid(p));
    EXPECT_EQ((uint32_t)(script::OUTPUT_VOLUME | script::OUTPUT_SDF), script::output_flags(p, MD_STR_LIT("s")));
    EXPECT_EQ((uint32_t)script::OUTPUT_VOLUME, script::output_flags(p, MD_STR_LIT("d")));
    script::Eval* e = script::eval_create(nullptr, p, eval_desc(), alloc);
    ASSERT_TRUE(e != nullptr);
    ASSERT_TRUE(script::eval_frames(e, p, &peptide()->sys, MD_STR_LIT("run"), nullptr, 0));
    const md_attributes_t* attrs = script::attributes(e);
    const md_attribute_t* as = md_attributes_find(attrs, MD_STR_LIT("script/s"));
    ASSERT_TRUE(as != nullptr);
    ASSERT_EQ(3u, as->format.rank);
    EXPECT_EQ(48u, as->format.shape[0]);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/s/origin")) != nullptr);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/s/basis")) != nullptr);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/d/bulk")) != nullptr);
    const size_t nv = 48 * 48 * 48;
    std::vector<float> full((const float*)as->data, (const float*)as->data + nv);

    // Over every frame, by a job: pending until it is done, then the volume published
    const md_bitset_t all = md_bitset_range(0, (uint32_t)nf);
    ASSERT_TRUE(script::aggregate_request(e, 9, MD_STR_LIT("s"), &all));
    script::AggregateInfo info;
    for (int i = 0; i < 2000; ++i) {
        ASSERT_TRUE(script::aggregate_info(&info, e, 9, MD_STR_LIT("s")));
        if (!info.num_pending) break;
        md_thread_sleep(5);
    }
    ASSERT_EQ(0u, info.num_pending);
    ASSERT_EQ(nv, info.num_values);
    std::vector<float> got(nv), w(nv);
    ASSERT_TRUE(script::aggregate_read(&info, got.data(), w.data(), nv, e, 9, MD_STR_LIT("s")));
    EXPECT_EQ(0, memcmp(full.data(), got.data(), nv * sizeof(float)));
    // Another period while one is made: asked again, and the evaluation freed with a job on its way
    const md_bitset_t half = md_bitset_range(0, (uint32_t)(nf / 2));
    ASSERT_TRUE(script::aggregate_request(e, 9, MD_STR_LIT("s"), &half));
    ASSERT_TRUE(script::aggregate_info(&info, e, 9, MD_STR_LIT("s")));
    EXPECT_GT(info.num_pending, 0u);

    // Drawn: its four structures, each with its pose
    script::DrawCtx ctx;
    ctx.prog = p;
    ctx.sys = &peptide()->sys;
    ctx.state = &peptide()->state;
    ctx.run = MD_STR_LIT("run");
    ctx.eval = e;
    script::Drawing dr;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    script::drawing_init(&dr, arena);
    ASSERT_TRUE(script::draw_ref(&dr, script::output_ref(p, MD_STR_LIT("s")), 0, ctx, script::DRAW_SDF));
    EXPECT_EQ(4u, md_array_size(dr.sdf.matrices));
    ASSERT_EQ(4u, md_array_size(dr.sdf.structures));
    EXPECT_EQ(10u, (unsigned)md_bitset_count(dr.sdf.structures[0]));
    md_arena_allocator_destroy(arena);

    script::eval_free(e);
    script::program_free(p);
}
