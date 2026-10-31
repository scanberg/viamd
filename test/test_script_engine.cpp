#include "utest.h"

#include <script_engine.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_writer.h>

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

namespace {
// The frames one after another: enough to show that what the evaluation spreads over threads is all done
void serial_for(void*, size_t count, void (*fn)(void* arg, size_t i), void* arg) {
    for (size_t i = 0; i < count; ++i) fn(arg, i);
}

script::EvalDesc eval_desc() {
    script::EvalDesc d;
    d.num_frames = peptide()->num_frames;
    d.sys = &peptide()->sys;
    d.run = STR_LIT("run");
    d.num_threads = 3;
    d.parallel_for = serial_for;
    return d;
}
}  // namespace

UTEST(script_v2, evaluate) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t nf = peptide()->num_frames;
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n";
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
    EXPECT_TRUE(script::eval_frames(e, p, &peptide()->sys, STR_LIT("run"), nullptr, 0));
    EXPECT_EQ(1.0f, script::eval_progress(e));
    pending = script::eval_pending_frames(e, alloc);
    EXPECT_EQ(0u, (size_t)md_bitset_count(pending));
    md_bitset_free(&pending, alloc);
    const md_attributes_t* attrs = script::attributes(e);
    ASSERT_TRUE(attrs != nullptr);
    const md_attribute_t* ad = md_attributes_find(attrs, STR_LIT("script/d"));
    ASSERT_TRUE(ad != nullptr);
    EXPECT_EQ(nf, (size_t)ad->format.shape[0]);
    EXPECT_TRUE(md_attributes_find(attrs, STR_LIT("script/h/range")) != nullptr);
    md_bitset_t held = script::eval_output_frames(e, STR_LIT("d"), alloc);
    EXPECT_EQ(nf, (size_t)md_bitset_count(held));
    md_bitset_free(&held, alloc);

    // A period: the histogram over its frames, from the values kept per frame
    const md_bitset_t first_half = md_bitset_range(0, (uint32_t)(nf / 2));
    ASSERT_TRUE(script::aggregate_request(e, 7, STR_LIT("h"), &first_half));
    EXPECT_FALSE(script::aggregate_request(e, 7, STR_LIT("d"), &first_half));     // not a distribution
    script::AggregateInfo info;
    ASSERT_TRUE(script::aggregate_info(&info, e, 7, STR_LIT("h")));
    EXPECT_EQ(10u, info.num_values);
    EXPECT_EQ((uint32_t)(nf / 2), info.num_frames);
    float values[10], weights[10];
    ASSERT_TRUE(script::aggregate_read(&info, values, weights, 10, e, 7, STR_LIT("h")));
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
    ASSERT_TRUE(script::aggregate_request(e, 7, STR_LIT("h"), &all));
    ASSERT_TRUE(script::aggregate_info(&info, e, 7, STR_LIT("h")));
    EXPECT_NE(v0, info.version);
    script::aggregate_release(e, 7, STR_LIT("h"));
    EXPECT_FALSE(script::aggregate_info(&info, e, 7, STR_LIT("h")));

    // An edit: renamed, d is taken over and nothing is pending; its values are there at once
    script::Program* q = compile_v2("ca = particles(name='CA')\ndist = distance(ca[0], ca[1])\nh = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\n", alloc);
    script::Plan plans[4];
    md_bitset_t would = md_bitset_empty_set();
    ASSERT_TRUE(script::plan(plans, 4, &would, e, q, nf, alloc));
    EXPECT_EQ((int)script::Status::Held, (int)plans[0].status);
    EXPECT_TRUE(str_eq(STR_LIT("d"), plans[0].from));
    EXPECT_EQ(0u, (size_t)md_bitset_count(would));
    script::Eval* e2 = script::eval_create(e, q, eval_desc(), alloc);
    ASSERT_TRUE(e2 != nullptr);
    pending = script::eval_pending_frames(e2, alloc);
    EXPECT_EQ(0u, (size_t)md_bitset_count(pending));
    md_bitset_free(&pending, alloc);
    const md_attribute_t* a2 = md_attributes_find(script::attributes(e2), STR_LIT("script/dist"));
    ASSERT_TRUE(a2 != nullptr);
    EXPECT_EQ(0, memcmp(ad->data, a2->data, nf * sizeof(float)));
    script::eval_free(e2);
    script::eval_free(e);
    script::program_free(q);
    script::program_free(p);
}

UTEST(script_v2, drawing) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "a  = angle(ca[0], ca[1], ca[2])\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(script::valid(p));
    script::DrawCtx ctx;
    ctx.prog = p;
    ctx.sys = &peptide()->sys;
    ctx.state = &peptide()->state;
    ctx.run = STR_LIT("run");
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
    ASSERT_TRUE(script::draw_ref(&dr, mark_of("particles(name='CA')"), -1, ctx, script::DRAW_ATOMS));
    ASSERT_EQ(1u, md_array_size(dr.groups));
    const uint64_t nca = md_bitset_count(dr.atom_mask);
    EXPECT_GT(nca, 2u);
    script::drawing_clear(&dr);
    ASSERT_TRUE(script::draw_text(&dr, STR_LIT("residues()"), ctx, script::DRAW_DEFAULT));
    EXPECT_EQ(md_system_residue_count(&peptide()->sys), md_array_size(dr.groups));
    script::drawing_clear(&dr);
    EXPECT_FALSE(script::draw_ref(&dr, script::Ref{1, 1}, -1, ctx, script::DRAW_DEFAULT));   // another compilation's
    script::program_free(p);
    md_arena_allocator_destroy(alloc);
}

UTEST(script_v2, filters) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const script::Language V2 = script::Language::V2;
    script::Program* p = compile_v2("ca = particles(name='CA')\n", alloc);
    script::FilterCtx ctx;
    ctx.sys = &peptide()->sys;
    ctx.state = &peptide()->state;
    ctx.prog = p;
    ctx.run = STR_LIT("run");

    md_bitset_t atoms = md_bitset_empty_set();
    script::FilterResult r = script::filter(V2, &atoms, STR_LIT("protein"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_FALSE(r.is_dynamic);
    EXPECT_EQ((uint64_t)peptide()->sys.atom.count, md_bitset_count(atoms));
    // The program's selection, and the language's set algebra
    r = script::filter(V2, &atoms, STR_LIT("ca"), ctx, alloc);
    ASSERT_TRUE(r.success);
    const uint64_t nca = md_bitset_count(atoms);
    EXPECT_GT(nca, 2u);
    r = script::filter(V2, &atoms, STR_LIT("protein and not ca"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ((uint64_t)peptide()->sys.atom.count - nca, md_bitset_count(atoms));

    md_array(md_bitset_t) sets = nullptr;
    r = script::filter_sets(V2, &sets, STR_LIT("residues()"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_EQ(md_system_residue_count(&peptide()->sys), md_array_size(sets));

    // Bonds: the pairs bonds() gives, or the bonds among atoms; the same for one residue
    md_array(md_atom_pair_t) b1 = nullptr;
    md_array(md_atom_pair_t) b2 = nullptr;
    ASSERT_TRUE(script::filter_bonds(V2, &b1, STR_LIT("bonds(residues()[0])"), ctx, alloc).success);
    ASSERT_TRUE(script::filter_bonds(V2, &b2, STR_LIT("residues()[0]"), ctx, alloc).success);
    EXPECT_GT(md_array_size(b1), 0u);
    EXPECT_EQ(md_array_size(b1), md_array_size(b2));

    // Dynamic, within a scope, and an error with its range in the expression
    r = script::filter(V2, &atoms, STR_LIT("within(3.0, of=ca)"), ctx, alloc);
    ASSERT_TRUE(r.success);
    EXPECT_TRUE(r.is_dynamic);
    md_bitset_t scope = md_bitset_empty_set();
    ASSERT_TRUE(script::filter(V2, &scope, STR_LIT("residues()[0]"), ctx, alloc).success);
    script::FilterCtx sc = ctx;
    sc.atoms = &scope;
    ASSERT_TRUE(script::filter(V2, &atoms, STR_LIT("protein"), sc, alloc).success);
    EXPECT_EQ(md_bitset_count(scope), md_bitset_count(atoms));
    r = script::filter(V2, &atoms, STR_LIT("particles(nme='CA')"), ctx, alloc);
    EXPECT_FALSE(r.success);
    EXPECT_EQ('E', r.error[0]);
    EXPECT_GE(r.error_range.beg, 0);
    EXPECT_LE(r.error_range.end, 19);

    script::Completions c = script::filter_complete(V2, STR_LIT("particles(res"), 13, ctx, alloc);
    bool found = false;
    for (size_t i = 0; i < c.count; ++i) found |= str_eq(c.items[i].label, STR_LIT("resname"));
    EXPECT_TRUE(found);
    script::program_free(p);
    md_arena_allocator_destroy(alloc);
}

UTEST(script_v2, runtime_errors) {
    ASSERT_TRUE(peptide() != nullptr);
    md_allocator_i* alloc = md_get_heap_allocator();
    // A position past the CA particles from some frame on: compiles, and is missing there (R0101), a warning
    const char* src = "ca = particles(name='CA')\nx = distance(ca[0], ca[frame()])\nd = distance(ca[0], ca[1])\n";
    script::Program* p = compile_v2(src, alloc);
    ASSERT_TRUE(script::valid(p));
    script::Eval* e = script::eval_create(nullptr, p, eval_desc(), alloc);
    ASSERT_TRUE(e != nullptr);
    EXPECT_EQ(0u, script::eval_errors(e).num_warnings);
    EXPECT_TRUE(script::eval_frames(e, p, &peptide()->sys, STR_LIT("run"), nullptr, 0));
    const script::EvalErrors er = script::eval_errors(e);
    EXPECT_NE(0u, er.version);
    EXPECT_EQ(0u, er.count);
    ASSERT_EQ(1u, er.num_warnings);
    EXPECT_EQ(0, strncmp(er.warnings[0].text.ptr, "R0101", 5));
    // About the position, within the call it is in
    const int at = (int)(strstr(src, "ca[frame()]") - src);
    EXPECT_GE(er.warnings[0].range.beg, at);
    EXPECT_LE(er.warnings[0].range.end, at + 11);
    EXPECT_TRUE(md_attributes_find(script::attributes(e), STR_LIT("script/d")) != nullptr);
    script::eval_free(e);
    script::program_free(p);
}
