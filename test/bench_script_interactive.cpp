#include "utest.h"

#include <script_engine.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <chrono>
#include <string>
#include <vector>

/* What the script costs while it is being written, v1 against VIAL (v2), through viamd's own interface: what runs as
 * one types (compile, completion), as one hovers (drawing an expression and its value text), and as a frame changes
 * (filters). Run on demand:
 *
 *   VIAMD_BENCH=1 viamd_test --filter=script_bench.*
 *
 * VIAMD_BENCH_TOP and VIAMD_BENCH_TRAJ name another structure and trajectory, VIAMD_BENCH_REPS the repetitions (the
 * median is reported). Every line starts with the version below, bumped with every change to what is measured or how,
 * so a measurement of a stale build is told apart. */

#define SCRIPT_BENCH_VERSION "script-bench r4"

namespace {

using clock_type = std::chrono::steady_clock;

double ms_since(clock_type::time_point t0) {
    return std::chrono::duration<double, std::milli>(clock_type::now() - t0).count();
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

bool ends_with(const char* s, const char* suffix) {
    const size_t a = strlen(s), b = strlen(suffix);
    return a >= b && strcmp(s + a - b, suffix) == 0;
}

struct Bench {
    md_allocator_i*   alloc = nullptr;
    md_system_t       sys = {};
    md_system_state_t state = {};
    size_t            num_frames = 0;
};

bool load(Bench* b, const char* top, const char* traj) {
    b->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    b->sys = {.alloc = b->alloc};
    b->state = {.alloc = b->alloc};
    const str_t t = str_from_cstr(top);
    bool ok = ends_with(top, ".gro") ? md_gro_system_init_from_file(&b->sys, &b->state, t)
                                      : md_pdb_system_init_from_file(&b->sys, &b->state, t, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE);
    ok = ok && md_util_system_infer(&b->sys, &b->state, MD_UTIL_INFER_ALL);
    if (!ok) return false;
    const str_t r = str_from_cstr(traj ? traj : top);
    bool run = false;
    if (!traj) run = md_pdb_system_publish_run(&b->sys, r, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    else if (ends_with(traj, ".xtc")) run = md_xtc_system_publish_run(&b->sys, r, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    else if (ends_with(traj, ".trr")) run = md_trr_system_publish_run(&b->sys, r, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    if (run) {
        const md_attribute_t* tm = md_attributes_find(&b->sys.attributes, STR_LIT("run/time"));
        b->num_frames = tm ? (size_t)tm->format.shape[0] : 0;
    }
    return true;
}

void serial_for(void*, size_t count, void (*fn)(void* arg, size_t i), void* arg) {
    for (size_t i = 0; i < count; ++i) fn(arg, i);
}

// The same analysis in each language, with the filters a user of it would write
struct Case {
    const char* name;
    const char* v1;
    const char* v2;
    const char* filters_v1[4];
    const char* filters_v2[4];
};

// Written for any system: its first residues, its waters and its oxygens
const Case cases[] = {
    {"distances",
     "d1 = distance(com(residue(1)), com(residue(2)));\n"
     "dm = distance_min(residue(1), residue(2));\n"
     "a = angle(1, 2, 3);\n",
     "d1 = distance(residue[0].com, residue[1].com)\n"
     "dm = distance(residue[0], residue[1], between='closest')\n"
     "a = angle(particle[0], particle[1], particle[2])\n",
     {"residue(1:3)", "element(\"O\") and water()", nullptr},
     {"particle(residue[0:3])", "particle(element='O') and water", nullptr}},
    {"neighbourhood",
     "nw = count(within(3.5, residue(1)) and water(), \"residue\");\n"
     "ns = count(within(3.0:5.0, residue(1)) and element(\"O\"));\n",
     "r0 = residue[0]\n"
     "nw = len(residue(within(3.5, of=r0) and water and not r0))\n"
     "ns = count(within(5.0, of=r0) and not within(3.0, of=r0) and particle(element='O') and not r0)\n",
     {"within(3.5, residue(1)) and water()", "residue(within(5.0, residue(1)))", nullptr},
     {"within(3.5, of=residue[0]) and water and not residue[0]",
      "particle(residue(within(5.0, of=residue[0]) and not residue[0]))", nullptr}},
    {"per residue",
     "c = distance(1, 2) in residue(1:3);\n"
     "s = count(within_z(20:25));\n",
     "c = [distance(particle(r)[0], particle(r)[1]) for r in residue[0:3]]\n"
     "s = count(particle[particle.z in 20.0..25.0])\n",
     {"within_z(20:25)", "residue(2:3)", nullptr},
     {"particle[particle.z in 20.0..25.0]", "particle(residue[1:3])", nullptr}},
};

struct Times {
    double compile = 0, hover = 0, hover_max = 0, complete = 0, filter = 0;
    size_t marks = 0, cursors = 0, filters = 0;
    std::string slowest;    // the mark that hovers slowest
    bool ok = false;
};

Times measure(script::Language lang, const Case& c, const Bench& b, int reps) {
    Times t;
    const char* src = lang == script::Language::V1 ? c.v1 : c.v2;
    const char* const* filters = lang == script::Language::V1 ? c.filters_v1 : c.filters_v2;
    const str_t s = str_from_cstr(src);
    script::CompileDesc cd;
    cd.sys = &b.sys;
    cd.run = STR_LIT("run");

    // Compile: once per edit
    std::vector<double> times;
    script::Program* p = nullptr;
    md_allocator_i* palloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    for (int r = 0; r < reps; ++r) {
        if (p) script::program_free(p);
        const auto t0 = clock_type::now();
        p = script::compile(lang, s, cd, palloc);
        times.push_back(ms_since(t0));
    }
    t.compile = median(times);
    if (!script::valid(p)) {
        printf("%s: %s does not compile in %s\n", SCRIPT_BENCH_VERSION, c.name, lang == script::Language::V1 ? "v1" : "v2");
        return t;
    }

    // An evaluation, which hover reads values over the trajectory from
    script::EvalDesc ed;
    ed.num_frames = b.num_frames;
    ed.sys = &b.sys;
    ed.run = STR_LIT("run");
    ed.num_threads = 1;
    ed.parallel_for = serial_for;
    script::Eval* e = script::eval_create(nullptr, p, ed, palloc);
    if (e && b.num_frames) {
        std::vector<uint32_t> frames(b.num_frames);
        for (size_t i = 0; i < frames.size(); ++i) frames[i] = (uint32_t)i;
        script::eval_frames(e, p, &b.sys, STR_LIT("run"), frames.data(), frames.size());
    }

    // Hover: every mark drawn and its value written, as moving the mouse over it does
    script::DrawCtx ctx;
    ctx.prog = p;
    ctx.sys = &b.sys;
    ctx.state = &b.state;
    ctx.eval = e;
    ctx.run = STR_LIT("run");
    const size_t nm = script::num_marks(p);
    const script::Mark* marks = script::marks(p);
    t.marks = nm;
    std::vector<double> per_mark(nm, 0.0);
    md_allocator_i* dalloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    for (size_t m = 0; m < nm; ++m) {
        std::vector<double> mt;
        for (int r = 0; r < reps; ++r) {
            script::Drawing d;
            script::drawing_init(&d, dalloc);
            char buf[256];
            const auto t0 = clock_type::now();
            script::draw_ref(&d, marks[m].ref, -1, ctx, script::DRAW_DEFAULT);
            script::value_text(buf, sizeof(buf), marks[m].ref, -1, ctx);
            mt.push_back(ms_since(t0));
            script::drawing_clear(&d);
            md_arena_allocator_reset(dalloc);
        }
        per_mark[m] = median(mt);
    }
    for (size_t m = 0; m < nm; ++m) {
        t.hover += per_mark[m];
        if (per_mark[m] > t.hover_max) {
            t.hover_max = per_mark[m];
            t.slowest = std::string(src + marks[m].range.beg, (size_t)(marks[m].range.end - marks[m].range.beg));
        }
    }
    if (nm) t.hover /= (double)nm;
    md_arena_allocator_destroy(dalloc);

    // Completion: at the end of every line, and inside the first call, as typing does
    std::vector<int> cursors;
    for (int i = 0; src[i]; ++i) if (src[i] == '\n') cursors.push_back(i);
    if (const char* paren = strchr(src, '(')) cursors.push_back((int)(paren - src) + 1);
    t.cursors = cursors.size();
    md_allocator_i* calloc_ = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    double sum = 0;
    for (int cur : cursors) {
        std::vector<double> ct;
        for (int r = 0; r < reps; ++r) {
            const auto t0 = clock_type::now();
            script::complete(lang, s, cur, &b.sys, calloc_);
            ct.push_back(ms_since(t0));
            md_arena_allocator_reset(calloc_);
        }
        sum += median(ct);
    }
    t.complete = cursors.empty() ? 0 : sum / (double)cursors.size();
    md_arena_allocator_destroy(calloc_);

    // Filters: evaluated again whenever the frame changes, for every representation that has one
    script::FilterCtx fc;
    fc.sys = &b.sys;
    fc.state = &b.state;
    fc.prog = p;
    fc.run = STR_LIT("run");
    md_allocator_i* falloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    sum = 0;
    for (int i = 0; i < 4 && filters[i]; ++i) {
        std::vector<double> ft;
        for (int r = 0; r < reps; ++r) {
            md_bitset_t out = md_bitset_empty_set();
            const auto t0 = clock_type::now();
            const script::FilterResult fr = script::filter(lang, &out, str_from_cstr(filters[i]), fc, falloc);
            ft.push_back(ms_since(t0));
            if (!fr.success) printf("%s: %s: filter '%s' fails: %s\n", SCRIPT_BENCH_VERSION, c.name, filters[i], fr.error);
            md_arena_allocator_reset(falloc);
        }
        sum += median(ft);
        t.filters++;
    }
    t.filter = t.filters ? sum / (double)t.filters : 0;
    md_arena_allocator_destroy(falloc);

    script::eval_free(e);
    script::program_free(p);
    md_arena_allocator_destroy(palloc);
    t.ok = true;
    return t;
}

}  // namespace

UTEST(script_bench, interactive) {
    if (!getenv("VIAMD_BENCH")) return;
    const char* top = getenv("VIAMD_BENCH_TOP");
    const char* traj = getenv("VIAMD_BENCH_TRAJ");
    const std::string def_top = std::string(VIAMD_TEST_DATA_DIR) + "/catalyst.gro";
    const std::string def_traj = std::string(VIAMD_TEST_DATA_DIR) + "/catalyst.xtc";
    if (!top) { top = def_top.c_str(); traj = def_traj.c_str(); }
    const int reps = getenv("VIAMD_BENCH_REPS") ? std::max(1, atoi(getenv("VIAMD_BENCH_REPS"))) : 15;
    Bench b;
    ASSERT_TRUE(load(&b, top, traj));
    printf("%s: %s, %zu particles, %zu frames; median of %d; ms per call\n", SCRIPT_BENCH_VERSION, top, b.sys.atom.count, b.num_frames, reps);
    printf("%s: %-14s %-4s %9s %9s %9s %9s %9s   %s\n", SCRIPT_BENCH_VERSION, "case", "lang", "compile", "hover", "hover max", "complete", "filter",
           "(marks, cursors, filters)");
    for (const Case& c : cases) {
        for (int l = 0; l < 2; ++l) {
            const script::Language lang = l ? script::Language::V2 : script::Language::V1;
            const Times t = measure(lang, c, b, reps);
            if (!t.ok) continue;
            printf("%s: %-14s %-4s %9.3f %9.3f %9.3f %9.3f %9.3f   (%zu, %zu, %zu) slowest hover: %.48s\n", SCRIPT_BENCH_VERSION, c.name,
                   l ? "v2" : "v1", t.compile, t.hover, t.hover_max, t.complete, t.filter, t.marks, t.cursors, t.filters, t.slowest.c_str());
        }
    }
    md_arena_allocator_destroy(b.alloc);
}
