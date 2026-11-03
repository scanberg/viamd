#include "utest.h"
#include "test_executor.h"

#include <script_engine.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_edr.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_writer.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/* The acceptance gate of VIAL (script v2) against v1: pairs of scripts, one in each language, which are to compute the
 * same thing, evaluated through viamd's own interface (script_engine.h) on the same data, and compared output by
 * output. Every output v1 publishes has to be published by the VIAL script under the same name, with the same number
 * of values, and the same values within a tolerance, compared in v1's unit (v1 gives angles in radians, VIAL in
 * degrees, which the units convert). A periodic output is compared around its circle. What VIAL publishes beyond that
 * is not compared.
 *
 * Selections are compared by filters: each line of a .filters file is a v1 filter and a VIAL filter, and the two
 * have to pick the same atoms.
 *
 * The corpus is test/script_pairs: .pair files of two scripts, and .filters files of filter pairs. Its header says what
 * it runs on:
 *   # data: 1ALA-560ns.pdb                  a structure; a PDB of several models is its own trajectory
 *   # data: catalyst.gro catalyst.xtc       a structure and a trajectory (xtc, trr) or energies (edr)
 *   # tol: 1e-3 1e-4                        absolute and relative tolerance, in v1's unit (the default)
 * and, in a .pair file, the two scripts follow `=== v1` and `=== vial`. Comment lines (#) say what the pair checks and
 * why the two read as they do where the languages differ in what they mean. */

namespace {

namespace fs = std::filesystem;

struct Data {
    std::string     key;
    md_allocator_i* alloc = nullptr;
    md_system_t     sys = {};
    md_system_state_t state = {};
    size_t          num_frames = 0;
    bool            ok = false;
};

std::vector<Data*>& datasets() {
    static std::vector<Data*> d;
    return d;
}

bool ends_with(const std::string& s, const char* suffix) {
    const size_t n = strlen(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

// The structure and what is loaded along it, once per run of the tests
const Data* load(const std::string& files) {
    for (Data* d : datasets()) if (d->key == files) return d->ok ? d : nullptr;
    Data* d = new Data;
    d->key = files;
    datasets().push_back(d);
    std::istringstream in(files);
    std::vector<std::string> names;
    for (std::string f; in >> f;) names.push_back(std::string(VIAMD_TEST_DATA_DIR) + "/" + f);
    if (names.empty()) return nullptr;

    d->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    d->sys = {.alloc = d->alloc};
    d->state = {.alloc = d->alloc};
    const str_t top = {names[0].c_str(), names[0].size()};
    bool ok = ends_with(names[0], ".gro") ? md_gro_system_init_from_file(&d->sys, &d->state, top, NULL)
                                           : md_pdb_system_init_from_file(&d->sys, &d->state, top, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL);
    ok = ok && md_util_system_infer(&d->sys, &d->state, MD_UTIL_INFER_ALL, NULL);
    if (!ok) return nullptr;
    bool run = false;
    if (ends_with(names[0], ".pdb") && names.size() == 1) run = md_pdb_system_publish_run(&d->sys, top, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    for (size_t i = 1; i < names.size(); ++i) {
        const str_t f = {names[i].c_str(), names[i].size()};
        if (ends_with(names[i], ".xtc"))      run = md_xtc_system_publish_run(&d->sys, f, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
        else if (ends_with(names[i], ".trr")) run = md_trr_system_publish_run(&d->sys, f, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
        else if (ends_with(names[i], ".edr")) run = md_edr_system_supplement_from_file(&d->sys, f, STR_LIT("run"));
        else return nullptr;
    }
    if (run) {
        const md_attribute_t* t = md_attributes_find(&d->sys.attributes, STR_LIT("run/time"));
        d->num_frames = t ? (size_t)t->format.shape[0] : 0;
    }
    d->ok = true;
    return d;
}

struct Header {
    std::string data;
    double abs_tol = 1e-3, rel_tol = 1e-4;
};

void parse_header_line(Header* h, const std::string& line) {
    if (line.rfind("# data:", 0) == 0) h->data = line.substr(7);
    if (line.rfind("# tol:", 0) == 0) { std::istringstream in(line.substr(6)); in >> h->abs_tol >> h->rel_tol; }
}

std::string read_file(const fs::path& p) {
    std::ifstream in(p);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// Compiles and evaluates over every frame; NULL, said why, when either fails
script::Eval* run_script(script::Language lang, const std::string& src, const Data& d, script::Program** out_prog, md_allocator_i* alloc) {
    script::CompileDesc cd;
    cd.sys = &d.sys;
    cd.run = STR_LIT("run");
    script::Program* p = script::compile(lang, str_t{src.c_str(), src.size()}, cd, alloc);
    *out_prog = p;
    if (!script::valid(p)) {
        printf("  %s does not compile:\n", lang == script::Language::V1 ? "v1" : "VIAL");
        script::write_diagnostics(md_writer_stdout(), p);
        return nullptr;
    }
    script::EvalDesc ed;
    ed.num_frames = d.num_frames;
    ed.sys = &d.sys;
    ed.run = STR_LIT("run");
    static test_executor_t shuffled = test_executor_shuffled(2);   // one thread, two thread_idx, pieces out of order
    ed.executor = test_executor(&shuffled);
    script::Eval* e = script::eval_create(nullptr, p, ed, alloc);
    if (!e) return nullptr;
    if (d.num_frames) {
        std::vector<uint32_t> frames(d.num_frames);
        for (size_t i = 0; i < frames.size(); ++i) frames[i] = (uint32_t)i;
        script::eval_frames(e, p, &d.sys, STR_LIT("run"), frames.data(), frames.size());
    }
    return e;
}

const md_attribute_t* output_attr(const md_attributes_t* attrs, const std::string& name, const char* leaf = "") {
    const std::string path = "script/" + name + leaf;
    return md_attributes_find(attrs, str_t{path.c_str(), path.size()});
}

std::vector<double> values_of(const md_attribute_t* a, md_unit_t unit) {
    size_t n = 1;
    for (uint32_t i = 0; i < a->format.rank; ++i) n *= a->format.shape[i];
    n *= a->format.components ? a->format.components : 1;
    std::vector<double> v(n);
    const size_t got = md_attribute_extract_f64(v.data(), n, a, md_attribute_slice_all(), unit);
    v.resize(got);
    return v;
}

// Every output of v1's, in VIAL's, alike. Says what differs.
bool compare_outputs(const std::string& pair, const script::Program* p1, const script::Eval* e1, const script::Program* p2,
                     const script::Eval* e2, const Header& h) {
    const md_attributes_t* a1 = script::attributes(e1);
    const md_attributes_t* a2 = script::attributes(e2);
    const size_t n1 = script::num_outputs(p1);
    const str_t* names1 = script::output_names(p1);
    bool ok = true;
    size_t compared = 0;
    for (size_t i = 0; i < n1; ++i) {
        const std::string name(names1[i].ptr, names1[i].len);
        const md_attribute_t* x = output_attr(a1, name);
        const md_attribute_t* y = output_attr(a2, name);
        if (!x) continue;   // v1 names it, publishes nothing (a selection)
        if (!y) {
            printf("  %s: VIAL does not publish %s\n", pair.c_str(), name.c_str());
            ok = false;
            continue;
        }
        // In v1's unit; where either has none (v1 drops the unit of some arithmetic), each in its own
        const bool unitless = md_unit_is_none(x->unit) || md_unit_is_none(y->unit);
        const md_unit_t unit = x->unit;
        const std::vector<double> u = values_of(x, x->unit), v = values_of(y, unitless ? y->unit : unit);
        ++compared;
        if (u.size() != v.size()) {
            printf("  %s: %s has %zu values in v1 and %zu in VIAL\n", pair.c_str(), name.c_str(), u.size(), v.size());
            ok = false;
            continue;
        }
        // Around its circle where it is periodic
        double period = 0;
        if (const md_attribute_t* per = output_attr(a2, name, "/period")) {
            const std::vector<double> pr = values_of(per, unit);
            if (pr.size() == 2) period = pr[1] - pr[0];
        }
        size_t bad = 0;
        double worst = 0;
        size_t worst_i = 0;
        for (size_t k = 0; k < u.size(); ++k) {
            if (std::isnan(u[k]) && std::isnan(v[k])) continue;
            double d = fabs(u[k] - v[k]);
            if (period > 0) d = fmin(d, fabs(period - fmod(d, period)));
            const double tol = fmax(h.abs_tol, h.rel_tol * fabs(u[k]));
            if (!(d <= tol)) {
                if (!bad || d > worst) { worst = d; worst_i = k; }
                ++bad;
            }
        }
        if (bad) {
            printf("  %s: %s differs in %zu of %zu values; the most at %zu: v1 %.6g, VIAL %.6g\n", pair.c_str(), name.c_str(), bad,
                   u.size(), worst_i, u[worst_i], v[worst_i]);
            ok = false;
        }
    }
    if (!compared) {
        printf("  %s: v1 publishes nothing to compare\n", pair.c_str());
        ok = false;
    }
    return ok;
}

fs::path pairs_dir() { return fs::path(VIAMD_SCRIPT_PAIRS_DIR); }

std::vector<fs::path> files_of(const char* ext) {
    std::vector<fs::path> out;
    for (const auto& e : fs::directory_iterator(pairs_dir())) if (e.path().extension() == ext) out.push_back(e.path());
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

UTEST(script_pairs, outputs_agree) {
    const std::vector<fs::path> files = files_of(".pair");
    ASSERT_GT(files.size(), 0u);
    size_t passed = 0;
    for (const fs::path& f : files) {
        Header h;
        std::string v1, vial, *cur = nullptr;
        std::istringstream in(read_file(f));
        for (std::string line; std::getline(in, line);) {
            if (line == "=== v1") { cur = &v1; continue; }
            if (line == "=== vial") { cur = &vial; continue; }
            if (!cur) { parse_header_line(&h, line); continue; }
            *cur += line + "\n";
        }
        const std::string name = f.stem().string();
        const Data* d = load(h.data);
        if (!d) { printf("  %s: data %s does not load\n", name.c_str(), h.data.c_str()); EXPECT_TRUE(false); continue; }
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
        script::Program* p1 = nullptr;
        script::Program* p2 = nullptr;
        script::Eval* e1 = run_script(script::Language::V1, v1, *d, &p1, alloc);
        script::Eval* e2 = run_script(script::Language::V2, vial, *d, &p2, alloc);
        const bool ok = e1 && e2 && compare_outputs(name, p1, e1, p2, e2, h);
        if (!ok) printf("  %s: FAILS\n", name.c_str());
        passed += ok;
        EXPECT_TRUE(ok);
        script::eval_free(e1);
        script::eval_free(e2);
        script::program_free(p1);
        script::program_free(p2);
        md_arena_allocator_destroy(alloc);
    }
    printf("  script pairs: %zu of %zu agree\n", passed, files.size());
}

UTEST(script_pairs, filters_agree) {
    const std::vector<fs::path> files = files_of(".filters");
    ASSERT_GT(files.size(), 0u);
    size_t lines = 0, agree = 0;
    for (const fs::path& f : files) {
        Header h;
        std::vector<std::pair<std::string, std::string>> pairs;
        std::istringstream in(read_file(f));
        for (std::string line; std::getline(in, line);) {
            if (line.empty()) continue;
            if (line[0] == '#') { parse_header_line(&h, line); continue; }
            const size_t sep = line.find("||");
            if (sep == std::string::npos) continue;
            auto trim = [](std::string s) {
                s.erase(0, s.find_first_not_of(" \t"));
                s.erase(s.find_last_not_of(" \t") + 1);
                return s;
            };
            pairs.emplace_back(trim(line.substr(0, sep)), trim(line.substr(sep + 2)));
        }
        const Data* d = load(h.data);
        if (!d) { printf("  %s: data %s does not load\n", f.stem().string().c_str(), h.data.c_str()); EXPECT_TRUE(false); continue; }
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
        script::FilterCtx ctx;
        ctx.sys = &d->sys;
        ctx.state = &d->state;
        ctx.run = STR_LIT("run");
        for (const auto& pr : pairs) {
            ++lines;
            md_bitset_t a = md_bitset_empty_set(), b = md_bitset_empty_set();
            const script::FilterResult r1 = script::filter(script::Language::V1, &a, str_t{pr.first.c_str(), pr.first.size()}, ctx, alloc);
            const script::FilterResult r2 = script::filter(script::Language::V2, &b, str_t{pr.second.c_str(), pr.second.size()}, ctx, alloc);
            bool ok = r1.success && r2.success && md_bitset_equal(a, b);
            if (!r1.success) printf("  %s: v1 '%s' fails: %s\n", f.stem().string().c_str(), pr.first.c_str(), r1.error);
            if (!r2.success) printf("  %s: VIAL '%s' fails: %s\n", f.stem().string().c_str(), pr.second.c_str(), r2.error);
            if (r1.success && r2.success && !ok) {
                printf("  %s: '%s' picks %llu atoms, '%s' %llu\n", f.stem().string().c_str(), pr.first.c_str(),
                       (unsigned long long)md_bitset_count(a), pr.second.c_str(), (unsigned long long)md_bitset_count(b));
            }
            agree += ok;
            EXPECT_TRUE(ok);
        }
        md_arena_allocator_destroy(alloc);
    }
    printf("  filter pairs: %zu of %zu agree\n", agree, lines);
}
