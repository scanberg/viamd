// The lines of the Script menu, written in VIAL (script v2). What they offer is what script_suggestions.h says; how
// they are written is VIAL's: selections are particles(...) with keywords, which all have to hold, and the predefined
// sets (protein, water, ...); a particle by its position is particles[i]; and the same measurement in every residue of a
// name is one call over arrays of sets, angle(residues('ALA') and particles(name='N'), ...).
//
// Every line is checked as v1's are: a selection is evaluated (script::filter_sets) and has to pick exactly the selected
// particles, each argument of a measurement has to pick its particle alone, and the line has to compile.

#include "script_suggestions.h"

#include <md_system.h>
#include <md_util.h>
#include <script_engine.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>

#include <stdlib.h>
#include <string.h>

namespace {

constexpr size_t MAX_EXPR_LEN = 240;            // as v1: past this a line lists rather than describes
constexpr size_t MAX_LISTED_EXPR_LEN = 4096;
constexpr size_t MAX_SUGGESTIONS = 6;
constexpr size_t MAX_DESCRIPTIONS = 3;          // descriptions by keywords offered, the shortest
constexpr size_t MAX_TERMS = 12;                // terms joined by `or` in a listing that is still offered
constexpr size_t MAX_RESID_LIST = 8;            // resid=[...] of single ids, beyond which it is not a description

// The keywords a description is made of, in the order they are written: where, then what
enum Key { Key_Chain, Key_Resid, Key_Resname, Key_Name, Key_Element, Key_Count };
constexpr size_t key_max[Key_Count] = {4, 0, 3, 4, 3};     // the most values each takes; resid by its own rule

const char* const predefined[] = {"protein", "nucleic", "water", "ion", "lipid", "backbone", "sidechain"};
constexpr size_t num_predefined = sizeof(predefined) / sizeof(predefined[0]);

str_t take(md_str_writer_t* sb, md_allocator_i* alloc) {
    str_t s = str_copy(md_str_writer_to_str(sb), alloc);
    md_str_writer_reset(sb);
    return s;
}

// A string is written in single quotes, or double where it holds a single one; one which holds both, a backslash,
// something that cannot be typed or a glob character (a name is matched as a glob) is not written
bool quotable(str_t s) {
    if (str_empty(s)) return false;
    bool dq = false, sq = false;
    for (size_t i = 0; i < s.len; ++i) {
        const unsigned char ch = (unsigned char)s.ptr[i];
        if (ch < 0x20 || ch == 0x7F || ch == '\\' || ch == '*' || ch == '?') return false;
        dq |= ch == '"';
        sq |= ch == '\'';
    }
    return !(dq && sq);
}

void write_quoted(md_str_writer_t* sb, str_t s) {
    size_t loc;
    const char q = str_find_char(&loc, s, '\'') ? '"' : '\'';
    *sb += q;
    *sb += s;
    *sb += q;
}

// 'A', or ['A', 'B']
void write_strings(md_str_writer_t* sb, const str_t* s, size_t n) {
    if (n > 1) *sb += '[';
    for (size_t i = 0; i < n; ++i) {
        if (i) *sb += ", ";
        write_quoted(sb, s[i]);
    }
    if (n > 1) *sb += ']';
}

int cmp_int32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

size_t sort_unique(int32_t* v, size_t n) {
    if (n < 2) return n;
    qsort(v, n, sizeof(int32_t), cmp_int32);
    size_t k = 1;
    for (size_t i = 1; i < n; ++i) if (v[i] != v[k - 1]) v[k++] = v[i];
    return k;
}

// Sorted unique ids as one filter value: 12, 10..20, or [3, 7, 9]; a list of ranges is not a value, so false when
// they are neither one run nor few enough to list
bool write_ids(md_str_writer_t* sb, const int32_t* v, size_t n) {
    if (n == 0) return false;
    if (n == 1) { md_str_writer_fmt(sb, "%d", v[0]); return true; }
    if ((int64_t)v[n - 1] - v[0] == (int64_t)n - 1) { md_str_writer_fmt(sb, "%d..%d", v[0], v[n - 1]); return true; }
    if (n > MAX_RESID_LIST) return false;
    *sb += '[';
    for (size_t i = 0; i < n; ++i) md_str_writer_fmt(sb, i ? ", %d" : "%d", v[i]);
    *sb += ']';
    return true;
}

// Terms joined by `or`. `and` binds tighter, so a term with one needs nothing around it.
void write_or(md_str_writer_t* sb, const str_t* terms, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (i) *sb += " or ";
        *sb += terms[i];
    }
}

void push_unique_str(md_array(str_t)* arr, str_t s, md_allocator_i* alloc) {
    for (size_t i = 0; i < md_array_size(*arr); ++i) if (str_eq((*arr)[i], s)) return;
    md_array_push(*arr, s, alloc);
}

// ### THE SYSTEM ###

struct Ctx {
    const md_system_t* sys = nullptr;
    md_allocator_i* temp = nullptr;
    size_t limit = MAX_EXPR_LEN;
    size_t max_terms = MAX_TERMS;
};

int32_t atom_res(const Ctx& c, size_t atom) { return md_system_residue_find_by_atom_idx(c.sys, atom); }
int32_t res_chain(const Ctx& c, int32_t res) { return res >= 0 ? md_system_chain_find_by_residue_idx(c.sys, (size_t)res) : -1; }

// The value of a keyword for a particle, as a string or an int; false where it has none (a particle of no residue has
// no resid, one of no chain no chain)
bool key_str(const Ctx& c, Key k, size_t atom, str_t* out) {
    const md_system_t& sys = *c.sys;
    switch (k) {
    case Key_Name: *out = md_atom_name(&sys.atom, atom); return quotable(*out);
    case Key_Element: {
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom);
        if (z <= 0) return false;
        *out = md_util_element_symbol((md_element_t)z);
        return quotable(*out);
    }
    case Key_Resname: {
        const int32_t r = atom_res(c, atom);
        if (r < 0) return false;
        *out = md_system_residue_name(&sys, (size_t)r);
        return quotable(*out);
    }
    case Key_Chain: {
        const int32_t ch = res_chain(c, atom_res(c, atom));
        if (ch < 0) return false;
        *out = md_system_chain_id(&sys, (size_t)ch);
        return quotable(*out);
    }
    default: return false;
    }
}

bool key_int(const Ctx& c, size_t atom, int32_t* out) {
    const int32_t r = atom_res(c, atom);
    if (r < 0 || !c.sys->residue.seq_id) return false;
    *out = (int32_t)md_system_residue_seq_id(c.sys, (size_t)r);
    return true;
}

const char* key_word(Key k) {
    static const char* words[Key_Count] = {"chain", "resid", "resname", "name", "element"};
    return words[k];
}

// The values the selected particles have for a keyword, in the form a filter takes them. Not usable where one of them
// has none, or there are too many.
struct KeyValues {
    bool usable = false;
    md_array(str_t) strs = 0;       // the strings, in the order first met
    md_array(int32_t) ids = 0;      // resid: sorted, unique
    str_t text = {};                // as written: 'CA', ['N', 'CA'], 10..20
};

void key_values(KeyValues* kv, const Ctx& c, Key k, md_bitset_t sel) {
    md_bitset_iter_t it = md_bitset_iter(sel);
    while (md_bitset_iter_next(&it)) {
        if (k == Key_Resid) {
            int32_t id;
            if (!key_int(c, it.idx, &id)) return;
            md_array_push(kv->ids, id, c.temp);
        } else {
            str_t s;
            if (!key_str(c, k, it.idx, &s)) return;
            push_unique_str(&kv->strs, s, c.temp);
            if (md_array_size(kv->strs) > key_max[k]) return;
        }
    }
    md_str_writer_t sb = md_str_writer_create(c.temp);
    if (k == Key_Resid) {
        const size_t n = sort_unique(kv->ids, md_array_size(kv->ids));
        md_array_shrink(kv->ids, n);
        if (!write_ids(&sb, kv->ids, n)) return;
    } else {
        if (md_array_size(kv->strs) == 0) return;
        write_strings(&sb, kv->strs, md_array_size(kv->strs));
    }
    kv->text = take(&sb, c.temp);
    kv->usable = true;
}

// Whether a particle has one of the values of a keyword: the cheap test that keeps the evaluation to the candidates
// that may pick the selection. Names and residue names exactly, elements ignoring case, as the filter does.
bool key_holds(const Ctx& c, Key k, const KeyValues& kv, size_t atom) {
    if (k == Key_Resid) {
        int32_t id;
        if (!key_int(c, atom, &id)) return false;
        for (size_t i = 0; i < md_array_size(kv.ids); ++i) if (kv.ids[i] == id) return true;
        return false;
    }
    str_t s;
    if (!key_str(c, k, atom, &s)) return false;
    for (size_t i = 0; i < md_array_size(kv.strs); ++i) {
        if (k == Key_Element ? str_eq_ignore_case(kv.strs[i], s) : str_eq(kv.strs[i], s)) return true;
    }
    return false;
}

// ### VALIDATION ###

bool filter_sets(md_array(md_bitset_t)* sets, str_t expr, const md_system_t& sys, md_allocator_i* temp) {
    script::FilterCtx ctx = {};
    ctx.sys = &sys;
    *sets = 0;
    return script::filter_sets(script::Language::V2, sets, expr, ctx, temp).success;
}

// The particles the expression picks, all its sets as one
bool filter_union(md_bitmask_t* out, str_t expr, const md_system_t& sys, md_allocator_i* temp, size_t* num_sets = nullptr) {
    md_array(md_bitset_t) sets = 0;
    if (!filter_sets(&sets, expr, sys, temp)) return false;
    *out = md_bitmask_create(0, (uint32_t)sys.atom.count, temp);
    for (size_t i = 0; i < md_array_size(sets); ++i) md_bitmask_or(out, sets[i]);
    if (num_sets) *num_sets = md_array_size(sets);
    return true;
}

bool picks_exactly(str_t expr, md_bitset_t sel, const md_system_t& sys, md_allocator_i* temp, size_t* num_sets = nullptr) {
    md_bitmask_t m;
    return filter_union(&m, expr, sys, temp, num_sets) && md_bitset_equal(md_bitmask_view(&m), sel);
}

bool compiles(str_t src, const md_system_t& sys, md_allocator_i* temp) {
    script::Program* ir = script::compile(script::Language::V2, src, script::CompileDesc{&sys}, temp);
    const bool ok = script::valid(ir) && script::num_errors(ir) == 0;
    script::program_free(ir);
    return ok;
}

// ### SELECTION: WHAT THE PARTICLES ARE ###

// particles(k=v, ...) of every combination of up to three keywords, optionally within a predefined set that holds the
// selection, which picks it: tested cheaply first, then evaluated. The shortest few.
void describe(md_array(str_t)* out, const Ctx& c, md_bitset_t sel) {
    const md_system_t& sys = *c.sys;
    const uint64_t count = md_bitset_count(sel);

    KeyValues kv[Key_Count];
    for (int k = 0; k < Key_Count; ++k) key_values(&kv[k], c, (Key)k, sel);

    // The predefined sets which hold the selection, each a context the keywords can narrow
    struct Context { str_t name; md_bitmask_t atoms; };
    md_array(Context) contexts = 0;
    md_array_push(contexts, Context{}, c.temp);     // none: the system
    for (size_t i = 0; i < num_predefined; ++i) {
        Context ctx = {str_from_cstr(predefined[i]), {}};
        if (!filter_union(&ctx.atoms, ctx.name, sys, c.temp)) continue;
        const md_bitset_t view = md_bitmask_view(&ctx.atoms);
        if (md_bitset_count(view) < count) continue;
        md_bitmask_t both = md_bitmask_create(0, (uint32_t)sys.atom.count, c.temp);
        md_bitmask_or(&both, view);
        md_bitmask_and(&both, sel);
        if (md_bitset_count(md_bitmask_view(&both)) != count) continue;
        if (md_bitset_count(view) == count) {
            md_array_push(*out, ctx.name, c.temp);  // it is the selection
            continue;
        }
        md_array_push(contexts, ctx, c.temp);
    }

    // The particles of the system that have one of the selection's values, keyword by keyword: what a combination
    // picks is then the intersection of these, without evaluating it
    md_bitmask_t holds[Key_Count] = {};
    for (int k = 0; k < Key_Count; ++k) {
        if (!kv[k].usable) continue;
        holds[k] = md_bitmask_create(0, (uint32_t)sys.atom.count, c.temp);
        for (size_t a = 0; a < sys.atom.count; ++a) if (key_holds(c, (Key)k, kv[k], a)) md_bitmask_set_bit(&holds[k], (uint32_t)a);
    }

    // Every combination which picks as many particles as are selected, and so picks them: it holds them
    struct Found { uint32_t mask; size_t ctx; };
    md_array(Found) fits = 0;
    md_bitmask_t picked = md_bitmask_create(0, (uint32_t)sys.atom.count, c.temp);
    for (uint32_t mask = 1; mask < (1u << Key_Count); ++mask) {
        bool usable = __builtin_popcount(mask) <= 3;
        for (int k = 0; k < Key_Count && usable; ++k) if (mask & (1u << k)) usable = kv[k].usable;
        if (!usable) continue;
        for (size_t ci = 0; ci < md_array_size(contexts); ++ci) {
            md_bitmask_clear(&picked);
            md_bitmask_set_range(&picked, 0, (uint32_t)sys.atom.count);
            if (ci) md_bitmask_and(&picked, md_bitmask_view(&contexts[ci].atoms));
            for (int k = 0; k < Key_Count; ++k) if (mask & (1u << k)) md_bitmask_and(&picked, md_bitmask_view(&holds[k]));
            if (md_bitmask_count(&picked) == count) md_array_push(fits, (Found{mask, ci}), c.temp);
        }
    }
    // Only those that say nothing more than they need: no other picks it with fewer of its keywords, in the system or
    // in its own context. particles(name='CA') is enough; protein and particles(name='CA') adds nothing.
    md_array(str_t) found = 0;
    md_str_writer_t sb = md_str_writer_create(c.temp);
    for (size_t i = 0; i < md_array_size(fits); ++i) {
        const Found& f = fits[i];
        bool minimal = true;
        for (size_t j = 0; j < md_array_size(fits) && minimal; ++j) {
            const Found& g = fits[j];
            if (j == i || (g.mask & ~f.mask)) continue;     // g's keywords are not a subset of f's
            if (g.ctx == 0 && (g.mask != f.mask || f.ctx != 0)) minimal = false;
            else if (g.ctx == f.ctx && g.mask != f.mask) minimal = false;
        }
        if (!minimal) continue;
        if (f.ctx) { sb += contexts[f.ctx].name; sb += " and "; }
        sb += "particles(";
        bool first = true;
        for (int k = 0; k < Key_Count; ++k) {
            if (!(f.mask & (1u << k))) continue;
            if (!first) sb += ", ";
            first = false;
            sb += key_word((Key)k);
            sb += '=';
            sb += kv[k].text;
        }
        sb += ')';
        if (md_str_writer_len(&sb) <= c.limit) md_array_push(found, take(&sb, c.temp), c.temp);
        md_str_writer_reset(&sb);
    }
    // The shortest first, of those which pick it
    for (size_t i = 1; i < md_array_size(found); ++i) {
        for (size_t k = i; k > 0 && found[k].len < found[k - 1].len; --k) { str_t t = found[k]; found[k] = found[k - 1]; found[k - 1] = t; }
    }
    size_t taken = 0;
    for (size_t i = 0; i < md_array_size(found) && taken < MAX_DESCRIPTIONS; ++i) {
        if (!picks_exactly(found[i], sel, sys, c.temp)) continue;
        md_array_push(*out, found[i], c.temp);
        ++taken;
    }
}

// ### SELECTION: WHERE THE PARTICLES ARE ###

// By the residues which hold them, chain by chain: whole residues as runs of ids, particles(chain='A', resid=10..20),
// and part of one by the names of its particles, particles(chain='A', resid=12, name=['N', 'CA']). False where a
// particle is in no residue, a name or chain id cannot be written, or there are more terms than make a line.
bool write_residues(md_str_writer_t* out, const Ctx& c, md_bitset_t sel) {
    const md_system_t& sys = *c.sys;
    if (!sys.residue.seq_id) return false;
    md_array(str_t) terms = 0;
    md_str_writer_t sb = md_str_writer_create(c.temp);
    md_str_writer_t names = md_str_writer_create(c.temp);

    // A run: residues one after another in the system, of one chain, with ids one after another, all whole or all with
    // the same names selected
    int32_t run_chain = -2, run_beg = 0, run_end = 0, run_last_res = -2;
    str_t run_names = {};
    auto flush = [&]() {
        if (run_chain == -2) return;
        sb += "particles(";
        if (run_chain >= 0) {
            sb += "chain=";
            write_quoted(&sb, md_system_chain_id(&sys, (size_t)run_chain));
            sb += ", ";
        }
        if (run_beg == run_end) md_str_writer_fmt(&sb, "resid=%d", run_beg);
        else md_str_writer_fmt(&sb, "resid=%d..%d", run_beg, run_end);
        if (!str_empty(run_names)) { sb += ", name="; sb += run_names; }
        sb += ')';
        md_array_push(terms, take(&sb, c.temp), c.temp);
        run_chain = -2;
    };

    uint64_t covered = 0;
    for (size_t r = 0; r < sys.residue.count; ++r) {
        const md_urange_t ar = md_system_residue_atom_range(&sys, r);
        const uint64_t n = md_bitset_count_range(sel, ar.beg, ar.end);
        if (!n) continue;
        covered += n;
        const int32_t chain = res_chain(c, (int32_t)r);
        if (chain >= 0 && !quotable(md_system_chain_id(&sys, (size_t)chain))) return false;
        str_t nm = {};
        if (n != ar.end - ar.beg) {
            md_array(str_t) list = 0;
            for (uint32_t a = ar.beg; a < ar.end; ++a) {
                if (!md_bitset_test(sel, a)) continue;
                const str_t s = md_atom_name(&sys.atom, a);
                if (!quotable(s)) return false;
                md_array_push(list, s, c.temp);
            }
            write_strings(&names, list, md_array_size(list));
            nm = take(&names, c.temp);
        }
        const int32_t id = (int32_t)md_system_residue_seq_id(&sys, r);
        const bool extends = run_chain == chain && run_last_res + 1 == (int32_t)r && run_end + 1 == id && str_eq(run_names, nm);
        if (!extends) {
            flush();
            if (md_array_size(terms) >= c.max_terms) return false;
            run_chain = chain;
            run_beg = id;
            run_names = nm;
        }
        run_end = id;
        run_last_res = (int32_t)r;
    }
    flush();
    if (covered != md_bitset_count(sel) || md_array_size(terms) == 0 || md_array_size(terms) > c.max_terms) return false;
    write_or(out, terms, md_array_size(terms));
    return md_str_writer_len(out) <= c.limit;
}

// By their positions in the system: particles[4], particles[10:20], particles[[3, 7, 9]], joined by `or`
bool write_positions(md_str_writer_t* out, const Ctx& c, md_bitset_t sel) {
    md_array(str_t) terms = 0;
    md_array(uint32_t) singles = 0;
    md_str_writer_t sb = md_str_writer_create(c.temp);
    md_bitset_iter_t it = md_bitset_iter(sel);
    int64_t beg = -1, end = -1;
    auto flush = [&]() {
        if (beg < 0) return;
        if (end - beg >= 3) {
            md_str_writer_fmt(&sb, "particles[%lld:%lld]", (long long)beg, (long long)end);
            md_array_push(terms, take(&sb, c.temp), c.temp);
        } else {
            for (int64_t i = beg; i < end; ++i) md_array_push(singles, (uint32_t)i, c.temp);
        }
    };
    while (md_bitset_iter_next(&it)) {
        const int64_t i = (int64_t)it.idx;
        if (i == end) { ++end; continue; }
        flush();
        beg = i;
        end = i + 1;
    }
    flush();
    const size_t ns = md_array_size(singles);
    if (ns == 1) md_str_writer_fmt(&sb, "particles[%u]", singles[0]);
    else if (ns > 1) {
        sb += "particles[[";
        for (size_t i = 0; i < ns && md_str_writer_len(&sb) <= c.limit; ++i) md_str_writer_fmt(&sb, i ? ", %u" : "%u", singles[i]);
        sb += "]]";
    }
    if (ns) md_array_push(terms, take(&sb, c.temp), c.temp);
    if (md_array_size(terms) == 0 || md_array_size(terms) > c.max_terms) return false;
    write_or(out, terms, md_array_size(terms));
    return md_str_writer_len(out) <= c.limit;
}

void push_suggestion(md_array(ScriptSuggestion)* out, str_t text, str_t hint, md_allocator_i* alloc) {
    for (size_t i = 0; i < md_array_size(*out); ++i) if (str_eq((*out)[i].text, text)) return;
    ScriptSuggestion s = {str_copy(text, alloc), str_copy(hint, alloc)};
    md_array_push(*out, s, alloc);
}

}  // namespace

md_array(ScriptSuggestion) script_suggest_selection_vial(str_t ident, md_bitset_t sel, const md_system_t& sys, md_allocator_i* alloc) {
    md_array(ScriptSuggestion) out = 0;
    if (md_bitset_empty(sel) || sys.atom.count == 0 || str_empty(ident)) return out;

    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    Ctx c = {&sys, temp};

    md_array(str_t) cand = 0;
    const uint64_t count = md_bitset_count(sel);
    if (count == sys.atom.count) md_array_push(cand, STR_LIT("particles"), temp);
    else describe(&cand, c, sel);

    md_str_writer_t sb = md_str_writer_create(temp);
    if (write_residues(&sb, c, sel)) md_array_push(cand, take(&sb, temp), temp);
    md_str_writer_reset(&sb);
    if (write_positions(&sb, c, sel)) md_array_push(cand, take(&sb, temp), temp);
    md_str_writer_reset(&sb);

    md_str_writer_t line = md_str_writer_create(temp);
    auto offer = [&](str_t expr, bool listed) {
        if (!picks_exactly(expr, sel, sys, temp)) return false;
        md_str_writer_fmt(&sb, "%llu particle%s", (unsigned long long)count, count == 1 ? "" : "s");
        if (listed) sb += ", listed: nothing shorter describes them";
        md_str_writer_fmt(&line, STR_FMT " = " STR_FMT, STR_ARG(ident), STR_ARG(expr));
        push_suggestion(&out, md_str_writer_to_str(&line), md_str_writer_to_str(&sb), alloc);
        md_str_writer_reset(&line);
        md_str_writer_reset(&sb);
        return true;
    };
    for (size_t i = 0; i < md_array_size(cand) && md_array_size(out) < MAX_SUGGESTIONS; ++i) offer(cand[i], false);

    // Nothing describes it: the shortest listing which picks it
    if (md_array_size(out) == 0) {
        c.limit = MAX_LISTED_EXPR_LEN;
        c.max_terms = MAX_LISTED_EXPR_LEN / 16;
        md_array(str_t) listed = 0;
        if (write_residues(&sb, c, sel)) md_array_push(listed, take(&sb, temp), temp);
        md_str_writer_reset(&sb);
        if (write_positions(&sb, c, sel)) md_array_push(listed, take(&sb, temp), temp);
        md_str_writer_reset(&sb);
        if (md_array_size(listed) == 2 && listed[1].len < listed[0].len) { str_t t = listed[0]; listed[0] = listed[1]; listed[1] = t; }
        for (size_t i = 0; i < md_array_size(listed) && !offer(listed[i], true); ++i) {}
    }

    md_arena_allocator_destroy(temp);
    return out;
}

md_array(ScriptSuggestion) script_suggest_measurement_vial(str_t ident, const int32_t* atoms, size_t count, const md_system_t& sys, md_allocator_i* alloc) {
    md_array(ScriptSuggestion) out = 0;
    if (count < 2 || count > 4 || str_empty(ident)) return out;
    for (size_t i = 0; i < count; ++i) {
        if (atoms[i] < 0 || (size_t)atoms[i] >= sys.atom.count) return out;
        for (size_t k = 0; k < i; ++k) if (atoms[k] == atoms[i]) return out;
    }
    static const char* procs[] = {nullptr, nullptr, "distance", "angle", "dihedral"};
    const char* proc = procs[count];

    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    Ctx c = {&sys, temp};
    struct Line { str_t text; str_t hint; };
    md_array(Line) lines = 0;
    md_str_writer_t sb = md_str_writer_create(temp);

    auto push_call = [&](const str_t* args, str_t hint) {
        md_str_writer_fmt(&sb, STR_FMT " = %s(", STR_ARG(ident), proc);
        for (size_t i = 0; i < count; ++i) {
            if (i) sb += ", ";
            sb += args[i];
        }
        sb += ')';
        Line l = {take(&sb, temp), hint};
        md_array_push(lines, l, temp);
    };
    auto picks_atom = [&](str_t expr, int32_t atom) {
        return picks_exactly(expr, md_bitset_range((uint32_t)atom, (uint32_t)atom + 1), sys, temp);
    };

    // By position: always right, whatever the particles are called
    str_t args[4] = {};
    for (size_t i = 0; i < count; ++i) {
        md_str_writer_fmt(&sb, "particles[%d]", atoms[i]);
        args[i] = take(&sb, temp);
    }
    push_call(args, {});

    // By the residue and name, as the file numbers them: particles(chain='A', resid=12, name='CA')
    int32_t res[4] = {-1, -1, -1, -1};
    str_t names[4] = {};
    bool named = sys.residue.seq_id != nullptr;
    bool same_res = true;
    for (size_t i = 0; i < count && named; ++i) {
        res[i] = atom_res(c, (size_t)atoms[i]);
        names[i] = md_atom_name(&sys.atom, (size_t)atoms[i]);
        named = res[i] >= 0 && quotable(names[i]);
        same_res &= res[i] == res[0];
    }
    bool by_res = named;
    for (size_t i = 0; i < count && by_res; ++i) {
        const int32_t chain = res_chain(c, res[i]);
        sb += "particles(";
        if (chain >= 0) {
            const str_t id = md_system_chain_id(&sys, (size_t)chain);
            if (!quotable(id)) { by_res = false; break; }
            sb += "chain=";
            write_quoted(&sb, id);
            sb += ", ";
        }
        md_str_writer_fmt(&sb, "resid=%d, name=", (int)md_system_residue_seq_id(&sys, (size_t)res[i]));
        write_quoted(&sb, names[i]);
        sb += ')';
        args[i] = take(&sb, temp);
        by_res = picks_atom(args[i], atoms[i]);
    }
    md_str_writer_reset(&sb);
    if (by_res) push_call(args, {});

    // The same measurement in every residue of the name, where each of them has each of the particles once
    const str_t resname = named && same_res ? md_system_residue_name(&sys, (size_t)res[0]) : str_t{};
    if (quotable(resname)) {
        size_t num_res = 0;
        bool ok = true;
        for (size_t i = 0; i < count && ok; ++i) {
            sb += "residues(";
            write_quoted(&sb, resname);
            sb += ") and particles(name=";
            write_quoted(&sb, names[i]);
            sb += ')';
            args[i] = take(&sb, temp);
            md_array(md_bitset_t) sets = 0;
            ok = filter_sets(&sets, args[i], sys, temp) && md_array_size(sets) > 0 && (i == 0 || md_array_size(sets) == num_res);
            num_res = md_array_size(sets);
            for (size_t k = 0; k < md_array_size(sets) && ok; ++k) ok = md_bitset_count(sets[k]) == 1;
        }
        if (ok && num_res > 1) {
            md_str_writer_fmt(&sb, "One value for each of the %zu " STR_FMT, num_res, STR_ARG(resname));
            const str_t hint = take(&sb, temp);
            push_call(args, hint);
        }
    }

    for (size_t i = 0; i < md_array_size(lines); ++i) {
        if (compiles(lines[i].text, sys, temp)) push_suggestion(&out, lines[i].text, lines[i].hint, alloc);
    }
    md_arena_allocator_destroy(temp);
    return out;
}
