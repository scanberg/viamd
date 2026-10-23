#include "script_suggestions.h"

#include <md_system.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>

#include <stdlib.h>
#include <string.h>

// Longer than this a line no longer reads at a glance in a menu: it has stopped describing the selection and lists it.
// The indices of a scattered selection are the case in point. Whatever grows past it is dropped while it is written.
static constexpr size_t MAX_EXPR_LEN = 240;
// A selection which nothing describes within that (the atoms within some distance of a residue) is still offered
// listed, by the structures which hold it or by its atoms, up to this. Past it, it is better made with a query.
static constexpr size_t MAX_LISTED_EXPR_LEN = 4096;
static constexpr size_t MAX_SELECTION_SUGGESTIONS = 6;

// A residue holding the selected atoms of other residues in many different ways is not described, it is listed:
// past this many ways (each one a term of its own) the description is given up.
static constexpr size_t MAX_GROUPS = 12;

// At most this many names, elements or residue names make a description of what the atoms are, rather than a list.
static constexpr size_t MAX_NAMES    = 4;
// The atoms of part of a residue are listed by name up to this many, and beyond it by their positions in it
static constexpr size_t MAX_PART_NAMES = 8;
static constexpr size_t MAX_ELEMENTS = 3;
static constexpr size_t MAX_RESNAMES = 3;

// ### WRITING ###

static bool too_long(const md_str_writer_t* sb, size_t limit) {
    return md_str_writer_len(sb) > limit;
}

static bool str_has(str_t s, const char* needle) {
    size_t loc;
    return str_find_str(&loc, s, str_from_cstr(needle));
}

static void write_int(md_str_writer_t* sb, int64_t v) {
    md_str_writer_fmt(sb, "%lld", (long long)v);
}

// A negative bound of a range has to be parenthesised: '-3:5' reads as a subtraction
static void write_bound(md_str_writer_t* sb, int64_t v) {
    md_str_writer_fmt(sb, v < 0 ? "(%lld)" : "%lld", (long long)v);
}

static int cmp_int32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

// Sorts the values and removes the duplicates, returns how many are left
static size_t sort_unique(int32_t* v, size_t n) {
    if (n < 2) return n;
    qsort(v, n, sizeof(int32_t), cmp_int32);
    size_t k = 1;
    for (size_t i = 1; i < n; ++i) {
        if (v[i] != v[k - 1]) v[k++] = v[i];
    }
    return k;
}

// Sorted and unique values as the script reads a list of them: "5", "3:7", "{1:4, 9, 12:13}"
static void write_ranges(md_str_writer_t* sb, const int32_t* v, size_t n) {
    size_t runs = 0;
    for (size_t i = 0; i < n; ++i) {
        if (i == 0 || (int64_t)v[i] != (int64_t)v[i - 1] + 1) ++runs;
    }
    if (runs > 1) *sb += '{';
    for (size_t i = 0; i < n && !too_long(sb, MAX_LISTED_EXPR_LEN);) {
        size_t j = i;
        while (j + 1 < n && (int64_t)v[j + 1] == (int64_t)v[j] + 1) ++j;
        if (i > 0) *sb += ", ";
        if (j > i) {
            write_bound(sb, v[i]);
            *sb += ':';
            write_bound(sb, v[j]);
        } else {
            write_int(sb, v[i]);
        }
        i = j + 1;
    }
    if (runs > 1) *sb += '}';
}

// A string can be written when one of the quotes is not in it, as the script has no escapes, and it is nothing but
// what can be typed: a name padded with zeros (or ending in a carriage return) would be written as it is matched
static bool quotable(str_t s) {
    if (str_empty(s)) return false;
    bool dq = false, sq = false;
    for (size_t i = 0; i < s.len; ++i) {
        const unsigned char ch = (unsigned char)s.ptr[i];
        if (ch < 0x20 || ch == 0x7F) return false;
        dq |= ch == '"';
        sq |= ch == '\'';
    }
    return !(dq && sq);
}

static void write_quoted(md_str_writer_t* sb, str_t s) {
    size_t loc;
    const char q = str_find_char(&loc, s, '"') ? '\'' : '"';
    *sb += q;
    *sb += s;
    *sb += q;
}

// "A" or {"A", "B"}
static void write_strings(md_str_writer_t* sb, const str_t* s, size_t n) {
    if (n > 1) *sb += '{';
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) *sb += ", ";
        write_quoted(sb, s[i]);
    }
    if (n > 1) *sb += '}';
}

static bool contains_str(const str_t* arr, size_t n, str_t s) {
    for (size_t i = 0; i < n; ++i) {
        if (str_eq(arr[i], s)) return true;
    }
    return false;
}

// Terms joined by 'or'. 'in' binds looser than 'or' - 'a in b or c' is 'a in (b or c)' - so a term with one is put in
// parentheses.
static void write_or(md_str_writer_t* sb, const str_t* terms, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) *sb += " or ";
        const bool paren = n > 1 && str_has(terms[i], " in ");
        if (paren) *sb += '(';
        *sb += terms[i];
        if (paren) *sb += ')';
    }
}

// A term within a context: 'resid(10:20) in chain("A")'. A term which is itself an 'in' nests, which reads from the
// inside out: 'name("CA") in resid(21) in chain("A")'.
static void write_in(md_str_writer_t* sb, str_t term, str_t context) {
    *sb += term;
    *sb += " in ";
    *sb += context;
}

// An operand of 'and': a term with an 'in' or an 'or' in parentheses, as both bind looser
static void write_and_operand(md_str_writer_t* sb, str_t term) {
    const bool paren = str_has(term, " in ") || str_has(term, " or ");
    if (paren) *sb += '(';
    *sb += term;
    if (paren) *sb += ')';
}

static str_t take(md_str_writer_t* sb, md_allocator_i* alloc) {
    str_t s = str_copy(md_str_writer_to_str(sb), alloc);
    md_str_writer_reset(sb);
    return s;
}

// ### THE SYSTEM ###

struct Ctx {
    const md_system_t* sys = nullptr;
    md_allocator_i* temp = nullptr;
    md_bitset_t sel = {};
    size_t limit = MAX_EXPR_LEN;        // the longest expression written
    size_t max_groups = MAX_GROUPS;     // the most ways of holding part of a residue written (Groups)
    // [chain] Whether the residues of the chain have sequence ids no other residue of it has, so that
    // resid(i) in chain("A") picks one: 0 not yet known, 1 they have, 2 they have not
    uint8_t* chain_seq = nullptr;
};

static void ctx_init(Ctx* c, const md_system_t& sys, md_bitset_t sel, md_allocator_i* temp) {
    c->sys  = &sys;
    c->temp = temp;
    c->sel  = sel;
    const size_t num_chains = md_system_chain_count(&sys);
    c->chain_seq = (uint8_t*)md_alloc(temp, num_chains + 1);
    MEMSET(c->chain_seq, 0, num_chains + 1);
}

// The chain of a residue (md_chain_data_t), -1 for none
static int32_t res_chain(const Ctx& c, int32_t res) {
    return res >= 0 ? md_system_chain_find_by_residue_idx(c.sys, (size_t)res) : -1;
}

// The chain of an instance: that of its first residue. The residues of an instance are of one chain, or of none,
// except in a molecule of a topology which holds several strands (a DNA duplex), whose chains are all inferred.
static int32_t inst_first_chain(const Ctx& c, int32_t inst) {
    if (inst < 0) return -1;
    const md_urange_t residues = md_system_instance_residue_range(c.sys, (size_t)inst);
    return residues.beg < residues.end ? res_chain(c, (int32_t)residues.beg) : -1;
}

// Every atom of the chain is selected
static bool chain_whole(const Ctx& c, int32_t chain) {
    const md_urange_t* end = md_system_chain_run_end(c.sys, (size_t)chain);
    for (const md_urange_t* it = md_system_chain_run_beg(c.sys, (size_t)chain); it && it < end; ++it) {
        const md_urange_t r = md_system_chain_run_atom_range(c.sys, *it);
        if (r.beg < r.end && !md_bitset_test_all_range(c.sel, r.beg, r.end)) return false;
    }
    return true;
}

// A chain is written by the id the data gives it: chain("A"), all of the atoms of the id. A chain inferred where the
// data gives none is written by its position: chain(2).
static bool chain_by_id(const Ctx& c, int32_t chain) {
    return chain >= 0 && quotable(md_system_chain_id(c.sys, (size_t)chain));
}

static void write_chain(md_str_writer_t* sb, const Ctx& c, int32_t chain) {
    if (chain_by_id(c, chain)) {
        *sb += "chain(";
        write_quoted(sb, md_system_chain_id(c.sys, (size_t)chain));
        *sb += ')';
    } else {
        md_str_writer_fmt(sb, "chain(%d)", chain + 1);
    }
}

// The sequence ids of a chain are only looked at when it holds part of what is described. Those of all of its
// residues: the waters the data files under a chain are numbered with it.
static bool chain_by_seq(Ctx* c, int32_t chain) {
    if (chain < 0 || !c->sys->residue.seq_id) return false;
    if (c->chain_seq[chain] == 0) {
        md_array(int32_t) ids = 0;
        const md_urange_t* end = md_system_chain_run_end(c->sys, (size_t)chain);
        for (const md_urange_t* it = md_system_chain_run_beg(c->sys, (size_t)chain); it && it < end; ++it) {
            for (uint32_t k = it->beg; k < it->end; ++k) md_array_push(ids, (int32_t)c->sys->residue.seq_id[k], c->temp);
        }
        const size_t n = md_array_size(ids);
        c->chain_seq[chain] = n > 0 && sort_unique(ids, n) == n ? 1 : 2;
    }
    return c->chain_seq[chain] == 1;
}

// What the selection holds: the residues with a selected atom, whole or in part, the instances which hold them,
// and the selected atoms which are in no residue.
struct Breakdown {
    md_array(int32_t) res       = 0;
    md_array(uint8_t) res_full = 0;     // [res] every atom of the residue is selected
    md_array(int32_t) res_inst = 0;     // [res] its instance, -1 for none
    md_array(int32_t) res_slot = 0;     // [res] where its instance is in inst, -1 for none
    md_array(int32_t) inst      = 0;    // the instances of the residues, in order
    md_array(uint8_t) inst_full = 0;    // [inst] every atom of the instance is selected
    md_array(int32_t) loose     = 0;    // selected atoms of no residue
    bool residues_full = true;          // the selection is whole residues (and nothing else)
    bool insts_full = true;             // the selection is whole instances (and nothing else)
};

static void breakdown_init(Breakdown* b, const Ctx& c) {
    const md_system_t& sys = *c.sys;
    const md_residue_data_t* res = &sys.residue;
    const md_bitset_t sel = c.sel;
    const uint64_t popcount = md_bitset_count(sel);

    uint64_t in_residues = 0;
    if (res->count && res->atom_offset && !md_bitset_empty(sel)) {
        // The first residue which ends past the first selected atom
        size_t lo = 0, hi = res->count;
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (res->atom_offset[mid + 1] <= sel.beg) lo = mid + 1;
            else hi = mid;
        }
        int32_t inst = -1;
        md_urange_t inst_residues = {0, 0};
        for (size_t i = lo; i < res->count && res->atom_offset[i] < sel.end; ++i) {
            const md_urange_t r = md_residue_atom_range(res, i);
            const uint64_t n = md_bitset_count_range(sel, r.beg, r.end);
            if (!n) continue;
            in_residues += n;

            // The residues of an instance are consecutive, so the instance of the previous one is the likely one
            if (!(inst_residues.beg <= i && i < inst_residues.end)) {
                inst = md_instance_find_by_residue_idx(&sys.instance, i);
                inst_residues = inst >= 0 ? md_system_instance_residue_range(&sys, inst) : md_urange_t{0, 0};
                if (inst >= 0) {
                    const md_urange_t ar = md_system_instance_atom_range(&sys, inst);
                    const bool full = md_bitset_test_all_range(sel, ar.beg, ar.end);
                    md_array_push(b->inst, inst, c.temp);
                    md_array_push(b->inst_full, (uint8_t)full, c.temp);
                    b->insts_full &= full;
                }
            }
            const bool full = n == (uint64_t)(r.end - r.beg);
            md_array_push(b->res, (int32_t)i, c.temp);
            md_array_push(b->res_full, (uint8_t)full, c.temp);
            md_array_push(b->res_inst, inst, c.temp);
            md_array_push(b->res_slot, inst >= 0 ? (int32_t)md_array_size(b->inst) - 1 : -1, c.temp);
            b->residues_full &= full;
            b->insts_full &= inst >= 0;
        }
    }

    if (in_residues < popcount) {
        md_bitset_iter_t it = md_bitset_iter(sel);
        while (md_bitset_iter_next(&it)) {
            if (md_residue_find_by_atom_idx(res, it.idx) == -1) {
                md_array_push(b->loose, (int32_t)it.idx, c.temp);
            }
        }
        b->residues_full = false;
        b->insts_full = false;
    }
}

// ### THE STRUCTURES WHICH HOLD THE SELECTION ###

// The molecule which an instance is all of, -1 when it is part of a larger one (a ligand bonded to a protein, a metal
// to the residues around it), of several (the waters of an mmCIF file), or there are none
static int32_t molecule_of_instance(const md_system_t& sys, int32_t inst) {
    const md_urange_t r = md_system_instance_atom_range(&sys, (size_t)inst);
    if (r.beg >= r.end) return -1;
    const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, r.beg);
    return mol >= 0 && md_system_molecule_atom_count(&sys, mol) == r.end - r.beg ? mol : -1;
}

// The molecules of an instance which is not one molecule: those inside it (the waters of a chain, which an mmCIF file
// gives one instance), or those it is part of (a protein and the metal coordinated to it are two instances of one
// molecule) where they are wholly among the instances (sorted). None where they are not: its residues are written then.
static void molecules_of_instance(md_array(int32_t)* mols, int32_t inst, const int32_t* insts, size_t num_insts, const Ctx& c) {
    const md_system_t& sys = *c.sys;
    const md_molecule_data_t& st = sys.molecule;
    const md_urange_t r = md_system_instance_atom_range(&sys, (size_t)inst);
    md_array(int32_t) touched = 0;
    int32_t prev = -1;
    for (uint32_t i = r.beg; i < r.end; ++i) {
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, i);
        if (mol < 0) return;
        if (mol != prev) md_array_push(touched, mol, c.temp);
        prev = mol;
    }
    const size_t n = sort_unique(touched, md_array_size(touched));
    uint64_t size = 0;
    for (size_t k = 0; k < n; ++k) size += md_system_molecule_atom_count(&sys, touched[k]);
    // Molecules inside the instance, which is whole
    if (size == r.end - r.beg) {
        for (size_t k = 0; k < n; ++k) md_array_push(*mols, touched[k], c.temp);
        return;
    }
    for (size_t k = 0; k < n; ++k) {
        for (uint32_t s = st.offset[touched[k]]; s < st.offset[touched[k] + 1]; ++s) {
            const md_instance_idx_t other = md_system_instance_find_by_atom_idx(&sys, (size_t)st.atom_idx[s]);
            if (other < 0 || !bsearch(&other, insts, num_insts, sizeof(int32_t), cmp_int32)) return;
        }
    }
    for (size_t k = 0; k < n; ++k) md_array_push(*mols, touched[k], c.temp);
}

// Whole chains and molecules as terms of an 'or'. Chains by their ids: chain("A"), chain({"A", "B"}); a chain inferred
// where the data gives none (or every chain, by_id false) by its position among the chains: chain(2). A molecule
// which is in no chain, or whose chain is not whole, by its position among the molecules: molecule(5:9). An instance
// which is part of a molecule not selected whole (a ligand bonded to a protein) is none of these: it is pushed to
// unwritten, and its residues are written instead.
static bool whole_terms(md_array(str_t)* terms, md_array(int32_t)* unwritten, const int32_t* chains, size_t num_chains, const int32_t* insts, size_t num_insts, const Ctx& c, bool by_id) {
    md_str_writer_t sb = md_str_writer_create(c.temp);
    md_array(str_t)   ids = 0;
    md_array(int32_t) pos = 0;
    for (size_t i = 0; i < num_chains; ++i) {
        if (by_id && chain_by_id(c, chains[i])) md_array_push(ids, md_system_chain_id(c.sys, (size_t)chains[i]), c.temp);
        else md_array_push(pos, chains[i] + 1, c.temp);
    }
    md_array(int32_t) mols = 0;
    md_array(int32_t) sorted = 0;
    size_t num_sorted = 0;
    for (size_t i = 0; i < num_insts; ++i) {
        const int32_t mol = molecule_of_instance(*c.sys, insts[i]);
        if (mol >= 0) {
            md_array_push(mols, mol + 1, c.temp);
            continue;
        }
        if (!sorted) {
            md_array_push_array(sorted, insts, num_insts, c.temp);
            num_sorted = sort_unique(sorted, num_insts);
        }
        const size_t before = md_array_size(mols);
        molecules_of_instance(&mols, insts[i], sorted, num_sorted, c);
        if (md_array_size(mols) == before) md_array_push(*unwritten, insts[i], c.temp);
        else for (size_t k = before; k < md_array_size(mols); ++k) mols[k] += 1;
    }

    if (md_array_size(ids)) {
        // Thousands of chains are not written by their ids
        if (md_array_size(ids) * 5 > c.limit) return false;
        sb += "chain(";
        write_strings(&sb, ids, md_array_size(ids));
        sb += ')';
        md_array_push(*terms, take(&sb, c.temp), c.temp);
    }
    if (md_array_size(pos)) {
        sb += "chain(";
        write_ranges(&sb, pos, sort_unique(pos, md_array_size(pos)));
        sb += ')';
        md_array_push(*terms, take(&sb, c.temp), c.temp);
    }
    if (md_array_size(mols)) {
        sb += "molecule(";
        write_ranges(&sb, mols, sort_unique(mols, md_array_size(mols)));
        sb += ')';
        md_array_push(*terms, take(&sb, c.temp), c.temp);
    }
    size_t len = 0;
    for (size_t i = 0; i < md_array_size(*terms); ++i) len += (*terms)[i].len;
    return len <= c.limit;
}

// The chains the selection holds whole, and the instances it holds whole which are not all of whole chains: in no
// chain, or in one which is not whole (a ligand of chain A, without the rest of A)
static void whole_chains_and_insts(md_array(int32_t)* chains, md_array(int32_t)* insts, const Breakdown& b, const Ctx& c) {
    const size_t num_chains = md_system_chain_count(c.sys);
    // [chain] 0 not seen, 1 whole, 2 not whole
    uint8_t* state = (uint8_t*)md_alloc(c.temp, num_chains + 1);
    MEMSET(state, 0, num_chains + 1);
    for (size_t k = 0; k < md_array_size(b.res); ++k) {
        const int32_t chain = res_chain(c, b.res[k]);
        if (chain < 0 || state[chain]) continue;
        state[chain] = chain_whole(c, chain) ? 1 : 2;
        if (state[chain] == 1) md_array_push(*chains, chain, c.temp);
    }
    for (size_t i = 0; i < md_array_size(b.inst); ++i) {
        if (!b.inst_full[i]) continue;
        const md_urange_t residues = md_system_instance_residue_range(c.sys, (size_t)b.inst[i]);
        bool in_whole = true;
        for (uint32_t k = residues.beg; k < residues.end && in_whole; ++k) {
            const int32_t chain = res_chain(c, (int32_t)k);
            in_whole = chain >= 0 && state[chain] == 1;
        }
        if (!in_whole) md_array_push(*insts, b.inst[i], c.temp);
    }
}

// The selected atoms of a residue which holds others as well, as the script picks them within it: by their names
// where those pick them and nothing else in it, name({"N", "CA", "C"}); else by their positions in it, atom(2:4).
static void part_term(md_str_writer_t* sb, int32_t res_idx, const Ctx& c, bool by_name) {
    const md_system_t& sys = *c.sys;
    const md_urange_t r = md_residue_atom_range(&sys.residue, res_idx);

    if (by_name) {
        md_array(str_t) names = 0;
        bool ok = true;
        for (uint32_t i = r.beg; i < r.end && ok; ++i) {
            if (!md_bitset_test(c.sel, i)) continue;
            const str_t name = md_atom_name(&sys.atom, i);
            ok = quotable(name);
            if (ok && !contains_str(names, md_array_size(names), name)) {
                md_array_push(names, name, c.temp);
                ok = md_array_size(names) <= MAX_PART_NAMES;
            }
        }
        // An atom left out which goes by one of the names would be picked as well
        for (uint32_t i = r.beg; i < r.end && ok; ++i) {
            if (!md_bitset_test(c.sel, i) && contains_str(names, md_array_size(names), md_atom_name(&sys.atom, i))) ok = false;
        }
        if (ok) {
            *sb += "name(";
            write_strings(sb, names, md_array_size(names));
            *sb += ')';
            return;
        }
    }

    md_array(int32_t) idx = 0;
    for (uint32_t i = r.beg; i < r.end; ++i) {
        if (md_bitset_test(c.sel, i)) md_array_push(idx, (int32_t)(i - r.beg + 1), c.temp);
    }
    *sb += "atom(";
    write_ranges(sb, idx, md_array_size(idx));
    *sb += ')';
}

// Residues which hold the selected atoms the same way, by what picks those atoms within them (part_term)
struct Groups {
    md_array(str_t) term = 0;
    md_array(md_array(int32_t)) vals = 0;   // what the residues are written by: their positions or sequence ids
};

static bool groups_add(Groups* g, str_t term, int32_t val, size_t max_groups, md_allocator_i* temp) {
    for (size_t i = 0; i < md_array_size(g->term); ++i) {
        if (str_eq(g->term[i], term)) {
            md_array_push(g->vals[i], val, temp);
            return true;
        }
    }
    if (md_array_size(g->term) == max_groups) return false;
    md_array(int32_t) vals = 0;
    md_array_push(vals, val, temp);
    md_array_push(g->term, term, temp);
    md_array_push(g->vals, vals, temp);
    return true;
}

// "term in selector(vals)" for each group, and "selector(vals)" for the whole residues, pushed to terms
static void group_terms(md_array(str_t)* terms, md_array(int32_t) whole, const Groups& g, const char* selector, md_allocator_i* temp) {
    md_str_writer_t sb = md_str_writer_create(temp);
    if (md_array_size(whole)) {
        md_str_writer_fmt(&sb, "%s(", selector);
        write_ranges(&sb, whole, sort_unique(whole, md_array_size(whole)));
        sb += ')';
        md_array_push(*terms, take(&sb, temp), temp);
    }
    for (size_t i = 0; i < md_array_size(g.term); ++i) {
        md_str_writer_fmt(&sb, STR_FMT " in %s(", STR_ARG(g.term[i]), selector);
        write_ranges(&sb, g.vals[i], sort_unique(g.vals[i], md_array_size(g.vals[i])));
        sb += ')';
        md_array_push(*terms, take(&sb, temp), temp);
    }
}

// The selection by the structures which hold it: the chains and molecules it holds whole, the parts of several residues
// it holds whole of the others (a glycan, a strand bonded to another), the residues it holds whole of the rest, and what
// it holds of the residues it holds a part of (part_term). By id: the chains by their ids, the
// residues of a chain by their sequence ids within it and the atoms by their names, where those pick them alone. Else
// everything by position. False when it does not make a line short enough to read.
static bool write_structures(md_str_writer_t* out, const Breakdown& b, Ctx* c, bool by_id) {
    md_allocator_i* temp = c->temp;
    md_array(str_t) terms = 0;

    md_array(int32_t) whole_chains = 0;
    md_array(int32_t) whole_insts = 0;
    md_array(int32_t) unwritten = 0;
    whole_chains_and_insts(&whole_chains, &whole_insts, b, *c);
    if (by_id) {
        // By the ids of the data, a molecule of a chain which has them is its residues within the chain, resid(802) in
        // chain("B"), rather than its position among the molecules
        size_t n = 0;
        for (size_t i = 0; i < md_array_size(whole_insts); ++i) {
            const int32_t chain = inst_first_chain(*c, whole_insts[i]);
            if (chain_by_id(*c, chain) && chain_by_seq(c, chain)) md_array_push(unwritten, whole_insts[i], temp);
            else whole_insts[n++] = whole_insts[i];
        }
        md_array_shrink(whole_insts, n);
    }
    if (!whole_terms(&terms, &unwritten, whole_chains, md_array_size(whole_chains), whole_insts, md_array_size(whole_insts), *c, by_id)) return false;
    const size_t num_unwritten = sort_unique(unwritten, md_array_size(unwritten));
    // What the terms add up to: a selection held in thousands of chains is given up on before all of them are written
    size_t len = 0;
    for (size_t i = 0; i < md_array_size(terms); ++i) len += terms[i].len;

    // Written whole already
    const size_t num_chains = md_system_chain_count(c->sys);
    uint8_t* chain_written = (uint8_t*)md_alloc(temp, num_chains + 1);
    MEMSET(chain_written, 0, num_chains + 1);
    for (size_t i = 0; i < md_array_size(whole_chains); ++i) chain_written[whole_chains[i]] = 1;

    // The residues of a chain by their sequence ids within it, gathered over all of its instances (they need not be
    // consecutive): a slot per chain
    int32_t* chain_slot = (int32_t*)md_alloc(temp, (num_chains + 1) * sizeof(int32_t));
    for (size_t i = 0; i <= num_chains; ++i) chain_slot[i] = -1;
    md_array(int32_t) slot_chain = 0;
    md_array(md_array(int32_t)) slot_whole = 0;
    md_array(Groups) slot_parts = 0;

    // The residues which are written by their position in the system, whichever instance they are of, and the parts of
    // several residues it holds whole of them
    md_array(int32_t) whole = 0;
    md_array(int32_t) whole_parts = 0;
    Groups parts = {};
    int32_t cur_part = -1;
    bool cur_part_whole = false;

    md_str_writer_t sb = md_str_writer_create(temp);
    const size_t num_residues = md_array_size(b.res);
    for (size_t k = 0; k < num_residues; ++k) {
        const int32_t inst = b.res_inst[k];
        const int32_t chain = res_chain(*c, b.res[k]);
        // Its instance is written whole, as a molecule or with its chain
        if (inst >= 0 && b.inst_full[b.res_slot[k]] && !bsearch(&inst, unwritten, num_unwritten, sizeof(int32_t), cmp_int32)) continue;
        if (chain >= 0 && chain_written[chain]) continue;

        if (by_id && chain_by_seq(c, chain)) {
            if (chain_slot[chain] < 0) {
                chain_slot[chain] = (int32_t)md_array_size(slot_chain);
                md_array_push(slot_chain, chain, temp);
                md_array_push(slot_whole, (md_array(int32_t))0, temp);
                md_array_push(slot_parts, Groups{}, temp);
            }
            const int32_t slot = chain_slot[chain];
            const int32_t seq = (int32_t)md_residue_seq_id(&c->sys->residue, b.res[k]);
            if (b.res_full[k]) {
                md_array_push(slot_whole[slot], seq, temp);
            } else {
                part_term(&sb, b.res[k], *c, true);
                if (too_long(&sb, c->limit) || !groups_add(&slot_parts[slot], take(&sb, temp), seq, c->max_groups, temp)) return false;
            }
        } else {
            const int32_t pos = b.res[k] + 1;
            if (b.res_full[k]) {
                const int32_t part = md_system_part_find_by_residue_idx(c->sys, (size_t)b.res[k]);
                if (part >= 0 && part != cur_part) {
                    cur_part = part;
                    const md_urange_t pc = md_system_part_residue_range(c->sys, (size_t)part);
                    const md_urange_t pa = md_system_part_atom_range(c->sys, (size_t)part);
                    cur_part_whole = pc.end - pc.beg > 1 && md_bitset_test_all_range(c->sel, pa.beg, pa.end);
                    if (cur_part_whole) md_array_push(whole_parts, part + 1, temp);
                }
                if (part >= 0 && cur_part_whole) continue;
                md_array_push(whole, pos, temp);
            } else {
                part_term(&sb, b.res[k], *c, by_id);
                if (too_long(&sb, c->limit) || !groups_add(&parts, take(&sb, temp), pos, c->max_groups, temp)) return false;
            }
        }
    }

    // Each term within its chain on its own, which reads more easily than one context around several:
    // '(resid(10:20) in chain("A")) or (name("CA") in resid(21) in chain("A"))'
    for (size_t slot = 0; slot < md_array_size(slot_chain); ++slot) {
        md_array(str_t) inner = 0;
        group_terms(&inner, slot_whole[slot], slot_parts[slot], "resid", temp);
        write_chain(&sb, *c, slot_chain[slot]);
        const str_t chain = take(&sb, temp);
        for (size_t k = 0; k < md_array_size(inner); ++k) {
            write_in(&sb, inner[k], chain);
            len += md_str_writer_len(&sb) + 6;
            if (len > c->limit) return false;
            md_array_push(terms, take(&sb, temp), temp);
        }
    }
    if (md_array_size(whole_parts)) {
        sb += "part(";
        write_ranges(&sb, whole_parts, md_array_size(whole_parts));
        sb += ')';
        md_array_push(terms, take(&sb, temp), temp);
    }
    group_terms(&terms, whole, parts, "residue", temp);

    if (md_array_size(b.loose)) {
        sb += "atom(";
        md_array(int32_t) idx = 0;
        for (size_t i = 0; i < md_array_size(b.loose); ++i) md_array_push(idx, b.loose[i] + 1, temp);
        write_ranges(&sb, idx, md_array_size(idx));
        sb += ')';
        md_array_push(terms, take(&sb, temp), temp);
    }

    write_or(out, terms, md_array_size(terms));
    return md_array_size(terms) > 0 && !too_long(out, c->limit);
}

static bool write_atoms(md_str_writer_t* out, const Ctx& c) {
    // A run of consecutive atoms of any length is written in one, so it is the number of runs which is bounded
    md_array(md_urange_t) runs = 0;
    md_bitset_iter_t it = md_bitset_iter(c.sel);
    while (md_bitset_iter_next(&it)) {
        const size_t n = md_array_size(runs);
        if (n && runs[n - 1].end == it.idx) {
            runs[n - 1].end += 1;
        } else {
            if (n * 3 > c.limit) return false;
            md_urange_t r = {it.idx, it.idx + 1};
            md_array_push(runs, r, c.temp);
        }
    }
    const size_t n = md_array_size(runs);
    *out += "atom(";
    if (n > 1) *out += '{';
    for (size_t i = 0; i < n; ++i) {
        if (i) *out += ", ";
        if (runs[i].end - runs[i].beg == 1) md_str_writer_fmt(out, "%u", runs[i].beg + 1);
        else md_str_writer_fmt(out, "%u:%u", runs[i].beg + 1, runs[i].end);
    }
    if (n > 1) *out += '}';
    *out += ')';
    return n > 0 && !too_long(out, c.limit);
}

// ### WHAT THE SELECTED ATOMS ARE ###
//
// A predicate holds for every selected atom: what they are by name or element, the part of a residue they are
// (backbone()), or the kind or name of the residues they make up whole (water(), resname("HEM")). Where it holds for
// nothing else in the system it is the selection by itself. Where it holds for other atoms as well, it is the
// selection within a scope which leaves those out: the residues of its kind, of its names, the instances which hold
// the selection, or the residues which do.

enum Scope : uint32_t {
    Scope_System    = 1,
    Scope_Class     = 2,    // the residues of the one kind the selected ones are of: protein(), water(), ...
    Scope_Resname   = 4,    // the residues of the names of the selected ones
    Scope_Inst      = 8,    // the chains which hold the selection (and the molecules in no chain which do)
    Scope_Res       = 16,   // the residues which hold the selection
    Scope_All       = 31,
};

struct Pred {
    str_t text = {};
    bool res_level = false;                 // of whole residues (true), or of atoms
    // Atoms: by their type (name, element), or by their flags
    uint8_t* type_take = nullptr;           // [particle types]
    uint32_t flags = 0;                     // MD_ATOM_FLAG_...
    // Residues: of a kind, or of a name
    int res_kind = -1;
    md_array(str_t) res_names = 0;
};

struct Scopes {
    uint32_t avail = 0;                     // the scopes there are for the selection (and which say something)
    str_t text[5] = {};                     // indexed by the bit of the scope
    int class_kind = -1;
    md_array(str_t) resnames = 0;
    uint8_t* res_touched = nullptr;         // [residue]
    uint8_t* res_in_scope = nullptr;        // [residue] of the chains which hold the selection, or of a molecule in no chain which does
};

static int scope_slot(uint32_t scope) {
    int i = 0;
    while (scope > 1) { scope >>= 1; ++i; }
    return i;
}

static bool pred_atom(const Pred& p, const md_system_t& sys, size_t i) {
    if (p.type_take) {
        const md_particle_type_idx_t t = sys.atom.type_idx ? sys.atom.type_idx[i] : 0;
        return (size_t)t < sys.atom.type.count && p.type_take[t];
    }
    return ((uint32_t)md_atom_flags(&sys.atom, i) & p.flags) != 0;
}

static bool pred_res(const Pred& p, const md_system_t& sys, size_t res) {
    if (p.res_kind >= 0) return (int)md_residue_kind(&sys.residue, res) == p.res_kind;
    return contains_str(p.res_names, md_array_size(p.res_names), md_residue_name(&sys.residue, res));
}

// The scopes in which something other than the selection fulfils the predicate
static uint32_t pred_violations(const Pred& p, const Ctx& c, const Scopes& s) {
    const md_system_t& sys = *c.sys;
    const md_residue_data_t* res = &sys.residue;
    uint32_t violated = 0;

    auto classify = [&](size_t ci) {
        uint32_t v = Scope_System;
        if (s.class_kind >= 0 && (int)md_residue_kind(res, ci) == s.class_kind) v |= Scope_Class;
        if (md_array_size(s.resnames) && contains_str(s.resnames, md_array_size(s.resnames), md_residue_name(res, ci))) v |= Scope_Resname;
        if (s.res_in_scope[ci]) v |= Scope_Inst;
        if (s.res_touched[ci]) v |= Scope_Res;
        return v;
    };

    // The atoms before, between and after the residues belong to no scope but the system
    auto loose_atoms = [&](uint32_t beg, uint32_t end) {
        if (p.res_level || (violated & Scope_System)) return;
        for (uint32_t i = beg; i < end; ++i) {
            if (!md_bitset_test(c.sel, i) && pred_atom(p, sys, i)) {
                violated |= Scope_System;
                return;
            }
        }
    };

    uint32_t prev_end = 0;
    for (size_t ci = 0; ci < res->count; ++ci) {
        const md_urange_t r = md_residue_atom_range(res, ci);
        loose_atoms(prev_end, r.beg);
        prev_end = r.end;

        bool violation = false;
        if (p.res_level) {
            // The selected residues are whole, so any other of the kind is a violation
            violation = !s.res_touched[ci] && pred_res(p, sys, ci);
        } else {
            for (uint32_t i = r.beg; i < r.end && !violation; ++i) {
                violation = !md_bitset_test(c.sel, i) && pred_atom(p, sys, i);
            }
        }
        if (violation) {
            violated |= classify(ci);
            // The residues which hold the selection are within every other scope: nothing is left
            if (violated & Scope_Res) return Scope_All;
        }
    }
    loose_atoms(prev_end, (uint32_t)sys.atom.count);
    return violated;
}

// The class of a kind of residue, as the script selects it
static const char* class_selector(int kind) {
    switch (kind) {
    case MD_RESIDUE_KIND_AMINO_ACID: return "protein()";
    case MD_RESIDUE_KIND_NUCLEOTIDE: return "nucleic()";
    case MD_RESIDUE_KIND_WATER:        return "water()";
    case MD_RESIDUE_KIND_ION:          return "ion()";
    case MD_RESIDUE_KIND_LIPID:        return "lipid()";
    default: return nullptr;
    }
}

static void scopes_init(Scopes* s, const Breakdown& b, Ctx* c) {
    const md_system_t& sys = *c->sys;
    md_allocator_i* temp = c->temp;
    const size_t num_residues = sys.residue.count;

    s->res_touched   = (uint8_t*)md_alloc(temp, num_residues + 1);
    s->res_in_scope = (uint8_t*)md_alloc(temp, num_residues + 1);
    MEMSET(s->res_touched, 0, num_residues + 1);
    MEMSET(s->res_in_scope, 0, num_residues + 1);
    for (size_t i = 0; i < md_array_size(b.res); ++i) s->res_touched[b.res[i]] = 1;

    // The chains which hold the selection, all of their residues, and the molecules in no chain which hold it
    md_array(int32_t) chains = 0;
    md_array(int32_t) loose_insts = 0;
    bool scope_whole = true;    // the selection is all of them: then the scope says nothing
    {
        const size_t num_chains = md_system_chain_count(&sys);
        uint8_t* chain_seen = (uint8_t*)md_alloc(temp, num_chains + 1);
        MEMSET(chain_seen, 0, num_chains + 1);
        for (size_t k = 0; k < md_array_size(b.res); ++k) {
            const int32_t chain = res_chain(*c, b.res[k]);
            if (chain < 0 || chain_seen[chain]) continue;
            chain_seen[chain] = 1;
            md_array_push(chains, chain, temp);
            scope_whole = scope_whole && chain_whole(*c, chain);
            const md_urange_t* end = md_system_chain_run_end(&sys, (size_t)chain);
            for (const md_urange_t* it = md_system_chain_run_beg(&sys, (size_t)chain); it && it < end; ++it) {
                for (uint32_t ci = it->beg; ci < it->end; ++ci) s->res_in_scope[ci] = 1;
            }
        }
        for (size_t i = 0; i < md_array_size(b.inst); ++i) {
            if (inst_first_chain(*c, b.inst[i]) >= 0) continue;
            md_array_push(loose_insts, b.inst[i], temp);
            scope_whole = scope_whole && b.inst_full[i] != 0;
            const md_urange_t residues = md_system_instance_residue_range(&sys, (size_t)b.inst[i]);
            for (uint32_t ci = residues.beg; ci < residues.end; ++ci) s->res_in_scope[ci] = 1;
        }
    }

    s->avail = Scope_System;
    // A scope is of the residues which hold the selection: selected atoms which are in none fit none of them
    if (md_array_size(b.loose) || md_array_size(b.res) == 0) return;

    // Of one kind
    int kind = -1;
    for (size_t i = 0; i < md_array_size(b.res); ++i) {
        const int k = (int)md_residue_kind(&sys.residue, b.res[i]);
        if (i == 0) kind = k;
        else if (k != kind) kind = -1;
    }
    if (class_selector(kind)) {
        s->class_kind = kind;
        s->avail |= Scope_Class;
        s->text[scope_slot(Scope_Class)] = str_from_cstr(class_selector(kind));
    }

    // Of a few names
    bool names_ok = true;
    for (size_t i = 0; i < md_array_size(b.res) && names_ok; ++i) {
        const str_t name = md_residue_name(&sys.residue, b.res[i]);
        names_ok = quotable(name);
        if (names_ok && !contains_str(s->resnames, md_array_size(s->resnames), name)) {
            md_array_push(s->resnames, name, temp);
            names_ok = md_array_size(s->resnames) <= MAX_RESNAMES;
        }
    }
    if (names_ok) {
        md_str_writer_t sb = md_str_writer_create(temp);
        sb += "resname(";
        write_strings(&sb, s->resnames, md_array_size(s->resnames));
        sb += ')';
        s->text[scope_slot(Scope_Resname)] = take(&sb, temp);
        s->avail |= Scope_Resname;
    } else {
        md_array_shrink(s->resnames, 0);
    }

    // The chains which hold the selection (and the molecules in no chain), when it is not the whole of them: then it is
    // written as them
    bool all_in_inst = true;
    for (size_t i = 0; i < md_array_size(b.res_inst); ++i) all_in_inst &= b.res_inst[i] >= 0;
    if (all_in_inst && !scope_whole) {
        md_array(str_t) terms = 0;
        md_array(int32_t) unwritten = 0;
        // An instance which is not a molecule of its own, nor whole molecules, leaves no term to be the scope
        const bool ok = whole_terms(&terms, &unwritten, chains, md_array_size(chains), loose_insts, md_array_size(loose_insts), *c, true) &&
            md_array_size(unwritten) == 0;
        md_str_writer_t sb = md_str_writer_create(temp);
        write_or(&sb, terms, md_array_size(terms));
        if (ok && !too_long(&sb, c->limit)) {
            s->text[scope_slot(Scope_Inst)] = take(&sb, temp);
            s->avail |= Scope_Inst;
        }
    }

    // The residues which hold the selection, when it is not the whole of them, written whole
    if (!b.residues_full) {
        Breakdown whole = b;
        whole.res_full = (uint8_t*)0;
        whole.inst_full = (uint8_t*)0;
        for (size_t i = 0; i < md_array_size(b.res); ++i) md_array_push(whole.res_full, (uint8_t)1, temp);
        // An instance is whole when all of its residues hold the selection
        md_array(uint32_t) held = 0;
        md_array_resize(held, md_array_size(b.inst), temp);
        for (size_t i = 0; i < md_array_size(b.inst); ++i) held[i] = 0;
        for (size_t k = 0; k < md_array_size(b.res_slot); ++k) {
            if (b.res_slot[k] >= 0) held[b.res_slot[k]] += 1;
        }
        for (size_t i = 0; i < md_array_size(b.inst); ++i) {
            md_array_push(whole.inst_full, (uint8_t)(held[i] == md_system_instance_residue_count(&sys, b.inst[i])), temp);
        }
        md_str_writer_t sb = md_str_writer_create(temp);
        if (write_structures(&sb, whole, c, true)) {
            s->text[scope_slot(Scope_Res)] = take(&sb, temp);
            s->avail |= Scope_Res;
        }
    }
}

// The predicates which hold for every selected atom
static md_array(Pred) predicates(const Breakdown& b, const Scopes& s, const Ctx& c) {
    const md_system_t& sys = *c.sys;
    md_allocator_i* temp = c.temp;
    md_array(Pred) preds = 0;
    md_str_writer_t sb = md_str_writer_create(temp);

    // Of whole residues
    if (b.residues_full && md_array_size(b.loose) == 0 && md_array_size(b.res)) {
        if (s.class_kind >= 0) {
            Pred p = {};
            p.text = str_from_cstr(class_selector(s.class_kind));
            p.res_level = true;
            p.res_kind = s.class_kind;
            md_array_push(preds, p, temp);
        }
        if (md_array_size(s.resnames)) {
            Pred p = {};
            p.text = s.text[scope_slot(Scope_Resname)];
            p.res_level = true;
            p.res_names = s.resnames;
            md_array_push(preds, p, temp);
        }
    }

    // Of atoms. Of whole residues they say nothing the residues do not say better: name({"OW", "HW1", "HW2"}) is
    // every atom of the waters.
    const size_t num_types = sys.atom.type.count;
    if (!sys.atom.type_idx || num_types == 0 || b.residues_full) return preds;

    // The types of the atoms, and their names in the order the atoms have them (N, CA, C rather than C, CA, N)
    uint8_t* seen = (uint8_t*)md_alloc(temp, num_types);
    MEMSET(seen, 0, num_types);
    md_array(str_t) names = 0;
    bool names_ok = true;
    uint32_t flags_all = ~0u;
    md_bitset_iter_t it = md_bitset_iter(c.sel);
    while (md_bitset_iter_next(&it)) {
        const md_particle_type_idx_t t = sys.atom.type_idx[it.idx];
        if ((size_t)t < num_types && !seen[t]) {
            seen[t] = 1;
            const str_t name = md_particle_type_name(&sys.atom.type, t);
            names_ok &= quotable(name);
            if (names_ok && !contains_str(names, md_array_size(names), name)) {
                md_array_push(names, name, temp);
                names_ok = md_array_size(names) <= MAX_NAMES;
            }
        }
        flags_all &= (uint32_t)md_atom_flags(&sys.atom, it.idx);
    }

    // What the atom flags say the atoms are
    static const struct { uint32_t flag; const char* text; } roles[] = {
        {MD_ATOM_FLAG_BACKBONE,   "backbone()"},
        {MD_ATOM_FLAG_SIDE_CHAIN, "sidechain()"},
        {MD_ATOM_FLAG_NUCLEOBASE, "nucleobase()"},
        {MD_ATOM_FLAG_NUCLEOSIDE, "nucleoside()"},
    };
    if (sys.atom.flags) {
        for (size_t i = 0; i < ARRAY_SIZE(roles); ++i) {
            if (flags_all & roles[i].flag) {
                Pred p = {};
                p.text = str_from_cstr(roles[i].text);
                p.flags = roles[i].flag;
                md_array_push(preds, p, temp);
            }
        }
    }

    // By name, before the elements: of two descriptions as short, the one by name is the more particular
    if (names_ok && md_array_size(names)) {
        Pred p = {};
        p.type_take = (uint8_t*)md_alloc(temp, num_types);
        for (size_t t = 0; t < num_types; ++t) p.type_take[t] = contains_str(names, md_array_size(names), md_particle_type_name(&sys.atom.type, t));
        sb += "name(";
        write_strings(&sb, names, md_array_size(names));
        sb += ')';
        p.text = take(&sb, temp);
        md_array_push(preds, p, temp);
    }

    // By element
    md_array(md_atomic_number_t) zs = 0;
    bool elements_ok = true;
    for (size_t t = 0; t < num_types && elements_ok; ++t) {
        if (!seen[t]) continue;
        const md_atomic_number_t z = sys.atom.type.z ? sys.atom.type.z[t] : 0;
        elements_ok = z > 0;
        bool dup = false;
        for (size_t k = 0; k < md_array_size(zs); ++k) dup |= zs[k] == z;
        if (!dup) md_array_push(zs, z, temp);
    }
    if (elements_ok && md_array_size(zs) && sys.atom.type.z) {
        // Sorted by atomic number, the order in which they are usually spoken of (C, N, O)
        for (size_t i = 1; i < md_array_size(zs); ++i) {
            for (size_t k = i; k > 0 && zs[k] < zs[k - 1]; --k) { md_atomic_number_t tmp = zs[k]; zs[k] = zs[k - 1]; zs[k - 1] = tmp; }
        }
        bool has_h = false;
        for (size_t k = 0; k < md_array_size(zs); ++k) has_h |= zs[k] == 1;
        bool sys_has_h = false;
        for (size_t t = 0; t < num_types; ++t) sys_has_h |= sys.atom.type.z[t] == 1;

        if (md_array_size(zs) <= MAX_ELEMENTS) {
            Pred p = {};
            p.type_take = (uint8_t*)md_alloc(temp, num_types);
            md_array(str_t) syms = 0;
            for (size_t k = 0; k < md_array_size(zs); ++k) md_array_push(syms, md_util_element_symbol(zs[k]), temp);
            for (size_t t = 0; t < num_types; ++t) {
                bool take_t = false;
                for (size_t k = 0; k < md_array_size(zs); ++k) take_t |= sys.atom.type.z[t] == zs[k];
                p.type_take[t] = take_t;
            }
            sb += "element(";
            write_strings(&sb, syms, md_array_size(syms));
            sb += ')';
            p.text = take(&sb, temp);
            md_array_push(preds, p, temp);
        }
        // The heavy atoms: what is left out is the hydrogens
        if (!has_h && sys_has_h && md_array_size(zs) > 1) {
            Pred p = {};
            p.type_take = (uint8_t*)md_alloc(temp, num_types);
            for (size_t t = 0; t < num_types; ++t) p.type_take[t] = sys.atom.type.z[t] != 1;
            p.text = STR_LIT("not element(\"H\")");
            md_array_push(preds, p, temp);
        }
    }

    return preds;
}

// Each predicate within the widest scope where it is the selection
static md_array(str_t) describe(const Breakdown& b, Ctx* c) {
    md_allocator_i* temp = c->temp;
    md_array(str_t) out = 0;

    Scopes s = {};
    scopes_init(&s, b, c);
    md_array(Pred) preds = predicates(b, s, *c);

    // One description in each scope, the shortest: the names of the atoms and their elements often both pick them
    str_t best[5] = {};

    md_str_writer_t sb = md_str_writer_create(temp);
    for (size_t i = 0; i < md_array_size(preds); ++i) {
        const Pred& p = preds[i];
        const uint32_t violated = pred_violations(p, *c, s);
        uint32_t avail = s.avail & ~violated;
        // Whole residues within the residues they are of, or within themselves, say nothing. Within the instances
        // which hold them they say something of several (the HIS of chain A), and of one only that it is the one of
        // its name there, which its residue number says better.
        if (p.res_level) avail &= md_array_size(b.res) > 1 ? (Scope_System | Scope_Inst) : Scope_System;

        for (uint32_t scope = Scope_System; scope <= Scope_Res; scope <<= 1) {
            if (!(avail & scope)) continue;
            sb += p.text;
            if (scope != Scope_System) {
                sb += " and ";
                write_and_operand(&sb, s.text[scope_slot(scope)]);
            }
            str_t& b = best[scope_slot(scope)];
            if (!too_long(&sb, c->limit) && (str_empty(b) || md_str_writer_len(&sb) < b.len)) b = take(&sb, temp);
            md_str_writer_reset(&sb);
            break;
        }
    }
    for (size_t i = 0; i < ARRAY_SIZE(best); ++i) {
        if (!str_empty(best[i])) md_array_push(out, best[i], temp);
    }
    return out;
}

// ### VALIDATION ###

// The sets the expression gives, as a script would have them: one per residue of residue(1:3)
static bool filter_sets(md_array(md_bitset_t)* sets, str_t expr, const md_system_t& sys, md_allocator_i* temp) {
    md_filter_context_t ctx = {};
    ctx.sys = &sys;
    *sets = 0;
    return md_filter_sets(sets, expr, &ctx, temp).success;
}

// The expression picks the atoms of sel and no others, in as many sets as it gives (num_sets)
static bool picks_exactly(str_t expr, md_bitset_t sel, const md_system_t& sys, md_allocator_i* temp, size_t* num_sets = nullptr) {
    md_array(md_bitset_t) sets = 0;
    if (!filter_sets(&sets, expr, sys, temp)) return false;
    const size_t n = md_array_size(sets);
    if (num_sets) *num_sets = n;
    if (n == 0) return md_bitset_empty(sel);
    if (n == 1) return md_bitset_equal(sets[0], sel);
    md_bitmask_t all = md_bitmask_create(0, (uint32_t)sys.atom.count, temp);
    for (size_t i = 0; i < n; ++i) md_bitmask_or(&all, sets[i]);
    return md_bitset_equal(md_bitmask_view(&all), sel);
}

static bool compiles(str_t src, const md_system_t& sys, md_allocator_i* temp) {
    md_script_ir_t* ir = md_script_ir_create(temp);
    md_script_ir_compile_from_source(ir, src, &sys, nullptr);
    const bool ok = md_script_ir_valid(ir) && md_script_ir_num_errors(ir) == 0;
    md_script_ir_free(ir);
    return ok;
}

// Whole glycans (glycan(): the parts all of whose residues are sugars), alone and within the chains which hold them:
// glycan(), glycan() in chain("A"). Whether they pick the selection is for the validation to tell.
static void glycan_terms(md_array(str_t)* out, const Breakdown& b, Ctx* c) {
    const md_system_t& sys = *c->sys;
    if (!b.residues_full || md_array_size(b.loose) || md_array_size(b.res) == 0 || md_system_part_count(&sys) == 0) return;

    md_array(int32_t) chains = 0;
    bool all_in_chain = true;
    int32_t cur_part = -1;
    for (size_t k = 0; k < md_array_size(b.res); ++k) {
        const md_residue_kind_t kind = md_residue_kind(&sys.residue, (size_t)b.res[k]);
        if (kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE || kind == MD_RESIDUE_KIND_WATER ||
            kind == MD_RESIDUE_KIND_ION || kind == MD_RESIDUE_KIND_LIPID) return;
        const int32_t part = md_system_part_find_by_residue_idx(&sys, (size_t)b.res[k]);
        if (part < 0) return;
        if (part != cur_part) {
            cur_part = part;
            const md_urange_t pa = md_system_part_atom_range(&sys, (size_t)part);
            if (!md_bitset_test_all_range(c->sel, pa.beg, pa.end)) return;
        }
        const int32_t chain = res_chain(*c, b.res[k]);
        if (chain < 0) all_in_chain = false;
        else md_array_push(chains, chain, c->temp);
    }
    // All of them, else those of the chains which hold the selection
    if (picks_exactly(STR_LIT("glycan()"), c->sel, sys, c->temp)) {
        md_array_push(*out, STR_LIT("glycan()"), c->temp);
        return;
    }

    const size_t num_chains = sort_unique(chains, md_array_size(chains));
    if (!all_in_chain || num_chains == 0) return;
    md_array(str_t) terms = 0;
    md_array(int32_t) unwritten = 0;
    if (!whole_terms(&terms, &unwritten, chains, num_chains, nullptr, 0, *c, true)) return;
    md_str_writer_t sb = md_str_writer_create(c->temp);
    md_str_writer_t ctx = md_str_writer_create(c->temp);
    write_or(&ctx, terms, md_array_size(terms));
    const bool paren = md_array_size(terms) > 1;
    sb += "glycan() in ";
    if (paren) sb += '(';
    sb += md_str_writer_to_str(&ctx);
    if (paren) sb += ')';
    if (!too_long(&sb, c->limit)) md_array_push(*out, take(&sb, c->temp), c->temp);
}

static void push_suggestion(md_array(ScriptSuggestion)* out, str_t text, str_t hint, md_allocator_i* alloc) {
    for (size_t i = 0; i < md_array_size(*out); ++i) {
        if (str_eq((*out)[i].text, text)) return;
    }
    ScriptSuggestion s = {str_copy(text, alloc), str_copy(hint, alloc)};
    md_array_push(*out, s, alloc);
}

// ### SELECTION ###

md_array(ScriptSuggestion) script_suggest_selection(str_t ident, md_bitset_t sel, const md_system_t& sys, md_allocator_i* alloc) {
    md_array(ScriptSuggestion) out = 0;
    if (md_bitset_empty(sel) || sys.atom.count == 0 || str_empty(ident)) return out;

    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    Ctx c = {};
    ctx_init(&c, sys, sel, temp);
    Breakdown b = {};
    breakdown_init(&b, c);

    md_array(str_t) cand = 0;
    md_str_writer_t sb = md_str_writer_create(temp);

    if (md_bitset_count(sel) == sys.atom.count) {
        md_array_push(cand, STR_LIT("all()"), temp);
    } else {
        // What they are, the shortest first
        md_array(str_t) desc = describe(b, &c);
        glycan_terms(&desc, b, &c);
        for (size_t i = 1; i < md_array_size(desc); ++i) {
            for (size_t k = i; k > 0 && desc[k].len < desc[k - 1].len; --k) { str_t tmp = desc[k]; desc[k] = desc[k - 1]; desc[k - 1] = tmp; }
        }
        for (size_t i = 0; i < md_array_size(desc); ++i) md_array_push(cand, desc[i], temp);
    }

    // Where they are: by id and name, then by position
    if (write_structures(&sb, b, &c, true))  md_array_push(cand, take(&sb, temp), temp);
    md_str_writer_reset(&sb);
    if (write_structures(&sb, b, &c, false)) md_array_push(cand, take(&sb, temp), temp);
    md_str_writer_reset(&sb);
    if (write_atoms(&sb, c)) md_array_push(cand, take(&sb, temp), temp);
    md_str_writer_reset(&sb);

    const uint64_t count = md_bitset_count(sel);
    md_str_writer_t line = md_str_writer_create(temp);
    auto offer = [&](str_t expr, bool listed) {
        size_t num_sets = 0;
        if (!picks_exactly(expr, sel, sys, temp, &num_sets)) return false;
        // An array of selections behaves differently from one (in a context, in rmsd()), so how many it gives is told
        md_str_writer_fmt(&sb, "%llu atom%s", (unsigned long long)count, count == 1 ? "" : "s");
        if (num_sets > 1) md_str_writer_fmt(&sb, ", as an array of %zu selections", num_sets);
        if (listed) sb += ", listed: nothing shorter describes them";
        md_str_writer_fmt(&line, STR_FMT " = " STR_FMT ";", STR_ARG(ident), STR_ARG(expr));
        push_suggestion(&out, md_str_writer_to_str(&line), md_str_writer_to_str(&sb), alloc);
        md_str_writer_reset(&line);
        md_str_writer_reset(&sb);
        return true;
    };

    for (size_t i = 0; i < md_array_size(cand) && md_array_size(out) < MAX_SELECTION_SUGGESTIONS; ++i) {
        bool dup = false;
        for (size_t k = 0; k < i; ++k) dup |= str_eq(cand[k], cand[i]);
        if (!dup) offer(cand[i], false);
    }

    // Nothing describes it: the shortest listing which picks it, if any is short enough to be of use
    if (md_array_size(out) == 0) {
        c.limit = MAX_LISTED_EXPR_LEN;
        c.max_groups = MAX_LISTED_EXPR_LEN / 16;
        md_array(str_t) listed = 0;
        if (write_structures(&sb, b, &c, true))  md_array_push(listed, take(&sb, temp), temp);
        md_str_writer_reset(&sb);
        if (write_structures(&sb, b, &c, false)) md_array_push(listed, take(&sb, temp), temp);
        md_str_writer_reset(&sb);
        if (write_atoms(&sb, c)) md_array_push(listed, take(&sb, temp), temp);
        md_str_writer_reset(&sb);
        for (size_t i = 1; i < md_array_size(listed); ++i) {
            for (size_t k = i; k > 0 && listed[k].len < listed[k - 1].len; --k) { str_t tmp = listed[k]; listed[k] = listed[k - 1]; listed[k - 1] = tmp; }
        }
        for (size_t i = 0; i < md_array_size(listed) && !offer(listed[i], true); ++i) {}
    }

    md_arena_allocator_destroy(temp);
    return out;
}

// ### MEASUREMENT ###

// The name of the atom when no other atom of its residue has it
static bool unique_name_in_res(str_t* out, const md_system_t& sys, int32_t atom, int32_t res) {
    const str_t name = md_atom_name(&sys.atom, atom);
    if (!quotable(name)) return false;
    const md_urange_t r = md_residue_atom_range(&sys.residue, res);
    for (uint32_t i = r.beg; i < r.end; ++i) {
        if ((int32_t)i != atom && str_eq(md_atom_name(&sys.atom, i), name)) return false;
    }
    *out = name;
    return true;
}

// The context which picks the residue alone, as the file numbers it: 'resid(42) in chain("A")' (or of an inferred
// chain, 'resid(42) in chain(2)'), else 'resid(42)'
static bool resid_context(md_str_writer_t* sb, Ctx* c, int32_t res) {
    const md_system_t& sys = *c->sys;
    if (!sys.residue.seq_id) return false;
    const int32_t chain = res_chain(*c, res);
    *sb += "resid(";
    write_int(sb, md_residue_seq_id(&sys.residue, res));
    *sb += ')';
    if (chain_by_seq(c, chain)) {
        *sb += " in ";
        write_chain(sb, *c, chain);
    }
    return true;
}

// The context picks one structure, and within it each argument picks its atom and nothing else
static bool args_pick_atoms(const str_t* args, const int32_t* atoms, size_t count, str_t context, const md_system_t& sys, md_allocator_i* temp) {
    md_array(md_bitset_t) sets = 0;
    if (!filter_sets(&sets, context, sys, temp) || md_array_size(sets) != 1) return false;
    md_str_writer_t sb = md_str_writer_create(temp);
    for (size_t i = 0; i < count; ++i) {
        md_str_writer_fmt(&sb, STR_FMT " in " STR_FMT, STR_ARG(args[i]), STR_ARG(context));
        const bool ok = picks_exactly(md_str_writer_to_str(&sb), md_bitset_range(atoms[i], atoms[i] + 1), sys, temp);
        md_str_writer_reset(&sb);
        if (!ok) return false;
    }
    return true;
}

md_array(ScriptSuggestion) script_suggest_measurement(str_t ident, const int32_t* atoms, size_t count, const md_system_t& sys, md_allocator_i* alloc) {
    md_array(ScriptSuggestion) out = 0;
    if (count < 2 || count > 4 || str_empty(ident)) return out;
    for (size_t i = 0; i < count; ++i) {
        if (atoms[i] < 0 || (size_t)atoms[i] >= sys.atom.count) return out;
        for (size_t k = 0; k < i; ++k) {
            if (atoms[k] == atoms[i]) return out;
        }
    }

    static const char* procs[] = {nullptr, nullptr, "distance", "angle", "dihedral"};
    const char* proc = procs[count];

    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    Ctx c = {};
    ctx_init(&c, sys, md_bitset_empty_set(), temp);

    struct Line { str_t text; str_t hint; };
    md_array(Line) lines = 0;
    md_str_writer_t sb = md_str_writer_create(temp);

    auto write_call = [&](const str_t* args) {
        md_str_writer_fmt(&sb, STR_FMT " = %s(", STR_ARG(ident), proc);
        for (size_t i = 0; i < count; ++i) {
            if (i) sb += ", ";
            sb += args[i];
        }
        sb += ')';
    };
    auto push_line = [&](str_t hint) {
        sb += ';';
        Line l = {take(&sb, temp), hint};
        md_array_push(lines, l, temp);
    };

    // By index in the system: always right, and what the atoms are called does not matter
    str_t args[4] = {};
    for (size_t i = 0; i < count; ++i) {
        md_str_writer_fmt(&sb, "%d", atoms[i] + 1);
        args[i] = take(&sb, temp);
    }
    write_call(args);
    push_line({});

    int32_t atom_res[4] = {-1, -1, -1, -1};
    str_t names[4] = {};
    bool named = true;
    bool same_res = true;
    for (size_t i = 0; i < count; ++i) {
        atom_res[i] = md_residue_find_by_atom_idx(&sys.residue, atoms[i]);
        same_res &= atom_res[i] >= 0 && atom_res[i] == atom_res[0];
        named &= atom_res[i] >= 0 && unique_name_in_res(&names[i], sys, atoms[i], atom_res[i]);
    }

    if (same_res) {
        // Within the residue: the atoms by name, else by their positions in it
        const int32_t res = atom_res[0];
        const md_urange_t r = md_residue_atom_range(&sys.residue, res);
        str_t pick[4] = {};     // what picks each atom within the residue, as a selection
        for (size_t i = 0; i < count; ++i) {
            if (named) {
                sb += "name(";
                write_quoted(&sb, names[i]);
                sb += ')';
                args[i] = take(&sb, temp);
                pick[i] = args[i];
            } else {
                md_str_writer_fmt(&sb, "%u", atoms[i] - r.beg + 1);
                args[i] = take(&sb, temp);
                md_str_writer_fmt(&sb, "atom(%u)", atoms[i] - r.beg + 1);
                pick[i] = take(&sb, temp);
            }
        }

        str_t contexts[2] = {};
        md_str_writer_fmt(&sb, "residue(%d)", res + 1);
        contexts[0] = take(&sb, temp);
        if (resid_context(&sb, &c, res)) contexts[1] = take(&sb, temp);
        md_str_writer_reset(&sb);

        for (size_t k = 0; k < ARRAY_SIZE(contexts); ++k) {
            if (str_empty(contexts[k]) || !args_pick_atoms(pick, atoms, count, contexts[k], sys, temp)) continue;
            write_call(args);
            md_str_writer_fmt(&sb, " in " STR_FMT, STR_ARG(contexts[k]));
            push_line({});
        }

        // The same measurement in every residue of the name: each of them has to have each of the atoms, once
        const str_t resname = md_residue_name(&sys.residue, res);
        if (named && quotable(resname)) {
            sb += "resname(";
            write_quoted(&sb, resname);
            sb += ')';
            const str_t context = take(&sb, temp);
            md_array(md_bitset_t) residues = 0;
            bool ok = filter_sets(&residues, context, sys, temp) && md_array_size(residues) > 0;
            const size_t num_res = md_array_size(residues);
            for (size_t i = 0; i < count && ok; ++i) {
                md_str_writer_fmt(&sb, STR_FMT " in " STR_FMT, STR_ARG(pick[i]), STR_ARG(context));
                md_array(md_bitset_t) sets = 0;
                ok = filter_sets(&sets, md_str_writer_to_str(&sb), sys, temp) && md_array_size(sets) == num_res;
                for (size_t k = 0; k < md_array_size(sets) && ok; ++k) ok = md_bitset_count(sets[k]) == 1;
                md_str_writer_reset(&sb);
            }
            if (ok) {
                str_t hint = {};
                if (num_res > 1) {
                    md_str_writer_fmt(&sb, "One value for each of the %zu " STR_FMT, num_res, STR_ARG(resname));
                    hint = take(&sb, temp);
                }
                write_call(args);
                md_str_writer_fmt(&sb, " in " STR_FMT, STR_ARG(context));
                push_line(hint);
            }
        }
    } else if (named) {
        // Each atom by name within its residue, as the system counts them and as the file numbers them
        for (int form = 0; form < 2; ++form) {
            bool ok = true;
            for (size_t i = 0; i < count && ok; ++i) {
                sb += "name(";
                write_quoted(&sb, names[i]);
                sb += ") in ";
                if (form == 0) md_str_writer_fmt(&sb, "residue(%d)", atom_res[i] + 1);
                else ok = resid_context(&sb, &c, atom_res[i]);
                args[i] = take(&sb, temp);
                ok = ok && picks_exactly(args[i], md_bitset_range(atoms[i], atoms[i] + 1), sys, temp);
            }
            if (!ok) continue;
            write_call(args);
            push_line({});
        }
    }

    for (size_t i = 0; i < md_array_size(lines); ++i) {
        if (compiles(lines[i].text, sys, temp)) push_suggestion(&out, lines[i].text, lines[i].hint, alloc);
    }

    md_arena_allocator_destroy(temp);
    return out;
}
