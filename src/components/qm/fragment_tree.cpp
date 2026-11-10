#include "fragment_tree.h"

#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

namespace {

struct Builder {
    const md_system_t* sys;
    md_allocator_i*    alloc;
    md_allocator_i*    temp;
    const int32_t*     atom_res;       // [atom count] the residue of each atom, -1 for none

    md_array(uint32_t) parent = nullptr;
    md_array(uint8_t)  kind = nullptr;
    md_array(uint8_t)  frag_kind = nullptr;
    md_array(char)     name = nullptr;
    md_array(uint32_t) beg = nullptr;
    md_array(uint32_t) end = nullptr;
    md_array(int32_t)  atoms = nullptr;
    uint32_t*          atom_leaf = nullptr;

    uint32_t push(uint32_t par, FragTreeKind k, md_fragment_kind_t fk, const char* nm) {
        const uint32_t idx = (uint32_t)md_array_size(parent);
        md_array_push(parent, par, alloc);
        md_array_push(kind, (uint8_t)k, alloc);
        md_array_push(frag_kind, (uint8_t)fk, alloc);
        const size_t base = md_array_size(name);
        md_array_resize(name, base + FRAG_TREE_NAME_CAP, alloc);
        snprintf(name + base, FRAG_TREE_NAME_CAP, "%s", nm);
        md_array_push(beg, (uint32_t)md_array_size(atoms), alloc);
        md_array_push(end, (uint32_t)md_array_size(atoms), alloc);
        return idx;
    }

    void add_atom(uint32_t leaf, int32_t a) {
        md_array_push(atoms, a, alloc);
        atom_leaf[a] = leaf;
    }

    // Ends the ranges of a node and all of its ancestors at the current end
    void close(uint32_t node) {
        for (uint32_t f = node; f != MD_FRAGMENT_NONE; f = parent[f]) end[f] = (uint32_t)md_array_size(atoms);
    }
};

// Composition of a set of atoms, carbon and hydrogen first and the rest alphabetically (Hill)
void hill_formula(char* buf, size_t cap, const md_system_t* sys, const int32_t* atoms, size_t count) {
    uint32_t num[128] = {};
    for (size_t i = 0; i < count; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, (size_t)atoms[i]);
        if (z < 128) num[z]++;
    }
    size_t len = 0;
    buf[0] = '\0';
    auto append = [&](int z) {
        const md_str_t sym = md_util_element_symbol((md_element_t)z);
        const int w = num[z] > 1 ? snprintf(buf + len, cap - len, "%.*s%u", (int)sym.len, sym.ptr, num[z])
                                 : snprintf(buf + len, cap - len, "%.*s", (int)sym.len, sym.ptr);
        if (w > 0) len = MD_MIN(cap - 1, len + (size_t)w);
    };
    if (num[MD_Z_C]) append(MD_Z_C);
    if (num[MD_Z_H]) append(MD_Z_H);
    bool done[128] = {};
    done[0] = done[MD_Z_C] = done[MD_Z_H] = true;
    for (;;) {
        int best = -1;
        md_str_t best_sym = {};
        for (int z = 1; z < 128; ++z) {
            if (done[z] || !num[z]) continue;
            const md_str_t sym = md_util_element_symbol((md_element_t)z);
            if (best < 0 || md_str_cmp_lex(sym, best_sym) < 0) { best = z; best_sym = sym; }
        }
        if (best < 0) break;
        done[best] = true;
        append(best);
    }
}

// "ALA 42", or "ALA 42 (chain B)" where the chain says which: in a molecule of several chains, where the numbers repeat
void residue_label(char* buf, size_t cap, const md_system_t* sys, int32_t r, bool with_chain = false) {
    const md_str_t rn = md_system_residue_name(sys, (size_t)r);
    const int seq = (int)md_system_residue_seq_id(sys, (size_t)r);
    const int64_t chain = with_chain ? (int64_t)md_system_chain_find_by_residue_idx(sys, (size_t)r) : -1;
    if (chain >= 0) {
        const md_str_t id = md_system_chain_id(sys, (size_t)chain);
        snprintf(buf, cap, "%.*s %d (chain %.*s)", (int)rn.len, rn.ptr, seq, (int)id.len, id.ptr);
    } else {
        snprintf(buf, cap, "%.*s %d", (int)rn.len, rn.ptr, seq);
    }
}

int64_t residue_chain(const md_system_t* sys, int32_t r) {
    return (int64_t)md_system_chain_find_by_residue_idx(sys, (size_t)r);
}

// Numbers the names repeated among the children [first, last) of one parent, in order
void number_repeats(Builder& b, const uint32_t* nodes, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        char* ni = b.name + (size_t)nodes[i] * FRAG_TREE_NAME_CAP;
        size_t total = 0;
        for (size_t j = 0; j < count; ++j) total += strcmp(ni, b.name + (size_t)nodes[j] * FRAG_TREE_NAME_CAP) == 0;
        if (total < 2) continue;
        char base[FRAG_TREE_NAME_CAP];
        snprintf(base, sizeof(base), "%s", ni);
        uint32_t k = 0;
        for (size_t j = i; j < count; ++j) {
            char* nj = b.name + (size_t)nodes[j] * FRAG_TREE_NAME_CAP;
            if (strcmp(nj, base) != 0) continue;
            char tmp[FRAG_TREE_NAME_CAP + 12];
            snprintf(tmp, sizeof(tmp), "%s %u", base, ++k);
            snprintf(nj, FRAG_TREE_NAME_CAP, "%s", tmp);
        }
    }
}

// The atoms of fragment f in residue r (r < 0: all of them), in its slot order
void fragment_atoms(md_array(int32_t)* out, const Builder& b, md_fragment_idx_t f, int32_t r) {
    md_array_shrink(*out, 0);
    for (md_fragment_atom_iter_t it = md_fragment_atom_iter(b.sys, f); md_fragment_atom_iter_next(&it);) {
        if (r < 0 || (b.atom_res && b.atom_res[it.atom] == r)) md_array_push(*out, it.atom, b.temp);
    }
}

// The residues of molecule m, ascending
md_array(int32_t) molecule_residues(const Builder& b, size_t m) {
    const md_molecule_t mol = md_system_molecule(b.sys, m);
    md_array(int32_t) residues = nullptr;
    if (b.atom_res) {
        for (size_t k = 0; k < mol.count; ++k) {
            const int32_t r = b.atom_res[mol.atom_idx[k]];
            if (r < 0) continue;
            bool have = false;
            for (size_t j = 0; j < md_array_size(residues) && !have; ++j) have = residues[j] == r;
            if (!have) md_array_push(residues, r, b.temp);
        }
        for (size_t i = 1; i < md_array_size(residues); ++i) {
            const int32_t v = residues[i];
            size_t j = i;
            while (j > 0 && residues[j - 1] > v) { residues[j] = residues[j - 1]; --j; }
            residues[j] = v;
        }
    }
    return residues;
}

// What molecule m is called on its own: its fragment when it is one, else its residue, "protein", "nucleic acid",
// or its formula
void molecule_name(char* nm, size_t cap, const Builder& b, size_t m, const int32_t* residues, size_t num_residues) {
    const md_system_t* sys = b.sys;
    const uint32_t f_beg = sys->fragment.frag_offset[m], f_end = sys->fragment.frag_offset[m + 1];
    if (f_end - f_beg == 1) {
        const md_str_t fn = md_system_fragment_name(sys, f_beg);
        snprintf(nm, cap, "%.*s", (int)fn.len, fn.ptr);
    } else if (num_residues == 1) {
        residue_label(nm, cap, sys, residues[0]);
    } else {
        bool protein = false, nucleic = false;
        for (size_t j = 0; j < num_residues; ++j) {
            const md_residue_kind_t rk = md_system_residue_kind(sys, (size_t)residues[j]);
            protein |= rk == MD_RESIDUE_KIND_AMINO_ACID;
            nucleic |= rk == MD_RESIDUE_KIND_NUCLEOTIDE;
        }
        const md_molecule_t mol = md_system_molecule(sys, m);
        if (protein)      snprintf(nm, cap, "protein");
        else if (nucleic) snprintf(nm, cap, "nucleic acid");
        else              hill_formula(nm, cap, sys, mol.atom_idx, mol.count);
    }
}

// One molecule, under 'par'; 'name_override' names it within a set
void emit_molecule(Builder& b, size_t m, uint32_t par, const char* name_override) {
    const md_system_t* sys = b.sys;
    const md_fragment_data_t& F = sys->fragment;
    const uint32_t f_beg = F.frag_offset[m], f_end = F.frag_offset[m + 1];
    md_array(int32_t) scratch = nullptr;
    if (f_end == f_beg) return;

    md_array(int32_t) residues = molecule_residues(b, m);
    char nm[FRAG_TREE_NAME_CAP];
    if (name_override) snprintf(nm, sizeof(nm), "%s", name_override);
    else               molecule_name(nm, sizeof(nm), b, m, residues, md_array_size(residues));

    if (f_end - f_beg == 1) {
        const uint32_t leaf = b.push(par, FRAG_TREE_FRAGMENT, md_system_fragment_kind(sys, f_beg), nm);
        fragment_atoms(&scratch, b, f_beg, -1);
        for (size_t i = 0; i < md_array_size(scratch); ++i) b.add_atom(leaf, scratch[i]);
        b.close(leaf);
        return;
    }

    const uint32_t mnode = b.push(par, FRAG_TREE_MOLECULE, MD_FRAGMENT_KIND_COUNT, nm);
    md_array(uint32_t) kids = nullptr;

    if (md_array_size(residues) < 2) {
        for (uint32_t f = f_beg; f < f_end; ++f) {
            const md_str_t fn = md_system_fragment_name(sys, f);
            snprintf(nm, sizeof(nm), "%.*s", (int)fn.len, fn.ptr);
            const uint32_t leaf = b.push(mnode, FRAG_TREE_FRAGMENT, md_system_fragment_kind(sys, f), nm);
            md_array_push(kids, leaf, b.temp);
            fragment_atoms(&scratch, b, f, -1);
            for (size_t i = 0; i < md_array_size(scratch); ++i) b.add_atom(leaf, scratch[i]);
            b.close(leaf);
        }
        number_repeats(b, kids, md_array_size(kids));
        return;
    }

    // A molecule of residues: each residue with the fragments it is the first of. A fragment of several residues (a
    // peptide group, a phosphate, a disulfide) is whole, under the first of them, its name saying the next ("peptide,
    // with ALA 4"): every atom is in one leaf, and the leaves are the chemistry, at the cost of a residue's node
    // holding the N of the next residue and not its own.
    md_array(md_fragment_idx_t) pieces = nullptr;
    md_array(md_fragment_idx_t) owned = nullptr;
    md_residue_idx_t of[4];
    // Of several chains (chains joined by a disulfide), the residues say their chain, as their numbers repeat
    bool chains = false;
    for (size_t j = 1; j < md_array_size(residues) && !chains; ++j) chains = residue_chain(sys, residues[j]) != residue_chain(sys, residues[0]);
    for (size_t j = 0; j < md_array_size(residues); ++j) {
        const int32_t r = residues[j];
        const size_t num = md_system_residue_fragments(nullptr, 0, sys, (size_t)r);
        md_array_resize(pieces, num, b.temp);
        md_system_residue_fragments(pieces, num, sys, (size_t)r);
        md_array_shrink(owned, 0);
        for (size_t p = 0; p < num; ++p) {
            const size_t no = md_system_fragment_residues(of, MD_COUNTOF(of), sys, pieces[p]);
            if (no > 0 && of[0] == r) md_array_push(owned, pieces[p], b.temp);
        }
        const size_t num_owned = md_array_size(owned);
        if (num_owned == 0) continue;
        // In the residue's own order: by its first atom, those reaching on to a later residue (the peptide group to
        // the next) last. Of a protein: the N-terminus, the alpha carbon, the side chain, the peptide group.
        {
            uint64_t* key = (uint64_t*)md_alloc(b.temp, sizeof(uint64_t) * num_owned);
            for (size_t p = 0; p < num_owned; ++p) {
                uint32_t first_atom = UINT32_MAX;
                for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, owned[p]); md_fragment_atom_iter_next(&it);) {
                    first_atom = MD_MIN(first_atom, (uint32_t)it.atom);
                }
                const bool reaches_on = md_system_fragment_residues(of, MD_COUNTOF(of), sys, owned[p]) > 1;
                key[p] = ((uint64_t)reaches_on << 32) | first_atom;
            }
            for (size_t p = 1; p < num_owned; ++p) {
                const uint64_t k = key[p];
                const md_fragment_idx_t f = owned[p];
                size_t q = p;
                for (; q > 0 && key[q - 1] > k; --q) { key[q] = key[q - 1]; owned[q] = owned[q - 1]; }
                key[q] = k;
                owned[q] = f;
            }
        }

        char rl[FRAG_TREE_NAME_CAP];
        residue_label(rl, sizeof(rl), sys, r, chains);
        if (num_owned == 1 && num == 1) {
            // The residue is one fragment, of it alone
            const uint32_t leaf = b.push(mnode, FRAG_TREE_FRAGMENT, md_system_fragment_kind(sys, owned[0]), rl);
            md_array_push(kids, leaf, b.temp);
            fragment_atoms(&scratch, b, owned[0], -1);
            for (size_t i = 0; i < md_array_size(scratch); ++i) b.add_atom(leaf, scratch[i]);
            b.close(leaf);
            continue;
        }
        const uint32_t rnode = b.push(mnode, FRAG_TREE_RESIDUE, MD_FRAGMENT_KIND_COUNT, rl);
        md_array_push(kids, rnode, b.temp);
        md_array(uint32_t) rkids = nullptr;
        for (size_t p = 0; p < num_owned; ++p) {
            const md_fragment_idx_t f = owned[p];
            const md_str_t fn = md_system_fragment_name(sys, f);
            const size_t no = md_system_fragment_residues(of, MD_COUNTOF(of), sys, f);
            if (no > 1) {
                char ol[FRAG_TREE_NAME_CAP];
                residue_label(ol, sizeof(ol), sys, of[1], chains && residue_chain(sys, of[1]) != residue_chain(sys, r));
                snprintf(nm, sizeof(nm), "%.*s, with %s", (int)fn.len, fn.ptr, ol);
            } else {
                snprintf(nm, sizeof(nm), "%.*s", (int)fn.len, fn.ptr);
            }
            const uint32_t leaf = b.push(rnode, FRAG_TREE_FRAGMENT, md_system_fragment_kind(sys, f), nm);
            md_array_push(rkids, leaf, b.temp);
            fragment_atoms(&scratch, b, f, -1);
            for (size_t i = 0; i < md_array_size(scratch); ++i) b.add_atom(leaf, scratch[i]);
            b.close(leaf);
        }
        number_repeats(b, rkids, md_array_size(rkids));
    }
    // Fragments of atoms in no residue, after
    for (uint32_t f = f_beg; f < f_end; ++f) {
        if (md_system_fragment_residues(of, MD_COUNTOF(of), sys, f) > 0) continue;
        const md_str_t fn = md_system_fragment_name(sys, f);
        snprintf(nm, sizeof(nm), "%.*s", (int)fn.len, fn.ptr);
        const uint32_t leaf = b.push(mnode, FRAG_TREE_FRAGMENT, md_system_fragment_kind(sys, f), nm);
        md_array_push(kids, leaf, b.temp);
        fragment_atoms(&scratch, b, f, -1);
        for (size_t i = 0; i < md_array_size(scratch); ++i) b.add_atom(leaf, scratch[i]);
        b.close(leaf);
    }
    number_repeats(b, kids, md_array_size(kids));
}

}  // namespace

bool frag_tree_build(FragmentTree* out, const md_system_t* sys, const bool* include, md_allocator_i* alloc) {
    if (!out || !sys || !alloc) return false;
    *out = FragmentTree{};
    const size_t n = md_atom_count(&sys->atom);
    out->num_atoms = n;
    out->atom_fragment = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (n ? n : 1));
    for (size_t i = 0; i < n; ++i) out->atom_fragment[i] = MD_FRAGMENT_NONE;
    const size_t num_mols = md_system_molecule_count(sys);
    if (n == 0 || num_mols == 0 || md_array_size(sys->fragment.frag_offset) != num_mols + 1) return true;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    Builder b = {};
    b.sys = sys;
    b.alloc = alloc;
    b.temp = md_temp_allocator(temp);
    b.atom_leaf = out->atom_fragment;

    int32_t* atom_res = nullptr;
    if (sys->residue.count && sys->residue.atom_offset) {
        atom_res = md_temp_alloc_array(temp, int32_t, n);
        for (size_t i = 0; i < n; ++i) atom_res[i] = -1;
        for (size_t r = 0; r < sys->residue.count; ++r) {
            const md_urange_t range = md_system_residue_atom_range(sys, r);
            for (uint32_t a = range.beg; a < range.end && a < n; ++a) atom_res[a] = (int32_t)r;
        }
    }
    b.atom_res = atom_res;

    // The molecules shown, and how many of each template
    bool* shown = md_temp_alloc_zero_array(temp, bool, num_mols);
    uint32_t* tmpl_count = md_temp_alloc_zero_array(temp, uint32_t, sys->fragment.template_count ? sys->fragment.template_count : 1);
    for (size_t m = 0; m < num_mols; ++m) {
        const md_molecule_t mol = md_system_molecule(sys, m);
        bool any = include == nullptr;
        for (size_t k = 0; k < mol.count && !any; ++k) any = include[mol.atom_idx[k]];
        if (!any) continue;
        shown[m] = true;
        tmpl_count[sys->fragment.mol_template[m]]++;
    }

    bool* done = md_temp_alloc_zero_array(temp, bool, num_mols);
    md_array(uint32_t) roots = nullptr;
    for (size_t m = 0; m < num_mols; ++m) {
        if (!shown[m] || done[m]) continue;
        const uint32_t t = sys->fragment.mol_template[m];
        if (tmpl_count[t] < 2) {
            const uint32_t first = (uint32_t)md_array_size(b.parent);
            emit_molecule(b, m, MD_FRAGMENT_NONE, nullptr);
            if (md_array_size(b.parent) > first) md_array_push(roots, first, b.temp);
            done[m] = true;
            continue;
        }
        // A set, named as its first member would be on its own
        if (sys->fragment.frag_offset[m + 1] == sys->fragment.frag_offset[m]) { done[m] = true; continue; }
        char member[FRAG_TREE_NAME_CAP];
        {
            md_array(int32_t) residues = molecule_residues(b, m);
            molecule_name(member, sizeof(member), b, m, residues, md_array_size(residues));
        }
        char set_name[FRAG_TREE_NAME_CAP];
        snprintf(set_name, sizeof(set_name), "%s \xc3\x97%u", member, tmpl_count[t]);
        const uint32_t snode = b.push(MD_FRAGMENT_NONE, FRAG_TREE_SET, MD_FRAGMENT_KIND_COUNT, set_name);
        md_array_push(roots, snode, b.temp);
        uint32_t k = 0;
        for (size_t p = m; p < num_mols; ++p) {
            if (!shown[p] || done[p] || sys->fragment.mol_template[p] != t) continue;
            char nm[FRAG_TREE_NAME_CAP + 12];
            snprintf(nm, sizeof(nm), "%s %u", member, ++k);
            char clipped[FRAG_TREE_NAME_CAP];
            snprintf(clipped, sizeof(clipped), "%s", nm);
            emit_molecule(b, p, snode, clipped);
            done[p] = true;
        }
    }
    number_repeats(b, roots, md_array_size(roots));

    out->count        = md_array_size(b.parent);
    out->parent       = b.parent;
    out->kind         = b.kind;
    out->frag_kind    = b.frag_kind;
    out->name         = b.name;
    out->atom_beg     = b.beg;
    out->atom_end     = b.end;
    out->atom_idx     = b.atoms;
    out->num_atom_idx = md_array_size(b.atoms);
    md_temp_end(temp);
    return true;
}
