#include "utest.h"

#include <script_suggestions.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_xyz.h>
#include <md_tpr.h>
#include <md_filter.h>
#include <md_script.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>

#include <string.h>
#include <initializer_list>     // The loops over braced lists: libc++ does not bring it in through the others

/* What the Script menu of the context menu offers for the selection. The promise is that a line does what it says once
 * it is in the script, so every test checks every line it gets by evaluating it again, apart from the generator: the
 * expression of a selection picks the selected atoms and no others, a measurement compiles. On top of that, the forms
 * which say what the atoms are have to be the ones offered where they apply. */

namespace {

struct Loaded {
    md_allocator_i* alloc = nullptr;
    md_system_t sys = {};
};

bool load(Loaded* l, const char* file) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s", VIAMD_TEST_DATA_DIR, file);
    l->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    l->sys = { .alloc = l->alloc };
    md_system_state_t st = { .alloc = l->alloc };
    const str_t p = str_from_cstr(path);
    bool ok = false;
    if (strstr(file, ".tpr"))      return md_tpr_system_init_from_file(&l->sys, &st, p) && md_util_system_infer(&l->sys, &st, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT);
    if (strstr(file, ".cif"))      ok = md_mmcif_system_init_from_file(&l->sys, &st, p);
    else if (strstr(file, ".gro")) ok = md_gro_system_init_from_file(&l->sys, &st, p);
    else if (strstr(file, ".xyz")) ok = md_xyz_system_init_from_file(&l->sys, &st, p, MD_XYZ_OPTION_NONE);
    else                           ok = md_pdb_system_init_from_file(&l->sys, &st, p, MD_PDB_OPTION_NONE);
    return ok && md_util_system_infer(&l->sys, &st, MD_UTIL_INFER_ALL);
}

void unload(Loaded* l) {
    md_arena_allocator_destroy(l->alloc);
}

md_bitset_t select(const Loaded& l, const char* expr) {
    md_filter_context_t ctx = {};
    ctx.sys = &l.sys;
    md_bitset_t sel = {};
    if (!md_filter(&sel, str_from_cstr(expr), &ctx, l.alloc).success) return {};
    return sel;
}

// The one atom the expression picks, -1 when it is not one
int32_t atom(const Loaded& l, const char* expr) {
    const md_bitset_t s = select(l, expr);
    return md_bitset_count(s) == 1 ? (int32_t)s.beg : -1;
}

bool compiles(const Loaded& l, str_t src) {
    md_script_ir_t* ir = md_script_ir_create(l.alloc);
    md_script_ir_compile_from_source(ir, src, &l.sys, nullptr);
    const bool ok = md_script_ir_valid(ir);
    if (!ok) md_script_ir_write_diagnostics(md_writer_stdout(), ir);
    md_script_ir_free(ir);
    return ok;
}

// Nothing which cannot be typed: a name padded with zeros would cut the line short where it is printed
bool printable(str_t s) {
    for (size_t i = 0; i < s.len; ++i) {
        const unsigned char c = (unsigned char)s.ptr[i];
        if (c < 0x20 || c == 0x7F) return false;
    }
    return true;
}

// "sel1 = <expr>;" -> <expr>
str_t rhs(str_t line, str_t ident) {
    const size_t pre = ident.len + 3;
    if (line.len <= pre + 1 || !str_eq(str_substr(line, 0, ident.len), ident) || !str_eq(str_substr(line, ident.len, 3), STR_LIT(" = ")) ||
        line.ptr[line.len - 1] != ';') {
        return {};
    }
    return str_substr(line, pre, line.len - pre - 1);
}

bool contains(const md_array(ScriptSuggestion) s, const char* text) {
    for (size_t i = 0; i < md_array_size(s); ++i) {
        if (str_eq_cstr(s[i].text, text)) return true;
    }
    return false;
}

bool any_contains(const md_array(ScriptSuggestion) s, const char* needle) {
    size_t loc;
    for (size_t i = 0; i < md_array_size(s); ++i) {
        if (str_find_str(&loc, s[i].text, str_from_cstr(needle))) return true;
    }
    return false;
}

void print(const md_array(ScriptSuggestion) s) {
    for (size_t i = 0; i < md_array_size(s); ++i) {
        printf("    " STR_FMT "    [" STR_FMT "]\n", STR_ARG(s[i].text), STR_ARG(s[i].hint));
    }
}

// Every line assigns to the name, compiles, and picks exactly the selection
bool all_pick(const Loaded& l, md_bitset_t sel, const md_array(ScriptSuggestion) s) {
    const str_t ident = STR_LIT("sel1");
    for (size_t i = 0; i < md_array_size(s); ++i) {
        const str_t expr = rhs(s[i].text, ident);
        md_bitset_t got = {};
        md_filter_context_t ctx = {};
        ctx.sys = &l.sys;
        const bool ok = !str_empty(expr) && printable(s[i].text) && compiles(l, s[i].text) &&
            md_filter(&got, expr, &ctx, l.alloc).success && md_bitset_equal(got, sel);
        if (!ok) {
            printf("  does not pick the selection: " STR_FMT "\n", STR_ARG(s[i].text));
            return false;
        }
    }
    return true;
}

md_array(ScriptSuggestion) suggest(const Loaded& l, md_bitset_t sel) {
    return script_suggest_selection(STR_LIT("sel1"), sel, l.sys, l.alloc);
}

md_array(ScriptSuggestion) measure(const Loaded& l, std::initializer_list<const char*> picks) {
    int32_t atoms[8];
    size_t n = 0;
    for (const char* p : picks) {
        if (n < ARRAY_SIZE(atoms)) atoms[n++] = atom(l, p);
    }
    return script_suggest_measurement(STR_LIT("m1"), atoms, n, l.sys, l.alloc);
}

}  // namespace

// The molecules which are not chains (a ligand, a water) are no chain: chain() of their position does not pick them,
// or picks another chain. This is what the menu used to write for a whole ligand.
UTEST(script_suggestions, molecule_which_is_not_a_chain) {
    Loaded l;
    ASSERT_TRUE(load(&l, "npt.gro"));
    const md_bitset_t sel = select(l, "residue(resname(\"LIG\"))");   // the ligand, after the protein
    ASSERT_EQ(10u, (unsigned)md_bitset_count(sel));
    md_array(ScriptSuggestion) s = suggest(l, sel);
    print(s);
    EXPECT_TRUE(all_pick(l, sel, s));
    EXPECT_FALSE(any_contains(s, "chain("));
    EXPECT_TRUE(contains(s, "sel1 = molecule(2);"));
    EXPECT_TRUE(contains(s, "sel1 = resname(\"LIG\");"));

    // Waters likewise, the protein with them. A GRO file gives no chain ids: the protein is a chain by its position
    const md_bitset_t mixed = select(l, "protein() or residue(resname(\"SOL\"))[1:5]");
    md_array(ScriptSuggestion) m = suggest(l, mixed);
    print(m);
    EXPECT_TRUE(all_pick(l, mixed, m));
    EXPECT_TRUE(contains(m, "sel1 = chain(1) or molecule(3:7);"));
    EXPECT_FALSE(any_contains(m, "chain(\""));
    unload(&l);
}

// A chain is all that the data gives its id: in an mmCIF file the protein of chain A and the ligands, ions and waters
// the author filed under A. The protein alone is the protein of the chain.
UTEST(script_suggestions, chains_by_author_id) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1fez.cif"));
    const md_bitset_t a = select(l, "chain('A')");
    md_array(ScriptSuggestion) s = suggest(l, a);
    print(s);
    EXPECT_TRUE(all_pick(l, a, s));
    ASSERT_GT(md_array_size(s), 0u);
    EXPECT_TRUE(str_eq_cstr(s[0].text, "sel1 = chain(\"A\");"));

    const md_bitset_t prot = select(l, "protein() and chain('A')");
    md_array(ScriptSuggestion) t = suggest(l, prot);
    print(t);
    EXPECT_TRUE(all_pick(l, prot, t));
    EXPECT_TRUE(contains(t, "sel1 = protein() and chain(\"A\");"));

    // Residues by the numbers of the chain
    const md_bitset_t res = select(l, "resid(10:20) in chain('B')");
    md_array(ScriptSuggestion) u = suggest(l, res);
    print(u);
    EXPECT_TRUE(all_pick(l, res, u));
    EXPECT_TRUE(contains(u, "sel1 = resid(10:20) in chain(\"B\");"));

    // The molecule of the protein of chain B (what the viewport selects at Molecule granularity): the Mg it
    // coordinates is not of it
    const md_bitset_t mol = select(l, "molecule(atom(2500))");
    md_array(ScriptSuggestion) v = suggest(l, mol);
    print(v);
    EXPECT_TRUE(all_pick(l, mol, v));
    EXPECT_TRUE(contains(v, "sel1 = protein() and chain(\"B\");"));
    EXPECT_TRUE(contains(v, "sel1 = molecule(2);"));

    // The Mg alone, a molecule of its own: by the ids of the file, and by its position among the molecules
    const md_bitset_t mg = select(l, "resname('MG') and chain('B')");
    ASSERT_EQ(1u, (unsigned)md_bitset_count(mg));
    md_array(ScriptSuggestion) w = suggest(l, mg);
    print(w);
    EXPECT_TRUE(all_pick(l, mg, w));
    EXPECT_TRUE(contains(w, "sel1 = resid(802) in chain(\"B\");"));
    EXPECT_TRUE(any_contains(w, "molecule("));
    EXPECT_FALSE(any_contains(w, "instance("));
    unload(&l);
}

// A DNA duplex which martinize-dna writes as one molecule type: one molecule, as GROMACS has it, of two chains, a strand
// each, which the topology gives no ids
UTEST(script_suggestions, strands_of_one_molecule) {
    Loaded l;
    ASSERT_TRUE(load(&l, "martini/martini2_dna.tpr"));
    const md_bitset_t strand = select(l, "chain(2)");
    md_array(ScriptSuggestion) s = suggest(l, strand);
    print(s);
    EXPECT_TRUE(all_pick(l, strand, s));
    EXPECT_TRUE(contains(s, "sel1 = chain(2);"));
    EXPECT_FALSE(any_contains(s, "molecule("));

    const md_bitset_t duplex = select(l, "molecule(1)");
    EXPECT_EQ(md_bitset_count(duplex), md_bitset_count(select(l, "chain(1:2)")));
    md_array(ScriptSuggestion) t = suggest(l, duplex);
    print(t);
    EXPECT_TRUE(all_pick(l, duplex, t));
    EXPECT_TRUE(contains(t, "sel1 = chain(1:2);") || contains(t, "sel1 = nucleic();"));
    unload(&l);
}

// A glycan bonded to a protein is of its molecule. Part granularity, part() and glycan() take the glycan alone: the
// threonine of chains A and B, each with a tree of three sugars (CHARMM-GUI names) on its side chain.
UTEST(script_suggestions, glycans) {
    static const char pdb[] =
        "CRYST1  246.108  246.108  398.534  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  N   THR A 323     220.790 133.490 398.290  0.00  0.00           N\n"
        "ATOM      2  CA  THR A 323     221.130 133.610   1.180  0.00  0.00           C\n"
        "ATOM      3  CB  THR A 323     221.340 132.260   1.860  0.00  0.00           C\n"
        "ATOM      4  OG1 THR A 323     221.250 132.440   3.290  0.00  0.00           O\n"
        "ATOM      5  CG2 THR A 323     222.640 131.430   1.500  0.00  0.00           C\n"
        "ATOM      6  C   THR A 323     222.370 134.490   1.350  0.00  0.00           C\n"
        "ATOM      7  O   THR A 323     222.780 134.850   2.480  0.00  0.00           O\n"
        "HETATM    8  C1  AGA A   1     219.860 132.350   3.670  0.00  0.00           C\n"
        "HETATM    9  C5  AGA A   1     220.490 130.470   4.970  0.00  0.00           C\n"
        "HETATM   10  O5  AGA A   1     219.510 131.000   4.050  0.00  0.00           O\n"
        "HETATM   11  C2  AGA A   1     219.580 133.160   4.900  0.00  0.00           C\n"
        "HETATM   12  N   AGA A   1     219.810 134.540   4.650  0.00  0.00           N\n"
        "HETATM   13  C   AGA A   1     218.940 135.570   4.510  0.00  0.00           C\n"
        "HETATM   14  O   AGA A   1     217.740 135.420   4.650  0.00  0.00           O\n"
        "HETATM   15  CT  AGA A   1     219.560 136.840   4.100  0.00  0.00           C\n"
        "HETATM   16  C3  AGA A   1     220.440 132.730   6.060  0.00  0.00           C\n"
        "HETATM   17  O3  AGA A   1     220.000 133.320   7.220  0.00  0.00           O\n"
        "HETATM   18  C4  AGA A   1     220.340 131.260   6.290  0.00  0.00           C\n"
        "HETATM   19  O4  AGA A   1     219.110 130.940   6.950  0.00  0.00           O\n"
        "HETATM   20  C6  AGA A   1     220.330 128.940   5.010  0.00  0.00           C\n"
        "HETATM   21  O6  AGA A   1     221.420 128.290   5.660  0.00  0.00           O\n"
        "HETATM   22  C1  BGA A   2     220.980 133.950   8.020  0.00  0.00           C\n"
        "HETATM   23  C5  BGA A   2     222.980 133.600   9.080  0.00  0.00           C\n"
        "HETATM   24  O5  BGA A   2     221.830 132.920   8.590  0.00  0.00           O\n"
        "HETATM   25  C2  BGA A   2     220.370 134.620   9.240  0.00  0.00           C\n"
        "HETATM   26  O2  BGA A   2     219.390 135.570   8.770  0.00  0.00           O\n"
        "HETATM   27  C3  BGA A   2     221.470 135.380  10.000  0.00  0.00           C\n"
        "HETATM   28  O3  BGA A   2     221.020 136.160  11.110  0.00  0.00           O\n"
        "HETATM   29  C4  BGA A   2     222.510 134.290  10.390  0.00  0.00           C\n"
        "HETATM   30  O4  BGA A   2     222.050 133.420  11.410  0.00  0.00           O\n"
        "HETATM   31  C6  BGA A   2     223.980 132.600   9.440  0.00  0.00           C\n"
        "HETATM   32  O6  BGA A   2     225.130 133.260   9.940  0.00  0.00           O\n"
        "HETATM   33  C1  ANE A   3     222.110 138.250  10.680  0.00  0.00           C\n"
        "HETATM   34  O11 ANE A   3     222.280 138.780   9.530  0.00  0.00           O\n"
        "HETATM   35  O12 ANE A   3     223.000 138.210  11.610  0.00  0.00           O\n"
        "HETATM   36  C2  ANE A   3     220.760 137.540  10.890  0.00  0.00           C\n"
        "HETATM   37  C6  ANE A   3     219.270 139.080   9.760  0.00  0.00           C\n"
        "HETATM   38  O6  ANE A   3     219.800 137.770   9.830  0.00  0.00           O\n"
        "HETATM   39  C3  ANE A   3     219.980 137.800  12.220  0.00  0.00           C\n"
        "HETATM   40  C4  ANE A   3     219.190 139.110  12.230  0.00  0.00           C\n"
        "HETATM   41  O4  ANE A   3     218.530 139.270  13.500  0.00  0.00           O\n"
        "HETATM   42  C5  ANE A   3     218.340 139.370  10.980  0.00  0.00           C\n"
        "HETATM   43  N   ANE A   3     217.770 140.720  10.910  0.00  0.00           N\n"
        "HETATM   44  C   ANE A   3     216.590 141.100  11.400  0.00  0.00           C\n"
        "HETATM   45  O   ANE A   3     215.770 140.350  11.840  0.00  0.00           O\n"
        "HETATM   46  CT  ANE A   3     216.270 142.580  11.320  0.00  0.00           C\n"
        "HETATM   47  C7  ANE A   3     218.670 139.400   8.410  0.00  0.00           C\n"
        "HETATM   48  O7  ANE A   3     217.770 138.390   8.060  0.00  0.00           O\n"
        "HETATM   49  C8  ANE A   3     219.730 139.520   7.270  0.00  0.00           C\n"
        "HETATM   50  O8  ANE A   3     220.710 138.450   7.350  0.00  0.00           O\n"
        "HETATM   51  C9  ANE A   3     220.360 140.880   7.030  0.00  0.00           C\n"
        "HETATM   52  O9  ANE A   3     221.100 141.240   8.210  0.00  0.00           O\n"
        "ATOM     53  N   THR B 323     180.790 133.490 398.290  0.00  0.00           N\n"
        "ATOM     54  CA  THR B 323     181.130 133.610   1.180  0.00  0.00           C\n"
        "ATOM     55  CB  THR B 323     181.340 132.260   1.860  0.00  0.00           C\n"
        "ATOM     56  OG1 THR B 323     181.250 132.440   3.290  0.00  0.00           O\n"
        "ATOM     57  CG2 THR B 323     182.640 131.430   1.500  0.00  0.00           C\n"
        "ATOM     58  C   THR B 323     182.370 134.490   1.350  0.00  0.00           C\n"
        "ATOM     59  O   THR B 323     182.780 134.850   2.480  0.00  0.00           O\n"
        "HETATM   60  C1  AGA B   1     179.860 132.350   3.670  0.00  0.00           C\n"
        "HETATM   61  C5  AGA B   1     180.490 130.470   4.970  0.00  0.00           C\n"
        "HETATM   62  O5  AGA B   1     179.510 131.000   4.050  0.00  0.00           O\n"
        "HETATM   63  C2  AGA B   1     179.580 133.160   4.900  0.00  0.00           C\n"
        "HETATM   64  N   AGA B   1     179.810 134.540   4.650  0.00  0.00           N\n"
        "HETATM   65  C   AGA B   1     178.940 135.570   4.510  0.00  0.00           C\n"
        "HETATM   66  O   AGA B   1     177.740 135.420   4.650  0.00  0.00           O\n"
        "HETATM   67  CT  AGA B   1     179.560 136.840   4.100  0.00  0.00           C\n"
        "HETATM   68  C3  AGA B   1     180.440 132.730   6.060  0.00  0.00           C\n"
        "HETATM   69  O3  AGA B   1     180.000 133.320   7.220  0.00  0.00           O\n"
        "HETATM   70  C4  AGA B   1     180.340 131.260   6.290  0.00  0.00           C\n"
        "HETATM   71  O4  AGA B   1     179.110 130.940   6.950  0.00  0.00           O\n"
        "HETATM   72  C6  AGA B   1     180.330 128.940   5.010  0.00  0.00           C\n"
        "HETATM   73  O6  AGA B   1     181.420 128.290   5.660  0.00  0.00           O\n"
        "HETATM   74  C1  BGA B   2     180.980 133.950   8.020  0.00  0.00           C\n"
        "HETATM   75  C5  BGA B   2     182.980 133.600   9.080  0.00  0.00           C\n"
        "HETATM   76  O5  BGA B   2     181.830 132.920   8.590  0.00  0.00           O\n"
        "HETATM   77  C2  BGA B   2     180.370 134.620   9.240  0.00  0.00           C\n"
        "HETATM   78  O2  BGA B   2     179.390 135.570   8.770  0.00  0.00           O\n"
        "HETATM   79  C3  BGA B   2     181.470 135.380  10.000  0.00  0.00           C\n"
        "HETATM   80  O3  BGA B   2     181.020 136.160  11.110  0.00  0.00           O\n"
        "HETATM   81  C4  BGA B   2     182.510 134.290  10.390  0.00  0.00           C\n"
        "HETATM   82  O4  BGA B   2     182.050 133.420  11.410  0.00  0.00           O\n"
        "HETATM   83  C6  BGA B   2     183.980 132.600   9.440  0.00  0.00           C\n"
        "HETATM   84  O6  BGA B   2     185.130 133.260   9.940  0.00  0.00           O\n"
        "HETATM   85  C1  ANE B   3     182.110 138.250  10.680  0.00  0.00           C\n"
        "HETATM   86  O11 ANE B   3     182.280 138.780   9.530  0.00  0.00           O\n"
        "HETATM   87  O12 ANE B   3     183.000 138.210  11.610  0.00  0.00           O\n"
        "HETATM   88  C2  ANE B   3     180.760 137.540  10.890  0.00  0.00           C\n"
        "HETATM   89  C6  ANE B   3     179.270 139.080   9.760  0.00  0.00           C\n"
        "HETATM   90  O6  ANE B   3     179.800 137.770   9.830  0.00  0.00           O\n"
        "HETATM   91  C3  ANE B   3     179.980 137.800  12.220  0.00  0.00           C\n"
        "HETATM   92  C4  ANE B   3     179.190 139.110  12.230  0.00  0.00           C\n"
        "HETATM   93  O4  ANE B   3     178.530 139.270  13.500  0.00  0.00           O\n"
        "HETATM   94  C5  ANE B   3     178.340 139.370  10.980  0.00  0.00           C\n"
        "HETATM   95  N   ANE B   3     177.770 140.720  10.910  0.00  0.00           N\n"
        "HETATM   96  C   ANE B   3     176.590 141.100  11.400  0.00  0.00           C\n"
        "HETATM   97  O   ANE B   3     175.770 140.350  11.840  0.00  0.00           O\n"
        "HETATM   98  CT  ANE B   3     176.270 142.580  11.320  0.00  0.00           C\n"
        "HETATM   99  C7  ANE B   3     178.670 139.400   8.410  0.00  0.00           C\n"
        "HETATM  100  O7  ANE B   3     177.770 138.390   8.060  0.00  0.00           O\n"
        "HETATM  101  C8  ANE B   3     179.730 139.520   7.270  0.00  0.00           C\n"
        "HETATM  102  O8  ANE B   3     180.710 138.450   7.350  0.00  0.00           O\n"
        "HETATM  103  C9  ANE B   3     180.360 140.880   7.030  0.00  0.00           C\n"
        "HETATM  104  O9  ANE B   3     181.100 141.240   8.210  0.00  0.00           O\n"
        "END\n";
    Loaded l;
    l.alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    l.sys = { .alloc = l.alloc };
    md_system_state_t st = { .alloc = l.alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&l.sys, &st, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&l.sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ(2u, (unsigned)md_molecule_count(&l.sys.molecule));
    ASSERT_EQ(4u, (unsigned)md_system_part_count(&l.sys));

    const md_bitset_t all = select(l, "glycan()");
    ASSERT_EQ(90u, (unsigned)md_bitset_count(all));
    md_array(ScriptSuggestion) s = suggest(l, all);
    print(s);
    EXPECT_TRUE(all_pick(l, all, s));
    EXPECT_TRUE(contains(s, "sel1 = glycan();"));

    // The glycan of chain A: what a click on it selects at Part granularity
    const md_bitset_t a = select(l, "part(atom(8))");
    ASSERT_EQ(45u, (unsigned)md_bitset_count(a));
    md_array(ScriptSuggestion) t = suggest(l, a);
    print(t);
    EXPECT_TRUE(all_pick(l, a, t));
    EXPECT_TRUE(contains(t, "sel1 = glycan() in chain(\"A\");"));
    EXPECT_TRUE(contains(t, "sel1 = part(2);"));
    EXPECT_FALSE(any_contains(t, "molecule("));
    unload(&l);
}

UTEST(script_suggestions, chains) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));
    const md_bitset_t b = select(l, "chain(2)");
    md_array(ScriptSuggestion) s = suggest(l, b);
    print(s);
    EXPECT_TRUE(all_pick(l, b, s));
    ASSERT_GT(md_array_size(s), 0u);
    EXPECT_TRUE(str_eq_cstr(s[0].text, "sel1 = chain(\"B\");"));
    EXPECT_TRUE(contains(s, "sel1 = chain(2);"));

    const md_bitset_t ac = select(l, "chain(1) or chain(3)");
    md_array(ScriptSuggestion) t = suggest(l, ac);
    print(t);
    EXPECT_TRUE(all_pick(l, ac, t));
    EXPECT_TRUE(contains(t, "sel1 = chain({\"A\", \"C\"});"));
    unload(&l);
}

// Residues of a chain by the numbers of the file, which is what the tooltip shows of them; by their position in the
// system as well
UTEST(script_suggestions, residues) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));
    const md_bitset_t sel = select(l, "resid(10:20) in chain(2)");
    md_array(ScriptSuggestion) s = suggest(l, sel);
    print(s);
    EXPECT_TRUE(all_pick(l, sel, s));
    ASSERT_GT(md_array_size(s), 0u);
    EXPECT_TRUE(str_eq_cstr(s[0].text, "sel1 = resid(10:20) in chain(\"B\");"));
    EXPECT_TRUE(contains(s, "sel1 = residue(405:415);"));
    // Eleven residues are eleven selections, which behaves differently in a context from one: the hint says so
    EXPECT_TRUE(str_eq_cstr(s[0].hint, "82 atoms, as an array of 11 selections"));
    unload(&l);
}

// What the atoms are, rather than where: their names within the chain, the backbone of residues
UTEST(script_suggestions, what_the_atoms_are) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));

    const md_bitset_t ca = select(l, "name(\"CA\") and chain(1)");
    md_array(ScriptSuggestion) s = suggest(l, ca);
    print(s);
    EXPECT_TRUE(all_pick(l, ca, s));
    ASSERT_GT(md_array_size(s), 0u);
    EXPECT_TRUE(str_eq_cstr(s[0].text, "sel1 = name(\"CA\") and chain(\"A\");"));

    const md_bitset_t bb = select(l, "backbone() and residue(10:12)");
    md_array(ScriptSuggestion) t = suggest(l, bb);
    print(t);
    EXPECT_TRUE(all_pick(l, bb, t));
    ASSERT_GT(md_array_size(t), 0u);
    EXPECT_TRUE(str_eq_cstr(t[0].text, "sel1 = backbone() and (resid(10:12) in chain(\"A\"));"));
    // The atoms of each residue by name, in the order the residue has them
    EXPECT_TRUE(contains(t, "sel1 = name({\"N\", \"CA\", \"C\"}) in resid(10:12) in chain(\"A\");"));
    unload(&l);

    ASSERT_TRUE(load(&l, "npt.gro"));
    const md_bitset_t w = select(l, "water()");
    md_array(ScriptSuggestion) u = suggest(l, w);
    print(u);
    EXPECT_TRUE(all_pick(l, w, u));
    ASSERT_GT(md_array_size(u), 0u);
    EXPECT_TRUE(str_eq_cstr(u[0].text, "sel1 = water();"));
    // The atoms of whole residues by name say nothing the residues do not
    EXPECT_FALSE(any_contains(u, "name(\"OW\""));

    const md_bitset_t o = select(l, "element(\"O\") and water()");
    md_array(ScriptSuggestion) v = suggest(l, o);
    print(v);
    EXPECT_TRUE(all_pick(l, o, v));
    EXPECT_TRUE(contains(v, "sel1 = name(\"OW\");"));
    EXPECT_TRUE(contains(v, "sel1 = element(\"O\") and water();"));
    unload(&l);
}

// 'in' binds looser than 'or': the terms of several residues are put in parentheses, or the first context would take
// all of what follows it
UTEST(script_suggestions, scattered) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));
    const md_bitset_t sel = select(l, "atom({5, 100, 2000, 7000})");
    md_array(ScriptSuggestion) s = suggest(l, sel);
    print(s);
    EXPECT_TRUE(all_pick(l, sel, s));
    EXPECT_TRUE(contains(s, "sel1 = atom({5, 100, 2000, 7000});"));
    EXPECT_TRUE(contains(s, "sel1 = (atom(5) in residue({1, 12})) or (atom(2) in residue(261)) or (atom(8) in residue(916));"));

    // Nothing short describes the atoms around a residue: they are listed, and the hint says so
    const md_bitset_t near = select(l, "within(6, residue(30)) and not residue(30)");
    ASSERT_GT(md_bitset_count(near), 20u);
    md_array(ScriptSuggestion) t = suggest(l, near);
    print(t);
    EXPECT_TRUE(all_pick(l, near, t));
    EXPECT_GT(md_array_size(t), 0u);
    unload(&l);
}

// Selections the user can make in any system, of any kind: whatever is offered for them picks them
UTEST(script_suggestions, every_line_picks_the_selection) {
    const char* files[] = {"1k4r.pdb", "npt.gro", "1LAF.pdb", "8g7u.cif", "2or2.cif", "nucl-dna.gro", "dppc64.pdb", "c60.pdb", "tryptophan.xyz", "1fez.cif", "tubulin-A-B.pdb"};
    const char* exprs[] = {
        "all()", "atom(1)", "atom(1:7)", "atom({3, 9, 27})", "residue(1)", "residue(2:4)", "chain(1)", "molecule(2)",
        "molecule(2:4)", "atom(2:3) in residue(2)", "element(\"C\")", "element(\"O\")", "protein()", "water()", "nucleic()",
        "backbone()", "not element(\"H\") and residue(1:3)", "residue(3) or atom(1)", "lipid()", "ion()",
        "atom(1:5:400)", "name(\"CA\")", "within(5, atom(1))",
    };
    for (size_t f = 0; f < ARRAY_SIZE(files); ++f) {
        Loaded l;
        ASSERT_TRUE(load(&l, files[f]));
        for (size_t e = 0; e < ARRAY_SIZE(exprs); ++e) {
            const md_bitset_t sel = select(l, exprs[e]);
            if (md_bitset_empty(sel)) continue;
            md_array(ScriptSuggestion) s = suggest(l, sel);
            if (!all_pick(l, sel, s)) {
                printf("  %s: %s\n", files[f], exprs[e]);
                print(s);
                EXPECT_TRUE(false);
            }
        }
        unload(&l);
    }
}

// Selections of no particular shape: some atoms of some residues, some residues whole, in any molecule
UTEST(script_suggestions, every_line_picks_a_random_selection) {
    const char* files[] = {"1k4r.pdb", "npt.gro", "1LAF.pdb", "8g7u.cif", "nucl-dna.gro", "dppc64.pdb", "tubulin-A-B.pdb"};
    uint64_t rng = 0x9E3779B97F4A7C15ull;
    auto next = [&rng](uint32_t n) {
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        return (uint32_t)(rng % n);
    };
    size_t offered = 0;
    for (size_t f = 0; f < ARRAY_SIZE(files); ++f) {
        Loaded l;
        ASSERT_TRUE(load(&l, files[f]));
        const md_system_t& sys = l.sys;
        for (int iter = 0; iter < 40; ++iter) {
            md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys.atom.count, l.alloc);
            const int parts = 1 + (int)next(4);
            for (int p = 0; p < parts; ++p) {
                const uint32_t c = next((uint32_t)sys.residue.count);
                const md_urange_t r = md_residue_atom_range(&sys.residue, c);
                switch (next(4)) {
                case 0: md_bitmask_set_range(&mask, r.beg, r.end); break;                            // a residue
                case 1: md_bitmask_set_bit(&mask, r.beg + next(r.end - r.beg)); break;               // an atom of one
                case 2: {                                                                            // a few residues
                    const md_urange_t r2 = md_residue_atom_range(&sys.residue, MIN(c + 1 + next(5), (uint32_t)sys.residue.count - 1));
                    md_bitmask_set_range(&mask, r.beg, r2.end);
                } break;
                default: {                                                                           // the same atom of several residues
                    const uint32_t k = next(r.end - r.beg);
                    for (uint32_t i = 0; i < 4 && c + i < sys.residue.count; ++i) {
                        const md_urange_t ri = md_residue_atom_range(&sys.residue, c + i);
                        if (ri.beg + k < ri.end) md_bitmask_set_bit(&mask, ri.beg + k);
                    }
                } break;
                }
            }
            const md_bitset_t sel = md_bitset_copy(md_bitmask_view(&mask), l.alloc);
            md_array(ScriptSuggestion) s = suggest(l, sel);
            offered += md_array_size(s) > 0;
            if (!all_pick(l, sel, s)) {
                printf("  %s, selection %d\n", files[f], iter);
                print(s);
                EXPECT_TRUE(false);
            }
        }
        unload(&l);
    }
    // Small selections of residues and atoms are listed at the very least
    EXPECT_EQ(ARRAY_SIZE(files) * 40, offered);
}

// The names of an XYZ file are the symbols it writes (mdlib used to pad them with zeros, which made them unwritable)
UTEST(script_suggestions, names_of_an_xyz_file) {
    Loaded l;
    ASSERT_TRUE(load(&l, "tryptophan.xyz"));
    const md_bitset_t sel = select(l, "element(\"C\")");
    md_array(ScriptSuggestion) s = suggest(l, sel);
    print(s);
    EXPECT_TRUE(all_pick(l, sel, s));
    EXPECT_TRUE(contains(s, "sel1 = name(\"C\");") || contains(s, "sel1 = element(\"C\");"));
    for (size_t i = 0; i < md_array_size(s); ++i) EXPECT_TRUE(printable(s[i].text));
    unload(&l);
}

UTEST(script_suggestions, measurement_within_a_residue) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));
    md_array(ScriptSuggestion) s = measure(l, {"name(\"N\") in residue(5)", "name(\"CA\") in residue(5)", "name(\"C\") in residue(5)"});
    print(s);
    for (size_t i = 0; i < md_array_size(s); ++i) EXPECT_TRUE(compiles(l, s[i].text));
    EXPECT_TRUE(contains(s, "m1 = angle(31, 32, 33);"));
    EXPECT_TRUE(contains(s, "m1 = angle(name(\"N\"), name(\"CA\"), name(\"C\")) in residue(5);"));
    EXPECT_TRUE(contains(s, "m1 = angle(name(\"N\"), name(\"CA\"), name(\"C\")) in resid(5) in chain(\"A\");"));
    // The same angle in every residue of the name, which the hint counts
    ASSERT_TRUE(contains(s, "m1 = angle(name(\"N\"), name(\"CA\"), name(\"C\")) in resname(\"HIS\");"));
    for (size_t i = 0; i < md_array_size(s); ++i) {
        size_t loc;
        const bool every = str_find_str(&loc, s[i].text, STR_LIT("resname("));
        EXPECT_TRUE(every == str_eq_cstr(s[i].hint, "One value for each of the 51 HIS"));
    }
    unload(&l);

    // A residue of which there is only one: no count to give
    ASSERT_TRUE(load(&l, "npt.gro"));
    md_array(ScriptSuggestion) t = measure(l, {"atom(1) in resname(\"LIG\")", "atom(2) in resname(\"LIG\")"});
    print(t);
    EXPECT_TRUE(contains(t, "m1 = distance(name(\"H00\"), name(\"C01\")) in resname(\"LIG\");"));
    for (size_t i = 0; i < md_array_size(t); ++i) EXPECT_TRUE(str_empty(t[i].hint));
    unload(&l);
}

// Atoms of several residues, each by name within its own; one which some residue of the name does not have is not
// measured in every residue of the name
UTEST(script_suggestions, measurement_across_residues) {
    Loaded l;
    ASSERT_TRUE(load(&l, "1k4r.pdb"));
    md_array(ScriptSuggestion) s = measure(l, {"name(\"CA\") in residue(5)", "name(\"CA\") in residue(404)"});
    print(s);
    for (size_t i = 0; i < md_array_size(s); ++i) EXPECT_TRUE(compiles(l, s[i].text));
    EXPECT_TRUE(contains(s, "m1 = distance(name(\"CA\") in residue(5), name(\"CA\") in residue(404));"));
    EXPECT_TRUE(contains(s, "m1 = distance(name(\"CA\") in resid(5) in chain(\"A\"), name(\"CA\") in resid(9) in chain(\"B\"));"));
    EXPECT_FALSE(any_contains(s, "resname("));

    // Two, three or four atoms, and each once
    EXPECT_EQ(0u, (unsigned)md_array_size(measure(l, {"atom(1)"})));
    EXPECT_EQ(0u, (unsigned)md_array_size(measure(l, {"atom(1)", "atom(1)"})));
    EXPECT_EQ(0u, (unsigned)md_array_size(measure(l, {"atom(1)", "atom(2)", "atom(3)", "atom(4)", "atom(5)"})));
    unload(&l);
}
