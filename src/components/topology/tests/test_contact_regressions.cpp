#include "utest.h"
#include <topology/critical_point_distribution_core.h>
#include <topology/nci_surface_core.h>
#include <md_system.h>
#include <md_topo.h>
#include <md_util.h>
#include <md_xyz.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_unit.h>

#include <vector>

namespace {
void check_contact(int* utest_result, const char* xyz, int a, int b, const char* expected, bool anion = false) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), 1 << 20);
    md_system_t sys = {};
    md_system_state_t state = {};
    sys.alloc = state.alloc = alloc;
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &state, str_from_cstr(xyz), MD_XYZ_OPTION_NONE));
    for (size_t i = 0; i < sys.atom.count; ++i) sys.atom.flags[i] |= MD_ATOM_FLAG_QM;
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    using namespace critical_point_distribution;
    std::vector<Position> atoms(sys.atom.count);
    std::vector<md_topo_vert_t> vertices(sys.atom.count + 1);
    std::vector<md_topo_critical_point_type_t> types(vertices.size(), MD_TOPO_MAXIMUM);
    const double scale = md_unit_scale(md_unit_angstrom()) / md_unit_scale(md_unit_bohr_radius());
    for (size_t i = 0; i < atoms.size(); ++i) {
        atoms[i] = {state.xyz[i].x, state.xyz[i].y, state.xyz[i].z};
        vertices[i].x = (float)(atoms[i][0] * scale);
        vertices[i].y = (float)(atoms[i][1] * scale);
        vertices[i].z = (float)(atoms[i][2] * scale);
        int hydrogens = 0;
        for (auto it = md_bond_iter(&sys.bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it))
            hydrogens += md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&it)) == 1;
        EXPECT_EQ(md_atom_hydrogen_count(&sys.atom, i), hydrogens);
    }
    const uint32_t cp = (uint32_t)sys.atom.count;
    types[cp] = MD_TOPO_SPLIT_SADDLE;
    md_topo_edge_t edges[] = {{cp, (uint32_t)a}, {cp, (uint32_t)b}};
    md_topo_extremum_graph_t graph = {};
    graph.vertices = vertices.data();
    graph.types = types.data();
    graph.num_vertices = (uint32_t)vertices.size();
    graph.edges = edges;
    graph.num_edges = 2;
    const auto labels = chemical_labels(graph, sys, atoms);
    EXPECT_TRUE(labels[cp].kind == Kind::Contact);
    EXPECT_STREQ(expected, labels[cp].interaction.c_str());
    if (anion) EXPECT_EQ(-1, group_charge(sys, b));
    if (hydrogen_bond_class(expected))
        EXPECT_EQ(nci_surface::reference_category(nci_surface::family_of(labels[cp])), (size_t)0);
    else if (std::string_view(expected) == "C-H\xC2\xB7\xC2\xB7\xC2\xB7X")
        EXPECT_EQ(nci_surface::reference_category(nci_surface::family_of(labels[cp])), (size_t)3);
    md_arena_allocator_destroy(alloc);
}
}

UTEST(critical_point_chemistry, reviewed_pyridine_dimer_ch_n) {
    const char* xyz = R"(22
S66x8 58, equilibrium geometry
N -0.941211224 0.790041387 0.0117189102
C -0.922755241 -0.552378118 0.0353787504
H 0.0572405085 -1.01558805 0.0513549112
C -2.07651901 -1.33301818 0.0392903499
H -1.99652898 -2.41058564 0.0588771999
C -3.31631303 -0.703339577 0.01759905
H -4.23157501 -1.27908432 0.0197937693
C -3.34889531 0.687018812 -0.00708595989
H -4.28544426 1.2261045 -0.0246589892
C -2.14310384 1.38263357 -0.00889005046
H -2.13809967 2.46565247 -0.0277829692
N 2.53321123 -0.950029314 0.0425178893
C 3.73499012 -1.5432055 0.04459773
H 3.72976637 -2.62616801 0.0664869025
C 4.94092655 -0.848246992 0.0205963496
H 5.87736464 -1.38778222 0.0236903597
C 4.90860891 0.542057455 -0.00715035992
H 5.82398367 1.1173085 -0.0263318699
C 3.66892838 1.17234361 -0.00962745957
H 3.58915567 2.24990225 -0.0307160299
C 2.51501489 0.392333984 0.0155662
H 1.53510439 0.855996549 0.0139033599
)";
    check_contact(utest_result, xyz, 2, 11, "C-H\xC2\xB7\xC2\xB7\xC2\xB7X");
    check_contact(utest_result, xyz, 21, 0, "C-H\xC2\xB7\xC2\xB7\xC2\xB7X");
}

UTEST(critical_point_chemistry, reviewed_pyridine_ethyne_ch_n) {
    check_contact(utest_result, R"(15
S66x8 65, equilibrium geometry
N -0.0830324888 0.000714590016 1.05519998
C -0.202853754 -1.1417259 0.364933699
H -0.0984856263 -2.05509806 0.937432647
C -0.446781427 -1.19176364 -1.00451231
H -0.533649206 -2.14585519 -1.50417149
C -0.574682117 0.00343953003 -1.70430946
H -0.763683915 0.00448009977 -2.76872659
C -0.453456759 1.1972425 -1.00091648
H -0.545630813 2.1522727 -1.4977951
C -0.209311113 1.14450765 0.368367314
H -0.110167071 2.05669737 0.943573952
C 0.47183603 -0.00605819 5.54171896
H 0.58724606 -0.00548399985 6.59673262
C 0.339766264 -0.00660791993 4.33547163
H 0.221618146 -0.00634548999 3.27096629
)", 14, 0, "C-H\xC2\xB7\xC2\xB7\xC2\xB7X");
}

UTEST(critical_point_chemistry, reviewed_ammonia_cyanide_ionic_hbond) {
    check_contact(utest_result, R"(6
IHB100 077
N -0.828252971 1.80791438 -1.20548964
H -0.418648124 0.952789009 -0.796174467
H -0.371406853 2.57512808 -0.72367084
H -1.78969753 1.82020807 -0.881549656
C 0.347015053 -0.726530194 0.355711043
N 0.716545939 -1.57022452 1.07349384
)", 1, 4, "Ionic hydrogen bond", true);
}

UTEST(critical_point_chemistry, reviewed_ammonium_methylazide_ionic_hbond) {
    check_contact(utest_result, R"(12
IHB100 017
N 0.111620858 -3.43490195 0
H 0.0825956166 -2.38472104 0
H -0.358757675 -3.79181409 -0.831749499
H -0.358757675 -3.79181409 0.831749499
H 1.08091617 -3.7524972 0
H 1.07208192 2.679039 -0.893579185
H 1.07208192 2.679039 0.893579185
H -0.12307781 3.64089203 0
C 0.452161878 2.72387719 0
N -0.524340093 1.61599612 0
N -0.174571916 0.47307539 0
N 0.0222487841 -0.648567259 0
)", 1, 11, "Ionic hydrogen bond");
}

UTEST(critical_point_chemistry, reviewed_furan_lone_pair_pi_unchanged) {
    check_contact(utest_result, R"(20
SH250 236
O -1.48000002 -0.147 0.326999992
C -3.44700003 0.796000004 -0.128000006
C -3.62700009 -0.621999979 -0.0209999997
C -2.13199997 1.02999997 0.0900000036
C -2.41000009 -1.14300001 0.254000008
H -4.19899988 1.53499997 -0.338
H -4.54400015 -1.17200005 -0.135000005
H -1.52900004 1.91700006 0.115000002
H -2.04500008 -2.13899994 0.416000009
F 0.476999998 -2.16100001 1.30700004
C 1.51400006 -1.324 1.46700001
F 1.83700001 -0.597000003 -2.68600011
F 0.238999993 -1.71599996 -1.73199999
F 1.36300004 -0.691999972 2.64199996
F 2.62700009 -2.0539999 1.53100002
F 2.30299997 -2.28999996 -1.41299999
N 4.53999996 0.510999978 -0.0670000017
C 1.48399997 -1.24699998 -1.56599998
C 3.41400003 0.270000011 -0.050999999
As 1.505 0.0909999982 -0.0179999992
)", 9, 0, "Lone pair\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80");
}

UTEST(critical_point_chemistry, reviewed_hs368_ch_n) {
    check_contact(utest_result, R"(39
VeloxTest probes/HS-368
C 0.673681021 -0.898754835 -0.379980028
C -2.03479218 -0.762090862 0.220945343
C -1.40735829 -1.89118075 -0.297186553
N -0.0979113057 -1.95742059 -0.599003732
H 1.72925341 -0.986414611 -0.655086279
H -3.10764551 -0.760805011 0.425221711
H -1.99071467 -2.79761982 -0.493849933
C 0.145672441 0.27986595 0.163186073
C -1.2411077 0.370866954 0.4555749
C -1.48666763 1.70596254 0.924790382
C -0.289262861 2.36644959 0.897821546
N 0.721433997 1.52061141 0.448082119
H -2.4385438 2.12565851 1.24155557
H -0.0530639179 3.39144969 1.16438127
C 4.92831802 2.51760292 0.127622694
C 2.08008265 1.85084474 0.363411218
C 3.08451581 0.959719718 0.686347306
C 4.45816994 1.28661919 0.561841786
C 6.35188484 2.81526971 0.0172995273
S 7.54641342 1.62332368 0.453100711
C 8.85678673 2.70866346 0.0906304047
C 8.37601376 3.92167354 -0.369841993
C 6.96394396 3.98416591 -0.413659602
H 9.00606537 4.74617386 -0.69466418
H 6.39679241 4.8497324 -0.746352255
C 10.2104092 2.11766648 0.236814573
O 11.3012428 2.91294694 0.27669552
O 10.3586912 0.922419429 0.31782186
C 11.2801695 4.31002188 0.517380893
H 12.2646723 4.5651288 0.933308125
H 11.1412859 4.87832212 -0.416908324
H 10.5002975 4.59460831 1.24008667
C 3.90955591 3.48254704 -0.202000886
C 2.49273324 3.15046406 -0.0856309235
N 1.6846993 4.15085936 -0.438353866
S 2.63295722 5.41126966 -0.880720794
N 4.10321665 4.72928143 -0.639525652
H 5.17715693 0.512846351 0.83924979
H 2.82023478 -0.026623698 1.06583178
)", 24, 36, "C-H\xC2\xB7\xC2\xB7\xC2\xB7X");
}
