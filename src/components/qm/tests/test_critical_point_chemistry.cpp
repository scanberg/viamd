#include "utest.h"
#include <qm/critical_point_distribution_core.h>
#include <md_system.h>
#include <md_topo.h>
#include <md_util.h>
#include <core/md_unit.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <md_gto.h>

#include <array>
#include <cmath>
#include <span>
#include <string>
#include <vector>

namespace {
const double ANGSTROM_TO_BOHR = md_unit_scale(md_unit_angstrom()) / md_unit_scale(md_unit_bohr_radius());
struct Fixture {
    md_system_t sys = {};
    md_particle_type_idx_t type_idx[3] = {0, 1, 2};
    md_atomic_number_t z[3] = {6, 8, 1};
    md_atom_flags_t flags[3] = {};
    md_topo_vert_t vertices[4] = {};
    md_topo_critical_point_type_t types[4] = {MD_TOPO_MAXIMUM, MD_TOPO_MAXIMUM, MD_TOPO_SPLIT_SADDLE, MD_TOPO_JOIN_SADDLE};
    md_topo_edge_t edges[2] = {{2, 0}, {2, 1}};
    md_topo_extremum_graph_t graph = {};
    critical_point_distribution::Position atoms[3] = {{0.0, 0.0, 0.0}, {1.2, 0.0, 0.0}, {4.0, 0.0, 0.0}};

    Fixture() {
        sys.atom.count = 3;
        sys.atom.type_idx = type_idx;
        sys.atom.type.count = 3;
        sys.atom.type.z = z;
        sys.atom.flags = flags;
        vertices[1].x = (float)(1.2 * ANGSTROM_TO_BOHR);
        graph.vertices = vertices;
        graph.types = types;
        graph.num_vertices = 4;
        graph.edges = edges;
        graph.num_edges = 2;
    }

    auto labels() { return critical_point_distribution::chemical_labels(graph, sys, atoms); }
};
}

UTEST(critical_point_chemistry, hierarchy_uses_chemistry_not_character_or_detail) {
    Fixture f;
    using namespace critical_point_distribution;
    ChemicalLabel bond;
    bond.kind = Kind::Bond;
    bond.atom[0] = 0;
    bond.atom[1] = 1;
    bond.group = "C carbonyl - O";
    auto path = environment_path(bond, f.sys);
    ASSERT_EQ(path.size(), (size_t)2);
    EXPECT_STREQ(path[0].c_str(), "Bonds");
    EXPECT_STREQ(path[1].c_str(), "C-O");
    std::swap(bond.atom[0], bond.atom[1]);
    bond.group = "Different detailed name";
    bond.finer = {"An even finer name"};
    EXPECT_TRUE(environment_path(bond, f.sys) == path);
    bond.kind = Kind::BondPath;
    path = environment_path(bond, f.sys);
    EXPECT_STREQ(path[0].c_str(), "Bond paths (connectivity unknown)");
    EXPECT_STREQ(path[1].c_str(), "C~O");
    bond.kind = Kind::Contact;
    bond.interaction = "Other contact";
    path = environment_path(bond, f.sys);
    ASSERT_EQ(path.size(), (size_t)3);
    EXPECT_STREQ(path[0].c_str(), "Contacts");
    EXPECT_STREQ(path[1].c_str(), "Other contacts");
    EXPECT_STREQ(path[2].c_str(), "C\xC2\xB7\xC2\xB7\xC2\xB7O");
    bond.interaction = "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    path = environment_path(bond, f.sys);
    EXPECT_STREQ(path[1].c_str(), "\xCF\x80-stacking");
}

UTEST(critical_point_chemistry, hydrogen_bond_hierarchy_has_donor_and_acceptor) {
    Fixture f;
    using namespace critical_point_distribution;
    uint32_t offsets[] = {0, 0, 1, 2};
    md_atom_idx_t others[] = {2, 1};
    md_bond_idx_t bonds[] = {0, 0};
    f.sys.bond.conn.offset = offsets;
    f.sys.bond.conn.offset_count = 4;
    f.sys.bond.conn.atom_idx = others;
    f.sys.bond.conn.bond_idx = bonds;
    f.z[0] = 7;
    ChemicalLabel label;
    label.kind = Kind::Contact;
    label.interaction = "Hydrogen bond";
    label.atom[0] = 2;
    label.atom[1] = 0;
    const auto path = environment_path(label, f.sys);
    ASSERT_EQ(path.size(), (size_t)3);
    EXPECT_STREQ(path[0].c_str(), "Contacts");
    EXPECT_STREQ(path[1].c_str(), "Hydrogen bonds");
    EXPECT_STREQ(path[2].c_str(), "O-H\xC2\xB7\xC2\xB7\xC2\xB7N");
    std::swap(label.atom[0], label.atom[1]);
    EXPECT_TRUE(environment_path(label, f.sys) == path);
    label.atom[0] = -1;
    label.atom[1] = 100;
    EXPECT_EQ(environment_path(label, f.sys).size(), (size_t)2);
}

UTEST(critical_point_chemistry, hierarchy_ring_and_unassigned_branches) {
    Fixture f;
    using namespace critical_point_distribution;
    ChemicalLabel label;
    label.kind = Kind::Ring;
    label.group = "Aromatic 6-ring C6";
    auto path = environment_path(label, f.sys);
    EXPECT_STREQ(path[1].c_str(), "Aromatic rings");
    label.group = "Closed by Hydrogen bond: 8-ring";
    path = environment_path(label, f.sys);
    EXPECT_STREQ(path[1].c_str(), "Contact-closed rings");
    label.kind = Kind::Unassigned;
    EXPECT_TRUE(environment_path(label, f.sys) == std::vector<std::string>{"Unassigned"});
}

UTEST(critical_point_chemistry, hierarchy_tree_preserves_all_descendants_and_repeated_leaves) {
    using namespace critical_point_distribution;
    const std::vector<std::string> paths[] = {
        {"BCPs", "Bonds", "C-C", "C aromatic - C aromatic"},
        {"BCPs", "Contacts", "Hydrogen bonds", "N-H...O", "H-bond environment"},
        {"BCPs", "Bonds", "C-C", "C aromatic - C aromatic"},
        {"RCPs", "Rings", "Aromatic rings", "Aromatic 6-ring C6"},
        {"BCPs", "Contacts", "Hydrogen bonds", "O-H...O", "Another environment"},
    };
    const auto tree = environment_tree(paths);
    ASSERT_EQ(tree[0].children.size(), (size_t)2);
    EXPECT_EQ(tree[0].environments.size(), (size_t)5);
    const auto& bcps = tree[tree[0].children[0]];
    EXPECT_STREQ(bcps.name.c_str(), "BCPs");
    EXPECT_TRUE(bcps.environments == (std::vector<uint32_t>{0, 1, 2, 4}));
    ASSERT_EQ(bcps.children.size(), (size_t)2);
    const auto& bonds = tree[bcps.children[0]];
    EXPECT_TRUE(bonds.environments == (std::vector<uint32_t>{0, 2}));
    const auto& cc = tree[bonds.children[0]];
    ASSERT_EQ(cc.children.size(), (size_t)2);
    const auto& a = tree[cc.children[0]];
    const auto& b = tree[cc.children[1]];
    EXPECT_STREQ(a.name.c_str(), b.name.c_str());
    EXPECT_EQ(a.environment, 0);
    EXPECT_EQ(b.environment, 2);
    const auto& contacts = tree[bcps.children[1]];
    const auto& hbonds = tree[contacts.children[0]];
    EXPECT_TRUE(hbonds.environments == (std::vector<uint32_t>{1, 4}));
    EXPECT_EQ(hbonds.children.size(), (size_t)2);
    EXPECT_TRUE(environment_tree({})[0].children.empty());
}

UTEST(critical_point_chemistry, unique_atom_and_unknown_connectivity) {
    Fixture f;
    f.flags[0] = MD_ATOM_FLAG_AROMATIC;
    const auto labels = f.labels();
    EXPECT_STREQ(labels[0].group.c_str(), "C aromatic");
    EXPECT_STREQ(labels[1].group.c_str(), "O");
    EXPECT_EQ(labels[2].atom[0], 0);
    EXPECT_EQ(labels[2].atom[1], 1);
    EXPECT_STREQ(labels[2].group.c_str(), "C aromatic ~ O");
    EXPECT_STREQ(labels[3].group.c_str(), "Unassigned ring CP");
    EXPECT_TRUE(labels[0].kind == critical_point_distribution::Kind::Atom);
    EXPECT_TRUE(labels[1].kind == critical_point_distribution::Kind::Atom);
    EXPECT_TRUE(labels[2].kind == critical_point_distribution::Kind::BondPath);
    EXPECT_TRUE(labels[3].kind == critical_point_distribution::Kind::Unassigned);
    const std::string atom_name = critical_point_distribution::display_name(labels[0]);
    const std::string ring_name = critical_point_distribution::display_name(labels[3]);
    EXPECT_STREQ(atom_name.c_str(), "Atom: C aromatic");
    EXPECT_STREQ(ring_name.c_str(), "Unassigned ring CP");
}

UTEST(critical_point_chemistry, carbonyl_and_topology_bond) {
    Fixture f;
    uint32_t offsets[] = {0, 1, 2, 2};
    md_atom_idx_t others[] = {1, 0};
    md_bond_idx_t bond_idx[] = {0, 0};
    md_bond_flags_t flags[] = {md_bond_flags_set_order(MD_BOND_FLAG_NONE, MD_BOND_ORDER_DOUBLE)};
    f.sys.bond.count = 1;
    f.sys.bond.flags = flags;
    f.sys.bond.conn.offset = offsets;
    f.sys.bond.conn.offset_count = 4;
    f.sys.bond.conn.atom_idx = others;
    f.sys.bond.conn.bond_idx = bond_idx;
    const auto labels = f.labels();
    EXPECT_STREQ(labels[0].group.c_str(), "C carbonyl");
    EXPECT_STREQ(labels[1].group.c_str(), "O carbonyl");
    EXPECT_STREQ(labels[2].group.c_str(), "C carbonyl \xE2\x80\x93 O carbonyl");
    EXPECT_TRUE(labels[2].kind == critical_point_distribution::Kind::Bond);
}

UTEST(critical_point_chemistry, nonbonded_is_not_inferred_without_connectivity) {
    Fixture f;
    uint32_t offsets[] = {0, 0, 0, 0};
    f.sys.bond.conn.offset = offsets;
    f.sys.bond.conn.offset_count = 4;
    const auto labels = f.labels();
    EXPECT_STREQ(labels[2].group.c_str(), "Other contact: C \xC2\xB7\xC2\xB7\xC2\xB7 O");
    EXPECT_TRUE(labels[2].kind == critical_point_distribution::Kind::Contact);
}

UTEST(critical_point_chemistry, hydrogen_environments_use_existing_metadata) {
    Fixture f;
    uint8_t hydrogen_count[] = {0, 1, 0};
    f.sys.atom.hydrogen_count = hydrogen_count;
    auto labels = f.labels();
    EXPECT_STREQ(labels[1].group.c_str(), "O (H-bearing)");
    f.z[1] = 7;
    labels = f.labels();
    EXPECT_STREQ(labels[1].group.c_str(), "N (H-bearing)");
}

UTEST(critical_point_chemistry, unmatched_and_ambiguous_attractors_stay_unassigned) {
    Fixture f;
    f.atoms[0] = {0.4, 0.0, 0.0};
    auto labels = f.labels();
    EXPECT_EQ(labels[0].atom[0], -1);
    EXPECT_TRUE(labels[0].kind == critical_point_distribution::Kind::Unassigned);
    EXPECT_STREQ(labels[2].group.c_str(), "Unassigned BCP");
    EXPECT_TRUE(labels[2].kind == critical_point_distribution::Kind::Unassigned);
    f.atoms[0] = {0.0, 0.0, 0.0};
    f.atoms[2] = {0.1, 0.0, 0.0};
    labels = f.labels();
    EXPECT_EQ(labels[0].atom[0], -1);
}

UTEST(critical_point_chemistry, duplicate_attractors_and_missing_paths_are_not_guessed) {
    Fixture f;
    f.vertices[1].x = 0.1f;
    auto labels = f.labels();
    EXPECT_EQ(labels[0].atom[0], -1);
    EXPECT_EQ(labels[1].atom[0], -1);
    f.vertices[1].x = (float)(1.2 * ANGSTROM_TO_BOHR);
    f.edges[1].to = 0;
    labels = f.labels();
    EXPECT_STREQ(labels[2].group.c_str(), "Unassigned BCP");
}

UTEST(critical_point_chemistry, double_precision_and_reversed_edges_are_supported) {
    Fixture f;
    md_topo_vert_f64_t vertices[4] = {};
    vertices[1].x = 1.2 * ANGSTROM_TO_BOHR;
    f.graph.vertices_f64 = vertices;
    f.vertices[0].x = 100.0f;
    f.edges[0] = {0, 2};
    const auto labels = f.labels();
    EXPECT_EQ(labels[0].atom[0], 0);
    EXPECT_EQ(labels[2].atom[1], 1);
}

UTEST(critical_point_chemistry, ring_cp_at_a_unique_ring_centroid) {
    Fixture f;
    md_allocator_i* alloc = md_get_heap_allocator();
    const int32_t ring[3] = {0, 1, 2};
    md_ring_data_push(&f.sys.ring, ring, 3, MD_RING_FLAG_AROMATIC, alloc);
    f.vertices[3].x = (float)((0.0 + 1.2 + 4.0) / 3.0 * ANGSTROM_TO_BOHR);
    auto labels = f.labels();
    EXPECT_TRUE(labels[3].kind == critical_point_distribution::Kind::Ring);
    EXPECT_STREQ(labels[3].group.c_str(), "Aromatic 3-ring CHO");

    f.vertices[3].x += (float)(0.6 * ANGSTROM_TO_BOHR);
    labels = f.labels();
    EXPECT_TRUE(labels[3].kind == critical_point_distribution::Kind::Unassigned);

    f.vertices[3].x -= (float)(0.6 * ANGSTROM_TO_BOHR);
    md_ring_data_push(&f.sys.ring, ring, 3, MD_RING_FLAG_NONE, alloc);
    labels = f.labels();
    EXPECT_TRUE(labels[3].kind == critical_point_distribution::Kind::Unassigned);
    md_ring_data_free(&f.sys.ring, alloc);
}

UTEST(critical_point_chemistry, cage_cps_stay_unassigned) {
    Fixture f;
    f.types[3] = MD_TOPO_MINIMUM;
    const auto labels = f.labels();
    EXPECT_TRUE(labels[3].kind == critical_point_distribution::Kind::Unassigned);
    EXPECT_STREQ(labels[3].group.c_str(), "Cage CP (unassigned)");
}

namespace {
using critical_point_distribution::Kind;
using critical_point_distribution::Position;

// Atoms with an attractor at each nucleus (vertex index = atom index), bonds, and further critical points.
struct Scene {
    md_system_t sys = {};
    std::vector<md_particle_type_idx_t> type_idx;
    std::vector<md_atomic_number_t> z;
    std::vector<md_atom_flags_t> flags;
    std::vector<Position> atoms;
    std::vector<std::pair<int, int>> bonds;
    std::vector<md_topo_vert_t> vertices;
    std::vector<md_topo_critical_point_type_t> types;
    std::vector<md_topo_edge_t> edges;
    std::vector<uint32_t> offsets;
    std::vector<md_atom_idx_t> others;
    std::vector<md_bond_idx_t> bond_idx;
    md_topo_extremum_graph_t graph = {};

    uint32_t point(md_topo_critical_point_type_t type, Position p, float rho) {
        md_topo_vert_t v = {};
        v.x = (float)(p[0] * ANGSTROM_TO_BOHR); v.y = (float)(p[1] * ANGSTROM_TO_BOHR); v.z = (float)(p[2] * ANGSTROM_TO_BOHR);
        v.value = rho;
        vertices.push_back(v);
        types.push_back(type);
        return (uint32_t)vertices.size() - 1;
    }
    int atom(md_atomic_number_t number, Position p, md_atom_flags_t f = MD_ATOM_FLAG_NONE) {
        type_idx.push_back((md_particle_type_idx_t)z.size());
        z.push_back(number);
        flags.push_back(f);
        atoms.push_back(p);
        point(MD_TOPO_MAXIMUM, p, 100.0f);
        return (int)atoms.size() - 1;
    }
    // A BCP at the midpoint with paths to both attractors
    uint32_t path(int a, int b, float rho) {
        Position mid;
        for (int k = 0; k < 3; ++k) mid[k] = 0.5 * (atoms[a][k] + atoms[b][k]);
        const uint32_t i = point(MD_TOPO_SPLIT_SADDLE, mid, rho);
        edges.push_back({i, (uint32_t)a});
        edges.push_back({i, (uint32_t)b});
        return i;
    }
    uint32_t bond(int a, int b) { bonds.push_back({a, b}); return path(a, b, 0.3f); }

    std::vector<critical_point_distribution::ChemicalLabel> labels(std::span<const std::vector<uint32_t>> ring_surfaces = {}) {
        const size_t n = atoms.size();
        offsets.assign(n + 1, 0);
        others.assign(2 * bonds.size(), 0);
        bond_idx.assign(2 * bonds.size(), 0);
        for (auto [a, b] : bonds) { offsets[a + 1]++; offsets[b + 1]++; }
        for (size_t i = 0; i < n; ++i) offsets[i + 1] += offsets[i];
        std::vector<uint32_t> fill(offsets.begin(), offsets.end() - 1);
        for (size_t k = 0; k < bonds.size(); ++k) {
            auto [a, b] = bonds[k];
            others[fill[a]] = b; bond_idx[fill[a]++] = (md_bond_idx_t)k;
            others[fill[b]] = a; bond_idx[fill[b]++] = (md_bond_idx_t)k;
        }
        sys.atom.count = n;
        sys.atom.type_idx = type_idx.data();
        sys.atom.type.count = n;
        sys.atom.type.z = z.data();
        sys.atom.flags = flags.data();
        sys.bond.count = bonds.size();
        sys.bond.conn.offset = offsets.data();
        sys.bond.conn.offset_count = n + 1;
        sys.bond.conn.atom_idx = others.data();
        sys.bond.conn.bond_idx = bond_idx.data();
        graph.vertices = vertices.data();
        graph.types = types.data();
        graph.num_vertices = (uint32_t)vertices.size();
        graph.edges = edges.data();
        graph.num_edges = (uint32_t)edges.size();
        return critical_point_distribution::chemical_labels(graph, sys, atoms, ring_surfaces);
    }
};
}

UTEST(critical_point_chemistry, hydrogen_bond_needs_a_polar_donor_and_its_angle) {
    Scene s;
    const int n = s.atom(7, {0.0, 0.0, 0.0});
    const int o = s.atom(8, {3.0, 0.0, 0.0});
    const int h = s.atom(1, {2.0, 0.0, 0.0});
    s.bond(o, h);
    const uint32_t contact = s.path(n, h, 0.02f);
    // Molecules {N} and {O, H}
    uint32_t mol_offset[] = {0, 1, 3};
    int32_t slot_atom[] = {0, 1, 2};
    s.sys.molecule.count = 2;
    s.sys.molecule.offset = mol_offset;
    s.sys.molecule.atom_idx = slot_atom;
    s.sys.molecule.atom_slot = slot_atom;
    auto labels = s.labels();
    EXPECT_TRUE(labels[contact].kind == Kind::Contact);
    std::string group = labels[contact].group;
    EXPECT_STREQ(group.c_str(), "Hydrogen bond (intermolecular): H (on O) \xC2\xB7\xC2\xB7\xC2\xB7 N");
    group = labels[h].group;
    EXPECT_STREQ(group.c_str(), "H (on O)");

    // O-H...N at 72 degrees is not a hydrogen bond
    s.atoms[h] = {3.0, 1.0, 0.0};
    s.vertices[h].x = (float)(3.0 * ANGSTROM_TO_BOHR);
    s.vertices[h].y = (float)(1.0 * ANGSTROM_TO_BOHR);
    labels = s.labels();
    group = labels[contact].group;
    EXPECT_STREQ(group.c_str(), "Other contact (intermolecular): H (on O) \xC2\xB7\xC2\xB7\xC2\xB7 N");
}

UTEST(critical_point_chemistry, pi_stacking_between_aromatic_atoms) {
    Scene s;
    const int a = s.atom(6, {0.0, 0.0, 0.0}, MD_ATOM_FLAG_AROMATIC);
    const int b = s.atom(6, {0.0, 0.0, 3.4}, MD_ATOM_FLAG_AROMATIC);
    const uint32_t contact = s.path(a, b, 0.01f);
    const auto labels = s.labels();
    const std::string group = labels[contact].group;
    EXPECT_STREQ(group.c_str(), "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80: C aromatic \xC2\xB7\xC2\xB7\xC2\xB7 C aromatic");
}

UTEST(critical_point_chemistry, ring_cp_inside_a_contact_closed_cycle) {
    Scene s;
    const int a = s.atom(6, {0.0, 0.0, 0.0});
    const int b = s.atom(6, {1.5, 0.0, 0.0});
    const int c = s.atom(6, {1.5, 1.5, 0.0});
    const int d = s.atom(6, {0.0, 1.5, 0.0});
    s.bond(a, b); s.bond(b, c); s.bond(c, d);
    s.path(d, a, 0.05f);
    const uint32_t rcp = s.point(MD_TOPO_JOIN_SADDLE, {0.75, 0.75, 0.0}, 0.03f);
    auto labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Ring);
    std::string group = labels[rcp].group;
    EXPECT_STREQ(group.c_str(), "Closed by Other contact");

    // A ring CP with the density of the ring's contact is flagged
    s.vertices[rcp].value = 0.0495f;
    labels = s.labels();
    group = labels[rcp].group;
    EXPECT_STREQ(group.c_str(), "Closed by Other contact (near-degenerate)");

    // Too far out of the ring's plane
    s.vertices[rcp].value = 0.03f;
    s.vertices[rcp].z = (float)(1.2 * ANGSTROM_TO_BOHR);
    labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Unassigned);
}

UTEST(critical_point_chemistry, ring_cp_of_a_twelve_ring) {
    Scene s;
    const double r = 1.5 / (2.0 * std::sin(3.14159265358979323846 / 12.0));
    for (int k = 0; k < 12; ++k) s.atom(6, {r * std::cos(k * 3.14159265358979323846 / 6.0), r * std::sin(k * 3.14159265358979323846 / 6.0), 0.0});
    for (int k = 0; k < 12; ++k) s.bond(k, (k + 1) % 12);
    const uint32_t rcp = s.point(MD_TOPO_JOIN_SADDLE, {0.0, 0.0, 0.0}, 0.01f);
    const auto labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Ring);
    const std::string group = labels[rcp].group;
    EXPECT_STREQ(group.c_str(), "Bonded 12-ring C12");
}

UTEST(critical_point_chemistry, ring_cp_denser_than_a_bcp_of_its_cycle) {
    Scene s;
    const int a = s.atom(6, {0.0, 0.0, 0.0});
    const int b = s.atom(6, {1.5, 0.0, 0.0});
    const int c = s.atom(6, {1.5, 1.5, 0.0});
    const int d = s.atom(6, {0.0, 1.5, 0.0});
    s.bond(a, b); s.bond(b, c); s.bond(c, d);
    s.path(d, a, 0.005f);
    const uint32_t rcp = s.point(MD_TOPO_JOIN_SADDLE, {0.75, 0.75, 0.0}, 0.01f);
    const auto labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Unassigned);
    EXPECT_NE(labels[rcp].detail.find("less dense than it"), std::string::npos);
}

UTEST(critical_point_chemistry, ring_cp_with_parallel_bond_paths) {
    // d-a joined by a weak path and a stronger parallel one: the stronger one satisfies the density rule
    Scene s;
    const int a = s.atom(6, {0.0, 0.0, 0.0});
    const int b = s.atom(6, {1.5, 0.0, 0.0});
    const int c = s.atom(6, {1.5, 1.5, 0.0});
    const int d = s.atom(6, {0.0, 1.5, 0.0});
    s.bond(a, b); s.bond(b, c); s.bond(c, d);
    s.path(d, a, 0.005f);
    s.path(d, a, 0.05f);
    const uint32_t rcp = s.point(MD_TOPO_JOIN_SADDLE, {0.75, 0.75, 0.0}, 0.01f);
    const auto labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Ring);
}

UTEST(critical_point_chemistry, ring_cp_between_cycles_of_the_same_size) {
    // Squares a-b-c-d and a-b-e-f share a-b; the larger one holds c and d in its face, so the smaller is the ring
    Scene s;
    const int a = s.atom(6, {-1.0, -1.0, 0.0});
    const int b = s.atom(6, {1.0, -1.0, 0.0});
    const int c = s.atom(6, {1.0, 1.0, 0.0});
    const int d = s.atom(6, {-1.0, 1.0, 0.0});
    const int e = s.atom(6, {1.5, 1.5, 0.0});
    const int f = s.atom(6, {-1.5, 1.5, 0.0});
    s.bond(a, b); s.bond(b, c); s.bond(c, d); s.bond(d, a);
    s.bond(b, e); s.bond(e, f); s.bond(f, a);
    const uint32_t rcp = s.point(MD_TOPO_JOIN_SADDLE, {0.0, 0.0, 0.0}, 0.01f);
    auto labels = s.labels();
    EXPECT_TRUE(labels[rcp].kind == Kind::Ring);
    EXPECT_NE(labels[rcp].detail.find("atoms 1, 2, 3, 4 "), std::string::npos);

    // Two mirror-image squares folded about a-b: neither face holds the other's atoms and their areas match
    Scene m;
    m.atom(6, {-1.0, -1.0, 0.0});
    m.atom(6, {1.0, -1.0, 0.0});
    m.atom(6, {1.0, 1.0, 1.0});
    m.atom(6, {-1.0, 1.0, 1.0});
    m.atom(6, {1.0, 1.0, -1.0});
    m.atom(6, {-1.0, 1.0, -1.0});
    m.bond(a, b); m.bond(b, c); m.bond(c, d); m.bond(d, a);
    m.bond(b, e); m.bond(e, f); m.bond(f, a);
    const uint32_t folded = m.point(MD_TOPO_JOIN_SADDLE, {0.0, 0.0, 0.0}, 0.01f);
    labels = m.labels();
    EXPECT_TRUE(labels[folded].kind == Kind::Unassigned);
    EXPECT_NE(labels[folded].detail.find("Equally small"), std::string::npos);
    EXPECT_TRUE(critical_point_distribution::ring_surface_may_decide(labels[folded]));

    // Its ring surface reaching c (vertex 2) decides for a-b-c-d
    std::vector<std::vector<uint32_t>> surfaces(m.vertices.size());
    surfaces[folded] = {0, 1, 2};
    labels = m.labels(surfaces);
    EXPECT_TRUE(labels[folded].kind == Kind::Ring);
    EXPECT_NE(labels[folded].detail.find("atoms 1, 2, 3, 4 "), std::string::npos);
    EXPECT_NE(labels[folded].detail.find("ring surface"), std::string::npos);

    // Reaching only the shared a and b, or nuclei of both squares, decides nothing
    surfaces[folded] = {0, 1};
    labels = m.labels(surfaces);
    EXPECT_TRUE(labels[folded].kind == Kind::Unassigned);
    EXPECT_NE(labels[folded].detail.find("ring surface"), std::string::npos);
    surfaces[folded] = {2, 0, 4};
    labels = m.labels(surfaces);
    EXPECT_TRUE(labels[folded].kind == Kind::Unassigned);
}

UTEST(critical_point_chemistry, cage_cp_from_its_assigned_faces) {
    Scene s;
    const Position p[4] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    for (const Position& q : p) s.atom(6, q);
    for (int i = 0; i < 4; ++i) for (int j = i + 1; j < 4; ++j) s.bond(i, j);
    const uint32_t ccp = s.point(MD_TOPO_MINIMUM, {0.0, 0.0, 0.0}, 0.01f);
    std::vector<uint32_t> faces;
    for (int skip = 0; skip < 4; ++skip) {
        Position c = {0, 0, 0};
        for (int k = 0; k < 4; ++k) if (k != skip) for (int x = 0; x < 3; ++x) c[x] += p[k][x] / 3.0;
        faces.push_back(s.point(MD_TOPO_JOIN_SADDLE, c, 0.05f));
        s.edges.push_back({faces.back(), ccp});
    }
    auto labels = s.labels();
    for (uint32_t f : faces) {
        EXPECT_TRUE(labels[f].kind == Kind::Ring);
        const std::string group = labels[f].group;
        EXPECT_STREQ(group.c_str(), "Bonded 3-ring C3");
    }
    EXPECT_TRUE(labels[ccp].kind == Kind::Cage);
    std::string group = labels[ccp].group;
    EXPECT_STREQ(group.c_str(), "Bonded cage, 4 faces");

    // With a face missing, it stays unassigned
    s.edges.pop_back();
    s.edges.pop_back();
    labels = s.labels();
    EXPECT_TRUE(labels[ccp].kind == Kind::Unassigned);
}

UTEST(critical_point_chemistry, non_nuclear_attractor_away_from_every_nucleus) {
    Fixture f;
    f.vertices[1].x = (float)(0.6 * ANGSTROM_TO_BOHR);
    const auto labels = f.labels();
    EXPECT_TRUE(labels[1].kind == Kind::NonNuclear);
    const std::string group = labels[1].group;
    EXPECT_STREQ(group.c_str(), "near C, O");
}

UTEST(critical_point_chemistry, ellipticity_of_negative_eigenvalues) {
    const double ev[3] = {-2.0, -1.0, 3.0};
    EXPECT_NEAR(critical_point_distribution::ellipticity(ev), 1.0, 1e-12);
}

UTEST(critical_point_chemistry, finer_detail_levels) {
    // Methyl carbon a bonded to aromatic-free carbon b with three H; b bonded to a and one H
    Scene s;
    const int a = s.atom(6, {0.0, 0.0, 0.0});
    const int b = s.atom(6, {1.5, 0.0, 0.0});
    const int h1 = s.atom(1, {-0.6, 0.9, 0.0});
    const int h2 = s.atom(1, {-0.6, -0.9, 0.0});
    const int h3 = s.atom(1, {-0.6, 0.0, 0.9});
    const int h4 = s.atom(1, {2.1, 0.9, 0.0});
    s.bond(a, b);
    s.bond(a, h1); s.bond(a, h2); s.bond(a, h3);
    const uint32_t ch = s.bond(b, h4);
    auto labels = s.labels();
    { const std::string name1 = critical_point_distribution::display_name(labels[a], 1); EXPECT_STREQ(name1.c_str(), "Atom: C [CH3]"); }
    { const std::string name2 = critical_point_distribution::display_name(labels[b], 1); EXPECT_STREQ(name2.c_str(), "Atom: C [CH]"); }
    { const std::string name3 = critical_point_distribution::display_name(labels[b], 0); EXPECT_STREQ(name3.c_str(), "Atom: C"); }
    { const std::string name4 = critical_point_distribution::display_name(labels[h1], 1); EXPECT_STREQ(name4.c_str(), "Atom: H (on C [CH3])"); }
    { const std::string name5 = critical_point_distribution::display_name(labels[ch], 1); EXPECT_STREQ(name5.c_str(), "Bond: C [CH] \xE2\x80\x93 H (on C [CH])"); }
    // No residues: no finer level than the neighbours
    { const std::string name6 = critical_point_distribution::display_name(labels[a], 2); EXPECT_STREQ(name6.c_str(), "Atom: C [CH3]"); }

    // The finest level keeping the carbons in at most the groups asked for
    const uint32_t carbons[2] = {(uint32_t)a, (uint32_t)b};
    EXPECT_EQ(critical_point_distribution::auto_detail_level(labels, carbons, 1), 0);
    EXPECT_EQ(critical_point_distribution::auto_detail_level(labels, carbons, 2), 2);

    // Residues: a and its hydrogens in ALA, b and h4 in GLY
    md_label_t names[2] = {};
    snprintf(names[0].buf, sizeof(names[0].buf), "ALA");
    snprintf(names[1].buf, sizeof(names[1].buf), "GLY");
    // Atoms are 0 a, 1 b, 2-4 h1-h3, 5 h4: residue ranges need contiguous atoms, so only a, b are placed: ALA [0,1), GLY [1,2)
    uint32_t offsets[3] = {0, 1, 2};
    s.sys.residue.count = 2;
    s.sys.residue.name = names;
    s.sys.residue.atom_offset = offsets;
    labels = s.labels();
    { const std::string name7 = critical_point_distribution::display_name(labels[a], 2); EXPECT_STREQ(name7.c_str(), "Atom: C [CH3] in ALA"); }
    const uint32_t cc = 6 + 0;   // the first bond's BCP follows the six attractors
    { const std::string name8 = critical_point_distribution::display_name(labels[cc], 2); EXPECT_STREQ(name8.c_str(), "Bond: C [CH3] \xE2\x80\x93 C [CH] between ALA and GLY"); }
    // An atom in no residue gets no residue level
    { const std::string name9 = critical_point_distribution::display_name(labels[h1], 2); EXPECT_STREQ(name9.c_str(), "Atom: H (on C [CH3])"); }
}

namespace {
// A curve sampled from a sum of Gaussians: (centre, height, width) each
critical_point_distribution::Curve gaussian_curve(std::initializer_list<std::array<double, 3>> peaks) {
    critical_point_distribution::Curve curve;
    curve.valid = true;
    const int n = critical_point_distribution::Curve::sample_count;
    for (int j = 0; j < n; ++j) {
        const double x = 10.0 * j / (n - 1);
        double y = 0.0;
        for (const auto& p : peaks) y += p[1] * std::exp(-0.5 * (x - p[0]) * (x - p[0]) / (p[2] * p[2]));
        curve.density[j] = x;
        curve.probability[j] = y;
    }
    return curve;
}
}

UTEST(critical_point_peaks, two_peaks_one_valley) {
    const auto curve = gaussian_curve({{2.0, 1.0, 0.4}, {7.0, 0.5, 0.4}});
    const auto valleys = critical_point_distribution::valleys(curve);
    ASSERT_EQ(valleys.size(), (size_t)1);
    EXPECT_GT(valleys[0].density, 2.5);
    EXPECT_LT(valleys[0].density, 6.5);
    EXPECT_NEAR(valleys[0].persistence, 0.5, 0.01);
}

UTEST(critical_point_peaks, shallow_ripples_are_ignored) {
    auto curve = gaussian_curve({{5.0, 1.0, 1.0}});
    for (int j = 0; j < critical_point_distribution::Curve::sample_count; ++j) curve.probability[j] += 0.002 * ((j % 7) - 3);
    EXPECT_EQ(critical_point_distribution::valleys(curve).size(), (size_t)0);
}

UTEST(critical_point_peaks, zoom_splits_within_the_overview) {
    // Three major peaks; the middle one has two shallow parts
    const auto curve = gaussian_curve({{1.5, 1.0, 0.3}, {4.6, 0.8, 0.3}, {5.4, 0.8, 0.3}, {8.5, 1.0, 0.3}});
    const auto valleys = critical_point_distribution::valleys(curve);
    ASSERT_EQ(valleys.size(), (size_t)3);
    // Overview with at most 3 peaks: the two deep valleys
    const auto overview = critical_point_distribution::peak_cuts(valleys, 1.0, 0.0, 3);
    ASSERT_EQ(overview.size(), (size_t)2);
    EXPECT_LT(overview[0], 4.0);
    EXPECT_GT(overview[1], 6.0);
    // Zoomed onto the middle peak: its shallow valley too, and the overview's cuts outside the view kept
    const auto zoomed = critical_point_distribution::peak_cuts(valleys, 4.0, 6.0, 3);
    ASSERT_EQ(zoomed.size(), (size_t)3);
    EXPECT_EQ(zoomed[0], overview[0]);
    EXPECT_NEAR(zoomed[1], 5.0, 0.1);
    EXPECT_EQ(zoomed[2], overview[1]);
    EXPECT_EQ(critical_point_distribution::peak_cuts(valleys, 1.0, 0.0, 1).size(), (size_t)0);
}

UTEST(critical_point_peaks, refined_cuts_nest_and_part_cps) {
    const std::vector<double> overview = {3.0};
    const std::vector<double> values = {1.0, 1.1, 4.0, 4.1, 5.0, 5.1, 6.0, 8.0};
    // View [3.5, 7], at most 3 peaks in it: the two most persistent valleys there; 2.0 lies outside the view
    const std::vector<critical_point_distribution::Valley> view = {{5.5, 0.5}, {4.5, 0.4}, {2.0, 0.9}, {4.55, 0.3}};
    const auto cuts = critical_point_distribution::refine_cuts(overview, view, values, 3.5, 7.0, 3);
    // The overview cut, outside the view, is kept
    ASSERT_EQ(cuts.size(), (size_t)3);
    EXPECT_EQ(cuts[0], 3.0);
    EXPECT_EQ(cuts[1], 4.5);
    EXPECT_EQ(cuts[2], 5.5);
    // A valley in the same gap as one already cut is skipped
    const auto same = critical_point_distribution::refine_cuts(overview, std::vector<critical_point_distribution::Valley>{{4.5, 0.5}, {4.6, 0.4}, {5.5, 0.3}}, values, 3.5, 7.0, 4);
    ASSERT_EQ(same.size(), (size_t)3);
    EXPECT_EQ(same[1], 4.5);
    EXPECT_EQ(same[2], 5.5);
}

UTEST(critical_point_energy, densities_and_hydrogen_bonds) {
    // A closed-shell contact: lap > 0, small G
    const auto e = critical_point_distribution::energy_densities(0.02, 0.1);
    EXPECT_NEAR(e.V, 0.025 - 0.04, 1e-15);
    EXPECT_NEAR(e.H, 0.02 - 0.015, 1e-15);
    EXPECT_NEAR(e.ratio, 0.75, 1e-12);
    EXPECT_TRUE(std::isnan(critical_point_distribution::energy_densities(NAN, 0.1).H));
    EXPECT_NEAR(critical_point_distribution::hydrogen_bond_energy(critical_point_distribution::HBondFormula::Emamian, 0.03, NAN), -5.9501, 1e-4);
    EXPECT_NEAR(critical_point_distribution::hydrogen_bond_energy(critical_point_distribution::HBondFormula::Espinosa, NAN, -0.02), -6.275, 1e-3);
    critical_point_distribution::ChemicalLabel hb;
    hb.kind = critical_point_distribution::Kind::Contact;
    hb.interaction = "Hydrogen bond";
    EXPECT_TRUE(critical_point_distribution::hydrogen_bond_energy_applies(hb));
    hb.interaction = "Other contact";
    EXPECT_FALSE(critical_point_distribution::hydrogen_bond_energy_applies(hb));
}

UTEST(critical_point_energy, mayer_bond_orders) {
    // H2 in an orthonormal minimal basis: sigma doubly occupied is a single bond; H2+ (one alpha electron) half of one
    const std::vector<double> S = {1, 0, 0, 1};
    const std::vector<double> P = {1, 1, 1, 1};
    const std::vector<double> Pa = {0.5, 0.5, 0.5, 0.5}, Pb = {0, 0, 0, 0};
    const std::vector<uint32_t> atom = {0, 1};
    const std::vector<std::array<int, 2>> pairs = {{0, 1}, {1, 0}, {0, 2}};
    const double* total[] = {P.data()};
    const auto b = critical_point_distribution::mayer_bond_orders(S, total, 1.0, 2, atom, pairs);
    EXPECT_NEAR(b[0], 1.0, 1e-12);
    EXPECT_NEAR(b[1], 1.0, 1e-12);
    EXPECT_TRUE(std::isnan(b[2]));
    const double* spin[] = {Pa.data(), Pb.data()};
    EXPECT_NEAR(critical_point_distribution::mayer_bond_orders(S, spin, 2.0, 2, atom, pairs)[0], 0.5, 1e-12);
    // Non-orthogonal: the same sigma orbital with overlap s, c = 1 / sqrt(2 (1 + s)) each; still order 1
    const double s = 0.6, c2 = 1.0 / (2.0 * (1.0 + s));
    const std::vector<double> S2 = {1, s, s, 1};
    const std::vector<double> P2 = {2 * c2, 2 * c2, 2 * c2, 2 * c2};
    const double* total2[] = {P2.data()};
    EXPECT_NEAR(critical_point_distribution::mayer_bond_orders(S2, total2, 1.0, 2, atom, pairs)[0], 1.0, 1e-12);
}

UTEST(critical_point_energy, approximate_interaction_energy) {
    constexpr double kcal = 627.5094740631;
    // Neutral atoms sharing a single bond 2 Bohr apart: -1/4 Hartree; like charges, no bond: repulsive
    EXPECT_NEAR(critical_point_distribution::approximate_interaction_energy(0.0, 0.0, 1.0, 2.0), -0.25 * kcal, 1e-9);
    EXPECT_NEAR(critical_point_distribution::approximate_interaction_energy(0.5, 0.5, 0.0, 1.0), 0.25 * kcal, 1e-9);
    EXPECT_LT(critical_point_distribution::approximate_interaction_energy(0.4, -0.4, 0.05, 4.0), 0.0);
    EXPECT_TRUE(std::isnan(critical_point_distribution::approximate_interaction_energy(NAN, 0.0, 1.0, 2.0)));
    EXPECT_TRUE(std::isnan(critical_point_distribution::approximate_interaction_energy(0.0, 0.0, 1.0, 0.0)));

    // Mulliken: the H2 sigma orbital with overlap s puts one electron on each atom
    const double s = 0.6, c2 = 1.0 / (2.0 * (1.0 + s));
    const std::vector<double> S = {1, s, s, 1};
    const std::vector<double> P = {2 * c2, 2 * c2, 2 * c2, 2 * c2};
    const std::vector<uint32_t> atom = {0, 1};
    const double* total[] = {P.data()};
    const auto N = critical_point_distribution::mulliken_populations(S, total, 2, atom, 3);
    EXPECT_NEAR(N[0], 1.0, 1e-12);
    EXPECT_NEAR(N[1], 1.0, 1e-12);
    EXPECT_TRUE(std::isnan(N[2]));
}

UTEST(critical_point_energy, populations_trace_and_cancellation) {
    using namespace critical_point_distribution;
    const std::vector<double> S = {1, 0.3, 0.3, 1};
    const std::vector<double> Pa = {0.7, 0.2, 0.2, 0.4}, Pb = {0.2, 0.1, 0.1, 0.8};
    const double* spin[] = {Pa.data(), Pb.data()};
    const std::vector<uint32_t> atoms = {0, 1};
    const std::vector<std::array<int, 2>> pairs = {{0, 1}, {1, 0}};
    const auto N = mulliken_populations(S, spin, 2, atoms, 2);
    double trace = 0.0;
    for (const double* P : spin) for (size_t i = 0; i < S.size(); ++i) trace += P[i] * S[i];
    EXPECT_NEAR(trace, N[0] + N[1], 1e-14);
    const auto cancelled = [] { return true; };
    const auto missing = mulliken_populations(S, spin, 2, atoms, 2, cancelled);
    for (double value : missing) EXPECT_TRUE(std::isnan(value));
    const auto stopped = mayer_bond_orders(S, spin, 2.0, 2, atoms, pairs, cancelled);
    for (double value : stopped) EXPECT_TRUE(std::isnan(value));
    const auto ordinary = mayer_bond_orders(S, spin, 2.0, 2, atoms, pairs);
    const auto resumed = mayer_bond_orders(S, spin, 2.0, 2, atoms, pairs, [] { return false; });
    for (size_t i = 0; i < pairs.size(); ++i) EXPECT_EQ(ordinary[i], resumed[i]);
    size_t calls = 0;
    const auto interrupted = mulliken_populations(S, spin, 2, atoms, 2, [&] { return ++calls > 1; });
    EXPECT_EQ(size_t(2), calls);
    for (double value : interrupted) EXPECT_TRUE(std::isnan(value));
}

namespace {
// A promolecular density of exponential atoms a exp(-b r), with its gradient and Hessian
struct ExpAtoms {
    std::vector<std::array<double, 3>> centre;
    double a = 0.5, b = 1.8;
    void eval(const double* x, double& rho, double g[3], double h[6]) const {
        rho = 0.0;
        for (int k = 0; k < 3; ++k) g[k] = 0.0;
        for (int k = 0; k < 6; ++k) h[k] = 0.0;
        for (const auto& c : centre) {
            const double d[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
            const double r = std::max(std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), 1e-12);
            const double f = a * std::exp(-b * r);
            rho += f;
            const double u[3] = {d[0] / r, d[1] / r, d[2] / r};
            for (int k = 0; k < 3; ++k) g[k] += -b * f * u[k];
            static const int idx[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
            for (int k = 0; k < 6; ++k) {
                const int i = idx[k][0], j = idx[k][1];
                h[k] += f * (b * b * u[i] * u[j] - b / r * ((i == j ? 1.0 : 0.0) - u[i] * u[j]));
            }
        }
    }
    critical_point_distribution::DensityEvaluator evaluator() const {
        return [this](const double* xyz, size_t n, double* rho, double* g, double* h) {
            for (size_t t = 0; t < n; ++t) eval(xyz + 3 * t, rho[t], g + 3 * t, h + 6 * t);
            return true;
        };
    }
    // The BCP near x0 (Newton on the gradient) as an NCI site
    critical_point_distribution::NciSite site(std::array<double, 3> x) const {
        for (int it = 0; it < 50; ++it) {
            double rho, g[3], h[6];
            eval(x.data(), rho, g, h);
            const double m[3][3] = {{h[0], h[1], h[2]}, {h[1], h[3], h[4]}, {h[2], h[4], h[5]}};
            const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            double dx[3];
            for (int c = 0; c < 3; ++c) {
                double mc[3][3];
                for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) mc[i][j] = j == c ? -g[i] : m[i][j];
                dx[c] = (mc[0][0] * (mc[1][1] * mc[2][2] - mc[1][2] * mc[2][1]) - mc[0][1] * (mc[1][0] * mc[2][2] - mc[1][2] * mc[2][0]) +
                         mc[0][2] * (mc[1][0] * mc[2][1] - mc[1][1] * mc[2][0])) / det;
            }
            for (int c = 0; c < 3; ++c) x[c] += dx[c];
        }
        critical_point_distribution::NciSite s = {};
        double g[3], h[6];
        eval(x.data(), s.rho, g, h);
        for (int c = 0; c < 3; ++c) s.center[c] = x[c];
        // Eigenvectors by inverse iteration on the known eigenvalues
        const auto l = critical_point_distribution::symmetric_eigenvalues(h);
        for (int k = 0; k < 3; ++k) {
            double v[3] = {0.3 + k, 0.7, 1.1 - 0.4 * k};
            for (int it = 0; it < 30; ++it) {
                const double shift = l[k] + 1e-7;
                const double m[3][3] = {{h[0] - shift, h[1], h[2]}, {h[1], h[3] - shift, h[4]}, {h[2], h[4], h[5] - shift}};
                const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                double w[3];
                for (int c = 0; c < 3; ++c) {
                    double mc[3][3];
                    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) mc[i][j] = j == c ? v[i] : m[i][j];
                    w[c] = (mc[0][0] * (mc[1][1] * mc[2][2] - mc[1][2] * mc[2][1]) - mc[0][1] * (mc[1][0] * mc[2][2] - mc[1][2] * mc[2][0]) +
                            mc[0][2] * (mc[1][0] * mc[2][1] - mc[1][1] * mc[2][0])) / det;
                }
                const double nw = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
                for (int c = 0; c < 3; ++c) v[c] = w[c] / nw;
            }
            for (int c = 0; c < 3; ++c) s.axes[k][c] = v[c];
        }
        return s;
    }
};
}

UTEST(critical_point_energy, symmetric_eigenvalues_and_reduced_gradient) {
    const double h[6] = {2.0, 1.0, 0.0, 2.0, 0.0, -1.0};   // eigenvalues -1, 1, 3
    const auto l = critical_point_distribution::symmetric_eigenvalues(h);
    EXPECT_NEAR(l[0], -1.0, 1e-12);
    EXPECT_NEAR(l[1], 1.0, 1e-12);
    EXPECT_NEAR(l[2], 3.0, 1e-12);
    const double d[6] = {0.5, 0.0, 0.0, 0.5, 0.0, 0.5};
    EXPECT_NEAR(critical_point_distribution::symmetric_eigenvalues(d)[1], 0.5, 1e-12);
    const double g[3] = {0.0, 0.0, 0.0}, g1[3] = {0.3, 0.0, 0.4};
    EXPECT_EQ(0.0, critical_point_distribution::reduced_density_gradient(0.01, g));
    EXPECT_NEAR(critical_point_distribution::reduced_density_gradient(1.0, g1), 0.5 / (2.0 * std::cbrt(3.0 * 9.869604401089358)), 1e-12);
}

UTEST(critical_point_energy, nci_attractive_and_crowded_contacts) {
    using critical_point_distribution::Character;
    // Two atoms 6 Bohr apart: an isolated contact, attractive
    ExpAtoms pair;
    pair.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    const auto site = pair.site({0.0, 0.0, 0.0});
    EXPECT_LT(site.rho, 0.05);
    // The test density is dilute: the weak band, by default, with an attractive lean; without it, the sign of lambda2 decides
    const auto weak = critical_point_distribution::nci_contacts(std::span(&site, 1), pair.evaluator());
    EXPECT_TRUE(weak[0].character == Character::WeakAttractive);
    EXPECT_LT(weak[0].share, critical_point_distribution::nci_weak_attractive_below);
    critical_point_distribution::NciParams signed_params;
    signed_params.rho_weak = 0.0;
    const auto isolated = critical_point_distribution::nci_contacts(std::span(&site, 1), pair.evaluator(), signed_params);
    ASSERT_EQ((size_t)1, isolated.size());
    EXPECT_TRUE(isolated[0].character == Character::Attractive);
    EXPECT_GT(isolated[0].samples, (size_t)10);
    EXPECT_LT(isolated[0].share, 0.05);

    // A third atom closing a crowded triangle: the ring CP's depletion reaches the contact's region
    ExpAtoms triangle;
    triangle.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {0.0, 4.5, 0.0}};   // ring CP about 1 Bohr from the BCP
    const auto crowded = triangle.site({0.0, 0.13, 0.0});
    EXPECT_NEAR(crowded.center[1], 0.13, 0.02);
    const auto c = critical_point_distribution::nci_contacts(std::span(&crowded, 1), triangle.evaluator(), signed_params);
    EXPECT_GT(c[0].share, isolated[0].share + 0.05);
    EXPECT_TRUE(c[0].character != Character::Attractive);
    EXPECT_NE(c[0].reason.find("repulsive"), std::string::npos);
    const auto c_weak = critical_point_distribution::nci_contacts(std::span(&crowded, 1), triangle.evaluator());
    EXPECT_GT(c_weak[0].share, weak[0].share);

    // Dense: no region sampled, covalent when bonded, strong attractive when not; a failing evaluator leaves the rest
    // undetermined
    ExpAtoms close;
    close.centre = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    critical_point_distribution::NciSite sites[3] = {close.site({0.0, 0.0, 0.0}), close.site({0.0, 0.0, 0.0}), site};
    sites[0].bonded = true;
    const auto all = critical_point_distribution::nci_contacts(sites, [](const double*, size_t, double*, double*, double*) { return false; });
    EXPECT_TRUE(all[0].character == Character::Covalent);
    EXPECT_TRUE(all[1].character == Character::StrongAttractive);
    EXPECT_TRUE(all[2].character == Character::Undetermined);
}

UTEST(critical_point_energy, pauli_projected_region_weighting) {
    using namespace critical_point_distribution;
    // The projected surface is x = 0.1 (y^2 + z^2), not the initial disk x = 0.
    const DensityEvaluator density = [](const double* xyz, size_t n, double* rho, double* g, double* h) {
        for (size_t t = 0; t < n; ++t) {
            const double x = xyz[3 * t], y = xyz[3 * t + 1], z = xyz[3 * t + 2];
            const double q = x - 0.1 * (y * y + z * z);
            rho[t] = 0.02 - 0.001 * (y * y + z * z) + 0.002 * q * q;
            g[3 * t] = 0.004 * q;
            g[3 * t + 1] = -0.002 * y - 0.0008 * y * q;
            g[3 * t + 2] = -0.002 * z - 0.0008 * z * q;
            h[6 * t] = 0.004;
            h[6 * t + 1] = -0.0008 * y;
            h[6 * t + 2] = -0.0008 * z;
            h[6 * t + 3] = -0.002 - 0.0008 * q + 0.00016 * y * y;
            h[6 * t + 4] = 0.00016 * y * z;
            h[6 * t + 5] = -0.002 - 0.0008 * q + 0.00016 * z * z;
        }
        return true;
    };
    NciSite site = {{0, 0, 0}, {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}, 0.02};
    NciParams params;
    params.radius = 0.6;
    params.step = 0.2;
    const auto baseline = nci_contacts(std::span(&site, 1), density, params);
    size_t calls = 0, sampled = 0;
    bool projected = true;
    const PauliEvaluator pauli = [&](const double* xyz, size_t n, double* P, double* rho) {
        ++calls;
        sampled = n;
        for (size_t t = 0; t < n; ++t) {
            const double x = xyz[3 * t], y = xyz[3 * t + 1], z = xyz[3 * t + 2];
            projected &= std::abs(x - 0.1 * (y * y + z * z)) < 1e-12;
            rho[t] = 0.02 - 0.001 * (y * y + z * z);
            P[t] = 0.003 + x * 0.01;
        }
        return true;
    };
    const auto result = nci_contacts(std::span(&site, 1), density, params, pauli);
    ASSERT_EQ(size_t(1), result.size());
    EXPECT_EQ(size_t(1), calls);
    EXPECT_TRUE(projected);
    EXPECT_EQ(result[0].samples, sampled);
    EXPECT_EQ(result[0].samples, result[0].pauli.samples);
    EXPECT_TRUE(result[0].pauli.reason.empty());
    EXPECT_TRUE(result[0].character == baseline[0].character);
    EXPECT_EQ(result[0].share, baseline[0].share);
    EXPECT_STREQ(result[0].reason.c_str(), baseline[0].reason.c_str());
    EXPECT_EQ(result[0].attractive, baseline[0].attractive);
    EXPECT_TRUE(std::isnan(baseline[0].pauli.mean));
    double weight = 0.0, sum = 0.0, alpha_sum = 0.0;
    size_t count = 0;
    const double tf = 0.3 * std::pow(3.0 * 9.869604401089358, 2.0 / 3.0);
    const int half = std::max(1, (int)std::floor(params.radius / params.step));
    for (int i = -half; i <= half; ++i) for (int j = -half; j <= half; ++j) {
        const double r2 = (i * i + j * j) * params.step * params.step;
        if (r2 > params.radius * params.radius * (1.0 + 1e-9)) continue;
        const double rho = 0.02 - 0.001 * r2, P = 0.003 + 0.001 * r2;
        weight += rho;
        sum += rho * P;
        alpha_sum += rho * P / (tf * std::pow(rho, 5.0 / 3.0));
        ++count;
    }
    EXPECT_EQ(count, sampled);
    EXPECT_NEAR(sum / weight, result[0].pauli.mean, 1e-14);
    EXPECT_NEAR(alpha_sum / weight, result[0].pauli.alpha, 1e-14);
    EXPECT_GT(result[0].pauli.mean, 0.003); // differs from the value at the CP

    const auto failed = nci_contacts(std::span(&site, 1), density, params,
        [](const double*, size_t, double*, double*) { return false; });
    EXPECT_TRUE(std::isnan(failed[0].pauli.mean));
    EXPECT_EQ(size_t(0), failed[0].pauli.samples);
    EXPECT_FALSE(failed[0].pauli.reason.empty());
    EXPECT_TRUE(failed[0].character == baseline[0].character);
    const auto partial = nci_contacts(std::span(&site, 1), density, params,
        [&](const double* xyz, size_t n, double* P, double* rho) {
            pauli(xyz, n, P, rho);
            P[0] = NAN;
            return true;
        });
    EXPECT_EQ(sampled - 1, partial[0].pauli.samples);
    EXPECT_TRUE(std::isnan(partial[0].pauli.mean));
    EXPECT_TRUE(std::isnan(partial[0].pauli.alpha));
    EXPECT_FALSE(partial[0].pauli.reason.empty());

    const NciSite two_sites[2] = {site, site};
    const auto separate = nci_contacts(two_sites, density, params,
        [&](const double* xyz, size_t n, double* P, double* rho) {
            pauli(xyz, n, P, rho);
            P[0] = NAN; // incomplete first region must not invalidate the second one
            return true;
        });
    EXPECT_EQ(result[0].samples - 1, separate[0].pauli.samples);
    EXPECT_TRUE(std::isnan(separate[0].pauli.mean));
    EXPECT_EQ(result[0].samples, separate[1].pauli.samples);
    EXPECT_NEAR(result[0].pauli.mean, separate[1].pauli.mean, 1e-14);
    EXPECT_NEAR(result[0].pauli.alpha, separate[1].pauli.alpha, 1e-14);
    const auto roundoff = nci_contacts(std::span(&site, 1), density, params,
        [&](const double* xyz, size_t n, double* P, double* rho) {
            pauli(xyz, n, P, rho);
            for (size_t t = 0; t < n; ++t) P[t] = -1e-16;
            return true;
        });
    EXPECT_NEAR(-1e-16, roundoff[0].pauli.mean, 1e-30);
    EXPECT_LT(roundoff[0].pauli.alpha, 0.0);
    calls = 0;
    const auto no_region = nci_contacts(std::span(&site, 1),
        [&](const double* xyz, size_t n, double* rho, double* g, double* h) {
            density(xyz, n, rho, g, h);
            for (size_t t = 0; t < n; ++t) rho[t] = 0.1;
            return true;
        }, params, pauli);
    EXPECT_EQ(size_t(0), calls);
    EXPECT_EQ(size_t(0), no_region[0].samples);
    EXPECT_TRUE(std::isnan(no_region[0].pauli.mean));
    EXPECT_FALSE(no_region[0].pauli.reason.empty());

    site.rho = 0.2;
    calls = 0;
    const auto dense = nci_contacts(std::span(&site, 1), density, params, pauli);
    EXPECT_EQ(size_t(0), calls);
    EXPECT_TRUE(std::isnan(dense[0].pauli.mean));
    EXPECT_FALSE(dense[0].pauli.reason.empty());
}

UTEST(critical_point_energy, nci_projection_reuses_unchanged_and_drops_inactive_points) {
    using namespace critical_point_distribution;
    NciSite site = {{0, 0, 0}, {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}, 0.02};
    NciParams params;
    params.radius = 0.6;
    std::vector<size_t> batches;
    const DensityEvaluator density = [&](const double* xyz, size_t n, double* rho, double* g, double* h) {
        batches.push_back(n);
        for (size_t t = 0; t < n; ++t) {
            const double target = xyz[3 * t + 1] > 0.1 ? 2.0 : 0.2;
            const double q = xyz[3 * t] - target;
            rho[t] = 0.02 + 0.002 * q * q;
            g[3 * t] = 0.004 * q;
            g[3 * t + 1] = g[3 * t + 2] = 0.0;
            std::fill_n(h + 6 * t, 6, 0.0);
            h[6 * t] = 0.004;
            h[6 * t + 3] = h[6 * t + 5] = -0.002;
        }
        return true;
    };
    size_t pauli_points = 0;
    bool projected = true;
    const auto out = nci_contacts(std::span(&site, 1), density, params,
        [&](const double* xyz, size_t n, double* P, double* rho) {
            pauli_points += n;
            for (size_t t = 0; t < n; ++t) {
                projected &= std::abs(xyz[3 * t] - 0.2) < 1e-14 && xyz[3 * t + 1] <= 0.1;
                P[t] = 0.003;
                rho[t] = 0.02;
            }
            return true;
        });
    ASSERT_EQ(size_t(2), batches.size());
    EXPECT_LT(batches[1], batches[0]);
    EXPECT_TRUE(projected);
    EXPECT_GT(pauli_points, size_t(1));
    EXPECT_EQ(batches[1], pauli_points);
    EXPECT_EQ(pauli_points, out[0].samples);
    EXPECT_EQ(pauli_points, out[0].pauli.samples);
    EXPECT_NEAR(0.003, out[0].pauli.mean, 1e-14);
}

UTEST(critical_point_energy, nci_cancels_between_batches) {
    using namespace critical_point_distribution;
    ExpAtoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto site = pair.site({0, 0, 0});
    bool cancel = true;
    size_t calls = 0, pauli_calls = 0;
    const DensityEvaluator density = [&](const double* x, size_t n, double* rho, double* g, double* h) {
        ++calls;
        pair.evaluator()(x, n, rho, g, h);
        cancel = true;
        return true;
    };
    const PauliEvaluator pauli = [&](const double*, size_t, double*, double*) { ++pauli_calls; return true; };
    auto out = nci_contacts(std::span(&site, 1), density, {}, pauli, [&] { return cancel; });
    EXPECT_EQ(size_t(0), calls);
    EXPECT_TRUE(out[0].character == Character::Undetermined);
    EXPECT_TRUE(std::isnan(out[0].pauli.mean));
    EXPECT_NE(out[0].reason.find("cancelled"), std::string::npos);
    cancel = false;
    out = nci_contacts(std::span(&site, 1), density, {}, pauli, [&] { return cancel; });
    EXPECT_EQ(size_t(1), calls);
    EXPECT_EQ(size_t(0), pauli_calls);
    EXPECT_TRUE(out[0].character == Character::Undetermined);
    EXPECT_TRUE(std::isnan(out[0].pauli.mean));
    EXPECT_NE(out[0].pauli.reason.find("cancelled"), std::string::npos);
}

UTEST(critical_point_energy, pauli_gaussian_contact_regions) {
    using namespace critical_point_distribution;
    md_gto_shell_t shells[2] = {{0, 0, 1, 0}, {1, 1, 1, 0}};
    float exponent[2] = {0.5f, 0.5f}, coeff[2] = {0.1f, 0.1f};
    float atoms[6] = {-1.5f, 0, 0, 1.5f, 0, 0};
    double D[4] = {1, 0, 0, 1};
    md_gto_basis_t basis = {2, 2, shells, exponent, coeff};
    md_topo_gto_desc_t desc = {};
    desc.basis = &basis;
    desc.atom_xyz = atoms;
    desc.density_matrix = D;
    desc.num_threads = 1;
    const DensityEvaluator density = [&](const double* x, size_t n, double* rho, double* g, double* h) {
        return md_topo_gto_density_derivatives(rho, g, h, x, n, &desc);
    };
    const PauliEvaluator pauli = [&](const double* x, size_t n, double* P, double* rho) {
        return md_topo_gto_pauli_kinetic_energy_density(P, nullptr, rho, nullptr, x, n, &desc);
    };
    NciSite site = {{0, 0, 0}, {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}, 0};
    double g[3], h[6];
    ASSERT_TRUE(density(site.center, 1, &site.rho, g, h));
    NciParams params;
    params.radius = 0.6;
    params.step = 0.2;
    // This is a numerical model, not an interacting wavefunction or a hydrogen-bond benchmark.
    params.s_max = 2.0;
    const auto independent = nci_contacts(std::span(&site, 1), density, params, pauli);
    ASSERT_GT(independent[0].samples, size_t(1));
    EXPECT_GT(independent[0].pauli.mean, 0.0);
    EXPECT_GT(independent[0].pauli.alpha, 0.0);
    EXPECT_EQ(independent[0].samples, independent[0].pauli.samples);
    // Hold the sampled disk fixed: dilution lowers t_P but raises t_P / tau_TF.
    NciParams fixed_region = params;
    fixed_region.s_max = 1e6;
    const auto full = nci_contacts(std::span(&site, 1), density, fixed_region, pauli);
    D[0] = D[3] = 0.125;
    ASSERT_TRUE(density(site.center, 1, &site.rho, g, h));
    const auto dilute = nci_contacts(std::span(&site, 1), density, fixed_region, pauli);
    ASSERT_GT(full[0].samples, size_t(1));
    ASSERT_EQ(full[0].samples, dilute[0].samples);
    EXPECT_EQ(dilute[0].samples, dilute[0].pauli.samples);
    EXPECT_NEAR(0.125 * full[0].pauli.mean, dilute[0].pauli.mean, 1e-14);
    EXPECT_NEAR(4.0 * full[0].pauli.alpha, dilute[0].pauli.alpha, 1e-10);
    D[0] = D[3] = 1.0;
    // A coherent single orbital has t_P = 0 everywhere, including away from the CP.
    D[1] = D[2] = 1.0;
    ASSERT_TRUE(density(site.center, 1, &site.rho, g, h));
    const auto coherent = nci_contacts(std::span(&site, 1), density, params, pauli);
    ASSERT_GT(coherent[0].samples, size_t(1));
    EXPECT_NEAR(0.0, coherent[0].pauli.mean, 1e-15);
    EXPECT_NEAR(0.0, coherent[0].pauli.alpha, 1e-12);
}
