#include "utest.h"
#include <qm/critical_point_distribution_core.h>
#include <md_system.h>
#include <md_topo.h>
#include <md_util.h>
#include <core/md_unit.h>
#include <core/md_allocator.h>

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
