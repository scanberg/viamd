#include "critical_point_distribution_core.h"

#include <core/md_log.h>
#include <core/md_unit.h>
#include <md_system.h>
#include <md_topo.h>
#include <md_util.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <vector>

namespace critical_point_distribution {

static std::string element(const md_system_t& sys, size_t i) {
    const str_t symbol = md_atomic_number_symbol(md_atom_atomic_number(&sys.atom, i));
    return symbol.len ? std::string(symbol.ptr, symbol.len) : "Unknown element";
}

// The atom's first bonded neighbour, -1 without connectivity or bonds
static int bonded_partner(const md_system_t& sys, size_t i) {
    if (!sys.bond.conn.offset || i + 1 >= sys.bond.conn.offset_count) return -1;
    for (auto it = md_bond_iter(&sys.bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        return (int)md_bond_iter_atom_index(&it);
    }
    return -1;
}

static std::string environment(const md_system_t& sys, size_t i) {
    const auto z = md_atom_atomic_number(&sys.atom, i);
    std::string label = element(sys, i);
    if (z == 1) {
        const int parent = bonded_partner(sys, i);
        return parent >= 0 ? label + " (on " + element(sys, parent) + ")" : label;
    }
    if (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_AROMATIC) return label + " aromatic";
    bool carbonyl = false, hydrogen = false;
    for (auto it = md_bond_iter(&sys.bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const auto other_z = md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&it));
        hydrogen |= other_z == 1;
        const auto flags = sys.bond.flags ? sys.bond.flags[md_bond_iter_bond_index(&it)] : MD_BOND_FLAG_NONE;
        carbonyl |= md_bond_order(flags) == 2 && ((z == 6 && other_z == 8) || (z == 8 && other_z == 6));
    }
    if (carbonyl) return label + " carbonyl";
    hydrogen |= sys.atom.hydrogen_count && sys.atom.hydrogen_count[i] > 0;
    if ((z == 7 || z == 8) && hydrogen) return label + " (H-bearing)";
    return label;
}

// Hill-order formula of the atoms: "C5N"
static std::string formula(const md_system_t& sys, const int32_t* atoms, size_t size) {
    std::vector<std::pair<std::string, int>> counts;
    for (size_t k = 0; k < size; ++k) {
        const std::string e = element(sys, atoms[k]);
        auto it = std::find_if(counts.begin(), counts.end(), [&](const auto& c) { return c.first == e; });
        if (it == counts.end()) counts.push_back({e, 1});
        else it->second++;
    }
    // Hill order: C first, then the others alphabetically
    std::sort(counts.begin(), counts.end(), [](const auto& a, const auto& b) {
        if ((a.first == "C") != (b.first == "C")) return a.first == "C";
        return a.first < b.first;
    });
    std::string result;
    for (const auto& c : counts) result += c.first + (c.second > 1 ? std::to_string(c.second) : "");
    return result;
}

// The ring's kind, size and composition: "Aromatic 6-ring C5N"
static std::string ring_name(const md_system_t& sys, size_t ring) {
    const md_ring_flags_t flags = md_ring_flags(&sys.ring, ring);
    const size_t size = md_ring_size(&sys.ring, ring);
    const std::string formula = critical_point_distribution::formula(sys, md_ring_atoms(&sys.ring, ring), size);
    const char* kind = (flags & MD_RING_FLAG_AROMATIC) ? "Aromatic" : (flags & MD_RING_FLAG_SUGAR) ? "Sugar" : "Aliphatic";
    return std::string(kind) + " " + std::to_string(size) + "-ring " + formula;
}

const char* kind_name(Kind kind) {
    switch (kind) {
    case Kind::Atom:       return "Atom";
    case Kind::Bond:       return "Bond";
    case Kind::Contact:    return "Contact";
    case Kind::BondPath:   return "Bond path";
    case Kind::Ring:       return "Ring";
    case Kind::Cage:       return "Cage";
    case Kind::NonNuclear: return "Non-nuclear attractor";
    case Kind::Unassigned: return "Unassigned";
    }
    return "Unassigned";
}

std::string display_name(const ChemicalLabel& label) {
    if (label.kind == Kind::Unassigned) return label.group;
    return std::string(kind_name(label.kind)) + ": " + label.group;
}

std::string display_name(const ChemicalLabel& label, int level) {
    const int finer = std::min(level, (int)label.finer.size());
    if (finer <= 0) return display_name(label);
    const std::string& group = label.finer[finer - 1];
    if (label.kind == Kind::Unassigned) return group;
    return std::string(kind_name(label.kind)) + ": " + group;
}

int auto_detail_level(std::span<const ChemicalLabel> labels, std::span<const uint32_t> members, size_t max_groups) {
    std::vector<std::string> names;
    for (int level = detail_levels - 1; level > 0; --level) {
        names.clear();
        for (uint32_t i : members) {
            if (i < labels.size()) names.push_back(display_name(labels[i], level));
        }
        std::sort(names.begin(), names.end());
        if ((size_t)(std::unique(names.begin(), names.end()) - names.begin()) <= max_groups) return level;
    }
    return 0;
}

// The finest name a label has so far
static const std::string& finest(const ChemicalLabel& label) { return label.finer.empty() ? label.group : label.finer.back(); }

// Adds a finer name, unless it says nothing more than the last
static void add_detail(ChemicalLabel& label, const std::string& name) {
    if (!name.empty() && name != finest(label)) label.finer.push_back(name);
}

// The atoms bonded to atom i, as a formula ("C2H"); empty without connectivity
static std::string neighbours(const md_system_t& sys, size_t i) {
    if (!sys.bond.conn.offset || i + 1 >= sys.bond.conn.offset_count) return {};
    std::vector<int32_t> atoms;
    for (auto it = md_bond_iter(&sys.bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) atoms.push_back((int32_t)md_bond_iter_atom_index(&it));
    return atoms.empty() ? std::string() : formula(sys, atoms.data(), atoms.size());
}

// The atom's environment with what it is bonded to: "C aromatic [C2H]"; a hydrogen by its parent: "H (on N [C2H])"
static std::string atom_detail(const md_system_t& sys, size_t i) {
    if (md_atom_atomic_number(&sys.atom, i) == 1) {
        const int parent = bonded_partner(sys, i);
        if (parent < 0) return environment(sys, i);
        std::string p = environment(sys, parent);
        const std::string implied = " (H-bearing)";
        if (p.size() > implied.size() && p.compare(p.size() - implied.size(), implied.size(), implied) == 0) p.resize(p.size() - implied.size());
        const std::string n = neighbours(sys, parent);
        return "H (on " + p + (n.empty() ? "" : " [" + n + "]") + ")";
    }
    const std::string n = neighbours(sys, i);
    return environment(sys, i) + (n.empty() ? "" : " [" + n + "]");
}

// The residues of the atoms by name: "in PHE", "between GLY and SER", "across DG\xC3\x974"; empty without residue names
static std::string residue_context(const md_system_t& sys, std::span<const int> atoms) {
    std::vector<md_residue_idx_t> residues;
    for (int a : atoms) {
        if (a < 0) return {};
        const md_residue_idx_t r = md_system_residue_find_by_atom_idx(&sys, (size_t)a);
        if (r < 0 || md_system_residue_name(&sys, r).len == 0) return {};
        residues.push_back(r);
    }
    std::sort(residues.begin(), residues.end());
    residues.erase(std::unique(residues.begin(), residues.end()), residues.end());
    if (residues.empty()) return {};
    std::vector<std::string> names;
    for (md_residue_idx_t r : residues) {
        const str_t name = md_system_residue_name(&sys, r);
        names.push_back(std::string(name.ptr, name.len));
    }
    std::sort(names.begin(), names.end());
    if (names.size() == 1) return "in " + names[0];
    if (names.size() == 2) return "between " + names[0] + " and " + names[1];
    std::string result = "across ";
    for (size_t k = 0; k < names.size();) {
        size_t n = 1;
        while (k + n < names.size() && names[k + n] == names[k]) ++n;
        result += (k ? ", " : "") + names[k] + (n > 1 ? "\xC3\x97" + std::to_string(n) : "");
        k += n;
    }
    return result;
}

// Appends a residue context as the finest level
static void add_residue_detail(ChemicalLabel& label, const md_system_t& sys, std::span<const int> atoms) {
    const std::string context = residue_context(sys, atoms);
    if (!context.empty()) add_detail(label, finest(label) + " " + context);
}

static double distance_squared(const Position& a, const Position& b) {
    double sum = 0.0;
    for (int k = 0; k < 3; ++k) sum += (a[k] - b[k]) * (a[k] - b[k]);
    return sum;
}

double ellipticity(const double ev[3]) {
    return ev[1] != 0.0 ? ev[0] / ev[1] - 1.0 : 0.0;
}

static double angle_degrees(const Position& a, const Position& vertex, const Position& b) {
    double dot = 0.0, la = 0.0, lb = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double u = a[k] - vertex[k], v = b[k] - vertex[k];
        dot += u * v; la += u * u; lb += v * v;
    }
    if (la <= 0.0 || lb <= 0.0) return 0.0;
    return std::acos(std::clamp(dot / std::sqrt(la * lb), -1.0, 1.0)) * 180.0 / 3.14159265358979323846;
}

static bool share_ring(const md_system_t& sys, int a, int b) {
    for (size_t r = 0; r < md_ring_data_count(&sys.ring); ++r) {
        const int32_t* atoms = md_ring_atoms(&sys.ring, r);
        const size_t size = md_ring_size(&sys.ring, r);
        bool has_a = false, has_b = false;
        for (size_t k = 0; k < size; ++k) { has_a |= atoms[k] == a; has_b |= atoms[k] == b; }
        if (has_a && has_b) return true;
    }
    return false;
}

// The class of a contact bond path between atoms a and b, from elements, connectivity and geometry. The angle
// criteria are the usual ones (IUPAC 2011 for hydrogen bonds: D-H...A above 110 degrees; halogen bonds: C-X...Y
// above 140 degrees); a contact failing them is "Other contact". A class names the geometry and topology only.
static std::string contact_class(const md_system_t& sys, std::span<const Position> xyz, int a, int b, std::string& why) {
    auto z = [&](int i) { return (int)md_atom_atomic_number(&sys.atom, i); };
    auto ring_aromatic = [&](int i) { return (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_AROMATIC) != 0; };
    // Of an aromatic system: in an aromatic ring, or sp2/sp and bonded to one (exocyclic C=O, NH2 of a nucleobase)
    auto aromatic = [&](int i) {
        if (ring_aromatic(i)) return true;
        const md_hybridization_t h = md_atom_hybridization(&sys.atom, i);
        if (z(i) == 1 || (h != MD_HYBRIDIZATION_SP2 && h != MD_HYBRIDIZATION_SP) || !sys.bond.conn.offset) return false;
        for (auto it = md_bond_iter(&sys.bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            if (ring_aromatic((int)md_bond_iter_atom_index(&it))) return true;
        }
        return false;
    };
    auto acceptor = [](int n) { return n == 7 || n == 8 || n == 9 || n == 16 || n == 17 || n == 35 || n == 53; };
    auto donor = [](int n) { return n == 7 || n == 8 || n == 9 || n == 16; };
    auto halogen = [](int n) { return n == 17 || n == 35 || n == 53; };
    char buf[128];

    if (z(a) == 1 && z(b) == 1) { why = "Both ends are hydrogens."; return "H\xC2\xB7\xC2\xB7\xC2\xB7H"; }
    if (z(a) == 1 || z(b) == 1) {
        const int h = z(a) == 1 ? a : b, other = h == a ? b : a;
        const int parent = bonded_partner(sys, h);
        if (parent < 0) { why = "The hydrogen has no bonded atom."; return "Other contact"; }
        const double angle = angle_degrees(xyz[parent], xyz[h], xyz[other]);
        snprintf(buf, sizeof(buf), "%s-H\xC2\xB7\xC2\xB7\xC2\xB7%s angle %.0f degrees", element(sys, parent).c_str(), element(sys, other).c_str(), angle);
        why = buf;
        if (angle < 110.0) { why += " (below 110)."; return "Other contact"; }
        why += ".";
        const bool polar = donor(z(parent));
        if (polar && acceptor(z(other))) return "Hydrogen bond";
        if (ring_aromatic(other)) return polar ? "X-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80" : z(parent) == 6 ? "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80" : "Other contact";
        if (z(parent) == 6 && acceptor(z(other))) return "C-H\xC2\xB7\xC2\xB7\xC2\xB7X";
        return "Other contact";
    }
    for (int pass = 0; pass < 2; ++pass) {
        const int x = pass ? b : a, y = pass ? a : b;
        if (!halogen(z(x)) || !(z(y) == 7 || z(y) == 8 || z(y) == 16)) continue;
        const int parent = bonded_partner(sys, x);
        if (parent < 0 || z(parent) != 6) continue;
        const double angle = angle_degrees(xyz[parent], xyz[x], xyz[y]);
        snprintf(buf, sizeof(buf), "C-%s\xC2\xB7\xC2\xB7\xC2\xB7%s angle %.0f degrees", element(sys, x).c_str(), element(sys, y).c_str(), angle);
        why = buf;
        if (angle >= 140.0) { why += "."; return "Halogen bond"; }
    }
    if (aromatic(a) && aromatic(b) && !share_ring(sys, a, b)) { why = "Both ends belong to aromatic systems (ring atoms, or conjugated substituents) of different rings."; return "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80"; }
    if (why.empty()) why = "No class criterion is met.";
    return "Other contact";
}

// Whether p, projected on the plane of the closed polygon, lies inside it; out_of_plane is its distance from it.
static bool inside_polygon(const Position& p, const std::vector<Position>& polygon, double& out_of_plane, double& area) {
    const size_t n = polygon.size();
    Position c = {0, 0, 0}, normal = {0, 0, 0};
    for (const auto& q : polygon) for (int k = 0; k < 3; ++k) c[k] += q[k] / n;
    for (size_t i = 0; i < n; ++i) {
        const Position& u = polygon[i];
        const Position& v = polygon[(i + 1) % n];
        normal[0] += (u[1] - v[1]) * (u[2] + v[2]);
        normal[1] += (u[2] - v[2]) * (u[0] + v[0]);
        normal[2] += (u[0] - v[0]) * (u[1] + v[1]);
    }
    const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (length < 1e-9) return false;
    area = 0.5 * length;
    for (double& x : normal) x /= length;
    // An orthonormal basis (e1, e2) of the plane
    Position e1 = std::abs(normal[0]) < 0.9 ? Position{1, 0, 0} : Position{0, 1, 0};
    double d = e1[0] * normal[0] + e1[1] * normal[1] + e1[2] * normal[2];
    for (int k = 0; k < 3; ++k) e1[k] -= d * normal[k];
    const double l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (double& x : e1) x /= l1;
    const Position e2 = {normal[1] * e1[2] - normal[2] * e1[1], normal[2] * e1[0] - normal[0] * e1[2], normal[0] * e1[1] - normal[1] * e1[0]};
    auto project = [&](const Position& q, double& x, double& y) {
        const double r[3] = {q[0] - c[0], q[1] - c[1], q[2] - c[2]};
        x = r[0] * e1[0] + r[1] * e1[1] + r[2] * e1[2];
        y = r[0] * e2[0] + r[1] * e2[1] + r[2] * e2[2];
    };
    out_of_plane = std::abs((p[0] - c[0]) * normal[0] + (p[1] - c[1]) * normal[1] + (p[2] - c[2]) * normal[2]);
    double px, py;
    project(p, px, py);
    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi, yi, xj, yj;
        project(polygon[i], xi, yi);
        project(polygon[j], xj, yj);
        if ((yi > py) != (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
}

struct PathEdge {
    int atom;
    Kind kind;
    uint32_t bcp;
};

// The simple cycles of 3 to max_size atoms of the bond-path graph among the candidate atoms, each once. False when
// the search exceeds its step budget (dense contact networks), in which case nothing is concluded.
static bool bond_path_cycles(const std::vector<std::vector<PathEdge>>& adjacency, const std::vector<int>& candidates,
                             size_t max_size, std::vector<std::vector<int>>& cycles) {
    std::vector<char> is_candidate(adjacency.size(), 0), on_path(adjacency.size(), 0);
    for (int a : candidates) is_candidate[a] = 1;
    std::vector<int> path;
    size_t steps = 0;
    const size_t budget = 200000;
    bool complete = true;
    for (int start : candidates) {
        std::function<void(int)> visit = [&](int atom) {
            if (++steps > budget) { complete = false; return; }
            for (const PathEdge& e : adjacency[atom]) {
                if (!complete) return;
                if (e.atom == start && path.size() >= 3 && path[1] < path.back()) cycles.push_back(path);
                if (e.atom <= start || !is_candidate[e.atom] || on_path[e.atom] || path.size() >= max_size) continue;
                on_path[e.atom] = 1; path.push_back(e.atom);
                visit(e.atom);
                on_path[e.atom] = 0; path.pop_back();
            }
        };
        path = {start};
        on_path[start] = 1;
        visit(start);
        on_path[start] = 0;
        if (!complete) return false;
    }
    return true;
}

static const char tied_ring_detail[] = "Equally small bond-path cycles of about the same area enclose it: ";

bool ring_surface_may_decide(const ChemicalLabel& label) {
    return label.kind == Kind::Unassigned && label.detail.rfind(tied_ring_detail, 0) == 0;
}

std::vector<ChemicalLabel> chemical_labels(const md_topo_extremum_graph_t& graph,
                                         const md_system_t& system,
                                         std::span<const Position> atoms_angstrom,
                                         std::span<const std::vector<uint32_t>> ring_surfaces) {
    std::vector<ChemicalLabel> labels(graph.num_vertices);
    const double bohr_to_angstrom = md_unit_scale(md_unit_bohr_radius()) / md_unit_scale(md_unit_angstrom());
    const size_t num_atoms = std::min(system.atom.count, atoms_angstrom.size());
    auto position = [&](uint32_t i) {
        const auto& v = graph.vertices[i];
        const auto* d = graph.vertices_f64 ? graph.vertices_f64 + i : nullptr;
        return Position{(d ? d->x : v.x) * bohr_to_angstrom, (d ? d->y : v.y) * bohr_to_angstrom, (d ? d->z : v.z) * bohr_to_angstrom};
    };

    // Rings whose atoms all have coordinates, by their centroids
    std::vector<std::pair<size_t, Position>> rings;
    for (size_t r = 0; r < md_ring_data_count(&system.ring); ++r) {
        const int32_t* atoms = md_ring_atoms(&system.ring, r);
        const size_t size = md_ring_size(&system.ring, r);
        Position centroid = {0.0, 0.0, 0.0};
        bool complete = size > 0;
        for (size_t k = 0; k < size && complete; ++k) {
            if (atoms[k] < 0 || (size_t)atoms[k] >= num_atoms) { complete = false; break; }
            for (int c = 0; c < 3; ++c) centroid[c] += atoms_angstrom[atoms[k]][c] / size;
        }
        if (complete) rings.push_back({r, centroid});
    }

    std::vector<int> assigned_atom(system.atom.count, -1);
    std::vector<int> assigned_ring(md_ring_data_count(&system.ring), -1);
    // The atoms of each assigned ring CP, for its cage
    std::vector<std::vector<int>> ring_atoms(graph.num_vertices);
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        auto& label = labels[i];
        const int type = md_topo_vertex_type(&graph, i);
        if (type == MD_TOPO_MINIMUM) {
            label.group = "Cage CP (unassigned)";
            label.detail = "Cage critical points are not assigned to a cage.";
            continue;
        }
        if (type == MD_TOPO_JOIN_SADDLE) {
            label.group = "Unassigned ring CP";
            label.detail = "No unique perceived ring (3-7 atoms) with its centroid within 0.50 Angstrom.";
            const Position p = position(i);
            int match = -1;
            bool ambiguous = false;
            for (const auto& ring : rings) {
                if (distance_squared(p, ring.second) <= 0.5 * 0.5) {
                    ambiguous |= match >= 0;
                    match = (int)ring.first;
                }
            }
            if (match < 0 || ambiguous) continue;
            if (assigned_ring[match] != -1) {
                const int previous = assigned_ring[match];
                if (previous >= 0) {
                    labels[previous] = {Kind::Unassigned, "Unassigned ring CP", "Multiple ring CPs at the centroid of the same ring."};
                }
                assigned_ring[match] = -2;
                label.detail = "Multiple ring CPs at the centroid of the same ring.";
                continue;
            }
            assigned_ring[match] = (int)i;
            ring_atoms[i].assign(md_ring_atoms(&system.ring, match), md_ring_atoms(&system.ring, match) + md_ring_size(&system.ring, match));
            label.kind = Kind::Ring;
            label.group = ring_name(system, match);
            const int32_t* atoms = md_ring_atoms(&system.ring, match);
            std::string list;
            for (size_t k = 0; k < md_ring_size(&system.ring, match); ++k) list += (k ? ", " : "") + std::to_string(atoms[k] + 1);
            label.detail = "Ring of atoms " + list + "; its centroid is the only one within 0.50 Angstrom.";
            add_residue_detail(label, system, ring_atoms[i]);
            continue;
        }
        if (type != MD_TOPO_MAXIMUM) continue;
        const Position p = position(i);
        int match = -1;
        bool ambiguous = false;
        for (size_t j = 0; j < num_atoms; ++j) {
            if (distance_squared(p, atoms_angstrom[j]) <= 0.3 * 0.3) {
                if (match >= 0) ambiguous = true;
                match = (int)j;
            }
        }
        label.group = "Unassigned attractor";
        label.detail = "No unique nucleus within 0.30 Angstrom; may be non-nuclear or unresolved.";
        if (match < 0 && num_atoms >= 2) {
            // A non-nuclear attractor only when clearly away from every nucleus (0.30 to 0.50 Angstrom stays unresolved)
            int first = -1, second = -1;
            for (size_t j = 0; j < num_atoms; ++j) {
                const double d = distance_squared(p, atoms_angstrom[j]);
                if (first < 0 || d < distance_squared(p, atoms_angstrom[first])) { second = first; first = (int)j; }
                else if (second < 0 || d < distance_squared(p, atoms_angstrom[second])) second = (int)j;
            }
            const double d1 = std::sqrt(distance_squared(p, atoms_angstrom[first]));
            const double d2 = std::sqrt(distance_squared(p, atoms_angstrom[second]));
            if (d1 > 0.5) {
                std::string ea = element(system, first), eb = element(system, second);
                if (ea > eb) std::swap(ea, eb);
                char buf[160];
                snprintf(buf, sizeof(buf), "Nearest nuclei: atom %d at %.2f Angstrom and atom %d at %.2f Angstrom.", first + 1, d1, second + 1, d2);
                label.kind = Kind::NonNuclear;
                label.group = "near " + ea + ", " + eb;
                label.detail = buf;
            }
            continue;
        }
        if (match < 0 || ambiguous) continue;
        // Multiple attractors near one nucleus are not multiple confident atom assignments.
        if (assigned_atom[match] != -1) {
            const int previous = assigned_atom[match];
            if (previous >= 0) {
                labels[previous] = {Kind::Unassigned, "Unassigned attractor", "Multiple attractors within 0.30 Angstrom of the same nucleus."};
            }
            assigned_atom[match] = -2;
            label.detail = "Multiple attractors within 0.30 Angstrom of the same nucleus.";
            continue;
        }
        assigned_atom[match] = (int)i;
        label.kind = Kind::Atom;
        label.atom[0] = match;
        label.group = environment(system, match);
        label.detail = "Atom " + std::to_string(match + 1) + " (" + label.group + "); unique nucleus within 0.30 Angstrom.";
        add_detail(label, atom_detail(system, match));
        const int self[1] = {match};
        add_residue_detail(label, system, self);
    }
    std::vector<std::vector<uint32_t>> endpoints(graph.num_vertices);
    for (uint32_t e = 0; e < graph.num_edges; ++e) {
        const auto edge = graph.edges[e];
        if (edge.from >= graph.num_vertices || edge.to >= graph.num_vertices) continue;
        uint32_t bcp = edge.from, attractor = edge.to;
        if (md_topo_vertex_type(&graph, bcp) == MD_TOPO_MAXIMUM) std::swap(bcp, attractor);
        if (md_topo_vertex_type(&graph, bcp) != MD_TOPO_SPLIT_SADDLE ||
            md_topo_vertex_type(&graph, attractor) != MD_TOPO_MAXIMUM) continue;
        auto& ends = endpoints[bcp];
        if (std::find(ends.begin(), ends.end(), attractor) == ends.end()) ends.push_back(attractor);
    }
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (md_topo_vertex_type(&graph, i) != MD_TOPO_SPLIT_SADDLE) continue;
        auto& label = labels[i];
        label.group = "Unassigned BCP";
        label.detail = "Requires two distinct bond-path attractors, each assigned to a nucleus.";
        const auto& ends = endpoints[i];
        if (ends.size() != 2) continue;
        int a = labels[ends[0]].atom[0], b = labels[ends[1]].atom[0];
        if (a < 0 || b < 0 || a == b) continue;
        std::string ea = environment(system, a), eb = environment(system, b);
        if (ea > eb) { std::swap(a, b); std::swap(ea, eb); }
        const bool connectivity = system.bond.conn.offset != nullptr;
        const bool bonded = connectivity && md_bond_find(&system.bond, a, b) >= 0;
        label.kind = !connectivity ? Kind::BondPath : bonded ? Kind::Bond : Kind::Contact;
        const char* relation = !connectivity ? "connectivity unknown" : bonded ? "bonded in the topology" : "not bonded in the topology";
        label.atom[0] = a; label.atom[1] = b;
        label.group = ea + (label.kind == Kind::Bond ? " \xE2\x80\x93 " : label.kind == Kind::Contact ? " \xC2\xB7\xC2\xB7\xC2\xB7 " : " ~ ") + eb;
        label.detail = "Bond paths to atoms " + std::to_string(a + 1) + " and " + std::to_string(b + 1) +
                       ", " + relation + ".";
        if (label.kind == Kind::Contact && (size_t)std::max(a, b) < num_atoms) {
            std::string why;
            std::string cls = contact_class(system, atoms_angstrom, a, b, why);
            label.interaction = cls;
            const md_molecule_idx_t ma = md_system_molecule_find_by_atom_idx(&system, a);
            const md_molecule_idx_t mb = md_system_molecule_find_by_atom_idx(&system, b);
            if (ma >= 0 && mb >= 0) cls += ma == mb ? " (intramolecular)" : " (intermolecular)";
            label.group = cls + ": " + label.group;
            label.detail += " " + why;
        }
        {
            const char* sep = label.kind == Kind::Bond ? " \xE2\x80\x93 " : label.kind == Kind::Contact ? " \xC2\xB7\xC2\xB7\xC2\xB7 " : " ~ ";
            const std::string coarse = ea + sep + eb;
            const std::string prefix = label.group.substr(0, label.group.size() - coarse.size());
            std::string da = atom_detail(system, a), db = atom_detail(system, b);
            if (ea == eb && da > db) std::swap(da, db);
            add_detail(label, prefix + da + sep + db);
            const int ends_atoms[2] = {a, b};
            add_residue_detail(label, system, ends_atoms);
        }
        if (graph.vertices_f64) {
            const double* ev = graph.vertices_f64[i].hessian_eigenvalues;
            const double laplacian = ev[0] + ev[1] + ev[2];
            char buf[160];
            snprintf(buf, sizeof(buf), " Laplacian %.4g au (%s), ellipticity %.3f.", laplacian,
                     laplacian < 0.0 ? "negative: shared-shell" : "positive: closed-shell", ellipticity(ev));
            label.detail += buf;
        }
    }

    // Ring CPs not at a perceived ring's centroid: inside the smallest cycle of bond paths (bonds and contacts) around
    // them. Bond paths are taken as straight lines between the nuclei.
    auto density = [&](uint32_t i) { return graph.vertices_f64 ? graph.vertices_f64[i].value : (double)graph.vertices[i].value; };
    // One edge per pair of nuclei: of parallel bond paths, the densest BCP, which the density rule below should see
    std::vector<std::vector<PathEdge>> adjacency(num_atoms);
    auto add_edge = [&](int a, int b, Kind kind, uint32_t bcp) {
        for (PathEdge& e : adjacency[a]) {
            if (e.atom != b) continue;
            if (density(bcp) > density(e.bcp)) e = {b, kind, bcp};
            return;
        }
        adjacency[a].push_back({b, kind, bcp});
    };
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        const auto& l = labels[i];
        if (md_topo_vertex_type(&graph, i) != MD_TOPO_SPLIT_SADDLE || l.atom[0] < 0 || l.atom[1] < 0) continue;
        if ((size_t)std::max(l.atom[0], l.atom[1]) >= num_atoms) continue;
        add_edge(l.atom[0], l.atom[1], l.kind, i);
        add_edge(l.atom[1], l.atom[0], l.kind, i);
    }
    auto edge = [&](int a, int b) {
        for (const PathEdge& e : adjacency[a]) if (e.atom == b) return e;
        return PathEdge{b, Kind::Unassigned, UINT32_MAX};
    };
    // Every enclosing cycle of every remaining ring CP is a candidate, scored by its size and the CP's distance from
    // its plane (0.25 Angstrom counts as one atom). The best-scoring pairs are matched first and each cycle takes one
    // ring CP, so a CP whose best cycle is taken falls back to its next one rather than both being dropped.
    std::vector<std::vector<int>> cycle_atoms;
    std::map<std::vector<int>, int> cycle_index;   // by sorted atoms
    std::vector<char> cycle_claimed;
    auto cycle_id = [&](const std::vector<int>& cycle) {
        std::vector<int> key = cycle;
        std::sort(key.begin(), key.end());
        auto [it, inserted] = cycle_index.try_emplace(key, (int)cycle_atoms.size());
        if (inserted) { cycle_atoms.push_back(cycle); cycle_claimed.push_back(0); }
        return it->second;
    };
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (md_topo_vertex_type(&graph, i) == MD_TOPO_JOIN_SADDLE && !ring_atoms[i].empty()) cycle_claimed[cycle_id(ring_atoms[i])] = 1;
    }
    struct Candidate { uint32_t rcp; int cycle; double offset, area, score; };
    std::vector<Candidate> candidates;
    std::vector<char> pending(graph.num_vertices, 0);
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (md_topo_vertex_type(&graph, i) != MD_TOPO_JOIN_SADDLE || labels[i].kind != Kind::Unassigned) continue;
        if (labels[i].detail.rfind("Multiple", 0) == 0) continue;
        const Position p = position(i);
        // Rings of up to 10 bond paths first; only without one, the larger hydrogen-bonded rings of proteins, up to 14
        bool thinner = false, crowded = false;
        for (const auto [max_size, radius] : {std::pair<size_t, double>{10, 5.0}, {14, 6.0}}) {
            std::vector<int> near;
            for (size_t j = 0; j < num_atoms; ++j) {
                if (!adjacency[j].empty() && distance_squared(p, atoms_angstrom[j]) <= radius * radius) near.push_back((int)j);
            }
            std::vector<std::vector<int>> cycles;
            if (!bond_path_cycles(adjacency, near, max_size, cycles)) {
                crowded = true;
                break;
            }
            for (const auto& cycle : cycles) {
                std::vector<Position> polygon;
                for (int a : cycle) polygon.push_back(atoms_angstrom[a]);
                double offset = 0.0, area = 0.0;
                if (!inside_polygon(p, polygon, offset, area) || offset > 1.0) continue;
                // On its ring's surface a ring CP is the density minimum: every bond path around it is denser at its BCP
                bool denser = true;
                for (size_t k = 0; k < cycle.size() && denser; ++k) {
                    const PathEdge e = edge(cycle[k], cycle[(k + 1) % cycle.size()]);
                    denser = e.bcp != UINT32_MAX && density(e.bcp) > density(i);
                }
                if (!denser) { thinner = true; continue; }
                const int id = cycle_id(cycle);
                if (cycle_claimed[id]) continue;
                candidates.push_back({i, id, offset, area, (double)cycle.size() + offset / 0.25});
                pending[i] = 1;
            }
            if (pending[i]) break;
        }
        labels[i].detail = crowded && !pending[i] ? "Too many bond-path cycles around it to decide its ring."
                         : thinner ? "The cycles of 3-14 bond paths around it each have a BCP less dense than it, so none can be its ring."
                         : "No cycle of 3-14 bond paths encloses it within 1.0 Angstrom of its plane.";
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        return a.score != b.score ? a.score < b.score : a.rcp < b.rcp;
    });
    std::vector<std::pair<int, double>> chosen(graph.num_vertices, {-1, 0.0});
    std::vector<std::string> surface_note(graph.num_vertices);
    auto atom_list = [](const std::vector<int>& atoms) {
        std::string list;
        for (int a : atoms) list += (list.empty() ? "" : ", ") + std::to_string(a + 1);
        return list;
    };
    for (const Candidate& c : candidates) {
        if (!pending[c.rcp] || cycle_claimed[c.cycle]) continue;
        // Equally good cycles of the same size: one nested in another shares most of its path and detours around
        // more, so the ring CP sits in the face of the smallest. Undecided when the two smallest are about as large.
        std::vector<const Candidate*> rivals = {&c};
        for (const Candidate& o : candidates) {
            if (o.rcp != c.rcp || o.cycle == c.cycle || cycle_claimed[o.cycle]) continue;
            if (cycle_atoms[o.cycle].size() == cycle_atoms[c.cycle].size() && std::abs(o.offset - c.offset) < 0.1) rivals.push_back(&o);
        }
        std::sort(rivals.begin(), rivals.end(), [](const Candidate* a, const Candidate* b) { return a->area < b->area; });
        pending[c.rcp] = 0;
        // A ring's face holds no other bond path: a cycle with another's atom inside it encloses that one
        if (rivals.size() > 1) {
            std::vector<const Candidate*> empty;
            for (const Candidate* x : rivals) {
                const std::vector<int>& own = cycle_atoms[x->cycle];
                std::vector<Position> polygon;
                for (int a : own) polygon.push_back(atoms_angstrom[a]);
                bool holds = false;
                for (const Candidate* y : rivals) {
                    for (int a : cycle_atoms[y->cycle]) {
                        if (holds || std::find(own.begin(), own.end(), a) != own.end()) continue;
                        double off = 0.0, area = 0.0;
                        holds = inside_polygon(atoms_angstrom[a], polygon, off, area) && off < 1.0;
                    }
                }
                if (!holds) empty.push_back(x);
            }
            if (empty.size() == 1) rivals = empty;
        }
        // The ring surface of a ring CP ends at the nuclei of its ring: of the cycles, only one may hold all those it reaches.
        // Its paths can stray to a neighbouring nucleus (the angles to some nuclei are very narrow), so the nuclei in no
        // cycle are left out.
        if (rivals.size() > 1 && rivals[1]->area < 1.15 * rivals[0]->area && c.rcp < ring_surfaces.size() && !ring_surfaces[c.rcp].empty()) {
            auto holds = [&](const Candidate* x, int a) {
                const std::vector<int>& own = cycle_atoms[x->cycle];
                return std::find(own.begin(), own.end(), a) != own.end();
            };
            std::vector<int> reached;
            for (uint32_t v : ring_surfaces[c.rcp]) {
                if (v >= graph.num_vertices || labels[v].kind != Kind::Atom) continue;
                const int a = labels[v].atom[0];
                const bool in_cycle = std::any_of(rivals.begin(), rivals.end(), [&](const Candidate* x) { return holds(x, a); });
                if (in_cycle && std::find(reached.begin(), reached.end(), a) == reached.end()) reached.push_back(a);
            }
            std::vector<const Candidate*> holding;
            for (const Candidate* x : rivals) {
                if (std::all_of(reached.begin(), reached.end(), [&](int a) { return holds(x, a); })) holding.push_back(x);
            }
            if (!reached.empty() && holding.size() == 1) {
                surface_note[c.rcp] = " Of " + std::to_string(rivals.size()) + " equally small cycles around it, the only one with every nucleus its ring surface reaches (atoms " + atom_list(reached) + ").";
                rivals = holding;
            } else if (!reached.empty()) {
                surface_note[c.rcp] = " Its ring surface reaches atoms " + atom_list(reached) + ", in more than one of them.";
            }
        }
        if (rivals.size() > 1 && rivals[1]->area < 1.15 * rivals[0]->area) {
            auto describe = [&](const Candidate& x) {
                std::string list;
                for (int a : cycle_atoms[x.cycle]) list += (list.empty() ? "" : ", ") + std::to_string(a + 1);
                char buf[64];
                snprintf(buf, sizeof(buf), " (%.2f Angstrom off, %.1f Angstrom\xC2\xB2)", x.offset, x.area);
                return "atoms " + list + buf;
            };
            std::string text = tied_ring_detail;
            for (size_t k = 0; k < rivals.size(); ++k) text += (k ? "; " : "") + describe(*rivals[k]);
            labels[c.rcp].detail = text + "." + surface_note[c.rcp];
            continue;
        }
        const Candidate& best = *rivals[0];
        cycle_claimed[best.cycle] = 1;
        chosen[best.rcp] = {best.cycle, best.offset};
    }
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (chosen[i].first < 0) {
            if (pending[i]) labels[i].detail = "Every bond-path cycle around it is taken by a nearer ring CP.";
            continue;
        }
        const std::vector<int>& best = cycle_atoms[chosen[i].first];
        const double best_offset = chosen[i].second;
        int contacts = 0, unknown = 0;
        uint32_t coalescing = UINT32_MAX;
        std::vector<std::string> classes;
        for (size_t k = 0; k < best.size(); ++k) {
            const PathEdge e = edge(best[k], best[(k + 1) % best.size()]);
            contacts += e.kind == Kind::Contact;
            unknown += e.kind == Kind::BondPath;
            if (e.kind == Kind::Contact && e.bcp != UINT32_MAX) {
                const std::string& cls = labels[e.bcp].interaction.empty() ? std::string("Other contact") : labels[e.bcp].interaction;
                if (std::find(classes.begin(), classes.end(), cls) == classes.end()) classes.push_back(cls);
            }
            // A ring CP with nearly the density of a BCP on its ring is close to annihilating with it
            if (e.bcp != UINT32_MAX && std::abs(density(e.bcp) - density(i)) <= 0.02 * std::abs(density(e.bcp))) coalescing = e.bcp;
        }
        // The ring is grouped by the most specific of its closing contacts; all of them are in its detail
        static const char* priority[] = {"Hydrogen bond", "Halogen bond", "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80",
                                         "X-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80",
                                         "C-H\xC2\xB7\xC2\xB7\xC2\xB7X", "H\xC2\xB7\xC2\xB7\xC2\xB7H", "Other contact"};
        auto rank = [](const std::string& cls) {
            for (size_t k = 0; k < sizeof(priority) / sizeof(priority[0]); ++k) if (cls == priority[k]) return k;
            return sizeof(priority) / sizeof(priority[0]);
        };
        std::sort(classes.begin(), classes.end(), [&](const auto& a, const auto& b) { return rank(a) < rank(b); });
        std::string closed_by;
        for (const auto& cls : classes) closed_by += (closed_by.empty() ? "" : ", ") + cls;
        const std::string size = std::to_string(best.size()) + "-ring";
        const std::string atoms_formula = formula(system, best.data(), best.size());
        auto& label = labels[i];
        label.kind = Kind::Ring;
        // Contact rings are named by what closes them; their sizes and compositions vary too much to group by
        label.group = contacts ? "Closed by " + classes.front()
                    : unknown ? "Bond-path " + size + " " + atoms_formula
                    : "Bonded " + size + " " + atoms_formula;
        std::string list;
        for (size_t k = 0; k < best.size(); ++k) list += (k ? ", " : "") + std::to_string(best[k] + 1);
        char buf[160];
        snprintf(buf, sizeof(buf), "; %.2f Angstrom from its mean plane", best_offset);
        label.detail = "Cycle of bond paths around it: atoms " + list + " (" + atoms_formula + ")" + buf;
        if (contacts) label.detail += "; " + size + ", " + std::to_string(contacts) + (contacts == 1 ? " contact: " : " contacts: ") + closed_by + ".";
        else label.detail += ".";
        label.detail += surface_note[i];
        if (contacts) add_detail(label, label.group + ", " + size + " " + atoms_formula);
        add_residue_detail(label, system, best);
        if (coalescing != UINT32_MAX) {
            label.group += " (near-degenerate)";
            for (auto& name : label.finer) name += " (near-degenerate)";
            snprintf(buf, sizeof(buf), " Its density is within 2%% of BCP %u on the ring: the pair is close to annihilating.", coalescing);
            label.detail += buf;
        }
        ring_atoms[i] = best;
    }

    // Cage CPs: the ring CPs whose paths end at them are the faces of the cage
    std::vector<std::vector<uint32_t>> faces(graph.num_vertices);
    for (uint32_t e = 0; e < graph.num_edges; ++e) {
        uint32_t rcp = graph.edges[e].from, ccp = graph.edges[e].to;
        if (rcp >= graph.num_vertices || ccp >= graph.num_vertices) continue;
        if (md_topo_vertex_type(&graph, rcp) == MD_TOPO_MINIMUM) std::swap(rcp, ccp);
        if (md_topo_vertex_type(&graph, rcp) != MD_TOPO_JOIN_SADDLE || md_topo_vertex_type(&graph, ccp) != MD_TOPO_MINIMUM) continue;
        if (std::find(faces[ccp].begin(), faces[ccp].end(), rcp) == faces[ccp].end()) faces[ccp].push_back(rcp);
    }
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (md_topo_vertex_type(&graph, i) != MD_TOPO_MINIMUM) continue;
        const auto& f = faces[i];
        if (f.size() < 3) {
            labels[i].detail = "Cage critical points need ring-CP paths to at least 3 faces; " + std::to_string(f.size()) + " found.";
            continue;
        }
        std::vector<int> atoms;
        std::string list;
        bool complete = true;
        for (uint32_t r : f) {
            if (labels[r].kind != Kind::Ring || ring_atoms[r].empty()) { complete = false; break; }
            atoms.insert(atoms.end(), ring_atoms[r].begin(), ring_atoms[r].end());
            list += (list.empty() ? "" : "; ") + labels[r].group;
        }
        if (!complete) {
            labels[i].detail = "Not every ring CP whose path ends here is assigned to a ring.";
            continue;
        }
        std::sort(atoms.begin(), atoms.end());
        atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
        bool contact_face = false;
        for (uint32_t r : f) contact_face |= labels[r].group.rfind("Closed by ", 0) == 0;
        labels[i].kind = Kind::Cage;
        labels[i].group = std::string(contact_face ? "Contact cage" : "Bonded cage") + ", " + std::to_string(f.size()) + " faces";
        labels[i].detail = "Atoms " + formula(system, atoms.data(), atoms.size()) + "; faces: " + list + ".";
        add_detail(labels[i], labels[i].group + " " + formula(system, atoms.data(), atoms.size()));
        add_residue_detail(labels[i], system, atoms);
    }
    return labels;
}

bool within_focus(const Position& position_angstrom, std::span<const Position> centers_angstrom) {
    for (const Position& center : centers_angstrom) {
        if (distance_squared(position_angstrom, center) <= 3.0 * 3.0) return true;
    }
    return false;
}

bool Filter::accepts(double density) const {
    if (!std::isfinite(density)) return false;
    if (!enabled) return true;
    for (const Range& range : ranges) {
        if (density >= range.min && density <= range.max) return true;
    }
    return false;
}

Curve estimate(std::span<const double> density, std::span<const uint32_t> component, size_t component_count) {
    Curve curve;
    if (density.empty()) return curve;
    for (double value : density) {
        if (!std::isfinite(value)) {
            MD_LOG_ERROR("Critical-point distribution: non-finite density");
            return curve;
        }
    }
    if (!component.empty()) {
        if (component.size() != density.size()) {
            MD_LOG_ERROR("Critical-point distribution: component count does not match density count");
            return curve;
        }
        for (uint32_t c : component) {
            if (c >= component_count) {
                MD_LOG_ERROR("Critical-point distribution: component index out of range");
                return curve;
            }
        }
    }

    std::vector<size_t> order(density.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return density[a] < density[b]; });
    std::vector<double> values(density.size());
    std::vector<uint32_t> member(density.size(), 0);
    for (size_t i = 0; i < order.size(); ++i) {
        values[i] = density[order[i]];
        if (!component.empty()) member[i] = component[order[i]];
    }
    curve.components.assign(component.empty() ? 0 : component_count, {});
    curve.logarithmic = values.front() > 0.0;
    if (curve.logarithmic) {
        for (double& value : values) value = std::log10(value);
    }

    const size_t count = values.size();
    const size_t neighbors = std::min(count - 1, std::min(size_t(32), std::max(size_t(2), size_t(std::sqrt(count)))));
    const double fallback = curve.logarithmic ? 0.05 : std::max(std::abs(values.front()) * 0.05, 1.0e-12);
    std::vector<double> bandwidth(count);
    for (size_t i = 0; i < count; ++i) {
        size_t lo = i, hi = i;
        for (size_t k = 0; k < neighbors; ++k) {
            if (lo > 0 && (hi + 1 == count || values[i] - values[lo - 1] < values[hi + 1] - values[i])) --lo;
            else ++hi;
        }
        bandwidth[i] = std::max(values[i] - values[lo], values[hi] - values[i]) * 0.5;
    }
    const double max_bandwidth = std::max(*std::max_element(bandwidth.begin(), bandwidth.end()), fallback);
    const double begin = values.front() - 4.0 * max_bandwidth;
    const double end = values.back() + 4.0 * max_bandwidth;
    const double step = (end - begin) / (Curve::sample_count - 1);
    // Keep even coincident samples resolved by the curve's sampling grid.
    for (double& width : bandwidth) width = std::max(width, 2.0 * step);

    constexpr double inv_sqrt_two_pi = 0.3989422804014327;
    for (int j = 0; j < Curve::sample_count; ++j) {
        const double x = begin + j * step;
        curve.density[j] = curve.logarithmic ? std::pow(10.0, x) : x;
        double sum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            const double z = (x - values[i]) / bandwidth[i];
            const double kernel = std::exp(-0.5 * z * z) / bandwidth[i];
            sum += kernel;
            if (!curve.components.empty()) curve.components[member[i]][j] += kernel * inv_sqrt_two_pi / count;
        }
        curve.probability[j] = sum * inv_sqrt_two_pi / count;
    }
    curve.valid = true;
    return curve;
}


std::vector<Valley> valleys(const Curve& curve, double min_persistence) {
    std::vector<Valley> result;
    if (!curve.valid) return result;
    const auto& p = curve.probability;
    const int n = (int)p.size();
    double top = 0.0;
    for (double v : p) top = std::max(top, v);
    if (!(top > 0.0)) return result;
    // Superlevel sets from the top down: a sample joining two peaks is a valley, and the lower peak dies there
    std::vector<int> order(n);
    for (int k = 0; k < n; ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return p[a] > p[b]; });
    std::vector<int> parent(n, -1);
    std::vector<double> peak(n, 0.0);
    std::function<int(int)> root = [&](int k) { return parent[k] == k ? k : parent[k] = root(parent[k]); };
    for (int k : order) {
        parent[k] = k;
        peak[k] = p[k];
        int joined = -1;
        for (int nb : {k - 1, k + 1}) {
            if (nb < 0 || nb >= n || parent[nb] < 0) continue;
            const int r = root(nb);
            if (joined < 0) {
                parent[k] = r;
                joined = r;
                continue;
            }
            const int a = root(k);
            if (a == r) continue;
            const int younger = peak[a] < peak[r] ? a : r, elder = younger == a ? r : a;
            const double persistence = (peak[younger] - p[k]) / top;
            if (persistence >= min_persistence) result.push_back({curve.density[k], persistence});
            parent[younger] = elder;
        }
    }
    std::sort(result.begin(), result.end(), [](const Valley& a, const Valley& b) { return a.density < b.density; });
    return result;
}

std::vector<double> peak_cuts(std::span<const Valley> valleys, double view_min, double view_max, size_t max_peaks) {
    std::vector<double> cuts;
    if (max_peaks < 2) return cuts;
    auto deepest = [&](auto&& accept) {
        std::vector<Valley> chosen;
        for (const Valley& v : valleys) if (accept(v)) chosen.push_back(v);
        std::stable_sort(chosen.begin(), chosen.end(), [](const Valley& a, const Valley& b) { return a.persistence > b.persistence; });
        if (chosen.size() > max_peaks - 1) chosen.resize(max_peaks - 1);
        return chosen;
    };
    // Outside the view only the overview's cuts
    for (const Valley& v : deepest([](const Valley&) { return true; })) {
        if (v.density < view_min || v.density > view_max) cuts.push_back(v.density);
    }
    for (const Valley& v : deepest([&](const Valley& v) { return v.density >= view_min && v.density <= view_max; })) cuts.push_back(v.density);
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
}

std::vector<double> refine_cuts(std::span<const double> overview, std::span<const Valley> view_valleys,
                                std::span<const double> sorted_values, double view_min, double view_max, size_t max_peaks) {
    std::vector<double> cuts(overview.begin(), overview.end());
    std::sort(cuts.begin(), cuts.end());
    auto in_view = [&](double d) { return d >= view_min && d <= view_max; };
    size_t peaks_in_view = 1;
    for (double c : cuts) peaks_in_view += in_view(c);
    std::vector<Valley> candidates;
    for (const Valley& v : view_valleys) if (in_view(v.density)) candidates.push_back(v);
    std::stable_sort(candidates.begin(), candidates.end(), [](const Valley& a, const Valley& b) { return a.persistence > b.persistence; });
    auto rank = [&](double d) { return std::upper_bound(sorted_values.begin(), sorted_values.end(), d) - sorted_values.begin(); };
    for (const Valley& v : candidates) {
        if (peaks_in_view >= max_peaks) break;
        // Both sides must keep CPs: the values between the neighbouring cuts and this one
        const auto it = std::lower_bound(cuts.begin(), cuts.end(), v.density);
        const auto r = rank(v.density);
        const auto below = it == cuts.begin() ? 0 : rank(*(it - 1));
        const auto above = it == cuts.end() ? (std::ptrdiff_t)sorted_values.size() : rank(*it);
        if (r == below || r == above) continue;
        cuts.insert(it, v.density);
        ++peaks_in_view;
    }
    return cuts;
}

}
