#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

struct md_system_t;
struct md_topo_extremum_graph_t;

namespace critical_point_distribution {

using Position = std::array<double, 3>;
bool within_focus(const Position& position_angstrom, std::span<const Position> centers_angstrom);

// What a critical point is assigned to. Unassigned keeps the CP's topology in its group name.
enum class Kind : uint8_t {
    Atom,       // an attractor at a unique nucleus
    Bond,       // a BCP whose bond path joins two bonded atoms
    Contact,    // a BCP whose bond path joins two atoms without a bond between them
    BondPath,   // a BCP joining two atoms when the connectivity is unknown
    Ring,       // an RCP inside a unique ring: a perceived ring, or a cycle of bond paths
    Cage,       // a CCP whose ring-CP faces are all assigned rings
    NonNuclear, // an attractor more than 0.50 Angstrom from every nucleus
    Unassigned,
};
inline constexpr int kind_count = 8;
const char* kind_name(Kind kind);

// lambda1 / lambda2 - 1 of the two negative Hessian eigenvalues (ascending) of a BCP
double ellipticity(const double hessian_eigenvalues[3]);

struct ChemicalLabel {
    Kind kind = Kind::Unassigned;
    std::string group;
    std::string detail;
    std::string interaction;   // a contact's class ("Hydrogen bond", "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", ...); empty otherwise
    int atom[2] = {-1, -1};
};

// "Bond: C – O"; unassigned groups already name their topology ("Unassigned BCP")
std::string display_name(const ChemicalLabel& label);

// ring_surfaces, by vertex, is optional: for a ring CP, the maxima its ring surface reaches
// (md_topo_gto_ring_surface_maxima), empty when not traced. They decide between equally small cycles.
std::vector<ChemicalLabel> chemical_labels(const md_topo_extremum_graph_t& graph,
                                         const md_system_t& system,
                                         std::span<const Position> atoms_angstrom,
                                         std::span<const std::vector<uint32_t>> ring_surfaces = {});

// A ring CP left between equally small cycles, which its ring surface may decide
bool ring_surface_may_decide(const ChemicalLabel& label);

struct Range {
    double min = 0.0;
    double max = 1.0;
};

// An enabled filter accepts the union of its inclusive ranges; a disabled one accepts every finite density.
struct Filter {
    bool enabled = false;
    std::vector<Range> ranges;

    bool accepts(double density) const;
};

struct Curve {
    static constexpr int sample_count = 256;
    std::array<double, sample_count> density = {};
    std::array<double, sample_count> probability = {};
    // Per component, the share of probability from its members; the components sum to probability.
    std::vector<std::array<double, sample_count>> components;
    bool logarithmic = false;
    bool valid = false;
};

// Probability per log10-density unit when all values are positive, otherwise per density unit.
// component[i] in [0, component_count) assigns density[i] to a component of the returned curve.
Curve estimate(std::span<const double> density, std::span<const uint32_t> component = {}, size_t component_count = 0);

}
