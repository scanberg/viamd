#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
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
    // The group at increasing detail, each a refinement of the one before: the bonded neighbours of the atoms
    // ("C aromatic [C2H]"), a contact ring's size and composition, then the residues ("in PHE"). May be empty.
    std::vector<std::string> finer;
};

// Levels of detail of a group: 0 is ChemicalLabel::group, then its finer names
inline constexpr int detail_levels = 3;

// "Bond: C – O"; unassigned groups already name their topology ("Unassigned BCP")
std::string display_name(const ChemicalLabel& label);
// The same at a level of detail, the finest the label has up to it
std::string display_name(const ChemicalLabel& label, int level);
// The finest level of detail at which the members' labels fall into at most max_groups groups; 0 if none does
int auto_detail_level(std::span<const ChemicalLabel> labels, std::span<const uint32_t> members, size_t max_groups);

// Stable chemistry branches above an environment leaf, independent of detail, character and density peaks.
std::vector<std::string> environment_path(const ChemicalLabel& label, const md_system_t& system);

struct EnvironmentTreeNode {
    std::string name;
    std::vector<uint32_t> children;
    std::vector<uint32_t> environments;  // descendant input indices, including collapsed leaves
    int environment = -1;               // input index for a leaf, -1 for a branch
};

// Root at index 0. Every path ends in a leaf; identical leaf names remain separate rows.
std::vector<EnvironmentTreeNode> environment_tree(std::span<const std::vector<std::string>> paths);

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

// Probability per log10-density unit when all values are positive and allow_log, otherwise per density unit.
// component[i] in [0, component_count) assigns density[i] to a component of the returned curve.
Curve estimate(std::span<const double> density, std::span<const uint32_t> component = {}, size_t component_count = 0,
               bool allow_log = true);

// The energy densities at a CP (Hartree / Bohr^3) from the kinetic energy density G and the Laplacian, by the local
// virial theorem: V = lap / 4 - 2 G, H = G + V. ratio is |V| / G: below 1 closed-shell, above 2 shared-shell.
struct EnergyDensities {
    double G = NAN, V = NAN, H = NAN, ratio = NAN;
};
EnergyDensities energy_densities(double G, double laplacian);

// The energy of a hydrogen bond from its BCP, kcal/mol, negative for binding. Emamian et al. (J. Comput. Chem. 2019)
// fit the binding energies of neutral complexes to rho alone: -223.08 rho + 0.7423. Espinosa, Molins and Lecomte
// (Chem. Phys. Lett. 1998): V / 2, known to overestimate strong ones.
enum class HBondFormula : uint8_t { Emamian, Espinosa };
double hydrogen_bond_energy(HBondFormula formula, double rho, double V);
// Only for what the formulas were fitted on: a BCP whose contact is classified a hydrogen bond
bool hydrogen_bond_energy_applies(const ChemicalLabel& label);

// Whether the interaction at a BCP is attractive or repulsive, from the density alone, by a non-covalent
// interaction (NCI) analysis (Johnson et al., JACS 2010; Contreras-Garcia et al., JCTC 2011). Around the BCP, the
// low reduced density gradient region s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)) <= s_max at rho <= rho_max is
// split by the sign of the Hessian's second eigenvalue: lambda2 < 0 where the density accumulates between the
// atoms (attractive), lambda2 > 0 where it is depleted (repulsive, steric crowding, as towards a ring CP). Where rho
// is below rho_weak, lambda2 is near 0: NCI's green, van der Waals (dispersion) band, mostly dispersion, which the
// density does not capture. A region mostly of such density is Weak, with a lean only when its rho-weighted repulsive
// share is clearly to one side; otherwise the repulsive share of its denser part gives the character. A BCP denser
// than rho_max is outside NCI's range: Covalent when its atoms are bonded, StrongAttractive when not (a short, strong
// hydrogen bond). Density shape, not an energy: see approximate_interaction_energy for a sign of that.
enum class Character : uint8_t { None, Covalent, StrongAttractive, Attractive, Mixed, Repulsive, WeakAttractive, Weak, WeakRepulsive, Undetermined };
inline constexpr int character_count = 10;
inline constexpr double nci_mixed_from = 0.25;      // repulsive share
inline constexpr double nci_repulsive_from = 0.5;
inline constexpr double nci_weak_attractive_below = 0.35;   // repulsive share of a weak region: a lean only outside
inline constexpr double nci_weak_repulsive_from = 0.65;     // this wide neutral band
const char* character_name(Character character);

struct NciParams {
    double radius     = 2.0;    // Bohr: of the disk sampled across the bond path
    double step       = 0.2;    // Bohr: its grid spacing
    double s_max      = 0.5;
    double rho_max    = 0.05;   // au
    double rho_weak   = 0.01;   // au: below it, the van der Waals band
    int    projections = 2;     // Newton steps onto the density minimum along the bond path
    double max_offset = 1.0;    // Bohr: a sample projected further is dropped
};
// A BCP: its position (Bohr), density and Hessian eigenvectors, ascending (axes[2] along the bond path)
struct NciSite {
    double center[3];
    double axes[3][3];
    double rho;
    bool   bonded = false;   // its atoms are bonded (or the connectivity is unknown): Covalent rather than StrongAttractive
};
struct PauliRegion {
    double mean = NAN;        // rho-weighted mean t_P, Hartree / Bohr^3
    double alpha = NAN;       // rho-weighted mean of t_P / tau_TF; tau_TF = 3/10 (3 pi^2)^(2/3) rho^(5/3)
    size_t samples = 0;       // valid samples; means require every sample in the NCI region to be valid
    std::string reason = "Pauli sampling was not requested.";
};
struct NciContact {
    Character   character = Character::None;
    double      attractive = NAN, repulsive = NAN;   // sum of rho dA over the region at rho >= rho_weak, by the sign of lambda2
    double      weak = NAN;                          // ... and at rho < rho_weak
    double      weak_repulsive = NAN;                // the part of weak where lambda2 > 0
    double      share = NAN;                         // the repulsive share of the part that decided: repulsive / (attractive
                                                     // + repulsive), or weak_repulsive / weak for a Weak region
    size_t      samples = 0;                         // in the region
    double      repulsive_peak = NAN;                // the largest rho where lambda2 > 0 (0 when none)
    std::string reason;
    PauliRegion pauli = {};
};
// Evaluates rho, its gradient (3 per point) and Hessian (xx xy xz yy yz zz) at num points xyz (Bohr); false if it fails
using DensityEvaluator = std::function<bool(const double* xyz, size_t num, double* rho, double* grad, double* hess)>;
// Total-density t_P and rho at arbitrary points (Bohr); false on failure or cancellation.
using PauliEvaluator = std::function<bool(const double* xyz, size_t num, double* pauli, double* rho)>;
// The NCI character of each site. Its samples are a disk of the given radius across the bond path, projected along it
// onto the density minimum: the surface between the atoms, where the region lies. The region is the samples
// connected to the BCP within the thresholds. All sites are batched together; after each projection only moved,
// still-active points are reevaluated (at most 1 + projections density batches).
// Optional Pauli sampling uses only the accepted, projected region in one additional batch. These are local,
// density-weighted means over a projected disk, not a 3D interaction energy or a repulsion classification.
// Optional cancellation is checked between evaluation batches and during regional processing. Unfinished contacts
// become Undetermined with a reason; regional Pauli values are invalidated. Evaluators must also honour cancellation.
std::vector<NciContact> nci_contacts(std::span<const NciSite> sites, const DensityEvaluator& evaluate, const NciParams& params = {},
                                     const PauliEvaluator& evaluate_pauli = {}, const std::function<bool()>& cancelled = {});
double reduced_density_gradient(double rho, const double grad[3]);
// Of a symmetric 3 x 3 matrix given as its upper triangle (xx xy xz yy yz zz), ascending
std::array<double, 3> symmetric_eigenvalues(const double h[6]);

// Mulliken populations of num_atoms atoms, N_A = sum_s sum_{mu in A} (P_s S)_mu,mu, of the densities given (the
// total alone, or alpha and beta). The charge of A is its nuclear charge less N_A.
std::vector<double> mulliken_populations(std::span<const double> S, std::span<const double* const> P, size_t num_ao,
                                         std::span<const uint32_t> ao_atom, size_t num_atoms, const std::function<bool()>& cancelled = {});

// An approximate interaction energy of two atoms, kcal/mol, negative favourable: the Coulomb energy of their
// charges plus an exchange term from their bond order, q_A q_B / R - B_AB / (2 R), R in Bohr. A crude stand-in
// for the IQA pair energy: point charges, basis dependent Mulliken charges, and no deformation energy.
double approximate_interaction_energy(double qa, double qb, double bond_order, double distance_bohr);

// Mayer bond orders of atom pairs, B_AB = w sum_s sum_{mu in A, nu in B} (P_s S)_mu,nu (P_s S)_nu,mu: the total density
// alone with w = 1, or the alpha and beta densities with w = 2. S and every P_s are num_ao x num_ao, row major, in one
// AO basis (Cartesian or spherical: the order is invariant); ao_atom is each AO's atom. A pair with an atom that has
// no AOs is NaN. Pairs are computed in parallel; the optional cancellation callback must be thread-safe.
// Unfinished pairs remain NaN. Cancelled Mulliken evaluation above leaves all populations NaN.
std::vector<double> mayer_bond_orders(std::span<const double> S, std::span<const double* const> P, double w, size_t num_ao,
                                      std::span<const uint32_t> ao_atom, std::span<const std::array<int, 2>> pairs,
                                      const std::function<bool()>& cancelled = {});

// A local minimum of a curve's total probability, where it splits into peaks. Its persistence is how far the
// smaller of the two peaks it separates rises above it, as a fraction of the curve's maximum: a deep valley parts
// two major peaks, a shallow one two parts of one. Valleys below the threshold are noise and not reported.
struct Valley {
    double density = 0.0;
    double persistence = 0.0;
};
std::vector<Valley> valleys(const Curve& curve, double min_persistence = 0.02);

// The densities at which a curve is cut into peaks, ascending. The overview is the max_peaks - 1 most persistent
// valleys; within the view [view_min, view_max] the max_peaks - 1 most persistent valleys there. As the view
// narrows, its cuts are the overview's and then ever shallower ones: each peak a part of one seen before.
std::vector<double> peak_cuts(std::span<const Valley> valleys, double view_min, double view_max, size_t max_peaks);

// A zoomed view's cuts: every overview cut, and within [view_min, view_max] the most persistent valleys of the curve
// estimated on the view's CPs alone (view_valleys), until the view holds max_peaks peaks. sorted_values are all the
// densities; a valley that parts no CPs from the cuts already made is the same gap and skipped.
std::vector<double> refine_cuts(std::span<const double> overview, std::span<const Valley> view_valleys,
                                std::span<const double> sorted_values, double view_min, double view_max, size_t max_peaks);

}
