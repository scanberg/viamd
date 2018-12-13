#include "mol/hydrogen_bond.h"
#include <core/common.h>
#include <core/log.h>
#include <mol/trajectory_utils.h>
#include <mol/spatial_hash.h>

namespace hydrogen_bond {

// Computes the potential donors given a set of atom labels.
// OH and NH atoms are assumed to be donors if the concecutive atom is marked with 'H' for Hydrogen.
int32 compute_donors(DynamicArray<HydrogenBondDonor>* donors, Array<const Label> labels) {
    ASSERT(donors);
    int32 pre_count = (int32)donors->count;
    const int32 num_labels = (int32)labels.count;
    for (int32 i = 0; i < num_labels; i++) {
        if (compare_n(labels[i], "OH", 2) || compare_n(labels[i], "NH", 2)) {
            if (i + 1 < num_labels && compare_n(labels[i + 1], "H", 1)) {
                donors->push_back({i, i + 1});
            }
            if (i + 2 < num_labels && compare_n(labels[i + 2], "H", 1)) {
                donors->push_back({i, i + 2});
            }
        }
    }
    return (int32)donors->count - pre_count;
}

DynamicArray<HydrogenBondDonor> compute_donors(Array<const Label> labels) {
    DynamicArray<HydrogenBondDonor> donors;
    compute_donors(&donors, labels);
    return donors;
}

// Computes the potential acceptors given a set of atom elements.
// This essentially just a filter on atom element which extracts Oxygen and Nitrogen
int32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements) {
    ASSERT(acceptors);
    const int32 pre_count = (int32)acceptors->count;
    for (int32 i = 0; i < (int32)elements.count; i++) {
        if (elements[i] == Element::O || elements[i] == Element::N) {
            acceptors->push_back(i);
        }
    }
    return (int32)acceptors->count - pre_count;
}

DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements) {
    DynamicArray<HydrogenBondAcceptor> acceptors;
    compute_acceptors(&acceptors, elements);
    return acceptors;
}

// Computes hydrogen bonds given a certain set of potential donors, acceptors and atomic positions from a frame.
// The distance cutoff sets the distance from bonds to potential acceptors.
//

int32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                    Array<const vec3> atom_positions, float dist_cutoff, float angle_cutoff) {
    const int32 num_acceptors = (int32)acceptors.count;
    if (!num_acceptors) return 0;

    DynamicArray<vec3> acceptor_pos(num_acceptors);
    DynamicArray<AtomIdx> acceptor_idx(num_acceptors);
    for (int32 i = 0; i < num_acceptors; i++) {
        acceptor_pos[i] = atom_positions[acceptors[i]];
        acceptor_idx[i] = acceptors[i];
    }

    int32 pre_count = (int32)bonds->count;
    spatialhash::Frame frame = spatialhash::compute_frame(acceptor_pos, vec3(dist_cutoff));
    for (const auto& don : donors) {
        vec3 donor_pos = atom_positions[don.donor_idx];
        vec3 hydro_pos = atom_positions[don.hydro_idx];
        spatialhash::for_each_within(frame, hydro_pos, dist_cutoff,
                                     [bonds, &donor_pos, &hydro_pos, &acceptor_idx, &don, angle_cutoff](int32 idx, const vec3& pos) {
                                         AtomIdx g_idx = acceptor_idx[idx];
                                         if (g_idx == don.donor_idx) return;
                                         const vec3 a = hydro_pos - donor_pos;
                                         const vec3 b = pos - hydro_pos;
                                         if (math::angle(a, b) < angle_cutoff) {
                                             bonds->push_back({g_idx, don.donor_idx, don.hydro_idx});
                                         }
                                     });
    }
    return (int32)bonds->count - pre_count;
}

DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                                         Array<const vec3> atom_positions, float dist_cutoff, float angle_cutoff) {
    DynamicArray<HydrogenBond> bonds;
    compute_bonds(&bonds, donors, acceptors, atom_positions, dist_cutoff, angle_cutoff);
    return bonds;
}

void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    ASSERT(hbt);
    for (int32 i = 0; i < dyn.trajectory.num_frames; i++) {
        Array<HydrogenBond> frame_bonds{(HydrogenBond*)(hbt->bond_data.end()), int64(0)};
        Array<const vec3> atom_positions = get_trajectory_positions(dyn.trajectory, i);
        frame_bonds.count = compute_bonds(&hbt->bond_data, dyn.molecule.hydrogen_bond.donors, dyn.molecule.hydrogen_bond.acceptors, atom_positions,
                                          max_dist, max_angle);
    }
}

HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    HydrogenBondTrajectory hbt;
    compute_bonds_trajectory(&hbt, dyn, max_dist, max_angle);
    return hbt;
}
}  // namespace hydrogen_bond