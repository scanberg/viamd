#pragma once

#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>

struct HydrogenBond {
    AtomIdx acc_idx = 0;
    AtomIdx don_idx = 0;
    AtomIdx hyd_idx = 0;
};

struct HydrogenBondTrajectory {
    DynamicArray<HydrogenBond> bond_data{};
    DynamicArray<Array<HydrogenBond>> frame_bonds{};
};

namespace hydrogen_bond {
int32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements);
DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements);

int32 compute_donors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Label> labels);
DynamicArray<HydrogenBondDonor> compute_donors(Array<const Label> labels);

int32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                    Array<const vec3> atom_positions, float dist_cutoff = 3.f, float angle_cutoff = 20.f * math::DEG_TO_RAD);
DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                                         Array<const vec3> atom_positions, float dist_cutoff = 3.f, float angle_cutoff = 20.f * math::DEG_TO_RAD);

void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
}  // namespace hydrogen_bond