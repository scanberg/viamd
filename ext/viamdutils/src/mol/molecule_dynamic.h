#pragma once

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>

struct MoleculeDynamic {
    MoleculeStructure molecule{};
    MoleculeTrajectory trajectory{};

    operator bool() const {
        bool mol_ok = molecule.operator bool();
        bool traj_ok = trajectory.operator bool();
        return mol_ok && traj_ok;
    }
};
