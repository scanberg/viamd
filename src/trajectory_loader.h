#pragma once

#include <core/string_types.h>
#include <task_system.h>

struct MoleculeStructure;
struct MoleculeTrajectory;

namespace load {
namespace mol {
    bool load_molecule(MoleculeStructure* mol, CStringView filename);
    task::TaskID load_molecule_async(MoleculeStructure* mol, CStringView filename);
}

namespace traj {
    bool load_trajectory(MoleculeTrajectory* traj, CStringView filename);
    task::TaskID load_trajectory_async(MoleculeTrajectory* traj, CStringView filename);
}

}  // namespace load