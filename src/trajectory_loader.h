#pragma once

#include <core/types.h>
#include <core/string_types.h>

struct MoleculeStructure;
struct MoleculeTrajectory;

namespace load {
namespace mol {
    bool is_extension_supported(CStringView filename);
    bool load_molecule(MoleculeStructure* mol, CStringView filename);
}

namespace traj {
    bool is_extension_supported(CStringView filename);
    bool read_num_frames(i32* num_frames, CStringView filename);
    bool load_trajectory(MoleculeTrajectory* traj, i32 num_atoms, CStringView filename);
}

}  // namespace load