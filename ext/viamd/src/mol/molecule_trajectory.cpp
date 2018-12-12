#include "molecule_trajectory.h"

bool init_trajectory(MoleculeTrajectory* traj, int32 num_atoms, int32 num_frames) {
    ASSERT(traj);

    traj->num_atoms = num_atoms;
    traj->num_frames = num_frames;
    traj->total_simulation_time = 0;
    traj->simulation_type = MoleculeTrajectory::NVT;
    traj->path_to_file = {};
    traj->file_handle = nullptr;

    traj->frame_offsets = {};
    traj->position_data = {(vec3*)CALLOC(num_frames * num_atoms, sizeof(vec3)), num_frames * num_atoms};
    traj->frame_buffer = {(TrajectoryFrame*)CALLOC(num_frames, sizeof(TrajectoryFrame)), num_frames};

    return true;
}

void free_trajectory(MoleculeTrajectory* traj) {
    ASSERT(traj);

    free_string(&traj->path_to_file);
    if (traj->frame_offsets.data) FREE(traj->frame_offsets.data);
    if (traj->position_data.data) FREE(traj->position_data.data);
    if (traj->frame_buffer.data) FREE(traj->frame_buffer.data);

    *traj = {};
}
