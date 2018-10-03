#pragma once

#include <mol/molecule_trajectory.h>

// Reads the header info of a trajectory and allocates space for it
// Only supports XTC currently
bool load_and_allocate_trajectory(MoleculeTrajectory* traj, CString path);

// Reads the actual trajectory position information... Necessary?
bool read_trajectory_data(MoleculeTrajectory* traj);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool all_trajectory_frames_read(MoleculeTrajectory* traj);
bool close_file_handle(MoleculeTrajectory* traj);

void copy_trajectory_frame(TrajectoryFrame* dst, const MoleculeTrajectory& src_traj, int frame_index);
void copy_trajectory_positions(Array<vec3> dst_array, const MoleculeTrajectory& traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], const MoleculeTrajectory& traj, int frame_index);

inline TrajectoryFrame get_trajectory_frame(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(frame_index < traj.num_frames);
    return traj.frame_buffer.data[frame_index];
}

inline Array<vec3> get_trajectory_positions(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(frame_index < traj.num_frames);
    return traj.frame_buffer.data[frame_index].atom_positions;
}
