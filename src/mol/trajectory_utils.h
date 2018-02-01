#pragma once

#include <mol/trajectory.h>

Trajectory read_trajectory(const char* file);
void destroy_trajectory(Trajectory* traj);

TrajectoryFrame get_trajectory_frame(Trajectory traj, int frame_index);
void read_trajectory_positions(Array<vec3> atom_positions, Trajectory traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], Trajectory traj, int frame_index);