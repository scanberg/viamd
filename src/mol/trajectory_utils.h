#pragma once

#include <mol/trajectory.h>
#include <core/allocator.h>

// Reads the header info of a trajectory and allocates space for it
Trajectory* allocate_trajectory(const char* file);
void free_trajectory(Trajectory* traj);

// Reads the actual trajectory position information
void read_trajectory_async(Trajectory* traj, void(*on_finish) = nullptr);

TrajectoryFrame allocate_trajectory_frame(int num_atoms, Allocator* alloc = nullptr);
void free_trajectory_frame(TrajectoryFrame* frame);

void copy_trajectory_frame(TrajectoryFrame* dst, const Trajectory& src_traj, int frame_index);
void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index);

TrajectoryFrame get_trajectory_frame(const Trajectory& traj, int frame_index);
Array<vec3> get_trajectory_positions(const Trajectory& traj, int frame_index);