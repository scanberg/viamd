#pragma once

#include <mol/trajectory.h>
#include <core/allocator.h>

Trajectory read_and_allocate_trajectory(const char* file, Allocator& alloc = default_alloc);
void free_trajectory(Trajectory* traj);

TrajectoryFrame allocate_trajectory_frame(int num_atoms, Allocator& alloc = default_alloc);
void free_trajectory_frame(TrajectoryFrame* frame);

void copy_trajectory_frame(TrajectoryFrame* dst, const Trajectory& src_traj, int frame_index);
void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index);