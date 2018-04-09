#pragma once

#include <mol/trajectory.h>
#include <thread>

// Reads the header info of a trajectory and allocates space for it
bool init_trajectory(Trajectory* traj, CString path);

// Allocates space for trajectory
bool init_trajectory(Trajectory* traj, int32 num_atoms, int32 num_frames);

// Frees memory allocated by trajectory
void free_trajectory(Trajectory* traj);

// Reads the actual trajectory position information

void read_trajectory(Trajectory* traj);

struct TrajAsyncDefaultFunctor {
	void operator()() const {};
};

template<typename OnFinishFunctor = TrajAsyncDefaultFunctor>
void read_trajectory_async(Trajectory* traj, OnFinishFunctor on_finish = TrajAsyncDefaultFunctor()) {
	ASSERT(traj);
	if (traj->path_to_file) {
		std::thread([traj, on_finish]() {
			read_trajectory(traj);
			on_finish();
		}).detach();
	}
}

TrajectoryFrame allocate_trajectory_frame(int num_atoms);
void free_trajectory_frame(TrajectoryFrame* frame);

void copy_trajectory_frame(TrajectoryFrame* dst, const Trajectory& src_traj, int frame_index);
void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index);

TrajectoryFrame get_trajectory_frame(const Trajectory& traj, int frame_index);
Array<vec3> get_trajectory_positions(const Trajectory& traj, int frame_index);