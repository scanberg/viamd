#pragma once

#include <mol/trajectory.h>
#include <thread>

// Reads the header info of a trajectory and allocates space for it
Trajectory* allocate_trajectory(const char* file);
void free_trajectory(Trajectory* traj);

// Reads the actual trajectory position information

void read_trajectory(Trajectory* traj);

struct TrajOnFinishFunctor {
	void operator()() const {};
};
template<typename Functor = TrajOnFinishFunctor>
void read_trajectory_async(Trajectory* traj, Functor on_finish = TrajOnFinishFunctor()) {
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