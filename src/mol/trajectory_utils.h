#pragma once

#include <mol/trajectory.h>
#include <thread>

// Reads the header info of a trajectory and allocates space for it
bool load_and_allocate_trajectory(Trajectory* traj, CString path);

// Reads the actual trajectory position information... Necessary?
bool read_trajectory_data(Trajectory* traj);

struct TrajAsyncDefaultFunctor {
	void operator()() const {};
};

template<typename OnFinishFunctor = TrajAsyncDefaultFunctor>
void read_trajectory_async(Trajectory* traj, OnFinishFunctor on_finish = TrajAsyncDefaultFunctor()) {
	ASSERT(traj);
	if (traj->path_to_file) {
		std::thread([traj, on_finish]() {
			if (read_trajectory_data(traj)) {
				on_finish();
			}
		}).detach();
	}
}

void copy_trajectory_frame(TrajectoryFrame* dst, const Trajectory& src_traj, int frame_index);
void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index);
void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index);

TrajectoryFrame get_trajectory_frame(const Trajectory& traj, int frame_index);
Array<vec3> get_trajectory_positions(const Trajectory& traj, int frame_index);