#include "trajectory_utils.h"'
#include <core/string_utils.h>

#include <stdio.h>
#include <xdrfile_xtc.h>
#include <string>

// @TODO: Remove dependency of string

Trajectory read_and_allocate_trajectory(const char* path, Allocator& alloc) {
	std::string url(path);

	auto pos = url.find_last_of("\\/");
	std::string dir = url.substr(0, pos);
	std::string file = url.substr(pos + 1);
	
	printf("'%s' '%s'", dir.c_str(), file.c_str());

	XDRFILE* file_handle = xdrfile_open(path, "r");
	int num_atoms, step;
	float time;
	read_xtc_header(file_handle, &num_atoms, &step, &time);

	return { num_atoms, step, time, Trajectory::NVT, file_handle, {} };
}

void free_trajectory(Trajectory* traj) {
	ASSERT(traj);
	if (traj->file_handle) xdrfile_close((XDRFILE*)traj->file_handle);
	traj->file_handle = nullptr;
}

TrajectoryFrame copy_trajectory_frame(Trajectory traj, int frame_index, Allocator& alloc) {
	return {};
}

void read_trajectory_positions(Array<vec3> atom_positions, Trajectory traj, int frame_index) {

}

void read_trajectory_box_vectors(vec3 box_vectors[3], Trajectory traj, int frame_index) {

}