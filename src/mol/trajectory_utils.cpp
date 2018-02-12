#include "trajectory_utils.h"
#include <core/string_utils.h>

#include <stdio.h>
#include <xdrfile_xtc.h>


Trajectory read_trajectory(const char* file) {
	CString url(file);
	char file_buffer[256];
	char dir_buffer[256];
	String file_str(file_buffer, 256);
	String dir_str(dir_buffer, 256);

	copy(dir_str, get_directory(url));
	copy(file_str, get_file_without_extension(url));

	char trajectory_cache[256];
	String cache_str(trajectory_cache, 256);

	XDRFILE* file_handle = xdrfile_open(file, "r");
	int num_atoms, step;
	float time;
	read_xtc_header(file_handle, &num_atoms, &step, &time);

	return {};
}

void destroy_trajectory(Trajectory* traj) {
	ASSERT(traj);
	if (traj->file_handle) xdrfile_close((XDRFILE*)traj->file_handle);
}

TrajectoryFrame copy_trajectory_frame(Trajectory traj, int frame_index, Allocator& alloc) {
	return {};
}

void read_trajectory_positions(Array<vec3> atom_positions, Trajectory traj, int frame_index) {

}

void read_trajectory_box_vectors(vec3 box_vectors[3], Trajectory traj, int frame_index) {

}