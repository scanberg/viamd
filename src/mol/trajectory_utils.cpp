#include "trajectory_utils.h"
#include <core/string_utils.h>

#include <stdio.h>
#include <xdrfile_xtc.h>
#include <string>

// @TODO: Remove dependency of string

Trajectory read_and_allocate_trajectory(const char* path, Allocator& alloc) {
    std::string url(path);

    auto dot_pos = url.find_last_of(".");
    std::string file_without_ext = url.substr(0, dot_pos);
    std::string cache_file = file_without_ext + ".cache";

    XDRFILE* file_handle = xdrfile_open(path, "r");
	if (!file_handle) {
		ASSERT(false, "Failed to open xdrfile");
	}

    int num_atoms;
    int num_frames;
    int64* offsets;
	if (read_xtc_natoms(path, &num_atoms) != exdrOK) {
		ASSERT(false, "Failed to extract num atoms");
	}
	if (read_xtc_frame_offsets(path, &num_frames, &offsets) != exdrOK) {
		ASSERT(false, "Failed to extract frames and offsets");
	}

    Trajectory traj{num_atoms, num_frames, 0.f, Trajectory::NVT, file_handle};
    traj.frame_offsets.resize(num_frames);
    memcpy(traj.frame_offsets.data, offsets, num_frames * sizeof(int64));
    free(offsets);

    printf("num_atoms: %i, num_frames: %i\n", num_atoms, num_frames);

    // @TODO: Only read in data if it fits into memory
    traj.position_data.resize(num_frames * num_atoms);
    traj.frame_buffer.resize(num_frames);

    for (int i = 0; i < num_frames; i++) {
        vec3* pos_data = traj.position_data.data + (i * num_atoms);
        TrajectoryFrame* frame = traj.frame_buffer.data + i;
        frame->atom_positions.data = pos_data;
        frame->atom_positions.count = num_atoms;
		frame->index = i;
		int step;
        float precision;
        read_xtc(file_handle, num_atoms, &step, &frame->time, (float(*)[3])&frame->box, (float(*)[3])pos_data, &precision);
		for (int j = 0; j < num_atoms; j++) {
			pos_data[j] *= 10.f;
		}
    }

    return traj;
}

void free_trajectory(Trajectory* traj) {
    ASSERT(traj);
    if (traj->file_handle) xdrfile_close((XDRFILE*)traj->file_handle);
    traj->file_handle = nullptr;
}

TrajectoryFrame copy_trajectory_frame(const Trajectory& traj, int frame_index, Allocator& alloc) { return {}; }

void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index) {
	ASSERT(dst_array);
	ASSERT(dst_array.count >= traj.num_atoms);
	ASSERT(frame_index < traj.num_frames);
	memcpy(dst_array.data, traj.frame_buffer.data[frame_index].atom_positions.data, traj.num_atoms * sizeof(vec3));
}

void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index) {}