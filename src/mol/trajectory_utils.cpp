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

    int num_atoms;
    int num_frames;
    int64* offsets;
    read_xtc_natoms(path, &num_atoms);
    read_xtc_frame_offsets(path, &num_frames, &offsets);

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
        float precision;
        read_xtc(file_handle, num_atoms, &frame->index, &frame->time, (float(*)[3])&frame->box, (float(*)[3])pos_data, &precision);
    }

    return traj;
}

void free_trajectory(Trajectory* traj) {
    ASSERT(traj);
    if (traj->file_handle) xdrfile_close((XDRFILE*)traj->file_handle);
    traj->file_handle = nullptr;
}

TrajectoryFrame copy_trajectory_frame(Trajectory traj, int frame_index, Allocator& alloc) { return {}; }

void read_trajectory_positions(Array<vec3> atom_positions, Trajectory traj, int frame_index) {}

void read_trajectory_box_vectors(vec3 box_vectors[3], Trajectory traj, int frame_index) {}