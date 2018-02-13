#pragma once

#include <core/types.h>
#include <core/array.h>
#include <core/common.h>

constexpr int MAX_TRAJECTORY_FRAME_BUFFER_SIZE = GIGABYTES(2);

struct TrajectoryFrame {
    int frame_index;
    Array<vec3> atom_positions;
    vec3 box_vectors[3];
};

struct Trajectory {
	enum Type { NVT, NPT };

	int32	num_atoms = 0;
	int32	num_frames = 0;
	float32 total_simulation_time = 0;
	Type	simulation_type = NVT;

    // This is optionally used since the data may fit into memory.
    void* file_handle = nullptr;

    // @NOTE: The frame_buffer may not contain all frames in trajectory.
    // If the trajectory is large, frame_buffer will be used as a cache towards the trajectory streamed from disk.
	DynamicArray<TrajectoryFrame> frame_buffer{};
};
