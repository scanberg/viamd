#pragma once

#include <core/types.h>
#include <core/array.h>
#include <core/string_utils.h>
#include <core/common.h>

struct TrajectoryFrame {
    int index;
    float time;
    mat3 box;
    Array<vec3> atom_positions;
};

struct Trajectory {
	enum Type { NVT, NPT };

	int32	num_atoms = 0;
	int32	num_frames = 0;
	float32 total_simulation_time = 0;
	Type	simulation_type = NVT;

	CString path_to_file = "";

    // @NOTE: The frame_buffer may not contain all frames in trajectory.
    // If the trajectory is large, frame_buffer will be used as a cache towards the trajectory streamed from disk.
	Array<TrajectoryFrame> frame_buffer{};

    // This is the position data of the trajectories
	Array<vec3> position_data{};

    // These are the offsets for each frame within the compressed blob of XTC data
	Array<int64> frame_offsets{};

	// This is for synchronization
	volatile bool is_loading = false;
};
