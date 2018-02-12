#pragma once

#include <core/types.h>
#include <core/array.h>
#include <core/common.h>

constexpr int MAX_TRAJECTORY_FRAME_BUFFER_SIZE = GIGABYTES(2);

enum class SimulationType { NVT, NPT };

struct TrajectoryFrame {
    int frame_index;
    Array<vec3> atom_positions;
    vec3 box_vectors[3];
};

struct Trajectory {
    SimulationType type = SimulationType::NVT;

    // This is optionally used since the data may fit into memory.
    void* file_handle = nullptr;

    // @NOTE: The frame_buffer may not contain all frames in trajectory.
    // If the trajectory is large, frame_buffer will be used as a cache towards the trajectory streamed from disk.
	Array<TrajectoryFrame> frame_buffer{};
};
