#pragma once

#include <core/types.h>
#include <core/common.h>

struct TrajectoryFileHandle;

enum class SimulationType { NVT, NVP; };

struct TrajectoryFrame {
    int frame_index;
    Array<vec3> atom_positions;
    vec3 box_vectors[3];
};

struct Trajectory {
    SimulationType type;
    TrajectoryFileHandle file;
    // Note that frame_buffer may not contain all frames in trajectory.
    // If the trajectory is large, frame_buffer will be used as a cache towards the trajectory streamed from disk.
    Array<TrajectoryFrame> frame_buffer;
};

constexpr int MaxTrajectoryFrameBufferSize = (GIGABYTE(2) / sizeof(TrajectoryFrame)) * sizeof(TrajectoryFrame);
