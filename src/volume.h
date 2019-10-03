#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>

struct Volume {
    using VoxelDataType = float;
    ivec3 dim = {0, 0, 0};
    ArrayView<VoxelDataType> voxel_data{};
    Range<float32> voxel_range = {0, 0};
};

void init_volume(Volume* vol, ivec3 dim);
void free_volume(Volume* vol);
void clear_volume(Volume* vol);
