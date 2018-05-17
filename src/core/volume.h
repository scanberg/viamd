#pragma once

#include <core/types.h>
#include <core/array.h>

struct Volume {
    typedef float VoxelDataType;
    ivec3 dim = {0, 0, 0};
    Range voxel_range = {0, 0};
    Array<VoxelDataType> voxel_data{};
    vec3 min_box = {0, 0, 0};
    vec3 max_box = {0, 0, 0};
};

void init_volume(Volume* vol, ivec3 dim);
void free_volume(Volume* vol);
void clear_volume(Volume* vol);
