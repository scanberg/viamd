#pragma once

#include <core/types.h>
#include <core/array.h>

struct DensityVolume {
	typedef float VoxelDataType;
    ivec3 dim;
	Range voxel_range;
    Array<VoxelDataType> voxel_data;
};

void init_volume(DensityVolume* vol, ivec3 dim);
void free_volume(DensityVolume* vol);
void clear_volume(DensityVolume* vol);