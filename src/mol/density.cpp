#include "density.h"

void init_volume(DensityVolume* vol, ivec3 dim) {
	ASSERT(vol);
	if (vol->voxel_data) {
		free_array(&vol->voxel_data);
	}
	int32 count = dim.x * dim.y * dim.z;
	ASSERT(count > 0);
	vol->voxel_data = allocate_array<DensityVolume::VoxelDataType>(count);
	vol->dim = dim;
	vol->voxel_range = { 0,0 };
}

void free_volume(DensityVolume* vol) {
	ASSERT(vol);
	if (vol->voxel_data) {
		free_array(&vol->voxel_data);
	}
	vol->dim = { 0 };
	vol->voxel_range = { 0,0 };
}

void clear_volume(DensityVolume* vol) {
	ASSERT(vol);
	memset(vol->voxel_data.data, 0, vol->voxel_data.count * sizeof(DensityVolume::VoxelDataType));
}