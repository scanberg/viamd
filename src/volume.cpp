#include "volume.h"

void init_volume(Volume* vol, ivec3 dim) {
    ASSERT(vol);
    if (vol->voxel_data) {
        free_array(&vol->voxel_data);
    }
    int32 count = dim.x * dim.y * dim.z;
    ASSERT(count > 0);
    vol->voxel_data = allocate_array<Volume::VoxelDataType>(count);
    vol->dim = dim;
    vol->voxel_range = {0, 0};
}

void free_volume(Volume* vol) {
    ASSERT(vol);
    if (vol->voxel_data) {
        free_array(&vol->voxel_data);
    }
    vol->dim = {0, 0, 0};
    vol->voxel_range = {0, 0};
}

void clear_volume(Volume* vol) {
    ASSERT(vol);
    memset(vol->voxel_data.ptr, 0, vol->voxel_data.count * sizeof(Volume::VoxelDataType));
    vol->voxel_range = {0, 0};
}
