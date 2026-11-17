#include "volume.h"

void init_volume(Volume* vol, ivec3 dim) {
    MD_ASSERT(vol);
    if (vol->voxel_data) {
        free_array(&vol->voxel_data);
    }
    i32 count = dim.x * dim.y * dim.z;
    MD_ASSERT(count > 0);
    vol->voxel_data = allocate_array<Volume::VoxelDataType>(count);
    vol->dim = dim;
    vol->voxel_range = {0, 0};
}

void free_volume(Volume* vol) {
    MD_ASSERT(vol);
    if (vol->voxel_data) {
        free_array(&vol->voxel_data);
    }
    vol->dim = {0, 0, 0};
    vol->voxel_range = {0, 0};
}

void clear_volume(Volume* vol) {
    MD_ASSERT(vol);
    memset(vol->voxel_data.ptr, 0, vol->voxel_data.count * sizeof(Volume::VoxelDataType));
    vol->voxel_range = {0, 0};
}
