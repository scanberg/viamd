#include "trajectory_utils.h"
#include <core/string_utils.h>
#include <core/log.h>

#include <stdio.h>
#include <xdrfile_xtc.h>

bool load_and_allocate_trajectory(MoleculeTrajectory* traj, CString path) {
    ASSERT(traj);
    free_trajectory(traj);

    CString directory = get_directory(path);
    CString file = get_file_without_extension(path);
    StringBuffer<512> cache_file = directory;
    cache_file += "/";
    cache_file += file;
    cache_file += ".cache";

    XDRFILE* file_handle = xdrfile_open(path, "r");
    if (!file_handle) {
        return false;
    }

    int32 num_atoms = 0;
    int32 num_frames = 0;
    int64* offsets = nullptr;
    if (read_xtc_natoms(path, &num_atoms) != exdrOK) {
        return false;
    }

    FILE* offset_cache_handle = fopen(cache_file, "rb");
    if (offset_cache_handle) {
        fseek(offset_cache_handle, 0, SEEK_END);
        int64 byte_size = ftell(offset_cache_handle);
        offsets = (int64*)malloc(byte_size);
        num_frames = (int32)(byte_size / sizeof(int64));
        fread(offsets, sizeof(int64), num_frames, offset_cache_handle);
        fclose(offset_cache_handle);
    } else {
        if (read_xtc_frame_offsets(path, &num_frames, &offsets) != exdrOK) {
            return false;
        }
        FILE* write_cache_handle = fopen(cache_file, "wb");
        if (write_cache_handle) {
            fwrite(offsets, sizeof(int64), num_frames, write_cache_handle);
            fclose(write_cache_handle);
        }
    }

    if (!offsets) {
        return false;
    }

    traj->num_atoms = num_atoms;
    traj->num_frames = 0;
    traj->total_simulation_time = 0;
    traj->simulation_type = MoleculeTrajectory::NVT;
    traj->path_to_file = allocate_string(path);
    traj->file_handle = file_handle;

    traj->frame_offsets = {offsets, num_frames};
    traj->position_data = {(vec3*)MALLOC(num_frames * num_atoms * sizeof(vec3)), num_frames * num_atoms};
    traj->frame_buffer = {(TrajectoryFrame*)MALLOC(num_frames * sizeof(TrajectoryFrame)), num_frames};

    return true;
}

bool read_trajectory_data(MoleculeTrajectory* traj) {
    ASSERT(traj);
    auto num_frames = traj->frame_offsets.count;
    XDRFILE* file = xdrfile_open(traj->path_to_file, "r");
    if (!file) {
        LOG_ERROR("Could not open file %s\n", traj->path_to_file.beg());
        return false;
    }

    for (int i = 0; i < num_frames; i++) {
        vec3* pos_data = traj->position_data.data + (i * traj->num_atoms);
        TrajectoryFrame* frame = traj->frame_buffer.data + i;
        frame->atom_positions.data = pos_data;
        frame->atom_positions.count = traj->num_atoms;
        frame->index = i;
        int step;
        float precision;
        read_xtc(file, traj->num_atoms, &step, &frame->time, (float(*)[3]) & frame->box, (float(*)[3])pos_data, &precision);
        for (int j = 0; j < traj->num_atoms; j++) {
            pos_data[j] *= 10.f;
        }
        frame->box *= 10.f;
        traj->num_frames++;
    }

    return true;
}

bool read_next_trajectory_frame(MoleculeTrajectory* traj) {
    ASSERT(traj);
    if (!traj->file_handle) return false;
    auto num_frames = traj->frame_offsets.count;
    if (traj->num_frames == num_frames) return false;

    // Next index to be loaded
    int i = traj->num_frames;
    vec3* pos_data = traj->position_data.data + (i * traj->num_atoms);
    TrajectoryFrame* frame = traj->frame_buffer.data + i;
    frame->atom_positions.data = pos_data;
    frame->atom_positions.count = traj->num_atoms;
    frame->index = i;
    int step;
    float precision;
    float matrix[3][3];
    float* pos_buf = (float*)TMP_MALLOC(traj->num_atoms * 3 * sizeof(float));

    read_xtc((XDRFILE*)traj->file_handle, traj->num_atoms, &step, &frame->time, matrix, (float(*)[3])pos_buf, &precision);
    frame->box = mat3(matrix[0][0], matrix[0][1], matrix[0][2], matrix[1][0], matrix[1][1], matrix[1][2], matrix[2][0], matrix[2][1], matrix[2][2]);

    for (int j = 0; j < traj->num_atoms; j++) {
        pos_data[j] = vec3(10.f) * vec3(pos_buf[j * 3 + 0], pos_buf[j * 3 + 1], pos_buf[j * 3 + 2]);
    }
    frame->box *= 10.f;
    traj->num_frames++;

    TMP_FREE(pos_buf);

    return true;
}

bool all_trajectory_frames_read(MoleculeTrajectory* traj) {
    ASSERT(traj);
    return (traj->num_frames == (int32)traj->frame_offsets.count);
}

bool close_file_handle(MoleculeTrajectory* traj) {
    ASSERT(traj);
    if (traj->file_handle) {
        xdrfile_close((XDRFILE*)traj->file_handle);
        traj->file_handle = nullptr;
        return true;
    }
    return false;
}

void copy_trajectory_frame(TrajectoryFrame* frame, const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(frame);
    ASSERT(frame_index < traj.num_frames);
    memcpy(frame, &traj.frame_buffer[frame_index], sizeof(TrajectoryFrame));
}

void copy_trajectory_positions(Array<vec3> dst_array, const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(dst_array);
    ASSERT(dst_array.count >= traj.num_atoms);
    ASSERT(frame_index < traj.num_frames);
    constexpr auto size = sizeof(vec3);
    memcpy(dst_array.data, traj.frame_buffer.data[frame_index].atom_positions.data, traj.num_atoms * sizeof(vec3));
}

void read_trajectory_box_vectors(vec3 box_vectors[3], const MoleculeTrajectory& traj, int frame_index) {
    (void)box_vectors;
    (void)traj;
    (void)frame_index;
}
