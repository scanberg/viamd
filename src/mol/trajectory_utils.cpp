#include "trajectory_utils.h"
#include <core/string_utils.h>
#include <core/log.h>

#include <stdio.h>
#include <xdrfile_xtc.h>
#include <string>
#include <fstream>

// @TODO: Remove dependency of string

bool load_and_allocate_trajectory(Trajectory* traj, CString path) {
	ASSERT(traj);
	free_trajectory(traj);

	// @TODO: Remove this shit
    std::string url(path);
    auto dot_pos = url.find_last_of(".");
    std::string file_without_ext = url.substr(0, dot_pos);
    std::string cache_file = file_without_ext + ".cache";

    XDRFILE* file_handle = xdrfile_open(path, "r");
    if (!file_handle) {
        return false;
    }

    int num_atoms = 0;
    int num_frames = 0;
    int64* offsets = nullptr;
    if (read_xtc_natoms(path, &num_atoms) != exdrOK) {
        return false;
    }

    std::ifstream offset_cache_stream(cache_file, std::ios::binary);
    if (offset_cache_stream) {
        offset_cache_stream.seekg(0, std::ios::end);
        int64 byte_size = offset_cache_stream.tellg();
        offset_cache_stream.seekg(0, std::ios::beg);
        offsets = (int64*)malloc(byte_size);
        offset_cache_stream.read((char*)offsets, byte_size);
        num_frames = (int)(byte_size / sizeof(int64));
    } else {
        if (read_xtc_frame_offsets(path, &num_frames, &offsets) != exdrOK) {
            return false;
        }
        std::ofstream write_offset_cache_stream(cache_file, std::ios::binary);
        if (write_offset_cache_stream) {
            write_offset_cache_stream.write((char*)offsets, num_frames * sizeof(int64));
        }
    }

    if (!offsets) {
        return false;
    }

	traj->num_atoms = num_atoms;
	traj->num_frames = 0;
	traj->total_simulation_time = 0;
	traj->simulation_type = Trajectory::NVT;
	traj->path_to_file = allocate_string(path);
	traj->file_handle = file_handle;

	traj->frame_offsets = { offsets, num_frames };
	traj->position_data = { (vec3*)MALLOC(num_frames * num_atoms * sizeof(vec3)), num_frames * num_atoms };
	traj->frame_buffer = { (TrajectoryFrame*)MALLOC(num_frames * sizeof(TrajectoryFrame)), num_frames };

    //traj->is_loading = false;
    //traj->signal_stop = false;

    return true;
}

bool read_trajectory_data(Trajectory* traj) {
	ASSERT(traj);
	auto num_frames = traj->frame_offsets.count;
	XDRFILE* file = xdrfile_open(traj->path_to_file, "r");
	if (!file) {
		LOG_ERROR("Could not open file %s\n", traj->path_to_file.beg());
		return false;
	}
	//traj->is_loading = true;
	for (int i = 0; i < num_frames; i++) {
		//if (traj->signal_stop) {
		//	traj->is_loading = false;
		//	return false;
		//}
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
	//traj->is_loading = false;
	return true;
}

bool read_next_trajectory_frame(Trajectory* traj) {
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
	read_xtc((XDRFILE*)traj->file_handle, traj->num_atoms, &step, &frame->time, (float(*)[3]) & frame->box, (float(*)[3])pos_data, &precision);
	for (int j = 0; j < traj->num_atoms; j++) {
		pos_data[j] *= 10.f;
	}
	frame->box *= 10.f;
	traj->num_frames++;

	return true;
}

bool all_trajectory_frames_read(Trajectory* traj) {
	ASSERT(traj);
	return (traj->num_frames == (int32)traj->frame_offsets.count);
}

bool close_file_handle(Trajectory* traj) {
	ASSERT(traj);
	if (traj->file_handle) {
		xdrfile_close((XDRFILE*)traj->file_handle);
		traj->file_handle = nullptr;
		return true;
	}
	return false;
}

void copy_trajectory_frame(TrajectoryFrame* frame, const Trajectory& traj, int frame_index) {
	ASSERT(frame);
	ASSERT(frame_index < traj.num_frames);
	memcpy(frame, &traj.frame_buffer[frame_index], sizeof(TrajectoryFrame));
}

void copy_trajectory_positions(Array<vec3> dst_array, const Trajectory& traj, int frame_index) {
    ASSERT(dst_array);
    ASSERT(dst_array.count >= traj.num_atoms);
    ASSERT(frame_index < traj.num_frames);
    memcpy(dst_array.data, traj.frame_buffer.data[frame_index].atom_positions.data, traj.num_atoms * sizeof(vec3));
}

void read_trajectory_box_vectors(vec3 box_vectors[3], const Trajectory& traj, int frame_index) {
	(void)box_vectors;
	(void)traj;
	(void)frame_index;
}

TrajectoryFrame get_trajectory_frame(const Trajectory& traj, int frame_index) {
    ASSERT(frame_index < traj.num_frames);
    return traj.frame_buffer.data[frame_index];
}

Array<vec3> get_trajectory_positions(const Trajectory& traj, int frame_index) {
    ASSERT(frame_index < traj.num_frames);
    return traj.frame_buffer.data[frame_index].atom_positions;
}