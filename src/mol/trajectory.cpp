#include "trajectory.h"

bool init_trajectory(Trajectory* traj, int32 num_atoms, int32 num_frames) {
	ASSERT(traj);

	traj->num_atoms = num_atoms;
	traj->num_frames = num_frames;
	traj->total_simulation_time = 0;
	traj->simulation_type = Trajectory::NVT;
	traj->path_to_file = {};
	traj->file_handle = nullptr;

	traj->frame_offsets = {};
	traj->position_data = { (vec3*)CALLOC(num_frames * num_atoms, sizeof(vec3)), num_frames * num_atoms };
	traj->frame_buffer = { (TrajectoryFrame*)CALLOC(num_frames, sizeof(TrajectoryFrame)), num_frames };

	//traj->is_loading = false;
	//traj->signal_stop = false;

	return true;
}

void free_trajectory(Trajectory* traj) {
	ASSERT(traj);

	// @TODO: Fix this using condition variables (platform specific???)
	//if (traj->is_loading) {
	//	traj->signal_stop = true;
	//	while (traj->is_loading) {
	//		// SLEEP?
	//	}
	//}

	free_string(&traj->path_to_file);
	if (traj->frame_offsets.data) FREE(traj->frame_offsets.data);
	if (traj->position_data.data) FREE(traj->position_data.data);
	if (traj->frame_buffer.data) FREE(traj->frame_buffer.data);

	*traj = {};
}
