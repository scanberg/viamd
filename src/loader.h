#pragma once

#include <core/md_str.h>

struct md_allocator_i;
struct md_system_t;
struct md_molecule_loader_i;
struct md_trajectory_i;
struct md_trajectory_loader_i;
struct md_bitfield_t;

// @NOTE(Robin): This API is currently a mess.

enum LoaderStateFlag_ {
    LoaderStateFlag_None = 0,
    LoaderStateFlag_RequiresDialogue = 1,
};

enum LoadTrajectoryFlag_ {
    LoadTrajectoryFlag_None = 0,
    LoadTrajectoryFlag_DisableCacheWrite = 1,
};

typedef uint32_t LoaderStateFlags;
typedef uint32_t LoadTrajectoryFlags;

namespace load {
    // This represents a loader state with arguments to load a molecule or trajectory from a file
    struct LoaderState {		
		md_molecule_loader_i*   mol_loader = 0;
		md_trajectory_loader_i* traj_loader = 0;
        const void*             mol_loader_arg = 0;
        LoaderStateFlags 		flags = LoaderStateFlag_None;

        size_t data_size = 0;
        void* data_ptr = 0;
	};

    bool init_loader_state(LoaderState* state, str_t filepath, md_allocator_i* alloc);
    void free_loader_state(LoaderState* state, md_allocator_i* alloc);

    size_t       loader_count();
    const str_t* loader_names();
    const str_t* loader_extensions();

namespace mol {
    md_molecule_loader_i* loader_from_ext(str_t ext);
}

namespace traj {
    md_trajectory_loader_i* loader_from_ext(str_t ext);

    md_trajectory_i* open_file(str_t filename, md_trajectory_loader_i* loader, const md_system_t* mol, md_allocator_i* alloc, LoadTrajectoryFlags flags = LoadTrajectoryFlag_None);
    bool close(md_trajectory_i* traj);

	// Get the internal trajectory, this can be used to access custom loader functionality
	// This is the internal trajectory without any form of caching or recentering applied
	md_trajectory_i* get_raw_trajectory(md_trajectory_i* traj);

    bool has_recenter_target(md_trajectory_i* traj);
    bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask);

    bool clear_cache(md_trajectory_i* traj);
    size_t num_cache_frames(md_trajectory_i* traj);
}

}  // namespace load
