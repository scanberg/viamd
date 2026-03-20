#pragma once

#include <stdint.h>
#include <core/md_str.h>

struct md_allocator_i;
struct md_system_t;
struct md_system_loader_i;
struct md_trajectory_i;
struct md_bitfield_t;

typedef bool (*md_traj_attach_fn)(md_system_t* sys, str_t filename, uint32_t flags);

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
        // This represents the create callbacks and arguments needed to load a molecule or trajectory from a file.
    struct LoaderState {		
		md_system_loader_i*     sys_loader = 0;
		md_traj_attach_fn       traj_attach = 0;
        const void*             sys_loader_arg = 0;
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
    md_system_loader_i* loader_from_ext(str_t ext);
}

namespace traj {
    md_traj_attach_fn attach_fn_from_ext(str_t ext);

    bool attach_file(md_system_t* sys, str_t filename, md_traj_attach_fn fn, LoadTrajectoryFlags flags = LoadTrajectoryFlag_None);

    bool has_recenter_target(md_trajectory_i* traj);
    bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask);
}

}  // namespace load
