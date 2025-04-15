#pragma once

#include <core/md_str.h>

struct md_allocator_i;
struct md_molecule_t;
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

typedef struct FrameLock;

typedef struct FrameData {
    // Header
    size_t  num_atoms;
    int64_t index;
    double  timestamp;
    md_unit_cell_t unit_cell;

    // Coordinates
    float* x;
    float* y;
    float* z;
};

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
    typedef uint64_t Handle;
    md_trajectory_loader_i* loader_from_ext(str_t ext);

    md_trajectory_i* open_file(str_t filename, md_trajectory_loader_i* loader, const md_molecule_t* mol, md_allocator_i* alloc, LoadTrajectoryFlags flags = LoadTrajectoryFlag_None);
    bool close(md_trajectory_i* traj);

    bool has_recenter_target(md_trajectory_i* traj);
    bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask);

    bool clear_cache(md_trajectory_i* traj);
    size_t num_cache_frames(md_trajectory_i* traj);


    // New interface

    size_t set_cache_memory_budget(size_t megabytes);

    Handle open(str_t filename);
    void   close(Handle traj);

    md_bitfield_t* get_recenter_target(Handle traj);
    bool set_recenter_target(Handle traj, const md_bitfield_t* atom_mask);

    bool evict_from_cache(Handle traj);

    // This operation pins a frame in the cache, which will ensures that it cannot be evicted by other threads.
    // If the frames data is not present in the cache, it will load the data first.
    bool pin_frame    (Handle traj, FrameLock* lock, FrameData* out_frame_data, size_t frame_idx);

    // This function attempts to pin a frame in place, but will only succeeed if it is already in the cache.
    bool try_pin_frame(Handle traj, FrameLock* lock, FrameData* out_frame_data, size_t frame_idx);

    void unpin_frame  (Handle traj, FrameLock* lock);
}

}  // namespace load
