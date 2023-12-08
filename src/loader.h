#pragma once

#include <core/md_str.h>

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_loader_i;
struct md_trajectory_i;
struct md_trajectory_loader_i;
struct md_bitfield_t;

namespace load {
    int64_t      supported_extension_count();
    const str_t* supported_extensions();

namespace mol {
    md_molecule_loader_i* get_loader_from_ext(str_t ext);
}

namespace traj {
    md_trajectory_loader_i* get_loader_from_ext(str_t ext);

    // loader is optional, the default loader (determined from file extension will be used) if NULL
    md_trajectory_i* open_file(str_t filename, md_trajectory_loader_i* loader, const md_molecule_t* mol, md_allocator_i* alloc, bool deperiodize_on_load);
    bool close(md_trajectory_i* traj);

    bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask);
    bool set_deperiodize(md_trajectory_i* traj, bool deperiodize);

    bool clear_cache(md_trajectory_i* traj);
    int64_t num_cache_frames(md_trajectory_i* traj);
}

}  // namespace load
