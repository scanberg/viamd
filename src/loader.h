#pragma once

#include <stdint.h>
#include <core/md_str.h>

struct md_allocator_i;
struct md_molecule_t;
struct md_trajectory_i;

namespace load {
namespace mol {
    bool is_extension_supported(str_t filename);

    bool load_string_pdb(md_molecule_t* mol, str_t string, md_allocator_i* alloc);
    bool load_string_gro(md_molecule_t* mol, str_t string, md_allocator_i* alloc);
    bool load_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc);
    bool free(md_molecule_t* mol);
}

namespace traj {
    bool is_extension_supported(str_t filename);

    bool open_file(md_trajectory_i* traj, str_t filename, const md_molecule_t* mol, md_allocator_i* alloc);
    bool close(md_trajectory_i* traj);

    // Loads data from trajectory frames.
    // If the frame does not reside in the cache, then the data is loaded from disk and
    bool load_trajectory_frame_box(const md_trajectory_i* traj, float* box[3][3], int64_t frame_idx);
    bool load_trajectory_frame_coords(const md_trajectory_i* traj, float* x, float* y, float* z, int64_t atom_count, int64_t frame_idx);
}

}  // namespace load