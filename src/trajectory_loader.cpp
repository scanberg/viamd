#include "trajectory_loader.h"
#include <core/common.h>
#include <core/log.h>
#include <core/file.h>
#include <core/math_utils.h>
#include <core/string_types.h>
#include <mol/xtc_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/trajectory_utils.h>

#include <task_system.h>

struct MoleculeLoader {
    CStringView id;
    bool (*read_molecule)(MoleculeStructure* mol, CStringView filename);
};

static constexpr MoleculeLoader mol_loaders[]{
    {"pdb", pdb::load_molecule_from_file},
    {"gro", gro::load_molecule_from_file},
};

struct TrajectoryLoader {
    CStringView id;
    bool (*read_num_frames)(i32* num_frames, CStringView filename);
    bool (*read_simulation_box)(mat3* sim_box, CStringView filename);
    bool (*read_frame_bytes)(FrameBytes* frame_bytes, CStringView filename);
    bool (*extract_frame)(TrajectoryFrame* frame, i32 num_atoms, Array<u8> data);
    // void (*post_init)(MoleculeTrajectory* traj, CStringView filename) = 0;
};

static constexpr TrajectoryLoader traj_loaders[]{
    {"xtc", xtc::read_trajectory_num_frames, NULL, xtc::read_trajectory_frame_bytes, xtc::decompress_trajectory_frame},
    {"pdb", pdb::read_trajectory_num_frames, pdb::read_trajectory_simulation_box, pdb::read_trajectory_frame_bytes, pdb::extract_trajectory_frame},
};

namespace load {

namespace mol {
static const MoleculeLoader* get_loader(CStringView extension) {
    for (const auto& loader : mol_loaders) {
        if (compare_ignore_case(extension, loader.id)) {
            return &loader;
        }
    }
    return NULL;
}

bool is_extension_supported(CStringView filename) { return get_loader(get_file_extension(filename)) != NULL; }

bool load_molecule(MoleculeStructure* mol, CStringView filename) {
    const MoleculeLoader* loader = get_loader(get_file_extension(filename));
    if (!loader) {
        LOG_ERROR("Could not find loader for file '%.*s'", filename.length(), filename.beg());
        return false;
    }
    return loader->read_molecule(mol, filename);
}

}  // namespace mol

namespace traj {
static const TrajectoryLoader* get_loader(CStringView extension) {
    for (const auto& loader : traj_loaders) {
        if (compare_ignore_case(extension, loader.id)) {
            return &loader;
        }
    }
    return NULL;
}

bool is_extension_supported(CStringView filename) { return get_loader(get_file_extension(filename)) != NULL; }

bool read_num_frames(i32* num_frames, CStringView filename) {
    const TrajectoryLoader* loader = get_loader(get_file_extension(filename));
    if (!loader) {
        LOG_ERROR("Could not find loader for file '%.*s'", filename.length(), filename.beg());
        return false;
    }
    return loader->read_num_frames(num_frames, filename);
}

bool load_trajectory(MoleculeTrajectory* traj, i32 num_atoms, CStringView filename) {
    const TrajectoryLoader* loader = get_loader(get_file_extension(filename));
    if (!loader) {
        LOG_ERROR("Could not find loader for file '%.*s'", filename.length(), filename.beg());
        return false;
    }

    if (num_atoms <= 0) {
        LOG_ERROR("A trajectory needs a positive number of atoms supplied to be matched against molecule structure");
        return false;
    }

    FILE* file = fopen(filename, "rb");
    defer { fclose(file); };
    if (!file) {
        LOG_ERROR("Could not open file '%.*s'", filename.length(), filename.beg());
        return false;
    }

    i32 num_frames = 0;
    if (!loader->read_num_frames(&num_frames, filename)) {
        LOG_ERROR("Could not read number of frames");
        return false;
    }

    mat3 sim_box = {};
    if (loader->read_simulation_box && !loader->read_simulation_box(&sim_box, filename)) {
        LOG_ERROR("Could not read simulation box");
        return false;
    }
    ASSERT(num_frames > 0);
    init_trajectory(traj, num_atoms, num_frames, 1.0f, sim_box);

    FrameBytes* frame_bytes = (FrameBytes*)TMP_MALLOC(num_frames * sizeof(frame_bytes));
    defer { TMP_FREE(frame_bytes); };

    if (!loader->read_frame_bytes(frame_bytes, filename)) {
        LOG_ERROR("Could not read frame bytes");
        return false;
    }

    constexpr u64 read_buffer_size = MEGABYTES(128);
    void* mem = TMP_MALLOC(read_buffer_size);
    defer { TMP_FREE(mem); };

    u64 max_frame_bytes = 0;
    for (i32 i = 0; i < num_frames; i++) {
        max_frame_bytes = math::max(max_frame_bytes, frame_bytes[i].extent);
    }

    const i32 batch_size = (i32)(read_buffer_size / max_frame_bytes);

    for (i32 i = 0; i < num_frames; i += batch_size) {
        int batch_beg = i;
        int batch_end = math::min(num_frames, i + batch_size);
        int batch_ext = batch_end - batch_beg;

        const i64 batch_bytes = frame_bytes[batch_end - 1].offset + frame_bytes[batch_end - 1].extent - frame_bytes[batch_beg].offset;
        fseeki64(file, frame_bytes[batch_beg].offset, SEEK_SET);
        fread(mem, 1, batch_bytes, file);

        task_system::ID id = task_system::enqueue_pool(
            "Extracting Trajectory Frames", batch_ext,
            [traj, loader, num_atoms, frame_bytes, mem, batch_offset = batch_beg](task_system::TaskSetRange range) {
                for (u32 i = batch_offset + range.beg; i < batch_offset + range.end; i++) {
                    const i64 mem_offset = frame_bytes[i].offset - frame_bytes[batch_offset].offset;
                    auto& frame = get_trajectory_frame(*traj, i);
                    loader->extract_frame(&frame, num_atoms, {(u8*)mem + mem_offset, (i64)frame_bytes[i].extent});
                }
            });
        task_system::wait_for_task(id);
        //task_system::clear_completed_tasks();  // @TODO: REMOVE THIS WHEN CALLBACKS ARE AVAILABLE SO TASKS CAN BE FREED AND PUT BACK INTO THE QUEUE
    }
    return true;
}

}  // namespace traj

}  // namespace load