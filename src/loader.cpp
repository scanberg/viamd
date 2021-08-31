#include "loader.h"
//#include <core/common.h>
//#include <core/log.h>
//#include <core/file.h>
#include <core/string_types.h>
#include <core/string_utils.h>
#include <core/sync.h>

#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_util.h>

#include <task_system.h>

#include <thread>

#include <core/lru_cache.h>

struct LoadedMolecule {
    uint64_t key;
    StringBuffer<8> extension;
    md_allocator_i* alloc;
};

struct FrameData {
    int64_t num_atoms;
    float* x;
    float* y;
    float* z;
    float box[3][3];
    double timestamp;
};

struct LoadedTrajectory {
    uint64_t key;
    StringBuffer<8> extension;
    const md_molecule_t* mol;
    md_allocator_i* alloc;
    FrameData frame_data[8];
    LRU_Cache_8<int64_t, int32_t> frame_cache;
};


static LoadedMolecule loaded_molecules[8] = {};
static int64_t num_loaded_molecules = 0;

static inline LoadedMolecule* find_loaded_molecule(uint64_t key) {
    for (uint64_t i = 0; i < num_loaded_molecules; ++i) {
        if (loaded_molecules[i].key == key) return &loaded_molecules[i];
    }
    return nullptr;
}

static inline void add_loaded_molecule(LoadedMolecule obj) {
    ASSERT(!find_loaded_molecule(obj.key));
    ASSERT(num_loaded_molecules < ARRAY_SIZE(loaded_molecules));
    loaded_molecules[num_loaded_molecules++] = obj;
}

static inline void remove_loaded_molecule(uint64_t key) {
    for (uint64_t i = 0; i < num_loaded_molecules; ++i) {
        if (loaded_molecules[i].key == key) {
            loaded_molecules[i] = loaded_molecules[--num_loaded_molecules];
            return;
        }
    }
    ASSERT(false);
}

static LoadedTrajectory loaded_trajectories[8] = {};
static int64_t num_loaded_trajectories = 0;

static inline LoadedTrajectory* find_loaded_trajectory(uint64_t key) {
    for (uint64_t i = 0; i < num_loaded_trajectories; ++i) {
        if (loaded_trajectories[i].key == key) return &loaded_trajectories[i];
    }
    return nullptr;
}

static inline void add_loaded_trajectory(const LoadedTrajectory *traj) {
    ASSERT(!find_loaded_trajectory(traj->key));
    ASSERT(num_loaded_trajectories < ARRAY_SIZE(loaded_trajectories));
    memcpy(&loaded_trajectories[num_loaded_trajectories++], traj, sizeof(LoadedTrajectory));
}

static inline void remove_loaded_trajectory(uint64_t key) {
    for (uint64_t i = 0; i < num_loaded_trajectories; ++i) {
        if (loaded_trajectories[i].key == key) {
            memcpy(&loaded_trajectories[i], &loaded_trajectories[--num_loaded_trajectories], sizeof(LoadedTrajectory));
            return;
        }
    }
    ASSERT(false);
}

namespace load {

namespace mol {

bool is_extension_supported(str_t filename) {
    
    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) return true;
    if (compare_str_cstr(ext, "gro")) return true;

    return false;
}


bool load_string_pdb(md_molecule_t* mol, str_t string, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    md_pdb_data_t data = {0};
    defer { md_pdb_data_free(&data, alloc); };
    if (md_pdb_data_parse_str(string, &data, alloc) && md_pdb_molecule_init(mol, &data, alloc)) {
        LoadedMolecule obj = {
            .key = (uint64_t)mol,
            .extension = "pdb",
            .alloc = alloc
        };

        add_loaded_molecule(obj);
        return true;
    }
    return false;
}

bool load_string_gro(md_molecule_t* mol, str_t string, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    md_gro_data_t data = {0};
    defer { md_gro_data_free(&data, alloc); };
    if (md_gro_data_parse_str(string, &data, alloc) && md_gro_molecule_init(mol, &data, alloc)) {
        LoadedMolecule obj = {
            .key = (uint64_t)mol,
            .extension = "gro",
            .alloc = alloc
        };

        add_loaded_molecule(obj);
        return true;
    }
    return false;
}

bool load_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) {
        md_pdb_data_t data = {0};
        defer { md_pdb_data_free(&data, alloc); };
        if (md_pdb_data_parse_file(filename, &data, alloc) && md_pdb_molecule_init(mol, &data, alloc)) {
            goto success;
        }
    }
    else if (compare_str_cstr(ext, "gro")) {
        md_gro_data_t data = {0};
        defer { md_gro_data_free(&data, alloc); };
        if (md_gro_data_parse_file(filename, &data, alloc) && md_gro_molecule_init(mol, &data, alloc)) {
            goto success;
        }
    }

    md_printf(MD_LOG_TYPE_ERROR, "Unsupported file extension: '%.*s'", filename.len, filename.ptr);
    return false;

success:
    LoadedMolecule obj = {
        .key = (uint64_t)mol,
        .extension = {ext.ptr, ext.len},
        .alloc = alloc
    };
    add_loaded_molecule(obj);
    return true;
}

bool free(md_molecule_t* mol) {
    LoadedMolecule* loaded_mol = find_loaded_molecule((uint64_t)mol);
    if (loaded_mol) {
        CStringView ext = loaded_mol->extension;
        md_allocator_i* alloc = loaded_mol->alloc;
        remove_loaded_molecule((uint64_t)mol);
        if (ext == "pdb") {
            return md_pdb_molecule_free(mol, alloc);
        }
        else if (ext == "gro") {
            return md_gro_molecule_free(mol, alloc);
        }
        ASSERT(false);
    }
    md_print(MD_LOG_TYPE_ERROR, "Attempting to free molecule which was not loaded with loader");
    return false;
}

}  // namespace mol

namespace traj {

bool is_extension_supported(str_t filename) {
    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) return true;
    if (compare_str_cstr(ext, "xtc")) return true;

    return false;
}

void init_frame_data(FrameData* frame_data, int64_t num_frames, int64_t num_atoms, md_allocator_i* alloc) {
    for (int64_t i = 0; i < num_frames; ++i) {
        frame_data[i].num_atoms = num_atoms;
        frame_data[i].x = (float*)md_alloc(alloc, num_atoms * sizeof(float));
        frame_data[i].y = (float*)md_alloc(alloc, num_atoms * sizeof(float));
        frame_data[i].z = (float*)md_alloc(alloc, num_atoms * sizeof(float));
    }   
}

void free_frame_data(FrameData* frame_data, int64_t num_frames, md_allocator_i* alloc) {
    for (int64_t i = 0; i < num_frames; ++i) {
        md_free(alloc, frame_data[i].x, frame_data[i].num_atoms * sizeof(float));
        md_free(alloc, frame_data[i].y, frame_data[i].num_atoms * sizeof(float));
        md_free(alloc, frame_data[i].z, frame_data[i].num_atoms * sizeof(float));
    }
    memset(frame_data, 0, sizeof(FrameData) * num_frames);
}

bool open_file(md_trajectory_i* traj, str_t filename, const md_molecule_t* mol, md_allocator_i* alloc) {
    ASSERT(traj);
    ASSERT(mol);
    ASSERT(alloc);

    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) {
        if (md_pdb_trajectory_open(traj, filename, alloc)) {
            goto success;
        }
    }
    else if (compare_str_cstr(ext, "xtc")) {
        if (md_xtc_trajectory_open(traj, filename, alloc)) {
            goto success;
        }
    }

    md_printf(MD_LOG_TYPE_ERROR, "Unsupported file extension: '%.*s'", filename.len, filename.ptr);
    return false;

success:
    LoadedTrajectory obj = {
        .key = (uint64_t)traj,
        .extension = {ext.ptr, ext.len},
        .mol = mol,
        .alloc = alloc,
        .frame_data = {0},
        .frame_cache = {
            .data = {0, 1, 2, 3, 4, 5, 6, 7},
            .keys = {-1, -1, -1, -1, -1, -1, -1, -1},
        },
    };
    
    init_frame_data(obj.frame_data, ARRAY_SIZE(obj.frame_data), traj->num_atoms, default_allocator);
    add_loaded_trajectory(&obj);
    return true;
}

bool close(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        CStringView ext = loaded_traj->extension;
        free_frame_data(loaded_traj->frame_data, ARRAY_SIZE(loaded_traj->frame_data), default_allocator);
        remove_loaded_trajectory(loaded_traj->key);

        if (ext == "pdb") {
            return md_pdb_trajectory_close(traj);
        }
        else if (ext == "xtc") {
            return md_xtc_trajectory_close(traj);
        }
        ASSERT(false);
    }
    md_print(MD_LOG_TYPE_ERROR, "Attempting to free trajectory which was not loaded with loader");
    return false;
}

bool load_frame_data(const md_molecule_t* mol, const md_trajectory_i* traj, int64_t frame_idx, FrameData* dst) {
    ASSERT(mol);
    ASSERT(traj);
    ASSERT(dst);

    md_trajectory_frame_header_t header = {0};
    if (md_trajectory_load_frame(traj, frame_idx, &header, dst->x, dst->y, dst->z, dst->num_atoms)) {
        memcpy(&dst->box, header.box, sizeof(dst->box));
        dst->timestamp = header.timestamp;

        md_util_apply_pbc_args_t args = {
            .atom = {
                .count = traj->num_atoms,
                .x = dst->x,
                .y = dst->y,
                .z = dst->z,
            },
            .residue = {
                .count = mol->residue.count,
                .atom_range = mol->residue.atom_range,
            },
            .chain = {
                .count = mol->chain.count,
                .residue_range = mol->chain.residue_range,
            },
        };
        memcpy(&args.pbc.box, &dst->box, sizeof(args.pbc.box));
        md_util_apply_pbc(dst->x, dst->y, dst->z, dst->num_atoms, args);

        return true;
    }
    return false;
}

bool fetch_trajectory_frame_data(const md_trajectory_i* traj, i64 frame_idx, float* x, float* y, float* z, int64_t num_coords, float* box[3][3]) {
    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    ASSERT(loaded_traj);

    // Use some kind of synchronization mechanism here so contenting threads don't push cached frames out from each other
    //loaded_traj->semaphore.acquire();
    //defer { loaded_traj->semaphore.release(); };

    int* slot = loaded_traj->frame_cache.get(frame_idx);
    if (!slot) {
        // not in cache
        slot = loaded_traj->frame_cache.reserve(frame_idx);
        ASSERT(slot);
        FrameData& fd = loaded_traj->frame_data[*slot];
        if (!load_frame_data(loaded_traj->mol, traj, frame_idx, &fd)) {
            md_print(MD_LOG_TYPE_ERROR, "Severe error occured, failed to load frame_data");
            return false;
        }
    }

    const FrameData* data = &loaded_traj->frame_data[*slot];

    if (num_coords > 0) {
        ASSERT(x);
        ASSERT(y);
        ASSERT(z);

        if (num_coords != data->num_atoms) {
            md_print(MD_LOG_TYPE_ERROR, "Number of coordinates requested differs from number of coordinates stored in trajectory frame");
        }
        const int64_t coord_copy_count = MIN(num_coords, data->num_atoms);
        memcpy(x, data->x, coord_copy_count * sizeof(float));
        memcpy(y, data->y, coord_copy_count * sizeof(float));
        memcpy(z, data->z, coord_copy_count * sizeof(float));
    }

    if (box) memcpy(box, &data->box, sizeof(data->box));

    return true;
}

bool load_trajectory_frame_box(const md_trajectory_i* traj, float* box[3][3], i64 frame_idx) {
    ASSERT(traj);
    ASSERT(traj->inst);
    return fetch_trajectory_frame_data(traj, frame_idx, NULL, NULL, NULL, 0, box);
}

bool load_trajectory_frame_coords(const md_trajectory_i* traj, float* x, float* y, float* z, i64 atom_count, i64 frame_idx) {
    ASSERT(traj);
    ASSERT(traj->inst);
    return fetch_trajectory_frame_data(traj, frame_idx, x, y, z, atom_count, NULL);
}


}  // namespace traj

}  // namespace load