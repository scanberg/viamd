#include "loader.h"
//#include <core/common.h>
//#include <core/log.h>
//#include <core/file.h>
//#include <core/string_types.h>
//#include <core/string_utils.h>
//#include <core/sync.h>

#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_frame_cache.h>
#include <md_util.h>

#include <string_util.h>
#include <string.h>

#include "task_system.h"

struct LoadedMolecule {
    uint64_t key;
    StrBuf<8> extension;
    md_allocator_i* alloc;
};

struct LoadedTrajectory {
    uint64_t key;
    StrBuf<8> extension;
    const md_molecule_t* mol;
    md_trajectory_i traj;
    md_frame_cache_t cache;
    md_allocator_i* alloc;
};

static LoadedMolecule loaded_molecules[8] = {};
static int64_t num_loaded_molecules = 0;

static LoadedTrajectory loaded_trajectories[8] = {};
static int64_t num_loaded_trajectories = 0;

static inline LoadedMolecule* find_loaded_molecule(uint64_t key) {
    for (int64_t i = 0; i < num_loaded_molecules; ++i) {
        if (loaded_molecules[i].key == key) return &loaded_molecules[i];
    }
    return nullptr;
}

static inline void add_loaded_molecule(LoadedMolecule obj) {
    ASSERT(!find_loaded_molecule(obj.key));
    ASSERT(num_loaded_molecules < (int64_t)ARRAY_SIZE(loaded_molecules));
    loaded_molecules[num_loaded_molecules++] = obj;
}

static inline void remove_loaded_molecule(uint64_t key) {
    for (int64_t i = 0; i < num_loaded_molecules; ++i) {
        if (loaded_molecules[i].key == key) {
            loaded_molecules[i] = loaded_molecules[--num_loaded_molecules];
            return;
        }
    }
    ASSERT(false);
}

static inline LoadedTrajectory* find_loaded_trajectory(uint64_t key) {
    for (int64_t i = 0; i < num_loaded_trajectories; ++i) {
        if (loaded_trajectories[i].key == key) return &loaded_trajectories[i];
    }
    return nullptr;
}

static inline void add_loaded_trajectory(LoadedTrajectory traj) {
    ASSERT(!find_loaded_trajectory(traj.key));
    ASSERT(num_loaded_trajectories < (int64_t)ARRAY_SIZE(loaded_trajectories));
    loaded_trajectories[num_loaded_trajectories++] = traj;
}

static inline void remove_loaded_trajectory(uint64_t key) {
    for (int64_t i = 0; i < num_loaded_trajectories; ++i) {
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
    md_pdb_data_t data = {};
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
    md_gro_data_t data = {};
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
        md_pdb_data_t data = {};
        defer { md_pdb_data_free(&data, alloc); };
        if (md_pdb_data_parse_file(filename, &data, alloc) && md_pdb_molecule_init(mol, &data, alloc)) {
            goto success;
        }
    }
    else if (compare_str_cstr(ext, "gro")) {
        md_gro_data_t data = {};
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
        str_t ext = loaded_mol->extension;
        md_allocator_i* alloc = loaded_mol->alloc;
        remove_loaded_molecule((uint64_t)mol);
        if (compare_str_cstr(ext, "pdb")) {
            return md_pdb_molecule_free(mol, alloc);
        }
        else if (compare_str_cstr(ext, "gro")) {
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

int64_t fetch_frame_data(struct md_trajectory_o* inst, int64_t idx, void* data_ptr) {
    if (data_ptr) {
        *((int64_t*)data_ptr) = idx;
    }
    return sizeof(int64_t);
}

bool decode_frame_data(struct md_trajectory_o* inst, const void* data_ptr, int64_t data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    ASSERT(loaded_traj);
    ASSERT(data_size == sizeof(int64_t));

    int64_t idx = (int64_t)data_ptr;
    ASSERT(0 <= idx && idx < md_trajectory_num_frames(&loaded_traj->traj));

    md_frame_data_t* frame_data;
    md_frame_cache_lock_t* lock;
    bool result = true;
    if (md_frame_cache_reserve_frame(&loaded_traj->cache, idx, &frame_data, &lock)) {
        const int64_t frame_data_size = md_trajectory_fetch_frame_data(&loaded_traj->traj, idx, 0);
        void* frame_data_ptr = md_alloc(default_temp_allocator, frame_data_size);
        md_trajectory_fetch_frame_data(&loaded_traj->traj, idx, frame_data_ptr);

        result = md_trajectory_decode_frame_data(&loaded_traj->traj, frame_data_ptr, frame_data_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

        if (result) {
            // Deperiodize
            const md_molecule_t& mol = *loaded_traj->mol;
            md_util_apply_pbc_args_t args = {
                .atom = {
                    .count = mol.atom.count,
                    .x = frame_data->x,
                    .y = frame_data->y,
                    .z = frame_data->z,
                },
                .residue = {
                        .count = mol.residue.count,
                        .atom_range = mol.residue.atom_range,
                },
                .chain = {
                        .count = mol.chain.count,
                        .residue_range = mol.chain.residue_range,
                }
            };
            memcpy(args.pbc.box, frame_data->header.box, sizeof(args.pbc.box));
            md_util_apply_pbc(frame_data->x, frame_data->y, frame_data->z, mol.atom.count, args);
        }
    }

    if (result) {
        const int64_t num_atoms = frame_data->header.num_atoms;
        if (header) memcpy(header, &frame_data->header, sizeof(md_trajectory_header_t));
        if (x) memcpy(x, frame_data->x, sizeof(float) * num_atoms);
        if (y) memcpy(y, frame_data->y, sizeof(float) * num_atoms);
        if (z) memcpy(z, frame_data->z, sizeof(float) * num_atoms);
    }

    md_frame_cache_release_frame_lock(lock);

    return result;
}

bool load_frame(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    void* frame_data = &idx;
    return decode_frame_data(inst, frame_data, sizeof(int64_t), header, x, y, z);
}

void launch_prefetch_job() {
    #if 0
#define NUM_SLOTS 64

    // This is the number of slots we have to work with in parallel.
    // This number should ideally be more than the number of cores available.
    // We pre-allocate the number of slots * max frame data size, so don't go bananas here if you want to save some on memory.
    task_system::enqueue_pool("Preloading frames", 1, [data, num_frames](task_system::TaskSetRange)
        {
            timestamp_t t0 = md_os_time_current();
            const int64_t slot_size = md_trajectory_max_frame_data_size(&data->mold.traj);
            void* slot_mem = md_alloc(default_allocator, slot_size * NUM_SLOTS);

            void* slots[NUM_SLOTS];
            atomic_queue::AtomicQueue<uint32_t, NUM_SLOTS, 0xFFFFFFFF> slot_queue;

            for (uint32_t i = 0; i < NUM_SLOTS; ++i) {
                slots[i] = (char*)slot_mem + i * slot_size;
                slot_queue.push(i);
            }

            // Iterate over all frames and load the raw data, then spawn a task for each frame
            for (uint32_t i = 0; i < num_frames; ++i) {
                uint32_t slot_idx = slot_queue.pop();
                ASSERT(slot_idx < NUM_SLOTS);

                void* frame_mem = slots[slot_idx];
                const int64_t frame_size = data->mold.traj.fetch_frame_data(data->mold.traj.inst, i, NULL);
                data->mold.traj.fetch_frame_data(data->mold.traj.inst, i, frame_mem);

                // Spawn task: Load, Decode and postprocess
                //printf("Spawning task to decode frame %i\n", i);
                auto id = task_system::enqueue_pool("##Decode frame", 1, [slot_idx, frame_mem, frame_size, frame_idx = i, &slot_queue, data](task_system::TaskSetRange)
                    {
                        md_frame_data_t* frame_data;
                        md_frame_cache_lock_t* lock;
                        if (md_frame_cache_reserve_frame(&data->mold.frame_cache, frame_idx, &frame_data, &lock)) {
                            data->mold.traj.decode_frame_data(data->mold.traj.inst, frame_mem, frame_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

                            // Free the data slot directly here
                            slot_queue.push(slot_idx);
                            
                            // deperiodize
                            const md_molecule_t& mol = data->mold.mol;
                            md_util_apply_pbc_args_t args = {
                                .atom = {
                                    .count = mol.atom.count,
                                    .x = frame_data->x,
                                    .y = frame_data->y,
                                    .z = frame_data->z,
                                },
                                .residue = {
                                        .count = mol.residue.count,
                                        .atom_range = mol.residue.atom_range,
                                },
                                .chain = {
                                        .count = mol.chain.count,
                                        .residue_range = mol.chain.residue_range,
                                }
                            };
                            memcpy(args.pbc.box, frame_data->header.box, sizeof(args.pbc.box));
                            md_util_apply_pbc(frame_data->x, frame_data->y, frame_data->z, mol.atom.count, args);

                            if (mol.backbone.count > 0) {
                                md_util_backbone_angle_args_t bb_args = {
                                    .atom = {
                                        .count = mol.atom.count,
                                        .x = frame_data->x,
                                        .y = frame_data->y,
                                        .z = frame_data->z,
                                    },
                                    .backbone = {
                                            .count = mol.backbone.count,
                                            .atoms = mol.backbone.atoms,
                                    },
                                    .chain = {
                                            .count = mol.chain.count,
                                            .backbone_range = mol.chain.backbone_range,
                                    }
                                };
                                md_util_compute_backbone_angles(data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx, data->trajectory_data.backbone_angles.stride, &bb_args);

                                md_util_secondary_structure_args_t ss_args = {
                                    .atom = {
                                        .count = data->mold.mol.atom.count,
                                        .x = frame_data->x,
                                        .y = frame_data->y,
                                        .z = frame_data->z,
                                    },
                                    .backbone = {
                                            .count = mol.backbone.count,
                                            .atoms = mol.backbone.atoms,
                                    },
                                    .chain = {
                                            .count = data->mold.mol.chain.count,
                                            .backbone_range = data->mold.mol.chain.backbone_range,
                                    }
                                };
                                md_util_compute_secondary_structure(data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx, data->trajectory_data.secondary_structure.stride, &ss_args);
                            }

                            md_frame_cache_release_frame_lock(lock);
                        } else {
                            slot_queue.push(slot_idx);
                        }
                        //printf("Finished decoding frame %i\n", frame_idx);
                    }
                );
            }
            
            //printf("Sitting back waiting for tasks to complete...\n");
            while (!slot_queue.was_full()) {
                _mm_pause(); // Back off for a bit
            }

            md_free(default_allocator, slot_mem, slot_size * NUM_SLOTS);

            timestamp_t t1 = md_os_time_current();
            md_printf(MD_LOG_TYPE_INFO, "Frame preload took %.2f seconds", md_os_time_delta_in_s(t0, t1));
        });
#undef NUM_SLOTS
#endif
}

bool open_file(md_trajectory_i* traj, str_t filename, const md_molecule_t* mol, md_allocator_i* alloc) {
    ASSERT(traj);
    ASSERT(mol);
    ASSERT(alloc);

    md_trajectory_i internal_traj = {0};

    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) {
        if (md_pdb_trajectory_open(&internal_traj, filename, alloc)) {
            goto success;
        }
    }
    else if (compare_str_cstr(ext, "xtc")) {
        if (md_xtc_trajectory_open(&internal_traj, filename, alloc)) {
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
        .traj = internal_traj,
        .cache = {0},
        .alloc = alloc,
    };
    
    md_frame_cache_init(&obj.cache, &obj.traj, alloc, MEGABYTES(512));

    // We only overload load frame and decode frame data to apply PBC upon loading data
    *traj = obj.traj;
    traj->load_frame = load_frame;
    traj->fetch_frame_data = fetch_frame_data;
    traj->decode_frame_data = decode_frame_data;
    
    add_loaded_trajectory(obj);
    return true;
}

bool close(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        str_t ext = loaded_traj->extension;
        remove_loaded_trajectory(loaded_traj->key);

        if (compare_str_cstr(ext, "pdb")) {
            return md_pdb_trajectory_close(traj);
        }
        else if (compare_str_cstr(ext, "xtc")) {
            return md_xtc_trajectory_close(traj);
        }
        ASSERT(false);
    }
    md_print(MD_LOG_TYPE_ERROR, "Attempting to free trajectory which was not loaded with loader");
    return false;
}

}  // namespace traj

}  // namespace load
