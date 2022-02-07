#include "loader.h"

#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_bitfield.h>
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
    md_exp_bitfield_t recenter_target;
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

static inline LoadedTrajectory* alloc_loaded_trajectory(uint64_t key) {
    ASSERT(find_loaded_trajectory(key) == NULL);
    ASSERT(num_loaded_trajectories < (int64_t)ARRAY_SIZE(loaded_trajectories));
    LoadedTrajectory* traj = &loaded_trajectories[num_loaded_trajectories++];
    *traj = {0}; // Clear
    traj->key = key;
    return traj;
}

static inline void remove_loaded_trajectory(uint64_t key) {
    for (int64_t i = 0; i < num_loaded_trajectories; ++i) {
        if (loaded_trajectories[i].key == key) {
            md_frame_cache_free(&loaded_trajectories[i].cache);
            md_trajectory_free(&loaded_trajectories[i].traj);
            loaded_trajectories[i] = loaded_trajectories[--num_loaded_trajectories];
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
    if (md_gro_data_parse_str(&data, string, alloc) && md_gro_molecule_init(mol, &data, alloc)) {
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
        if (md_gro_data_parse_file(&data, filename, alloc) && md_gro_molecule_init(mol, &data, alloc)) {
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

md_allocator_i* get_internal_allocator(md_molecule_t* mol) {
    LoadedMolecule* loaded_mol = find_loaded_molecule((uint64_t)mol);
    if (loaded_mol) {
        str_t ext = loaded_mol->extension;
        if (compare_str_cstr(ext, "pdb")) {
            return md_pdb_molecule_internal_allocator(mol);
        }
        else if (compare_str_cstr(ext, "gro")) {
            return md_gro_molecule_internal_allocator(mol);
        }
        ASSERT(false);
    }
    md_print(MD_LOG_TYPE_ERROR, "Attempting to get allocator for molecule which was not loaded with loader");
    return NULL;
}

}  // namespace mol

namespace traj {

bool is_extension_supported(str_t filename) {
    str_t ext = extract_ext(filename);
    if (compare_str_cstr(ext, "pdb")) return true;
    if (compare_str_cstr(ext, "xtc")) return true;

    return false;
}

bool get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    return md_trajectory_get_header(&loaded_traj->traj, header);
}

int64_t fetch_frame_data(struct md_trajectory_o*, int64_t idx, void* data_ptr) {
    if (data_ptr) {
        *((int64_t*)data_ptr) = idx;
    }
    return sizeof(int64_t);
}

bool decode_frame_data(struct md_trajectory_o* inst, const void* data_ptr, int64_t data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    ASSERT(loaded_traj);
    ASSERT(data_size == sizeof(int64_t));

    int64_t idx = *((int64_t*)data_ptr);
    ASSERT(0 <= idx && idx < md_trajectory_num_frames(&loaded_traj->traj));

    md_frame_data_t* frame_data;
    md_frame_cache_lock_t* lock = 0;
    bool result = true;
    bool in_cache = md_frame_cache_find_or_reserve(&loaded_traj->cache, idx, &frame_data, &lock);
    if (!in_cache) {
        const int64_t frame_data_size = md_trajectory_fetch_frame_data(&loaded_traj->traj, idx, 0);
        void* frame_data_ptr = md_alloc(default_temp_allocator, frame_data_size);
        md_trajectory_fetch_frame_data(&loaded_traj->traj, idx, frame_data_ptr);

        result = md_trajectory_decode_frame_data(&loaded_traj->traj, frame_data_ptr, frame_data_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

        bool have_box = (frame_data->header.box[0][0] + frame_data->header.box[1][1] + frame_data->header.box[2][2]) > 0;

        if (result) {
            // If we have a recenter target, then compute and apply that transformation
            if (!md_bitfield_empty(&loaded_traj->recenter_target)) {
                int64_t count = md_bitfield_popcount(&loaded_traj->recenter_target);
                if (count > 0) {
                    // Allocate data for substructure
                    int64_t stride = ROUND_UP(count, md_simd_widthf);
                    float* mem = (float*)md_alloc(default_temp_allocator, stride * 4 * sizeof(float));
                    float* tmp_x = mem + 0 * stride;
                    float* tmp_y = mem + 1 * stride;
                    float* tmp_z = mem + 2 * stride;
                    float* tmp_w = mem + 3 * stride;

                    // Extract fields for substructure
                    int64_t dst_idx = 0;
                    int64_t beg_bit = loaded_traj->recenter_target.beg_bit;
                    int64_t end_bit = loaded_traj->recenter_target.end_bit;
                    while ((beg_bit = md_bitfield_scan(&loaded_traj->recenter_target, beg_bit, end_bit)) != 0) {
                        int64_t src_idx = beg_bit - 1;
                        tmp_x[dst_idx] = frame_data->x[src_idx];
                        tmp_y[dst_idx] = frame_data->y[src_idx];
                        tmp_z[dst_idx] = frame_data->z[src_idx];
                        tmp_w[dst_idx] = loaded_traj->mol->atom.mass[src_idx];
                        dst_idx += 1;
                    }

                    // Compute deperiodized com for substructure
                    vec3_t com;
                    if (have_box)
                        com = md_util_compute_periodic_com(tmp_x, tmp_y, tmp_z, tmp_w, count, frame_data->header.box);
                    else {
                        com = md_util_compute_com(tmp_x, tmp_y, tmp_z, tmp_w, count);
                    }

                    vec3_t ext = {0,0,0};
                    if (have_box) {
                        ext = {frame_data->header.box[0][0], frame_data->header.box[1][1], frame_data->header.box[2][2]};
                    }

                    // Translate all
                    vec3_t trans = ext * 0.5f - com;
                    for (int64_t i = 0; i < frame_data->header.num_atoms; ++i) {
                        frame_data->x[i] += trans.x;
                        frame_data->y[i] += trans.y;
                        frame_data->z[i] += trans.z;
                    }
                }
            }

            // Deperiodize
            if (have_box) {
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
    }

    if (result) {
        const int64_t num_atoms = frame_data->header.num_atoms;
        if (header) *header = frame_data->header;
        if (x) memcpy(x, frame_data->x, sizeof(float) * num_atoms);
        if (y) memcpy(y, frame_data->y, sizeof(float) * num_atoms);
        if (z) memcpy(z, frame_data->z, sizeof(float) * num_atoms);
    }

    if (lock) {
        md_frame_cache_release_frame_lock(lock);
    }

    return result;
}

bool load_frame(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    void* frame_data = &idx;
    return decode_frame_data(inst, frame_data, sizeof(int64_t), header, x, y, z);
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
    LoadedTrajectory* inst = alloc_loaded_trajectory((uint64_t)traj);
    inst->extension = ext;
    inst->mol = mol;
    inst->traj = internal_traj;
    inst->cache = {0};
    inst->recenter_target = {0};
    inst->alloc = alloc;
    
    md_frame_cache_init(&inst->cache, &inst->traj, alloc, md_trajectory_num_frames(&internal_traj));
    md_bitfield_init(&inst->recenter_target, alloc);

    // We only overload load frame and decode frame data to apply PBC upon loading data
    traj->inst = (struct md_trajectory_o*)inst;
    traj->get_header = get_header;
    traj->load_frame = load_frame;
    traj->fetch_frame_data = fetch_frame_data;
    traj->decode_frame_data = decode_frame_data;
    
    return true;
}

bool close(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        remove_loaded_trajectory(loaded_traj->key);
        memset(traj, 0, sizeof(md_trajectory_i));
        return true;
    }
    md_print(MD_LOG_TYPE_ERROR, "Attempting to free trajectory which was not loaded with loader");
    ASSERT(false);
    return false;
}

bool set_recenter_target(md_trajectory_i* traj, const md_exp_bitfield_t* atom_mask) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        if (atom_mask) {
            md_bitfield_copy(&loaded_traj->recenter_target, atom_mask);
        }
        else {
            md_bitfield_clear(&loaded_traj->recenter_target);
        }
        return true;
    }
    md_print(MD_LOG_TYPE_ERROR, "Supplied trajectory was not loaded with loader");
    return false;
}

bool clear_cache(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        md_frame_cache_clear(&loaded_traj->cache);
        return true;
    }
    md_print(MD_LOG_TYPE_ERROR, "Supplied trajectory was not loaded with loader");
    return false;
}

}  // namespace traj

}  // namespace load
