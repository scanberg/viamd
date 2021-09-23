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
#include <md_util.h>

#include <string_util.h>

struct LoadedMolecule {
    uint64_t key;
    StrBuf<8> extension;
    md_allocator_i* alloc;
};

struct LoadedTrajectory {
    uint64_t key;
    StrBuf<8> extension;
    const md_molecule_t* mol;
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
    ASSERT(num_loaded_molecules < ARRAY_SIZE(loaded_molecules));
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

static inline void add_loaded_trajectory(const LoadedTrajectory *traj) {
    ASSERT(!find_loaded_trajectory(traj->key));
    ASSERT(num_loaded_trajectories < ARRAY_SIZE(loaded_trajectories));
    memcpy(&loaded_trajectories[num_loaded_trajectories++], traj, sizeof(LoadedTrajectory));
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
    };
    
    add_loaded_trajectory(&obj);
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