#include "loader.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_bitfield.h>
#include <core/md_os.h>
#include <core/md_parse.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_dcd.h>
#include <md_trajectory.h>
#include <md_util.h>
#if MD_VLX
#include <md_vlx.h>
#endif

#include <task_system.h>

enum sys_loader_t {
    SYS_LOADER_UNKNOWN,
    SYS_LOADER_PDB,
    SYS_LOADER_GRO,
    SYS_LOADER_XYZ,
    SYS_LOADER_CIF,
    SYS_LOADER_LAMMPS,
#if MD_VLX
    SYS_LOADER_VELOXCHEM,
#endif
    SYS_LOADER_COUNT
};

static const str_t sys_loader_name[] {
    STR_LIT("Unknown"),
    STR_LIT("Standard Protein Data Bank (pdb)"),
    STR_LIT("Gromacs Structure (gro)"),
    STR_LIT("XYZ"),
    STR_LIT("PDBx/mmCIF (cif)"),
    STR_LIT("LAMMPS (data)"),
#if MD_VLX
    STR_LIT("VeloxChem"),
#endif
};

static const str_t sys_loader_ext[] {
    {},
    STR_LIT("pdb"),
    STR_LIT("gro"),
    STR_LIT("xyz;xmol;arc"),
    STR_LIT("cif"),
    STR_LIT("data"),
#if MD_VLX
    STR_LIT("out;h5"),
#endif
};

static md_system_loader_i* sys_loader[] = {
    NULL,
    md_pdb_system_loader(),
    md_gro_system_loader(),
    md_xyz_system_loader(),
    md_mmcif_system_loader(),
    md_lammps_system_loader(),
#if MD_VLX
    md_vlx_system_loader(),
#endif
};

enum traj_loader_t {
    TRAJ_LOADER_UNKNOWN,
    TRAJ_LOADER_PDB,
    TRAJ_LOADER_XTC,
    TRAJ_LOADER_TRR,
    TRAJ_LOADER_XYZ,
    TRAJ_LOADER_LAMMPSTRJ,
    TRAJ_LOADER_DCD,
    TRAJ_LOADER_COUNT,
};

static str_t traj_loader_ext[] {
	{},
	STR_LIT("pdb"),
	STR_LIT("xtc"),
	STR_LIT("trr"),
	STR_LIT("xyz;xmol;arc"),
	STR_LIT("lammpstrj"),
    STR_LIT("dcd")
};

static md_traj_attach_fn traj_creator[] = {
	NULL,
    md_pdb_trajectory_create,
    md_xtc_trajectory_create,
    md_trr_trajectory_create,
    md_xyz_trajectory_create,
    md_lammps_trajectory_create,
    md_dcd_trajectory_create
};

struct LoadedMolecule {
    uint64_t key;
    md_allocator_i* alloc;
};

struct LoadedTrajectory {
    uint64_t key;
    const md_system_t* mol;
    md_trajectory_i* traj;
    md_allocator_i*  alloc;

    md_array(int32_t) recenter_indices;
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
            if (loaded_trajectories[i].traj && loaded_trajectories[i].traj->free) {
                loaded_trajectories[i].traj->free(loaded_trajectories[i].traj);
            }
            // Swap back and pop
            loaded_trajectories[i] = loaded_trajectories[--num_loaded_trajectories];
            return;
        }
    }
    ASSERT(false);
}

// In here each loader gets a chance to do a precheck with the file to be loaded
static void mol_loader_preload_check(load::LoaderState* state, sys_loader_t loader, str_t file_path, md_allocator_i* alloc) {
    switch (loader) {
    case SYS_LOADER_LAMMPS: {
        md_lammps_atom_format_t format = md_lammps_atom_format_from_file(file_path);
        if (format != MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
            // Encode this into the argument
            const char* format_str = md_lammps_atom_format_strings()[format];
            md_lammps_system_loader_arg_t arg = md_lammps_system_loader_arg(format_str);
            state->data_size = sizeof(arg);
            state->data_ptr = md_alloc(alloc, sizeof(arg));
            MEMCPY(state->data_ptr, &arg, sizeof(arg));
            state->sys_loader_arg = state->data_ptr;
        } else {
            state->flags |= LoaderStateFlag_RequiresDialogue;
        }
    }
        break;
    default:
        break;
    }
}

static void traj_loader_preload_check(load::LoaderState*, traj_loader_t, str_t, md_allocator_i*) {
    return;
}

namespace load {

#define NUM_ENTRIES 13
struct table_entry_t {
    str_t name[NUM_ENTRIES];
    str_t ext[NUM_ENTRIES];
    md_system_loader_i*     sys_loader[NUM_ENTRIES];
    md_trajectory_creator_fn traj_creator[NUM_ENTRIES];
    uint32_t flags[NUM_ENTRIES];
};

enum {
    FLAG_NONE = 0,
	FLAG_REQUIRES_DIALOGUE = 1 << 0,
};

static const table_entry_t table = {
    {
        STR_LIT("Standard Protein Data Bank (pdb)"),
        STR_LIT("Gromacs Structure (gro)"),
        STR_LIT("Gromacs Compressed Trajectory (xtc)"),
        STR_LIT("Gromacs Lossless Trajectory (trr)"),
        STR_LIT("xyz (xyz)"),
        STR_LIT("xyz (xmol)"),
        STR_LIT("xyz (arc)"),
        STR_LIT("PDBx/mmCIF (cif)"),
        STR_LIT("LAMMPS (data)"),
        STR_LIT("LAMMPS Trajectory (lammpstrj)"),
        STR_LIT("DCD Trajectory (dcd)"),
#if MD_VLX
        STR_LIT("VeloxChem (out)"),
        STR_LIT("VeloxChem (h5)")
#endif
    },
    {
        STR_LIT("pdb"),
        STR_LIT("gro"),
        STR_LIT("xtc"),
        STR_LIT("trr"),
        STR_LIT("xyz"),
        STR_LIT("xmol"),
        STR_LIT("arc"),
        STR_LIT("cif"),
        STR_LIT("data"),
        STR_LIT("lammpstrj"),
        STR_LIT("dcd"),
#if MD_VLX
        STR_LIT("out"),
        STR_LIT("h5")
#endif
    },
    { 
        md_pdb_system_loader(),
        md_gro_system_loader(),  
        NULL,
        NULL,
        md_xyz_system_loader(),
        md_xyz_system_loader(),
        md_xyz_system_loader(),
        md_mmcif_system_loader(),
        md_lammps_system_loader(),
        NULL,
        //NULL,
#if MD_VLX
        md_vlx_system_loader(),
        md_vlx_system_loader(),
#endif
    },
	{ 
        md_pdb_trajectory_create,
    	NULL,
        md_xtc_trajectory_create,
        md_trr_trajectory_create,
        md_xyz_trajectory_create,
        md_xyz_trajectory_create,
        md_xyz_trajectory_create,
    	NULL,
        NULL,
        md_lammps_trajectory_create,
        md_dcd_trajectory_create,
#if MD_VLX
        NULL,
#endif
    }
};

sys_loader_t mol_loader_from_ext(str_t ext) {
    str_t tok[16];
	for (size_t i = 1; i < SYS_LOADER_COUNT; ++i) {
        str_t exts = sys_loader_ext[i];
        size_t num_tok = extract_tokens_delim(tok, ARRAY_SIZE(tok), &exts, ';');
        for (size_t j = 0; j < num_tok; ++j) {
			if (str_eq_ignore_case(ext, tok[j])) {
				return (sys_loader_t)i;
			}
		}
	}
	return SYS_LOADER_UNKNOWN;
}

traj_loader_t traj_loader_from_ext(str_t ext) {
    str_t tok[16];
    for (size_t i = 1; i < TRAJ_LOADER_COUNT; ++i) {
        str_t exts = traj_loader_ext[i];
        size_t num_tok = extract_tokens_delim(tok, ARRAY_SIZE(tok), &exts, ';');
        for (size_t j = 0; j < num_tok; ++j) {
            if (str_eq_ignore_case(ext, tok[j])) {
                return (traj_loader_t)i;
            }
        }
    }
    return TRAJ_LOADER_UNKNOWN;
}

bool init_loader_state(LoaderState* state, str_t file_path, md_allocator_i* alloc) {
    ASSERT(state);
    MEMSET(state, 0, sizeof(LoaderState));
    sys_loader_t  sys   = SYS_LOADER_UNKNOWN;
    traj_loader_t traj  = TRAJ_LOADER_UNKNOWN;

    str_t ext = {0};
    if (extract_ext(&ext, file_path)) {
        sys = mol_loader_from_ext(ext);
        if (sys) {
            state->sys_loader = sys_loader[sys];
            mol_loader_preload_check(state, sys, file_path, alloc);
        }

        traj = traj_loader_from_ext(ext);
        if (traj) {
            state->traj_creator = traj_creator[traj];
            traj_loader_preload_check(state, traj, file_path, alloc);
        }
    }

    return (sys || traj);
}

void free_loader_state(LoaderState* state, md_allocator_i* alloc) {
    ASSERT(state);
    if (state->data_size) {
        md_free(alloc, state->data_ptr, state->data_size);
        state->data_size = 0;
        state->data_ptr  = 0;
    }
}

size_t loader_count() {
    return NUM_ENTRIES;
}

const str_t* loader_names() {
    return table.name;
}

const str_t* loader_extensions() {
	return table.ext;
}

namespace mol {

md_system_loader_i* loader_from_ext(str_t ext) {
    for (size_t i = 0; i < NUM_ENTRIES; ++i) {
    	if (str_eq(ext, table.ext[i])) {
            return table.sys_loader[i];
        }
    }
    return NULL;
}

bool loader_requires_dialogue(md_system_loader_i* loader) {
    if (loader) {
        for (size_t i = 0; i < NUM_ENTRIES; ++i) {
		    if (table.sys_loader[i] == loader) {
			    return table.flags[i] & FLAG_REQUIRES_DIALOGUE;
		    }
	    }
    }
    return false;
}

}  // namespace mol

namespace traj {

md_trajectory_creator_fn creator_from_ext(str_t ext) {
	for (size_t i = 0; i < NUM_ENTRIES; ++i) {
        if (str_eq(ext, table.ext[i])) {
            return table.traj_creator[i];
        }
    }
    return NULL;
}

md_trajectory_i* open_file(str_t filename, md_trajectory_creator_fn creator, const md_system_t* mol, md_allocator_i* alloc, LoadTrajectoryFlags flags) {
    ASSERT(mol);
    ASSERT(alloc);

    if (!creator) {
        str_t ext;
        if (extract_ext(&ext, filename)) {
            creator = creator_from_ext(ext);
        }
    }
    if (!creator) {
        MD_LOG_ERROR("Unsupported file extension: '%.*s'", filename.len, filename.ptr);
        return NULL;
    }

    return creator(filename, alloc, (md_trajectory_flags_t)flags);
}

bool close(md_trajectory_i* traj) {
    if (traj) {
        LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
        if (loaded_traj) {
            remove_loaded_trajectory(loaded_traj->key);
            return true;
        }
        if (traj->free) {
            traj->free(traj);
            return true;
        }
        MD_LOG_ERROR("Supplied trajectory cannot free itself");
    }
    return false;
}

md_trajectory_i* get_raw_trajectory(md_trajectory_i* traj) {
    if (traj) {
        LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
        if (loaded_traj) {
            return loaded_traj->traj;
        }
        return traj;
    }
	return nullptr;
}

bool has_recenter_target(md_trajectory_i* traj) {
    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        return md_array_size(loaded_traj->recenter_indices) > 0;
    }
    return false;
}

bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        if (atom_mask) {
            size_t count = md_bitfield_popcount(atom_mask);
            md_array_resize(loaded_traj->recenter_indices, count, loaded_traj->alloc);
            md_bitfield_iter_extract_indices(loaded_traj->recenter_indices, count, md_bitfield_iter_create(atom_mask));
        }
        else {
            md_array_shrink(loaded_traj->recenter_indices, 0);
        }
        return true;
    }
    MD_LOG_ERROR("Supplied trajectory is not tracked by the cache wrapper");
    return false;
}

}  // namespace traj

}  // namespace load
