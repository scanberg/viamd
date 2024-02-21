#include "loader.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_simd.h>
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
//#include <md_dcd.h>
#include <md_trajectory.h>
#include <md_frame_cache.h>
#include <md_util.h>
#if MD_VLX
#include <md_vlx.h>
#endif

#include <string.h>

#include "task_system.h"

enum mol_loader_t {
    MOL_LOADER_UNKNOWN,
    MOL_LOADER_PDB,
    MOL_LOADER_GRO,
    MOL_LOADER_XYZ,
    MOL_LOADER_CIF,
    MOL_LOADER_LAMMPS,
    MOL_LOADER_COUNT
};

static const str_t mol_loader_name[] {
    STR_LIT("Unknown"),
    STR_LIT("Standard Protein Data Bank (pdb)"),
    STR_LIT("Gromacs Structure (gro)"),
    STR_LIT("XYZ"),
    STR_LIT("PDBx/mmCIF (cif)"),
    STR_LIT("LAMMPS (data)"),
#if MD_VLX
    STR_LIT("Veloxchem (out)"),
#endif
};

static const str_t mol_loader_ext[] {
    {},
    STR_LIT("pdb"),
    STR_LIT("gro"),
    STR_LIT("xyz;xmol;arc"),
    STR_LIT("cif"),
    STR_LIT("data"),
#if MD_VLX
    STR_LIT("out"),
#endif
};

static md_molecule_loader_i* mol_loader_api[] = {
    NULL,
    md_pdb_molecule_api(),
    md_gro_molecule_api(),
    md_xyz_molecule_api(),
    md_mmcif_molecule_api(),
    md_lammps_molecule_api(),
    md_vlx_molecule_api(),
};

enum traj_loader_t {
    TRAJ_LOADER_UNKNOWN,
    TRAJ_LOADER_PDB,
    TRAJ_LOADER_XTC,
    TRAJ_LOADER_TRR,
    TRAJ_LOADER_XYZ,
    TRAJ_LOADER_LAMMPSTRJ,
    TRAJ_LOADER_COUNT,
};

static str_t traj_loader_name[] {
	STR_LIT("Unknown"),
	STR_LIT("Standard Protein Data Bank (pdb)"),
	STR_LIT("Gromacs Compressed Trajectory (xtc)"),
	STR_LIT("Gromacs Lossless Trajectory (trr)"),
	STR_LIT("XYZ"),
    STR_LIT("Lammps Trajectory [ASCII] (lammpstrj)"),
};

static str_t traj_loader_ext[] {
	{},
	STR_LIT("pdb"),
	STR_LIT("xtc"),
	STR_LIT("trr"),
	STR_LIT("xyz;xmol;arc"),
	STR_LIT("lammpstrj"),
};

static md_trajectory_loader_i* traj_loader_api[] = {
	NULL,
	md_pdb_trajectory_loader(),
	md_xtc_trajectory_loader(),
	md_trr_trajectory_loader(),
	md_xyz_trajectory_loader(),
    md_lammps_trajectory_loader(),
};

struct LoadedMolecule {
    uint64_t key;
    md_allocator_i* alloc;
};

struct LoadedTrajectory {
    uint64_t key;
    const md_molecule_t* mol;
    md_trajectory_loader_i* loader;
    md_trajectory_i* traj;
    md_frame_cache_t cache;
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
            md_frame_cache_free(&loaded_trajectories[i].cache);
            loaded_trajectories[i].loader->destroy(loaded_trajectories[i].traj);
            // Swap back and pop
            loaded_trajectories[i] = loaded_trajectories[--num_loaded_trajectories];
            return;
        }
    }
    ASSERT(false);
}

// In here each loader gets a chance to do a precheck with the file to be loaded
static void mol_loader_preload_check(load::LoaderState* state, mol_loader_t loader, str_t file_path, md_allocator_i* alloc) {
    switch (loader) {
    case MOL_LOADER_LAMMPS: {
        md_lammps_atom_format_t format = md_lammps_atom_format_from_file(file_path);
        if (format != MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
            // Encode this into the argument
            const char* format_str = md_lammps_atom_format_strings()[format];
            md_lammps_molecule_loader_arg_t arg = md_lammps_molecule_loader_arg(format_str);
            state->data_size = sizeof(arg);
            state->data_ptr = md_alloc(alloc, sizeof(arg));
            MEMCPY(state->data_ptr, &arg, sizeof(arg));
            state->mol_loader_arg = state->data_ptr;
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

#define NUM_ENTRIES 11
struct table_entry_t {
    str_t name[NUM_ENTRIES];
    str_t ext[NUM_ENTRIES];
    md_molecule_loader_i*   mol_loader[NUM_ENTRIES];
    md_trajectory_loader_i* traj_loader[NUM_ENTRIES];
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
        //STR_LIT("DCD Trajectory (dcd)"),
#if MD_VLX
        STR_LIT("Veloxchem (out)")
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
        //STR_LIT("dcd"),
#if MD_VLX
        STR_LIT("out")
#endif
    },
    { 
        md_pdb_molecule_api(),
        md_gro_molecule_api(),  
        NULL,
        NULL,
        md_xyz_molecule_api(),
        md_xyz_molecule_api(),
        md_xyz_molecule_api(),
        md_mmcif_molecule_api(),
        md_lammps_molecule_api(),
        NULL,
        //NULL,
#if MD_VLX
        md_vlx_molecule_api(),
#endif
    },
	{ 
        md_pdb_trajectory_loader(),
    	NULL,
    	md_xtc_trajectory_loader(),
    	md_trr_trajectory_loader(),
    	md_xyz_trajectory_loader(),
    	md_xyz_trajectory_loader(),
    	md_xyz_trajectory_loader(),
    	NULL,
        NULL,
        md_lammps_trajectory_loader(),
        //md_dcd_trajectory_loader(),
#if MD_VLX
        NULL,
#endif
    }
};

mol_loader_t mol_loader_from_ext(str_t ext) {
    str_t tok[16];
	for (size_t i = 1; i < MOL_LOADER_COUNT; ++i) {
        str_t exts = mol_loader_ext[i];
        size_t num_tok = extract_tokens_delim(tok, ARRAY_SIZE(tok), &exts, ';');
        for (size_t j = 0; j < num_tok; ++j) {
			if (str_eq_ignore_case(ext, tok[j])) {
				return (mol_loader_t)i;
			}
		}
	}
	return MOL_LOADER_UNKNOWN;
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
    mol_loader_t mol_loader   = MOL_LOADER_UNKNOWN;
    traj_loader_t traj_loader = TRAJ_LOADER_UNKNOWN;

    str_t ext = {0};
    if (extract_ext(&ext, file_path)) {
        mol_loader = mol_loader_from_ext(ext);
        if (mol_loader) {
            state->mol_loader = mol_loader_api[mol_loader];
            mol_loader_preload_check(state, mol_loader, file_path, alloc);
        }

        traj_loader = traj_loader_from_ext(ext);
        if (traj_loader) {
            state->traj_loader = traj_loader_api[traj_loader];
            traj_loader_preload_check(state, traj_loader, file_path, alloc);
        }
    }

    return (mol_loader || traj_loader);
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

md_molecule_loader_i* loader_from_ext(str_t ext) {
    for (size_t i = 0; i < NUM_ENTRIES; ++i) {
    	if (str_eq(ext, table.ext[i])) {
            return table.mol_loader[i];
        }
    }
    return NULL;
}

bool loader_requires_dialogue(md_molecule_loader_i* loader) {
    if (loader) {
        for (size_t i = 0; i < NUM_ENTRIES; ++i) {
		    if (table.mol_loader[i] == loader) {
			    return table.flags[i] & FLAG_REQUIRES_DIALOGUE;
		    }
	    }
    }
    return false;
}

}  // namespace mol

namespace traj {

md_trajectory_loader_i* loader_from_ext(str_t ext) {
	for (size_t i = 0; i < NUM_ENTRIES; ++i) {
        if (str_eq(ext, table.ext[i])) {
            return table.traj_loader[i];
        }
    }
    return NULL;
}

bool get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    return md_trajectory_get_header(loaded_traj->traj, header);
}

#if 0
size_t fetch_frame_data(struct md_trajectory_o*, int64_t idx, void* data_ptr) {
    if (data_ptr) {
        *((int64_t*)data_ptr) = idx;
    }
    return sizeof(int64_t);
}

bool decode_frame_data(struct md_trajectory_o* inst, const void* data_ptr, [[maybe_unused]] size_t data_size, md_trajectory_frame_header_t* header, float* out_x, float* out_y, float* out_z) {
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    ASSERT(loaded_traj);
    ASSERT(data_size == sizeof(int64_t));

    int64_t idx = *((int64_t*)data_ptr);
    ASSERT(0 <= idx && idx < (int64_t)md_trajectory_num_frames(loaded_traj->traj));

    md_frame_data_t* frame_data;
    md_frame_cache_lock_t* lock = 0;
    bool result = true;
    bool in_cache = md_frame_cache_find_or_reserve(&loaded_traj->cache, idx, &frame_data, &lock);
    if (!in_cache) {
        md_allocator_i* alloc = md_heap_allocator;
        size_t frame_data_size = md_trajectory_fetch_frame_data(loaded_traj->traj, idx, 0);
        void*  frame_data_ptr  = md_alloc(alloc, frame_data_size);
        md_trajectory_fetch_frame_data(loaded_traj->traj, idx, frame_data_ptr);
        result = md_trajectory_decode_frame_data(loaded_traj->traj, frame_data_ptr, frame_data_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

        if (result) {
            const md_unit_cell_t* cell = &frame_data->header.unit_cell;
            const md_molecule_t* mol = loaded_traj->mol;
            float* x = frame_data->x;
            float* y = frame_data->y;
            float* z = frame_data->z;
            const size_t num_atoms = frame_data->header.num_atoms;

            // If we have a recenter target, then compute the com and apply that transformation
            if (md_array_size(loaded_traj->recenter_indices) > 0) {
                size_t count = md_array_size(loaded_traj->recenter_indices);
                const int32_t* indices = loaded_traj->recenter_indices;

                vec3_t com = {0};
                if (count == 1) {
                    const int32_t i = indices[0];
					com = vec3_set(x[i], y[i], z[i]);
				} else {
					com = md_util_com_compute(x, y, z, mol->atom.mass, indices, count, &mol->unit_cell);
                    md_util_pbc(&com.x, &com.y, &com.z, 0, 1, cell);
                }

                // Translate all
                const vec3_t center = cell->flags ? cell->basis * vec3_set1(0.5f) : vec3_zero();
                const vec3_t trans  = center - com;
                vec3_batch_translate_inplace(x, y, z, num_atoms, trans);
            }
        }

        md_free(alloc, frame_data_ptr, frame_data_size);
    }

    if (result) {
        const int64_t num_atoms = frame_data->header.num_atoms;
        if (header) *header = frame_data->header;
        if (out_x) MEMCPY(out_x, frame_data->x, sizeof(float) * num_atoms);
        if (out_y) MEMCPY(out_y, frame_data->y, sizeof(float) * num_atoms);
        if (out_z) MEMCPY(out_z, frame_data->z, sizeof(float) * num_atoms);
    }

    if (lock) {
        md_frame_cache_frame_lock_release(lock);
    }

    return result;
}
#endif

bool load_frame(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* out_header, float* out_x, float* out_y, float* out_z) {
    ASSERT(inst);
    LoadedTrajectory* loaded_traj = (LoadedTrajectory*)inst;
    ASSERT(0 <= idx && idx < (int64_t)md_trajectory_num_frames(loaded_traj->traj));

    if ((out_x || out_y || out_z) && !(out_x && out_y && out_z))  {
        MD_LOG_ERROR("One coordinate stream (x,y,z) was null, when attempting to read out coordinates");
        return false;
    }

    md_frame_data_t* frame_data;
    md_frame_cache_lock_t* lock = 0;
    bool result = true;
    bool in_cache = md_frame_cache_find_or_reserve(&loaded_traj->cache, idx, &frame_data, &lock);
    if (!in_cache) {
        //md_allocator_i* alloc = md_heap_allocator;
        //size_t frame_data_size = md_trajectory_fetch_frame_data(loaded_traj->traj, idx, 0);
        //void*  frame_data_ptr  = md_alloc(alloc, frame_data_size);
        //md_trajectory_fetch_frame_data(loaded_traj->traj, idx, frame_data_ptr);
        //result = md_trajectory_decode_frame_data(loaded_traj->traj, frame_data_ptr, frame_data_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);
        result = md_trajectory_load_frame(loaded_traj->traj, idx, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

        if (result) {
            const md_unit_cell_t* cell = &frame_data->header.unit_cell;
            const md_molecule_t* mol = loaded_traj->mol;
            float* x = frame_data->x;
            float* y = frame_data->y;
            float* z = frame_data->z;
            const size_t num_atoms = frame_data->header.num_atoms;

            // If we have a recenter target, then compute the com and apply that transformation
            if (md_array_size(loaded_traj->recenter_indices) > 0) {
                size_t count = md_array_size(loaded_traj->recenter_indices);
                const int32_t* indices = loaded_traj->recenter_indices;

                vec3_t com = {0};
                if (count == 1) {
                    const int32_t i = indices[0];
                    com = vec3_set(x[i], y[i], z[i]);
                } else {
                    com = md_util_com_compute(x, y, z, mol->atom.mass, indices, count, &mol->unit_cell);
                    md_util_pbc(&com.x, &com.y, &com.z, 0, 1, cell);
                }

                // Translate all
                const vec3_t center = cell->flags ? cell->basis * vec3_set1(0.5f) : vec3_zero();
                const vec3_t trans  = center - com;
                vec3_batch_translate_inplace(x, y, z, num_atoms, trans);
            }
        }

        //md_free(alloc, frame_data_ptr, frame_data_size);
    }

    if (result) {
        const size_t num_bytes = frame_data->header.num_atoms * sizeof(float);
        if (out_header) *out_header = frame_data->header;
        if (out_x && out_y && out_z) {
            MEMCPY(out_x, frame_data->x, num_bytes);
            MEMCPY(out_y, frame_data->y, num_bytes);
            MEMCPY(out_z, frame_data->z, num_bytes);
        }
    }

    if (lock) {
        md_frame_cache_frame_lock_release(lock);
    }

    return result;
}

md_trajectory_i* open_file(str_t filename, md_trajectory_loader_i* loader, const md_molecule_t* mol, md_allocator_i* alloc, LoadTrajectoryFlags flags) {
    ASSERT(mol);
    ASSERT(alloc);

    if (!loader) {
        str_t ext;
        if (extract_ext(&ext, filename)) {
            loader = loader_from_ext(ext);
        }
    }
    if (!loader) {
        MD_LOG_ERROR("Unsupported file extension: '%.*s'", filename.len, filename.ptr);
        return NULL;
    }

    md_trajectory_i* internal_traj = loader->create(filename, alloc, flags);
    if (!internal_traj) {
        return NULL;
    }
    
    if (md_trajectory_num_atoms(internal_traj) != mol->atom.count) {
        MD_LOG_ERROR("Trajectory is not compatible with the loaded molecule.");
        loader->destroy(internal_traj);
        return NULL;
    }

    md_trajectory_i* traj = (md_trajectory_i*)md_alloc(alloc, sizeof(md_trajectory_i));
    memset(traj, 0, sizeof(md_trajectory_i));

    LoadedTrajectory* inst = alloc_loaded_trajectory((uint64_t)traj);
    inst->mol = mol;
    inst->loader = loader;
    inst->traj = internal_traj;
    inst->cache = {0};
    inst->recenter_indices = 0;
    inst->alloc = alloc;
    
    const size_t num_traj_frames      = md_trajectory_num_frames(internal_traj);
    const size_t frame_cache_size     = CLAMP(MEGABYTES(VIAMD_FRAME_CACHE_SIZE), MEGABYTES(4), md_os_physical_ram() / 4);
    const size_t approx_frame_size    = mol->atom.count * 3 * sizeof(float);
    const size_t max_num_cache_frames = frame_cache_size / approx_frame_size;

    const size_t  num_cache_frames    = MIN(num_traj_frames, max_num_cache_frames);
    
    MD_LOG_DEBUG("Initializing frame cache with %i frames.", (int)num_cache_frames);
    md_frame_cache_init(&inst->cache, inst->traj, alloc, num_cache_frames);

    // We only overload load frame and decode frame data to apply PBC upon loading data
    traj->inst = (md_trajectory_o*)inst;
    traj->get_header = get_header;
    traj->load_frame = load_frame;
    //traj->fetch_frame_data = 0;
    //traj->decode_frame_data = 0;
    
    return traj;
}

bool close(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        remove_loaded_trajectory(loaded_traj->key);
        memset(traj, 0, sizeof(md_trajectory_i));
        return true;
    }
    MD_LOG_ERROR("Attempting to free trajectory which was not loaded with loader");
    ASSERT(false);
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
    MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    return false;
}

bool clear_cache(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        md_frame_cache_clear(&loaded_traj->cache);
        return true;
    }
    MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    return false;
}

size_t num_cache_frames(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = find_loaded_trajectory((uint64_t)traj);
    if (loaded_traj) {
        return md_frame_cache_num_frames(&loaded_traj->cache);
    }
    MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    return 0;
}

}  // namespace traj

}  // namespace load
