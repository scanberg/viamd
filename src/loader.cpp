#include "loader.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_bitfield.h>
#include <core/md_parse.h>
#include <core/md_lru_cache.inl>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_lammps.h>
//#include <md_dcd.h>
#include <md_trajectory.h>
#include <md_util.h>
#if MD_VLX
#include <md_vlx.h>
#endif

#include <string.h>

#include "task_system.h"

#ifndef md_lru_cache_validate
#define md_lru_cache_validate md_lru_cache8_validate
#endif

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

static md_trajectory_loader_i* traj_loader[] = {
	NULL,
	md_pdb_trajectory_loader(),
	md_xtc_trajectory_loader(),
	md_trr_trajectory_loader(),
	md_xyz_trajectory_loader(),
    md_lammps_trajectory_loader(),
};

struct LoadedTrajectory {
    const md_system_t* mol;
    md_trajectory_loader_i* loader;
    md_trajectory_i* traj;
    load::traj::FrameCache* cache;
    md_allocator_i*  alloc;
    size_t           num_cache_frames;

    md_array(int32_t) recenter_indices;
};

namespace load::traj {
bool get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header);
bool load_frame(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* out_header, float* out_x, float* out_y, float* out_z);
}

static inline LoadedTrajectory* get_loaded_trajectory(md_trajectory_i* traj) {
    if (traj && traj->inst && traj->get_header == load::traj::get_header && traj->load_frame == load::traj::load_frame) {
        return (LoadedTrajectory*)traj->inst;
    }
    return nullptr;
}

static inline void frame_cache_touch(load::traj::FrameCache* cache, int32_t slot_idx) {
    ASSERT(cache);
    md_lru_cache8_t lru = (md_lru_cache8_t)cache->lru;
    md_lru_cache8_set_mru(lru, slot_idx);
    cache->lru = (uint64_t)lru;
}

static inline int32_t frame_cache_lru_slot(const load::traj::FrameCache* cache) {
    ASSERT(cache);
    md_lru_cache8_t lru = (md_lru_cache8_t)cache->lru;
    return (int32_t)md_lru_cache8_get_lru(lru);
}

static void frame_cache_clear(load::traj::FrameCache* cache) {
    if (!cache) return;

    md_lru_cache8_t lru = 0;
    md_lru_cache8_init(lru);
    cache->lru = (uint64_t)lru;

    for (size_t i = 0; i < load::traj::FrameCacheCapacity; ++i) {
        cache->slots[i].frame_idx = -1;
        cache->slots[i].header = {};
    }
}

static void frame_cache_free(load::traj::FrameCache* cache) {
    if (!cache) return;

    if (cache->coord_buffer) {
        const size_t coord_count = load::traj::FrameCacheCapacity * cache->coord_stride * 3;
        md_free(cache->alloc, cache->coord_buffer, coord_count * sizeof(float));
    }

    *cache = {};
}

static bool frame_cache_init(load::traj::FrameCache* cache, md_trajectory_i* traj, const md_system_t* mol, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(traj);
    ASSERT(mol);
    ASSERT(alloc);

    frame_cache_free(cache);

    cache->traj = traj;
    cache->alloc = alloc;
    cache->coord_stride = ALIGN_TO(mol->atom.count, 16);
    cache->num_frames = MIN(md_trajectory_num_frames(traj), (size_t)load::traj::FrameCacheCapacity);

    const size_t coord_count = load::traj::FrameCacheCapacity * cache->coord_stride * 3;
    if (coord_count) {
        cache->coord_buffer = (float*)md_alloc(alloc, coord_count * sizeof(float));
        if (!cache->coord_buffer) {
            MD_LOG_ERROR("Failed to allocate trajectory frame cache");
            frame_cache_free(cache);
            return false;
        }
    }

    for (size_t i = 0; i < load::traj::FrameCacheCapacity; ++i) {
        float* ptr = cache->coord_buffer ? cache->coord_buffer + i * cache->coord_stride * 3 : nullptr;
        cache->slots[i].frame_idx = -1;
        cache->slots[i].header = {};
        cache->slots[i].x = ptr ? ptr + 0 * cache->coord_stride : nullptr;
        cache->slots[i].y = ptr ? ptr + 1 * cache->coord_stride : nullptr;
        cache->slots[i].z = ptr ? ptr + 2 * cache->coord_stride : nullptr;
    }

    frame_cache_clear(cache);
    return true;
}

static bool ensure_frame_cache(LoadedTrajectory* loaded_traj) {
    ASSERT(loaded_traj);

    if (!loaded_traj->cache) {
        return false;
    }

    if (loaded_traj->cache->traj != loaded_traj->traj) {
        return frame_cache_init(loaded_traj->cache, loaded_traj->traj, loaded_traj->mol, loaded_traj->alloc);
    }

    return true;
}

static void apply_recenter(LoadedTrajectory* loaded_traj, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(loaded_traj);
    ASSERT(header);
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);

    if (md_array_size(loaded_traj->recenter_indices) == 0) {
        return;
    }

    const md_unitcell_t* cell = &header->unitcell;
    const md_system_t* mol = loaded_traj->mol;
    const size_t num_atoms = header->num_atoms;
    const size_t count = md_array_size(loaded_traj->recenter_indices);
    const int32_t* indices = loaded_traj->recenter_indices;

    vec3_t com = {0};
    if (count == 1) {
        const int32_t i = indices[0];
        com = vec3_set(x[i], y[i], z[i]);
    } else {
        com = md_util_com_compute(x, y, z, NULL, indices, count, &mol->unitcell);
        md_util_pbc(&com.x, &com.y, &com.z, 0, 1, cell);
    }

    const mat3_t basis = md_unitcell_basis_mat3(cell);
    const vec3_t center = basis * vec3_set1(0.5f);
    const vec3_t trans = center - com;
    vec3_batch_translate_inplace(x, y, z, num_atoms, trans);
}

// In here each loader gets a chance to do a precheck with the file to be loaded
static void mol_loader_preload_check(load::LoaderState* state, sys_loader_t loader, str_t file_path, md_allocator_i* alloc) {
    switch (loader) {
    case SYS_LOADER_LAMMPS: {
        md_lammps_atom_format_t format = md_lammps_atom_format_from_file(file_path);
        if (format != MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
            // Encode this into the argument
            const char* format_str = md_lammps_atom_format_strings()[format];
            md_lammps_molecule_loader_arg_t arg = md_lammps_molecule_loader_arg(format_str);
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

#define NUM_ENTRIES 12
struct table_entry_t {
    str_t name[NUM_ENTRIES];
    str_t ext[NUM_ENTRIES];
    md_system_loader_i*     sys_loader[NUM_ENTRIES];
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
        //STR_LIT("dcd"),
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
            state->traj_loader = traj_loader[traj];
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
        md_allocator_i* alloc = md_get_heap_allocator();
        size_t frame_data_size = md_trajectory_fetch_frame_data(loaded_traj->traj, idx, 0);
        void*  frame_data_ptr  = md_alloc(alloc, frame_data_size);
        md_trajectory_fetch_frame_data(loaded_traj->traj, idx, frame_data_ptr);
        result = md_trajectory_decode_frame_data(loaded_traj->traj, frame_data_ptr, frame_data_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

        if (result) {
            const md_unit_cell_t* cell = &frame_data->header.unit_cell;
            const md_system_t* mol = loaded_traj->mol;
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

    if (!ensure_frame_cache(loaded_traj)) {
        md_trajectory_frame_header_t header = {};
        md_trajectory_frame_header_t* header_ptr = out_header ? out_header : ((out_x && out_y && out_z) ? &header : nullptr);
        const bool result = md_trajectory_load_frame(loaded_traj->traj, idx, header_ptr, out_x, out_y, out_z);
        if (result && out_x && out_y && out_z) {
            apply_recenter(loaded_traj, header_ptr, out_x, out_y, out_z);
        }
        return result;
    }

    load::traj::FrameCache* cache = loaded_traj->cache;
    load::traj::FrameCacheSlot* slot = nullptr;
    int32_t slot_idx = -1;

    for (int32_t i = 0; i < (int32_t)load::traj::FrameCacheCapacity; ++i) {
        if (cache->slots[i].frame_idx == idx) {
            slot = &cache->slots[i];
            slot_idx = i;
            break;
        }
    }

    bool result = true;
    if (!slot) {
        for (int32_t i = 0; i < (int32_t)load::traj::FrameCacheCapacity; ++i) {
            if (cache->slots[i].frame_idx < 0) {
                slot = &cache->slots[i];
                slot_idx = i;
                break;
            }
        }

        if (!slot) {
            slot_idx = frame_cache_lru_slot(cache);
            slot = &cache->slots[slot_idx];
        }

        result = md_trajectory_load_frame(loaded_traj->traj, idx, &slot->header, slot->x, slot->y, slot->z);
        if (result) {
            slot->frame_idx = idx;
            apply_recenter(loaded_traj, &slot->header, slot->x, slot->y, slot->z);
        } else {
            slot->frame_idx = -1;
            slot->header = {};
        }
    }

    if (result) {
        frame_cache_touch(cache, slot_idx);
        const size_t num_bytes = slot->header.num_atoms * sizeof(float);
        if (out_header) *out_header = slot->header;
        if (out_x && out_y && out_z) {
            MEMCPY(out_x, slot->x, num_bytes);
            MEMCPY(out_y, slot->y, num_bytes);
            MEMCPY(out_z, slot->z, num_bytes);
        }
    }

    return result;
}

md_trajectory_i* open_file(str_t filename, md_trajectory_loader_i* loader, const md_system_t* mol, md_allocator_i* alloc, load::traj::FrameCache* frame_cache, LoadTrajectoryFlags flags) {
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
    
    size_t traj_atom_count = md_trajectory_num_atoms(internal_traj);
    if (traj_atom_count != mol->atom.count) {
        MD_LOG_ERROR("Trajectory is not compatible with the loaded molecule.");
        loader->destroy(internal_traj);
        return NULL;
    }

    md_trajectory_i* traj = (md_trajectory_i*)md_alloc(alloc, sizeof(md_trajectory_i));
    MEMSET(traj, 0, sizeof(md_trajectory_i));

    LoadedTrajectory* inst = (LoadedTrajectory*)md_alloc(alloc, sizeof(LoadedTrajectory));
    MEMSET(inst, 0, sizeof(LoadedTrajectory));
    inst->mol = mol;
    inst->loader = loader;
    inst->traj = internal_traj;
    inst->cache = frame_cache;
    inst->recenter_indices = 0;
    inst->alloc = alloc;
    inst->num_cache_frames = MIN(md_trajectory_num_frames(internal_traj), (size_t)load::traj::FrameCacheCapacity);

    // We only overload load frame and decode frame data to apply PBC upon loading data
    traj->inst = (md_trajectory_o*)inst;
    traj->get_header = get_header;
    traj->load_frame = load_frame;
    
    return traj;
}

bool close(md_trajectory_i* traj) {
    if (traj) {
        LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
        if (loaded_traj) {
            frame_cache_free(loaded_traj->cache);
            md_array_free(loaded_traj->recenter_indices, loaded_traj->alloc);
            loaded_traj->loader->destroy(loaded_traj->traj);
            md_free(loaded_traj->alloc, loaded_traj, sizeof(LoadedTrajectory));
            MEMSET(traj, 0, sizeof(md_trajectory_i));
            return true;
        }
        MD_LOG_ERROR("Attempting to free trajectory which was not loaded with loader");
    }
    return false;
}

md_trajectory_i* get_raw_trajectory(md_trajectory_i* traj) {
    if (traj) {
        LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
        if (loaded_traj) {
            return loaded_traj->traj;
        }
        MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    }
	return nullptr;
}

bool has_recenter_target(md_trajectory_i* traj) {
    LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
    if (loaded_traj) {
        return md_array_size(loaded_traj->recenter_indices) > 0;
    }
    return false;
}

bool set_recenter_target(md_trajectory_i* traj, const md_bitfield_t* atom_mask) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
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

    LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
    if (loaded_traj) {
        frame_cache_clear(loaded_traj->cache);
        return true;
    }
    MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    return false;
}

size_t num_cache_frames(md_trajectory_i* traj) {
    ASSERT(traj);

    LoadedTrajectory* loaded_traj = get_loaded_trajectory(traj);
    if (loaded_traj) {
        return loaded_traj->cache ? loaded_traj->num_cache_frames : 0;
    }
    MD_LOG_ERROR("Supplied trajectory was not loaded with loader");
    return 0;
}

}  // namespace traj

}  // namespace load
