#include "loader.h"

#include <core/md_log.h>

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

namespace loader {

static const str_t loader_name[LoaderType_COUNT] = {
        STR_LIT("Undefined"),
        STR_LIT("Standard Protein Data Bank (pdb)"),
        STR_LIT("Gromacs Structure (gro)"),
        STR_LIT("xyz (xyz)"),
        STR_LIT("xyz (xmol)"),
        STR_LIT("xyz (arc)"),
        STR_LIT("PDBx/mmCIF (cif)"),
        STR_LIT("LAMMPS (data)"),
        STR_LIT("LAMMPS Trajectory (lammpstrj)"),
        STR_LIT("Gromacs Compressed Trajectory (xtc)"),
        STR_LIT("Gromacs Lossless Trajectory (trr)"),
        STR_LIT("DCD Trajectory (dcd)"),
#if MD_VLX
        STR_LIT("VeloxChem (h5)")
#endif
};

static const str_t loader_ext[LoaderType_COUNT] = {
        STR_LIT(""),
        STR_LIT("pdb"),
        STR_LIT("gro"),
        STR_LIT("xyz"),
        STR_LIT("xmol"),
        STR_LIT("arc"),
        STR_LIT("cif"),
        STR_LIT("data"),
        STR_LIT("lammpstrj"),
        STR_LIT("xtc"),
        STR_LIT("trr"),
        STR_LIT("dcd"),
#if MD_VLX
        STR_LIT("h5")
#endif
};

static const LoaderFlags loader_flags[LoaderType_COUNT] = {
        LoaderFlag_None,                            // Unknown    
        LoaderFlag_System | LoaderFlag_Trajectory,  // PDB
        LoaderFlag_System,                          // GRO
        LoaderFlag_System | LoaderFlag_Trajectory,  // XYZ
        LoaderFlag_System | LoaderFlag_Trajectory,  // XMOL
        LoaderFlag_System | LoaderFlag_Trajectory,  // ARC
        LoaderFlag_System,                          // CIF
        LoaderFlag_System,                          // LAMMPS DATA
        LoaderFlag_Trajectory,                      // LAMMPS Trajectory
        LoaderFlag_Trajectory,                      // XTC
        LoaderFlag_Trajectory,                      // TRR
        LoaderFlag_Trajectory,                      // DCD
#if MD_VLX
        LoaderFlag_System | LoaderFlag_Trajectory,  // Veloxchem (h5)
#endif
};

void init(State* state, str_t filepath) {
    ASSERT(state);
	*state = { 0 };

    str_t ext = {0};
    if (extract_ext(&ext, filepath)) {
        state->type = type_from_ext(ext);
        if (state->type != LoaderType_Undefined) {
            state->flags = loader_flags[state->type];

            // Perform special check if LAMMPS to see if we can identify the format
            if (state->type == LoaderType_LAMMPSDATA) {
                md_lammps_atom_format_t format = md_lammps_atom_format_from_file(filepath);
                if (format) {
                    state->arg = md_lammps_atom_format_strings()[format];
                } else {
                    MD_LOG_INFO("Could not determine LAMMPS atom format for file '" STR_FMT "'", STR_ARG(filepath));
                    state->flags |= LoaderFlag_RequiresDialogue;
                }
            }
            return;
        }
    }
    MD_LOG_INFO("Could not determine loader type from file extension '" STR_FMT "'", STR_ARG(ext));
    state->flags |= LoaderFlag_RequiresDialogue;
}

bool load(md_system_t* sys, str_t filepath, const State& state) {
    ASSERT(sys);    

    md_trajectory_flags_t traj_flags = MD_TRAJECTORY_FLAG_NONE;
    if (state.flags & LoaderFlag_DisableCacheWrite) {
        traj_flags |= MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE;
    }

    switch (state.type) {
        case LoaderType_PDB: {
            md_pdb_options_t options = MD_PDB_OPTION_NONE;
            if (state.flags & LoaderFlag_DisableCacheWrite) {
                options |= MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE;
            }
            return md_pdb_system_init_from_file(sys, filepath, options);
        }
        case LoaderType_GRO:
            return md_gro_system_init_from_file(sys, filepath);
        case LoaderType_XYZ:
        case LoaderType_XMOL:
        case LoaderType_ARC: {
            md_xyz_options_t options = MD_XYZ_OPTION_NONE;
            if (state.flags & LoaderFlag_DisableCacheWrite) {
                options |= MD_XYZ_OPTION_DISABLE_CACHE_WRITE;
            }
            return md_xyz_system_init_from_file(sys, filepath, options);
        }
        case LoaderType_CIF:
            return md_mmcif_system_init_from_file(sys, filepath);
        case LoaderType_LAMMPSDATA: {
            const char* format = (const char*)state.arg;
            return md_lammps_system_init_from_file(sys, filepath, format);
        }
        case LoaderType_LAMMPSTRJ:
            return md_lammps_trajectory_attach_from_file(sys, filepath, traj_flags);
        case LoaderType_XTC:
            return md_xtc_attach_from_file(sys, filepath, traj_flags);
        case LoaderType_TRR:
            return md_trr_attach_from_file(sys, filepath, traj_flags);
        case LoaderType_DCD:
            return md_dcd_attach_from_file(sys, filepath, traj_flags);
#if MD_VLX
        case LoaderType_VLX_H5:
            return md_vlx_system_init_from_file(sys, filepath);
#endif
        default:
            MD_LOG_ERROR("Undefined loader type in loader state");
            return false;
    }
}

str_t type_name(LoaderType type) {
    if (type < LoaderType_COUNT) {
        return loader_name[type];
    }
    return loader_name[LoaderType_Undefined];
}

str_t type_ext(LoaderType type) {
    if (type < LoaderType_COUNT) {
        return loader_ext[type];
    }
    return loader_ext[LoaderType_Undefined];
}

LoaderFlags type_flags(LoaderType type) {
    if (type < LoaderType_COUNT) {
        return loader_flags[type];
    }
    return loader_flags[LoaderType_Undefined];
}

LoaderType type_from_ext(str_t ext) {
    for (size_t i = 1; i < LoaderType_COUNT; ++i) {
        if (str_eq_ignore_case(ext, loader_ext[i])) {
            return (LoaderType)i;
        }
    }
    return LoaderType_Undefined;
}

}  // namespace load
