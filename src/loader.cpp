#include "loader.h"

#include <core/md_log.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_tpr.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_dcd.h>
#include <md_util.h>

#include <md_molden.h>
#include <md_itp.h>
#include <md_edr.h>
#include <md_xvg.h>
#include <md_csv.h>

#if MD_VLX
#include <md_vlx.h>
#endif

#if MD_TREXIO
#include <md_trexio.h>
#endif

namespace loader {

static const str_t loader_name[LoaderType_COUNT] = {
        STR_INIT("Undefined"),
        STR_INIT("Standard Protein Data Bank (pdb)"),
        STR_INIT("Gromacs Structure (gro)"),
        STR_INIT("xyz (xyz)"),
        STR_INIT("xyz (xmol)"),
        STR_INIT("xyz (arc)"),
        STR_INIT("PDBx/mmCIF (cif)"),
        STR_INIT("LAMMPS (data)"),
        STR_INIT("LAMMPS Trajectory (lammpstrj)"),
        STR_INIT("Gromacs Compressed Trajectory (xtc)"),
        STR_INIT("Gromacs Lossless Trajectory (trr)"),
        STR_INIT("DCD Trajectory (dcd)"),
#if MD_VLX
        STR_LIT("VeloxChem (h5)"),
#endif
        STR_LIT("Molden (molden)"),
#if MD_TREXIO
        STR_LIT("TREXIO (trexio)"),
#endif
        STR_LIT("Gromacs Topology (itp/top)"),
        STR_INIT("Gromacs Run Input (tpr)"),
        STR_INIT("Gromacs Energy (edr)"),
        STR_INIT("xmgrace columns (xvg)"),
        STR_INIT("Comma separated values (csv)"),
};

static const str_t loader_ext[LoaderType_COUNT] = {
        STR_INIT(""),
        STR_INIT("pdb"),
        STR_INIT("gro"),
        STR_INIT("xyz"),
        STR_INIT("xmol"),
        STR_INIT("arc"),
        STR_INIT("cif"),
        STR_INIT("data"),
        STR_INIT("lammpstrj"),
        STR_INIT("xtc"),
        STR_INIT("trr"),
        STR_INIT("dcd"),
#if MD_VLX
        STR_LIT("h5"),
#endif
        STR_LIT("molden"),
#if MD_TREXIO
        STR_LIT("trexio"),
#endif
        STR_LIT("itp"),
        STR_INIT("tpr"),
        STR_INIT("edr"),
        STR_INIT("xvg"),
        STR_INIT("csv"),
};

static const LoaderFlags loader_flags[LoaderType_COUNT] = {
        LoaderFlag_None,                                            // Unknown    
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM,  // PDB
        LoaderFlag_System | LoaderFlag_MM,                          // GRO
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM,  // XYZ
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM,  // XMOL
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM,  // ARC
        LoaderFlag_System | LoaderFlag_MM,                          // CIF
        LoaderFlag_System | LoaderFlag_MM,                          // LAMMPS DATA
        LoaderFlag_Trajectory | LoaderFlag_MM,                      // LAMMPS Trajectory
        LoaderFlag_Trajectory | LoaderFlag_MM,                      // XTC
        LoaderFlag_Trajectory | LoaderFlag_MM,                      // TRR
        LoaderFlag_Trajectory | LoaderFlag_MM,                      // DCD
#if MD_VLX
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM | LoaderFlag_QM,  // Veloxchem (h5)
#endif
        LoaderFlag_System | LoaderFlag_QM,                          // Molden
#if MD_TREXIO
        LoaderFlag_System | LoaderFlag_QM,                          // TREXIO
#endif
        LoaderFlag_Supplemental | LoaderFlag_MM,                    // GROMACS topology
        LoaderFlag_System | LoaderFlag_MM | LoaderFlag_Topology,    // GROMACS run input
        LoaderFlag_Supplemental | LoaderFlag_Temporal | LoaderFlag_MM, // GROMACS energy
        LoaderFlag_Supplemental | LoaderFlag_Temporal,                 // XVG
        LoaderFlag_Supplemental | LoaderFlag_Temporal,                 // CSV
};

void init(LoaderState* state, str_t filepath, const md_system_t* sys) {
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
#if MD_TREXIO
            // .h5 is not one format. A TREXIO file has a nucleus group and a VeloxChem one does
            // not, so the extension picks the reader and the CONTENT corrects it - which is also
            // why md_trexio_file_is_trexio exists.
            if (state->type == LoaderType_TREXIO || (state->type == LoaderType_VLX_H5 && md_trexio_file_is_trexio(filepath))) {
                state->type  = LoaderType_TREXIO;
                state->flags = loader_flags[LoaderType_TREXIO];
                return;
            }
#endif
#if MD_VLX
            if (state->type == LoaderType_VLX_H5 && sys) {
                // Send check to vlx to see if we can supplement the existing system with qm data
                if (md_vlx_system_is_file_supplemental(sys, filepath)) {
                    state->flags |= LoaderFlag_Supplemental;
                }
            }
#endif
            return;
        }
    }

    // A Molden file is commonly called .molden.input or .mold, and .input is not an extension worth
    // claiming, so the one case where looking inside the file is cheaper than guessing.
    if (md_molden_file_is_molden(filepath)) {
        state->type  = LoaderType_MOLDEN;
        state->flags = loader_flags[LoaderType_MOLDEN];
        return;
    }

    MD_LOG_INFO("Could not determine loader type from file extension '" STR_FMT "'", STR_ARG(ext));
    state->flags |= LoaderFlag_RequiresDialogue;
}

bool load(md_system_t* out_sys, md_system_state_t* out_state, str_t filepath, const LoaderState& state) {
    ASSERT(out_sys);
    ASSERT(out_state);

    switch (state.type) {
        case LoaderType_PDB: {
            md_pdb_options_t options = MD_PDB_OPTION_NONE;
            if (state.flags & LoaderFlag_DisableCacheWrite) {
                options |= MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE;
            }
            return md_pdb_system_init_from_file(out_sys, out_state, filepath, options);
        }
        case LoaderType_GRO:
            return md_gro_system_init_from_file(out_sys, out_state, filepath);
        case LoaderType_TPR:
            return md_tpr_system_init_from_file(out_sys, out_state, filepath);
        case LoaderType_XYZ:
        case LoaderType_XMOL:
        case LoaderType_ARC: {
            md_xyz_options_t options = MD_XYZ_OPTION_NONE;
            if (state.flags & LoaderFlag_DisableCacheWrite) {
                options |= MD_XYZ_OPTION_DISABLE_CACHE_WRITE;
            }
            return md_xyz_system_init_from_file(out_sys, out_state, filepath, options);
        }
        case LoaderType_CIF:
            return md_mmcif_system_init_from_file(out_sys, out_state, filepath);
        case LoaderType_LAMMPSDATA: {
            const char* format = (const char*)state.arg;
            return md_lammps_system_init_from_file(out_sys, out_state, filepath, format);
        }
        case LoaderType_LAMMPSTRJ:
        case LoaderType_XTC:
        case LoaderType_TRR:
        case LoaderType_DCD:
            // A trajectory is not loaded into the system: it is published as a run (publish_run).
            MD_LOG_ERROR("'" STR_FMT "' is a trajectory; it is opened with publish_run, not loaded", STR_ARG(filepath));
            return false;
#if MD_VLX
        case LoaderType_VLX_H5:
            return md_vlx_system_init_from_file(out_sys, out_state, filepath);
#endif
        case LoaderType_MOLDEN:
            return md_molden_system_init_from_file(out_sys, out_state, filepath);
#if MD_TREXIO
        case LoaderType_TREXIO:
            return md_trexio_system_init_from_file(out_sys, out_state, filepath);
#endif
        default:
            return false;
    }
}

bool load_supplemental(md_system_t* out_sys, str_t filepath, const LoaderState& state, str_t run) {
    ASSERT(out_sys);

    switch (state.type) {
#if MD_VLX
        case LoaderType_VLX_H5:
            return md_vlx_system_supplement_from_file(out_sys, filepath);
#endif
        case LoaderType_ITP:
            return md_itp_system_supplement_from_file(out_sys, filepath);
        case LoaderType_EDR:
            return md_edr_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_XVG:
            return md_xvg_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_CSV:
            return md_csv_system_supplement_from_file(out_sys, filepath, run);
        default:
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

bool publish_run(md_system_t* sys, str_t filepath, str_t run, uint32_t flags) {
    str_t ext = {0};
    if (!extract_ext(&ext, filepath)) {
        return false;
    }
    switch (type_from_ext(ext)) {
    case LoaderType_XTC:       return md_xtc_system_publish_run(sys, filepath, run, flags);
    case LoaderType_TRR:       return md_trr_system_publish_run(sys, filepath, run, flags);
    case LoaderType_DCD:       return md_dcd_system_publish_run(sys, filepath, run, flags);
    case LoaderType_PDB:       return md_pdb_system_publish_run(sys, filepath, run, flags);
    case LoaderType_XYZ:
    case LoaderType_XMOL:
    case LoaderType_ARC:       return md_xyz_system_publish_run(sys, filepath, run, flags);
    case LoaderType_LAMMPSTRJ: return md_lammps_system_publish_run(sys, filepath, run, flags);
    default:                   return false;
    }
}

LoaderType type_from_ext(str_t ext) {
    // One loader, two extensions: a .top is an .itp with a [ molecules ] section
    if (str_eq_ignore_case(ext, STR_LIT("top"))) {
        return LoaderType_ITP;
    }
    for (size_t i = 1; i < LoaderType_COUNT; ++i) {
        if (str_eq_ignore_case(ext, loader_ext[i])) {
            return (LoaderType)i;
        }
    }
    return LoaderType_Undefined;
}

}  // namespace load
