#include "loader.h"

#include <core/md_log.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_tpr.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_xyz.h>
#include <md_ase_traj.h>
#include <md_mmcif.h>
#include <task_system.h>
#include <md_lammps.h>
#include <md_dcd.h>
#include <md_util.h>
#include <md_system.h>

#include <md_molden.h>
#include <md_itp.h>
#include <md_ndx.h>
#include <md_edr.h>
#include <md_xvg.h>
#include <md_csv.h>
#include <md_amber_prmtop.h>
#include <md_amber_rst7.h>
#include <md_amber_netcdf.h>
#include <md_amber_energy.h>
#include <md_psf.h>
#include <md_namdbin.h>

#include <core/md_arena_allocator.h>

#if MD_HDF5
#include <md_vlx.h>
#include <md_trexio.h>
#include <md_h5md.h>
#endif

namespace loader {

static const md_str_t loader_name[LoaderType_COUNT] = {
        MD_STR_INIT("Undefined"),
        MD_STR_INIT("Standard Protein Data Bank (pdb)"),
        MD_STR_INIT("Gromacs Structure (gro)"),
        MD_STR_INIT("xyz (xyz)"),
        MD_STR_INIT("xyz (xmol)"),
        MD_STR_INIT("xyz (arc)"),
        MD_STR_INIT("PDBx/mmCIF (cif)"),
        MD_STR_INIT("LAMMPS (data)"),
        MD_STR_INIT("LAMMPS Trajectory (lammpstrj)"),
        MD_STR_INIT("Gromacs Compressed Trajectory (xtc)"),
        MD_STR_INIT("Gromacs Lossless Trajectory (trr)"),
        MD_STR_INIT("DCD Trajectory (dcd)"),
#if MD_HDF5
        MD_STR_LIT("VeloxChem (h5)"),
#endif
        MD_STR_LIT("Molden (molden)"),
#if MD_HDF5
        MD_STR_LIT("TREXIO (trexio)"),
        MD_STR_INIT("H5MD (h5md)"),
#endif
        MD_STR_LIT("Gromacs Topology (itp/top)"),
        MD_STR_INIT("Gromacs Run Input (tpr)"),
        MD_STR_INIT("Gromacs Index Groups (ndx)"),
        MD_STR_INIT("Gromacs Energy (edr)"),
        MD_STR_INIT("xmgrace columns (xvg)"),
        MD_STR_INIT("Comma separated values (csv)"),
        MD_STR_INIT("ASE Trajectory (traj)"),
        MD_STR_INIT("Amber Topology (prmtop/parm7)"),
        MD_STR_INIT("Amber Restart (rst7/inpcrd/ncrst)"),
        MD_STR_INIT("Amber NetCDF Trajectory (nc)"),
        MD_STR_INIT("Amber Energies (mdout/mden)"),
        MD_STR_INIT("CHARMM/NAMD Topology (psf)"),
        MD_STR_INIT("NAMD Binary Coordinates (coor)"),
};

static const md_str_t loader_ext[LoaderType_COUNT] = {
        MD_STR_INIT(""),
        MD_STR_INIT("pdb"),
        MD_STR_INIT("gro"),
        MD_STR_INIT("xyz"),
        MD_STR_INIT("xmol"),
        MD_STR_INIT("arc"),
        MD_STR_INIT("cif"),
        MD_STR_INIT("data"),
        MD_STR_INIT("lammpstrj"),
        MD_STR_INIT("xtc"),
        MD_STR_INIT("trr"),
        MD_STR_INIT("dcd"),
#if MD_HDF5
        MD_STR_LIT("h5"),
#endif
        MD_STR_LIT("molden"),
#if MD_HDF5
        MD_STR_LIT("trexio"),
        MD_STR_INIT("h5md"),
#endif
        MD_STR_LIT("itp"),
        MD_STR_INIT("tpr"),
        MD_STR_INIT("ndx"),
        MD_STR_INIT("edr"),
        MD_STR_INIT("xvg"),
        MD_STR_INIT("csv"),
        MD_STR_INIT("traj"),
        MD_STR_INIT("prmtop"),
        MD_STR_INIT("rst7"),
        MD_STR_INIT("nc"),
        MD_STR_INIT("mdout"),
        MD_STR_INIT("psf"),
        MD_STR_INIT("coor"),
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
#if MD_HDF5
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM | LoaderFlag_QM,  // Veloxchem (h5)
#endif
        LoaderFlag_System | LoaderFlag_QM,                          // Molden
#if MD_HDF5
        LoaderFlag_System | LoaderFlag_QM,                          // TREXIO
        // The system from the file, then its positions as the run (init_trajectory_data publishes the
        // structure file's own frames). Topology: GROMACS writes its bonds; a file without any has
        // them inferred, see load_data_from_file
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM | LoaderFlag_Topology,  // H5MD
#endif
        LoaderFlag_Supplemental | LoaderFlag_MM,                    // GROMACS topology
        LoaderFlag_System | LoaderFlag_MM | LoaderFlag_Topology,    // GROMACS run input
        LoaderFlag_Supplemental | LoaderFlag_MM,                    // GROMACS index groups
        LoaderFlag_Supplemental | LoaderFlag_Temporal | LoaderFlag_MM, // GROMACS energy
        LoaderFlag_Supplemental | LoaderFlag_Temporal,                 // XVG
        LoaderFlag_Supplemental | LoaderFlag_Temporal,                 // CSV
        LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM,     // ASE trajectory
        // A topology builds the system once its coordinates are there; it fills in a loaded system of its atoms
        LoaderFlag_System | LoaderFlag_MM | LoaderFlag_Topology | LoaderFlag_NeedsCoordinates,  // Amber topology
        LoaderFlag_Coordinates | LoaderFlag_MM,                        // Amber restart
        LoaderFlag_Trajectory | LoaderFlag_MM,                         // Amber NetCDF trajectory
        LoaderFlag_Supplemental | LoaderFlag_Temporal | LoaderFlag_MM, // Amber energies
        LoaderFlag_System | LoaderFlag_MM | LoaderFlag_Topology | LoaderFlag_NeedsCoordinates,  // PSF
        LoaderFlag_Coordinates | LoaderFlag_MM,                        // NAMD binary coordinates
};

// Where the extension does not tell, or is shared by several formats: Amber names its files without one (prmtop, inpcrd,
// mdout, mden), .top is a GROMACS topology or an Amber one, .out and .log are anything, .coor is binary or a PDB file
// and .nc is NetCDF of any convention. The contents decide.
static LoaderType type_from_content(md_str_t filepath, md_str_t ext, LoaderType from_ext) {
    if (md_str_empty(filepath)) return LoaderType_Undefined;
    switch (from_ext) {
    case LoaderType_ITP:
        return md_amber_prmtop_file_is_topology(filepath) ? LoaderType_AMBER_PRMTOP : LoaderType_ITP;
    case LoaderType_AMBER_NETCDF:
    case LoaderType_AMBER_RESTART:
        switch (md_amber_netcdf_file_kind(filepath)) {
        case MD_AMBER_NETCDF_TRAJECTORY: return LoaderType_AMBER_NETCDF;
        case MD_AMBER_NETCDF_RESTART:    return LoaderType_AMBER_RESTART;
        default:                         return from_ext;   // the reader says what is wrong with it
        }
    case LoaderType_NAMD_COOR:
        // binaryoutput off: NAMD writes its coordinates as a PDB file of the same name
        return md_namdbin_file_is_coor(filepath) ? LoaderType_NAMD_COOR : LoaderType_PDB;
    case LoaderType_Undefined:
        break;
    default:
        return from_ext;
    }
    if (!md_str_empty(ext) && !md_str_eq_ignore_case(ext, MD_STR_LIT("out")) && !md_str_eq_ignore_case(ext, MD_STR_LIT("log")) &&
        !md_str_eq_ignore_case(ext, MD_STR_LIT("txt")) && !md_str_eq_ignore_case(ext, MD_STR_LIT("dat"))) {
        return LoaderType_Undefined;
    }
    if (md_amber_energy_file_is_mdout(filepath) || md_amber_energy_file_is_mden(filepath)) return LoaderType_AMBER_ENERGY;
    if (!md_str_empty(ext)) return LoaderType_Undefined;
    if (md_amber_prmtop_file_is_topology(filepath)) return LoaderType_AMBER_PRMTOP;
    switch (md_amber_netcdf_file_kind(filepath)) {
    case MD_AMBER_NETCDF_TRAJECTORY: return LoaderType_AMBER_NETCDF;
    case MD_AMBER_NETCDF_RESTART:    return LoaderType_AMBER_RESTART;
    default: break;
    }
    if (md_amber_rst7_file_is_restart(filepath)) return LoaderType_AMBER_RESTART;
    if (md_psf_file_is_psf(filepath)) return LoaderType_PSF;
    return LoaderType_Undefined;
}

void init(LoaderState* state, md_str_t filepath, const md_system_t* sys) {
    MD_ASSERT(state);
	*state = { 0 };

    md_str_t ext = {0};
    md_extract_ext(&ext, filepath);
    {
        state->type = type_from_content(filepath, ext, md_str_empty(ext) ? LoaderType_Undefined : type_from_ext(ext));
        if (state->type != LoaderType_Undefined) {
            state->flags = loader_flags[state->type];

            // Perform special check if LAMMPS to see if we can identify the format
            if (state->type == LoaderType_LAMMPSDATA) {
                md_lammps_atom_format_t format = md_lammps_atom_format_from_file(filepath);
                if (format) {
                    state->arg = md_lammps_atom_format_strings()[format];
                } else {
                    MD_LOG_INFO("Could not determine LAMMPS atom format for file '" MD_STR_FMT "'", MD_STR_ARG(filepath));
                    state->flags |= LoaderFlag_RequiresDialogue;
                }
            }
#if MD_HDF5
            // .h5 is not one format. A TREXIO file has a nucleus group and a VeloxChem one does
            // not, so the extension picks the reader and the CONTENT corrects it - which is also
            // why md_trexio_file_is_trexio exists.
            if (state->type == LoaderType_TREXIO || (state->type == LoaderType_VLX_H5 && md_trexio_file_is_trexio(filepath))) {
                state->type  = LoaderType_TREXIO;
                state->flags = loader_flags[LoaderType_TREXIO];
                return;
            }
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

    MD_LOG_INFO("Could not determine loader type from file extension '" MD_STR_FMT "'", MD_STR_ARG(ext));
    state->flags |= LoaderFlag_RequiresDialogue;
}

bool load(md_system_t* out_sys, md_system_state_t* out_state, md_str_t filepath, const LoaderState& state) {
    MD_ASSERT(out_sys);
    MD_ASSERT(out_state);

    switch (state.type) {
        case LoaderType_PDB: {
            md_pdb_options_t options = MD_PDB_OPTION_NONE;
            if (state.flags & LoaderFlag_DisableCacheWrite) {
                options |= MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE;
            }
            return md_pdb_system_init_from_file(out_sys, out_state, filepath, options, task_system::executor());
        }
        case LoaderType_GRO:
            return md_gro_system_init_from_file(out_sys, out_state, filepath, task_system::executor());
        case LoaderType_TPR:
            return md_tpr_system_init_from_file(out_sys, out_state, filepath);
        case LoaderType_ASE_TRAJ:
            return md_ase_traj_system_init_from_file(out_sys, out_state, filepath);
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
            return md_mmcif_system_init_from_file(out_sys, out_state, filepath, task_system::executor());
        case LoaderType_LAMMPSDATA: {
            const char* format = state.arg;
            return md_lammps_system_init_from_file(out_sys, out_state, filepath, format);
        }
        case LoaderType_AMBER_PRMTOP:
        case LoaderType_PSF: {
            const bool ok = (state.type == LoaderType_AMBER_PRMTOP) ? md_amber_prmtop_system_init_from_file(out_sys, out_state, filepath)
                                                                   : md_psf_system_init_from_file(out_sys, out_state, filepath);
            if (!ok) return false;
            if (!state.coordinates || !state.coordinates[0]) {
                MD_LOG_ERROR("'" MD_STR_FMT "' is a topology without coordinates; load a restart, a structure file or a trajectory of its atoms with it", MD_STR_ARG(filepath));
                return false;
            }
            // The topology's own box stands where the coordinates bring none
            const md_unitcell_t topology_cell = out_state->unitcell;
            if (!load_coordinates(out_state, out_sys, md_str_from_cstr(state.coordinates))) return false;
            if (!md_system_state_has_unitcell(out_state)) out_state->unitcell = topology_cell;
            return true;
        }
        case LoaderType_AMBER_RESTART:
        case LoaderType_NAMD_COOR:
            MD_LOG_ERROR("'" MD_STR_FMT "' holds coordinates only; they are loaded for a topology of their atoms (prmtop, psf)", MD_STR_ARG(filepath));
            return false;
        case LoaderType_LAMMPSTRJ:
        case LoaderType_XTC:
        case LoaderType_TRR:
        case LoaderType_DCD:
        case LoaderType_AMBER_NETCDF:
            // A trajectory is not loaded into the system: it is published as a run (publish_run).
            MD_LOG_ERROR("'" MD_STR_FMT "' is a trajectory; it is opened with publish_run, not loaded", MD_STR_ARG(filepath));
            return false;
#if MD_HDF5
        case LoaderType_VLX_H5:
            return md_vlx_system_init_from_file(out_sys, out_state, filepath, task_system::executor());
#endif
        case LoaderType_MOLDEN:
            return md_molden_system_init_from_file(out_sys, out_state, filepath);
#if MD_HDF5
        case LoaderType_TREXIO:
            return md_trexio_system_init_from_file(out_sys, out_state, filepath);
        case LoaderType_H5MD:
            return md_h5md_system_init_from_file(out_sys, out_state, filepath);
#endif
        default:
            return false;
    }
}

bool load_supplemental(md_system_t* out_sys, md_str_t filepath, const LoaderState& state, md_str_t run) {
    MD_ASSERT(out_sys);

    switch (state.type) {
#if MD_HDF5
        case LoaderType_VLX_H5:
            return md_vlx_system_supplement_from_file(out_sys, filepath, task_system::executor());
#endif
        case LoaderType_ITP:
            return md_itp_system_supplement_from_file(out_sys, filepath);
        case LoaderType_NDX:
            return md_ndx_system_supplement_from_file(out_sys, filepath);
        case LoaderType_EDR:
            return md_edr_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_XVG:
            return md_xvg_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_CSV:
            return md_csv_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_AMBER_ENERGY:
            return md_amber_energy_system_supplement_from_file(out_sys, filepath, run);
        case LoaderType_AMBER_PRMTOP:
            return md_amber_prmtop_system_supplement_from_file(out_sys, filepath);
        case LoaderType_PSF:
            return md_psf_system_supplement_from_file(out_sys, filepath);
        default:
            return false;
    }
}

md_str_t type_name(LoaderType type) {
    if (type < LoaderType_COUNT) {
        return loader_name[type];
    }
    return loader_name[LoaderType_Undefined];
}

md_str_t type_ext(LoaderType type) {
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

bool publish_run(md_system_t* sys, md_str_t filepath, md_str_t run, uint32_t flags) {
    md_str_t ext = {0};
    md_extract_ext(&ext, filepath);
    switch (type_from_content(filepath, ext, md_str_empty(ext) ? LoaderType_Undefined : type_from_ext(ext))) {
    case LoaderType_AMBER_NETCDF: return md_amber_netcdf_system_publish_run(sys, filepath, run, flags);
    case LoaderType_XTC:       return md_xtc_system_publish_run(sys, filepath, run, flags);
    case LoaderType_TRR:       return md_trr_system_publish_run(sys, filepath, run, flags);
    case LoaderType_DCD:       return md_dcd_system_publish_run(sys, filepath, run, flags);
    case LoaderType_PDB:       return md_pdb_system_publish_run(sys, filepath, run, flags);
    case LoaderType_ASE_TRAJ:  return md_ase_traj_system_publish_run(sys, filepath, run, flags);
    case LoaderType_XYZ:
    case LoaderType_XMOL:
    case LoaderType_ARC:       return md_xyz_system_publish_run(sys, filepath, run, flags);
    case LoaderType_LAMMPSTRJ: return md_lammps_system_publish_run(sys, filepath, run, flags);
#if MD_HDF5
    case LoaderType_H5MD:      return md_h5md_system_publish_run(sys, filepath, run, flags);
#endif
    default:                   return false;
    }
}

LoaderType type_from_ext(md_str_t ext) {
    // One loader, two extensions: a .top is an .itp with a [ molecules ] section (or an Amber topology, which
    // init tells from the contents)
    if (md_str_eq_ignore_case(ext, MD_STR_LIT("top"))) {
        return LoaderType_ITP;
    }
    static const struct { const char* ext; LoaderType type; } aliases[] = {
        { "parm7", LoaderType_AMBER_PRMTOP },
        { "inpcrd", LoaderType_AMBER_RESTART }, { "restrt", LoaderType_AMBER_RESTART }, { "rst", LoaderType_AMBER_RESTART },
        { "ncrst", LoaderType_AMBER_RESTART },
        { "ncdf", LoaderType_AMBER_NETCDF }, { "netcdf", LoaderType_AMBER_NETCDF },
        { "mden", LoaderType_AMBER_ENERGY },
    };
    for (size_t i = 0; i < MD_COUNTOF(aliases); ++i) {
        if (md_str_eq_cstr_ignore_case(ext, aliases[i].ext)) return aliases[i].type;
    }
    for (size_t i = 1; i < LoaderType_COUNT; ++i) {
        if (md_str_eq_ignore_case(ext, loader_ext[i])) {
            return (LoaderType)i;
        }
    }
    return LoaderType_Undefined;
}

size_t atom_count(md_str_t filepath, const LoaderState& state) {
    switch (state.type) {
    case LoaderType_AMBER_PRMTOP:  return md_amber_prmtop_file_atom_count(filepath);
    case LoaderType_PSF:           return md_psf_file_atom_count(filepath);
    case LoaderType_AMBER_NETCDF:  return md_amber_netcdf_file_atom_count(filepath);
    case LoaderType_AMBER_RESTART: {
        const size_t n = md_amber_netcdf_file_atom_count(filepath);
        return n ? n : md_amber_rst7_file_atom_count(filepath);
    }
    case LoaderType_NAMD_COOR:     return md_namdbin_file_atom_count(filepath);
    case LoaderType_DCD:           return md_dcd_file_atom_count(filepath);
    case LoaderType_XTC:           return md_xtc_file_atom_count(filepath);
    case LoaderType_TRR:           return md_trr_file_atom_count(filepath);
    default:                       return 0;
    }
}

bool load_coordinates(md_system_state_t* out_state, md_system_t* sys, md_str_t filepath) {
    MD_ASSERT(out_state && out_state->alloc);
    MD_ASSERT(sys);
    const size_t N = sys->atom.count;
    LoaderState cs;
    init(&cs, filepath, nullptr);

    switch (cs.type) {
    case LoaderType_AMBER_RESTART:
        if (md_amber_netcdf_file_kind(filepath) == MD_AMBER_NETCDF_RESTART) {
            return md_amber_netcdf_state_init_from_file(out_state, filepath, N);
        }
        return md_amber_rst7_state_init_from_file(out_state, filepath, N);
    case LoaderType_AMBER_NETCDF:
        return md_amber_netcdf_state_init_from_file(out_state, filepath, N);
    case LoaderType_NAMD_COOR:
        return md_namdbin_state_init_from_file(out_state, filepath, N);
    default:
        break;
    }

    if ((cs.flags & LoaderFlag_System) && !(cs.flags & LoaderFlag_NeedsCoordinates)) {
        // A structure file of the same atoms: its coordinates. The atoms are the topology's, which names them.
        md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
        md_system_t tmp_sys = {};
        tmp_sys.alloc = arena;
        md_system_state_t tmp_state = {};
        tmp_state.alloc = arena;
        bool ok = load(&tmp_sys, &tmp_state, filepath, cs) && md_system_state_has_coords(&tmp_state);
        if (ok && tmp_state.num_atoms != N) {
            MD_LOG_ERROR("'" MD_STR_FMT "' holds %zu atoms, the topology %zu", MD_STR_ARG(filepath), tmp_state.num_atoms, N);
            ok = false;
        }
        ok = ok && md_system_state_copy(out_state, &tmp_state);
        md_system_free(&tmp_sys);
        md_arena_allocator_destroy(arena);
        return ok;
    }

    if (cs.flags & LoaderFlag_Trajectory) {
        // The first frame of a trajectory, read through a run of its own that goes again after
        const md_str_t run = MD_STR_LIT("run/__coordinates");
        if (!publish_run(sys, filepath, run, 0)) return false;
        const bool has_cell = md_attributes_find(&sys->attributes, MD_STR_LIT("run/__coordinates/unitcell")) != nullptr;
        const md_str_t paths[] = { MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell") };
        bool ok = md_system_state_init(out_state, N);
        md_system_extract_t* ex = ok ? md_system_extract_begin(sys, run, paths, has_cell ? 2 : 1, out_state->alloc) : nullptr;
        ok = ex && md_system_extract_frame(ex, 0, out_state);
        if (ex) md_system_extract_end(ex);
        md_attributes_remove_prefix(&sys->attributes, run);
        return ok;
    }

    MD_LOG_ERROR("'" MD_STR_FMT "' holds no coordinates that are known here", MD_STR_ARG(filepath));
    return false;
}

}  // namespace loader
