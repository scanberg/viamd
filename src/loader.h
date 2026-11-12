#pragma once

#include <core/md_str.h>
#include <md_io.h>
#include <stdint.h>

struct md_system_t;
struct md_system_state_t;

// The loaders of the formats in mdlib's table (md_io.h), for the application: what is shown, asked for and stored

// The formats and what they hold are mdlib's (md_io.h); the flags of the application's own are above them.
enum LoaderFlag_ {
    LoaderFlag_None             = 0,
    LoaderFlag_System           = MD_IO_FLAG_SYSTEM,
    LoaderFlag_Trajectory       = MD_IO_FLAG_RUN,
    LoaderFlag_MM               = MD_IO_FLAG_MM,                // Molecular Mechanics data
    LoaderFlag_QM               = MD_IO_FLAG_QM,                // Quantum Mechanics data
    LoaderFlag_Supplemental     = MD_IO_FLAG_SUPPLEMENT,        // May contain supplemental data to existing system
    LoaderFlag_Topology         = MD_IO_FLAG_TOPOLOGY,          // The format supplies a complete topology: its bonds are not to be inferred
    LoaderFlag_Temporal         = MD_IO_FLAG_TEMPORAL,          // Supplemental data sampled along a trajectory: it joins the trajectory's run, so that has to be loaded first
    LoaderFlag_NeedsCoordinates = MD_IO_FLAG_NEEDS_COORDINATES, // A topology without coordinates (Amber prmtop, PSF): the system is built when its coordinates are there
    LoaderFlag_Coordinates      = MD_IO_FLAG_COORDINATES,       // Coordinates without atoms (a restart, NAMD binary coordinates): they are for a topology
    // The application's own
    LoaderFlag_RequiresDialogue  = 1u << 16,
    LoaderFlag_DisableCacheWrite = 1u << 17,
    LoaderFlag_ExplicitType      = 1u << 18,    // The loader was chosen by the user rather than told from the file: a workspace has to store it
};

enum LoaderType_ {
    LoaderType_Undefined       = MD_IO_FORMAT_UNKNOWN,
    LoaderType_PDB             = MD_IO_FORMAT_PDB,
    LoaderType_GRO             = MD_IO_FORMAT_GRO,
    LoaderType_XYZ             = MD_IO_FORMAT_XYZ,
    LoaderType_XMOL            = MD_IO_FORMAT_XMOL,
    LoaderType_ARC             = MD_IO_FORMAT_ARC,
    LoaderType_CIF             = MD_IO_FORMAT_MMCIF,
    LoaderType_LAMMPSDATA      = MD_IO_FORMAT_LAMMPS_DATA,
    LoaderType_LAMMPSTRJ       = MD_IO_FORMAT_LAMMPS_TRJ,
    LoaderType_XTC             = MD_IO_FORMAT_XTC,
    LoaderType_TRR             = MD_IO_FORMAT_TRR,
    LoaderType_DCD             = MD_IO_FORMAT_DCD,
    LoaderType_VLX_H5          = MD_IO_FORMAT_VLX,
    LoaderType_MOLDEN          = MD_IO_FORMAT_MOLDEN,
    LoaderType_TREXIO          = MD_IO_FORMAT_TREXIO,
    LoaderType_H5MD            = MD_IO_FORMAT_H5MD,          // H5MD (GROMACS mdrun -o traj.h5md, and others): a system and its trajectory in one file
    LoaderType_ITP             = MD_IO_FORMAT_ITP,           // GROMACS topology (.itp / .top), supplements a loaded system
    LoaderType_TPR             = MD_IO_FORMAT_TPR,           // GROMACS run input, topology and coordinates
    LoaderType_NDX             = MD_IO_FORMAT_NDX,           // GROMACS index groups (.ndx), supplements a loaded system with its groups
    LoaderType_EDR             = MD_IO_FORMAT_EDR,           // GROMACS energy file, supplements a loaded trajectory
    LoaderType_XVG             = MD_IO_FORMAT_XVG,           // xmgrace columns (GROMACS analysis output), a series along the loaded trajectory
    LoaderType_CSV             = MD_IO_FORMAT_CSV,           // comma separated columns, a series along the loaded trajectory
    LoaderType_ASE_TRAJ        = MD_IO_FORMAT_ASE_TRAJ,      // ASE trajectory (.traj): a system and its trajectory in one file
    LoaderType_AMBER_PRMTOP    = MD_IO_FORMAT_AMBER_PRMTOP,  // Amber topology (.prmtop / .parm7, CHAMBER too): atoms and bonds, no coordinates
    LoaderType_AMBER_RESTART   = MD_IO_FORMAT_AMBER_RESTART, // Amber restart (.rst7 / .inpcrd / .ncrst), ASCII or NetCDF: coordinates for a topology
    LoaderType_AMBER_NETCDF    = MD_IO_FORMAT_AMBER_NETCDF,  // Amber NetCDF trajectory (.nc)
    LoaderType_AMBER_ENERGY    = MD_IO_FORMAT_AMBER_ENERGY,  // Amber mdout or mden energies, a series along the loaded trajectory
    LoaderType_PSF             = MD_IO_FORMAT_PSF,           // CHARMM / NAMD topology (.psf): atoms and bonds, no coordinates
    LoaderType_NAMD_COOR       = MD_IO_FORMAT_NAMD_COOR,     // NAMD binary coordinates (.coor): coordinates for a topology
    LoaderType_COUNT           = MD_IO_FORMAT_COUNT
};

typedef uint32_t LoaderFlags;
typedef uint32_t LoaderType;

namespace loader {

    struct LoaderState {
        LoaderType  type  = LoaderType_Undefined;
        // An argument for the loader, as text so that a workspace can store it with the file: the
        // LAMMPS atom format. Not owned: it has to outlive the load.
        const char* arg   = 0;
        // For a topology without coordinates (LoaderFlag_NeedsCoordinates): the file its coordinates are read from,
        // a restart, a structure file of the same atoms or a trajectory (its first frame). Not owned.
        const char* coordinates = 0;
        LoaderFlags flags = LoaderFlag_None;
    };

    // The reason here why we don't directly provide prepackaged loaders based on extensions
    // Is to get a chance to glance into the file and see if we recognize it first.
    // And perhaps there is also some arguments or options that need to be supplied for the loader.
    void init(LoaderState* state, md_str_t filepath, const md_system_t* sys = nullptr);

    // out_state receives the loaded coordinates. Set out_state->alloc before calling; it may
    // differ from sys->alloc (a temp arena for the state, a persistent one for the system).
    bool load(md_system_t* sys, md_system_state_t* out_state, md_str_t filepath, const LoaderState& state);

    // For a file that ADDS to an already loaded system rather than replacing it (LoaderFlag_
    // Supplemental, set by init above). Its data lands in the system's attribute table; the atoms
    // and the coordinates are left alone, which is why there is no out_state. Returns false for a
    // loader type with nothing supplemental to contribute.
    //
    // run is the "run/<name>" of the loaded trajectory. A LoaderFlag_Temporal file publishes into
    // it and fails without one; every other loader ignores it.
    bool load_supplemental(md_system_t* sys, md_str_t filepath, const LoaderState& state, md_str_t run);

    // Publish the trajectory in filepath as the run "run/<name>" (see md_run_publish in mdlib's
    // md_system.h): frame axis, cells and positions streamed from the file. filepath is the file the
    // frames are in - the structure file itself for a multi model PDB or a multi frame XYZ. False for
    // a format that does not publish runs yet, and for a file with a single frame.
    bool publish_run(md_system_t* sys, md_str_t filepath, md_str_t run, uint32_t flags);

    // A topology without coordinates is loaded (load above) with the coordinates its LoaderState names. These are the
    // pieces of that, for the application which puts the two together in whatever order they come.
    //
    // The atom count of a topology, a coordinate file or a trajectory (DCD, XTC, TRR, Amber NetCDF), as its header says;
    // 0 when the format does not say up front.
    size_t atom_count(md_str_t filepath, const LoaderState& state);

    // The coordinates for the atoms of sys from filepath: a restart, NAMD binary coordinates, a structure file of as
    // many atoms, or the first frame of a trajectory. The file must hold as many atoms as sys. out_state->alloc set.
    bool load_coordinates(md_system_state_t* out_state, md_system_t* sys, md_str_t filepath);

    // To help enlist supported loader type. type_from_ext gives Undefined for a format which is not built in (those of
    // HDF5 without it), and type_available says so of a type.
    md_str_t       type_name(LoaderType type);
    md_str_t       type_ext(LoaderType type);
    LoaderFlags type_flags(LoaderType type);

    LoaderType type_from_ext(md_str_t ext);
    bool       type_available(LoaderType type);

}  // namespace loader
