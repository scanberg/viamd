#pragma once

#include <core/md_str.h>
#include <stdint.h>

struct md_system_t;
struct md_system_state_t;

// This is a stupid dispatch wrapper to select the appropriate loaders for system and trajectories

enum LoaderFlag_ {
    LoaderFlag_None = 0,
    LoaderFlag_RequiresDialogue = 1,
    LoaderFlag_System = 2,
    LoaderFlag_Trajectory = 4,
    LoaderFlag_DisableCacheWrite = 8,
    LoaderFlag_CoarseGrained = 16,
    LoaderFlag_MM = 32, // Molecular Mechanics data
    LoaderFlag_QM = 64, // Quantum Mechanics data
    LoaderFlag_Supplemental = 128, // May contain supplemental data to existing system
};

enum LoaderType_ {
    LoaderType_Undefined = 0,
    LoaderType_PDB,
    LoaderType_GRO,
    LoaderType_XYZ,
    LoaderType_XMOL,
    LoaderType_ARC,
    LoaderType_CIF,
    LoaderType_LAMMPSDATA,
    LoaderType_LAMMPSTRJ,
    LoaderType_XTC,
    LoaderType_TRR,
    LoaderType_DCD,
#if MD_VLX
    LoaderType_VLX_H5,
#endif
    LoaderType_COUNT
};

typedef uint32_t LoaderFlags;
typedef uint32_t LoaderType;

namespace loader {

    struct LoaderState {		
		LoaderType  type  = LoaderType_Undefined;
        const void* arg   = 0;
        LoaderFlags flags = LoaderFlag_None;
	};

    // The reason here why we don't directly provide prepackaged loaders based on extensions
    // Is to get a chance to glance into the file and see if we recognize it first.
    // And perhaps there is also some arguments or options that need to be supplied for the loader.
    void init(LoaderState* state, str_t filepath, const md_system_t* sys = nullptr);

    // out_state receives the loaded coordinates. Set out_state->alloc before calling; it may
    // differ from sys->alloc (a temp arena for the state, a persistent one for the system).
    bool load(md_system_t* sys, md_system_state_t* out_state, str_t filepath, const LoaderState& state);

    // For a file that ADDS to an already loaded system rather than replacing it (LoaderFlag_
    // Supplemental, set by init above). Its data lands in the system's attribute table; the atoms
    // and the coordinates are left alone, which is why there is no out_state. Returns false for a
    // loader type with nothing supplemental to contribute.
    bool load_supplemental(md_system_t* sys, str_t filepath, const LoaderState& state);

    // To help enlist supported loader type
    str_t       type_name(LoaderType type);
    str_t       type_ext(LoaderType type);
    LoaderFlags type_flags(LoaderType type);

    LoaderType type_from_ext(str_t ext);

}  // namespace loader
