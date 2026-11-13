#include "loader.h"

#include <core/md_log.h>
#include <core/md_path.h>
#include <core/md_arena_allocator.h>

#include <md_io.h>
#include <md_lammps.h>
#include <md_system.h>
#include <md_run.h>
#include <md_attributes.h>
#include <task_system.h>

#if MD_HDF5
#include <md_vlx.h>
#endif

namespace loader {

static const md_io_format_info_t* info(LoaderType type) {
    return md_io_format_info((md_io_format_t)type);
}

static md_load_desc_t load_desc(const LoaderState& state, md_str_t run) {
    md_load_desc_t desc = {};
    desc.executor = task_system::executor();
    desc.format   = state.arg ? md_str_from_cstr(state.arg) : md_str_t{};
    desc.run      = run;
    return desc;
}

bool type_available(LoaderType type) {
    return md_io_format_available((md_io_format_t)type);
}

void init(LoaderState* state, md_str_t filepath, const md_system_t* sys) {
    MD_ASSERT(state);
    *state = { 0 };
    (void)sys;

    const md_io_format_t format = md_io_format_from_file(filepath);
    if (format == MD_IO_FORMAT_UNKNOWN || !md_io_format_available(format)) {
        MD_LOG_INFO("Could not determine loader type of file '" MD_STR_FMT "'", MD_STR_ARG(filepath));
        state->flags |= LoaderFlag_RequiresDialogue;
        return;
    }
    state->type  = (LoaderType)format;
    state->flags = info(state->type)->flags;

    // Perform special check if LAMMPS to see if we can identify the format
    if (state->type == LoaderType_LAMMPSDATA) {
        md_lammps_atom_format_t lammps_format = md_lammps_atom_format_from_file(filepath);
        if (lammps_format) {
            state->arg = md_lammps_atom_format_strings()[lammps_format];
        } else {
            MD_LOG_INFO("Could not determine LAMMPS atom format for file '" MD_STR_FMT "'", MD_STR_ARG(filepath));
            state->flags |= LoaderFlag_RequiresDialogue;
        }
    }
#if MD_HDF5
    // Send check to vlx to see if we can supplement the existing system with qm data
    if (state->type == LoaderType_VLX_H5 && sys && md_vlx_system_is_file_supplemental(sys, filepath)) {
        state->flags |= LoaderFlag_Supplemental;
    }
#endif
}

bool load(md_system_t* out_sys, md_system_state_t* out_state, md_str_t filepath, const LoaderState& state) {
    MD_ASSERT(out_sys);
    MD_ASSERT(out_state);

    const md_io_format_info_t* f = info(state.type);
    if (f->flags & LoaderFlag_Coordinates) {
        MD_LOG_ERROR("'" MD_STR_FMT "' holds coordinates only; they are loaded for a topology of their atoms (prmtop, psf)", MD_STR_ARG(filepath));
        return false;
    }
    if (!f->system_init) {
        if (f->publish_run) {
            // A trajectory is not loaded into the system: it is published as a run (publish_run).
            MD_LOG_ERROR("'" MD_STR_FMT "' is a trajectory; it is opened with publish_run, not loaded", MD_STR_ARG(filepath));
        }
        return false;
    }

    const md_load_desc_t desc = load_desc(state, {});
    if (!f->system_init(out_sys, out_state, filepath, &desc)) return false;

    if (f->flags & LoaderFlag_NeedsCoordinates) {
        if (!state.coordinates || !state.coordinates[0]) {
            MD_LOG_ERROR("'" MD_STR_FMT "' is a topology without coordinates; load a restart, a structure file or a trajectory of its atoms with it", MD_STR_ARG(filepath));
            return false;
        }
        // The topology's own box stands where the coordinates bring none
        const md_unitcell_t topology_cell = out_state->unitcell;
        if (!load_coordinates(out_state, out_sys, md_str_from_cstr(state.coordinates))) return false;
        if (!md_system_state_has_unitcell(out_state)) out_state->unitcell = topology_cell;
    }
    return true;
}

bool load_supplemental(md_system_t* out_sys, md_str_t filepath, const LoaderState& state, md_str_t run) {
    MD_ASSERT(out_sys);
    const md_io_format_info_t* f = info(state.type);
    if (!f->supplement) return false;
    const md_load_desc_t desc = load_desc(state, run);
    return f->supplement(out_sys, filepath, &desc);
}

md_str_t type_name(LoaderType type) {
    return info(type)->name;
}

md_str_t type_ext(LoaderType type) {
    return info(type)->ext;
}

LoaderFlags type_flags(LoaderType type) {
    return info(type)->flags;
}

bool publish_run(md_system_t* sys, md_str_t filepath, md_str_t run, uint32_t flags) {
    const md_io_format_t format = md_io_format_from_file(filepath);
    const md_io_format_info_t* f = md_io_format_info(format);
    return f->publish_run ? f->publish_run(sys, filepath, run, flags) : false;
}

LoaderType type_from_ext(md_str_t ext) {
    const md_io_format_t format = md_io_format_from_ext(ext);
    return md_io_format_available(format) ? (LoaderType)format : LoaderType_Undefined;
}

size_t atom_count(md_str_t filepath, const LoaderState& state) {
    const md_io_format_info_t* f = info(state.type);
    return f->atom_count ? f->atom_count(filepath) : 0;
}

bool load_coordinates(md_system_state_t* out_state, md_system_t* sys, md_str_t filepath) {
    MD_ASSERT(out_state && out_state->alloc);
    MD_ASSERT(sys);
    const size_t N = sys->atom.count;
    LoaderState cs;
    init(&cs, filepath, nullptr);

    // A restart, NAMD binary coordinates, or a trajectory whose reader gives its first frame
    const md_io_format_info_t* f = info(cs.type);
    if (f->state_init) {
        return f->state_init(out_state, filepath, N);
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
