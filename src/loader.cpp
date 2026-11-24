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

void init(LoaderState* state, md_str_t filepath) {
    MD_ASSERT(state);
    *state = { 0 };

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
}

md_load_result_t load(md_system_t* sys, md_system_state_t* out_state, md_str_t filepath, const LoaderState& state, md_str_t run) {
    MD_ASSERT(sys);
    MD_ASSERT(out_state);
    const md_io_format_info_t* f = info(state.type);
    if (!f->load) {
        MD_LOG_ERROR("'" MD_STR_FMT "': no loader for its format in this build", MD_STR_ARG(filepath));
        return MD_LOAD_FAILED;
    }
    const md_load_desc_t desc = load_desc(state, run);
    return f->load(sys, out_state, filepath, &desc);
}

md_load_result_t check(const md_system_t* sys, md_str_t filepath, const LoaderState& state, md_str_t run) {
    MD_ASSERT(sys);
    const md_io_format_info_t* f = info(state.type);
    if (!f->check) return MD_LOAD_MISMATCH;
    const md_load_desc_t desc = load_desc(state, run);
    return f->check(sys, filepath, &desc);
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

}  // namespace loader
