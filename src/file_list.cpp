#include <file_list.h>
#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>

FileRole file_role(LoaderFlags flags) {
    if (flags & LoaderFlag_Supplemental) return (flags & LoaderFlag_Temporal) ? FileRole_Series : FileRole_Supplement;
    if (flags & LoaderFlag_System)       return FileRole_System;
    if (flags & LoaderFlag_Trajectory)   return FileRole_Trajectory;
    return FileRole_None;
}

void file_list_free(FileList* list, md_allocator_i* alloc) {
    ASSERT(list);
    md_array_free(list->entries, alloc);
    md_array_free(list->runs, alloc);
    *list = {};
}

void file_list_clear(FileList* list) {
    ASSERT(list);
    md_array_shrink(list->entries, 0);
    md_array_shrink(list->runs, 0);
    list->atom_count = 0;
}

FileRole file_list_role(const FileList* list, size_t i) {
    ASSERT(list);
    if (i >= md_array_size(list->entries)) return FileRole_None;
    const int run = list->entries[i].run;
    for (size_t j = 0; j < i; ++j) {
        if (list->entries[j].run == run) {
            // Not the first of its place
            return run < 0 ? FileRole_Supplement : FileRole_Series;
        }
    }
    return run < 0 ? FileRole_System : FileRole_Trajectory;
}

static bool same_file(const FileEntry& a, const FileEntry& b) {
    return a.type == b.type && strcmp(a.path, b.path) == 0 && strcmp(a.arg, b.arg) == 0;
}

bool file_list_record(FileList* list, FileEntry entry, FileRole role, md_allocator_i* alloc) {
    ASSERT(list);
    ASSERT(alloc);

    switch (role) {
    case FileRole_System:
        file_list_clear(list);
        entry.run = -1;
        break;
    case FileRole_Supplement:
        entry.run = -1;
        break;
    case FileRole_Trajectory:
        md_array_push(list->runs, RunFiles{}, alloc);
        entry.run = (int)md_array_size(list->runs) - 1;
        break;
    case FileRole_Series:
        if (md_array_size(list->runs) == 0) return false;
        entry.run = (int)md_array_size(list->runs) - 1;
        break;
    default:
        return false;
    }

    if (role == FileRole_Supplement || role == FileRole_Series) {
        for (size_t i = 0; i < md_array_size(list->entries); ++i) {
            FileEntry& e = list->entries[i];
            if (e.run == entry.run && same_file(e, entry)) {
                // Loaded again: it is where it was, and loaded now if it was missing before
                e.missing &= entry.missing;
                return true;
            }
        }
    }

    md_array_push(list->entries, entry, alloc);
    return true;
}

void file_list_drop_runs(FileList* list) {
    ASSERT(list);
    size_t n = 0;
    for (size_t i = 0; i < md_array_size(list->entries); ++i) {
        if (list->entries[i].run < 0) {
            list->entries[n++] = list->entries[i];
        }
    }
    md_array_shrink(list->entries, n);
    md_array_shrink(list->runs, 0);
}

const FileEntry* file_list_system_file(const FileList* list) {
    ASSERT(list);
    for (size_t i = 0; i < md_array_size(list->entries); ++i) {
        if (list->entries[i].run < 0) return &list->entries[i];
    }
    return NULL;
}

const FileEntry* file_list_run_file(const FileList* list) {
    ASSERT(list);
    const int last_run = (int)md_array_size(list->runs) - 1;
    if (last_run >= 0) {
        for (size_t i = 0; i < md_array_size(list->entries); ++i) {
            if (list->entries[i].run == last_run) return &list->entries[i];
        }
    }
    return file_list_system_file(list);
}

// ## Workspace

static void write_input(viamd::serialization_state_t& state, const FileEntry& e, md_allocator_i* alloc) {
    const str_t path = str_from_cstr(e.path);
    str_t rel = md_path_make_relative(state.filename, path, alloc);
    if (str_empty(rel)) {
        rel = path;
    }
    viamd::write_str(state, STR_LIT("Path"), rel);
    if (e.explicit_type) {
        // By the extension the loader is known by: a name that is the same in every build, which the enum is not
        viamd::write_str(state, STR_LIT("Loader"), loader::type_ext(e.type));
    }
    if (e.arg[0] != '\0') {
        viamd::write_str(state, STR_LIT("Argument"), str_from_cstr(e.arg));
    }
}

int file_list_serialize(viamd::serialization_state_t& state, const FileList* list, int system_id, int first_run_id) {
    ASSERT(list);
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_allocator_i* alloc = md_temp_allocator(temp);

    viamd::write_section_header(state, STR_LIT("System"));
    viamd::write_int(state, STR_LIT("Id"), system_id);
    if (list->atom_count > 0) {
        viamd::write_int(state, STR_LIT("AtomCount"), (int64_t)list->atom_count);
    }

    const size_t num_runs = md_array_size(list->runs);
    for (size_t r = 0; r < num_runs; ++r) {
        viamd::write_section_header(state, STR_LIT("Run"));
        viamd::write_int(state, STR_LIT("Id"), first_run_id + (int)r);
        viamd::write_int(state, STR_LIT("System"), system_id);
        if (list->runs[r].frame_count > 0) {
            viamd::write_int(state, STR_LIT("FrameCount"), (int64_t)list->runs[r].frame_count);
        }
    }

    for (size_t i = 0; i < md_array_size(list->entries); ++i) {
        const FileEntry& e = list->entries[i];
        viamd::write_section_header(state, STR_LIT("Input"));
        if (e.run < 0) {
            viamd::write_int(state, STR_LIT("System"), system_id);
        } else {
            viamd::write_int(state, STR_LIT("Run"), first_run_id + e.run);
        }
        write_input(state, e, alloc);
    }

    return first_run_id + (int)num_runs;
}

// A file named in the workspace: relative to it, unless it was absolute when stored. Kept when it is not
// there: a missing file is reported, and saved again, rather than forgotten.
static bool resolve_path(char* buf, size_t cap, str_t folder, str_t arg) {
    str_t file;
    viamd::extract_str(file, arg);
    if (str_empty(file)) return false;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_str_writer_t joined = md_str_writer_create(alloc);
    if (!md_path_is_absolute(file)) {
        joined += folder;
    }
    joined += file;
    const str_t joined_str = md_str_writer_to_str(&joined);
    const str_t canonical  = md_path_make_canonical(joined_str, alloc);
    str_copy_to_char_buf(buf, cap, str_empty(canonical) ? joined_str : canonical);
    return true;
}

static bool extract_count(size_t& out, str_t arg) {
    int64_t v;
    if (viamd::extract_int64(v, arg) && v >= 0) {
        out = (size_t)v;
        return true;
    }
    return false;
}

// -1 is an id not given. The legacy [Files] section gets one which no version 2 section can carry (those are
// read as non negative, see workspace_files_read).
static const int NO_ID     = -1;
static const int LEGACY_ID = -2;
static inline bool given(int id) { return id != NO_ID; }

static void read_legacy(WorkspaceFiles* files, viamd::deserialization_state_t& state, str_t folder, md_allocator_i* alloc) {
    // The keys may come in any order; the files load in this one
    str_t molecule = {}, trajectory = {}, energy = {};
    md_array(str_t) series = 0;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    defer { md_temp_end(temp); };

    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        str_t file;
        viamd::extract_str(file, arg);
        if (str_empty(file)) continue;
        if      (str_eq_cstr(ident, "MoleculeFile"))   molecule   = arg;
        else if (str_eq_cstr(ident, "TrajectoryFile")) trajectory = arg;
        else if (str_eq_cstr(ident, "EnergyFile"))     energy     = arg;
        else if (str_eq_cstr(ident, "SeriesFile"))     md_array_push(series, arg, md_temp_allocator(temp));
        // CoarseGrained, which older workspaces have, is not read: what is coarse grained is inferred
    }
    if (str_empty(molecule)) return;

    files->legacy = true;
    WorkspaceFiles::System sys = { LEGACY_ID, 0 };
    md_array_push(files->systems, sys, alloc);
    WorkspaceFiles::Run run = { LEGACY_ID, LEGACY_ID, 0 };
    md_array_push(files->runs, run, alloc);

    auto add = [&](str_t a, bool of_run, FileRole role, bool opens_own_run) {
        WorkspaceFiles::Input in = {};
        in.system = of_run ? NO_ID : LEGACY_ID;
        in.run    = of_run ? LEGACY_ID : NO_ID;
        in.role   = role;
        in.opens_own_run = opens_own_run;
        if (resolve_path(in.entry.path, sizeof(in.entry.path), folder, a)) {
            md_array_push(files->inputs, in, alloc);
        }
    };
    // Without a trajectory of its own, the run was the structure file's own frames (a multi model PDB)
    add(molecule, false, FileRole_System, str_empty(trajectory));
    if (!str_empty(trajectory)) add(trajectory, true, FileRole_Trajectory, false);
    if (!str_empty(energy))     add(energy, true, FileRole_Series, false);
    for (size_t i = 0; i < md_array_size(series); ++i) {
        add(series[i], true, FileRole_Series, false);
    }
}

bool workspace_files_read(WorkspaceFiles* files, viamd::deserialization_state_t& state, str_t folder, md_allocator_i* alloc) {
    ASSERT(files);
    const str_t section = viamd::section_header(state);
    str_t ident, arg;

    if (str_eq_cstr(section, "System")) {
        WorkspaceFiles::System sys = { NO_ID, 0 };
        while (viamd::next_entry(ident, arg, state)) {
            if      (str_eq_cstr(ident, "Id"))        viamd::extract_int(sys.id, arg);
            else if (str_eq_cstr(ident, "AtomCount")) extract_count(sys.atom_count, arg);
        }
        if (sys.id < 0) {
            MD_LOG_INFO("Workspace: a [System] without an Id is left out");
            return true;
        }
        md_array_push(files->systems, sys, alloc);
        return true;
    }
    if (str_eq_cstr(section, "Run")) {
        WorkspaceFiles::Run run = { NO_ID, NO_ID, 0 };
        while (viamd::next_entry(ident, arg, state)) {
            if      (str_eq_cstr(ident, "Id"))         viamd::extract_int(run.id, arg);
            else if (str_eq_cstr(ident, "System"))     viamd::extract_int(run.system, arg);
            else if (str_eq_cstr(ident, "FrameCount")) extract_count(run.frame_count, arg);
        }
        if (run.id < 0 || run.system < 0) {
            MD_LOG_INFO("Workspace: a [Run] without an Id or a System is left out");
            return true;
        }
        md_array_push(files->runs, run, alloc);
        return true;
    }
    if (str_eq_cstr(section, "Input")) {
        WorkspaceFiles::Input in = {};
        in.system = NO_ID;
        in.run = NO_ID;
        bool has_path = false;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Path")) {
                has_path = resolve_path(in.entry.path, sizeof(in.entry.path), folder, arg);
            } else if (str_eq_cstr(ident, "System")) {
                if (viamd::extract_int(in.system, arg) && in.system < 0) in.system = LEGACY_ID - 1;   // given, and names nothing
            } else if (str_eq_cstr(ident, "Run")) {
                if (viamd::extract_int(in.run, arg) && in.run < 0) in.run = LEGACY_ID - 1;
            } else if (str_eq_cstr(ident, "Loader")) {
                str_t ext;
                viamd::extract_str(ext, arg);
                const LoaderType type = loader::type_from_ext(ext);
                if (type != LoaderType_Undefined) {
                    in.entry.type = type;
                    in.entry.explicit_type = true;
                } else {
                    MD_LOG_INFO("Workspace: no loader '" STR_FMT "' in this build, the file is read with the one its extension picks", STR_ARG(ext));
                }
            } else if (str_eq_cstr(ident, "Argument")) {
                viamd::extract_to_char_buf(in.entry.arg, sizeof(in.entry.arg), arg);
            }
        }
        if (has_path) {
            md_array_push(files->inputs, in, alloc);
        }
        return true;
    }
    if (str_eq_cstr(section, "Files") || str_eq_cstr(section, "File")) {
        // DEPRECATED: one slot per kind of file, written before version 2. "File" is an older name still.
        read_legacy(files, state, folder, alloc);
        return true;
    }
    return false;
}

void workspace_files_plan(SystemPlan* plan, const WorkspaceFiles* files, size_t max_runs, md_allocator_i* alloc) {
    ASSERT(plan);
    ASSERT(files);
    *plan = {};

    const size_t num_systems = md_array_size(files->systems);
    if (num_systems == 0) {
        plan->skipped_inputs = md_array_size(files->inputs);
        return;
    }
    const WorkspaceFiles::System& sys = files->systems[0];
    plan->atom_count = sys.atom_count;
    plan->skipped_systems = num_systems - 1;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    defer { md_temp_end(temp); };

    // The system's runs in the order the workspace names them, as many as are loaded
    md_array(int) run_ids = 0;
    for (size_t r = 0; r < md_array_size(files->runs); ++r) {
        if (files->runs[r].system != sys.id) continue;
        if (md_array_size(run_ids) < max_runs) {
            md_array_push(run_ids, files->runs[r].id, md_temp_allocator(temp));
            md_array_push(plan->frame_counts, files->runs[r].frame_count, alloc);
        } else {
            plan->skipped_runs += 1;
        }
    }
    auto run_index = [&](int id) -> int {
        for (size_t r = 0; r < md_array_size(run_ids); ++r) {
            if (run_ids[r] == id) return (int)r;
        }
        return -1;
    };
    // Owned by exactly one system or run the workspace names
    auto owned = [&](const WorkspaceFiles::Input& in) {
        if (given(in.system) == given(in.run)) return false;
        for (size_t s = 0; s < num_systems; ++s) if (given(in.system) && files->systems[s].id == in.system) return true;
        for (size_t r = 0; r < md_array_size(files->runs); ++r) if (given(in.run) && files->runs[r].id == in.run) return true;
        return false;
    };

    // The system's own files first, then each of its runs: a run is opened on the complete system
    const int num_places = 1 + (int)md_array_size(run_ids);
    for (int place = -1; place < num_places - 1; ++place) {
        bool first = true;
        for (size_t i = 0; i < md_array_size(files->inputs); ++i) {
            const WorkspaceFiles::Input& in = files->inputs[i];
            if (!owned(in)) continue;
            const bool here = place < 0 ? (given(in.system) && in.system == sys.id) : (given(in.run) && run_index(in.run) == place);
            if (!here) continue;
            SystemPlan::Input p = {};
            p.entry = in.entry;
            p.run   = place;
            p.role  = in.role != FileRole_None ? in.role : place < 0 ? (first ? FileRole_System : FileRole_Supplement) : (first ? FileRole_Trajectory : FileRole_Series);
            p.opens_own_run = in.opens_own_run;
            md_array_push(plan->inputs, p, alloc);
            first = false;
        }
    }

    // What belongs nowhere the workspace names: a dangling id, both owners or none
    for (size_t i = 0; i < md_array_size(files->inputs); ++i) {
        if (!owned(files->inputs[i])) {
            plan->skipped_inputs += 1;
        }
    }
}
