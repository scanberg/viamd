#include <drop_files.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>

// The run a trajectory of the drop opens on the system it is, for the series of the drop to be matched along
static const md_str_t DROP_RUN = MD_STR_LIT("run/drop");

PlacementKind drop_file_kind(LoaderFlags flags) {
    if (flags & LoaderFlag_Temporal)                                        return PlacementKind_Series;
    if ((flags & LoaderFlag_Supplemental) && !(flags & LoaderFlag_System))  return PlacementKind_Data;
    if ((flags & LoaderFlag_Trajectory) && !(flags & LoaderFlag_System))    return PlacementKind_Trajectory;
    if ((flags & LoaderFlag_Topology) && !(flags & LoaderFlag_Trajectory))  return PlacementKind_Topology;
    return PlacementKind_System;
}

void drop_files_init(DropFiles* drop, const DropSystem* loaded, size_t num_loaded) {
    MD_ASSERT(drop);
    drop->num_files = 0;
    drop->num_skipped = 0;
    drop->loaded = loaded;
    drop->num_loaded = num_loaded;
}

void drop_files_add(DropFiles* drop, md_str_t path, const loader::LoaderState& ls, uint32_t target) {
    MD_ASSERT(drop);
    if (drop->num_files == PLACEMENT_MAX_FILES) {
        drop->num_skipped += 1;
        return;
    }
    DropFile& f = drop->files[drop->num_files++];
    f = DropFile{};
    f.entry  = file_entry_from(path, ls);
    f.ls     = ls;
    f.kind   = drop_file_kind(ls.flags);
    f.target = target;
}

static bool load_file(DropFile& f) {
    if (f.loaded) return f.loaded > 0;
    if (f.kind == PlacementKind_Data || f.kind == PlacementKind_Series) {
        f.loaded = -1;
        return false;
    }
    f.arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    f.sys.alloc = f.arena;
    f.state.alloc = f.arena;
    const bool ok = loader::load(&f.sys, &f.state, md_str_from_cstr(f.entry.path), f.ls) == MD_LOAD_OK;
    f.loaded = ok ? 1 : -1;
    return ok;
}

PlacementMatch drop_files_match(DropFiles* drop, size_t file, PlacementRef onto) {
    MD_ASSERT(drop && file < drop->num_files);
    DropFile& f = drop->files[file];
    const md_str_t path = md_str_from_cstr(f.entry.path);

    // What it is matched onto: a system of the drop (or the run its frames open), or a loaded system (or its run)
    const md_system_t* dst = nullptr;
    const md_system_state_t* dst_state = nullptr;
    md_str_t run = {};
    if (onto.file >= 0) {
        if ((size_t)onto.file >= drop->num_files || (size_t)onto.file == file) return PlacementMatch_No;
        DropFile& g = drop->files[onto.file];
        if (!load_file(g)) return PlacementMatch_No;
        dst = &g.sys;
        dst_state = &g.state;
        if (f.kind == PlacementKind_Series) {
            if (!(g.ls.flags & LoaderFlag_Trajectory)) return PlacementMatch_No;
            if (!g.run_published) {
                g.run_published = loader::publish_run(&g.sys, md_str_from_cstr(g.entry.path), DROP_RUN, 0);
                if (!g.run_published) return PlacementMatch_No;
            }
            run = DROP_RUN;
        }
    } else {
        const DropSystem* s = nullptr;
        for (size_t i = 0; i < drop->num_loaded && !s; ++i) {
            if (drop->loaded[i].id == onto.system) s = &drop->loaded[i];
        }
        if (!s || !s->sys) return PlacementMatch_No;
        dst = s->sys;
        dst_state = s->state;
        if (f.kind == PlacementKind_Series) {
            if (onto.run < 0 || (size_t)onto.run >= s->num_runs) return PlacementMatch_No;
            run = s->runs[onto.run];
        }
    }

    switch (f.kind) {
    case PlacementKind_Data:
    case PlacementKind_Series:
        return loader::check(dst, path, f.ls, run) == MD_LOAD_OK ? PlacementMatch_Join : PlacementMatch_No;
    default: {
        // A file which is data of a system's atoms as well (an itp's molecules, the QM region of a VeloxChem file): of
        // those atoms when its check says, whether or not it is all of them
        if ((f.ls.flags & LoaderFlag_Supplemental) && dst->atom.count > 0 && loader::check(dst, path, f.ls) == MD_LOAD_OK) {
            return PlacementMatch_Join;
        }
        if (!load_file(f)) return PlacementMatch_No;
        uint32_t override = 0;
        if (!md_system_match(&override, dst, dst_state, &f.sys, &f.state)) return PlacementMatch_No;
        // Frames add a run: the coordinates of a trajectory's first frame replace nothing
        if (f.kind == PlacementKind_Trajectory) override &= ~(uint32_t)MD_SYSTEM_KIND_COORDINATES;
        return override ? PlacementMatch_Override : PlacementMatch_Join;
    }
    }
}

static PlacementMatch query_match(void* user, size_t file, PlacementRef onto) {
    return drop_files_match((DropFiles*)user, file, onto);
}

void drop_files_plan(PlacementPlan* plan, DropFiles* drop, size_t max_systems) {
    MD_ASSERT(plan && drop);
    PlacementFile files[PLACEMENT_MAX_FILES];
    for (size_t i = 0; i < drop->num_files; ++i) {
        files[i].kind = drop->files[i].kind;
        files[i].target = drop->files[i].target;
    }
    PlacementLoaded loaded[PLACEMENT_MAX_CANDIDATES];
    const size_t num_loaded = MD_MIN(drop->num_loaded, MD_COUNTOF(loaded));
    for (size_t i = 0; i < num_loaded; ++i) {
        loaded[i].id = drop->loaded[i].id;
        loaded[i].num_runs = drop->loaded[i].num_runs;
    }
    const PlacementQuery query = { drop, query_match };
    file_placement_plan(plan, files, drop->num_files, loaded, num_loaded, max_systems, query);
}

const DropFile* drop_files_system(DropFiles* drop, size_t file) {
    MD_ASSERT(drop && file < drop->num_files);
    DropFile& f = drop->files[file];
    return load_file(f) ? &f : nullptr;
}

void drop_files_free(DropFiles* drop) {
    MD_ASSERT(drop);
    for (size_t i = 0; i < drop->num_files; ++i) {
        DropFile& f = drop->files[i];
        if (f.arena) {
            md_system_free(&f.sys);
            md_arena_allocator_destroy(f.arena);
        }
        f = DropFile{};
    }
    drop->num_files = 0;
}
