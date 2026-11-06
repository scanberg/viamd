#include <load_order.h>

#include <core/md_log.h>

#include <stdio.h>
#include <initializer_list>

static bool pending_has(const FileEntry& e) { return e.path[0] != '\0'; }

static size_t count_of(const LoadOps& ops, str_t path, const loader::LoaderState& ls) {
    return ops.atom_count ? ops.atom_count(ops.user, path, ls) : loader::atom_count(path, ls);
}

loader::LoaderState load_order_loader_state(const FileEntry& e) {
    loader::LoaderState ls = {};
    if (e.explicit_type) {
        ls.type  = e.type;
        ls.flags = loader::type_flags(e.type) | LoaderFlag_ExplicitType;
    } else {
        loader::init(&ls, str_from_cstr(e.path), nullptr);
    }
    if (e.arg[0] != '\0') ls.arg = e.arg;
    return ls;
}

// The system from a topology and the file its coordinates come from. When that is a trajectory, it opens as the run.
// Files whose headers say they are of different atom counts are refused before anything is freed.
static LoadResult build_from_topology(const LoadOps& ops, const FileEntry& topology_in, const FileEntry& coordinates_in) {
    // Copies: the entries may be slots of what waits, which the load clears
    const FileEntry topology = topology_in;
    const FileEntry coordinates = coordinates_in;
    loader::LoaderState ls = load_order_loader_state(topology);
    const loader::LoaderState cs = load_order_loader_state(coordinates);
    const size_t a = count_of(ops, str_from_cstr(topology.path), ls);
    const size_t b = count_of(ops, str_from_cstr(coordinates.path), cs);
    if (a && b && a != b) {
        MD_LOG_ERROR("'%s' has %zu atoms and '%s' %zu: they are not of the same system", topology.path, a, coordinates.path, b);
        return LoadResult_Failed;
    }
    ls.coordinates = coordinates.path;
    const LoadResult result = ops.load(ops.user, str_from_cstr(topology.path), ls, FileRole_System, false);
    if (result != LoadResult_System) return result;

    if ((cs.flags & LoaderFlag_Trajectory) && !(cs.flags & LoaderFlag_System)) {
        ops.load(ops.user, str_from_cstr(coordinates.path), cs, FileRole_Trajectory, false);
    }
    return result;
}

// After a system is built: what waited for it. A trajectory opens as its run, unless the system has one; the rest is
// dropped, as it was for a system that is no longer there.
static void pending_resolve(PendingFiles* pending, const LoadOps& ops) {
    const FileEntry trajectory = pending->trajectory;
    *pending = {};
    if (pending_has(trajectory) && ops.run_frames(ops.user) == 0) {
        const loader::LoaderState ts = load_order_loader_state(trajectory);
        const size_t n = count_of(ops, str_from_cstr(trajectory.path), ts);
        if (n == 0 || n == ops.system_atoms(ops.user)) {
            ops.load(ops.user, str_from_cstr(trajectory.path), ts, FileRole_Trajectory, false);
        } else {
            MD_LOG_INFO("'%s' does not have the atoms of the loaded system and no longer waits", trajectory.path);
        }
    }
}

static LoadResult pending_wait(FileEntry* slot, const LoadOps& ops, str_t path, const loader::LoaderState& load_state, const char* what) {
    *slot = file_entry_from(path, load_state);
    const size_t n = count_of(ops, path, load_state);
    char count[64] = "";
    if (n) snprintf(count, sizeof(count), " of %zu atoms", n);
    MD_LOG_INFO("'" STR_FMT "'%s is waiting for %s", STR_ARG(path), count, what);
    return LoadResult_Pending;
}

LoadResult load_in_order(PendingFiles* pending, const LoadOps& ops, str_t filepath, const loader::LoaderState& load_state) {
    ASSERT(pending);
    ASSERT(ops.load && ops.system_atoms && ops.run_frames && ops.files);

    // On its own the role follows from the loader flags, and a structure file opens its own frames as the run:
    // a multi model PDB or an XYZ trajectory is a system and a run in one file
    const FileRole role = file_role(load_state.flags);
    const size_t num_atoms = ops.system_atoms(ops.user);

    // ## A topology without coordinates
    if ((load_state.flags & LoaderFlag_NeedsCoordinates) && !(load_state.coordinates && load_state.coordinates[0])) {
        if (load_state.flags & LoaderFlag_Supplemental) {
            return ops.load(ops.user, filepath, load_state, FileRole_Supplement, false);
        }
        const size_t n = count_of(ops, filepath, load_state);
        if (n && n == num_atoms) {
            // A system of as many atoms is loaded (a PDB file of the simulation): the topology fills it in
            loader::LoaderState ls = load_state;
            ls.flags |= LoaderFlag_Supplemental;
            return ops.load(ops.user, filepath, ls, FileRole_Supplement, false);
        }
        // Its coordinates may be waiting: those of a restart first, else the first frame of a trajectory
        const FileEntry topology = file_entry_from(filepath, load_state);
        for (FileEntry* coords : { &pending->coordinates, &pending->trajectory }) {
            if (!pending_has(*coords)) continue;
            const FileEntry c = *coords;
            const LoadResult result = build_from_topology(ops, topology, c);
            if (result == LoadResult_System) {
                // The trajectory it was built with is open as the run already
                if (coords == &pending->trajectory) pending->trajectory = {};
                pending_resolve(pending, ops);
                return result;
            }
        }
        return pending_wait(&pending->topology, ops, filepath, load_state, "its coordinates: load a restart, a structure file or a trajectory of its atoms");
    }

    // ## Coordinates without atoms
    if (load_state.flags & LoaderFlag_Coordinates) {
        const FileEntry coords = file_entry_from(filepath, load_state);
        if (pending_has(pending->topology)) {
            const LoadResult result = build_from_topology(ops, pending->topology, coords);
            if (result == LoadResult_System) {
                pending_resolve(pending, ops);
                return result;
            }
        }
        // New coordinates for a system built from a topology: built again with them, its run opened again
        const FileList* list = ops.files(ops.user);
        const FileEntry* system = file_list_system_file(list);
        if (system && system->coordinates[0] != '\0' && num_atoms) {
            const FileEntry topology = *system;
            FileEntry run = {};
            for (size_t i = 0; i < file_list_count(list); ++i) {
                if (file_list_role(list, i) == FileRole_Trajectory) run = *file_list_get(list, i);
            }
            const LoadResult result = build_from_topology(ops, topology, coords);
            if (result == LoadResult_System && pending_has(run)) {
                ops.load(ops.user, str_from_cstr(run.path), load_order_loader_state(run), FileRole_Trajectory, false);
            }
            return result;
        }
        return pending_wait(&pending->coordinates, ops, filepath, load_state, "a topology of its atoms (prmtop, psf)");
    }

    // ## A trajectory for a waiting topology, or before anything it could be the run of
    if (role == FileRole_Trajectory) {
        const FileEntry trajectory = file_entry_from(filepath, load_state);
        if (pending_has(pending->topology)) {
            const LoadResult result = build_from_topology(ops, pending->topology, trajectory);
            if (result == LoadResult_System) pending_resolve(pending, ops);
            if (result != LoadResult_Failed) return result;
        }
        // Not of the loaded system's atoms, as far as its header tells
        const size_t n = count_of(ops, filepath, load_state);
        if (num_atoms == 0 || (n && n != num_atoms)) {
            return pending_wait(&pending->trajectory, ops, filepath, load_state, "a topology or a structure file of its atoms");
        }
    }

    // A structure file for a system built from a topology: the structure builds it, the topology fills it in and the
    // run opens again - as when the structure came first
    FileEntry topology_of_system = {};
    FileEntry run_of_system = {};
    if (role == FileRole_System) {
        const FileList* list = ops.files(ops.user);
        const FileEntry* system = file_list_system_file(list);
        if (system && system->coordinates[0] != '\0' && !pending_has(pending->topology)) {
            topology_of_system = *system;
            topology_of_system.coordinates[0] = '\0';
            for (size_t i = 0; i < file_list_count(list); ++i) {
                if (file_list_role(list, i) == FileRole_Trajectory) run_of_system = *file_list_get(list, i);
            }
        }
    }

    const LoadResult result = ops.load(ops.user, filepath, load_state, role, role == FileRole_System);
    if (result == LoadResult_System) {
        if (pending_has(topology_of_system)) {
            pending->topology = topology_of_system;
            if (pending_has(run_of_system) && !pending_has(pending->trajectory)) pending->trajectory = run_of_system;
        }
        // A structure file of the waiting topology's atoms: the topology fills it in, as when it comes second
        if (pending_has(pending->topology)) {
            const FileEntry topology = pending->topology;
            pending->topology = {};
            loader::LoaderState ts = load_order_loader_state(topology);
            if (count_of(ops, str_from_cstr(topology.path), ts) == ops.system_atoms(ops.user)) {
                ts.flags |= LoaderFlag_Supplemental;
                ops.load(ops.user, str_from_cstr(topology.path), ts, FileRole_Supplement, false);
            } else {
                MD_LOG_INFO("'%s' does not have the atoms of the loaded system and no longer waits", topology.path);
            }
        }
        pending_resolve(pending, ops);
    }
    return result;
}
