#include <file_placement.h>

// See file_placement.h for the rules. The passes below are in the order the header gives them.

namespace {

struct Context {
    PlacementPlan*         plan;
    const PlacementFile*   files;
    size_t                 num_files;
    const PlacementLoaded* loaded;
    size_t                 num_loaded;
    PlacementQuery         query;

    // Of each system of the plan: made new on purpose (PLACEMENT_NEW), and the system of the plan it is a copy of
    // (a file of the drop which replaces what another file of it gives), -1 for none
    bool explicit_new[PLACEMENT_MAX_FILES];
    int  copy_of[PLACEMENT_MAX_FILES];
};

PlacementRef loaded_ref(uint32_t system, int32_t run = -1) {
    PlacementRef r;
    r.system = system;
    r.run = run;
    return r;
}

PlacementRef file_ref(size_t file) {
    PlacementRef r;
    r.file = (int32_t)file;
    return r;
}

PlacementMatch match(const Context& c, size_t file, PlacementRef onto) {
    return c.query.match ? c.query.match(c.query.user, file, onto) : PlacementMatch_No;
}

bool is_system(PlacementKind k)     { return k == PlacementKind_System || k == PlacementKind_Topology; }
bool is_standalone(PlacementKind k) { return is_system(k) || k == PlacementKind_Trajectory; }
bool is_loaded_id(uint32_t t)       { return t != PLACEMENT_AUTO && t != PLACEMENT_NEW && t != PLACEMENT_NONE; }

void refuse(Context& c, size_t file, PlacementReason reason) {
    PlacementRefusal& r = c.plan->refused[c.plan->num_refused++];
    r.file = file;
    r.reason = reason;
}

size_t add_system(Context& c, PlacementAction action, uint32_t system) {
    const size_t i = c.plan->num_systems++;
    PlacementSystem& s = c.plan->systems[i];
    s = PlacementSystem{};
    s.action = action;
    s.system = system;
    c.explicit_new[i] = false;
    c.copy_of[i] = -1;
    return i;
}

void add_entry(Context& c, size_t sys, size_t file, PlacementRef run = PlacementRef{}) {
    PlacementSystem& s = c.plan->systems[sys];
    if (s.num_entries == PLACEMENT_MAX_FILES) return;
    PlacementEntry& e = s.entries[s.num_entries++];
    e.file = file;
    e.run = run;
}

void add_candidate(PlacementSystem& s, PlacementRef ref) {
    if (s.num_candidates < PLACEMENT_MAX_CANDIDATES) s.candidates[s.num_candidates++] = ref;
}

// The system of the plan which joins the loaded system, made when there is none
size_t joined(Context& c, uint32_t system) {
    for (size_t i = 0; i < c.plan->num_systems; ++i) {
        const PlacementSystem& s = c.plan->systems[i];
        if (s.action == PlacementAction_Join && s.system == system) return i;
    }
    return add_system(c, PlacementAction_Join, system);
}

// A system the drop makes: not one it joins to a loaded system on purpose, nor a choice of one file
bool made_by_drop(const Context& c, size_t sys) {
    const PlacementSystem& s = c.plan->systems[sys];
    return s.action == PlacementAction_New || s.action == PlacementAction_Copy || (s.action == PlacementAction_Ask && c.copy_of[sys] < 0 && s.num_entries && is_standalone(c.files[s.entries[0].file].kind));
}

// What a file of the drop is to the system another part of the drop makes: matched against each of its files which
// are systems of their own, or its trajectories where it has none
PlacementMatch match_drop_system(const Context& c, size_t file, size_t sys) {
    const PlacementSystem& s = c.plan->systems[sys];
    bool any_system = false;
    for (size_t i = 0; i < s.num_entries; ++i) any_system |= is_system(c.files[s.entries[i].file].kind);
    PlacementMatch result = PlacementMatch_Join;
    bool compared = false;
    for (size_t i = 0; i < s.num_entries; ++i) {
        const size_t other = s.entries[i].file;
        const PlacementKind k = c.files[other].kind;
        if (any_system ? !is_system(k) : k != PlacementKind_Trajectory) continue;
        compared = true;
        const PlacementMatch m = match(c, file, file_ref(other));
        if (m == PlacementMatch_No) return PlacementMatch_No;
        if (m == PlacementMatch_Override) result = PlacementMatch_Override;
    }
    return compared ? result : PlacementMatch_No;
}

// What the system a part of the drop makes is to a loaded one: each of its files which stand alone matched against it
PlacementMatch match_loaded(const Context& c, size_t sys, uint32_t system) {
    const PlacementSystem& s = c.plan->systems[sys];
    PlacementMatch result = PlacementMatch_Join;
    bool compared = false;
    for (size_t i = 0; i < s.num_entries; ++i) {
        const size_t file = s.entries[i].file;
        if (!is_standalone(c.files[file].kind)) continue;
        compared = true;
        const PlacementMatch m = match(c, file, loaded_ref(system));
        if (m == PlacementMatch_No) return PlacementMatch_No;
        if (m == PlacementMatch_Override) result = PlacementMatch_Override;
    }
    return compared ? result : PlacementMatch_No;
}

PlacementRef run_of(size_t file, PlacementKind kind) {
    return kind == PlacementKind_Trajectory ? file_ref(file) : PlacementRef{};
}

// ## Pass 1: the files which are systems of their own make the drop's systems, among themselves
void place_standalone(Context& c, size_t file) {
    const PlacementFile& f = c.files[file];
    if (is_loaded_id(f.target)) {
        if (match(c, file, loaded_ref(f.target)) == PlacementMatch_No) {
            refuse(c, file, PlacementReason_NotOfTarget);
            return;
        }
        add_entry(c, joined(c, f.target), file, run_of(file, f.kind));
        return;
    }
    const bool explicit_new = f.target == PLACEMENT_NEW;

    size_t joins[PLACEMENT_MAX_FILES], num_joins = 0;
    size_t overrides[PLACEMENT_MAX_FILES], num_overrides = 0;
    for (size_t i = 0; i < c.plan->num_systems; ++i) {
        if (!made_by_drop(c, i) || c.explicit_new[i] != explicit_new || c.copy_of[i] >= 0) continue;
        switch (match_drop_system(c, file, i)) {
        case PlacementMatch_Join:     joins[num_joins++] = i; break;
        case PlacementMatch_Override: overrides[num_overrides++] = i; break;
        default: break;
        }
    }

    if (num_joins == 1) {
        add_entry(c, joins[0], file, run_of(file, f.kind));
        return;
    }
    size_t sys;
    if (num_joins > 1) {
        sys = add_system(c, PlacementAction_Ask, PLACEMENT_NONE);
        for (size_t i = 0; i < num_joins; ++i) add_candidate(c.plan->systems[sys], file_ref(c.plan->systems[joins[i]].entries[0].file));
    } else {
        sys = add_system(c, PlacementAction_New, PLACEMENT_NONE);
        // Of the drop's own system, but replacing what it gives: a copy of it, the file applied after
        if (num_overrides) c.copy_of[sys] = (int)overrides[0];
    }
    c.explicit_new[sys] = explicit_new;
    add_entry(c, sys, file, run_of(file, f.kind));
}

// The systems the drop makes, against the loaded systems
void place_drop_systems(Context& c) {
    const size_t num = c.plan->num_systems;
    for (size_t i = 0; i < num; ++i) {
        PlacementSystem& s = c.plan->systems[i];
        if (s.action != PlacementAction_New || c.explicit_new[i] || c.copy_of[i] >= 0) continue;
        uint32_t joins[PLACEMENT_MAX_CANDIDATES], overrides[PLACEMENT_MAX_CANDIDATES];
        size_t num_joins = 0, num_overrides = 0;
        for (size_t l = 0; l < c.num_loaded; ++l) {
            const uint32_t id = c.loaded[l].id;
            switch (match_loaded(c, i, id)) {
            case PlacementMatch_Join:     if (num_joins < PLACEMENT_MAX_CANDIDATES) joins[num_joins++] = id; break;
            case PlacementMatch_Override: if (num_overrides < PLACEMENT_MAX_CANDIDATES) overrides[num_overrides++] = id; break;
            default: break;
            }
        }
        if (num_joins == 1) {
            // Into the system of the plan which joins it already, if there is one
            size_t into = PLACEMENT_MAX_FILES;
            for (size_t j = 0; j < c.plan->num_systems; ++j) {
                if (j != i && c.plan->systems[j].action == PlacementAction_Join && c.plan->systems[j].system == joins[0]) into = j;
            }
            s.action = PlacementAction_Join;
            s.system = joins[0];
            if (into != PLACEMENT_MAX_FILES) {
                // Emptied, and dropped when the plan is finished; still the join a copy of it resolves to
                for (size_t e = 0; e < s.num_entries; ++e) add_entry(c, into, s.entries[e].file, s.entries[e].run);
                s.num_entries = 0;
            }
        } else if (num_joins > 1) {
            s.action = PlacementAction_Ask;
            for (size_t k = 0; k < num_joins; ++k) add_candidate(s, loaded_ref(joins[k]));
        } else if (num_overrides == 1) {
            s.action = PlacementAction_Copy;
            s.system = overrides[0];
        }
        // Several it would replace something of: which to copy is not for placement to guess, and nothing of theirs
        // is replaced in a system of its own
    }
}

// A copy of a system the drop makes: of the loaded system that one joined, or of its files where it is new
void resolve_drop_copies(Context& c) {
    for (size_t i = 0; i < c.plan->num_systems; ++i) {
        const int of = c.copy_of[i];
        if (of < 0) continue;
        PlacementSystem& s = c.plan->systems[i];
        const PlacementSystem& base = c.plan->systems[of];
        if (base.action == PlacementAction_Join || base.action == PlacementAction_Copy) {
            // The plan's systems load in order, so the copy is made of what base joined or copied
            s.action = PlacementAction_Copy;
            s.system = base.system;
            if (base.action == PlacementAction_Copy) {
                // Of a copy, which is not loaded yet: the same copy, with base's files before its own
                PlacementEntry own[PLACEMENT_MAX_FILES];
                const size_t num_own = s.num_entries;
                for (size_t e = 0; e < num_own; ++e) own[e] = s.entries[e];
                s.num_entries = 0;
                for (size_t e = 0; e < base.num_entries; ++e) {
                    if (is_system(c.files[base.entries[e].file].kind)) add_entry(c, i, base.entries[e].file);
                }
                for (size_t e = 0; e < num_own; ++e) add_entry(c, i, own[e].file, own[e].run);
            }
        } else {
            // A new system, or one the user is asked about: the same files, the file applied after
            PlacementEntry own[PLACEMENT_MAX_FILES];
            const size_t num_own = s.num_entries;
            for (size_t e = 0; e < num_own; ++e) own[e] = s.entries[e];
            s.num_entries = 0;
            s.action = PlacementAction_New;
            for (size_t e = 0; e < base.num_entries; ++e) {
                if (is_system(c.files[base.entries[e].file].kind)) add_entry(c, i, base.entries[e].file);
            }
            for (size_t e = 0; e < num_own; ++e) add_entry(c, i, own[e].file, own[e].run);
        }
    }
}

// ## Pass 2: data of a system's atoms, onto the drop's systems or the loaded ones
void place_data(Context& c, size_t file) {
    const PlacementFile& f = c.files[file];
    if (f.target == PLACEMENT_NEW) {
        refuse(c, file, PlacementReason_NeedsSystem);
        return;
    }
    if (is_loaded_id(f.target)) {
        if (match(c, file, loaded_ref(f.target)) == PlacementMatch_No) {
            refuse(c, file, PlacementReason_NotOfTarget);
            return;
        }
        add_entry(c, joined(c, f.target), file);
        return;
    }
    // Candidates as systems of the plan, or loaded systems not in it yet (PLACEMENT_MAX_FILES + their index)
    size_t cands[PLACEMENT_MAX_CANDIDATES];
    size_t num = 0;
    for (size_t i = 0; i < c.plan->num_systems && num < PLACEMENT_MAX_CANDIDATES; ++i) {
        const PlacementSystem& s = c.plan->systems[i];
        PlacementMatch m = PlacementMatch_No;
        if (s.action == PlacementAction_Join) {
            m = match(c, file, loaded_ref(s.system));
        } else if (made_by_drop(c, i)) {
            m = match_drop_system(c, file, i);
        }
        // What data of the atoms brings, it adds: a match is a candidate, replacing or not
        if (m != PlacementMatch_No) cands[num++] = i;
    }
    for (size_t l = 0; l < c.num_loaded && num < PLACEMENT_MAX_CANDIDATES; ++l) {
        bool in_plan = false;
        for (size_t i = 0; i < c.plan->num_systems; ++i) {
            in_plan |= c.plan->systems[i].action == PlacementAction_Join && c.plan->systems[i].system == c.loaded[l].id;
        }
        if (!in_plan && match(c, file, loaded_ref(c.loaded[l].id)) != PlacementMatch_No) cands[num++] = PLACEMENT_MAX_FILES + l;
    }
    if (num == 0) {
        refuse(c, file, PlacementReason_NeedsSystem);
        return;
    }
    if (num == 1) {
        const size_t sys = cands[0] < PLACEMENT_MAX_FILES ? cands[0] : joined(c, c.loaded[cands[0] - PLACEMENT_MAX_FILES].id);
        add_entry(c, sys, file);
        return;
    }
    const size_t sys = add_system(c, PlacementAction_Ask, PLACEMENT_NONE);
    for (size_t k = 0; k < num; ++k) {
        if (cands[k] >= PLACEMENT_MAX_FILES) {
            add_candidate(c.plan->systems[sys], loaded_ref(c.loaded[cands[k] - PLACEMENT_MAX_FILES].id));
        } else {
            const PlacementSystem& s = c.plan->systems[cands[k]];
            add_candidate(c.plan->systems[sys], s.action == PlacementAction_Join ? loaded_ref(s.system) : file_ref(s.entries[0].file));
        }
    }
    add_entry(c, sys, file);
}

// ## Pass 3: series, onto a run the drop opens or one loaded
void place_series(Context& c, size_t file) {
    const PlacementFile& f = c.files[file];
    if (f.target == PLACEMENT_NEW) {
        refuse(c, file, PlacementReason_NeedsRun);
        return;
    }
    struct Cand { size_t sys; PlacementRef run; } cands[PLACEMENT_MAX_CANDIDATES];
    size_t num = 0;
    // The runs the drop opens
    for (size_t i = 0; i < c.plan->num_systems; ++i) {
        const PlacementSystem& s = c.plan->systems[i];
        if (is_loaded_id(f.target) && !(s.action == PlacementAction_Join && s.system == f.target)) continue;
        for (size_t e = 0; e < s.num_entries && num < PLACEMENT_MAX_CANDIDATES; ++e) {
            const size_t t = s.entries[e].file;
            if (c.files[t].kind != PlacementKind_Trajectory || s.entries[e].run.file != (int32_t)t) continue;
            if (match(c, file, file_ref(t)) != PlacementMatch_No) cands[num++] = { i, file_ref(t) };
        }
    }
    // The runs loaded (PLACEMENT_MAX_FILES + the loaded system's index: no system of the plan yet)
    for (size_t l = 0; l < c.num_loaded; ++l) {
        const uint32_t id = c.loaded[l].id;
        if (is_loaded_id(f.target) && id != f.target) continue;
        for (size_t k = 0; k < c.loaded[l].num_runs && num < PLACEMENT_MAX_CANDIDATES; ++k) {
            if (match(c, file, loaded_ref(id, (int32_t)k)) != PlacementMatch_No) cands[num++] = { PLACEMENT_MAX_FILES + l, loaded_ref(id, (int32_t)k) };
        }
    }
    if (num == 0) {
        refuse(c, file, is_loaded_id(f.target) ? PlacementReason_NotOfTarget : PlacementReason_NeedsRun);
        return;
    }
    if (num == 1) {
        const size_t sys = cands[0].sys < PLACEMENT_MAX_FILES ? cands[0].sys : joined(c, c.loaded[cands[0].sys - PLACEMENT_MAX_FILES].id);
        add_entry(c, sys, file, cands[0].run);
        return;
    }
    const size_t sys = add_system(c, PlacementAction_Ask, PLACEMENT_NONE);
    for (size_t k = 0; k < num; ++k) add_candidate(c.plan->systems[sys], cands[k].run);
    add_entry(c, sys, file);
}

// ## The cap, and what was emptied
void finish(Context& c, size_t max_systems) {
    size_t room = max_systems > c.num_loaded ? max_systems - c.num_loaded : 0;
    size_t n = 0;
    for (size_t i = 0; i < c.plan->num_systems; ++i) {
        PlacementSystem& s = c.plan->systems[i];
        if (s.num_entries == 0) continue;
        if (s.action == PlacementAction_New || s.action == PlacementAction_Copy) {
            if (room == 0) {
                for (size_t e = 0; e < s.num_entries; ++e) refuse(c, s.entries[e].file, PlacementReason_Cap);
                continue;
            }
            room -= 1;
        } else if (s.action == PlacementAction_Ask) {
            // A file of the drop, or a run, is no new system to offer; a system the drop makes is
            s.new_allowed = room > 0 && is_standalone(c.files[s.entries[0].file].kind);
        }
        if (n != i) c.plan->systems[n] = s;
        n += 1;
    }
    c.plan->num_systems = n;
}

}  // namespace

void file_placement_plan(PlacementPlan* plan, const PlacementFile* files, size_t num_files, const PlacementLoaded* loaded, size_t num_loaded, size_t max_systems, PlacementQuery query) {
    *plan = PlacementPlan{};
    Context c = {};
    c.plan = plan;
    c.files = files;
    c.num_files = num_files < PLACEMENT_MAX_FILES ? num_files : PLACEMENT_MAX_FILES;
    c.loaded = loaded;
    c.num_loaded = num_loaded;
    c.query = query;

    // Pass 1: topologies, so that their names stand, then the other systems, then trajectories
    static const PlacementKind standalone[] = { PlacementKind_Topology, PlacementKind_System, PlacementKind_Trajectory };
    for (PlacementKind kind : standalone) {
        for (size_t i = 0; i < c.num_files; ++i) {
            if (files[i].kind == kind) place_standalone(c, i);
        }
    }
    place_drop_systems(c);
    resolve_drop_copies(c);
    // Pass 2 and 3
    for (size_t i = 0; i < c.num_files; ++i) {
        if (files[i].kind == PlacementKind_Data) place_data(c, i);
    }
    for (size_t i = 0; i < c.num_files; ++i) {
        if (files[i].kind == PlacementKind_Series) place_series(c, i);
    }
    finish(c, max_systems);
    for (size_t i = c.num_files; i < num_files; ++i) {
        if (plan->num_refused < PLACEMENT_MAX_FILES) {
            plan->refused[plan->num_refused].file = i;
            plan->refused[plan->num_refused].reason = PlacementReason_TooMany;
            plan->num_refused += 1;
        }
    }
}
