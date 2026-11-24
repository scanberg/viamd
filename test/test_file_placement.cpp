#include "utest.h"

#include <file_placement.h>

#include <string>
#include <vector>

/* Where the files of one drop go (file_placement.h). Whether a file is of a system's atoms is the application's to
 * answer, from the files; here it is a table, so that each scenario says only what matches what. A plan is written
 * as text: a system per '|', its action, then its files in load order with the run each is of. */

namespace {

struct Rule {
    size_t         file;
    PlacementRef   onto;
    PlacementMatch match;
};

struct Drop {
    std::vector<PlacementFile>   files;
    std::vector<PlacementLoaded> loaded;
    std::vector<Rule>            rules;
    size_t                       max_systems = 8;
    PlacementPlan                plan;

    size_t file(PlacementKind kind, uint32_t target = PLACEMENT_AUTO) {
        PlacementFile f;
        f.kind = kind;
        f.target = target;
        files.push_back(f);
        return files.size() - 1;
    }
    void system(uint32_t id, size_t num_runs = 0) {
        PlacementLoaded l;
        l.id = id;
        l.num_runs = num_runs;
        loaded.push_back(l);
    }
    // a is of the system b builds (or the run trajectory b opens), both files of the drop
    void of_file(size_t a, size_t b, PlacementMatch m = PlacementMatch_Join) {
        PlacementRef r; r.file = (int32_t)b; rules.push_back({ a, r, m });
        // Of each other: whichever of them comes first
        if (m == PlacementMatch_Join) { PlacementRef q; q.file = (int32_t)a; rules.push_back({ b, q, m }); }
    }
    void of_system(size_t a, uint32_t system, PlacementMatch m = PlacementMatch_Join) {
        PlacementRef r; r.system = system; rules.push_back({ a, r, m });
    }
    void of_run(size_t a, uint32_t system, int32_t run) {
        PlacementRef r; r.system = system; r.run = run; rules.push_back({ a, r, PlacementMatch_Join });
    }

    static PlacementMatch query(void* user, size_t file, PlacementRef onto) {
        const Drop* d = (const Drop*)user;
        for (const Rule& r : d->rules) {
            if (r.file == file && r.onto.system == onto.system && r.onto.run == onto.run && r.onto.file == onto.file) return r.match;
        }
        return PlacementMatch_No;
    }

    std::string place() {
        PlacementQuery q;
        q.user = this;
        q.match = query;
        file_placement_plan(&plan, files.data(), files.size(), loaded.data(), loaded.size(), max_systems, q);
        return str();
    }

    static std::string ref(PlacementRef r) {
        if (r.file >= 0) return "f" + std::to_string(r.file);
        return "s" + std::to_string(r.system) + (r.run >= 0 ? "." + std::to_string(r.run) : "");
    }

    std::string str() const {
        static const char* actions[] = { "Join", "New", "Copy", "Ask" };
        static const char* reasons[] = { "NeedsSystem", "NeedsRun", "NotOfTarget", "Cap", "TooMany" };
        std::string s;
        for (size_t i = 0; i < plan.num_systems; ++i) {
            const PlacementSystem& sys = plan.systems[i];
            if (!s.empty()) s += " | ";
            s += actions[sys.action];
            if (sys.action == PlacementAction_Join || sys.action == PlacementAction_Copy) s += "(s" + std::to_string(sys.system) + ")";
            if (sys.action == PlacementAction_Ask) {
                s += "[";
                for (size_t k = 0; k < sys.num_candidates; ++k) s += (k ? " " : "") + ref(sys.candidates[k]);
                s += sys.new_allowed ? " +new]" : "]";
            }
            s += ":";
            for (size_t e = 0; e < sys.num_entries; ++e) {
                const PlacementEntry& en = sys.entries[e];
                s += " " + std::to_string(en.file);
                if (en.run.file >= 0 || en.run.system != PLACEMENT_NONE) s += "@" + ref(en.run);
            }
        }
        for (size_t i = 0; i < plan.num_refused; ++i) {
            s += " ! " + std::to_string(plan.refused[i].file) + ":" + reasons[plan.refused[i].reason];
        }
        return s;
    }
};

}  // namespace

// The plan of a drop against its text, the string held while it is compared
#define EXPECT_PLAN(text, drop) do { const std::string plan_ = (drop).place(); EXPECT_STREQ(text, plan_.c_str()); } while (0)

// A topology alone is a system without coordinates; its coordinates, dropped later, join it
UTEST(viamd_file_placement, a_topology_first_then_its_coordinates) {
    Drop first;
    first.file(PlacementKind_Topology);
    EXPECT_PLAN("New: 0", first);

    Drop then;
    then.system(1);
    const size_t dcd = then.file(PlacementKind_Trajectory);
    then.of_system(dcd, 1);
    EXPECT_PLAN("Join(s1): 0@f0", then);
}

// Two systems dropped together with their trajectories: each trajectory joins the structure of its atoms, whatever
// order the files came in
UTEST(viamd_file_placement, two_systems_with_their_trajectories) {
    Drop d;
    const size_t a_gro = d.file(PlacementKind_System);
    const size_t a_xtc = d.file(PlacementKind_Trajectory);
    const size_t b_gro = d.file(PlacementKind_System);
    const size_t b_xtc = d.file(PlacementKind_Trajectory);
    d.of_file(a_xtc, a_gro);
    d.of_file(b_xtc, b_gro);
    EXPECT_PLAN("New: 0 1@f1 | New: 2 3@f3", d);

    // With a system loaded of neither
    d.system(7);
    EXPECT_PLAN("New: 0 1@f1 | New: 2 3@f3", d);
}

// A topology, a structure of its atoms and a trajectory, together: one system, the topology loaded first so that its
// names stand
UTEST(viamd_file_placement, a_topology_its_structure_and_trajectory) {
    Drop d;
    const size_t pdb = d.file(PlacementKind_System);
    const size_t dcd = d.file(PlacementKind_Trajectory);
    const size_t psf = d.file(PlacementKind_Topology);
    d.of_file(pdb, psf);
    d.of_file(dcd, psf);
    d.of_file(dcd, pdb);
    EXPECT_PLAN("New: 2 0 1@f1", d);
}

// A structure of other atoms than the one system loaded: a system of its own, never a replacement
UTEST(viamd_file_placement, another_structure_is_another_system) {
    Drop d;
    d.system(1);
    d.file(PlacementKind_System);
    EXPECT_PLAN("New: 0", d);
}

// Of the atoms of several loaded systems: the user picks
UTEST(viamd_file_placement, several_candidates_ask) {
    Drop d;
    d.system(1);
    d.system(2);
    d.system(3);
    const size_t dcd = d.file(PlacementKind_Trajectory);
    d.of_system(dcd, 1);
    d.of_system(dcd, 3);
    EXPECT_PLAN("Ask[s1 s3 +new]: 0@f0", d);

    // No room for a new system: the choice is among the two
    d.max_systems = 3;
    EXPECT_PLAN("Ask[s1 s3]: 0@f0", d);
}

// A restart for a system which has coordinates: never replaces them, but makes a copy of the system with the restart's
UTEST(viamd_file_placement, coordinates_for_a_complete_system_copy_it) {
    Drop d;
    d.system(1);
    const size_t rst = d.file(PlacementKind_System);
    d.of_system(rst, 1, PlacementMatch_Override);
    EXPECT_PLAN("Copy(s1): 0", d);

    // Of two it would replace something of: which to copy is not guessed
    d.system(2);
    d.of_system(rst, 2, PlacementMatch_Override);
    EXPECT_PLAN("New: 0", d);

    // A candidate it replaces nothing of is chosen over a copy
    d.system(3);
    d.of_system(rst, 3);
    EXPECT_PLAN("Join(s3): 0", d);
}

// Two structures of the same atoms in one drop: two systems, the second a copy of the first with its own coordinates
UTEST(viamd_file_placement, two_snapshots_together) {
    Drop d;
    const size_t a = d.file(PlacementKind_System);
    const size_t b = d.file(PlacementKind_System);
    d.of_file(b, a, PlacementMatch_Override);
    EXPECT_PLAN("New: 0 | New: 0 1", d);

    // The first joins a loaded system it fills in: the second is a copy of that system
    d.system(5);
    d.of_system(a, 5);
    EXPECT_PLAN("Join(s5): 0 | Copy(s5): 1", d);
}

// Points first: a trajectory alone is a system of points, which a topology of as many atoms names later
UTEST(viamd_file_placement, points_first_then_a_topology) {
    Drop first;
    first.file(PlacementKind_Trajectory);
    EXPECT_PLAN("New: 0@f0", first);

    Drop then;
    then.system(1, 1);
    const size_t psf = then.file(PlacementKind_Topology);
    then.of_system(psf, 1);
    EXPECT_PLAN("Join(s1): 0", then);
}

// A series joins the run it is along: one of the drop's own, or a loaded one; of none it is refused
UTEST(viamd_file_placement, series_and_runs) {
    Drop d;
    const size_t gro = d.file(PlacementKind_System);
    const size_t xtc = d.file(PlacementKind_Trajectory);
    const size_t edr = d.file(PlacementKind_Series);
    d.of_file(xtc, gro);
    d.of_file(edr, xtc);
    EXPECT_PLAN("New: 0 1@f1 2@f1", d);

    // A loaded system of two runs, the series along the second
    Drop two;
    two.system(1, 2);
    const size_t e = two.file(PlacementKind_Series);
    two.of_run(e, 1, 1);
    EXPECT_PLAN("Join(s1): 0@s1.1", two);
    // Along both: asked, and a series is no new system
    two.of_run(e, 1, 0);
    EXPECT_PLAN("Ask[s1.0 s1.1]: 0", two);

    // Along none
    Drop none;
    none.system(1, 1);
    none.file(PlacementKind_Series);
    EXPECT_PLAN(" ! 0:NeedsRun", none);
}

// What cannot stand alone joins the system it is of, and is refused without one
UTEST(viamd_file_placement, data_needs_a_system) {
    Drop d;
    const size_t ndx = d.file(PlacementKind_Data);
    EXPECT_PLAN(" ! 0:NeedsSystem", d);

    d.system(1);
    d.system(2);
    d.of_system(ndx, 2);
    EXPECT_PLAN("Join(s2): 0", d);

    // With the structure it is of, dropped with it
    Drop with;
    const size_t gro = with.file(PlacementKind_System);
    const size_t idx = with.file(PlacementKind_Data);
    with.of_file(idx, gro);
    EXPECT_PLAN("New: 0 1", with);
}

// Dropped onto a system: into it if of its atoms, replacing what it does, and refused otherwise. Dropped below the
// systems: a new system, whatever it matches.
UTEST(viamd_file_placement, explicit_targets) {
    Drop onto;
    onto.system(1);
    onto.system(2);
    const size_t pdb = onto.file(PlacementKind_System, 1);
    const size_t xtc = onto.file(PlacementKind_Trajectory, 1);
    onto.of_system(pdb, 1, PlacementMatch_Override);
    EXPECT_PLAN("Join(s1): 0 ! 1:NotOfTarget", onto);

    // One onto the system, another of its atoms with it on its own: both into it, and a copy of what replaces
    Drop both;
    both.system(1);
    const size_t top = both.file(PlacementKind_Topology, 1);
    const size_t a = both.file(PlacementKind_System);
    const size_t b = both.file(PlacementKind_System);
    both.of_system(top, 1);
    both.of_system(a, 1);
    both.of_file(b, a, PlacementMatch_Override);
    EXPECT_PLAN("Join(s1): 0 1 | Copy(s1): 2", both);

    Drop below;
    below.system(1);
    const size_t gro = below.file(PlacementKind_System, PLACEMENT_NEW);
    below.of_system(gro, 1);
    EXPECT_PLAN("New: 0", below);
    (void)xtc;
}

// At most max_systems: a new system past them is refused, one that joins is not
UTEST(viamd_file_placement, the_cap) {
    Drop d;
    for (uint32_t i = 1; i <= 8; ++i) d.system(i);
    const size_t pdb = d.file(PlacementKind_System);
    const size_t xtc = d.file(PlacementKind_Trajectory);
    d.of_system(xtc, 4);
    (void)pdb;
    EXPECT_PLAN("Join(s4): 1@f1 ! 0:Cap", d);
}
