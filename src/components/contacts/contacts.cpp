// Contacts: a graphical front end for md_contact applied to the current state of the system.
//
// The system is divided into groups (residues or molecules). For the current frame every particle pair between
// two different groups is streamed once (md_contact_pairs), as far as either the contact criterion or the force
// field can reach. Two particles are in contact when the gap between their van der Waals spheres, d - r_i - r_j, is
// below the cutoff: the cutoff is a distance between surfaces, 0 means touching and a negative cutoff asks for
// overlap. The smallest gap over the particles of two groups is the distance between their van der Waals surfaces.
// Per group it reports the groups it is in contact with, the particle pairs in contact, and its non-bonded energy
// with every other group, as the
// simulation computed it (md_nonbonded). Pairs within a group are not contacts and are skipped.
//
// The energy is there when the system carries its non-bonded force field (md_system_t::nonbonded), which the loader
// provides when the file has it (a tpr). Otherwise there is no energy.
//
// A prototype: computed on demand in a background task from a copy of the coordinates.
//
// The groups in contact are shown as a map too, of the frame computed. Over the frames, the same is contacts() in the
// script (VIAL), with the same criterion: Over Frames adds the line, whose mark drags into the Heatmaps window.

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <md_system.h>
#include <md_contact.h>
#include <md_nonbonded.h>

#include <viamd_event.h>
#include <window_menu.h>
#include <viamd.h>
#include <script_engine.h>
#include <script_editor.h>
#include <event.h>
#include <task_system.h>
#include <serialization_utils.h>

#include <imgui.h>
#include <implot.h>

#include <atomic>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace {

// What the particles are grouped by: residues, molecules (md_molecule_data_t: atoms joined by bonds, each residue and
// polymer chain whole), chains (md_chain_data_t: all that the data gives a chain id) or parts (md_part_data_t: the
// covalent pieces of a molecule, a strand with its caps, a glycan, a ligand). A particle in no group of the kind (an
// ion in no chain) is left out. The values are those a workspace stores.
enum GroupBy {
    GroupBy_Residue,
    GroupBy_Molecule,
    GroupBy_Chain,
    GroupBy_Part,
    GroupBy_Count
};

static const char* group_by_lbl[] = { "Residue", "Molecule", "Chain", "Part" };

struct GroupResult {
    uint32_t contacts = 0;      // Groups in contact with
    uint64_t pairs = 0;         // Particle pairs within the contact cutoff
    double   e_lj = 0;
    double   e_coul = 0;
};

struct Result {
    bool   valid = false;
    bool   has_energy = false;
    double cutoff = 0;          // Å, between van der Waals surfaces
    double energy_cutoff = 0;   // Å
    GroupBy group_by = GroupBy_Residue;
    std::vector<GroupResult> groups;
    std::vector<uint32_t> group_off;    // Particles of group g: group_atom[group_off[g], group_off[g + 1])
    std::vector<uint32_t> group_atom;
    uint64_t num_group_pairs = 0;       // Pairs of groups in contact
    std::vector<uint64_t> group_pairs;  // Each of them, (gi << 32) | gj with gi < gj, ascending
    uint64_t num_pairs = 0;
    double   e_lj = 0, e_coul = 0;      // Between groups, each pair once
    double   seconds = 0;
};

struct StreamCtx {
    const md_nb_forcefield_t* ff;
    const uint32_t* label;
    const float* radius;        // Per particle, Å
    float cutoff;
    Result* res;
    std::vector<uint64_t>* keys;
};

static void stream_pairs(const uint32_t* a, const uint32_t* b, const float* r, size_t count, void* user) {
    StreamCtx* s = (StreamCtx*)user;
    Result& res = *s->res;
    for (size_t k = 0; k < count; ++k) {
        const uint32_t i = a[k], j = b[k];
        const uint32_t gi = s->label[i], gj = s->label[j];
        if (r[k] - s->radius[i] - s->radius[j] < s->cutoff) {
            res.groups[gi].pairs += 1;
            res.groups[gj].pairs += 1;
            res.num_pairs += 1;
            s->keys->push_back(((uint64_t)std::min(gi, gj) << 32) | std::max(gi, gj));
        }
        if (s->ff) {
            double elj, ec;
            md_nb_forcefield_pair_energy(s->ff, i, j, (double)r[k] * r[k] * 0.01, &elj, &ec);   // Å² -> nm²
            res.groups[gi].e_lj += elj;   res.groups[gj].e_lj += elj;
            res.groups[gi].e_coul += ec;  res.groups[gj].e_coul += ec;
            res.e_lj += elj;
            res.e_coul += ec;
        }
    }
    // Keep the keys compact
    if (s->keys->size() > (1u << 24)) {
        std::sort(s->keys->begin(), s->keys->end());
        s->keys->erase(std::unique(s->keys->begin(), s->keys->end()), s->keys->end());
    }
}

struct Contacts : viamd::EventHandler {
    bool show_window = false;
    ApplicationState* app_state = nullptr;

    // Settings
    GroupBy group_by = GroupBy_Residue;
    float cutoff = 1.0f;                // Å, gap between van der Waals surfaces

    // Computation
    task_system::Job task;
    Result result;
    Result pending;
    std::atomic<bool> pending_ready = false;
    md_system_state_t frame = {};       // Private copy of the frame the task works on

    int hovered_group = -1;
    int hovered_pair[2] = {-1, -1};     // The cell of the map under the pointer: two groups
    bool highlighting = false;          // The highlight is ours, to clear when nothing is hovered

    Contacts() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app_state = (ApplicationState*)e.payload;
                workspace_register_window("Contacts", &show_window);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                group_by = GroupBy_Residue;
                cutoff = 1.0f;
                break;
            case viamd::EventType_ViamdShutdown:
                task.cancel_and_wait();
                md_system_state_free(&frame);
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowMenuRequest:
                viamd::window_menu_add(e, "Analysis/Contacts", &show_window, {.order = 30});
                break;
            case viamd::EventType_ViamdSystemInit:
                // A new system: whatever was computed or loaded belongs to the old one
                task.cancel_and_wait();
                pending_ready = false;
                result = Result{};
                break;
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, MD_STR_LIT("Contacts"));
                viamd::write_int(state, MD_STR_LIT("GroupBy"), (int)group_by);
                viamd::write_flt(state, MD_STR_LIT("Cutoff"), cutoff);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (md_str_eq(viamd::section_header(state), MD_STR_LIT("Contacts"))) {
                    md_str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (md_str_eq(ident, MD_STR_LIT("GroupBy"))) {
                            viamd::extract_enum(group_by, arg, GroupBy_Count);
                        } else if (md_str_eq(ident, MD_STR_LIT("Cutoff"))) {
                            viamd::extract_flt(cutoff, arg);
                        }
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }

    void compute() {
        if (task.running()) return;
        const md_system_t& sys = app_state->mold.sys;
        const md_system_state_t& cur = app_state->mold.state;
        const size_t N = sys.atom.count;
        if (N == 0) return;

        // The frame as it is now: the task works on a copy
        frame.alloc = md_get_heap_allocator();
        if (!md_system_state_copy(&frame, &cur)) return;

        pending = Result{};
        pending.group_by = group_by;
        pending.cutoff = cutoff;
        std::vector<uint32_t>& off = pending.group_off;
        std::vector<uint32_t>& atoms = pending.group_atom;
        auto push_range = [&](md_urange_t r) { for (uint32_t k = r.beg; k < r.end && k < N; ++k) atoms.push_back(k); };
        off.push_back(0);
        if (group_by == GroupBy_Residue) {
            for (size_t c = 0; c < md_system_residue_count(&sys); ++c) {
                push_range(md_system_residue_atom_range(&sys, c));
                off.push_back((uint32_t)atoms.size());
            }
        } else if (group_by == GroupBy_Part) {
            for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
                push_range(md_system_part_atom_range(&sys, p));
                off.push_back((uint32_t)atoms.size());
            }
        } else if (group_by == GroupBy_Chain) {
            for (size_t c = 0; c < md_system_chain_count(&sys); ++c) {
                const md_urange_t* end = md_system_chain_run_end(&sys, c);
                for (const md_urange_t* it = md_system_chain_run_beg(&sys, c); it && it < end; ++it) {
                    push_range(md_system_chain_run_atom_range(&sys, *it));
                }
                off.push_back((uint32_t)atoms.size());
            }
        } else {
            for (size_t m = 0; m < md_system_molecule_count(&sys); ++m) {
                const md_molecule_t mol = md_system_molecule(&sys, m);
                for (size_t k = 0; k < mol.count; ++k) atoms.push_back((uint32_t)mol.atom_idx[k]);
                off.push_back((uint32_t)atoms.size());
            }
        }
        pending.groups.resize(off.size() - 1);

        task.start(MD_STR_LIT("Contacts"), [this]() {
            const md_system_t& sys = app_state->mold.sys;
            const size_t N = sys.atom.count;
            Result& res = pending;
            const md_tick_t t0 = md_tick_now();
            const md_nb_forcefield_t* ff = sys.nonbonded;
            res.has_energy = ff != nullptr;
            res.energy_cutoff = ff ? md_nb_potential_cutoff(&ff->potential) * 10.0 : 0.0;

            // Van der Waals radii: of the elements, and for coarse grained beads what their force field implies
            std::vector<float> radius(N);
            if (N) md_atom_extract_radii(radius.data(), 0, N, &sys.atom);
            float max_radius = 0.0f;
            for (float r : radius) max_radius = std::max(max_radius, r);

            std::vector<uint32_t> label(N);
            md_bitmask_t set_m = md_bitmask_create(0, (uint32_t)N, md_get_heap_allocator());
            // Particles in no group get a label of their own and are left out of the set
            const size_t num_groups = res.groups.size();
            for (size_t k = 0; k < N; ++k) label[k] = (uint32_t)(num_groups + k);
            for (size_t g = 0; g < num_groups; ++g) {
                for (uint32_t k = res.group_off[g]; k < res.group_off[g + 1]; ++k) {
                    label[res.group_atom[k]] = (uint32_t)g;
                    md_bitmask_set_bit(&set_m, res.group_atom[k]);
                }
            }
            const md_bitset_t set = md_bitmask_view(&set_m);

            md_contact_pairs_desc_t desc = {};
            desc.set_a = &set;
            desc.radius = std::max(res.cutoff + 2.0 * max_radius, res.energy_cutoff);
            desc.particle_label = label.data();
            md_contact_pairs_t pairs = {};
            md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
            std::vector<uint64_t> keys;
            if (md_contact_pairs_init(&pairs, &desc, &sys, arena)) {
                StreamCtx ctx = { ff, label.data(), radius.data(), (float)res.cutoff, &res, &keys };
                if (md_contact_pairs_for_each(&pairs, &frame, stream_pairs, &ctx)) {
                    std::sort(keys.begin(), keys.end());
                    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
                    for (uint64_t key : keys) {
                        res.groups[key >> 32].contacts += 1;
                        res.groups[key & 0xFFFFFFFFu].contacts += 1;
                    }
                    res.num_group_pairs = keys.size();
                    res.group_pairs = std::move(keys);
                    res.valid = true;
                }
                md_contact_pairs_free(&pairs);
            }
            md_arena_allocator_destroy(arena);
            md_bitmask_free(&set_m);
            res.seconds = md_tick_to_seconds(md_tick_now() - t0);
            // A cancelled one (a new system) is for the old system: not taken
            if (!task_system::cancelled()) pending_ready = true;
        });
    }

    void draw_window() {
        if (pending_ready.exchange(false)) {
            result = std::move(pending);
        }
        hovered_group = -1;
        hovered_pair[0] = hovered_pair[1] = -1;
        if (show_window && app_state) draw_contents();
        if (app_state) highlight();
    }

    // What is hovered lit up in the view: a group of the table, or the two of a cell of the map
    void highlight() {
        const Result& res = result;
        const int n = (int)res.groups.size();
        md_bitmask_t* mask = &app_state->selection.highlight_mask;
        auto light = [&](int g) {
            for (uint32_t k = res.group_off[g]; k < res.group_off[g + 1]; ++k) md_bitmask_set_bit(mask, res.group_atom[k]);
        };
        const bool group = hovered_group >= 0 && hovered_group < n;
        const bool pair = hovered_pair[0] >= 0 && hovered_pair[0] < n && hovered_pair[1] >= 0 && hovered_pair[1] < n;
        if (group || pair) {
            md_bitmask_clear(mask);
            if (group) light(hovered_group);
            if (pair) { light(hovered_pair[0]); light(hovered_pair[1]); }
            highlighting = true;
        } else if (highlighting) {
            md_bitmask_clear(mask);
            highlighting = false;
        }
    }

    // The groups there would be now, of the system as it is
    size_t num_groups_now() const {
        const md_system_t& sys = app_state->mold.sys;
        switch (group_by) {
        case GroupBy_Residue:  return md_system_residue_count(&sys);
        case GroupBy_Molecule: return md_system_molecule_count(&sys);
        case GroupBy_Chain:    return md_system_chain_count(&sys);
        case GroupBy_Part:     return md_system_part_count(&sys);
        default:               return 0;
        }
    }

    // A name the script binds nothing to: base1, base2, ...
    std::string unique_identifier(const char* base) const {
        const script::Program* ir = app_state->script.ir;
        const size_t n = script::num_identifiers(ir);
        const md_str_t* idents = script::identifiers(ir);
        char buf[64];
        for (int i = 1; i < 1000; ++i) {
            snprintf(buf, sizeof(buf), "%s%d", base, i);
            bool taken = false;
            for (size_t k = 0; k < n && !taken; ++k) taken = md_str_eq(idents[k], md_str_from_cstr(buf));
            if (!taken) return buf;
        }
        return base;
    }

    // The map over the frames: contacts() in the script, by the same criterion and groups, added as a line. Of the
    // residues every pair counts, as here: min_separation=0 (contacts() leaves out a residue's two neighbours on
    // either side by default).
    void draw_over_frames() {
        static const char* selector[GroupBy_Count] = { "residue", "molecule", "chain", "part" };
        const bool v2 = app_state->script.language == script::Language::V2;
        const size_t n = num_groups_now();
        const bool too_many = n > 5000;
        ImGui::BeginDisabled(!v2 || too_many || n < 2);
        if (ImGui::Button("Over frames")) {
            char args[128];
            int len = snprintf(args, sizeof(args), "%s", selector[group_by]);
            if (cutoff != 1.0f && len < (int)sizeof(args)) len += snprintf(args + len, sizeof(args) - len, ", surface_distance=%.3g", cutoff);
            if (group_by == GroupBy_Residue && len < (int)sizeof(args)) snprintf(args + len, sizeof(args) - len, ", min_separation=0");
            char line[256];
            snprintf(line, sizeof(line), "%s = contacts(%s)", unique_identifier("cmap").c_str(), args);
            script_editor::append_line(app_state->editor, md_str_from_cstr(line));
        }
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
            if (!v2) {
                ImGui::SetTooltip("Over the frames, the map is contacts() of a VIAL script: the script is in v1");
            } else if (too_many) {
                ImGui::SetTooltip("%zu groups: a map of %zu × %zu over the frames is too large.\n"
                                  "Write contacts() of fewer in the script, such as contacts(protein.residue).", n, n, n);
            } else {
                ImGui::SetTooltip("Adds contacts() of these groups by this surface distance to the script, as a line: evaluated, it is\n"
                                  "the map over every frame. Drag its mark from the gutter into the Heatmaps window to see it over\n"
                                  "all frames or a period.");
            }
        }
    }

    void group_label(char* buf, size_t cap, int g) const {
        const md_system_t& sys = app_state->mold.sys;
        switch (result.group_by) {
        case GroupBy_Residue: {
            const md_str_t name = md_system_residue_name(&sys, (size_t)g);
            snprintf(buf, cap, MD_STR_FMT "%d", MD_STR_ARG(name), (int)md_system_residue_seq_id(&sys, (size_t)g));
            break;
        }
        case GroupBy_Chain: {
            const md_str_t id = md_system_chain_id(&sys, (size_t)g);
            if (id.len) snprintf(buf, cap, "Chain " MD_STR_FMT, MD_STR_ARG(id));
            else snprintf(buf, cap, "Chain %d", g + 1);
            break;
        }
        default:
            snprintf(buf, cap, "%s %d", group_by_lbl[result.group_by], g + 1);
            break;
        }
    }

    // The groups in contact at the frame computed, group 0 at the top left. At most 256 cells a side of what is in view
    // are drawn, a cell lit where any pair of its groups is in contact, so zooming in shows every group.
    void draw_map() {
        const Result& res = result;
        const uint32_t n = (uint32_t)res.groups.size();
        if (n == 0) return;
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        if (ImPlot::BeginPlot("##contact_map", ImVec2(-1, MD_MAX(avail.y, 200.0f)), ImPlotFlags_NoLegend | ImPlotFlags_NoTitle | ImPlotFlags_NoMouseText | ImPlotFlags_Equal)) {
            ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoGridLines);
            ImPlot::SetupAxisFormat(ImAxis_X1, [](double v, char* buf, int size, void*) { return snprintf(buf, size, "%d", (int)floor(v + 1e-6) + 1); }, nullptr);
            ImPlot::SetupAxisFormat(ImAxis_Y1, [](double v, char* buf, int size, void* d) { return snprintf(buf, size, "%d", (int)floor(*(const uint32_t*)d - v + 1e-6) + 1); }, (void*)&n);
            ImPlot::SetupAxesLimits(0, n, 0, n, ImPlotCond_Once);
            const ImPlotRect lim = ImPlot::GetPlotLimits();
            const uint32_t c0 = (uint32_t)MD_CLAMP(floor(lim.X.Min), 0.0, (double)n), c1 = (uint32_t)MD_CLAMP(ceil(lim.X.Max), 0.0, (double)n);
            const uint32_t r0 = (uint32_t)MD_CLAMP(floor((double)n - lim.Y.Max), 0.0, (double)n), r1 = (uint32_t)MD_CLAMP(ceil((double)n - lim.Y.Min), 0.0, (double)n);
            if (c1 > c0 && r1 > r0) {
                const uint32_t bc = (c1 - c0 + 255) / 256, br = (r1 - r0 + 255) / 256;
                const uint32_t nc = (c1 - c0 + bc - 1) / bc, nr = (r1 - r0 + br - 1) / br;
                std::vector<float> cells((size_t)nr * nc, 0.0f);
                auto mark = [&](uint32_t r, uint32_t c) {
                    if (r < r0 || r >= r1 || c < c0 || c >= c1) return;
                    cells[(size_t)((r - r0) / br) * nc + (c - c0) / bc] = 1.0f;
                };
                for (uint64_t key : res.group_pairs) {
                    const uint32_t a = (uint32_t)(key >> 32), b = (uint32_t)(key & 0xFFFFFFFFu);
                    mark(a, b);
                    mark(b, a);
                }
                ImPlot::PlotHeatmap("##cells", cells.data(), (int)nr, (int)nc, 0.0, 1.0, nullptr,
                                    ImPlotPoint(c0, (double)n - (r0 + (double)nr * br)), ImPlotPoint(c0 + (double)nc * bc, (double)n - r0));
            }
            if (ImPlot::IsPlotHovered()) {
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                const long c = (long)floor(mp.x), r = (long)floor((double)n - mp.y);
                if (c >= 0 && r >= 0 && (uint32_t)c < n && (uint32_t)r < n && r != c) {
                    hovered_pair[0] = (int)r;
                    hovered_pair[1] = (int)c;
                    const uint64_t key = ((uint64_t)MD_MIN(r, c) << 32) | (uint64_t)MD_MAX(r, c);
                    const bool in = std::binary_search(res.group_pairs.begin(), res.group_pairs.end(), key);
                    const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint((double)c, (double)n - r));
                    const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint((double)c + 1, (double)n - r - 1));
                    ImPlot::GetPlotDrawList()->AddRect(p0, p1, IM_COL32(255, 255, 255, 220), 0.0f, 0, 1.5f);
                    char la[64], lb[64];
                    group_label(la, sizeof(la), (int)r);
                    group_label(lb, sizeof(lb), (int)c);
                    ImGui::SetTooltip("%s and %s\n%s", la, lb, in ? "In contact" : "Not in contact");
                }
            }
            ImPlot::EndPlot();
        }
        ImPlot::PopColormap();
    }

    void draw_contents() {
        ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Contacts", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            // Settings
            ImGui::PushItemWidth(200);
            int gb = (int)group_by;
            if (ImGui::Combo("Groups", &gb, group_by_lbl, GroupBy_Count)) group_by = (GroupBy)gb;
            ImGui::SetItemTooltip("Residue: each residue.\nMolecule: atoms joined by bonds, each residue and polymer chain whole.\n"
                                  "Chain: all that the file gives one chain id (none where it gives no ids: see Part).\n"
                                  "Part: a covalent piece of a molecule: a strand of a polymer with its caps, a glycan, a ligand.\n"
                                  "Particles in no group of the kind are left out.");
            ImGui::InputFloat("Surface distance (Å)", &cutoff, 0.1f, 1.0f, "%.2f");
            cutoff = MD_CLAMP(cutoff, -10.0f, 100.0f);
            ImGui::SetItemTooltip("Contact when the van der Waals surfaces of two particles are closer than this.\n0 is touching, a negative value asks for overlap.");
            ImGui::PopItemWidth();

            if (const md_nb_forcefield_t* ff = app_state->mold.sys.nonbonded) {
                static const char* mod_lbl[]  = { "plain cut-off", "potential-shift", "potential-switch", "force-switch" };
                static const char* coul_lbl[] = { "none", "cut-off", "reaction field", "Ewald real space" };
                ImGui::TextDisabled("Energy: from the force field of the system: %zu particle types, LJ %s, Coulomb %s, pairs within %.2f nm",
                    ff->num_types, mod_lbl[ff->potential.lj_modifier], coul_lbl[ff->potential.coulomb], md_nb_potential_cutoff(&ff->potential));
            } else {
                ImGui::TextDisabled("Energy: not available, the system carries no non-bonded force field");
            }

            const bool busy = task.running();
            ImGui::BeginDisabled(busy);
            if (ImGui::Button(busy ? "Computing..." : "Compute for current frame")) compute();
            ImGui::EndDisabled();
            ImGui::SameLine();
            draw_over_frames();

            ImGui::Separator();
            draw_result();
        }
        ImGui::End();
    }

    void draw_result() {
        const Result& res = result;
        if (!res.valid) {
            ImGui::TextDisabled("No result");
            return;
        }
        ImGui::Text("%zu groups (%s), surfaces within %.2f Å: %llu group pairs in contact, %llu particle pairs (%.2f s)",
            res.groups.size(), group_by_lbl[res.group_by], res.cutoff, (unsigned long long)res.num_group_pairs, (unsigned long long)res.num_pairs, res.seconds);
        if (res.has_energy) {
            ImGui::Text("Energy between groups: LJ %.3f, Coulomb %.3f, total %.3f kJ/mol (pairs within %.2f nm)",
                res.e_lj, res.e_coul, res.e_lj + res.e_coul, res.energy_cutoff * 0.1);
        }

        if (!ImGui::BeginTabBar("##contact_views")) return;
        if (ImGui::BeginTabItem("Map")) {
            draw_map();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Groups")) {
            draw_groups();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }

    void draw_groups() {
        const Result& res = result;
        // Distributions over the groups: contacts, and energy when there is one. Side by side when both.
        std::vector<double> v(res.groups.size());
        const float plot_w = res.has_energy ? (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f : -1.0f;
        if (ImPlot::BeginPlot("##contacts_dist", ImVec2(plot_w, 180))) {
            for (size_t g = 0; g < v.size(); ++g) v[g] = (double)res.groups[g].contacts;
            ImPlot::SetupAxes("Contacts per group", "Groups", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotHistogram("##contacts", v.data(), (int)v.size(), ImPlotBin_Sturges);
            ImPlot::EndPlot();
        }
        if (res.has_energy) {
            ImGui::SameLine();
            if (ImPlot::BeginPlot("##energy_dist", ImVec2(plot_w, 180))) {
                for (size_t g = 0; g < v.size(); ++g) v[g] = res.groups[g].e_lj + res.groups[g].e_coul;
                ImPlot::SetupAxes("Energy per group (kJ/mol)", "Groups", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotHistogram("##energy", v.data(), (int)v.size(), ImPlotBin_Sturges);
                ImPlot::EndPlot();
            }
        }

        const int num_cols = res.has_energy ? 6 : 3;
        const ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_Sortable | ImGuiTableFlags_Resizable;
        if (ImGui::BeginTable("groups", num_cols, flags)) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Group", ImGuiTableColumnFlags_DefaultSort);
            ImGui::TableSetupColumn("Contacts");
            ImGui::TableSetupColumn("Pairs");
            if (res.has_energy) {
                ImGui::TableSetupColumn("LJ");
                ImGui::TableSetupColumn("Coulomb");
                ImGui::TableSetupColumn("Total");
            }
            ImGui::TableHeadersRow();

            static std::vector<int> order;
            order.resize(res.groups.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
            if (ImGuiTableSortSpecs* spec = ImGui::TableGetSortSpecs(); spec && spec->SpecsCount > 0) {
                const int col = spec->Specs[0].ColumnIndex;
                const bool asc = spec->Specs[0].SortDirection == ImGuiSortDirection_Ascending;
                auto key = [&](int g) -> double {
                    const GroupResult& r = res.groups[g];
                    switch (col) {
                    case 1: return r.contacts;
                    case 2: return (double)r.pairs;
                    case 3: return r.e_lj;
                    case 4: return r.e_coul;
                    case 5: return r.e_lj + r.e_coul;
                    default: return g;
                    }
                };
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return asc ? key(a) < key(b) : key(a) > key(b); });
            }

            ImGuiListClipper clipper;
            clipper.Begin((int)order.size());
            while (clipper.Step()) {
                for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                    const int g = order[row];
                    const GroupResult& r = res.groups[g];
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char lbl[32];
                    snprintf(lbl, sizeof(lbl), "%d", g + 1);
                    ImGui::Selectable(lbl, false, ImGuiSelectableFlags_SpanAllColumns);
                    if (ImGui::IsItemHovered()) hovered_group = g;
                    ImGui::TableNextColumn(); ImGui::Text("%u", r.contacts);
                    ImGui::TableNextColumn(); ImGui::Text("%llu", (unsigned long long)r.pairs);
                    if (res.has_energy) {
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", r.e_lj);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", r.e_coul);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", r.e_lj + r.e_coul);
                    }
                }
            }
            ImGui::EndTable();
        }
    }
};

static Contacts instance;

}  // namespace
