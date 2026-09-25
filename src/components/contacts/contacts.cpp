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

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <md_system.h>
#include <md_contact.h>
#include <md_nonbonded.h>

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <task_system.h>
#include <serialization_utils.h>

#include <imgui.h>
#include <implot.h>

#include <atomic>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace {

enum GroupBy {
    GroupBy_Residue,
    GroupBy_Molecule,
    GroupBy_Count
};

static const char* group_by_lbl[] = { "Residue", "Molecule / chain" };

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
    std::vector<md_urange_t> ranges;    // Particles of each group
    uint64_t num_group_pairs = 0;       // Pairs of groups in contact
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
    task_system::ID task = 0;
    std::atomic<bool> busy = false;
    Result result;
    Result pending;
    std::atomic<bool> pending_ready = false;
    md_system_state_t frame = {};       // Private copy of the frame the task works on

    int hovered_group = -1;

    Contacts() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app_state = (ApplicationState*)e.payload;
                break;
            case viamd::EventType_ViamdShutdown:
                if (task) task_system::task_interrupt_and_wait_for(task);
                md_system_state_free(&frame);
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Contacts", &show_window);
                break;
            case viamd::EventType_ViamdSystemInit:
                // A new system: whatever was computed or loaded belongs to the old one
                if (task) task_system::task_interrupt_and_wait_for(task);
                result = Result{};
                break;
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, STR_LIT("Contacts"));
                viamd::write_int(state, STR_LIT("GroupBy"), (int)group_by);
                viamd::write_flt(state, STR_LIT("Cutoff"), cutoff);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("Contacts"))) {
                    str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (str_eq(ident, STR_LIT("GroupBy"))) {
                            int v = 0;
                            viamd::extract_int(v, arg);
                            group_by = (GroupBy)CLAMP(v, 0, GroupBy_Count - 1);
                        } else if (str_eq(ident, STR_LIT("Cutoff"))) {
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
        if (busy) return;
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
        if (group_by == GroupBy_Residue) {
            for (size_t c = 0; c < md_system_component_count(&sys); ++c) pending.ranges.push_back(md_system_component_atom_range(&sys, c));
        } else {
            for (size_t c = 0; c < md_system_instance_count(&sys); ++c) pending.ranges.push_back(md_system_instance_atom_range(&sys, c));
        }
        pending.groups.resize(pending.ranges.size());

        busy = true;
        task = task_system::create_pool_task(STR_LIT("Contacts"), [this]() {
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
            md_bitfield_t set = md_bitfield_create(md_get_heap_allocator());
            // Particles in no group get a label of their own and are left out of the set
            for (size_t k = 0; k < N; ++k) label[k] = (uint32_t)(res.ranges.size() + k);
            for (size_t g = 0; g < res.ranges.size(); ++g) {
                for (uint32_t k = res.ranges[g].beg; k < res.ranges[g].end; ++k) label[k] = (uint32_t)g;
                md_bitfield_set_range(&set, res.ranges[g].beg, res.ranges[g].end);
            }

            md_contact_pairs_desc_t desc = {};
            desc.set_a = &set;
            desc.radius = std::max(res.cutoff + 2.0 * max_radius, res.energy_cutoff);
            desc.particle_label = label.data();
            md_contact_pairs_t pairs = {};
            md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
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
                    res.valid = true;
                }
                md_contact_pairs_free(&pairs);
            }
            md_arena_allocator_destroy(arena);
            md_bitfield_free(&set);
            res.seconds = md_tick_to_seconds(md_tick_now() - t0);
            pending_ready = true;
            busy = false;
        });
        task_system::enqueue_task(task);
    }

    void draw_window() {
        if (pending_ready.exchange(false)) {
            result = std::move(pending);
        }
        if (!show_window || !app_state) return;

        ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Contacts", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            // Settings
            ImGui::PushItemWidth(200);
            int gb = (int)group_by;
            if (ImGui::Combo("Groups", &gb, group_by_lbl, GroupBy_Count)) group_by = (GroupBy)gb;
            ImGui::InputFloat("Surface distance (Å)", &cutoff, 0.1f, 1.0f, "%.2f");
            cutoff = CLAMP(cutoff, -10.0f, 100.0f);
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

            ImGui::BeginDisabled(busy);
            if (ImGui::Button(busy ? "Computing..." : "Compute for current frame")) compute();
            ImGui::EndDisabled();

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

        hovered_group = -1;
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

        if (hovered_group >= 0 && hovered_group < (int)res.ranges.size()) {
            md_bitfield_clear(&app_state->selection.highlight_mask);
            md_bitfield_set_range(&app_state->selection.highlight_mask, res.ranges[hovered_group].beg, res.ranges[hovered_group].end);
        }
    }
};

static Contacts instance;

}  // namespace
