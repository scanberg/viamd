// Ring puckering: the shape of the rings of five to seven atoms of the system at the frame shown (md_pucker_ring).
//
// Each ring is put in the order the conventions are written for (md_pucker_ring_order) and given its Cremer-Pople
// coordinates, the Altona-Sundaralingam pseudorotation of a furanose, and the name of the nearest canonical
// conformation, the same as pucker() of the script and its hovers give. Two maps show where the rings are:
//   - Five membered rings on the pseudorotation wheel: the phase P = phi2 + 90 around it (0 at the top, clockwise),
//     the amplitude Q out from the centre. P of a furanose is its Altona-Sundaralingam phase to within a degree.
//   - Six membered rings on (phi, theta): the chairs at the top and the bottom, the boats and skews along the middle.
// The shading is the region of each kind of conformation, the labels the canonical conformations: by IUPAC with the
// locants of a sugar when every ring on the map is one, by their kind otherwise. Rings flatter than
// MD_PUCKER_PLANAR_AMPLITUDE have no shape to place and stay off the maps.
//
// The frame shown only. For how a ring moves over the trajectory use pucker() and its accessors in the script.

#include <core/md_allocator.h>
#include <core/md_bitset.h>
#include <core/md_str.h>

#include <md_system.h>
#include <md_pucker.h>

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <serialization_utils.h>

#include <imgui.h>
#include <implot.h>

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cfloat>

namespace {

enum Kind {
    Kind_All,
    Kind_NotAromatic,
    Kind_Sugar,
    Kind_Count
};
static const char* kind_lbl[] = { "All rings", "Not aromatic", "Sugars" };

enum Size {
    Size_Any,
    Size_5,
    Size_6,
    Size_7,
    Size_Count
};
static const char* size_lbl[] = { "Any size", "5 atoms", "6 atoms", "7 atoms" };

struct Row {
    md_pucker_ring_t r = {};
    uint32_t ring_idx = 0;
    int32_t  res_idx = -1;
    md_ring_flags_t flags = MD_RING_FLAG_NONE;
    bool     valid = false;     // Coordinates were computed (not a degenerate ring)
};

struct NameCount {
    std::string name;
    int count;
};

static constexpr double DEG = 3.14159265358979323846 / 180.0;

static inline float wheel_phase(const md_pucker_ring_t& r) {
    float p = r.cp.phi2 + 90.0f;
    return p >= 360.0f ? p - 360.0f : p;
}

static ImU32 kind_color(md_pucker_kind_t kind, float alpha) {
    ImVec4 c;
    switch (kind) {
    case MD_PUCKER_KIND_CHAIR:      c = ImVec4(0.30f, 0.55f, 0.90f, 1); break;
    case MD_PUCKER_KIND_HALF_CHAIR: c = ImVec4(0.35f, 0.75f, 0.75f, 1); break;
    case MD_PUCKER_KIND_ENVELOPE:   c = ImVec4(0.45f, 0.75f, 0.35f, 1); break;
    case MD_PUCKER_KIND_BOAT:       c = ImVec4(0.90f, 0.55f, 0.25f, 1); break;
    case MD_PUCKER_KIND_SKEW:       c = ImVec4(0.85f, 0.35f, 0.45f, 1); break;
    case MD_PUCKER_KIND_TWIST:      c = ImVec4(0.65f, 0.45f, 0.85f, 1); break;
    default:                        c = ImVec4(0.5f, 0.5f, 0.5f, 1); break;
    }
    c.w = alpha;
    return ImGui::ColorConvertFloat4ToU32(c);
}

static char kind_letter(md_pucker_kind_t kind) {
    switch (kind) {
    case MD_PUCKER_KIND_CHAIR:      return 'C';
    case MD_PUCKER_KIND_BOAT:       return 'B';
    case MD_PUCKER_KIND_SKEW:       return 'S';
    case MD_PUCKER_KIND_HALF_CHAIR: return 'H';
    case MD_PUCKER_KIND_ENVELOPE:   return 'E';
    case MD_PUCKER_KIND_TWIST:      return 'T';
    default:                        return '?';
    }
}

// A conformation as IUPAC writes it, centred on pos: the atoms above raised before the letter, those below lowered
// after it. Without locants the letter alone.
static void draw_conformer(ImDrawList* dl, ImVec2 pos, md_pucker_conformer_t conf, const char* const* locants, ImU32 col, float scale = 1.0f) {
    ImFont* font = ImGui::GetFont();
    const float fs = ImGui::GetFontSize() * scale;
    const float ss = fs * 0.72f;
    char letter[2] = { kind_letter(conf.kind), '\0' };
    char above[16] = "", below[16] = "";
    if (locants) {
        md_pucker_locants(above, sizeof(above), conf.above, locants);
        md_pucker_locants(below, sizeof(below), conf.below, locants);
    }
    const ImVec2 wa = font->CalcTextSizeA(ss, FLT_MAX, 0.0f, above);
    const ImVec2 wl = font->CalcTextSizeA(fs, FLT_MAX, 0.0f, letter);
    const ImVec2 wb = font->CalcTextSizeA(ss, FLT_MAX, 0.0f, below);
    float x = pos.x - (wa.x + wl.x + wb.x) * 0.5f;
    const float y = pos.y - wl.y * 0.5f;
    if (above[0]) dl->AddText(font, ss, ImVec2(x, y - fs * 0.30f), col, above);
    x += wa.x;
    dl->AddText(font, fs, ImVec2(x, y), col, letter);
    x += wl.x;
    if (below[0]) dl->AddText(font, ss, ImVec2(x, y + fs * 0.45f), col, below);
}

struct Pucker : viamd::EventHandler {
    bool show_window = false;
    ApplicationState* app_state = nullptr;

    // Settings
    int  kind = Kind_NotAromatic;
    int  size = Size_Any;
    bool only_selected = false;

    // The rings at the frame shown
    std::vector<Row> rows;
    bool   dirty = true;
    size_t num_rings = 0;           // Of the system when rows were computed: the topology can change

    // The rows shown, and what they were filtered with
    std::vector<int> shown;
    std::vector<NameCount> names;
    uint64_t shown_key = 0;
    bool     shown_dirty = true;
    int      num_planar = 0;

    int  hovered = -1;              // Row hovered in the table or on a map, this frame
    bool highlighting = false;      // The highlight is ours, to clear when nothing is hovered

    md_pucker_canonical_t canon5[20] = {};
    md_pucker_canonical_t canon6[38] = {};
    std::vector<md_pucker_kind_t> regions6;     // Kind of conformation per cell of the (phi, theta) map

    static constexpr int REGION_PHI = 120;      // 3 degree cells
    static constexpr int REGION_THETA = 60;

    Pucker() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app_state = (ApplicationState*)e.payload;
                workspace_register_window("Ring Pucker", &show_window);
                init_tables();
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                kind = Kind_NotAromatic;
                size = Size_Any;
                only_selected = false;
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Ring Pucker", &show_window);
                break;
            case viamd::EventType_ViamdSystemInit:
            case viamd::EventType_ViamdSystemFree:
                rows.clear();
                shown.clear();
                names.clear();
                hovered = -1;
                highlighting = false;
                dirty = true;
                break;
            case viamd::EventType_ViamdSystemStateChanged:
                dirty = true;
                break;
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, STR_LIT("RingPucker"));
                viamd::write_int(state, STR_LIT("Kind"), kind);
                viamd::write_int(state, STR_LIT("Size"), size);
                viamd::write_bool(state, STR_LIT("OnlySelected"), only_selected);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("RingPucker"))) {
                    str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (str_eq(ident, STR_LIT("Kind"))) {
                            viamd::extract_enum(kind, arg, Kind_Count);
                        } else if (str_eq(ident, STR_LIT("Size"))) {
                            viamd::extract_enum(size, arg, Size_Count);
                        } else if (str_eq(ident, STR_LIT("OnlySelected"))) {
                            viamd::extract_bool(only_selected, arg);
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

    void init_tables() {
        md_pucker_canonical(canon5, ARRAY_SIZE(canon5), 5);
        md_pucker_canonical(canon6, ARRAY_SIZE(canon6), 6);
        // The region of each kind on the (phi, theta) map: the kind nearest the centre of each cell, as the names are
        // given (md_pucker_conformer)
        regions6.resize(REGION_PHI * REGION_THETA);
        for (int t = 0; t < REGION_THETA; ++t) {
            for (int f = 0; f < REGION_PHI; ++f) {
                md_pucker_t p = {};
                p.size = 6;
                p.amplitude = 1.0f;
                p.theta = (t + 0.5f) * 180.0f / REGION_THETA;
                p.phi2  = (f + 0.5f) * 360.0f / REGION_PHI;
                regions6[t * REGION_PHI + f] = md_pucker_conformer(&p).kind;
            }
        }
    }

    void compute() {
        const md_system_t& sys = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;
        num_rings = md_index_data_num_ranges(&sys.ring);
        rows.clear();
        if (!state.xyz) return;
        for (size_t i = 0; i < num_rings; ++i) {
            const size_t n = md_index_range_size(&sys.ring, i);
            if (n < 5 || n > 7) continue;
            Row row;
            row.ring_idx = (uint32_t)i;
            row.flags = md_system_ring_flags(&sys, i);
            row.valid = md_pucker_ring(&row.r, &sys, &state, i);
            if (row.r.size == 0) continue;
            row.res_idx = md_system_residue_find_by_atom_idx(&sys, row.r.atoms[0]);
            rows.push_back(row);
        }
        dirty = false;
        shown_dirty = true;
    }

    bool passes(const Row& row) const {
        if (kind == Kind_NotAromatic && (row.flags & MD_RING_FLAG_AROMATIC)) return false;
        if (kind == Kind_Sugar && !(row.flags & MD_RING_FLAG_SUGAR)) return false;
        if (size != Size_Any && row.r.size != 4 + size) return false;
        if (only_selected) {
            for (int k = 0; k < row.r.size; ++k) {
                if (!md_bitmask_test_bit(&app_state->selection.selection_mask, (uint32_t)row.r.atoms[k])) return false;
            }
        }
        return true;
    }

    void update_shown() {
        uint64_t key = (uint64_t)kind | ((uint64_t)size << 4) | ((uint64_t)only_selected << 8);
        if (only_selected) key ^= md_bitset_hash64(md_bitmask_view(&app_state->selection.selection_mask), 0) << 9;
        if (!shown_dirty && key == shown_key) return;
        shown_key = key;
        shown_dirty = false;
        shown.clear();
        num_planar = 0;
        std::unordered_map<std::string, int> count;
        for (size_t i = 0; i < rows.size(); ++i) {
            if (!passes(rows[i])) continue;
            shown.push_back((int)i);
            const md_pucker_ring_t& r = rows[i].r;
            if (rows[i].valid && r.cp.amplitude < MD_PUCKER_PLANAR_AMPLITUDE) num_planar += 1;
            count[r.name[0] ? std::string(r.name) : std::string("not named")] += 1;
        }
        names.clear();
        for (auto& [name, c] : count) names.push_back({name, c});
        std::sort(names.begin(), names.end(), [](const NameCount& a, const NameCount& b) { return a.count != b.count ? a.count > b.count : a.name < b.name; });
    }

    bool is_highlighted(const Row& row) const {
        const md_bitmask_t* mask = &app_state->selection.highlight_mask;
        for (int k = 0; k < row.r.size; ++k) {
            if (md_bitmask_test_bit(mask, (uint32_t)row.r.atoms[k])) return true;
        }
        return false;
    }

    void residue_label(char* buf, size_t cap, const Row& row) const {
        const md_system_t& sys = app_state->mold.sys;
        buf[0] = '\0';
        if (row.res_idx < 0) return;
        const str_t name = md_residue_name(&sys.residue, row.res_idx);
        const int seq = md_residue_seq_id(&sys.residue, row.res_idx);
        // The chain id the data gives, none for an inferred chain or a molecule in none
        const md_chain_idx_t c = md_system_chain_find_by_atom_idx(&sys, row.r.atoms[0]);
        const str_t chain = c >= 0 ? md_system_chain_id(&sys, (size_t)c) : str_t{};
        if (chain.len) snprintf(buf, cap, "%.*s %.*s %d", (int)chain.len, chain.ptr, (int)name.len, name.ptr, seq);
        else           snprintf(buf, cap, "%.*s %d", (int)name.len, name.ptr, seq);
    }

    void atoms_label(char* buf, size_t cap, const Row& row) const {
        const md_system_t& sys = app_state->mold.sys;
        size_t len = 0;
        buf[0] = '\0';
        for (int k = 0; k < row.r.size && len < cap; ++k) {
            const str_t s = md_atom_name(&sys.atom, row.r.atoms[k]);
            len += (size_t)snprintf(buf + len, cap - len, "%s%.*s", k ? " " : "", (int)s.len, s.ptr);
        }
    }

    void row_tooltip(const Row& row) const {
        char res[64], atoms[96];
        residue_label(res, sizeof(res), row);
        atoms_label(atoms, sizeof(atoms), row);
        const md_pucker_ring_t& r = row.r;
        ImGui::BeginTooltip();
        ImGui::Text("Ring %u, %s", row.ring_idx + 1, res);
        ImGui::TextDisabled("%s", atoms);
        if (r.name[0]) ImGui::Text("%s", r.name);
        ImGui::Text("Q %.3f Å", r.cp.amplitude);
        if (r.size == 6)  ImGui::Text("θ %.1f°  φ %.1f°", r.cp.theta, r.cp.phi2);
        else              ImGui::Text("φ2 %.1f°", r.cp.phi2);
        if (r.has_pseudo) ImGui::Text("P %.1f°  τm %.1f°", r.pseudo_phase, r.pseudo_amplitude);
        if (!r.anchored && !r.sugar) ImGui::TextDisabled("The order of this ring is not fixed by its atoms:\nonly Q and the kind of shape mean something");
        ImGui::EndTooltip();
    }

    void draw_window() {
        if (!show_window || !app_state) {
            if (highlighting && app_state) {
                md_bitmask_clear(&app_state->selection.highlight_mask);
                highlighting = false;
            }
            return;
        }
        const md_system_t& sys = app_state->mold.sys;
        if (md_index_data_num_ranges(&sys.ring) != num_rings) dirty = true;
        if (dirty) compute();
        update_shown();

        const int prev_hovered = hovered;
        hovered = -1;

        ImGui::SetNextWindowSize(ImVec2(760, 720), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Ring Pucker", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::PushItemWidth(150);
            ImGui::Combo("##kind", &kind, kind_lbl, Kind_Count);
            ImGui::SameLine();
            ImGui::Combo("##size", &size, size_lbl, Size_Count);
            ImGui::PopItemWidth();
            ImGui::SameLine();
            ImGui::Checkbox("Only selected", &only_selected);
            ImGui::SetItemTooltip("Only rings with every atom selected");
            update_shown();

            ImGui::TextDisabled("%zu of %zu rings (5 to 7 atoms), at the frame shown", shown.size(), rows.size());
            ImGui::SetItemTooltip("Over the trajectory: pucker() in the script, with amplitude(), phase(), theta() and pseudorotation()");

            if (rows.empty()) {
                ImGui::TextDisabled("The system has no rings of 5 to 7 atoms");
            } else {
                draw_maps(prev_hovered);
                if (ImGui::BeginTabBar("##views")) {
                    if (ImGui::BeginTabItem("Rings")) {
                        draw_table(prev_hovered);
                        ImGui::EndTabItem();
                    }
                    char lbl[48];
                    snprintf(lbl, sizeof(lbl), "Conformations (%zu)###names", names.size());
                    if (ImGui::BeginTabItem(lbl)) {
                        draw_names();
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
        }
        ImGui::End();

        if (hovered >= 0 && hovered < (int)rows.size()) {
            const Row& row = rows[hovered];
            md_bitmask_clear(&app_state->selection.highlight_mask);
            for (int k = 0; k < row.r.size; ++k) md_bitmask_set_bit(&app_state->selection.highlight_mask, (uint32_t)row.r.atoms[k]);
            highlighting = true;
        } else if (highlighting) {
            md_bitmask_clear(&app_state->selection.highlight_mask);
            highlighting = false;
        }
    }

    // Whether every ring of the given size on a map is a sugar: the labels take the locants of a sugar
    bool all_sugars(int ring_size) const {
        bool any = false;
        for (int i : shown) {
            const Row& row = rows[i];
            if (row.r.size != ring_size) continue;
            if (!row.r.sugar) return false;
            any = true;
        }
        return any;
    }

    void draw_maps(int prev_hovered) {
        int n5 = 0, n6 = 0;
        for (int i : shown) {
            const Row& row = rows[i];
            if (!row.valid || row.r.cp.amplitude < MD_PUCKER_PLANAR_AMPLITUDE) continue;
            n5 += row.r.size == 5;
            n6 += row.r.size == 6;
        }
        const bool show5 = size == Size_Any || size == Size_5;
        const bool show6 = size == Size_Any || size == Size_6;
        if (!show5 && !show6) return;

        const float avail = ImGui::GetContentRegionAvail().x;
        const float spacing = ImGui::GetStyle().ItemSpacing.x;
        const float h = 300.0f;
        const float w5 = show6 ? h : avail;
        const float w6 = show5 ? std::max(avail - w5 - spacing, 200.0f) : avail;

        if (show5) {
            draw_wheel(ImVec2(w5, h), prev_hovered, n5);
            if (show6) ImGui::SameLine();
        }
        if (show6) draw_sphere(ImVec2(w6, h), prev_hovered, n6);
        if (num_planar > 0) {
            ImGui::TextDisabled("%d planar rings (Q under %.1f Å) are not on the maps", num_planar, MD_PUCKER_PLANAR_AMPLITUDE);
        }
    }

    // Nearest point to the mouse within a few pixels: the row hovered on a map
    int pick(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<int>& idx) const {
        if (!ImPlot::IsPlotHovered()) return -1;
        const ImVec2 m = ImGui::GetMousePos();
        float best = 8.0f * 8.0f;
        int hit = -1;
        for (size_t i = 0; i < xs.size(); ++i) {
            const ImVec2 p = ImPlot::PlotToPixels(xs[i], ys[i]);
            const float d = (p.x - m.x) * (p.x - m.x) + (p.y - m.y) * (p.y - m.y);
            if (d < best) { best = d; hit = idx[i]; }
        }
        return hit;
    }

    void plot_points(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<int>& idx, int prev_hovered) {
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3.0f, ImVec4(0.95f, 0.95f, 0.95f, 0.85f), 1.0f, ImVec4(0.1f, 0.1f, 0.1f, 0.9f));
        ImPlot::PlotScatter("##rings", xs.data(), ys.data(), (int)xs.size());
        std::vector<double> hx, hy;
        for (size_t i = 0; i < xs.size(); ++i) {
            if (idx[i] == prev_hovered || is_highlighted(rows[idx[i]])) {
                hx.push_back(xs[i]);
                hy.push_back(ys[i]);
            }
        }
        if (!hx.empty()) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 6.0f, ImVec4(1.0f, 0.75f, 0.2f, 1.0f), 1.5f, ImVec4(0, 0, 0, 1));
            ImPlot::PlotScatter("##highlighted", hx.data(), hy.data(), (int)hx.size());
        }
    }

    // Five membered rings: P around, Q out
    void draw_wheel(ImVec2 sz, int prev_hovered, int count) {
        std::vector<double> xs, ys;
        std::vector<int> idx;
        double qmax = 0.6;
        for (int i : shown) {
            const Row& row = rows[i];
            if (row.r.size != 5 || !row.valid || row.r.cp.amplitude < MD_PUCKER_PLANAR_AMPLITUDE) continue;
            const double q = row.r.cp.amplitude, p = wheel_phase(row.r) * DEG;
            xs.push_back(q * sin(p));
            ys.push_back(q * cos(p));
            idx.push_back(i);
            qmax = std::max(qmax, q);
        }
        qmax = std::ceil(qmax / 0.2) * 0.2;
        const double lim = qmax * 1.32;
        const bool sugars = all_sugars(5);
        static const char* loc5[] = { "O", "1", "2", "3", "4" };

        const ImPlotAxisFlags ax = ImPlotAxisFlags_NoDecorations | ImPlotAxisFlags_Lock;
        // Locked limits with the aspect of the plot, so that a unit is as many pixels across as up: the wheel is round
        const ImVec2 pad = ImPlot::GetStyle().PlotPadding;
        const double aspect = std::max(1.0f, sz.x - 2 * pad.x) / std::max(1.0f, sz.y - 2 * pad.y);
        const double xlim = lim * std::max(1.0, aspect), ylim = lim * std::max(1.0, 1.0 / aspect);
        if (ImPlot::BeginPlot("##wheel", sz, ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
            ImPlot::SetupAxes(nullptr, nullptr, ax, ax);
            ImPlot::SetupAxesLimits(-xlim, xlim, -ylim, ylim, ImPlotCond_Always);
            ImPlot::SetupFinish();

            ImDrawList* dl = ImPlot::GetPlotDrawList();
            const ImVec2 c = ImPlot::PlotToPixels(0.0, 0.0);
            const float r_px = c.y - ImPlot::PlotToPixels(0.0, qmax).y;
            ImPlot::PushPlotClipRect();
            // The sectors of the envelopes and twists, each centred on its canonical phase
            for (const md_pucker_canonical_t& k : canon5) {
                const float p = (float)((k.phi2 + 90.0) * DEG);
                const float a0 = p - (float)(9.0 * DEG) - (float)(90.0 * DEG);
                const float a1 = p + (float)(9.0 * DEG) - (float)(90.0 * DEG);
                dl->PathLineTo(c);
                dl->PathArcTo(c, r_px, a0, a1, 8);
                dl->PathFillConvex(kind_color(k.conf.kind, 0.18f));
            }
            const ImU32 grid = ImGui::GetColorU32(ImGuiCol_TextDisabled, 0.5f);
            for (double q = 0.2; q <= qmax + 1e-6; q += 0.2) {
                const float rr = c.y - ImPlot::PlotToPixels(0.0, q).y;
                dl->AddCircle(c, rr, grid, 64);
                char t[16];
                snprintf(t, sizeof(t), "%.1f", q);
                dl->AddText(ImVec2(c.x + 2, c.y - rr), grid, t);
            }
            for (int k = 0; k < 10; ++k) {
                const float a = (float)(k * 36.0 * DEG);
                dl->AddLine(c, ImVec2(c.x + r_px * sinf(a), c.y - r_px * cosf(a)), grid);
            }
            ImPlot::PopPlotClipRect();
            // The canonical conformations around the rim
            const ImU32 text = ImGui::GetColorU32(ImGuiCol_Text);
            for (const md_pucker_canonical_t& k : canon5) {
                const double p = (k.phi2 + 90.0) * DEG;
                const double rl = qmax * 1.16;
                draw_conformer(dl, ImPlot::PlotToPixels(rl * sin(p), rl * cos(p)), k.conf, sugars ? loc5 : nullptr, text);
            }
            plot_points(xs, ys, idx, prev_hovered);
            const int hit = pick(xs, ys, idx);
            if (hit >= 0) {
                hovered = hit;
                row_tooltip(rows[hit]);
            }
            ImPlot::EndPlot();
        }
        if (ImGui::IsItemHovered() && hovered < 0) {
            ImGui::SetTooltip("Five membered rings: %d\nP = φ2 + 90° around, from the top, clockwise\nQ (Å) out from the centre%s",
                count, sugars ? "\nLabels: furanose numbering, O for the ring O" : "\nE envelope, T twist");
        }
    }

    // Six membered rings: phi across, theta down
    void draw_sphere(ImVec2 sz, int prev_hovered, int count) {
        std::vector<double> xs, ys;
        std::vector<int> idx;
        for (int i : shown) {
            const Row& row = rows[i];
            if (row.r.size != 6 || !row.valid || row.r.cp.amplitude < MD_PUCKER_PLANAR_AMPLITUDE) continue;
            xs.push_back(row.r.cp.phi2);
            ys.push_back(row.r.cp.theta);
            idx.push_back(i);
        }
        const bool sugars = all_sugars(6);
        static const char* loc6[] = { "O", "1", "2", "3", "4", "5" };

        if (ImPlot::BeginPlot("##sphere", sz, ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
            ImPlot::SetupAxes("φ (°)", "θ (°)", ImPlotAxisFlags_Lock, ImPlotAxisFlags_Lock | ImPlotAxisFlags_Invert);
            ImPlot::SetupAxesLimits(0, 360, -8, 188, ImPlotCond_Always);
            static const double xt[] = { 0, 60, 120, 180, 240, 300, 360 };
            static const double yt[] = { 0, 45, 90, 135, 180 };
            ImPlot::SetupAxisTicks(ImAxis_X1, xt, ARRAY_SIZE(xt));
            ImPlot::SetupAxisTicks(ImAxis_Y1, yt, ARRAY_SIZE(yt));
            ImPlot::SetupFinish();

            ImDrawList* dl = ImPlot::GetPlotDrawList();
            ImPlot::PushPlotClipRect();
            // The regions, a run of cells of one kind at a time
            for (int t = 0; t < REGION_THETA; ++t) {
                const double t0 = t * 180.0 / REGION_THETA, t1 = (t + 1) * 180.0 / REGION_THETA;
                int f = 0;
                while (f < REGION_PHI) {
                    const md_pucker_kind_t k = regions6[t * REGION_PHI + f];
                    int e = f + 1;
                    while (e < REGION_PHI && regions6[t * REGION_PHI + e] == k) ++e;
                    const ImVec2 a = ImPlot::PlotToPixels(f * 360.0 / REGION_PHI, t0);
                    const ImVec2 b = ImPlot::PlotToPixels(e * 360.0 / REGION_PHI, t1);
                    dl->AddRectFilled(ImVec2(std::min(a.x, b.x), std::min(a.y, b.y)), ImVec2(std::max(a.x, b.x), std::max(a.y, b.y)), kind_color(k, 0.18f));
                    f = e;
                }
            }
            // The labels are 30 degrees apart along a row: smaller when the map is narrow, and every other one when
            // even that does not fit
            const ImU32 text = ImGui::GetColorU32(ImGuiCol_Text);
            const float step = ImPlot::PlotToPixels(30.0, 90.0).x - ImPlot::PlotToPixels(0.0, 90.0).x;
            const float need = ImGui::GetFontSize() * (sugars ? 1.9f : 1.2f);
            const float scale = std::clamp(step / need, 0.7f, 1.0f);
            const bool every_other = step < need * 0.7f;
            for (const md_pucker_canonical_t& k : canon6) {
                if (k.conf.kind == MD_PUCKER_KIND_CHAIR) {
                    // At a pole: along the whole edge
                    const double t = k.theta < 90.0 ? 8.0 : 172.0;
                    for (double f : { 60.0, 180.0, 300.0 }) draw_conformer(dl, ImPlot::PlotToPixels(f, t), k.conf, sugars ? loc6 : nullptr, text);
                    continue;
                }
                if (every_other && fmodf(k.phi2, 60.0f) != 0.0f) continue;
                draw_conformer(dl, ImPlot::PlotToPixels(k.phi2, k.theta), k.conf, sugars ? loc6 : nullptr, text, scale);
                if (k.phi2 == 0.0f) draw_conformer(dl, ImPlot::PlotToPixels(360.0, k.theta), k.conf, sugars ? loc6 : nullptr, text, scale);
            }
            ImPlot::PopPlotClipRect();
            plot_points(xs, ys, idx, prev_hovered);
            const int hit = pick(xs, ys, idx);
            if (hit >= 0) {
                hovered = hit;
                row_tooltip(rows[hit]);
            }
            ImPlot::EndPlot();
        }
        if (ImGui::IsItemHovered() && hovered < 0) {
            ImGui::SetTooltip("Six membered rings: %d\nChairs at the top and bottom (θ 0 and 180), boats and skews along θ 90%s",
                count, sugars ? "\nLabels: pyranose numbering, O for the ring O" : "\nC chair, H half chair, E envelope, B boat, S skew");
        }
    }

    // How many of the rings shown have each name
    void draw_names() {
        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;
        if (ImGui::BeginTable("names", 3, flags)) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Conformation", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Rings");
            ImGui::TableSetupColumn("%");
            ImGui::TableHeadersRow();
            for (const NameCount& n : names) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(n.name.c_str());
                ImGui::TableNextColumn(); ImGui::Text("%d", n.count);
                ImGui::TableNextColumn(); ImGui::Text("%.1f", 100.0 * n.count / std::max<size_t>(shown.size(), 1));
            }
            ImGui::EndTable();
        }
    }

    void draw_table(int prev_hovered) {
        enum Col { Col_Ring, Col_Residue, Col_Atoms, Col_Name, Col_Q, Col_Theta, Col_Phi, Col_P, Col_Tau, Col_Count };
        const ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_Sortable |
            ImGuiTableFlags_Resizable | ImGuiTableFlags_Hideable | ImGuiTableFlags_SizingFixedFit;
        if (!ImGui::BeginTable("rings", Col_Count, flags)) return;
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_DefaultSort);
        ImGui::TableSetupColumn("Residue");
        ImGui::TableSetupColumn("Atoms", ImGuiTableColumnFlags_DefaultHide);
        ImGui::TableSetupColumn("Conformation", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Q (Å)");
        ImGui::TableSetupColumn("θ (°)");
        ImGui::TableSetupColumn("φ2 (°)");
        ImGui::TableSetupColumn("P (°)");
        ImGui::TableSetupColumn("τm (°)");
        ImGui::TableHeadersRow();

        static std::vector<int> order;
        order = shown;
        if (ImGuiTableSortSpecs* spec = ImGui::TableGetSortSpecs(); spec && spec->SpecsCount > 0) {
            const int col = spec->Specs[0].ColumnIndex;
            const bool asc = spec->Specs[0].SortDirection == ImGuiSortDirection_Ascending;
            auto num = [&](int i) -> double {
                const md_pucker_ring_t& r = rows[i].r;
                switch (col) {
                case Col_Residue: return rows[i].res_idx;
                case Col_Q:     return r.cp.amplitude;
                case Col_Theta: return r.size == 6 ? r.cp.theta : -1.0;
                case Col_Phi:   return r.cp.phi2;
                case Col_P:     return r.has_pseudo ? r.pseudo_phase : -1.0;
                case Col_Tau:   return r.has_pseudo ? r.pseudo_amplitude : -1.0;
                default:        return rows[i].ring_idx;
                }
            };
            if (col == Col_Name) {
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                    const int c = strcmp(rows[a].r.name, rows[b].r.name);
                    return asc ? c < 0 : c > 0;
                });
            } else {
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return asc ? num(a) < num(b) : num(a) > num(b); });
            }
        }

        ImGuiListClipper clipper;
        clipper.Begin((int)order.size());
        while (clipper.Step()) {
            for (int line = clipper.DisplayStart; line < clipper.DisplayEnd; ++line) {
                const int i = order[line];
                const Row& row = rows[i];
                const md_pucker_ring_t& r = row.r;
                // Phases of a ring whose order is not fixed by its atoms depend on that order: shown, but dimmed
                const bool loose = !r.anchored && !r.sugar;
                char buf[96];
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                snprintf(buf, sizeof(buf), "%u", row.ring_idx + 1);
                ImGui::Selectable(buf, i == prev_hovered || is_highlighted(row), ImGuiSelectableFlags_SpanAllColumns);
                if (ImGui::IsItemHovered()) {
                    hovered = i;
                    row_tooltip(row);
                }
                ImGui::TableNextColumn(); residue_label(buf, sizeof(buf), row); ImGui::TextUnformatted(buf);
                ImGui::TableNextColumn(); atoms_label(buf, sizeof(buf), row);   ImGui::TextUnformatted(buf);
                ImGui::TableNextColumn(); ImGui::TextUnformatted(r.name[0] ? r.name : "-");
                ImGui::TableNextColumn();
                if (!row.valid) {
                    ImGui::TextDisabled("degenerate");
                    continue;
                }
                ImGui::Text("%.3f", r.cp.amplitude);
                ImGui::TableNextColumn();
                if (r.size == 6) {
                    if (loose) ImGui::TextDisabled("%.1f", r.cp.theta); else ImGui::Text("%.1f", r.cp.theta);
                }
                ImGui::TableNextColumn();
                if (loose) ImGui::TextDisabled("%.1f", r.cp.phi2); else ImGui::Text("%.1f", r.cp.phi2);
                ImGui::TableNextColumn();
                if (r.has_pseudo) ImGui::Text("%.1f", r.pseudo_phase);
                ImGui::TableNextColumn();
                if (r.has_pseudo) ImGui::Text("%.1f", r.pseudo_amplitude);
            }
        }
        ImGui::EndTable();
    }
};

static Pucker instance;

}  // namespace
