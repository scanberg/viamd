// HEATMAPS
//
// Matrices the script publishes (script::OUTPUT_MATRIX), dragged in from the script editor's gutter, one per tab: a
// contact map, contacts(residue), or any value of two dimensions.
//
// A matrix of truths per frame is kept as the frames each of its entries is true in (md_pair_frames.h). It is shown at
// the frame shown, or as the fraction of frames each entry is true in, over all of them or a timeline period's, as far
// as they are evaluated. A matrix of numbers per frame is shown at the frame shown, or as its mean over the frames
// evaluated. One over the trajectory, or static, is shown as it is.
//
//   hover              a cell: what its row and column are, its value, and for a fraction the frames it is true in;
//                      the row and the column are lit up in the view
//   a period dropped   from the timeline, the fractions or the mean over its frames
//   wheel, drag        zoom and pan; double-click fits
//
// What is drawn is at most 256 cells a side of what is in view, so zooming in shows every cell.

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_str.h>
#include <core/md_pair_frames.h>

#include <md_attributes.h>

#include <viamd_event.h>
#include <window_menu.h>
#include <viamd.h>
#include <event.h>
#include <serialization_utils.h>
#include <script_engine.h>
#include <plot_series.h>
#include <timeline_periods.h>
#include <timeline_periods_ui.h>
#include <display_units.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>

#include <vector>
#include <cmath>
#include <cstring>
#include <cstdio>

namespace {

enum Over : int {
    Over_Frame = 0,     // The frame shown
    Over_All,           // Every frame
    Over_Period,        // The frames of a timeline period, where they are now
    Over_Count,
};

static const md_str_t over_ini[Over_Count] = { MD_STR_LIT("frame"), MD_STR_LIT("all"), MD_STR_LIT("period") };

#define HEATMAP_MAX_CELLS 256

struct Heatmap {
    char     name[64] = "";
    Over     over = Over_All;
    uint32_t period = 0;            // Of Over_Period, by id
    bool     select = false;        // Its tab, next frame

    // What it shows: the whole matrix, row by row, NaN where nothing is known
    std::vector<float> values;
    uint32_t rows = 0, cols = 0;
    bool     pairs = false;         // A matrix of truths per frame, as pairs
    bool     temporal = false;
    bool     truth = false;         // Truths, of the frame shown
    bool     fraction = false;      // Fractions of frames, of truths
    size_t   num_frames = 0;        // Of a fraction or a mean: over how many frames
    md_unit_t unit = {};
    float    lo = 0.0f, hi = 1.0f;
    uint64_t key = 0;
    char     status[192] = "";
};

static float attr_value(const md_attribute_t* a, size_t i) {
    switch (a->format.type) {
    case MD_ATTRIBUTE_TYPE_F32: return ((const float*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_F64: return (float)((const double*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U8:  return (float)((const uint8_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I8:  return (float)((const int8_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U16: return (float)((const uint16_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I16: return (float)((const int16_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U32: return (float)((const uint32_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I32: return (float)((const int32_t*)a->data)[i];
    default: return NAN;
    }
}

// A tick's text: the label of the row or column it is the edge of, else its index
struct AxisLabels {
    const md_attributes_t* attrs = nullptr;
    const md_attribute_t*  labels = nullptr;
    uint32_t n = 0;
    bool     rows = false;          // The y axis, rows from the top
};

static int axis_label(double value, char* buf, int size, void* data) {
    const AxisLabels* a = (const AxisLabels*)data;
    const double at = a->rows ? (double)a->n - value : value;
    const long i = (long)floor(at + 1e-6);
    if (i < 0 || (uint32_t)i >= a->n) return snprintf(buf, size, "%s", "");
    if (a->labels) {
        const md_str_t s = md_attribute_str(a->attrs, a->labels, (size_t)i);
        if (s.len) return snprintf(buf, size, MD_STR_FMT, MD_STR_ARG(s));
    }
    return snprintf(buf, size, "%ld", i);
}

struct Heatmaps : viamd::EventHandler {
    bool show_window = false;
    ApplicationState* app = nullptr;
    std::vector<Heatmap> maps;

    // The cell under the pointer, this frame: whose rows and columns are lit up
    int  hovered_map = -1;
    int  hovered_row = -1, hovered_col = -1;
    bool highlighting = false;

    Heatmaps() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app = (ApplicationState*)e.payload;
                workspace_register_window("Heatmaps", &show_window);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                maps.clear();
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowMenuRequest:
                viamd::window_menu_add(e, "Analysis/Heatmaps", &show_window, {.order = 25});
                break;
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, MD_STR_LIT("Heatmaps"));
                for (const Heatmap& h : maps) {
                    viamd::write_str(state, MD_STR_LIT("Matrix"), md_str_from_cstr(h.name));
                    viamd::write_str(state, MD_STR_LIT("Over"), over_ini[h.over]);
                    if (h.over == Over_Period) viamd::write_int(state, MD_STR_LIT("Period"), (int)h.period);
                }
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (md_str_eq(viamd::section_header(state), MD_STR_LIT("Heatmaps"))) {
                    md_str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (md_str_eq(ident, MD_STR_LIT("Matrix"))) {
                            Heatmap h;
                            viamd::extract_to_char_buf(h.name, sizeof(h.name), arg);
                            maps.push_back(h);
                        } else if (maps.empty()) {
                            continue;
                        } else if (md_str_eq(ident, MD_STR_LIT("Over"))) {
                            md_str_t s = {};
                            viamd::extract_str(s, arg);
                            for (int o = 0; o < Over_Count; ++o) if (md_str_eq(s, over_ini[o])) maps.back().over = (Over)o;
                        } else if (md_str_eq(ident, MD_STR_LIT("Period"))) {
                            int id = 0;
                            viamd::extract_int(id, arg);
                            maps.back().period = (uint32_t)id;
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

    // A matrix of the script, by the path its mark drags ("script/<name>"); one already here is selected
    void add(const char* path) {
        const char* name = strncmp(path, "script/", 7) == 0 ? path + 7 : path;
        for (Heatmap& h : maps) {
            if (strcmp(h.name, name) == 0) { h.select = true; return; }
        }
        Heatmap h;
        snprintf(h.name, sizeof(h.name), "%s", name);
        h.select = true;
        maps.push_back(h);
        show_window = true;
    }

    // The frames [beg, end) of what h is over, false when it covers none
    bool frames_of(const Heatmap& h, uint32_t* beg, uint32_t* end) const {
        const int nf = (int)run_num_frames(app);
        if (nf <= 0) return false;
        switch (h.over) {
        case Over_Frame: {
            const int f = MD_CLAMP((int)(app->animation.frame + 0.5), 0, nf - 1);
            *beg = (uint32_t)f;
            *end = (uint32_t)f + 1;
            return true;
        }
        case Over_All:
            *beg = 0;
            *end = (uint32_t)nf;
            return true;
        case Over_Period: {
            const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, h.period);
            int b = 0, e = 0;
            if (!p || !timeline_period_frames(&b, &e, *p, app->animation.frame, nf)) return false;
            *beg = (uint32_t)b;
            *end = (uint32_t)e + 1;
            return true;
        }
        default:
            return false;
        }
    }

    // What h shows, made again when anything it is made from has changed, and every frame while evaluating
    void update(Heatmap& h) {
        const script::Eval* eval = app->script.eval;
        const script::Program* ir = app->script.eval_ir;
        const md_str_t name = md_str_from_cstr(h.name);
        const uint32_t flags = ir ? script::output_flags(ir, name) : 0;
        uint32_t beg = 0, end = 0;
        const bool have_frames = frames_of(h, &beg, &end);
        const bool running = app->tasks.evaluate.running();
        const uint64_t version = eval ? script::eval_errors(eval).version : 0;
        const int tick = running ? ImGui::GetFrameCount() : 0;
        uint64_t key = md_hash64(h.name, strlen(h.name), (uint64_t)(uintptr_t)eval);
        key = md_hash64(&flags, sizeof(flags), key);
        const uint32_t state[6] = { (uint32_t)h.over, h.period, have_frames ? beg : ~0u, end, (uint32_t)tick, (uint32_t)version };
        key = md_hash64(state, sizeof(state), key ^ version);
        if (key == h.key) return;
        h.key = key;

        h.values.clear();
        h.rows = h.cols = 0;
        h.pairs = h.temporal = h.truth = h.fraction = false;
        h.num_frames = 0;
        h.unit = {};
        h.status[0] = '\0';
        if (!(flags & script::OUTPUT_MATRIX)) {
            snprintf(h.status, sizeof(h.status), "The script publishes no matrix '%s'", h.name);
            return;
        }
        h.temporal = (flags & script::OUTPUT_TEMPORAL) != 0;
        if (!eval) {
            snprintf(h.status, sizeof(h.status), "Not evaluated: Evaluate in the script editor");
            return;
        }
        const md_attributes_t* attrs = script::attributes(eval);
        char path[SERIES_PATH_CAP];
        snprintf(path, sizeof(path), "script/%s", h.name);
        const md_attribute_t* attr = attrs ? md_attributes_find(attrs, md_str_from_cstr(path)) : nullptr;

        if (const md_pair_frames_t* pf = script::eval_pairs(eval, name)) {
            h.pairs = true;
            h.rows = md_pair_frames_rows(pf);
            h.cols = md_pair_frames_cols(pf);
            const size_t n = (size_t)h.rows * h.cols;
            h.values.assign(n, NAN);
            h.lo = 0.0f;
            h.hi = 1.0f;
            if (!have_frames) {
                snprintf(h.status, sizeof(h.status), "The period covers none of the frames");
            } else if (h.over == Over_Frame) {
                h.truth = true;
                std::vector<uint8_t> m(n);
                if (md_pair_frames_frame(m.data(), pf, beg)) {
                    size_t ones = 0;
                    for (size_t i = 0; i < n; ++i) { h.values[i] = m[i] ? 1.0f : 0.0f; ones += m[i]; }
                    h.num_frames = 1;
                    snprintf(h.status, sizeof(h.status), "Frame %u: %zu entries true", beg, ones);
                } else {
                    snprintf(h.status, sizeof(h.status), "Frame %u is not evaluated yet", beg);
                }
            } else {
                h.fraction = true;
                h.num_frames = md_pair_frames_occupancy(h.values.data(), pf, md_bitset_range(beg, end));
                if (h.num_frames == 0) snprintf(h.status, sizeof(h.status), "None of its %u frames is evaluated yet", end - beg);
                else snprintf(h.status, sizeof(h.status), "The fraction of %zu frames%s each entry is true in, of %zu entries ever true",
                              h.num_frames, h.num_frames < end - beg ? " evaluated so far" : "", md_pair_frames_count(pf));
            }
            return;
        }

        if (!attr || !attr->data || (attr->format.rank != 2 && attr->format.rank != 3)) {
            snprintf(h.status, sizeof(h.status), running ? "Being evaluated" : "Not evaluated: Evaluate in the script editor");
            return;
        }
        const bool per_frame = attr->format.rank == 3;
        h.rows = attr->format.shape[per_frame ? 1 : 0];
        h.cols = attr->format.shape[per_frame ? 2 : 1];
        h.unit = attr->unit;
        const size_t n = (size_t)h.rows * h.cols;
        h.values.assign(n, NAN);
        if (!per_frame) {
            for (size_t i = 0; i < n; ++i) h.values[i] = attr_value(attr, i);
        } else if (!have_frames) {
            snprintf(h.status, sizeof(h.status), "The period covers none of the frames");
        } else {
            const uint32_t nf = attr->format.shape[0];
            end = MD_MIN(end, nf);
            // The frames it holds: a frame not computed yet is NaN, or 0 for integers, and is left out
            md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
            const md_bitset_t held = script::eval_output_frames(eval, name, temp);
            const uint64_t num = beg < end ? md_bitset_count_range(held, beg, end) : 0;
            if (num == 0) {
                snprintf(h.status, sizeof(h.status), h.over == Over_Frame ? "Frame %u is not evaluated yet" : "None of its frames is evaluated yet", beg);
            } else if (n * num > ((size_t)1 << 27)) {
                snprintf(h.status, sizeof(h.status), "Too many values to average here (%zu × %llu frames): mean(%s, axis='frame') in the script",
                         n, (unsigned long long)num, h.name);
            } else {
                std::vector<uint32_t> count(n, 0);
                std::vector<double> sum(n, 0.0);
                md_bitset_iter_t it = md_bitset_iter_range(held, beg, end);
                while (md_bitset_iter_next(&it)) {
                    const size_t base = (size_t)it.idx * n;
                    for (size_t i = 0; i < n; ++i) {
                        const float v = attr_value(attr, base + i);
                        if (std::isnan(v)) continue;
                        sum[i] += v;
                        count[i] += 1;
                    }
                }
                for (size_t i = 0; i < n; ++i) h.values[i] = count[i] ? (float)(sum[i] / count[i]) : NAN;
                h.num_frames = (size_t)num;
                if (h.over == Over_Frame) snprintf(h.status, sizeof(h.status), "Frame %u", beg);
                else snprintf(h.status, sizeof(h.status), "The mean over %llu frames", (unsigned long long)num);
            }
            md_arena_allocator_destroy(temp);
        }
        // Its range, in the display unit
        const double f = display_units::factor(nullptr, h.unit);
        float lo = INFINITY, hi = -INFINITY;
        for (float& v : h.values) {
            if (std::isnan(v)) continue;
            v = (float)(v * f);
            lo = MD_MIN(lo, v);
            hi = MD_MAX(hi, v);
        }
        if (!(hi >= lo)) { lo = 0.0f; hi = 1.0f; }
        if (hi == lo) hi = lo + 1.0f;
        h.lo = lo;
        h.hi = hi;
    }

    // The value of a cell as the tooltip says it
    void cell_text(char* buf, size_t cap, const Heatmap& h, float v) const {
        if (std::isnan(v)) { snprintf(buf, cap, "Not known"); return; }
        if (h.truth) { snprintf(buf, cap, v > 0.5f ? "True" : "False"); return; }
        if (h.fraction) {
            const size_t k = (size_t)llround((double)v * (double)h.num_frames);
            snprintf(buf, cap, "True in %zu of %zu frames (%.0f%%)", k, h.num_frames, 100.0 * v);
            return;
        }
        char unit[64] = "";
        display_units::factor_print(unit, sizeof(unit), h.unit);
        snprintf(buf, cap, "%.4g%s%s", v, unit[0] ? " " : "", unit);
    }

    void draw_map(Heatmap& h, int idx) {
        const md_attributes_t* attrs = app->script.eval ? script::attributes(app->script.eval) : nullptr;
        // What it is over: the frame shown, every frame, a period
        if (h.temporal) {
            const TimelinePeriod* cur = h.over == Over_Period ? timeline_periods_find(&app->timeline.periods, h.period) : nullptr;
            const char* preview = h.over == Over_Frame ? "Frame shown" : h.over == Over_All ? "All frames" : cur ? cur->name : "No period";
            ImGui::SetNextItemWidth(200);
            if (ImGui::BeginCombo("Over", preview)) {
                if (ImGui::Selectable("Frame shown", h.over == Over_Frame)) h.over = Over_Frame;
                if (ImGui::Selectable("All frames", h.over == Over_All)) h.over = Over_All;
                const TimelinePeriods& ps = app->timeline.periods;
                if (ps.count) ImGui::Separator();
                for (int i = 0; i < ps.count; ++i) {
                    ImGui::PushID(i);
                    if (ImGui::Selectable(ps.items[i].name, h.over == Over_Period && h.period == ps.items[i].id)) {
                        h.over = Over_Period;
                        h.period = ps.items[i].id;
                    }
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip(h.pairs ? "Truths at the frame shown, or the fraction of frames each entry is true in.\nA period dropped here from the timeline: over its frames."
                                          : "The values at the frame shown, or their mean over frames.\nA period dropped here from the timeline: over its frames.");
            ImGui::SameLine();
        }
        ImGui::TextDisabled("%s", h.status);
        if (h.values.empty() || !h.rows || !h.cols) return;

        char row_path[SERIES_PATH_CAP + 16], col_path[SERIES_PATH_CAP + 16];
        snprintf(row_path, sizeof(row_path), "script/%s/labels", h.name);
        snprintf(col_path, sizeof(col_path), "script/%s/column_labels", h.name);
        const md_attribute_t* rl = attrs ? md_attributes_find(attrs, md_str_from_cstr(row_path)) : nullptr;
        const md_attribute_t* cl = attrs ? md_attributes_find(attrs, md_str_from_cstr(col_path)) : nullptr;
        if (rl && (rl->format.type != MD_ATTRIBUTE_TYPE_STR || rl->format.shape[0] != h.rows)) rl = nullptr;
        if (cl && (cl->format.type != MD_ATTRIBUTE_TYPE_STR || cl->format.shape[0] != h.cols)) cl = nullptr;
        AxisLabels ax = { attrs, cl, h.cols, false }, ay = { attrs, rl, h.rows, true };

        const float scale_w = 70.0f;
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const ImVec2 size(MD_MAX(avail.x - scale_w - ImGui::GetStyle().ItemSpacing.x, 100.0f), MD_MAX(avail.y, 100.0f));
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        char id[96];
        snprintf(id, sizeof(id), "##heatmap_%s", h.name);
        if (ImPlot::BeginPlot(id, size, ImPlotFlags_NoLegend | ImPlotFlags_NoTitle | ImPlotFlags_NoMouseText | ImPlotFlags_Equal)) {
            ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoGridLines);
            ImPlot::SetupAxisFormat(ImAxis_X1, axis_label, &ax);
            ImPlot::SetupAxisFormat(ImAxis_Y1, axis_label, &ay);
            ImPlot::SetupAxesLimits(0, h.cols, 0, h.rows, ImPlotCond_Once);
            // What is in view, at most HEATMAP_MAX_CELLS cells a side: blocks of cells as one, a truth where any of
            // them is true, else their mean
            const ImPlotRect lim = ImPlot::GetPlotLimits();
            const uint32_t c0 = (uint32_t)MD_CLAMP(floor(lim.X.Min), 0.0, (double)h.cols), c1 = (uint32_t)MD_CLAMP(ceil(lim.X.Max), 0.0, (double)h.cols);
            const uint32_t r0 = (uint32_t)MD_CLAMP(floor((double)h.rows - lim.Y.Max), 0.0, (double)h.rows), r1 = (uint32_t)MD_CLAMP(ceil((double)h.rows - lim.Y.Min), 0.0, (double)h.rows);
            if (c1 > c0 && r1 > r0) {
                const uint32_t bc = (c1 - c0 + HEATMAP_MAX_CELLS - 1) / HEATMAP_MAX_CELLS, br = (r1 - r0 + HEATMAP_MAX_CELLS - 1) / HEATMAP_MAX_CELLS;
                const uint32_t nc = (c1 - c0 + bc - 1) / bc, nr = (r1 - r0 + br - 1) / br;
                std::vector<float> cells((size_t)nr * nc);
                for (uint32_t r = 0; r < nr; ++r) {
                    for (uint32_t c = 0; c < nc; ++c) {
                        double sum = 0.0, top = -INFINITY;
                        uint32_t count = 0;
                        for (uint32_t i = r0 + r * br; i < MD_MIN(r0 + (r + 1) * br, r1); ++i) {
                            for (uint32_t j = c0 + c * bc; j < MD_MIN(c0 + (c + 1) * bc, c1); ++j) {
                                const float v = h.values[(size_t)i * h.cols + j];
                                if (std::isnan(v)) continue;
                                sum += v;
                                top = MD_MAX(top, (double)v);
                                ++count;
                            }
                        }
                        // Nothing known: below the scale, drawn as its least
                        cells[(size_t)r * nc + c] = count == 0 ? h.lo : (float)(h.truth ? top : sum / count);
                    }
                }
                ImPlot::PlotHeatmap("##cells", cells.data(), (int)nr, (int)nc, h.lo, h.hi, nullptr,
                                    ImPlotPoint(c0, (double)h.rows - (r0 + (double)nr * br)), ImPlotPoint(c0 + (double)nc * bc, (double)h.rows - r0));
            }

            if (h.temporal && ImPlot::BeginDragDropTargetPlot()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                    if (payload->DataSize == (int)sizeof(uint32_t)) {
                        h.over = Over_Period;
                        h.period = *(const uint32_t*)payload->Data;
                    }
                }
                ImPlot::EndDragDropTarget();
            }

            if (ImPlot::IsPlotHovered()) {
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                const long c = (long)floor(mp.x), r = (long)floor((double)h.rows - mp.y);
                if (c >= 0 && r >= 0 && (uint32_t)c < h.cols && (uint32_t)r < h.rows) {
                    hovered_map = idx;
                    hovered_row = (int)r;
                    hovered_col = (int)c;
                    const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint((double)c, (double)h.rows - r));
                    const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint((double)c + 1, (double)h.rows - r - 1));
                    ImPlot::GetPlotDrawList()->AddRect(p0, p1, IM_COL32(255, 255, 255, 220), 0.0f, 0, 1.5f);
                    draw_cell_tooltip(h, (uint32_t)r, (uint32_t)c, rl, cl, attrs);
                }
            }
            ImPlot::EndPlot();
        }
        ImGui::SameLine();
        ImPlot::ColormapScale("##scale", h.lo, h.hi, ImVec2(scale_w, size.y), h.truth || h.fraction ? "%.2f" : "%.3g");
        ImPlot::PopColormap();
    }

    void draw_cell_tooltip(const Heatmap& h, uint32_t r, uint32_t c, const md_attribute_t* rl, const md_attribute_t* cl, const md_attributes_t* attrs) {
        const md_str_t rs = rl ? md_attribute_str(attrs, rl, r) : md_str_t{};
        const md_str_t cs = cl ? md_attribute_str(attrs, cl, c) : md_str_t{};
        char text[128];
        cell_text(text, sizeof(text), h, h.values[(size_t)r * h.cols + c]);
        ImGui::BeginTooltip();
        if (rs.len) ImGui::Text(MD_STR_FMT, MD_STR_ARG(rs)); else ImGui::Text("Row %u", r);
        ImGui::SameLine();
        ImGui::TextDisabled("and");
        ImGui::SameLine();
        if (cs.len) ImGui::Text(MD_STR_FMT, MD_STR_ARG(cs)); else ImGui::Text("column %u", c);
        ImGui::TextUnformatted(text);
        // A fraction: when, as a strip over the frames, the run's start at the left
        const md_pair_frames_t* pf = h.fraction && app->script.eval ? script::eval_pairs(app->script.eval, md_str_from_cstr(h.name)) : nullptr;
        if (pf) {
            md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
            const md_bitset_t when = md_pair_frames_entry(pf, r, c, temp);
            const uint32_t nf = md_pair_frames_num_frames(pf);
            const float w = 18.0f * ImGui::GetFontSize(), ht = 0.8f * ImGui::GetFontSize();
            const ImVec2 o = ImGui::GetCursorScreenPos();
            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddRectFilled(o, ImVec2(o.x + w, o.y + ht), IM_COL32(60, 60, 60, 255));
            const int px = (int)w;
            for (int x = 0; x < px; ++x) {
                const uint32_t b0 = (uint32_t)((uint64_t)x * nf / px), b1 = MD_MAX((uint32_t)((uint64_t)(x + 1) * nf / px), b0 + 1);
                const uint64_t k = md_bitset_count_range(when, b0, b1);
                if (!k) continue;
                const float a = (float)k / (float)(b1 - b0);
                dl->AddRectFilled(ImVec2(o.x + x, o.y), ImVec2(o.x + x + 1, o.y + ht), ImGui::GetColorU32(ImVec4(0.99f, 0.91f, 0.14f, 0.25f + 0.75f * a)));
            }
            ImGui::Dummy(ImVec2(w, ht));
            ImGui::TextDisabled("When, over the %u frames", nf);
            md_arena_allocator_destroy(temp);
        }
        ImGui::EndTooltip();
    }

    void draw_window() {
        hovered_map = hovered_row = hovered_col = -1;
        if (show_window && app) {
            ImGui::SetNextWindowSize(ImVec2(640, 640), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Heatmaps", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
                if (maps.empty()) {
                    ImGui::TextDisabled("Drag a matrix here from the gutter of the script editor:\ncontacts(residue), or any value of two dimensions");
                } else if (ImGui::BeginTabBar("##heatmaps", ImGuiTabBarFlags_Reorderable | ImGuiTabBarFlags_FittingPolicyScroll)) {
                    int remove = -1;
                    for (int i = 0; i < (int)maps.size(); ++i) {
                        Heatmap& h = maps[i];
                        bool open = true;
                        const ImGuiTabItemFlags tf = h.select ? ImGuiTabItemFlags_SetSelected : 0;
                        h.select = false;
                        if (ImGui::BeginTabItem(h.name, &open, tf)) {
                            update(h);
                            draw_map(h, i);
                            ImGui::EndTabItem();
                        }
                        if (!open) remove = i;
                    }
                    ImGui::EndTabBar();
                    if (remove >= 0) maps.erase(maps.begin() + remove);
                }
                // A matrix dropped anywhere in the window
                ImGuiWindow* win = ImGui::GetCurrentWindow();
                if (ImGui::BeginDragDropTargetCustom(win->InnerRect, win->ID)) {
                    if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(HEATMAP_SERIES_DND)) {
                        if (payload->DataSize == (int)sizeof(SeriesDragPayload)) add(((const SeriesDragPayload*)payload->Data)->key.path);
                    }
                    ImGui::EndDragDropTarget();
                }
            }
            ImGui::End();
        }

        // The row and the column of the cell under the pointer, lit up in the view
        if (hovered_map >= 0 && hovered_map < (int)maps.size() && app->script.eval_ir) {
            const Heatmap& h = maps[hovered_map];
            const md_str_t name = md_str_from_cstr(h.name);
            md_bitmask_t* mask = &app->selection.highlight_mask;
            md_bitmask_clear(mask);
            std::vector<uint32_t> atoms;
            for (int d = 0; d < 2; ++d) {
                const uint32_t index = (uint32_t)(d ? hovered_col : hovered_row);
                const size_t n = script::output_atoms(app->script.eval_ir, name, d, index, nullptr, 0);
                atoms.resize(n);
                script::output_atoms(app->script.eval_ir, name, d, index, atoms.data(), n);
                for (uint32_t a : atoms) md_bitmask_set_bit(mask, a);
            }
            highlighting = true;
        } else if (highlighting) {
            md_bitmask_clear(&app->selection.highlight_mask);
            highlighting = false;
        }
    }
};

static Heatmaps instance;

}  // namespace
