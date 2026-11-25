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
// It is drawn in tiles of 128 × 128 cells, a cell 2^level entries a side, at the level that makes a cell a pixel or a
// little more, so what is drawn is never more than the screen shows. A cell of many entries is the highest of them for
// truths and fractions, so a contact is not lost from far away, and their mean for other numbers.
//
// Nothing of it is made on the main thread, which only draws and uploads. What a heatmap shows is a View, made by a job
// whenever what it is over changes: for a matrix of truths, the frames counted and an index of where its entries are
// kept, so that a tile costs what is true in it (heatmaps_core.h); for any other, the matrix. Tiles in view that are
// not of the view shown are then made by another job, over the threads, and a large tile's cell rows over them too;
// each is uploaded as a texture as soon as it is done, and kept, those drawn least recently going when there are many.
// One job runs at a time, and a View or its index is replaced only when none does. Until a tile is made it is drawn as
// it was, or from a coarser tile, or as a blank, with a spinner once it has waited a while.

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_pair_frames.h>

#include <md_attributes.h>

#include <viamd_event.h>
#include <window_menu.h>
#include <viamd.h>
#include <system_scene.h>
#include <event.h>
#include <serialization_utils.h>
#include <script_engine.h>
#include <plot_series.h>
#include <timeline_periods.h>
#include <timeline_periods_ui.h>
#include <display_units.h>
#include <task_system.h>
#include "heatmaps_core.h"
#include <gfx/gl.h>
#include <gfx/gl_utils.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>
#include <plot_input.h>
#include <implot_internal.h>

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
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

using heatmaps::TILE;
using heatmaps::Index;
using heatmaps::View;
using heatmaps::ViewSource;

#define TILES_PER_JOB   256     // Tiles made by one job, those nearest the middle of the view first
#define TILE_KEEP       512     // Tiles kept, or twice those in view when more
#define SPINNER_AFTER   0.15    // Seconds a tile waits before it shows that it does

struct Tile {
    GLuint   tex = 0;               // TILE × TILE, RGBA8, its first row at the top
    uint64_t gen = 0;               // Of the view it shows
    int      used = 0;              // The frame it was last drawn in
    double   waiting = 0;           // Since when it has been in view and not of the view shown, 0 when it is
};

// A tile a job makes: its cells, colored, then done
struct TileWork {
    uint64_t key = 0;
    uint32_t level = 0, tx = 0, ty = 0;
    std::vector<uint32_t> rgba;
    std::atomic<int> done{0};
};

enum Work : int { Work_None = 0, Work_View, Work_Tiles };

struct Heatmap {
    char     name[64] = "";
    // Of the evaluation of the run of which instance of the first system (by id, instance_resolve): its frames, the
    // frame it shows. The script is of the first system.
    uint32_t instance = 0;
    const script::Eval* eval = nullptr;     // That evaluation, as last read: what its view and labels are of
    Over     over = Over_All;
    uint32_t period = 0;            // Of Over_Period, by id
    bool     select = false;        // Its tab, next frame
    bool     temporal = false;
    bool     blank = false;         // Nothing to show now: status says why
    char     status[192] = "";      // Why, or what is being made
    uint64_t want = 0;              // What it should show, as a key
    uint64_t made = 0;              // What the view shown, or the one being made, is of
    uint64_t gen = 0;

    // Shown; and made by a view job, taken in when it has finished
    std::unique_ptr<View>  view, next_view;
    std::unique_ptr<Index> index, next_index;
    uint32_t rows = 0, cols = 0;    // Of the tiles kept
    uint32_t lut[256] = {};         // The colormap, for the tile job

    // Labels of rows and columns, read where the evaluation's attributes are not being changed
    std::vector<std::string> row_labels, col_labels;
    bool labels_stale = true;

    std::unordered_map<uint64_t, Tile> tiles;   // By tile_key
    std::vector<uint64_t> wanted;               // In view and not of the view shown, this frame, nearest the middle first
    std::unique_ptr<TileWork[]> work;
    size_t   num_work = 0;
    Work     doing = Work_None;     // What job runs
    bool     from_attr = false;     // A view job reading the evaluation's attributes
    Work     last = Work_None;      // And the one that ran before: a view and its tiles take turns, so that a view
                                    // made again at every frame (the frame shown, while playing) has its tiles made
    task_system::Job job;           // Last: it goes first, and nothing it reads goes before it

    Heatmap() = default;
    Heatmap(const Heatmap&) = delete;
    Heatmap& operator=(const Heatmap&) = delete;
    ~Heatmap() {
        job.cancel_and_wait();
        for (auto& kv : tiles) gl::free_texture(&kv.second.tex);
    }
};

static inline uint64_t tile_key(uint32_t level, uint32_t tx, uint32_t ty) {
    return (uint64_t)level << 58 | (uint64_t)ty << 29 | tx;
}
static inline uint32_t key_level(uint64_t k) { return (uint32_t)(k >> 58); }
static inline uint32_t key_ty(uint64_t k) { return (uint32_t)(k >> 29) & ((1u << 29) - 1); }
static inline uint32_t key_tx(uint64_t k) { return (uint32_t)k & ((1u << 29) - 1); }

static void free_tiles(Heatmap& h) {
    for (auto& kv : h.tiles) gl::free_texture(&kv.second.tex);
    h.tiles.clear();
}

// A tick's text: the label of the row or column it is the edge of, else its index
struct AxisLabels {
    const std::vector<std::string>* labels = nullptr;
    uint32_t n = 0;
    bool     rows = false;          // The y axis, rows from the top
};

static int axis_label(double value, char* buf, int size, void* data) {
    const AxisLabels* a = (const AxisLabels*)data;
    const double at = a->rows ? (double)a->n - value : value;
    const long i = (long)floor(at + 1e-6);
    if (i < 0 || (uint32_t)i >= a->n) return snprintf(buf, size, "%s", "");
    if (a->labels && (size_t)i < a->labels->size() && !(*a->labels)[i].empty()) return snprintf(buf, size, "%s", (*a->labels)[i].c_str());
    return snprintf(buf, size, "%ld", i);
}

// Waiting, drawn over a tile
static void spinner(ImDrawList* dl, ImVec2 c, float r, double t) {
    const float a = (float)fmod(t * 5.0, 2.0 * 3.14159265358979);
    dl->PathArcTo(c, r, a, a + 4.5f, 20);
    dl->PathStroke(IM_COL32(235, 235, 235, 210), 0, MD_MAX(1.5f, r * 0.22f));
}

struct Heatmaps : viamd::EventHandler {
    bool show_window = false;
    ApplicationState* app = nullptr;
    std::vector<std::unique_ptr<Heatmap>> maps;

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
            case viamd::EventType_ViamdShutdown:
            case viamd::EventType_ViamdDeserializeBegin:
                maps.clear();
                break;
            case viamd::EventType_ViamdScriptEvalFree:
                // Nothing of it is read from now on: its jobs stopped, and the views of it gone. Tiles stay, drawn as
                // they were until those of the next evaluation are made.
                for (auto& h : maps) drop_views(*h);
                break;
            case viamd::EventType_ViamdScriptEvalBegin:
                // A view made from its attributes stops: they change while it runs
                for (auto& h : maps) {
                    if (h->doing != Work_View || !h->from_attr || h->eval != app->script.evaluating) continue;
                    h->job.cancel_and_wait();
                    h->made = 0;
                }
                break;
            case viamd::EventType_ViamdScriptEvalInit:
                for (auto& h : maps) {
                    h->made = 0;
                    h->labels_stale = true;
                }
                break;
            case viamd::EventType_ViamdScriptEvalEnd:
                for (auto& h : maps) {
                    if (h->eval != app->script.evaluating) continue;
                    h->made = 0;
                    h->labels_stale = true;
                }
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
                for (const auto& h : maps) {
                    viamd::write_str(state, MD_STR_LIT("Matrix"), md_str_from_cstr(h->name));
                    // Of another instance than the first
                    const SystemInstance& inst = instance_resolve(app->system_context, h->instance);
                    if (&inst != &system_instance(app->system_context)) viamd::write_int(state, MD_STR_LIT("Instance"), (int)inst.id);
                    viamd::write_str(state, MD_STR_LIT("Over"), over_ini[h->over]);
                    if (h->over == Over_Period) viamd::write_int(state, MD_STR_LIT("Period"), (int)h->period);
                }
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (md_str_eq(viamd::section_header(state), MD_STR_LIT("Heatmaps"))) {
                    md_str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (md_str_eq(ident, MD_STR_LIT("Matrix"))) {
                            maps.push_back(std::make_unique<Heatmap>());
                            viamd::extract_to_char_buf(maps.back()->name, sizeof(maps.back()->name), arg);
                        } else if (maps.empty()) {
                            continue;
                        } else if (md_str_eq(ident, MD_STR_LIT("Instance"))) {
                            int id = 0;
                            if (viamd::extract_int(id, arg)) maps.back()->instance = (uint32_t)MD_MAX(id, 0);
                        } else if (md_str_eq(ident, MD_STR_LIT("Over"))) {
                            md_str_t s = {};
                            viamd::extract_str(s, arg);
                            for (int o = 0; o < Over_Count; ++o) if (md_str_eq(s, over_ini[o])) maps.back()->over = (Over)o;
                        } else if (md_str_eq(ident, MD_STR_LIT("Period"))) {
                            int id = 0;
                            viamd::extract_int(id, arg);
                            maps.back()->period = (uint32_t)id;
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
        for (auto& h : maps) {
            if (strcmp(h->name, name) == 0) { h->select = true; return; }
        }
        maps.push_back(std::make_unique<Heatmap>());
        Heatmap& h = *maps.back();
        snprintf(h.name, sizeof(h.name), "%s", name);
        h.select = true;
        show_window = true;
    }

    // Its job stopped and waited for, and every view and index dropped: what they read is about to go
    static void drop_views(Heatmap& h) {
        h.job.cancel_and_wait();
        h.eval = nullptr;
        h.work.reset();
        h.num_work = 0;
        h.next_view.reset();
        h.next_index.reset();
        h.view.reset();
        h.index.reset();
        h.made = 0;
        h.labels_stale = true;
    }

    const SystemInstance& instance_of(const Heatmap& h) const { return instance_resolve(app->system_context, h.instance); }

    // The evaluation of the run of its instance: the first instance's run's for one without a run
    const script::Eval* eval_of(const Heatmap& h) const {
        const SystemRun* run = instance_run(instance_of(h));
        return script_eval_of_run(app, run ? md_str_from_cstr(run->path) : md_str_t{});
    }

    // The frames [beg, end) of what h is over, false when it covers none
    bool frames_of(const Heatmap& h, uint32_t* beg, uint32_t* end) const {
        const SystemInstance& inst = instance_of(h);
        const SystemRun* run = instance_run(inst);
        const int nf = (int)run_num_frames(app->system_context, run);
        if (nf <= 0) return false;
        switch (h.over) {
        case Over_Frame: {
            const int f = MD_CLAMP((int)(instance_current_frame(inst) + 0.5), 0, nf - 1);
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
            if (!p || !run_period_frames(&b, &e, app, app->system_context, run, *p) || b >= nf) return false;
            e = MD_MIN(e, nf - 1);
            *beg = (uint32_t)b;
            *end = (uint32_t)e + 1;
            return true;
        }
        default:
            return false;
        }
    }

    // The labels of its rows and columns, copied: read while the evaluation's attributes do not change
    void read_labels(Heatmap& h) {
        h.labels_stale = false;
        h.row_labels.clear();
        h.col_labels.clear();
        const md_attributes_t* attrs = h.eval ? script::attributes(h.eval) : nullptr;
        if (!attrs) return;
        for (int d = 0; d < 2; ++d) {
            char path[SERIES_PATH_CAP + 16];
            snprintf(path, sizeof(path), d ? "script/%s/column_labels" : "script/%s/labels", h.name);
            const md_attribute_t* a = md_attributes_find(attrs, md_str_from_cstr(path));
            if (!a || a->format.type != MD_ATTRIBUTE_TYPE_STR) continue;
            std::vector<std::string>& out = d ? h.col_labels : h.row_labels;
            out.resize(a->format.shape[0]);
            for (uint32_t k = 0; k < a->format.shape[0]; ++k) {
                const md_str_t s = md_attribute_str(attrs, a, k);
                out[k].assign(s.ptr ? s.ptr : "", s.len);
            }
        }
    }

    // What h should show, as a key, and a job to make its view when that is not what it shows. On the main thread
    // nothing is computed here: what is read is what a view is made from.
    void request(Heatmap& h) {
        const script::Eval* eval = eval_of(h);
        // Of another evaluation (another instance, or a run shown anew): its labels too
        if (eval != h.eval) {
            h.eval = eval;
            h.labels_stale = true;
        }
        const script::Program* ir = app->script.eval_ir;
        const md_str_t name = md_str_from_cstr(h.name);
        const uint32_t flags = ir ? script::output_flags(ir, name) : 0;
        const bool running = eval && app->script.evaluating == eval;
        h.temporal = (flags & script::OUTPUT_TEMPORAL) != 0;
        h.blank = false;
        h.status[0] = '\0';
        if (h.labels_stale && !running) read_labels(h);
        if (!(flags & script::OUTPUT_MATRIX)) {
            h.blank = true;
            snprintf(h.status, sizeof(h.status), "The script publishes no matrix '%s'", h.name);
            return;
        }
        if (!eval) {
            h.blank = true;
            snprintf(h.status, sizeof(h.status), "Not evaluated: Evaluate in the script editor");
            return;
        }
        uint32_t beg = 0, end = 0;
        const bool have_frames = frames_of(h, &beg, &end);
        ViewSource src;
        snprintf(src.name, sizeof(src.name), "%s", h.name);
        src.frame = h.over == Over_Frame;
        src.pf = script::eval_pairs(eval, name);
        if (src.pf) {
            if (!have_frames) {
                h.blank = true;
                snprintf(h.status, sizeof(h.status), "The period covers none of the frames");
                return;
            }
        } else {
            // From its attribute, which is not read while the evaluation runs: what was made stays
            if (running) {
                h.blank = !h.view;
                snprintf(h.status, sizeof(h.status), "Being evaluated");
                return;
            }
            char path[SERIES_PATH_CAP];
            snprintf(path, sizeof(path), "script/%s", h.name);
            const md_attributes_t* attrs = script::attributes(eval);
            const md_attribute_t* attr = attrs ? md_attributes_find(attrs, md_str_from_cstr(path)) : nullptr;
            if (!attr || !attr->data || (attr->format.rank != 2 && attr->format.rank != 3)) {
                h.blank = true;
                snprintf(h.status, sizeof(h.status), "Not evaluated: Evaluate in the script editor");
                return;
            }
            if (attr->format.rank == 3 && !have_frames) {
                h.blank = true;
                snprintf(h.status, sizeof(h.status), "The period covers none of the frames");
                return;
            }
            src.attr = *attr;
            src.eval = eval;
            src.factor = display_units::factor(nullptr, attr->unit);
        }
        src.beg = beg;
        src.end = end;

        // While evaluating, made again four times a second
        const uint64_t version = script::eval_errors(eval).version;
        const int tick = running ? (int)(ImGui::GetTime() * 4.0) : 0;
        uint64_t key = md_hash64(h.name, strlen(h.name), (uint64_t)(uintptr_t)eval ^ ((uint64_t)(uintptr_t)src.pf << 1));
        const uint32_t state[7] = { flags, (uint32_t)h.over, h.period, have_frames ? beg : ~0u, end, (uint32_t)tick, (uint32_t)version };
        key = md_hash64(state, sizeof(state), key ^ version);
        key = md_hash64(&src.factor, sizeof(src.factor), key);
        h.want = key ? key : 1;
        if (h.want != h.made && !h.job.running() && !(h.last == Work_View && !h.wanted.empty())) start_view(h, src);
    }

    void start_view(Heatmap& h, const ViewSource& src) {
        h.made = h.want;
        h.next_view.reset(new View);
        h.next_view->gen = ++h.gen;
        h.next_index.reset();
        h.doing = Work_View;
        h.from_attr = !src.pf;
        Heatmap* hp = &h;
        const Index* cur = h.index.get();
        h.job.start(MD_STR_LIT(""), [hp, cur, src]() { heatmaps::make_view(*hp->next_view, &hp->next_index, cur, src); },
                    [hp](task_system::Status st) { finish_view(*hp, st); });
    }

    // On the main thread, with no job running: the view made shown, and the index it is of kept
    static void finish_view(Heatmap& h, task_system::Status st) {
        h.doing = Work_None;
        h.last = Work_View;
        if (st != task_system::Status::Completed || !h.next_view) {
            h.next_view.reset();
            h.next_index.reset();
            h.made = 0;
            return;
        }
        if (h.next_view->rows != h.rows || h.next_view->cols != h.cols) {
            free_tiles(h);
            h.rows = h.next_view->rows;
            h.cols = h.next_view->cols;
        }
        h.view = std::move(h.next_view);
        if (h.next_index) h.index = std::move(h.next_index);
        else if (!h.view->pairs) h.index.reset();
    }

    // The tiles wanted, those nearest the middle first, made by a job of the view shown: unless another view is to be
    // made, and its turn has come
    void start_tiles(Heatmap& h) {
        if (h.job.running() || !h.view || h.wanted.empty() || (h.want != h.made && h.last != Work_View)) return;
        const size_t n = MD_MIN(h.wanted.size(), (size_t)TILES_PER_JOB);
        h.work.reset(new TileWork[n]);
        h.num_work = n;
        for (size_t k = 0; k < n; ++k) {
            TileWork& w = h.work[k];
            w.key = h.wanted[k];
            w.level = key_level(w.key);
            w.tx = key_tx(w.key);
            w.ty = key_ty(w.key);
        }
        for (int i = 0; i < 256; ++i) h.lut[i] = ImGui::ColorConvertFloat4ToU32(ImPlot::SampleColormap((float)i / 255.0f, ImPlotColormap_Viridis));
        h.doing = Work_Tiles;
        Heatmap* hp = &h;
        h.job.start(MD_STR_LIT(""), [hp]() {
            const View& v = *hp->view;
            task_system::parallel_for((uint32_t)hp->num_work, 1, [&](uint32_t beg, uint32_t end) {
                for (uint32_t k = beg; k < end && !task_system::cancelled(); ++k) {
                    TileWork& w = hp->work[k];
                    w.rgba.resize(TILE * TILE);
                    if (heatmaps::make_tile(w.rgba.data(), v, w.level, w.tx, w.ty, hp->lut, IM_COL32(90, 90, 90, 255))) w.done.store(1, std::memory_order_release);
                }
            });
        }, [hp](task_system::Status) {
            upload_done(*hp);
            hp->work.reset();
            hp->num_work = 0;
            hp->doing = Work_None;
            hp->last = Work_Tiles;
        });
    }

    // Tiles a job has made, uploaded: as soon as each is done
    static void upload_done(Heatmap& h) {
        if (!h.work || !h.view) return;
        for (size_t k = 0; k < h.num_work; ++k) {
            TileWork& w = h.work[k];
            if (w.done.load(std::memory_order_acquire) != 1) continue;
            w.done.store(2, std::memory_order_relaxed);
            Tile& t = h.tiles[w.key];
            if (!t.tex) {
                gl::init_texture_2D(&t.tex, TILE, TILE, GL_RGBA8);
                glBindTexture(GL_TEXTURE_2D, t.tex);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            }
            gl::set_texture_2D_data(t.tex, 0, w.rgba.data(), GL_RGBA8);
            glBindTexture(GL_TEXTURE_2D, 0);
            t.gen = h.view->gen;
            t.waiting = 0;
            w.rgba = {};
        }
    }

    // The value of a cell as the tooltip says it
    static void cell_text(char* buf, size_t cap, const View& v, float x) {
        if (std::isnan(x)) { snprintf(buf, cap, "Not known"); return; }
        if (v.truth) { snprintf(buf, cap, x > 0.5f ? "True" : "False"); return; }
        if (v.fraction) {
            const size_t k = (size_t)llround((double)x * (double)v.num_frames);
            snprintf(buf, cap, "True in %zu of %zu frames (%.0f%%)", k, v.num_frames, 100.0 * x);
            return;
        }
        char unit[64] = "";
        display_units::factor_print(unit, sizeof(unit), v.unit);
        snprintf(buf, cap, "%.4g%s%s", x, unit[0] ? " " : "", unit);
    }

    void draw_map(Heatmap& h, int idx) {
        const View* v = h.view.get();
        h.wanted.clear();       // Of what is drawn below
        // Of which instance's run, where the first system has several
        if (app->system_context.num_instances > 1) {
            ImGui::SetNextItemWidth(200);
            instance_combo("Instance", &h.instance, app->system_context);
            ImGui::SetItemTooltip("Of the evaluation of its run: its frames, and the frame it shows");
            ImGui::SameLine();
        }
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
            ImGui::SetItemTooltip(v && v->pairs ? "Truths at the frame shown, or the fraction of frames each entry is true in.\nA period dropped here from the timeline: over its frames."
                                                : "The values at the frame shown, or their mean over frames.\nA period dropped here from the timeline: over its frames.");
            ImGui::SameLine();
        }
        ImGui::TextDisabled("%s", h.status[0] ? h.status : v ? v->status : "Being made");
        if (h.blank || !h.rows || !h.cols) return;

        AxisLabels ax = { &h.col_labels, h.cols, false }, ay = { &h.row_labels, h.rows, true };
        const float scale_w = 70.0f;
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const ImVec2 size(MD_MAX(avail.x - scale_w - ImGui::GetStyle().ItemSpacing.x, 100.0f), MD_MAX(avail.y, 100.0f));
        const bool shown = v && v->shown();
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        char id[96];
        snprintf(id, sizeof(id), "##heatmap_%s", h.name);
        if (ImPlot::BeginPlot(id, size, ImPlotFlags_NoLegend | ImPlotFlags_NoTitle | ImPlotFlags_NoMouseText | ImPlotFlags_Equal)) {
            ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoGridLines);
            ImPlot::SetupAxisFormat(ImAxis_X1, axis_label, &ax);
            ImPlot::SetupAxisFormat(ImAxis_Y1, axis_label, &ay);
            ImPlot::SetupAxesLimits(0, h.cols, 0, h.rows, ImPlotCond_Once);
            upload_done(h);
            if (!v || shown) draw_tiles(h);
            // A fit, by double-click, is to the whole matrix: there is no item to fit to
            if (ImPlot::FitThisFrame()) {
                ImPlot::FitPoint(ImPlotPoint(0, 0));
                ImPlot::FitPoint(ImPlotPoint(h.cols, h.rows));
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

            if (plot_input_hovered() && shown) {
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                const long c = (long)floor(mp.x), r = (long)floor((double)h.rows - mp.y);
                if (c >= 0 && r >= 0 && (uint32_t)c < h.cols && (uint32_t)r < h.rows) {
                    hovered_map = idx;
                    hovered_row = (int)r;
                    hovered_col = (int)c;
                    const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint((double)c, (double)h.rows - r));
                    const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint((double)c + 1, (double)h.rows - r - 1));
                    ImPlot::GetPlotDrawList()->AddRect(p0, p1, IM_COL32(255, 255, 255, 220), 0.0f, 0, 1.5f);
                    draw_cell_tooltip(h, *v, (uint32_t)r, (uint32_t)c);
                }
            }
            ImPlot::EndPlot();
        }
        ImGui::SameLine();
        const float lo = v ? v->lo : 0.0f, hi = v ? v->hi : 1.0f;
        ImPlot::ColormapScale("##scale", lo, hi, ImVec2(scale_w, size.y), v && (v->truth || v->fraction) ? "%.2f" : "%.3g");
        ImPlot::PopColormap();
        start_tiles(h);
    }

    // The tiles in view, at the level that makes a cell a pixel or a little more. Those not of the view shown are
    // wanted, and drawn as they were, from a coarser tile, or blank; with a spinner once they have waited a while.
    void draw_tiles(Heatmap& h) {
        const View* v = h.view.get();
        const uint64_t gen = v ? v->gen : 0;
        const ImPlotRect lim = ImPlot::GetPlotLimits();
        const ImVec2 size = ImPlot::GetPlotSize();
        const ImVec2 fb = ImGui::GetIO().DisplayFramebufferScale;
        const double per_px = MD_MAX(lim.X.Size() / MD_MAX(size.x * fb.x, 1.0f), lim.Y.Size() / MD_MAX(size.y * fb.y, 1.0f));
        const uint32_t top = heatmaps::tile_top(h.rows, h.cols);
        const uint32_t level = heatmaps::tile_level(per_px, h.rows, h.cols);
        const uint64_t span = (uint64_t)TILE << level;
        h.wanted.clear();

        const uint32_t c0 = (uint32_t)MD_CLAMP(floor(lim.X.Min), 0.0, (double)h.cols), c1 = (uint32_t)MD_CLAMP(ceil(lim.X.Max), 0.0, (double)h.cols);
        const uint32_t r0 = (uint32_t)MD_CLAMP(floor((double)h.rows - lim.Y.Max), 0.0, (double)h.rows), r1 = (uint32_t)MD_CLAMP(ceil((double)h.rows - lim.Y.Min), 0.0, (double)h.rows);
        if (c1 <= c0 || r1 <= r0) return;
        const uint32_t tx0 = (uint32_t)(c0 / span), tx1 = (uint32_t)((c1 - 1) / span);
        const uint32_t ty0 = (uint32_t)(r0 / span), ty1 = (uint32_t)((r1 - 1) / span);
        const double mx = 0.5 * (lim.X.Min + lim.X.Max), my = (double)h.rows - 0.5 * (lim.Y.Min + lim.Y.Max);

        const int now = ImGui::GetFrameCount();
        const double t_now = ImGui::GetTime();
        std::vector<std::pair<double, uint64_t>> want;
        ImPlot::PushPlotClipRect();
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        const ImVec2 clip0 = ImPlot::GetPlotPos(), clip1(clip0.x + size.x, clip0.y + size.y);
        for (uint32_t ty = ty0; ty <= ty1; ++ty) {
            for (uint32_t tx = tx0; tx <= tx1; ++tx) {
                const uint64_t key = tile_key(level, tx, ty);
                Tile& t = h.tiles[key];
                t.used = now;
                const bool current = t.tex && v && t.gen == gen;
                if (current) {
                    t.waiting = 0;
                } else {
                    if (t.waiting == 0) t.waiting = t_now;
                    if (v) {
                        const double dx = ((double)tx + 0.5) * span - mx, dy = ((double)ty + 0.5) * span - my;
                        want.emplace_back(dx * dx + dy * dy, key);
                    }
                }
                if (t.tex) {
                    draw_tile(dl, h, t, level, tx, ty, level, tx, ty);
                } else {
                    bool drawn = false;
                    for (uint32_t l = level + 1; l <= top && !drawn; ++l) {
                        const uint32_t s = l - level;
                        auto p = h.tiles.find(tile_key(l, tx >> s, ty >> s));
                        if (p == h.tiles.end() || !p->second.tex) continue;
                        p->second.used = now;
                        draw_tile(dl, h, p->second, l, tx >> s, ty >> s, level, tx, ty);
                        drawn = true;
                    }
                    // Blank, where the matrix is
                    if (!drawn) {
                        const ImVec2 a = ImPlot::PlotToPixels(ImPlotPoint((double)tx * span, (double)h.rows - (double)ty * span));
                        const ImVec2 b = ImPlot::PlotToPixels(ImPlotPoint(MD_MIN((double)(tx + 1) * span, (double)h.cols), MD_MAX((double)h.rows - (double)(ty + 1) * span, 0.0)));
                        dl->AddRectFilled(a, b, IM_COL32(48, 48, 48, 255));
                    }
                }
                // Waited a while: a spinner in the middle of what is seen of it
                if (!current && t_now - t.waiting > SPINNER_AFTER) {
                    const ImVec2 a = ImPlot::PlotToPixels(ImPlotPoint((double)tx * span, (double)h.rows - (double)ty * span));
                    const ImVec2 b = ImPlot::PlotToPixels(ImPlotPoint(MD_MIN((double)(tx + 1) * span, (double)h.cols), MD_MAX((double)h.rows - (double)(ty + 1) * span, 0.0)));
                    const ImVec2 q0(MD_MAX(a.x, clip0.x), MD_MAX(a.y, clip0.y)), q1(MD_MIN(b.x, clip1.x), MD_MIN(b.y, clip1.y));
                    const float w = q1.x - q0.x, ht = q1.y - q0.y;
                    if (w >= 24.0f && ht >= 24.0f) spinner(dl, ImVec2(0.5f * (q0.x + q1.x), 0.5f * (q0.y + q1.y)), MD_MIN(10.0f, 0.2f * MD_MIN(w, ht)), t_now);
                }
            }
        }
        ImPlot::PopPlotClipRect();
        std::sort(want.begin(), want.end());
        for (const auto& w : want) h.wanted.push_back(w.second);

        // The least recently drawn go, none drawn now
        const size_t keep = MD_MAX((size_t)TILE_KEEP, 2 * (size_t)(tx1 - tx0 + 1) * (ty1 - ty0 + 1));
        if (h.tiles.size() > keep) {
            std::vector<std::pair<int, uint64_t>> old;
            for (const auto& kv : h.tiles) if (kv.second.used != now) old.emplace_back(kv.second.used, kv.first);
            std::sort(old.begin(), old.end());
            for (size_t i = 0; i < old.size() && h.tiles.size() > keep; ++i) {
                auto it = h.tiles.find(old[i].second);
                gl::free_texture(&it->second.tex);
                h.tiles.erase(it);
            }
        }
    }

    // Where tile (level, tx, ty) is, drawn from t, tile (tl, ttx, tty), the same or one coarser holding it
    static void draw_tile(ImDrawList* dl, const Heatmap& h, const Tile& t, uint32_t tl, uint32_t ttx, uint32_t tty, uint32_t level, uint32_t tx, uint32_t ty) {
        const double s = (double)((uint64_t)TILE << level), ts = (double)((uint64_t)TILE << tl);
        const double x = tx * s, y = ty * s;        // Its first column and row
        const double u = x - ttx * ts, w = y - tty * ts;
        const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint(x, (double)h.rows - y));
        const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint(x + s, (double)h.rows - y - s));
        dl->AddImage((ImTextureID)(intptr_t)t.tex, p0, p1, ImVec2((float)(u / ts), (float)(w / ts)), ImVec2((float)((u + s) / ts), (float)((w + s) / ts)));
    }

    void draw_cell_tooltip(const Heatmap& h, const View& v, uint32_t r, uint32_t c) {
        const char* rs = r < h.row_labels.size() ? h.row_labels[r].c_str() : "";
        const char* cs = c < h.col_labels.size() ? h.col_labels[c].c_str() : "";
        char text[128];
        cell_text(text, sizeof(text), v, v.value(r, c));
        ImGui::BeginTooltip();
        if (rs[0]) ImGui::TextUnformatted(rs); else ImGui::Text("Row %u", r);
        ImGui::SameLine();
        ImGui::TextDisabled("and");
        ImGui::SameLine();
        if (cs[0]) ImGui::TextUnformatted(cs); else ImGui::Text("column %u", c);
        ImGui::TextUnformatted(text);
        // A fraction: when, as a strip over the frames, the run's start at the left
        if (v.fraction && v.pf) {
            md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
            const md_bitset_t when = heatmaps::sparse_frames(v.sparse(), r, c, temp);
            const uint32_t nf = md_pair_frames_num_frames(v.pf);
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
                        Heatmap& h = *maps[i];
                        bool open = true;
                        const ImGuiTabItemFlags tf = h.select ? ImGuiTabItemFlags_SetSelected : 0;
                        h.select = false;
                        if (ImGui::BeginTabItem(h.name, &open, tf)) {
                            request(h);
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
            const Heatmap& h = *maps[hovered_map];
            const md_str_t name = md_str_from_cstr(h.name);
            md_bitmask_t* mask = &app->system_context.selection.highlight_mask;
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
            md_bitmask_clear(&app->system_context.selection.highlight_mask);
            highlighting = false;
        }
    }
};

static Heatmaps instance;

}  // namespace
