// SHAPESPACE
//
// The shapes of structures (sets of atoms) as the eigenvalues of their gyration tensors put them: how linear, planar or
// spherical each is, a point in the triangle with those three at its corners.
//
// What is plotted is a stack of LAYERS the user defines, as in the Ramachandran window: the structures a filter gives
// (one for each set of atoms it picks: residue(:) is a structure per residue), weighted by mass or not, over the current
// frame, all frames, or the frames of a timeline period, drawn as points (in the layer's colour, by structure, or by
// time), as contours of their density (at levels of the fraction of the samples they enclose) or as a heat map of it.
// The first layer is drawn on top.
//
// The shapes over every frame are an EVALUATION, made once for the structures and the weighting, off the main thread,
// and shared by every layer of the same: frames are read once whatever the number of layers over them. The points of
// the frames done are shown as it goes. The shapes of the frame shown are found from its state, as it changes.
//
// The densities are density_map.h's, over a bounded plane clipped to the triangle.
//
//   hover                  highlights the structure of the point under the pointer, and says what it is
//   click                  on a point of a frame, goes to that frame
//   shift-drag             selects the structures of the points in the box (right button: deselects them)
//   shift-click            on a point, selects its structure (right button: deselects it); shift-right click elsewhere clears
//   double-click           fits the view
//   a period dropped       from the timeline on the plot, a new layer of the density over it; on a layer, its frames

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_os.h>

#include <md_csv.h>
#include <md_xvg.h>
#include <md_util.h>
#include <md_filter.h>

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>

#include <serialization_utils.h>
#include <timeline_periods_ui.h>
#include <filter_input.h>
#include <property_grid.h>
#include <density_map.h>
#include <layer_groups.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <implot_internal.h>
#include <task_system.h>

#include <atomic>
#include <memory>
#include <vector>
#include <math.h>
#include <string.h>

namespace shapespace {

// ## Constants

// The corners of the triangle, in the plot's coordinates: an equilateral triangle of side 1
static const vec2_t tri_lin = { 0.0f, 0.0f };
static const vec2_t tri_pla = { 1.0f, 0.0f };
static const vec2_t tri_iso = { 0.5f, 0.86602540378f };

// The plane the densities are over: square, around the triangle, with room for the images of the samples near its
// edges that keep the density from leaking out of it (density_map.h), some 4 sigma
static const DensityMapDomain shape_domain = {
    -0.25f, 1.25f, 0.433f - 0.75f, 0.433f + 0.75f, false,     // 0.25 = 4 * MAX_BANDWIDTH, past the corners
    true, { tri_lin.x, tri_lin.y, tri_pla.x, tri_pla.y, tri_iso.x, tri_iso.y },
};

constexpr int MAX_LAYERS      = 16;
constexpr int MAX_LEVELS      = DENSITY_MAP_MAX_LEVELS;
constexpr size_t MAX_POINTS   = 100000;     // Of a layer of points over frames: frames are skipped beyond
constexpr float MAX_BANDWIDTH = 0.06f;      // The domain has room around the triangle for images within 4 sigma of it

static const ImVec4 layer_palette[] = {
    {0.95f, 0.95f, 0.95f, 1.0f},
    {1.00f, 0.62f, 0.20f, 1.0f},
    {0.35f, 0.80f, 1.00f, 1.0f},
    {0.95f, 0.42f, 0.80f, 1.0f},
    {0.55f, 0.90f, 0.35f, 1.0f},
    {1.00f, 0.88f, 0.30f, 1.0f},
    {0.68f, 0.58f, 1.00f, 1.0f},
    {1.00f, 0.42f, 0.42f, 1.0f},
};

static inline vec2_t shape_coord(vec3_t w) {
    return tri_lin * w.x + tri_pla * w.y + tri_iso * w.z;
}

// The weights back from the point: they sum to one (md_util_shape_weights), which the corners' barycentric coordinates
// are, so the point is all there is to keep
static inline vec3_t shape_weights_of(vec2_t c) {
    const float iso = c.y / tri_iso.y;
    const float pla = c.x - tri_iso.x * iso;
    return { 1.0f - pla - iso, pla, iso };
}

static inline bool coord_valid(vec2_t c) {
    return c.x == c.x && c.y == c.y;    // Not NaN: the frames not evaluated yet are
}

// The shape weights (linear, planar, isotropic) of the atoms of a structure in a state: the eigenvalues of the gyration
// tensor about their centre, placed as one where the cell is periodic. masses may be NULL. NaN where there is no shape
// (fewer than two atoms, or all in one place).
static vec3_t structure_shape(md_bitset_t atoms, const vec3_t* xyz, const md_unitcell_t* cell, const float* masses, md_temp_scope_t temp) {
    const size_t count = md_bitset_count(atoms);
    if (count < 2) return { NAN, NAN, NAN };   // No shape: not a point
    vec4_t* xyzw = (vec4_t*)md_temp_alloc_aligned(temp, sizeof(vec4_t) * count, 16);
    md_bitset_iter_t it = md_bitset_iter(atoms);
    size_t k = 0;
    while (md_bitset_iter_next(&it) && k < count) {
        xyzw[k++] = vec4_from_vec3(xyz[it.idx], masses ? masses[it.idx] : 1.0f);
    }
    // Covariance about the mean of the placed points, not the circular mean the placement is seeded from (their offset
    // enters the covariance squared)
    vec3_t com = vec3_zero();
    md_util_deperiodize_self_vec4(xyzw, k, cell, &com);
    const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, k, com);
    return md_util_shape_weights(&M);
}

// ## Evaluations

// The shapes of a set of structures over every frame of the run, weighted by mass or not: shared by the layers of the
// same structures and weighting, and freed when none uses it
struct Evaluation {
    uint64_t key = 0;                       // Of the structures and the weighting
    md_allocator_i* arena = nullptr;
    md_array(md_bitset_t) structures = 0;   // In arena
    bool     use_mass = true;
    size_t   num_frames = 0;
    size_t   num_structures = 0;
    vec2_t*  coords = nullptr;              // [structure * num_frames + frame], in the triangle: NaN for frames not done
    std::atomic<uint32_t> frames_done{0};
    std::atomic<int32_t>  cancel{0};
    task_system::ID task = task_system::INVALID_ID;

    bool complete() const { return num_frames > 0 && frames_done.load() >= num_frames; }
    // Stopped before all of its frames were done (interrupted from the task overlay, or by a load): never will be
    // (Asked in this order: a task that has stopped has counted all the frames it did.)
    bool failed() const { return !task_system::task_is_running(task) && !complete(); }

    // Waits for its task, which on the main thread may run pool work meanwhile: what is dropped in use goes through
    // Shapespace::retired first, freed once its task has stopped
    ~Evaluation() {
        cancel.store(1);
        task_system::task_interrupt_and_wait_for(task);
        if (arena) md_arena_allocator_destroy(arena);
    }
};

using EvaluationRef = std::shared_ptr<Evaluation>;

// ## Layers

enum LayerFrames : int {
    LayerFrames_Current = 0,
    LayerFrames_All,
    LayerFrames_Period,
    LayerFrames_Count,
};

enum LayerStyle : int {
    LayerStyle_Points = 0,
    LayerStyle_Contours,
    LayerStyle_Heatmap,
    LayerStyle_Count,
};

enum PointColoring : int {
    PointColoring_Layer = 0,    // The layer's colour
    PointColoring_Structure,    // A colour for each structure
    PointColoring_Time,         // The frame, through the layer's colour map
    PointColoring_Count,
};

static const char* style_names[LayerStyle_Count] = { "Points", "Contours", "Heatmap" };
static const char* point_coloring_names[PointColoring_Count] = { "Layer color", "Structure", "Time" };

struct Layer {
    uint32_t id = 0;
    RepresentationNode node = {};   // Where it is in the list: the group it is in (layer_groups.h). node.id is id.
    char     name[64] = "";
    bool     visible = true;

    AtomFilter  filter = {};        // The structures: a set of atoms each, within its group's subset. Empty: the
                                    // group's structures, or without a group all atoms, as one.
    bool        use_mass = true;    // Weighted by mass
    LayerFrames frames = LayerFrames_All;
    uint32_t    period = 0;

    LayerStyle style = LayerStyle_Points;
    ImVec4   color = {1, 1, 1, 1};
    float    opacity = 1.0f;

    bool     fill = false;
    bool     lines = true;
    float    line_width = 1.5f;
    int      num_levels = 3;
    float    levels[MAX_LEVELS] = { 0.5f, 0.8f, 0.95f };
    float    bandwidth = 0.0f;      // Of the gaussian the density is smoothed with, of the triangle's side. 0: automatic.

    ImPlotColormap colormap = ImPlotColormap_Viridis;
    bool     log_scale = false;

    PointColoring point_coloring = PointColoring_Layer;
    float    point_size = 4.0f;

    // ### Derived, not saved

    md_array(md_bitset_t) structures = 0;   // Heap allocated, as the filter gave them
    uint64_t structs_key = 0;       // Of what they were found from
    uint64_t structs_hash = 0;      // Of the structures themselves
    bool     structs_valid = false;

    EvaluationRef eval;             // Of frames: over the run

    // Of the frame shown
    md_array(vec2_t) current_coords = 0;
    uint64_t current_key = 0;

    // The density
    uint32_t tex = 0;
    DensityMapResult result = {};
    bool     has_result = false;
    uint64_t shown_hash = 0;
    uint64_t want_hash = 0;
    DensityMapJob* job = nullptr;
    uint64_t job_serial = 0;        // Of the job launched last: the callbacks of cancelled ones are not its

    // The bandwidth found for the samples (Scott's rule), and what of
    float    auto_bandwidth = 0.0f;
    uint64_t auto_key = 0;

    uint32_t point_stride = 1;
    uint64_t point_count = 0;
};

// Scott's rule for two dimensions, sigma n^(-1/6), of the spread of the samples (both axes alike, as the plane is
// isotropic): what smooths a single gaussian best, and much of anything else well enough. From at most some 100k of them.
static float scott_bandwidth(const Evaluation& e, uint32_t beg, uint32_t end) {
    const uint64_t total = (uint64_t)e.num_structures * (end - beg);
    const uint32_t step = (uint32_t)MAX((uint64_t)1, total / 100000);
    double n = 0, sx = 0, sy = 0, sxx = 0, syy = 0;
    for (size_t s = 0; s < e.num_structures; ++s) {
        const vec2_t* c = e.coords + s * e.num_frames;
        for (uint32_t f = beg; f < end; f += step) {
            if (!coord_valid(c[f])) continue;
            n += 1;
            sx += c[f].x;
            sy += c[f].y;
            sxx += (double)c[f].x * c[f].x;
            syy += (double)c[f].y * c[f].y;
        }
    }
    if (n < 2) return 0.02f;
    const double var = MAX(0.0, (sxx - sx * sx / n) + (syy - sy * sy / n)) / (2.0 * (n - 1));
    // n: of all the samples, not of those the spread was taken from
    const double h = sqrt(var) * pow((double)total, -1.0 / 6.0);
    return CLAMP((float)h, 0.005f, MAX_BANDWIDTH);
}

struct Shapespace : viamd::EventHandler {
    ApplicationState* app = nullptr;
    bool show_window = false;

    // In the order of the list (apply_tree): what is listed first is drawn on top
    Layer layers[MAX_LAYERS];
    int   num_layers = 0;
    uint32_t next_layer_id = 1;
    uint32_t selected_layer = 0;    // By id: the one the inspector shows, or
    uint32_t selected_group = 0;    // the group it shows

    LayerGroup groups[LAYER_MAX_GROUPS];
    int   num_groups = 0;
    uint32_t next_group_id = 1;
    bool  tree_dirty = false;       // The layers are to be put in the order of the list again
    static constexpr uint32_t LIST_ID = 2;  // What its rows are dragged as (LayerNodeRef)
    // The groups read from a workspace, by the ids the file has for them
    int      num_read_groups = 0;
    int      read_group_file_id[LAYER_MAX_GROUPS] = {};
    uint32_t read_group_id[LAYER_MAX_GROUPS] = {};

    std::vector<EvaluationRef> evaluations;
    std::vector<EvaluationRef> retired;     // Cancelled, freed once their tasks have stopped
    uint64_t system_generation = 1;
    uint64_t state_generation = 1;          // Of the state shown: its coordinates may change with nothing else

    bool show_sidebar = true;
    bool show_grid = true;

    // The highlight is whoever sets it to clear: what this window highlighted and no longer does is cleared once
    bool owns_highlight = false;
    bool highlights_now = false;
    bool plot_drawn = false;                // This frame: a box selection is of a plot that is drawn

    // The view: fitted to the triangle (and its labels) at the plot's shape, until it is panned or zoomed
    bool   fit_view = true;
    bool   view_fitted = false;
    ImVec2 plot_size = {};
    ImPlotRect fitted_rect;
    DensityMapCanvas canvas = {};

    // Box selection
    struct {
        bool active = false;
        SelectionOperator op = SelectionOperator::Or;
        ImPlotPoint p0, p1;
        int64_t click_structure = -1;   // Of the layer click_layer: the point pressed on
        int click_layer = -1;
    } box;

    // Of a workspace from before there were layers
    struct {
        bool  present = false;
        char  filter[256] = "all";
        float marker_size = 1.5f;
        bool  use_mass = true;
        int   period = 0;
    } legacy;
    bool read_layers = false;
    bool layers_declared = false;

    Shapespace() { viamd::event_system_register_handler(*this); }

    // ## Events

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app = (ApplicationState*)e.payload;
                workspace_register_window("ShapeSpace", &show_window);
                create_default_layers();
                break;
            case viamd::EventType_ViamdShutdown:
                release_all(true);
                clear_layers();
                density_map_canvas_free(&canvas);
                density_map_shutdown();
                break;
            case viamd::EventType_ViamdFrameTick:
                update();
                draw();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Shapespace", &show_window);
                break;
            case viamd::EventType_ViamdSystemInit:
                system_generation += 1;
                break;
            case viamd::EventType_ViamdSystemFree:
            case viamd::EventType_ViamdTrajectoryFree:
                // What the evaluations read goes with it
                release_all(true);
                system_generation += 1;
                break;
            case viamd::EventType_ViamdLoadData: {
                // A system or a run loaded in place of the one there: a run of another file may have as many frames
                // as the one before, and is opened without a SystemInit or SystemFree (the tasks reading the old one
                // were interrupted and waited for already)
                const LoadDataPayload* load = (const LoadDataPayload*)e.payload;
                if (load && (load->loader_state.flags & (LoaderFlag_System | LoaderFlag_Trajectory))) {
                    release_all(true);
                    system_generation += 1;
                }
                break;
            }
            case viamd::EventType_ViamdSystemStateChanged:
                // Coordinates written in place (a QM optimization step, a vibration), the frame the same
                state_generation += 1;
                break;
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                clear_layers();
                legacy = {};
                read_layers = false;
                layers_declared = false;
                num_read_groups = 0;
                show_sidebar = true;
                show_grid = true;
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                const str_t section = viamd::section_header(state);
                if (str_eq(section, STR_LIT("ShapeSpace"))) {
                    deserialize_settings(state);
                } else if (str_eq(section, STR_LIT("ShapeSpaceLayer"))) {
                    deserialize_layer(state);
                } else if (str_eq(section, STR_LIT("ShapeSpaceGroup"))) {
                    deserialize_group(state);
                }
                break;
            }
            case viamd::EventType_ViamdDeserializeEnd:
                if (!read_layers && !layers_declared) {
                    // Over the timeline filter of a workspace from before there were periods, as it was
                    if (legacy.present || (app && app->timeline.periods.legacy_filter_id != 0)) create_legacy_layers();
                    else create_default_layers();
                }
                break;
            default:
                break;
            }
        }
    }

    // Every evaluation and density let go: what they were of is gone
    void release_all(bool wait) {
        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            cancel_job(l, wait);
            l.eval.reset();
            l.has_result = false;
            l.shown_hash = 0;
            l.want_hash = 0;
            l.structs_key = 0;
            l.current_key = 0;
        }
        for (EvaluationRef& e : evaluations) e->cancel.store(1);
        for (EvaluationRef& e : retired) e->cancel.store(1);
        for (EvaluationRef& e : evaluations) task_system::task_interrupt_and_wait_for(e->task);
        for (EvaluationRef& e : retired) task_system::task_interrupt_and_wait_for(e->task);
        evaluations.clear();
        retired.clear();
        density_map_canvas_invalidate(&canvas);
    }

    // ## Layers

    Layer* find_layer(uint32_t id) {
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].id == id) return &layers[i];
        }
        return nullptr;
    }

    int layer_index(uint32_t id) const {
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].id == id) return i;
        }
        return -1;
    }

    ImVec4 next_layer_color() const {
        int best = 0, best_uses = INT32_MAX;
        for (int p = 0; p < (int)ARRAY_SIZE(layer_palette); ++p) {
            int uses = 0;
            for (int i = 0; i < num_layers; ++i) {
                const ImVec4 c = layers[i].color;
                if (c.x == layer_palette[p].x && c.y == layer_palette[p].y && c.z == layer_palette[p].z) ++uses;
            }
            if (uses < best_uses) {
                best = p;
                best_uses = uses;
            }
        }
        return layer_palette[best];
    }

    const Layer* find_layer(uint32_t id) const {
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].id == id) return &layers[i];
        }
        return nullptr;
    }

    // ### The list: layers and groups

    size_t tree_nodes(RepresentationNode* out[MAX_LAYERS + LAYER_MAX_GROUPS]) {
        size_t n = 0;
        for (int i = 0; i < num_layers; ++i) out[n++] = &layers[i].node;
        for (int i = 0; i < num_groups; ++i) out[n++] = &groups[i].node;
        return n;
    }

    // The tree made canonical, and the layers put in its order (moved: what points at them is no longer valid)
    void apply_tree() {
        RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
        representation_tree_normalize(nodes, tree_nodes(nodes));
        // Depth first, the orders are the positions in the list: a handful, sorted in place
        bool moved = false;
        for (int i = 1; i < num_layers; ++i) {
            for (int j = i; j > 0 && layers[j - 1].node.order > layers[j].node.order; --j) {
                Layer tmp;
                move_into(tmp, layers[j]);
                move_into(layers[j], layers[j - 1]);
                move_into(layers[j - 1], tmp);
                moved = true;
            }
        }
        if (moved) density_map_canvas_invalidate(&canvas);
        tree_dirty = false;
    }

    // The tree canonical again after a node was added or changed, the layers put in order later (apply_tree)
    void tree_changed() {
        RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
        representation_tree_normalize(nodes, tree_nodes(nodes));
        tree_dirty = true;
    }

    LayerGroup* find_group(uint32_t id) {
        for (int i = 0; i < num_groups; ++i) {
            if (groups[i].node.id == id) return &groups[i];
        }
        return nullptr;
    }

    const LayerGroup* find_group(uint32_t id) const {
        for (int i = 0; i < num_groups; ++i) {
            if (groups[i].node.id == id) return &groups[i];
        }
        return nullptr;
    }

    // The group a layer is in, NULL at the top level
    const LayerGroup* layer_group(const Layer& l) const { return l.node.parent ? find_group(l.node.parent) : nullptr; }

    // Its own eye, and its group's
    bool layer_shown(const Layer& l) const {
        if (!l.visible) return false;
        const LayerGroup* g = layer_group(l);
        return !g || g->visible;
    }

    // The group of what is selected, where what is added goes
    uint32_t current_group() const {
        if (selected_group) return selected_group;
        if (const Layer* l = find_layer(selected_layer)) return l->node.parent;
        return 0;
    }

    // Adds a layer: first in the group 'parent' (or at the top) for 'at' 0, last for -1. NULL when there is no room.
    // In place in the list once the tree is applied (apply_tree).
    Layer* add_layer(LayerFrames frames, LayerStyle style, int at = -1, uint32_t parent = 0) {
        if (num_layers >= MAX_LAYERS) return nullptr;
        Layer& l = layers[num_layers++];
        l.eval.reset();
        new (&l) Layer();
        l.id = next_layer_id++;
        l.node.id = l.id;
        l.node.parent = parent;
        if (at < 0) {
            l.node.order = LAYER_TREE_LAST;
        } else {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            l.node.order = layer_tree_first_order(nodes, tree_nodes(nodes), parent);
        }
        l.frames = frames;
        l.style = style;
        l.color = next_layer_color();
        tree_changed();
        return &l;
    }

    void free_layer(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        cancel_job(l, true);
        density_map_texture_free(&l.tex);
        for (size_t i = 0; i < md_array_size(l.structures); ++i) md_bitset_free(&l.structures[i], heap);
        md_array_free(l.structures, heap);
        md_array_free(l.current_coords, heap);
        l.eval.reset();
    }

    // Layers are moved as bytes but for their evaluation, which is counted: moved from one to the other
    static void move_into(Layer& dst, Layer& src) {
        EvaluationRef e = std::move(src.eval);
        dst.eval.reset();
        MEMCPY((void*)&dst, (const void*)&src, sizeof(Layer));
        new (&dst.eval) EvaluationRef(std::move(e));
        new (&src.eval) EvaluationRef();
    }

    void remove_layer(uint32_t id) {
        const int idx = layer_index(id);
        if (idx < 0) return;
        free_layer(layers[idx]);
        for (int i = idx; i + 1 < num_layers; ++i) move_into(layers[i], layers[i + 1]);
        num_layers -= 1;
        layers[num_layers].eval.reset();
        new (&layers[num_layers]) Layer();
        if (selected_layer == id) selected_layer = num_layers > 0 ? layers[MIN(idx, num_layers - 1)].id : 0;
        tree_changed();
        density_map_canvas_invalidate(&canvas);
    }

    // A copy of what is set, not what is derived (found again), at 'order' in the group 'parent'
    Layer* copy_layer(uint32_t src_id, uint32_t parent, int32_t order) {
        const int idx = layer_index(src_id);
        if (idx < 0 || num_layers >= MAX_LAYERS) return nullptr;
        Layer* l = add_layer(layers[idx].frames, layers[idx].style, -1, parent);
        if (!l) return nullptr;
        const Layer& src = layers[idx];
        STATIC_ASSERT(sizeof(l->name) == sizeof(src.name), "");
        MEMCPY(l->name, src.name, sizeof(l->name));
        l->visible = src.visible;
        l->filter = src.filter;
        l->use_mass = src.use_mass;
        l->frames = src.frames;
        l->period = src.period;
        l->style = src.style;
        l->color = src.color;
        l->opacity = src.opacity;
        l->fill = src.fill;
        l->lines = src.lines;
        l->line_width = src.line_width;
        l->num_levels = src.num_levels;
        MEMCPY(l->levels, src.levels, sizeof(l->levels));
        l->bandwidth = src.bandwidth;
        l->colormap = src.colormap;
        l->log_scale = src.log_scale;
        l->point_coloring = src.point_coloring;
        l->point_size = src.point_size;
        l->node.order = order;
        tree_changed();
        return l;
    }

    // Just above it, in its group
    void duplicate_layer(uint32_t id) {
        const Layer* src = find_layer(id);
        if (!src) return;
        if (Layer* l = copy_layer(id, src->node.parent, src->node.order - 1)) {
            selected_layer = l->id;
            selected_group = 0;
        }
    }

    void clear_layers() {
        for (int i = 0; i < num_layers; ++i) {
            free_layer(layers[i]);
            new (&layers[i]) Layer();
        }
        num_layers = 0;
        for (int i = 0; i < num_groups; ++i) {
            layer_group_free(&groups[i]);
            groups[i] = LayerGroup{};
        }
        num_groups = 0;
        selected_layer = 0;
        selected_group = 0;
        tree_dirty = false;
        density_map_canvas_invalidate(&canvas);
    }

    // ### Groups

    // A group at 'order' of the top level, NULL when there is no room
    LayerGroup* add_group(int32_t order) {
        if (num_groups >= LAYER_MAX_GROUPS) return nullptr;
        LayerGroup& g = groups[num_groups++];
        g = LayerGroup{};
        g.node = {next_group_id++, 0, order, true};
        tree_changed();
        return &g;
    }

    // A group of the layer, where it was, with the layer's structures as its subset: the layers added to it show them
    void group_layer(uint32_t id) {
        Layer* l = find_layer(id);
        if (!l || l->node.parent != 0) return;
        LayerGroup* g = add_group(l->node.order - 1);
        if (!g) return;
        MEMCPY(g->filter.text, l->filter.text, sizeof(g->filter.text));
        l->filter.text[0] = '\0';
        l->structs_key = 0;
        l->node.parent = g->node.id;
        l->node.order = LAYER_TREE_LAST;
        tree_changed();
        selected_group = g->node.id;
        selected_layer = 0;
    }

    // In a group now: a filter that is the group's own says nothing more (and within it would give its structures
    // again, or fewer), so it goes
    void entered_group(Layer& l) {
        l.structs_key = 0;
        const LayerGroup* g = layer_group(l);
        if (g && !g->follow_selection && g->filter.text[0] && strcmp(g->filter.text, l.filter.text) == 0) l.filter.text[0] = '\0';
    }

    // Into the group (0: out of any), last in it, with its filter as it is: what it shows is then within the group's
    void move_layer_to_group(uint32_t id, uint32_t group_id) {
        Layer* l = find_layer(id);
        if (!l || l->node.parent == group_id) return;
        l->node.parent = group_id;
        l->node.order = LAYER_TREE_LAST;
        entered_group(*l);
        tree_changed();
    }

    void move_node(LayerNodeRef src, LayerNodeRef target, RepresentationDrop drop) {
        RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t n = tree_nodes(nodes);
        const int si = representation_tree_find(nodes, n, src.group, src.id);
        const int ti = target.id ? representation_tree_find(nodes, n, target.group, target.id) : -1;    // None: the end of the list
        if (si < 0 || (target.id && ti < 0)) return;
        if (layer_tree_move(nodes, n, (size_t)si, ti, drop)) {
            if (!src.group) {
                if (Layer* l = find_layer(src.id)) entered_group(*l);
            }
            tree_dirty = true;
        }
    }

    // With what is in it: removed, or (ungroup) taken out to where the group was, each with a filter that keeps what it
    // shows. A copy: just below it, with copies of what is in it.
    void group_action(uint32_t id, LayerGroupAction action) {
        LayerGroup* g = find_group(id);
        if (!g || action == LayerGroupAction_None) return;
        if (action == LayerGroupAction_Duplicate) {
            char name[sizeof(g->name)];
            char filter[sizeof(g->filter.text)];
            MEMCPY(name, g->name, sizeof(name));
            MEMCPY(filter, g->filter.text, sizeof(filter));
            const bool visible = g->visible, open = g->open, follow = g->follow_selection;
            LayerGroup* copy = add_group(g->node.order + 1);
            if (!copy) return;
            MEMCPY(copy->name, name, sizeof(name));
            MEMCPY(copy->filter.text, filter, sizeof(filter));
            copy->visible = visible;
            copy->open = open;
            copy->follow_selection = follow;
            const uint32_t copy_id = copy->node.id;
            // In the order they are in, each last in the copy
            uint32_t ids[MAX_LAYERS];
            int n = 0;
            {
                RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
                const size_t nn = tree_nodes(nodes);
                int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
                const size_t nk = representation_tree_children(kids, ARRAY_SIZE(kids), nodes, nn, id);
                for (size_t k = 0; k < nk; ++k) {
                    if (kids[k] < num_layers) ids[n++] = layers[kids[k]].id;
                }
            }
            for (int k = 0; k < n; ++k) copy_layer(ids[k], copy_id, LAYER_TREE_LAST);
            selected_group = copy_id;
            selected_layer = 0;
            return;
        }
        if (action == LayerGroupAction_Remove) {
            for (int i = num_layers - 1; i >= 0; --i) {
                if (layers[i].node.parent == id) remove_layer(layers[i].id);
            }
            g = find_group(id);
        } else if (action == LayerGroupAction_Ungroup) {
            const int32_t at = g->node.order + 1;   // Where the group was, in the order they were in it
            for (int i = 0; i < num_layers; ++i) {
                Layer& l = layers[i];
                if (l.node.parent != id) continue;
                // What can be kept is (layer_group_keeps_filter); what cannot is left as it is, as the menu warned
                if (layer_group_keeps_filter(l.filter.text, false, *g) && !g->follow_selection) {
                    layer_group_compose_filter(l.filter.text, sizeof(l.filter.text), *g);
                }
                l.node.parent = 0;
                l.node.order = at;
                l.structs_key = 0;
            }
        }
        if (!g) return;
        const int gi = (int)(g - groups);
        layer_group_free(g);
        for (int i = gi; i + 1 < num_groups; ++i) groups[i] = groups[i + 1];
        num_groups -= 1;
        groups[num_groups] = LayerGroup{};
        if (selected_group == id) selected_group = 0;
        tree_changed();
    }

    // In a group, as an example of one: the shape of the frame shown, its trajectory as points coloured by time, and
    // the density over it. The group's subset is where to change what they are of all at once - residue(protein()) for
    // the residues of the protein, a structure each - and what is added to it is of the same. The whole system, as one
    // structure, to begin with: what every system has, and one structure to find over the frames.
    void create_default_layers() {
        clear_layers();
        uint32_t group = 0;
        if (LayerGroup* g = add_group(LAYER_TREE_LAST)) {
            snprintf(g->name, sizeof(g->name), "System");
            snprintf(g->filter.text, sizeof(g->filter.text), "all");
            group = g->node.id;
        }
        if (Layer* cur = add_layer(LayerFrames_Current, LayerStyle_Points, -1, group)) {
            cur->color = ImVec4(1, 1, 1, 1);
            cur->point_size = 7.0f;
            selected_layer = cur->id;
        }
        if (Layer* pts = add_layer(LayerFrames_All, LayerStyle_Points, -1, group)) {
            pts->point_coloring = PointColoring_Time;
            pts->point_size = 3.0f;
            pts->opacity = 0.8f;
        }
        if (Layer* den = add_layer(LayerFrames_All, LayerStyle_Contours, -1, group)) {
            den->color = ImVec4(1, 1, 1, 1);   // Lines that stand out on the points coloured by time
        }
    }

    // ## Jobs

    void cancel_job(Layer& l, bool wait) {
        if (!l.job) return;
        density_map_cancel(l.job, wait);
        l.job = nullptr;
        l.want_hash = 0;
    }

    // The density of the shapes of the layer's evaluation over frames [beg, end)
    void launch_job(Layer& l, uint64_t hash, uint32_t beg, uint32_t end, float sigma) {
        ASSERT(!l.job && l.eval && l.eval->complete());
        EvaluationRef eval = l.eval;
        const uint64_t num_samples = (uint64_t)eval->num_structures * (end - beg);
        const uint32_t num_chunks = (uint32_t)MIN((uint64_t)MIN(8u, end - beg), (num_samples >> 20) + 1);
        auto source = [eval, beg, end](uint32_t chunk, uint32_t chunks, DensityMapBins& bins, const std::atomic<int32_t>& cancel) {
            const uint32_t n = end - beg;
            const uint32_t f0 = beg + (uint32_t)((uint64_t)n * chunk / chunks);
            const uint32_t f1 = beg + (uint32_t)((uint64_t)n * (chunk + 1) / chunks);
            for (size_t s = 0; s < eval->num_structures; ++s) {
                if (cancel.load(std::memory_order_relaxed)) return;
                const vec2_t* c = eval->coords + s * eval->num_frames;
                for (uint32_t f = f0; f < f1; ++f) {
                    if (coord_valid(c[f])) bins.add(0, c[f].x, c[f].y);
                }
            }
        };
        const uint32_t id = l.id;
        const uint64_t serial = ++l.job_serial;
        auto done = [this, id, hash, serial](const DensityMapResult* result, const float* grid) {
            Layer* layer = find_layer(id);
            if (!layer || !layer->job || layer->job_serial != serial || layer->want_hash != hash) return;
            layer->job = nullptr;
            if (!result) {
                layer->want_hash = 0;
                return;
            }
            if (!layer->tex) layer->tex = density_map_texture_create(false);
            density_map_texture_upload(layer->tex, grid);
            layer->result = *result;
            layer->has_result = true;
            layer->shown_hash = hash;
        };
        l.want_hash = hash;
        // The bandwidth is of the triangle's side, which is 1
        l.job = density_map_launch(STR_LIT("Shapespace density"), shape_domain, sigma, num_chunks, source, done);
    }

    // ## Update

    // The structures of the layer: found again when what they are found from changes
    void update_structures(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        const LayerGroup* g = layer_group(l);
        const bool in_context = g && layer_group_has_context(*g);
        uint64_t key = md_hash64_combine(system_generation, md_hash64(l.filter.text, strlen(l.filter.text), app->script.ir_fingerprint));
        key = md_hash64_combine(key, app->mold.sys.topology_version);
        if (in_context) key = md_hash64_combine(key, g->hash);
        // What depends on the coordinates follows the frame shown as points of it, but is taken once for frames that
        // are many (all of them would have to be evaluated again at every frame)
        if (l.filter.text[0] && l.filter.result.is_dynamic && l.frames == LayerFrames_Current) {
            key = md_hash64_combine(key, md_hash64_combine((uint64_t)(app->animation.frame * 4096.0), state_generation));
        }
        if (key == l.structs_key) return;
        l.structs_key = key;

        for (size_t i = 0; i < md_array_size(l.structures); ++i) md_bitset_free(&l.structures[i], heap);
        md_array_shrink(l.structures, 0);
        l.structs_valid = false;

        const size_t num_atoms = app->mold.sys.atom.count;
        if (num_atoms == 0) return;
        if (in_context && !g->valid) {
            // Nothing within a group whose subset does not evaluate
        } else if (l.filter.text[0] == '\0' && in_context) {
            // The group's structures
            for (size_t i = 0; i < md_array_size(g->sets); ++i) md_array_push(l.structures, md_bitset_copy(g->sets[i], heap), heap);
            l.structs_valid = md_array_size(l.structures) > 0;
        } else if (l.filter.text[0] == '\0') {
            md_array_push(l.structures, md_bitset_range(0, (uint32_t)num_atoms), heap);
            l.structs_valid = true;
        } else {
            l.structs_valid = filter_evaluate_sets(&l.filter, &l.structures, layer_group_filter_context(app, in_context ? g : nullptr), heap) &&
                              md_array_size(l.structures) > 0;
        }
        uint64_t h = 0x5A9E;
        for (size_t i = 0; i < md_array_size(l.structures); ++i) h = md_bitset_hash64(l.structures[i], h + i);
        l.structs_hash = h;
    }

    // The evaluation of the layer's structures and weighting over the run: one shared, or a new one started
    void update_evaluation(Layer& l) {
        if (l.frames == LayerFrames_Current || !l.structs_valid) {
            l.eval.reset();
            return;
        }
        const size_t num_frames = run_num_frames(app);
        if (num_frames == 0) {
            l.eval.reset();
            return;
        }
        const uint64_t key = md_hash64_combine(md_hash64_combine(l.structs_hash, l.use_mass ? 0x11 : 0x22), md_hash64_combine(system_generation, num_frames));
        if (l.eval && l.eval->key == key && !l.eval->failed()) return;
        for (size_t i = 0; i < evaluations.size(); ++i) {
            if (evaluations[i]->key != key) continue;
            if (evaluations[i]->failed()) {
                // Started again, for every layer over it
                retire(i);
                break;
            }
            l.eval = evaluations[i];
            return;
        }
        l.eval.reset();
        EvaluationRef e = std::make_shared<Evaluation>();
        e->key = key;
        e->use_mass = l.use_mass;
        e->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
        for (size_t i = 0; i < md_array_size(l.structures); ++i) {
            md_array_push(e->structures, md_bitset_copy(l.structures[i], e->arena), e->arena);
        }
        e->num_structures = md_array_size(e->structures);
        e->num_frames = num_frames;
        const size_t n = e->num_structures * num_frames;
        e->coords = (vec2_t*)md_alloc(e->arena, sizeof(vec2_t) * n);
        for (size_t i = 0; i < n; ++i) e->coords[i] = { NAN, NAN };

        Evaluation* ev = e.get();
        ApplicationState* state = app;
        e->task = task_system::create_pool_task(STR_LIT("Shapespace"), (uint32_t)num_frames, [ev, state](uint32_t beg, uint32_t end, uint32_t) {
            const md_system_t& sys = state->mold.sys;
            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };
            vec3_t* xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * ALIGN_TO(sys.atom.count, 8));
            float* masses = nullptr;
            if (ev->use_mass) {
                masses = (float*)md_temp_alloc(temp, sizeof(float) * ALIGN_TO(sys.atom.count, 8));
                md_atom_extract_masses(masses, 0, sys.atom.count, &sys.atom);
            }
            // One context for the range, so the run's files stay open across its frames
            const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
            md_system_extract_t* ex = md_system_extract_begin(&sys, str_from_cstr(state->mold.run), paths, ARRAY_SIZE(paths), md_get_heap_allocator());
            if (!ex) {
                // Nothing to read them from: done, without shapes
                ev->frames_done.fetch_add(end - beg);
                return;
            }
            defer { md_system_extract_end(ex); };
            for (uint32_t f = beg; f < end; ++f) {
                if (ev->cancel.load(std::memory_order_relaxed)) return;
                md_system_state_t frame = {};
                frame.num_atoms = sys.atom.count;
                frame.xyz = xyz;
                if (!md_system_extract_frame(ex, f, &frame)) {
                    ev->frames_done.fetch_add(1);
                    continue;
                }
                for (size_t s = 0; s < ev->num_structures; ++s) {
                    md_temp_scope_t inner = md_temp_begin();
                    const vec3_t w = structure_shape(ev->structures[s], xyz, &frame.unitcell, masses, inner);
                    md_temp_end(inner);
                    // Degenerate (one atom, or all in one place): no shape, NaN
                    ev->coords[s * ev->num_frames + f] = shape_coord(w);
                }
                ev->frames_done.fetch_add(1);
            }
        });
        task_system::enqueue_task(e->task);
        evaluations.push_back(e);
        l.eval = e;
    }

    // The shapes of the structures in the state shown: found again when it changes
    void update_current(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        const md_system_state_t& st = app->mold.state;
        uint64_t key = md_hash64_combine(l.structs_hash, l.use_mass ? 0x11 : 0x22);
        key = md_hash64_combine(key, md_hash64(&st.frame, sizeof(st.frame), (uint64_t)(uintptr_t)st.xyz));
        key = md_hash64_combine(key, md_hash64_combine(system_generation, state_generation));
        if (key == l.current_key) return;
        l.current_key = key;
        const size_t n = l.structs_valid ? md_array_size(l.structures) : 0;
        md_array_resize(l.current_coords, n, heap);
        if (n == 0 || !st.xyz || st.num_atoms != app->mold.sys.atom.count) {
            md_array_shrink(l.current_coords, 0);
            return;
        }
        md_temp_scope_t temp = md_temp_begin();
        float* masses = nullptr;
        if (l.use_mass) {
            masses = (float*)md_temp_alloc(temp, sizeof(float) * ALIGN_TO(app->mold.sys.atom.count, 8));
            md_atom_extract_masses(masses, 0, app->mold.sys.atom.count, &app->mold.sys.atom);
        }
        for (size_t s = 0; s < n; ++s) {
            md_temp_scope_t inner = md_temp_begin();
            l.current_coords[s] = shape_coord(structure_shape(l.structures[s], st.xyz, &st.unitcell, masses, inner));
            md_temp_end(inner);
        }
        md_temp_end(temp);
    }

    // The frames [beg, end) of the layer, false when it covers none
    bool layer_frame_range(const Layer& l, uint32_t* beg, uint32_t* end) const {
        const int num_frames = (int)run_num_frames(app);
        if (num_frames <= 0) return false;
        switch (l.frames) {
        case LayerFrames_Current: {
            const int f = CLAMP((int)(app->animation.frame + 0.5), 0, num_frames - 1);
            *beg = (uint32_t)f;
            *end = (uint32_t)f + 1;
            return true;
        }
        case LayerFrames_All:
            *beg = 0;
            *end = (uint32_t)num_frames;
            return true;
        case LayerFrames_Period: {
            const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, l.period);
            int b = 0, e = 0;
            if (!p || !timeline_period_frames(&b, &e, *p, app->animation.frame, num_frames)) return false;
            *beg = (uint32_t)b;
            *end = (uint32_t)e + 1;
            return true;
        }
        default:
            return false;
        }
    }

    // A density over frames, with frames to find it over
    bool layer_draws_density(const Layer& l) const {
        return l.style != LayerStyle_Points && l.frames != LayerFrames_Current && run_num_frames(app) >= 2;
    }

    void update() {
        if (!app || !show_window) return;
        if (tree_dirty) apply_tree();

        // The groups' subsets, before what is within them. One that depends on the coordinates follows the frame only
        // when every layer in it is of the frame shown: what is over frames could not follow it.
        for (int gi = 0; gi < num_groups; ++gi) {
            LayerGroup& g = groups[gi];
            bool follow_frame = true;
            for (int i = 0; i < num_layers; ++i) {
                if (layers[i].node.parent == g.node.id && layers[i].frames != LayerFrames_Current) follow_frame = false;
            }
            const bool dynamic = !g.follow_selection && g.filter.result.is_dynamic;
            const uint64_t generation = md_hash64_combine(system_generation, follow_frame && dynamic ? state_generation : 0);
            layer_group_update(&g, app, generation, follow_frame, true);
        }

        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            update_structures(l);
            if (!layer_shown(l)) continue;
            if (l.frames == LayerFrames_Current || run_num_frames(app) < 2) {
                update_current(l);
            }
            update_evaluation(l);

            if (!layer_draws_density(l) || !l.eval || !l.eval->complete()) continue;
            uint32_t beg = 0, end = 0;
            if (!layer_frame_range(l, &beg, &end)) continue;
            uint64_t h = l.eval->key;
            const uint32_t frames[2] = { beg, end };
            h = md_hash64(frames, sizeof(frames), h);
            // Found once for the samples, not at every frame
            if (h != l.auto_key) {
                l.auto_key = h;
                l.auto_bandwidth = scott_bandwidth(*l.eval, beg, end);
            }
            const float sigma = l.bandwidth > 0.0f ? l.bandwidth : l.auto_bandwidth;
            h = md_hash64(&sigma, sizeof(sigma), h);
            h |= 1;
            if (h == l.shown_hash || h == l.want_hash || l.job) continue;
            launch_job(l, h, beg, end, sigma);
        }
        // Evaluations no layer holds any longer (nor a density job: they hold them too)
        for (size_t i = 0; i < evaluations.size();) {
            if (evaluations[i].use_count() == 1) retire(i);
            else ++i;
        }
        for (size_t i = 0; i < retired.size();) {
            if (retired[i].use_count() == 1 && !task_system::task_is_running(retired[i]->task)) {
                retired[i] = std::move(retired.back());
                retired.pop_back();
            } else {
                ++i;
            }
        }
    }

    // Out of evaluations[i]: cancelled, and freed once its task has stopped
    void retire(size_t i) {
        evaluations[i]->cancel.store(1);
        retired.push_back(std::move(evaluations[i]));
        evaluations[i] = std::move(evaluations.back());
        evaluations.pop_back();
    }

    // ## Workspace

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("ShapeSpace"));
        viamd::write_bool(state, STR_LIT("ShowLayers"), show_sidebar);
        viamd::write_bool(state, STR_LIT("ShowGrid"), show_grid);
        viamd::write_int(state, STR_LIT("Layers"), num_layers);
        // In the order of the list: a group, then what is in it
        RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t n = tree_nodes(nodes);
        int top[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t nt = representation_tree_children(top, ARRAY_SIZE(top), nodes, n, 0);
        for (size_t k = 0; k < nt; ++k) {
            if (top[k] < num_layers) {
                serialize_layer(state, layers[top[k]]);
                continue;
            }
            const LayerGroup& g = groups[top[k] - num_layers];
            layer_group_serialize(state, STR_LIT("ShapeSpaceGroup"), g);
            int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t nk = representation_tree_children(kids, ARRAY_SIZE(kids), nodes, n, g.node.id);
            for (size_t c = 0; c < nk; ++c) {
                if (kids[c] < num_layers) serialize_layer(state, layers[kids[c]]);
            }
        }
    }

    void serialize_layer(viamd::serialization_state_t& state, const Layer& l) {
        {
            viamd::write_section_header(state, STR_LIT("ShapeSpaceLayer"));
            if (l.node.parent) viamd::write_int(state, STR_LIT("Group"), (int64_t)l.node.parent);
            if (l.name[0]) viamd::write_str(state, STR_LIT("Name"), str_from_cstr(l.name));
            viamd::write_bool(state, STR_LIT("Visible"), l.visible);
            if (l.filter.text[0]) viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(l.filter.text));
            viamd::write_bool(state, STR_LIT("UseMass"), l.use_mass);
            static const str_t frames[LayerFrames_Count] = { STR_LIT("current"), STR_LIT("all"), STR_LIT("period") };
            viamd::write_str(state, STR_LIT("Frames"), frames[l.frames]);
            if (l.frames == LayerFrames_Period) viamd::write_int(state, STR_LIT("Period"), (int)l.period);
            static const str_t styles[LayerStyle_Count] = { STR_LIT("points"), STR_LIT("contours"), STR_LIT("heatmap") };
            viamd::write_str(state, STR_LIT("Style"), styles[l.style]);
            viamd::write_vec4(state, STR_LIT("Color"), vec4_set(l.color.x, l.color.y, l.color.z, l.color.w));
            viamd::write_flt(state, STR_LIT("Opacity"), l.opacity);
            viamd::write_bool(state, STR_LIT("Fill"), l.fill);
            viamd::write_bool(state, STR_LIT("Lines"), l.lines);
            viamd::write_flt(state, STR_LIT("LineWidth"), l.line_width);
            viamd::write_flt_vec(state, STR_LIT("Levels"), l.levels, (size_t)l.num_levels);
            if (l.bandwidth > 0.0f) viamd::write_flt(state, STR_LIT("Bandwidth"), l.bandwidth);
            else                    viamd::write_str(state, STR_LIT("Bandwidth"), STR_LIT("auto"));
            viamd::write_str(state, STR_LIT("Colormap"), str_from_cstr(ImPlot::GetColormapName(l.colormap)));
            viamd::write_bool(state, STR_LIT("LogScale"), l.log_scale);
            static const str_t colorings[PointColoring_Count] = { STR_LIT("layer"), STR_LIT("structure"), STR_LIT("time") };
            viamd::write_str(state, STR_LIT("PointColoring"), colorings[l.point_coloring]);
            viamd::write_flt(state, STR_LIT("PointSize"), l.point_size);
        }
    }

    void deserialize_settings(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "ShowLayers")) {
                viamd::extract_bool(show_sidebar, arg);
            } else if (str_eq_cstr(ident, "ShowGrid")) {
                viamd::extract_bool(show_grid, arg);
            } else if (str_eq_cstr(ident, "Layers")) {
                layers_declared = true;
            }
            // From before there were layers
            else if (str_eq_cstr(ident, "Filter")) {
                legacy.present = true;
                viamd::extract_to_char_buf(legacy.filter, sizeof(legacy.filter), arg);
            } else if (str_eq_cstr(ident, "MarkerSize")) {
                legacy.present = true;
                viamd::extract_flt(legacy.marker_size, arg);
            } else if (str_eq_cstr(ident, "UseMass")) {
                legacy.present = true;
                viamd::extract_bool(legacy.use_mass, arg);
            } else if (str_eq_cstr(ident, "Period")) {
                legacy.present = true;
                viamd::extract_int(legacy.period, arg);
            }
        }
    }

    // As the list is read: each last of where it goes
    void deserialize_group(viamd::deserialization_state_t& state) {
        if (!read_layers) {
            clear_layers();
            read_layers = true;
        }
        str_t ident, arg;
        LayerGroup* g = add_group(LAYER_TREE_LAST);
        if (!g) {
            while (viamd::next_entry(ident, arg, state)) {}
            return;
        }
        int file_id = 0;
        while (viamd::next_entry(ident, arg, state)) {
            layer_group_deserialize_entry(g, ident, arg, &file_id);
        }
        if (num_read_groups < LAYER_MAX_GROUPS) {
            read_group_file_id[num_read_groups] = file_id;
            read_group_id[num_read_groups] = g->node.id;
            num_read_groups += 1;
        }
    }

    void deserialize_layer(viamd::deserialization_state_t& state) {
        if (!read_layers) {
            clear_layers();
            read_layers = true;
        }
        Layer* l = add_layer(LayerFrames_All, LayerStyle_Points);
        str_t ident, arg;
        if (!l) {
            while (viamd::next_entry(ident, arg, state)) {}
            return;
        }
        while (viamd::next_entry(ident, arg, state)) {
            str_t s = {};
            if (str_eq_cstr(ident, "Name")) {
                viamd::extract_to_char_buf(l->name, sizeof(l->name), arg);
            } else if (str_eq_cstr(ident, "Group")) {
                int file_id = 0;
                if (viamd::extract_int(file_id, arg)) {
                    for (int k = 0; k < num_read_groups; ++k) {
                        if (read_group_file_id[k] == file_id) {
                            l->node.parent = read_group_id[k];
                            tree_changed();     // Last in it
                        }
                    }
                }
            } else if (str_eq_cstr(ident, "Visible")) {
                viamd::extract_bool(l->visible, arg);
            } else if (str_eq_cstr(ident, "Filter")) {
                viamd::extract_to_char_buf(l->filter.text, sizeof(l->filter.text), arg);
            } else if (str_eq_cstr(ident, "UseMass")) {
                viamd::extract_bool(l->use_mass, arg);
            } else if (str_eq_cstr(ident, "Frames") && viamd::extract_str(s, arg)) {
                l->frames = str_eq_cstr(s, "current") ? LayerFrames_Current : str_eq_cstr(s, "period") ? LayerFrames_Period : LayerFrames_All;
            } else if (str_eq_cstr(ident, "Period")) {
                int id = 0;
                if (viamd::extract_int(id, arg)) l->period = (uint32_t)MAX(id, 0);
            } else if (str_eq_cstr(ident, "Style") && viamd::extract_str(s, arg)) {
                l->style = str_eq_cstr(s, "contours") ? LayerStyle_Contours : str_eq_cstr(s, "heatmap") ? LayerStyle_Heatmap : LayerStyle_Points;
            } else if (str_eq_cstr(ident, "Color")) {
                vec4_t c;
                if (viamd::extract_vec4(c, arg)) l->color = ImVec4(c.x, c.y, c.z, c.w);
            } else if (str_eq_cstr(ident, "Opacity")) {
                viamd::extract_flt(l->opacity, arg);
            } else if (str_eq_cstr(ident, "Fill")) {
                viamd::extract_bool(l->fill, arg);
            } else if (str_eq_cstr(ident, "Lines")) {
                viamd::extract_bool(l->lines, arg);
            } else if (str_eq_cstr(ident, "LineWidth")) {
                viamd::extract_flt(l->line_width, arg);
            } else if (str_eq_cstr(ident, "Levels")) {
                float v[MAX_LEVELS];
                for (int k = MAX_LEVELS; k > 0; --k) {
                    if (viamd::extract_flt_vec(v, (size_t)k, arg)) {
                        l->num_levels = k;
                        MEMCPY(l->levels, v, sizeof(float) * k);
                        break;
                    }
                }
            } else if (str_eq_cstr(ident, "Bandwidth")) {
                if (!viamd::extract_flt(l->bandwidth, arg)) l->bandwidth = 0.0f;  // auto
            } else if (str_eq_cstr(ident, "Colormap") && viamd::extract_str(s, arg)) {
                char name[64];
                str_copy_to_char_buf(name, sizeof(name), s);
                const ImPlotColormap cm = ImPlot::GetColormapIndex(name);
                if (cm >= 0) l->colormap = cm;
            } else if (str_eq_cstr(ident, "LogScale")) {
                viamd::extract_bool(l->log_scale, arg);
            } else if (str_eq_cstr(ident, "PointColoring") && viamd::extract_str(s, arg)) {
                l->point_coloring = str_eq_cstr(s, "structure") ? PointColoring_Structure : str_eq_cstr(s, "time") ? PointColoring_Time : PointColoring_Layer;
            } else if (str_eq_cstr(ident, "PointSize")) {
                viamd::extract_flt(l->point_size, arg);
            }
        }
        sanitize_layer(*l);
        if (!selected_layer) selected_layer = l->id;
    }

    static void sanitize_layer(Layer& l) {
        l.frames = (LayerFrames)CLAMP((int)l.frames, 0, LayerFrames_Count - 1);
        l.style = (LayerStyle)CLAMP((int)l.style, 0, LayerStyle_Count - 1);
        l.opacity = CLAMP(l.opacity, 0.0f, 1.0f);
        l.line_width = CLAMP(l.line_width, 0.5f, 8.0f);
        l.bandwidth = l.bandwidth > 0.0f ? CLAMP(l.bandwidth, 0.005f, MAX_BANDWIDTH) : 0.0f;
        l.point_size = CLAMP(l.point_size, 1.0f, 16.0f);
        l.num_levels = CLAMP(l.num_levels, 1, MAX_LEVELS);
        for (int i = 0; i < l.num_levels; ++i) l.levels[i] = CLAMP(l.levels[i], 0.01f, 0.999f);
        if (l.colormap < 0 || l.colormap >= ImPlot::GetColormapCount()) l.colormap = ImPlotColormap_Viridis;
    }

    // From before there were layers: the structures of the filter over every frame (or the period), a colour each
    void create_legacy_layers() {
        clear_layers();
        Layer* pts = add_layer(LayerFrames_All, LayerStyle_Points);
        if (!pts) return;
        // "all" was the default: all atoms as one structure, which an empty filter is
        if (strcmp(legacy.filter, "all") != 0) {
            snprintf(pts->filter.text, sizeof(pts->filter.text), "%s", legacy.filter);
        }
        pts->use_mass = legacy.use_mass;
        pts->point_coloring = PointColoring_Structure;
        // ImPlot's marker size is a radius
        pts->point_size = CLAMP(legacy.marker_size * 2.0f, 1.0f, 16.0f);
        uint32_t period = (uint32_t)MAX(legacy.period, 0);
        if (period == 0 && app) period = app->timeline.periods.legacy_filter_id;
        if (period != 0) {
            pts->frames = LayerFrames_Period;
            pts->period = period;
        }
        selected_layer = pts->id;
    }

    // ## Drawing

    void layer_label(char* buf, size_t cap, const Layer& l) const {
        if (l.name[0]) {
            snprintf(buf, cap, "%s", l.name);
            return;
        }
        const char* frames = "Current frame";
        char period_name[TIMELINE_PERIOD_NAME_CAP + 8];
        if (l.frames == LayerFrames_All) {
            frames = "All frames";
        } else if (l.frames == LayerFrames_Period) {
            const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, l.period);
            snprintf(period_name, sizeof(period_name), "%s", p ? p->name : "No period");
            frames = period_name;
        }
        if (l.filter.text[0]) snprintf(buf, cap, "%s \xc2\xb7 %s", frames, l.filter.text);
        else                  snprintf(buf, cap, "%s", frames);
    }

    void layer_stats(char* buf, size_t cap, const Layer& l) const {
        if (!l.structs_valid) {
            const LayerGroup* g = layer_group(l);
            if (g && layer_group_has_context(*g) && !g->valid) snprintf(buf, cap, "the group's subset does not evaluate");
            else snprintf(buf, cap, "%s", l.filter.text[0] && l.filter.result.error[0] ? l.filter.result.error : "no structures");
            return;
        }
        const size_t n = md_array_size(l.structures);
        const char* what = n == 1 ? "structure" : "structures";
        if (l.frames != LayerFrames_Current && run_num_frames(app) >= 2) {
            if (l.eval && !l.eval->complete()) {
                snprintf(buf, cap, "%zu %s \xc2\xb7 %u%% evaluated", n, what, (uint32_t)(100.0 * l.eval->frames_done.load() / MAX((size_t)1, l.eval->num_frames)));
                return;
            }
            uint32_t beg = 0, end = 0;
            if (layer_frame_range(l, &beg, &end)) {
                snprintf(buf, cap, "%zu %s \xc2\xb7 %u frames", n, what, end - beg);
                return;
            }
            snprintf(buf, cap, "%zu %s \xc2\xb7 no frames", n, what);
            return;
        }
        snprintf(buf, cap, "%zu %s", n, what);
    }

    bool layer_draw(DensityMapDraw* d, const Layer& l) const {
        if (!layer_draws_density(l) || !l.has_result || !l.tex || !l.eval || !l.eval->complete()) return false;
        if (l.result.samples[0] == 0 || l.result.max[0] <= 0.0f) return false;
        density_map_draw_init(d);
        d->tex = l.tex;
        d->opacity = l.opacity;
        d->line_width = l.line_width;
        d->line_color[0] = l.color.x;
        d->line_color[1] = l.color.y;
        d->line_color[2] = l.color.z;
        d->line_color[3] = l.color.w;
        d->mode = l.style == LayerStyle_Heatmap ? 1 : 0;
        d->fill = l.fill;
        d->lines = l.lines;
        if (l.style == LayerStyle_Heatmap) density_map_draw_colormap(d, l.colormap, l.log_scale);
        density_map_draw_levels(d, l.result, 0, l.levels, l.num_levels);
        density_map_draw_fills_from(d, l.color);
        d->content = l.shown_hash;
        return true;
    }

    uint32_t render_densities() {
        DensityMapDraw draws[MAX_LAYERS];
        int num_draws = 0;
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!layer_shown(l) || l.style == LayerStyle_Points || l.opacity <= 0.0f) continue;
            if (layer_draw(&draws[num_draws], l)) ++num_draws;
        }
        int w = 0, h = 0;
        density_map_plot_canvas_size(&w, &h);
        if (num_draws == 0 || w <= 0 || h <= 0) return 0;
        float view_uv[4];
        density_map_plot_view_uv(view_uv, shape_domain);
        // Within the triangle only: the smoothing spreads the density past its edges, where there are no shapes
        const DensityMapDomain& D = shape_domain;
        auto u = [&](vec2_t p) { return (p.x - D.x_min) / (D.x_max - D.x_min); };
        auto v = [&](vec2_t p) { return (p.y - D.y_min) / (D.y_max - D.y_min); };
        const float clip[6] = { u(tri_lin), v(tri_lin), u(tri_pla), v(tri_pla), u(tri_iso), v(tri_iso) };
        return density_map_composite(&canvas, w, h, view_uv, draws, num_draws, clip);
    }

    // ### Points

    struct PointSource {
        const vec2_t* coords = nullptr;     // [structure * stride + frame]
        size_t   stride = 0;
        size_t   num_structures = 0;
        uint32_t frame_beg = 0, frame_end = 0, frame_step = 1;
        bool     current = false;
    };

    bool point_source(PointSource* ps, Layer& l) const {
        if (l.style != LayerStyle_Points || !l.structs_valid) return false;
        if (l.frames == LayerFrames_Current || run_num_frames(app) < 2) {
            if (md_array_size(l.current_coords) == 0) return false;
            ps->coords = l.current_coords;
            ps->stride = 1;
            ps->num_structures = md_array_size(l.current_coords);
            ps->frame_beg = 0;
            ps->frame_end = 1;
            ps->frame_step = 1;
            ps->current = true;
            l.point_stride = 1;
            return true;
        }
        if (!l.eval) return false;
        uint32_t beg = 0, end = 0;
        if (!layer_frame_range(l, &beg, &end)) return false;
        end = MIN(end, (uint32_t)l.eval->num_frames);
        if (end <= beg) return false;
        ps->coords = l.eval->coords;
        ps->stride = l.eval->num_frames;
        ps->num_structures = l.eval->num_structures;
        ps->frame_beg = beg;
        ps->frame_end = end;
        const uint64_t total = (uint64_t)ps->num_structures * (end - beg);
        ps->frame_step = (uint32_t)MAX((uint64_t)1, (total + MAX_POINTS - 1) / MAX_POINTS);
        ps->current = false;
        l.point_stride = ps->frame_step;
        return true;
    }

    template <typename F>
    static void for_each_point(const PointSource& ps, F&& f) {
        for (size_t s = 0; s < ps.num_structures; ++s) {
            const vec2_t* c = ps.coords + s * ps.stride;
            for (uint32_t frame = ps.frame_beg; frame < ps.frame_end; frame += ps.frame_step) {
                const vec2_t p = ps.current ? c[0] : c[frame];
                if (coord_valid(p)) f((uint32_t)s, frame, p);
            }
        }
    }

    ImVec4 point_color(const Layer& l, const PointSource& ps, uint32_t s, uint32_t frame) const {
        switch (l.point_coloring) {
        case PointColoring_Structure:
            return ImPlot::GetColormapColor((int)s % ImPlot::GetColormapSize(ImPlotColormap_Deep), ImPlotColormap_Deep);
        case PointColoring_Time:
            if (ps.current || ps.frame_end - ps.frame_beg <= 1) return l.color;
            return ImPlot::SampleColormap((float)(frame - ps.frame_beg) / (float)(ps.frame_end - ps.frame_beg - 1), l.colormap);
        default:
            return l.color;
        }
    }

    // The atoms of the structure s of layer l
    md_bitset_t structure_atoms(const Layer& l, const PointSource& ps, uint32_t s) const {
        if (!ps.current && l.eval && s < l.eval->num_structures) return l.eval->structures[s];
        if (s < md_array_size(l.structures)) return l.structures[s];
        return md_bitset_empty_set();
    }

    struct Hover {
        int      layer = -1;
        uint32_t structure = 0;
        uint32_t frame = 0;
        bool     current = true;
        vec2_t   coord = {};
    };

    void draw_tooltip(const Hover& h, const PointSource& ps) const {
        const Layer& l = layers[h.layer];
        const size_t ns = ps.num_structures;
        const vec3_t w = shape_weights_of(h.coord);
        const md_bitset_t atoms = structure_atoms(l, ps, h.structure);
        ImGui::BeginTooltip();
        if (ns > 1) ImGui::Text("Structure %u of %zu  \xc2\xb7  %u atoms", h.structure + 1, ns, (uint32_t)md_bitset_count(atoms));
        else        ImGui::Text("%u atoms", (uint32_t)md_bitset_count(atoms));
        ImGui::Text("linear %.2f   planar %.2f   spherical %.2f", w.x, w.y, w.z);
        if (h.current) ImGui::TextDisabled("frame shown");
        else           ImGui::TextDisabled("frame %u  \xc2\xb7  click to go there", h.frame);
        char label[128];
        layer_label(label, sizeof(label), l);
        ImGui::TextDisabled("%s", label);
        ImGui::EndTooltip();
    }

    void draw_points(bool hovered) {
        plot_drawn = true;
        md_bitmask_t* selection_mask = &app->selection.selection_mask;
        md_bitmask_t* highlight_mask = &app->selection.highlight_mask;

        const ImPlotRect lim = ImPlot::GetPlotLimits();
        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 size = ImPlot::GetPlotSize();
        const double sx = size.x / lim.Size().x;
        const double sy = size.y / lim.Size().y;
        auto to_pixels = [&](vec2_t c) {
            return ImVec2((float)(pos.x + (c.x - lim.X.Min) * sx), (float)(pos.y + (lim.Y.Max - c.y) * sy));
        };

        PointSource sources[MAX_LAYERS];
        bool has_source[MAX_LAYERS] = {};
        for (int i = 0; i < num_layers; ++i) {
            layers[i].point_count = 0;
            if (layer_shown(layers[i])) has_source[i] = point_source(&sources[i], layers[i]);
        }

        // What the pointer is on: the nearest point within reach, of the layer on top
        Hover hover;
        if (hovered && !box.active) {
            const ImVec2 mouse = ImGui::GetMousePos();
            for (int i = 0; i < num_layers && hover.layer < 0; ++i) {
                if (!has_source[i]) continue;
                const float r = MAX(8.0f, layers[i].point_size);
                float best = r * r;
                for_each_point(sources[i], [&](uint32_t s, uint32_t frame, vec2_t c) {
                    const ImVec2 p = to_pixels(c);
                    const float d2 = (p.x - mouse.x) * (p.x - mouse.x) + (p.y - mouse.y) * (p.y - mouse.y);
                    if (d2 < best) {
                        best = d2;
                        hover = Hover{ i, s, frame, sources[i].current, c };
                    }
                });
            }
        }

        // Not moved from where it was pressed: a click, on the point pressed on (if any)
        bool box_is_click = false;
        if (box.active) {
            const ImVec2 a = ImPlot::PlotToPixels(box.p0);
            const ImVec2 m = ImGui::GetMousePos();
            box_is_click = fabsf(a.x - m.x) < 3.0f && fabsf(a.y - m.y) < 3.0f;
        }

        // The highlight: what the box covers, or the structure under the pointer
        if (hovered || box.active) {
            md_bitmask_clear(highlight_mask);
            highlights_now = true;
        }
        if (box.active && box_is_click) {
            if (box.click_layer >= 0 && box.click_layer < num_layers && box.click_structure >= 0 && has_source[box.click_layer]) {
                md_bitmask_or(highlight_mask, structure_atoms(layers[box.click_layer], sources[box.click_layer], (uint32_t)box.click_structure));
            }
        } else if (box.active) {
            const double x0 = MIN(box.p0.x, box.p1.x), x1 = MAX(box.p0.x, box.p1.x);
            const double y0 = MIN(box.p0.y, box.p1.y), y1 = MAX(box.p0.y, box.p1.y);
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                // A structure once, however many of its points are in
                md_allocator_i* temp = md_get_heap_allocator();
                md_bitmask_t taken = md_bitmask_create(0, (uint32_t)sources[i].num_structures, temp);
                for_each_point(sources[i], [&](uint32_t s, uint32_t, vec2_t c) {
                    if (x0 <= c.x && c.x <= x1 && y0 <= c.y && c.y <= y1 && !md_bitmask_test_bit(&taken, s)) {
                        md_bitmask_set_bit(&taken, s);
                        md_bitmask_or(highlight_mask, structure_atoms(layers[i], sources[i], s));
                    }
                });
                md_bitmask_free(&taken);
            }
        } else if (hover.layer >= 0) {
            md_bitmask_or(highlight_mask, structure_atoms(layers[hover.layer], sources[hover.layer], hover.structure));
        }

        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        density_map_marks_begin(dl);
        for (int i = num_layers - 1; i >= 0; --i) {
            if (!has_source[i]) continue;
            Layer& l = layers[i];
            const PointSource& ps = sources[i];
            const float r = l.point_size * 0.5f;
            for_each_point(ps, [&](uint32_t s, uint32_t frame, vec2_t c) {
                ImVec4 col = point_color(l, ps, s, frame);
                col.w *= l.opacity;
                const ImU32 u = ImGui::ColorConvertFloat4ToU32(col);
                l.point_count += 1;
                if (ps.current) density_map_mark(dl, to_pixels(c), r + 1.0f, DensityMapMark_Disc, u);
                else            density_map_mark(dl, to_pixels(c), r, DensityMapMark_Square, u);
            });
        }
        // Rings around the points of the frame shown of what is selected and highlighted
        const vec4_t sel_col = app->selection.color.selection.outline;
        const vec4_t hl_col = app->selection.color.highlight.outline;
        const bool any_sel = !md_bitmask_empty(selection_mask);
        const bool any_hl = !md_bitmask_empty(highlight_mask);
        for (int i = num_layers - 1; i >= 0 && (any_sel || any_hl); --i) {
            if (!has_source[i] || !sources[i].current) continue;
            const Layer& l = layers[i];
            const float r = MAX(l.point_size * 0.5f + 3.0f, 5.0f);
            for_each_point(sources[i], [&](uint32_t s, uint32_t, vec2_t c) {
                const md_bitset_t atoms = structure_atoms(l, sources[i], s);
                const uint32_t first = md_bitset_first(atoms);
                if (md_bitset_empty(atoms)) return;
                const bool hl = any_hl && md_bitmask_test_bit(highlight_mask, first);
                const bool sel = any_sel && md_bitmask_test_bit(selection_mask, first);
                if (!hl && !sel) return;
                density_map_mark(dl, to_pixels(c), r, DensityMapMark_Ring, ImGui::ColorConvertFloat4ToU32(vec_cast(hl ? hl_col : sel_col)));
            });
        }
        // The point under the pointer
        if (hover.layer >= 0) {
            density_map_mark(dl, to_pixels(hover.coord), MAX(layers[hover.layer].point_size * 0.5f + 3.0f, 5.0f), DensityMapMark_Ring, ImGui::ColorConvertFloat4ToU32(vec_cast(hl_col)));
        }
        density_map_marks_end(dl);
        ImPlot::PopPlotClipRect();

        if (hover.layer >= 0) {
            draw_tooltip(hover, sources[hover.layer]);
            // A click on a point of a frame: that frame
            if (!hover.current && !ImGui::GetIO().KeyShift && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                app->animation.frame = hover.frame;
            }
        }

        // Selecting: shift-drag a box, shift-click a point; the right button takes away; a right click on nothing clears
        const bool shift = ImGui::GetIO().KeyShift;
        const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
        if (!box.active && hovered && shift && (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
            box.active = true;
            box.op = ImGui::IsMouseClicked(ImGuiMouseButton_Left) ? SelectionOperator::Or : SelectionOperator::AndNot;
            box.p0 = box.p1 = mouse;
            box.click_layer = hover.layer;
            box.click_structure = hover.layer >= 0 ? (int64_t)hover.structure : -1;
        }
        if (box.active) {
            box.p1 = mouse;
            if (!box_is_click) {
                const ImVec2 a = ImPlot::PlotToPixels(box.p0);
                const ImVec2 b = ImPlot::PlotToPixels(box.p1);
                const vec4_t c = box.op == SelectionOperator::Or ? sel_col : vec4_t{1.0f, 0.4f, 0.4f, 0.9f};
                ImPlot::PushPlotClipRect();
                dl->AddRectFilled(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 0.12f)));
                dl->AddRect(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(vec_cast(c)), 0.0f, 0, 1.0f);
                ImPlot::PopPlotClipRect();
            }
            const ImGuiMouseButton button = box.op == SelectionOperator::Or ? ImGuiMouseButton_Left : ImGuiMouseButton_Right;
            if (ImGui::IsMouseReleased(button) || !ImGui::IsMouseDown(button)) {
                if (!box_is_click || box.click_structure >= 0) {
                    grow_mask_by_selection_granularity(highlight_mask, app->selection.granularity, app->mold.sys);
                    modify_field(selection_mask, md_bitmask_view(highlight_mask), box.op);
                } else if (box.op == SelectionOperator::AndNot) {
                    md_bitmask_clear(selection_mask);
                }
                box.active = false;
                box.click_layer = -1;
                box.click_structure = -1;
            }
        }
    }

    // The triangle, its grid (lines of equal weight, a quarter apart) and its corners
    void draw_triangle() {
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        auto px = [](vec2_t p) { return ImPlot::PlotToPixels(ImPlotPoint(p.x, p.y)); };
        const ImVec2 a = px(tri_lin), b = px(tri_pla), c = px(tri_iso);
        ImPlot::PushPlotClipRect();
        dl->AddTriangleFilled(a, b, c, IM_COL32(255, 255, 255, 14));
        if (show_grid) {
            const ImU32 grid = IM_COL32(255, 255, 255, 28);
            for (int k = 1; k < 4; ++k) {
                const float t = k / 4.0f;
                // Of equal weight of each corner: parallel to the side opposite it
                dl->AddLine(ImLerp(a, b, t), ImLerp(a, c, t), grid);
                dl->AddLine(ImLerp(b, a, t), ImLerp(b, c, t), grid);
                dl->AddLine(ImLerp(c, a, t), ImLerp(c, b, t), grid);
            }
        }
        dl->AddTriangle(a, b, c, IM_COL32(255, 255, 255, 70), 1.0f);
        const ImU32 text = IM_COL32(255, 255, 255, 220);
        const char* t_lin = "Linear";
        const char* t_pla = "Planar";
        const char* t_iso = "Spherical";
        dl->AddText(c + ImGui::CalcTextSize(t_iso) * ImVec2(-0.5f, -1.2f), text, t_iso);
        dl->AddText(b + ImGui::CalcTextSize(t_pla) * ImVec2(-0.5f, 0.2f), text, t_pla);
        dl->AddText(a + ImGui::CalcTextSize(t_lin) * ImVec2(-0.5f, 0.2f), text, t_lin);
        ImPlot::PopPlotClipRect();
    }

    // A layer over a period, in the group of what is selected (or in 'group' when given)
    Layer* add_period_layer(uint32_t period_id, uint32_t group = UINT32_MAX) {
        const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, period_id);
        if (!p) return nullptr;
        Layer* l = add_like(LayerFrames_Period, LayerStyle_Contours, group == UINT32_MAX ? current_group() : group);
        if (!l) return nullptr;
        l->period = period_id;
        l->color = ImVec4(p->color.x, p->color.y, p->color.z, 1.0f);
        return l;
    }

    // A layer first in the group 'group', selected. In a group it shows the group's structures; at the top level it
    // is of the structures of the layer selected, if any. The weighting is that of the layer selected either way.
    Layer* add_like(LayerFrames frames, LayerStyle style, uint32_t group) {
        const Layer* like = find_layer(selected_layer);
        char filter[sizeof(AtomFilter::text)] = "";
        bool use_mass = true;
        if (like) {
            if (group == 0 && like->node.parent == 0) MEMCPY(filter, like->filter.text, sizeof(filter));
            use_mass = like->use_mass;
        }
        Layer* l = add_layer(frames, style, 0, group);
        if (!l) return nullptr;
        MEMCPY(l->filter.text, filter, sizeof(filter));
        l->use_mass = use_mass;
        selected_layer = l->id;
        selected_group = 0;
        return l;
    }

    // The view of the triangle and its corners' labels, centred, at the plot's shape (pixels)
    static ImPlotRect fit_rect(ImVec2 size) {
        const float line = ImGui::GetTextLineHeight();
        const float mx = ImGui::CalcTextSize("Spherical").x * 0.5f + 8.0f;    // The labels are centred on the corners
        const float mt = line * 1.6f + 6.0f;
        const float mb = line * 1.4f + 6.0f;
        const double h_tri = tri_iso.y;
        const double s = MAX(1.0, MIN((size.x - 2.0 * mx) / 1.0, (size.y - mt - mb) / h_tri));     // Pixels per unit
        const double w = size.x / s;
        const double h = size.y / s;
        const double extra = (size.y - s * h_tri - mt - mb) * 0.5;
        const double y0 = -(mb + MAX(0.0, extra)) / s;
        return ImPlotRect(0.5 - w * 0.5, 0.5 + w * 0.5, y0, y0 + h);
    }

    void draw_plot() {
        constexpr ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
        constexpr ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoDecorations;

        ImPlotInputMap& input = ImPlot::GetInputMap();
        const ImPlotInputMap saved_input = input;
        input.OverrideMod = ImGuiMod_Shift;     // Shift is for selecting
        ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(1, 1));
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));

        if (ImPlot::BeginPlot("##Shapespace", ImVec2(-1, -1), flags)) {
            // Fitted to the shape of the plot of the last frame: ImPlot's equal axes would otherwise keep the width
            // asked for and cut the triangle at the top and bottom of a wide plot
            if (fit_view && plot_size.x > 0 && plot_size.y > 0) {
                fitted_rect = fit_rect(plot_size);
                ImPlot::SetupAxesLimits(fitted_rect.X.Min, fitted_rect.X.Max, fitted_rect.Y.Min, fitted_rect.Y.Max, ImGuiCond_Always);
                fit_view = false;
                view_fitted = true;
            } else {
                ImPlot::SetupAxesLimits(-0.10, 1.10, -0.10, 0.98, ImGuiCond_Once);
            }
            ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, -4.0, 5.0);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, -4.0, 5.0);
            ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 0.02, 8.0);
            ImPlot::SetupAxisZoomConstraints(ImAxis_Y1, 0.02, 8.0);
            ImPlot::SetupAxes(nullptr, nullptr, axis_flags, axis_flags);
            ImPlot::SetupFinish();

            draw_triangle();
            density_map_plot_canvas(render_densities());
            draw_points(ImPlot::IsPlotHovered());

            if (ImPlot::BeginDragDropTargetPlot()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                    add_period_layer(*(const uint32_t*)payload->Data);
                }
                ImPlot::EndDragDropTarget();
            }
            // Double-click: fitted (ImPlot's own fit has no items to fit to)
            if (ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) fit_view = true;

            // Panned or zoomed: no longer fitted. Fitted and of another shape: fitted again.
            const ImVec2 size = ImPlot::GetPlotSize();
            const ImPlotRect lim = ImPlot::GetPlotLimits();
            const bool reshaped = fabsf(size.x - plot_size.x) > 0.5f || fabsf(size.y - plot_size.y) > 0.5f;
            if (view_fitted && !fit_view) {
                const double eps = 1e-6 * MAX(1e-6, fitted_rect.X.Size());
                const bool moved = fabs(lim.X.Min - fitted_rect.X.Min) > eps || fabs(lim.X.Max - fitted_rect.X.Max) > eps ||
                                   fabs(lim.Y.Min - fitted_rect.Y.Min) > eps || fabs(lim.Y.Max - fitted_rect.Y.Max) > eps;
                if (reshaped) fit_view = true;
                else if (moved) view_fitted = false;
            }
            if (plot_size.x <= 0.0f) fit_view = true;
            plot_size = size;
            ImPlot::EndPlot();
        }

        ImPlot::PopStyleColor();
        ImPlot::PopStyleVar();
        input = saved_input;
    }

    // ### The window

    void export_layer(const Layer& l, str_t ext) {
        const EvaluationRef& eval = l.eval;
        if (!eval || !eval->complete()) return;
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) return;
        const str_t path = { path_buf, strnlen(path_buf, sizeof(path_buf)) };

        md_temp_scope_t temp_scope = md_temp_begin();
        md_allocator_i* temp = md_temp_allocator(temp_scope);
        defer { md_temp_end(temp_scope); };

        // The time column is timeline.x_values, held in display units: its label is from the same place
        const md_unit_t time_unit = app->timeline.time_unit;
        str_t x_label = STR_INIT("Frame");
        if (!md_unit_is_none(time_unit)) {
            char unit_buf[64];
            md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
            x_label = str_printf(temp, "Time (%s)", unit_buf);
        }
        md_array(str_t)  labels = 0;
        md_array(float*) values = 0;
        md_array_push(labels, x_label, temp);
        md_array_push(values, NULL, temp);
        const size_t ns = eval->num_structures;
        const size_t nf = eval->num_frames;
        for (size_t s = 0; s < ns; ++s) {
            md_array_push(labels, ns > 1 ? str_printf(temp, "linear %zu", s + 1)    : STR_LIT("linear"),    temp);
            md_array_push(labels, ns > 1 ? str_printf(temp, "planar %zu", s + 1)    : STR_LIT("planar"),    temp);
            md_array_push(labels, ns > 1 ? str_printf(temp, "isotropic %zu", s + 1) : STR_LIT("isotropic"), temp);
            for (int k = 0; k < 3; ++k) md_array_push(values, NULL, temp);
        }
        const float* time = app->timeline.x_values;
        const bool has_time = time && md_array_size(app->timeline.x_values) >= nf;
        for (size_t f = 0; f < nf; ++f) {
            md_array_push(values[0], has_time ? time[f] : (float)f, temp);
            for (size_t s = 0; s < ns; ++s) {
                const vec3_t w = shape_weights_of(eval->coords[s * nf + f]);
                md_array_push(values[1 + s * 3 + 0], w.x, temp);
                md_array_push(values[1 + s * 3 + 1], w.y, temp);
                md_array_push(values[1 + s * 3 + 2], w.z, temp);
            }
        }
        const size_t num_cols = ns * 3 + 1;
        if (str_eq(ext, STR_LIT("csv"))) {
            md_csv_write_to_file(values, labels, num_cols, nf, path);
        } else {
            const str_t header = md_xvg_format_header(STR_LIT("Shape Space"), x_label, STR_LIT("Shape Weight"), num_cols - 1, labels + 1, temp);
            const str_t xvg = md_xvg_format(header, num_cols, nf, values, temp);
            md_file_t file = {0};
            if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                md_file_write(file, xvg.ptr, xvg.len);
                md_file_close(&file);
            }
        }
    }

    void draw_toolbar() {
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6, 2));
        if (ImGui::Button(ICON_FA_LAYER_GROUP "##layers")) show_sidebar = !show_sidebar;
        ImGui::SetItemTooltip(show_sidebar ? "Hide the layers" : "Show the layers");
        ImGui::SameLine();

        const Layer* sel = find_layer(selected_layer);
        const bool can_export = sel && sel->frames != LayerFrames_Current && sel->eval && sel->eval->complete();
        if (ImGui::Button("Export " ICON_FA_CARET_DOWN)) ImGui::OpenPopup("##export");
        if (ImGui::BeginPopup("##export")) {
            ImGui::TextDisabled("The shapes of the selected layer's structures, every frame");
            if (ImGui::MenuItem("CSV", nullptr, false, can_export)) export_layer(*sel, STR_LIT("csv"));
            if (ImGui::MenuItem("XVG", nullptr, false, can_export)) export_layer(*sel, STR_LIT("xvg"));
            if (!can_export) ImGui::TextDisabled("Select a layer over frames, once evaluated");
            ImGui::EndPopup();
        }

        ImGui::SameLine();
        const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + ImGui::CalcTextSize(ICON_FA_GEAR).x + ImGui::GetStyle().FramePadding.x * 4 + ImGui::GetStyle().ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
        if (ImGui::Button(ICON_FA_EXPAND)) fit_view = true;
        ImGui::SetItemTooltip("Fit the view (or double-click the plot)");
        ImGui::SameLine();
        if (ImGui::Button(ICON_FA_GEAR)) ImGui::OpenPopup("##settings");
        if (ImGui::BeginPopup("##settings")) {
            ImGui::Checkbox("Grid of equal weights", &show_grid);
            ImGui::Separator();
            ImGui::TextDisabled("Click a point of a frame to go there.");
            ImGui::TextDisabled("Shift-drag to select structures, with the right button to deselect them.");
            ImGui::TextDisabled("Drop a period from the timeline on the plot to add its density.");
            ImGui::EndPopup();
        }
        ImGui::PopStyleVar();
    }

    // What the rows asked for, done once the list is drawn
    struct SidebarActions {
        uint32_t remove = 0;
        uint32_t duplicate = 0;
        uint32_t group = 0;             // A new group of this layer
        uint32_t to_group_layer = 0;    // This layer into to_group (0: out of any)
        uint32_t to_group = 0;
        uint32_t group_id = 0;
        LayerGroupAction group_action = LayerGroupAction_None;
        bool     move = false;
        LayerNodeRef move_src, move_target;
        RepresentationDrop move_drop = RepresentationDrop_Before;
        uint32_t period = 0;            // A layer over it in period_group
        uint32_t period_group = 0;
    };

    void draw_layer_row(Layer& l, float indent, RepresentationNode* const nodes[], size_t num_nodes, SidebarActions* act) {
        ImGui::PushID("layer");
        ImGui::PushID((int)l.id);
        const float line_h = ImGui::GetTextLineHeight();
        const float row_h = line_h * 2.0f + 8.0f;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + indent);
        const ImVec2 row_min = ImGui::GetCursorScreenPos();
        const float row_w = ImGui::GetContentRegionAvail().x;
        const ImRect row_rect(row_min, ImVec2(row_min.x + row_w, row_min.y + row_h));

        char label[128], stats[128];
        layer_label(label, sizeof(label), l);
        layer_stats(stats, sizeof(stats), l);

        ImGui::SetNextItemAllowOverlap();
        if (ImGui::Selectable("##row", selected_layer == l.id && !selected_group, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h))) {
            selected_layer = l.id;
            selected_group = 0;
        }
        const bool row_hovered = ImGui::IsItemHovered();
        layer_tree_drag_source({LIST_ID, l.id, false}, label);
        if (ImGui::BeginPopupContextItem("##menu")) {
            selected_layer = l.id;
            selected_group = 0;
            if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false, num_layers < MAX_LAYERS)) act->duplicate = l.id;
            if (l.node.parent == 0) {
                if (ImGui::MenuItem(ICON_FA_FOLDER " Group", nullptr, false, num_groups < LAYER_MAX_GROUPS)) act->group = l.id;
                ImGui::SetItemTooltip("A new group of it, with its structures as the group's: the layers added to the group show them");
            }
            if (ImGui::BeginMenu("Move to", num_groups > 0)) {
                if (l.node.parent != 0 && ImGui::MenuItem("Out of the group")) {
                    act->to_group_layer = l.id;
                    act->to_group = 0;
                }
                for (int gi = 0; gi < num_groups; ++gi) {
                    if (groups[gi].node.id == l.node.parent) continue;
                    char gl[96];
                    layer_group_label(gl, sizeof(gl), groups[gi]);
                    ImGui::PushID(gi);
                    if (ImGui::MenuItem(gl)) {
                        act->to_group_layer = l.id;
                        act->to_group = groups[gi].node.id;
                    }
                    ImGui::PopID();
                }
                ImGui::EndMenu();
            }
            ImGui::Separator();
            if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) act->remove = l.id;
            ImGui::EndPopup();
        }
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                l.frames = LayerFrames_Period;
                l.period = *(const uint32_t*)payload->Data;
                selected_layer = l.id;
                selected_group = 0;
            }
            ImGui::EndDragDropTarget();
        }
        LayerNodeRef src;
        RepresentationDrop drop;
        if (layer_tree_drop_target(row_rect, row_min.x, {LIST_ID, l.id, false}, false, nodes, num_nodes, &src, &drop)) {
            act->move = true;
            act->move_src = src;
            act->move_target = {LIST_ID, l.id, false};
            act->move_drop = drop;
        }
        if (row_hovered && l.frames == LayerFrames_Period) timeline_periods_highlight(l.period);
        // Hovering a layer lights its structures up
        if (row_hovered && l.structs_valid) {
            highlights_now = true;
            md_bitmask_clear(&app->selection.highlight_mask);
            for (size_t i = 0; i < md_array_size(l.structures); ++i) md_bitmask_or(&app->selection.highlight_mask, l.structures[i]);
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float pad = 4.0f;
        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, row_min.y + (row_h - ImGui::GetFrameHeight()) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        if (ImGui::Button(l.visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) l.visible = !l.visible;
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip(l.visible ? "Hide" : "Show");

        const ImVec2 sw_min = ImVec2(row_min.x + pad + line_h + 8, row_min.y + (row_h - line_h) * 0.5f);
        const ImVec2 sw_max = ImVec2(sw_min.x + line_h, sw_min.y + line_h);
        const DensityMapSwatch kind = l.style == LayerStyle_Points ? DensityMapSwatch_Points : l.style == LayerStyle_Contours ? DensityMapSwatch_Contours : DensityMapSwatch_Heatmap;
        const bool by_map = l.style == LayerStyle_Points && l.point_coloring == PointColoring_Time;
        density_map_swatch(dl, sw_min, sw_max, by_map ? DensityMapSwatch_Heatmap : kind, l.color, l.fill, l.colormap, l.visible ? 1.0f : 0.35f);

        const ImVec2 text_min = ImVec2(sw_max.x + 6, row_min.y + 4);
        const ImVec4 clip = ImVec4(text_min.x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), text_min, ImGui::GetColorU32(l.visible ? ImGuiCol_Text : ImGuiCol_TextDisabled), label, nullptr, 0.0f, &clip);
        // On the selection's colour the disabled one is too faint to read
        const bool selected = selected_layer == l.id && !selected_group;
        const ImU32 dim_col = selected ? ImGui::GetColorU32(ImGuiCol_Text, 0.7f) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(text_min.x, text_min.y + line_h + 1),
                    !l.structs_valid ? IM_COL32(255, 110, 110, 255) : dim_col, stats, nullptr, 0.0f, &clip);
        // Working on it: its evaluation or its density
        if (l.job || (l.eval && !l.eval->complete())) {
            const float t = (float)ImGui::GetTime();
            const ImVec2 c = ImVec2(row_min.x + row_w - pad - 6, row_min.y + 9);
            dl->PathArcTo(c, 4.5f, t * 6.0f, t * 6.0f + 4.0f, 10);
            dl->PathStroke(ImGui::GetColorU32(ImGuiCol_TextDisabled), 0, 1.5f);
        }

        ImGui::SetCursorScreenPos(ImVec2(row_min.x - indent, row_min.y + row_h + 2));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
        ImGui::PopID();
    }

    void draw_group(LayerGroup& g, RepresentationNode* const nodes[], size_t num_nodes, SidebarActions* act) {
        int num_in = 0;
        for (int i = 0; i < num_layers; ++i) num_in += layers[i].node.parent == g.node.id ? 1 : 0;
        const bool context = layer_group_has_context(g);
        const bool invalid = context && !g.valid;
        char summary[64];
        if (invalid) {
            snprintf(summary, sizeof(summary), "does not evaluate");
        } else if (context) {
            const size_t ns = md_array_size(g.sets);
            snprintf(summary, sizeof(summary), "%zu structure%s", ns, ns == 1 ? "" : "s");
        } else {
            snprintf(summary, sizeof(summary), "%d layer%s", num_in, num_in == 1 ? "" : "s");
        }

        // What ungrouping would not keep: a layer has no selection of its own to take a group's on
        bool exact = !g.follow_selection;
        for (int i = 0; i < num_layers && exact; ++i) {
            const Layer& l = layers[i];
            if (l.node.parent == g.node.id) exact = layer_group_keeps_filter(l.filter.text, false, g);
        }
        const char* warning = exact ? nullptr :
            "Some of its layers cannot be given a filter that keeps what they show: those within the selection, or with "
            "one of their own within a filter that depends on the coordinates. Those keep their own, over the whole system.";
        const LayerGroupRow row = layer_group_row(g, LIST_ID, selected_group == g.node.id, summary, invalid, warning);
        if (row.clicked) {
            selected_group = g.node.id;
            selected_layer = 0;
        }
        if (row.action != LayerGroupAction_None) {
            act->group_id = g.node.id;
            act->group_action = row.action;
        }
        // Its structures lit up
        if (row.hovered && context && g.valid) {
            highlights_now = true;
            md_bitmask_clear(&app->selection.highlight_mask);
            md_bitmask_or(&app->selection.highlight_mask, g.atoms);
        }

        ImGui::PushID("group_drop");
        ImGui::PushID((int)g.node.id);
        LayerNodeRef src;
        RepresentationDrop drop;
        if (layer_tree_drop_target(row.rect, row.x0, {LIST_ID, g.node.id, true}, g.open, nodes, num_nodes, &src, &drop)) {
            act->move = true;
            act->move_src = src;
            act->move_target = {LIST_ID, g.node.id, true};
            act->move_drop = drop;
        }
        // A period from the timeline: a layer over it, in the group
        if (ImGui::BeginDragDropTargetCustom(row.rect, ImGui::GetID("##period"))) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                act->period = *(const uint32_t*)payload->Data;
                act->period_group = g.node.id;
            }
            ImGui::EndDragDropTarget();
        }
        ImGui::PopID();
        ImGui::PopID();

        if (!g.open) return;
        // What is within a hidden group is dimmed: hidden too, whatever its own eye says
        if (!g.visible) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t nk = representation_tree_children(kids, ARRAY_SIZE(kids), nodes, num_nodes, g.node.id);
        for (size_t k = 0; k < nk; ++k) {
            if (kids[k] < num_layers) draw_layer_row(layers[kids[k]], ImGui::GetStyle().IndentSpacing, nodes, num_nodes, act);
        }
        if (nk == 0) {
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + ImGui::GetStyle().IndentSpacing);
            ImGui::TextDisabled("Empty: drag layers here, or add with +");
            const ImRect r(ImGui::GetItemRectMin(), ImVec2(ImGui::GetItemRectMin().x + ImGui::GetContentRegionAvail().x, ImGui::GetItemRectMax().y));
            ImGui::PushID((int)g.node.id);
            LayerNodeRef src;
            if (layer_tree_drop_target_end(r, r.Min.x, LIST_ID, g.node.id, nodes, num_nodes, &src)) {
                act->move = true;
                act->move_src = src;
                act->move_target = {LIST_ID, g.node.id, true};
                act->move_drop = RepresentationDrop_IntoLast;
            }
            ImGui::PopID();
        }
        if (!g.visible) ImGui::PopStyleVar();
    }

    // What is added goes first in the group of what is selected: a layer there shows the group's structures
    void draw_add_menu() {
        const uint32_t group = current_group();
        if (const LayerGroup* g = find_group(group)) {
            char gl[96];
            layer_group_label(gl, sizeof(gl), *g);
            ImGui::TextDisabled("In %s", gl);
        }
        const bool room = num_layers < MAX_LAYERS;
        if (ImGui::MenuItem("Current frame", nullptr, false, room)) {
            if (Layer* l = add_like(LayerFrames_Current, LayerStyle_Points, group)) l->point_size = 7.0f;
        }
        if (ImGui::MenuItem("All frames", nullptr, false, room)) add_like(LayerFrames_All, LayerStyle_Contours, group);
        const TimelinePeriods& ps = app->timeline.periods;
        if (ImGui::BeginMenu("Period", ps.count > 0 && room)) {
            for (int i = 0; i < ps.count; ++i) {
                char label[TIMELINE_PERIOD_NAME_CAP + 16];
                snprintf(label, sizeof(label), "%s###p%u", ps.items[i].name, ps.items[i].id);
                if (ImGui::MenuItem(label)) add_period_layer(ps.items[i].id, group);
                if (ImGui::IsItemHovered()) timeline_periods_highlight(ps.items[i].id);
            }
            ImGui::EndMenu();
        }
        if (ps.count == 0) ImGui::SetItemTooltip("Mark periods in the timeline (shift-drag) to show the shapes over them");
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_FA_FOLDER " Group", nullptr, false, num_groups < LAYER_MAX_GROUPS)) {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            if (LayerGroup* ng = add_group(layer_tree_first_order(nodes, tree_nodes(nodes), 0))) {
                selected_group = ng->node.id;
                selected_layer = 0;
            }
        }
        ImGui::SetItemTooltip("A group of layers that show the same structures: its subset, which the layers in it are within");
    }

    void draw_inspector(Layer& l) {
        prop::begin("shapespace");
        ImGui::PushID((int)l.id);

        char auto_name[128];
        {
            char saved[sizeof(l.name)];
            MEMCPY(saved, l.name, sizeof(saved));
            l.name[0] = '\0';
            layer_label(auto_name, sizeof(auto_name), l);
            MEMCPY(l.name, saved, sizeof(saved));
        }
        prop::row("name");
        ImGui::InputTextWithHint("##name", auto_name, l.name, sizeof(l.name));

        prop::row("structures", "A structure for each set of atoms the filter gives: residue(:) a structure per residue. Empty: all atoms, as one.");
        {
            const bool invalid = !l.structs_valid && l.filter.text[0];
            if (invalid) ImGui::PushInvalid();
            const LayerGroup* g = layer_group(l);
            const bool in_context = g && layer_group_has_context(*g) && g->valid;
            char none[96];
            if (in_context) {
                char gl[64];
                layer_group_label(gl, sizeof(gl), *g);
                snprintf(none, sizeof(none), "the group's (%s)", gl);
            } else {
                snprintf(none, sizeof(none), "all atoms, as one");
            }
            filter_input("##filter", &l.filter, layer_group_filter_context(app, in_context ? g : nullptr), 0, nullptr, none);
            if (invalid) ImGui::PopInvalid();
            if (ImGui::IsItemHovered() && l.structs_valid) {
                highlights_now = true;
                md_bitmask_clear(&app->selection.highlight_mask);
                for (size_t i = 0; i < md_array_size(l.structures); ++i) md_bitmask_or(&app->selection.highlight_mask, l.structures[i]);
            }
            if (l.filter.result.is_dynamic && l.frames != LayerFrames_Current) {
                prop::hint("Depends on the coordinates: taken in the frame it was evaluated in");
            }
        }
        prop::row("by mass", "Atoms weighted by their mass, or alike");
        ImGui::Checkbox("##mass", &l.use_mass);

        prop::row("frames");
        {
            char preview[TIMELINE_PERIOD_NAME_CAP + 16];
            const TimelinePeriods& ps = app->timeline.periods;
            const TimelinePeriod* cur = timeline_periods_find(&ps, l.period);
            if (l.frames == LayerFrames_Current)  snprintf(preview, sizeof(preview), "Current frame");
            else if (l.frames == LayerFrames_All) snprintf(preview, sizeof(preview), "All frames");
            else snprintf(preview, sizeof(preview), "%s", cur ? cur->name : "Period (missing)");
            if (ImGui::BeginCombo("##frames", preview)) {
                if (ImGui::Selectable("Current frame", l.frames == LayerFrames_Current)) l.frames = LayerFrames_Current;
                if (ImGui::Selectable("All frames", l.frames == LayerFrames_All)) l.frames = LayerFrames_All;
                if (ps.count > 0) ImGui::SeparatorText("Periods");
                for (int i = 0; i < ps.count; ++i) {
                    char label[TIMELINE_PERIOD_NAME_CAP + 16];
                    snprintf(label, sizeof(label), "%s###p%u", ps.items[i].name, ps.items[i].id);
                    if (ImGui::Selectable(label, l.frames == LayerFrames_Period && l.period == ps.items[i].id)) {
                        l.frames = LayerFrames_Period;
                        l.period = ps.items[i].id;
                    }
                    if (ImGui::IsItemHovered()) timeline_periods_highlight(ps.items[i].id);
                }
                ImGui::EndCombo();
            }
            if (l.frames != LayerFrames_Current && run_num_frames(app) < 2) {
                prop::hint(l.style == LayerStyle_Points ? "No trajectory: the frame there is" : "No trajectory to find a density over");
            }
        }

        prop::row("style");
        {
            int s = (int)l.style;
            static const int values[] = { LayerStyle_Points, LayerStyle_Contours, LayerStyle_Heatmap };
            density_map_segmented("##style", style_names, values, LayerStyle_Count, &s, nullptr, true);
            l.style = (LayerStyle)s;
        }
        if (l.style != LayerStyle_Heatmap && !(l.style == LayerStyle_Points && l.point_coloring != PointColoring_Layer)) {
            prop::row("color");
            prop::color("##color", &l.color.x, true);
        }
        prop::row("opacity");
        ImGui::SliderFloat("##opacity", &l.opacity, 0.0f, 1.0f, "%.2f");

        if (l.style == LayerStyle_Contours) {
            prop::row("draw");
            ImGui::Checkbox("fill", &l.fill);
            ImGui::SameLine();
            ImGui::Checkbox("lines", &l.lines);
            if (l.lines) {
                prop::row("line width");
                ImGui::SliderFloat("##line_width", &l.line_width, 0.5f, 6.0f, "%.1f px");
            }
            prop::row("levels", "The contours enclose these fractions of the samples (highest density regions)");
            density_map_levels_editor(l.levels, &l.num_levels);
        }
        if (l.style == LayerStyle_Heatmap || (l.style == LayerStyle_Points && l.point_coloring == PointColoring_Time)) {
            prop::row("colormap");
            ImPlot::ColormapSelection("##cmap", &l.colormap);
        }
        if (l.style == LayerStyle_Heatmap) {
            prop::row("log scale", "Four decades below the maximum");
            ImGui::Checkbox("##log", &l.log_scale);
        }
        if (l.style != LayerStyle_Points) {
            prop::row("bandwidth", "The standard deviation of the gaussian the density is smoothed with, of the triangle's side. "
                                   "Auto: from the spread and the number of the samples (Scott's rule).");
            bool automatic = l.bandwidth <= 0.0f;
            if (ImGui::Checkbox("auto", &automatic)) {
                l.bandwidth = automatic ? 0.0f : (l.auto_bandwidth > 0.0f ? l.auto_bandwidth : 0.02f);
            }
            ImGui::SameLine();
            float pct = (automatic ? l.auto_bandwidth : l.bandwidth) * 100.0f;
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (automatic) ImGui::BeginDisabled();
            if (ImGui::SliderFloat("##bandwidth", &pct, 0.5f, MAX_BANDWIDTH * 100.0f, "%.1f%%", ImGuiSliderFlags_AlwaysClamp)) l.bandwidth = pct / 100.0f;
            if (automatic) ImGui::EndDisabled();
        }
        if (l.style == LayerStyle_Points) {
            prop::row("point size");
            ImGui::SliderFloat("##point_size", &l.point_size, 1.0f, 12.0f, "%.0f px");
            prop::row("color by");
            int c = (int)l.point_coloring;
            if (ImGui::Combo("##coloring", &c, point_coloring_names, PointColoring_Count)) l.point_coloring = (PointColoring)c;
            if (l.frames != LayerFrames_Current && l.point_stride > 1) {
                prop::hint("Every %u frames, to keep to %u points", l.point_stride, (uint32_t)MAX_POINTS);
            }
        }

        ImGui::PopID();
        prop::end();
    }

    void draw_sidebar() {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("LAYERS");
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - ImGui::GetFrameHeight()));
        if (ImGui::Button(ICON_FA_PLUS "##add", ImVec2(ImGui::GetFrameHeight(), 0))) ImGui::OpenPopup("##add_menu");
        ImGui::SetItemTooltip("Add a layer, in the group of what is selected, or a group");
        if (ImGui::BeginPopup("##add_menu")) {
            draw_add_menu();
            ImGui::EndPopup();
        }

        // The list: what is at the top level, and in the groups open what is in them
        SidebarActions act;
        {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t n = tree_nodes(nodes);
            int top[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t nt = representation_tree_children(top, ARRAY_SIZE(top), nodes, n, 0);
            for (size_t k = 0; k < nt; ++k) {
                if (top[k] < num_layers) draw_layer_row(layers[top[k]], 0.0f, nodes, n, &act);
                else                     draw_group(groups[top[k] - num_layers], nodes, n, &act);
            }
        }
        if (num_layers == 0 && num_groups == 0) ImGui::TextDisabled("No layers: add one with +");
        // Below the list: the end of it, out of any group
        {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t n = tree_nodes(nodes);
            const ImVec2 p = ImGui::GetCursorScreenPos();
            const ImRect r(p, ImVec2(p.x + ImGui::GetContentRegionAvail().x, p.y + ImGui::GetFrameHeight()));
            ImGui::Dummy(ImVec2(0, ImGui::GetFrameHeight() * 0.25f));
            LayerNodeRef src;
            if (layer_tree_drop_target_end(r, r.Min.x, LIST_ID, 0, nodes, n, &src)) {
                act.move = true;
                act.move_src = src;
                act.move_target = {};
                act.move_drop = RepresentationDrop_IntoLast;
            }
        }

        if (act.move) move_node(act.move_src, act.move_target, act.move_drop);
        if (act.duplicate) duplicate_layer(act.duplicate);
        if (act.group) group_layer(act.group);
        if (act.to_group_layer) move_layer_to_group(act.to_group_layer, act.to_group);
        if (act.period) add_period_layer(act.period, act.period_group);
        if (act.group_action != LayerGroupAction_None) group_action(act.group_id, act.group_action);
        if (act.remove) remove_layer(act.remove);

        if (LayerGroup* g = find_group(selected_group)) {
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();
            prop::begin("shapespace");
            layer_group_inspector(*g, app, "structures");
            prop::end();
        } else if (Layer* l = find_layer(selected_layer)) {
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();
            draw_inspector(*l);
        }

        // Delete takes the selected layer away, or the selected group with what is in it
        if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Delete)) {
            if (selected_group) group_action(selected_group, LayerGroupAction_Remove);
            else if (selected_layer) remove_layer(selected_layer);
        }
        if (tree_dirty) apply_tree();
    }

    void draw() {
        highlights_now = false;
        plot_drawn = false;
        draw_window();
        // Collapsed, hidden or closed in the middle of a box: no selection from it once it is back
        if (!plot_drawn && box.active) {
            box.active = false;
            box.click_layer = -1;
            box.click_structure = -1;
        }
        if (app && owns_highlight && !highlights_now) md_bitmask_clear(&app->selection.highlight_mask);
        owns_highlight = highlights_now;
    }

    void draw_window() {
        if (!show_window || !app) return;
        ImGui::SetNextWindowSize(ImVec2(720, 480), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Shapespace", &show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar)) {
            draw_toolbar();
            if (show_sidebar) {
                ImGui::BeginChild("##sidebar", ImVec2(320, 0), ImGuiChildFlags_ResizeX | ImGuiChildFlags_Borders);
                draw_sidebar();
                ImGui::EndChild();
                ImGui::SameLine();
            }
            ImGui::BeginChild("##plot", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
            draw_plot();
            ImGui::EndChild();
        }
        ImGui::End();
    }
};

static Shapespace instance = {};

}  // namespace shapespace
