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
//
// A selection is of points, and a point is a structure in a frame: what is selected is those structures (the atom
// selection) and those frames (the frame selection, ApplicationState::selection.frame_mask, which the timeline shows),
// or either alone (the settings). A point drawn in place of several frames (where frames are skipped to keep to a number
// of points) stands for them all. The points of the frame shown say nothing of time and select no frames. A layer may be
// over the frames selected.
//
// In 3D the triangle is stood up across the time axis (the timeline's time): each structure of a layer over frames is a
// line through its shapes, smoothed over a window of frames or not, the frame shown a slice across it, the timeline's
// periods bands along the side of its time axis. As in every 3D view (plot_input.h), left-drag turns it. The triangle stays as it is
// (equilateral, all of it, drawn at either end of the time shown); what moves is time: the wheel and middle-drag zoom
// it, right-drag pans it, double-click shows all of it again. The time shown is the timeline's visible range (linked,
// both ways), a window of its own, or a tail following the frame shown, and the paths are cut to it. The box fills the
// plot. Hover, click and shift-drag select as in 2D, what is drawn on the screen.

#include <core/md_log.h>
#include <md_run.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_file.h>

#include <md_csv.h>
#include <md_xvg.h>
#include <md_system.h>
#include <core/md_geom.h>

#include <viamd_event.h>
#include <window_menu.h>
#include <viamd.h>
#include <display_units.h>
#include <frame_axis.h>
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
#include <implot3d.h>
#include <plot_input.h>
#include <task_system.h>

#include <atomic>
#include <memory>
#include <vector>
#include <math.h>
#include <string.h>

namespace shapespace {

// ## Constants

// The corners of the triangle, in the plot's coordinates: an equilateral triangle of side 1
static const md_float2_t tri_lin = { 0.0f, 0.0f };
static const md_float2_t tri_pla = { 1.0f, 0.0f };
static const md_float2_t tri_iso = { 0.5f, 0.86602540378f };

// The plane the densities are over: square, around the triangle, with room for the images of the samples near its
// edges that keep the density from leaking out of it (density_map.h), some 4 sigma
static const DensityMapDomain shape_domain = {
    -0.25f, 1.25f, 0.433f - 0.75f, 0.433f + 0.75f, false,     // 0.25 = 4 * MAX_BANDWIDTH, past the corners
    true, { tri_lin.x, tri_lin.y, tri_pla.x, tri_pla.y, tri_iso.x, tri_iso.y },
};

constexpr int MAX_LAYERS      = 16;
constexpr int MAX_LEVELS      = DENSITY_MAP_MAX_LEVELS;
constexpr size_t MAX_POINTS   = 100000;     // Of a layer of points over frames: frames are skipped beyond
constexpr size_t MAX_LINE_POINTS = 60000;   // Of a layer's lines in 3D, whose segments are sorted by depth every frame
constexpr float PERSPECTIVE_FOV = 0.5236f;   // Radians, 30 degrees: of the 3D view, in perspective
constexpr int MAX_SMOOTHING   = 501;        // Frames, of the window the shapes are smoothed over (Shapespace::smoothing)
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

static inline md_float2_t shape_coord(md_float3_t w) {
    return tri_lin * w.x + tri_pla * w.y + tri_iso * w.z;
}

// The weights back from the point: they sum to one (md_geom_shape_weights), which the corners' barycentric coordinates
// are, so the point is all there is to keep
static inline md_float3_t shape_weights_of(md_float2_t c) {
    const float iso = c.y / tri_iso.y;
    const float pla = c.x - tri_iso.x * iso;
    return { 1.0f - pla - iso, pla, iso };
}

static inline bool coord_valid(md_float2_t c) {
    return c.x == c.x && c.y == c.y;    // Not NaN: the frames not evaluated yet are
}

// The shape weights (linear, planar, isotropic) of the atoms of a structure in a state: the eigenvalues of the gyration
// tensor about their centre, placed as one where the cell is periodic. masses may be NULL. NaN where there is no shape
// (fewer than two atoms, or all in one place).
static md_float3_t structure_shape(md_bitset_t atoms, const md_float3_t* xyz, const md_unitcell_t* cell, const float* masses, md_temp_scope_t temp) {
    const size_t count = md_bitset_count(atoms);
    if (count < 2) return { NAN, NAN, NAN };   // No shape: not a point
    md_float4_t* xyzw = (md_float4_t*)md_temp_alloc_aligned(temp, sizeof(md_float4_t) * count, 16);
    md_bitset_iter_t it = md_bitset_iter(atoms);
    size_t k = 0;
    while (md_bitset_iter_next(&it) && k < count) {
        xyzw[k++] = md_float4_from_float3(xyz[it.idx], masses ? masses[it.idx] : 1.0f);
    }
    // Covariance about the mean of the placed points, not the circular mean the placement is seeded from (their offset
    // enters the covariance squared)
    md_float3_t com = md_float3_zero();
    md_geom_deperiodize_self(xyzw, k, cell, &com);
    const md_float3x3_t M = md_float3x3_covariance_matrix_float4(xyzw, 0, k, com);
    return md_geom_shape_weights(&M);
}

// ## Evaluations

// The shapes of a set of structures over every frame of the run, weighted by mass or not: shared by the layers of the
// same structures and weighting, and freed when none uses it
struct Evaluation {
    uint64_t key = 0;                       // Of the structures and the weighting
    md_allocator_i* arena = nullptr;
    md_array(md_bitset_t) structures = 0;   // In arena
    bool     use_mass = true;
    uint32_t first_frame = 0;               // Of the run: its frames are [first_frame, first_frame + num_frames)
    size_t   num_frames = 0;
    size_t   num_structures = 0;
    md_float2_t*  coords = nullptr;              // [structure * num_frames + frame - first_frame], in the triangle: NaN for frames not done
    std::atomic<uint32_t> frames_done{0};

    // Its shapes smoothed over smoothed_window frames (Shapespace::smoothing), laid out as coords: made on the main
    // thread once it is complete, and held by what reads them elsewhere (a density job) while the window changes
    std::shared_ptr<const std::vector<md_float2_t>> smoothed;
    int smoothed_window = 0;
    task_system::Job task;

    bool complete() const { return num_frames > 0 && frames_done.load() >= num_frames; }
    // Stopped before all of its frames were done (interrupted from the task overlay, or by a load): never will be
    // (Asked in this order: a task that has stopped has counted all the frames it did.)
    bool failed() const { return !task.running() && !complete(); }

    // Waits for its task, which on the main thread may run pool work meanwhile: what is dropped in use goes through
    // Shapespace::retired first, freed once its task has stopped. Before the arena goes, which the task reads.
    ~Evaluation() {
        task.cancel_and_wait();
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

// How a layer over frames is drawn in 3D, over time
enum Trajectory : int {
    Trajectory_Lines = 0,       // A line through the frames of each structure
    Trajectory_Points,
    Trajectory_Both,
    Trajectory_Count,
};

enum ViewMode : int {
    ViewMode_2D = 0,
    ViewMode_3D,
    ViewMode_Split,             // Both, side by side
    ViewMode_Count,
};

// What a selection in the plot selects
enum SelectTarget : int {
    SelectTarget_Both = 0,      // The structures (atoms) and the frames of the points
    SelectTarget_Structures,
    SelectTarget_Frames,
    SelectTarget_Count,
};

// What the time shown in 3D is (Shapespace::time_win)
enum TimeWindow : int {
    TimeWindow_Timeline = 0,    // The timeline's visible range, both ways
    TimeWindow_Own,             // A window of its own
    TimeWindow_Follow,          // Following the frame shown
    TimeWindow_Count,
};

static const char* style_names[LayerStyle_Count] = { "Points", "Contours", "Heatmap" };
static const char* trajectory_names[Trajectory_Count] = { "Lines", "Points", "Both" };
static const char* view_mode_names[ViewMode_Count] = { "2D", "3D", "Split" };
static const char* select_target_names[SelectTarget_Count] = { "Structures and frames", "Structures", "Frames" };
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

    ImPlotColormap colormap = ImPlotColormap_Spectral;     // Of time: bright against the plots' dark background
    ImPlotColormap density_colormap = ImPlotColormap_Viridis;  // Of a heatmap: sequential, as a density is
    ImPlotColormap structure_colormap = ImPlotColormap_Deep;    // Of points coloured by structure: in turn, by default
    bool     log_scale = false;

    PointColoring point_coloring = PointColoring_Layer;
    float    point_size = 4.0f;

    // In 3D: points over frames as lines, as points or both
    Trajectory trajectory = Trajectory_Lines;
    float    trajectory_width = 1.0f;

    // ### Derived, not saved

    md_array(md_bitset_t) structures = 0;   // Heap allocated, as the filter gave them
    uint64_t structs_key = 0;       // Of what they were found from
    uint64_t structs_hash = 0;      // Of the structures themselves
    bool     structs_valid = false;

    EvaluationRef eval;             // Of frames: over the run

    // Of the frame shown
    md_array(md_float2_t) current_coords = 0;
    uint64_t current_key = 0;

    // Of the frame shown, smoothed, when its structures follow it: theirs over the frames of the window about it
    EvaluationRef window;
    md_array(md_float2_t) window_coords = 0;     // Smoothed shapes of window_structs, the last that were done
    md_array(md_bitset_t) window_structs = 0;    // Heap allocated
    uint64_t window_key = 0;                     // Of what window_coords are of

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
    uint32_t line_stride = 1;       // Of the lines in 3D

};

// Scott's rule for two dimensions, sigma n^(-1/6), of the spread of the samples (both axes alike, as the plane is
// isotropic): what smooths a single gaussian best, and much of anything else well enough. From at most some 100k of them.
static float scott_bandwidth(const Evaluation& e, const md_float2_t* coords, md_bitset_t frames) {
    const uint64_t total = (uint64_t)e.num_structures * md_bitset_count(frames);
    const uint32_t step = (uint32_t)MD_MAX((uint64_t)1, total / 100000);
    double n = 0, sx = 0, sy = 0, sxx = 0, syy = 0;
    for (size_t s = 0; s < e.num_structures; ++s) {
        const md_float2_t* c = coords + s * e.num_frames;
        md_bitset_iter_t it = md_bitset_iter_range(frames, 0, (uint32_t)e.num_frames);
        uint32_t k = 0;
        while (md_bitset_iter_next(&it)) {
            const uint32_t f = it.idx;
            if (k++ % step) continue;
            if (!coord_valid(c[f])) continue;
            n += 1;
            sx += c[f].x;
            sy += c[f].y;
            sxx += (double)c[f].x * c[f].x;
            syy += (double)c[f].y * c[f].y;
        }
    }
    if (n < 2) return 0.02f;
    const double var = MD_MAX(0.0, (sxx - sx * sx / n) + (syy - sy * sy / n)) / (2.0 * (n - 1));
    // n: of all the samples, not of those the spread was taken from
    const double h = sqrt(var) * pow((double)total, -1.0 / 6.0);
    return MD_CLAMP((float)h, 0.005f, MAX_BANDWIDTH);
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
        uint32_t click_frame = 0;       // and the frames it stands for, [click_frame, click_next)
        uint32_t click_next = 0;
    } box;

    ViewMode     view_mode = ViewMode_2D;

    // The shapes over frames are smoothed over a window of this many frames (a moving average; 1 is not smoothed), of
    // every layer and view alike: the points and densities in 2D, the paths in 3D, the points of the frame shown (as the
    // smoothed shapes there) and what is selected. Once a layer's evaluation is complete; not of what is exported.
    int smoothing = 1;
    bool perspective = false;       // Of the 3D view: seen from a point (PERSPECTIVE_FOV), not from afar

    // The window of time the 3D plot shows and cuts the paths to (### The time window of the 3D plot), in the run's own
    // unit: of its own, the times [beg, end] (all of the run when end <= beg); following the frame shown, offsets from
    // its time ([-length, 0] is a tail behind it); linked, the timeline's visible range. Not of the 2D plot.
    struct {
        TimeWindow mode = TimeWindow_Timeline;
        double beg = 0.0;
        double end = 0.0;
        bool   show_periods = true;     // The timeline's periods along time
    } time_win;
    SelectTarget select_target = SelectTarget_Both;

    // The frame highlight (selection.frame_highlight_mask) is ours to clear, as the atoms' highlight
    bool owns_frame_highlight = false;
    bool frame_highlights_now = false;

    // The 3D view
    struct {
        float  view_scale = 1.0f;       // Of the box on the screen (ImPlot3D's ViewScaleFactor): to fill the plot
        bool   panning = false;         // Right-drag, along time
        bool   dollying = false;        // Middle-drag: zooming in time, about the time where it began (dolly_u of the axis)
        bool   scrubbing = false;       // Ctrl + left-drag: the frame shown is the one at the pointer's time
        bool   over_box = false;        // The pointer, last frame: over one of the box's planes
        float  dolly_u = 0.5f;
        int    preset = -1;             // A view to turn to: 0 along the time axis, 1 from the side, 2 the default
        bool   drawn = false;           // This frame
        // Box selection, on the screen
        bool   active = false;
        SelectionOperator op = SelectionOperator::Or;
        ImVec2 p0, p1;
        int    click_layer = -1;
        int64_t click_structure = -1;
        uint32_t click_frame = 0, click_next = 0;
        // A press of the left button without a drag, on a point of a frame: goes there when let go
        bool   pressed = false;
        ImVec2 press_pos;
    } v3;

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
            case viamd::EventType_ViamdWindowMenuRequest:
                viamd::window_menu_add(e, "Analysis/Shapespace", &show_window, {.order = 50});
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
                view_mode = ViewMode_2D;
                select_target = SelectTarget_Both;
                time_win = {};
                smoothing = 1;
                perspective = false;
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                const md_str_t section = viamd::section_header(state);
                if (md_str_eq(section, MD_STR_LIT("ShapeSpace"))) {
                    deserialize_settings(state);
                } else if (md_str_eq(section, MD_STR_LIT("ShapeSpaceLayer"))) {
                    deserialize_layer(state);
                } else if (md_str_eq(section, MD_STR_LIT("ShapeSpaceGroup"))) {
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
            window_clear(l);
            l.has_result = false;
            l.shown_hash = 0;
            l.want_hash = 0;
            l.structs_key = 0;
            l.current_key = 0;
        }
        for (EvaluationRef& e : evaluations) e->task.cancel();
        for (EvaluationRef& e : retired) e->task.cancel();
        for (EvaluationRef& e : evaluations) e->task.wait();
        for (EvaluationRef& e : retired) e->task.wait();
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
        for (int p = 0; p < (int)MD_COUNTOF(layer_palette); ++p) {
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
        l.window.reset();
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
        window_clear(l);
        md_array_free(l.window_coords, heap);
        md_array_free(l.window_structs, heap);
        l.eval.reset();
        l.window.reset();
    }

    // Layers are moved as bytes but for their evaluation, which is counted: moved from one to the other
    static void move_into(Layer& dst, Layer& src) {
        EvaluationRef e = std::move(src.eval), w = std::move(src.window);
        dst.eval.reset();
        dst.window.reset();
        md_memcpy((void*)&dst, (const void*)&src, sizeof(Layer));
        new (&dst.eval) EvaluationRef(std::move(e));
        new (&dst.window) EvaluationRef(std::move(w));
        new (&src.eval) EvaluationRef();
        new (&src.window) EvaluationRef();
    }

    void remove_layer(uint32_t id) {
        const int idx = layer_index(id);
        if (idx < 0) return;
        free_layer(layers[idx]);
        for (int i = idx; i + 1 < num_layers; ++i) move_into(layers[i], layers[i + 1]);
        num_layers -= 1;
        layers[num_layers].eval.reset();
        layers[num_layers].window.reset();
        new (&layers[num_layers]) Layer();
        if (selected_layer == id) selected_layer = num_layers > 0 ? layers[MD_MIN(idx, num_layers - 1)].id : 0;
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
        MD_STATIC_ASSERT(sizeof(l->name) == sizeof(src.name), "");
        md_memcpy(l->name, src.name, sizeof(l->name));
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
        md_memcpy(l->levels, src.levels, sizeof(l->levels));
        l->bandwidth = src.bandwidth;
        l->colormap = src.colormap;
        l->density_colormap = src.density_colormap;
        l->structure_colormap = src.structure_colormap;
        l->log_scale = src.log_scale;
        l->point_coloring = src.point_coloring;
        l->point_size = src.point_size;
        l->trajectory = src.trajectory;
        l->trajectory_width = src.trajectory_width;
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
        md_memcpy(g->filter.text, l->filter.text, sizeof(g->filter.text));
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
            md_memcpy(name, g->name, sizeof(name));
            md_memcpy(filter, g->filter.text, sizeof(filter));
            const bool visible = g->visible, open = g->open, follow = g->follow_selection;
            LayerGroup* copy = add_group(g->node.order + 1);
            if (!copy) return;
            md_memcpy(copy->name, name, sizeof(name));
            md_memcpy(copy->filter.text, filter, sizeof(filter));
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
                const size_t nk = representation_tree_children(kids, MD_COUNTOF(kids), nodes, nn, id);
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
            // Everything, as one structure, in the language of the workspace
            const bool v2 = app && app->script.language == script::Language::V2;
            snprintf(g->filter.text, sizeof(g->filter.text), "%s", v2 ? "particle" : "all");
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

    // The density of the shapes of the layer's evaluation over a set of its frames
    void launch_job(Layer& l, uint64_t hash, md_bitset_t frames, float sigma) {
        MD_ASSERT(!l.job && l.eval && l.eval->complete());
        EvaluationRef eval = l.eval;
        // The frames, as the job's own: the set may be a view of the selection, which changes
        auto list = std::make_shared<std::vector<uint32_t>>();
        list->reserve((size_t)md_bitset_count(frames));
        md_bitset_iter_t it = md_bitset_iter_range(frames, 0, (uint32_t)eval->num_frames);
        while (md_bitset_iter_next(&it)) list->push_back(it.idx);
        const uint32_t nf = (uint32_t)list->size();
        if (nf == 0) return;
        const uint64_t num_samples = (uint64_t)eval->num_structures * nf;
        const uint32_t num_chunks = (uint32_t)MD_MIN((uint64_t)MD_MIN(8u, nf), (num_samples >> 20) + 1);
        // The shapes as shown: the smoothed ones held by the job, as the window may change while it runs
        std::shared_ptr<const std::vector<md_float2_t>> smoothed = eval_shapes(*eval) == eval->coords ? nullptr : eval->smoothed;
        auto source = [eval, list, smoothed](uint32_t chunk, uint32_t chunks, DensityMapBins& bins, const std::atomic<int32_t>& cancel) {
            const uint32_t n = (uint32_t)list->size();
            const uint32_t k0 = (uint32_t)((uint64_t)n * chunk / chunks);
            const uint32_t k1 = (uint32_t)((uint64_t)n * (chunk + 1) / chunks);
            const md_float2_t* coords = smoothed ? smoothed->data() : eval->coords;
            for (size_t s = 0; s < eval->num_structures; ++s) {
                if (cancel.load(std::memory_order_relaxed)) return;
                const md_float2_t* c = coords + s * eval->num_frames;
                for (uint32_t k = k0; k < k1; ++k) {
                    const uint32_t f = (*list)[k];
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
        l.job = density_map_launch(MD_STR_LIT("Shapespace density"), shape_domain, sigma, num_chunks, source, done);
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
            l.structs_valid = filter_evaluate_sets(&l.filter, &l.structures, app->script.language, layer_group_filter_context(app, in_context ? g : nullptr), heap) &&
                              md_array_size(l.structures) > 0;
        }
        uint64_t h = 0x5A9E;
        for (size_t i = 0; i < md_array_size(l.structures); ++i) h = md_bitset_hash64(l.structures[i], h + i);
        l.structs_hash = h;
    }

    // Whether the layer's structures follow the frame shown (what picks them depends on the coordinates): then they are
    // not the same over frames, and their smoothed shapes are of a window about the frame shown (update_window)
    bool structures_follow_frame(const Layer& l) const {
        if (l.filter.text[0] && l.filter.result.is_dynamic) return true;
        const LayerGroup* g = layer_group(l);
        return l.filter.text[0] == '\0' && g && !g->follow_selection && g->filter.result.is_dynamic;
    }

    // Whether the points of the frame shown are the smoothed shapes there, of an evaluation over the run, rather than the
    // shapes of the state shown
    bool current_smoothed(const Layer& l) const {
        return l.frames == LayerFrames_Current && smoothing > 1 && run_num_frames(app) >= 2 && !structures_follow_frame(l);
    }

    // The smoothed shape at frame f of the run (of num_frames) from the shapes c of frames [first, ...) about it: the
    // mean of those valid in a centred window of `smoothing` frames that shrinks at the ends, so that it stays centred.
    // As VIAL's smooth(x, n), so that smooth(shape_weights(s), n) gives the same shapes. NaN when none is valid.
    md_float2_t smoothed_at(const md_float2_t* c, uint32_t first, size_t num_frames, size_t f) const {
        const size_t w = MD_MIN((size_t)(smoothing / 2), MD_MIN(f, num_frames - 1 - f));
        double sx = 0, sy = 0;
        uint32_t n = 0;
        for (size_t g = f - w; g <= f + w; ++g) {
            const md_float2_t p = c[g - first];
            if (!coord_valid(p)) continue;
            sx += p.x;
            sy += p.y;
            n += 1;
        }
        return n ? md_float2_t{ (float)(sx / n), (float)(sy / n) } : md_float2_t{ NAN, NAN };
    }

    // The evaluation's shapes smoothed (smoothed_at), made once it is complete and again when the window changes
    void update_smoothed(Evaluation& e) {
        if (smoothing <= 1 || !e.complete()) return;
        if (e.smoothed && e.smoothed_window == smoothing) return;
        const size_t nf = e.num_frames;
        const size_t half = (size_t)(smoothing / 2);
        auto out = std::make_shared<std::vector<md_float2_t>>(e.num_structures * nf);
        for (size_t s = 0; s < e.num_structures; ++s) {
            const md_float2_t* c = e.coords + s * nf;
            md_float2_t* o = out->data() + s * nf;
            // A running sum over the window [f - w, f + w]: both ends only move forward, by two where it grows or shrinks
            double sx = 0, sy = 0;
            int64_t n = 0, lo = 0, hi = -1;
            for (size_t f = 0; f < nf; ++f) {
                const size_t w = MD_MIN(half, MD_MIN(f, nf - 1 - f));
                const int64_t want_lo = (int64_t)(f - w), want_hi = (int64_t)(f + w);
                while (hi < want_hi) {
                    ++hi;
                    if (coord_valid(c[hi])) { sx += c[hi].x; sy += c[hi].y; ++n; }
                }
                while (lo < want_lo) {
                    if (coord_valid(c[lo])) { sx -= c[lo].x; sy -= c[lo].y; --n; }
                    ++lo;
                }
                o[f] = n > 0 ? md_float2_t{ (float)(sx / n), (float)(sy / n) } : md_float2_t{ NAN, NAN };
            }
        }
        e.smoothed = std::move(out);
        e.smoothed_window = smoothing;
    }

    // The shapes of an evaluation as they are shown: smoothed, when they are and the smoothing is ready; else its own.
    // Of the window last smoothed over until this one's are made (update_smoothed): not its own shapes meanwhile, which
    // would flicker in and out of the smoothed ones as the window is changed.
    const md_float2_t* eval_shapes(const Evaluation& e) const {
        if (smoothing > 1 && e.smoothed) return e.smoothed->data();
        return e.coords;
    }

    // The evaluation of the layer's structures and weighting over the run: one shared, or a new one started. Of a layer
    // of the frame shown too, when its points there are the smoothed shapes.
    void update_evaluation(Layer& l) {
        const size_t num_frames = run_num_frames(app);
        if ((l.frames == LayerFrames_Current && !current_smoothed(l)) || !l.structs_valid || num_frames == 0) {
            l.eval.reset();
            return;
        }
        const uint64_t key = md_hash64_combine(md_hash64_combine(l.structs_hash, l.use_mass ? 0x11 : 0x22), md_hash64_combine(system_generation, num_frames));
        if (l.eval && l.eval->key == key && !l.eval->failed()) return;
        l.eval = evaluation(key, l.structures, l.use_mass, 0, num_frames);
    }

    // The evaluation of the key: one shared, or a new one started of the structures and weighting over the frames
    // [first, first + count) of the run
    EvaluationRef evaluation(uint64_t key, const md_array(md_bitset_t) structures, bool use_mass, uint32_t first, size_t count) {
        for (size_t i = 0; i < evaluations.size(); ++i) {
            if (evaluations[i]->key != key) continue;
            if (evaluations[i]->failed()) {
                // Started again, for every layer over it
                retire(i);
                break;
            }
            return evaluations[i];
        }
        EvaluationRef e = std::make_shared<Evaluation>();
        e->key = key;
        e->use_mass = use_mass;
        e->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
        for (size_t i = 0; i < md_array_size(structures); ++i) {
            md_array_push(e->structures, md_bitset_copy(structures[i], e->arena), e->arena);
        }
        e->num_structures = md_array_size(e->structures);
        e->first_frame = first;
        e->num_frames = count;
        const size_t n = e->num_structures * count;
        e->coords = (md_float2_t*)md_alloc(e->arena, sizeof(md_float2_t) * MD_MAX(n, (size_t)1));
        for (size_t i = 0; i < n; ++i) e->coords[i] = { NAN, NAN };

        Evaluation* ev = e.get();
        ApplicationState* state = app;
        e->task.start_for(MD_STR_LIT("Shapespace"), (uint32_t)count, 1, [ev, state](uint32_t beg, uint32_t end) {
            const md_system_t& sys = state->mold.sys;
            md_temp_scope_t temp = md_temp_begin();
            MD_DEFER { md_temp_end(temp); };
            md_float3_t* xyz = (md_float3_t*)md_temp_alloc(temp, sizeof(md_float3_t) * MD_ALIGN_TO(sys.atom.count, 8));
            float* masses = nullptr;
            if (ev->use_mass) {
                masses = (float*)md_temp_alloc(temp, sizeof(float) * MD_ALIGN_TO(sys.atom.count, 8));
                md_atom_extract_masses(masses, 0, sys.atom.count, &sys.atom);
            }
            // One context for the range, so the run's files stay open across its frames
            const md_str_t paths[] = { MD_STR_INIT("atom/position"), MD_STR_INIT("unitcell") };
            md_system_extract_t* ex = md_system_extract_begin(&sys, md_str_from_cstr(state->mold.run), paths, MD_COUNTOF(paths), md_get_heap_allocator());
            if (!ex) {
                // Nothing to read them from: done, without shapes
                ev->frames_done.fetch_add(end - beg);
                return;
            }
            MD_DEFER { md_system_extract_end(ex); };
            for (uint32_t f = beg; f < end; ++f) {
                if (task_system::cancelled()) return;
                md_system_state_t frame = {};
                frame.num_atoms = sys.atom.count;
                frame.xyz = xyz;
                if (!md_system_extract_frame(ex, ev->first_frame + f, &frame)) {
                    ev->frames_done.fetch_add(1);
                    continue;
                }
                for (size_t s = 0; s < ev->num_structures; ++s) {
                    md_temp_scope_t inner = md_temp_begin();
                    const md_float3_t w = structure_shape(ev->structures[s], xyz, &frame.unitcell, masses, inner);
                    md_temp_end(inner);
                    // Degenerate (one atom, or all in one place): no shape, NaN
                    ev->coords[s * ev->num_frames + f] = shape_coord(w);
                }
                ev->frames_done.fetch_add(1);
            }
        });
        evaluations.push_back(e);
        return e;
    }

    // Whether the points of the frame shown are smoothed shapes of structures that follow the frame: of those picked in
    // the state shown, over the frames of the window about it (they are not the same structures in other frames)
    bool window_smoothed(const Layer& l) const {
        return l.frames == LayerFrames_Current && smoothing > 1 && run_num_frames(app) >= 2 && structures_follow_frame(l);
    }

    void window_clear(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        l.window.reset();
        l.window_key = 0;
        md_array_shrink(l.window_coords, 0);
        for (size_t i = 0; i < md_array_size(l.window_structs); ++i) md_bitset_free(&l.window_structs[i], heap);
        md_array_shrink(l.window_structs, 0);
    }

    // The smoothed shapes of the frame shown, of structures that follow it: an evaluation of them over the frames of the
    // windows about the frames on either side of it, whose means there are blended as the paths are drawn. Until it is
    // done, those of the frame before stay shown (with their structures): a few frames late while playing, not raw.
    void update_window(Layer& l) {
        if (!window_smoothed(l) || !l.structs_valid) {
            window_clear(l);
            return;
        }
        const size_t nf = run_num_frames(app);
        const double fr = MD_CLAMP(app->animation.frame, 0.0, (double)(nf - 1));
        const size_t f0 = (size_t)fr, f1 = MD_MIN(f0 + 1, nf - 1);
        const size_t half = (size_t)(smoothing / 2);
        const size_t lo = f0 - MD_MIN(half, MD_MIN(f0, nf - 1 - f0));
        const size_t hi = f1 + MD_MIN(half, MD_MIN(f1, nf - 1 - f1));     // Of the windows of f0 and f1, the frames
        uint64_t key = md_hash64_combine(md_hash64_combine(l.structs_hash, l.use_mass ? 0x11 : 0x22), md_hash64_combine(system_generation, nf));
        key = md_hash64_combine(key, md_hash64_combine(0x3170 + lo, hi));
        if (!l.window || l.window->key != key || l.window->failed()) l.window = evaluation(key, l.structures, l.use_mass, (uint32_t)lo, hi - lo + 1);
        if (!l.window->complete()) return;

        const uint64_t shown = md_hash64(&fr, sizeof(fr), md_hash64_combine(key, (uint64_t)smoothing));
        if (shown == l.window_key) return;
        l.window_key = shown;
        md_allocator_i* heap = md_get_heap_allocator();
        const Evaluation& e = *l.window;
        const float a = (float)(fr - (double)f0);
        md_array_resize(l.window_coords, e.num_structures, heap);
        for (size_t s = 0; s < e.num_structures; ++s) {
            const md_float2_t* c = e.coords + s * e.num_frames;
            const md_float2_t p0 = smoothed_at(c, e.first_frame, nf, f0), p1 = smoothed_at(c, e.first_frame, nf, f1);
            l.window_coords[s] = !coord_valid(p1) ? p0 : !coord_valid(p0) ? p1 : md_float2_t{ p0.x + (p1.x - p0.x) * a, p0.y + (p1.y - p0.y) * a };
        }
        for (size_t i = 0; i < md_array_size(l.window_structs); ++i) md_bitset_free(&l.window_structs[i], heap);
        md_array_shrink(l.window_structs, 0);
        for (size_t i = 0; i < e.num_structures; ++i) md_array_push(l.window_structs, md_bitset_copy(e.structures[i], heap), heap);
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
            masses = (float*)md_temp_alloc(temp, sizeof(float) * MD_ALIGN_TO(app->mold.sys.atom.count, 8));
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
            const int f = MD_CLAMP((int)(app->animation.frame + 0.5), 0, num_frames - 1);
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

    // The frames of the layer as a set: a run of them. False when it covers none.
    bool layer_frames(const Layer& l, md_bitset_t* out) const {
        uint32_t beg = 0, end = 0;
        if (!layer_frame_range(l, &beg, &end)) return false;
        *out = md_bitset_range(beg, end);
        return true;
    }

    // What the frames of the layer are, as far as what is found over them goes
    uint64_t layer_frames_hash(const Layer& l, md_bitset_t frames) const {
        const uint32_t span[2] = { md_bitset_first(frames), md_bitset_last_end(frames) };
        return md_hash64(span, sizeof(span), 0);
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
            if (l.eval) update_smoothed(*l.eval);
            update_window(l);

            if (!layer_draws_density(l) || !l.eval || !l.eval->complete()) continue;
            md_bitset_t frames = {};
            if (!layer_frames(l, &frames)) continue;
            const int window = eval_shapes(*l.eval) == l.eval->coords ? 1 : smoothing;
            uint64_t h = md_hash64_combine(md_hash64_combine(l.eval->key, layer_frames_hash(l, frames)), (uint64_t)window);
            // Found once for the samples, not at every frame
            if (h != l.auto_key) {
                l.auto_key = h;
                l.auto_bandwidth = scott_bandwidth(*l.eval, eval_shapes(*l.eval), frames);
            }
            const float sigma = l.bandwidth > 0.0f ? l.bandwidth : l.auto_bandwidth;
            h = md_hash64(&sigma, sizeof(sigma), h);
            h |= 1;
            if (h == l.shown_hash || h == l.want_hash || l.job) continue;
            launch_job(l, h, frames, sigma);
        }
        // Evaluations no layer holds any longer (nor a density job: they hold them too)
        for (size_t i = 0; i < evaluations.size();) {
            if (evaluations[i].use_count() == 1) retire(i);
            else ++i;
        }
        for (size_t i = 0; i < retired.size();) {
            if (retired[i].use_count() == 1 && !retired[i]->task.running()) {
                retired[i] = std::move(retired.back());
                retired.pop_back();
            } else {
                ++i;
            }
        }
    }

    // Out of evaluations[i]: cancelled, and freed once its task has stopped
    void retire(size_t i) {
        evaluations[i]->task.cancel();
        retired.push_back(std::move(evaluations[i]));
        evaluations[i] = std::move(evaluations.back());
        evaluations.pop_back();
    }

    // ## Workspace

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, MD_STR_LIT("ShapeSpace"));
        viamd::write_bool(state, MD_STR_LIT("ShowLayers"), show_sidebar);
        viamd::write_bool(state, MD_STR_LIT("ShowGrid"), show_grid);
        viamd::write_str(state, MD_STR_LIT("View"), md_str_from_cstr(view_mode_names[view_mode]));
        static const md_str_t targets[SelectTarget_Count] = { MD_STR_LIT("both"), MD_STR_LIT("structures"), MD_STR_LIT("frames") };
        viamd::write_str(state, MD_STR_LIT("Selects"), targets[select_target]);
        static const md_str_t windows[TimeWindow_Count] = { MD_STR_LIT("timeline"), MD_STR_LIT("own"), MD_STR_LIT("follow") };
        viamd::write_str(state, MD_STR_LIT("TimeWindow"), windows[time_win.mode]);
        const double tw[2] = { time_win.beg, time_win.end };
        viamd::write_dbl_vec(state, MD_STR_LIT("TimeWindowRange"), tw, 2);
        viamd::write_bool(state, MD_STR_LIT("ShowPeriods"), time_win.show_periods);
        viamd::write_int(state, MD_STR_LIT("Smoothing"), smoothing);
        viamd::write_bool(state, MD_STR_LIT("Perspective"), perspective);
        viamd::write_int(state, MD_STR_LIT("Layers"), num_layers);
        // In the order of the list: a group, then what is in it
        RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t n = tree_nodes(nodes);
        int top[MAX_LAYERS + LAYER_MAX_GROUPS];
        const size_t nt = representation_tree_children(top, MD_COUNTOF(top), nodes, n, 0);
        for (size_t k = 0; k < nt; ++k) {
            if (top[k] < num_layers) {
                serialize_layer(state, layers[top[k]]);
                continue;
            }
            const LayerGroup& g = groups[top[k] - num_layers];
            layer_group_serialize(state, MD_STR_LIT("ShapeSpaceGroup"), g);
            int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t nk = representation_tree_children(kids, MD_COUNTOF(kids), nodes, n, g.node.id);
            for (size_t c = 0; c < nk; ++c) {
                if (kids[c] < num_layers) serialize_layer(state, layers[kids[c]]);
            }
        }
    }

    void serialize_layer(viamd::serialization_state_t& state, const Layer& l) {
        {
            viamd::write_section_header(state, MD_STR_LIT("ShapeSpaceLayer"));
            if (l.node.parent) viamd::write_int(state, MD_STR_LIT("Group"), (int64_t)l.node.parent);
            if (l.name[0]) viamd::write_str(state, MD_STR_LIT("Name"), md_str_from_cstr(l.name));
            viamd::write_bool(state, MD_STR_LIT("Visible"), l.visible);
            if (l.filter.text[0]) viamd::write_str(state, MD_STR_LIT("Filter"), md_str_from_cstr(l.filter.text));
            viamd::write_bool(state, MD_STR_LIT("UseMass"), l.use_mass);
            static const md_str_t frames[LayerFrames_Count] = { MD_STR_LIT("current"), MD_STR_LIT("all"), MD_STR_LIT("period") };
            viamd::write_str(state, MD_STR_LIT("Frames"), frames[l.frames]);
            if (l.frames == LayerFrames_Period) viamd::write_int(state, MD_STR_LIT("Period"), (int)l.period);
            static const md_str_t styles[LayerStyle_Count] = { MD_STR_LIT("points"), MD_STR_LIT("contours"), MD_STR_LIT("heatmap") };
            viamd::write_str(state, MD_STR_LIT("Style"), styles[l.style]);
            viamd::write_vec4(state, MD_STR_LIT("Color"), md_float4_set(l.color.x, l.color.y, l.color.z, l.color.w));
            viamd::write_flt(state, MD_STR_LIT("Opacity"), l.opacity);
            viamd::write_bool(state, MD_STR_LIT("Fill"), l.fill);
            viamd::write_bool(state, MD_STR_LIT("Lines"), l.lines);
            viamd::write_flt(state, MD_STR_LIT("LineWidth"), l.line_width);
            viamd::write_flt_vec(state, MD_STR_LIT("Levels"), l.levels, (size_t)l.num_levels);
            if (l.bandwidth > 0.0f) viamd::write_flt(state, MD_STR_LIT("Bandwidth"), l.bandwidth);
            else                    viamd::write_str(state, MD_STR_LIT("Bandwidth"), MD_STR_LIT("auto"));
            viamd::write_str(state, MD_STR_LIT("Colormap"), md_str_from_cstr(ImPlot::GetColormapName(l.colormap)));
            viamd::write_str(state, MD_STR_LIT("DensityColormap"), md_str_from_cstr(ImPlot::GetColormapName(l.density_colormap)));
            viamd::write_str(state, MD_STR_LIT("StructureColormap"), md_str_from_cstr(ImPlot::GetColormapName(l.structure_colormap)));
            viamd::write_bool(state, MD_STR_LIT("LogScale"), l.log_scale);
            static const md_str_t colorings[PointColoring_Count] = { MD_STR_LIT("layer"), MD_STR_LIT("structure"), MD_STR_LIT("time") };
            viamd::write_str(state, MD_STR_LIT("PointColoring"), colorings[l.point_coloring]);
            viamd::write_flt(state, MD_STR_LIT("PointSize"), l.point_size);
            static const md_str_t trajectories[Trajectory_Count] = { MD_STR_LIT("lines"), MD_STR_LIT("points"), MD_STR_LIT("both") };
            viamd::write_str(state, MD_STR_LIT("Trajectory"), trajectories[l.trajectory]);
            viamd::write_flt(state, MD_STR_LIT("TrajectoryWidth"), l.trajectory_width);
        }
    }

    void deserialize_settings(viamd::deserialization_state_t& state) {
        md_str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (md_str_eq_cstr(ident, "ShowLayers")) {
                viamd::extract_bool(show_sidebar, arg);
            } else if (md_str_eq_cstr(ident, "ShowGrid")) {
                viamd::extract_bool(show_grid, arg);
            } else if (md_str_eq_cstr(ident, "View")) {
                md_str_t v = {};
                if (viamd::extract_str(v, arg)) {
                    for (int k = 0; k < ViewMode_Count; ++k) {
                        if (md_str_eq_cstr(v, view_mode_names[k])) view_mode = (ViewMode)k;
                    }
                }
            } else if (md_str_eq_cstr(ident, "TimeWindow")) {
                md_str_t v = {};
                if (viamd::extract_str(v, arg)) {
                    time_win.mode = md_str_eq_cstr(v, "own") ? TimeWindow_Own : md_str_eq_cstr(v, "follow") ? TimeWindow_Follow : TimeWindow_Timeline;
                }
            } else if (md_str_eq_cstr(ident, "TimeWindowRange")) {
                double tw[2];
                if (viamd::extract_dbl_vec(tw, 2, arg) && tw[1] >= tw[0]) {
                    time_win.beg = tw[0];
                    time_win.end = tw[1];
                }
            } else if (md_str_eq_cstr(ident, "ShowPeriods")) {
                viamd::extract_bool(time_win.show_periods, arg);
            } else if (md_str_eq_cstr(ident, "Perspective")) {
                viamd::extract_bool(perspective, arg);
            } else if (md_str_eq_cstr(ident, "Smoothing")) {
                if (viamd::extract_int(smoothing, arg)) smoothing = MD_CLAMP(smoothing | 1, 1, MAX_SMOOTHING);
            } else if (md_str_eq_cstr(ident, "Selects")) {
                md_str_t v = {};
                if (viamd::extract_str(v, arg)) {
                    select_target = md_str_eq_cstr(v, "structures") ? SelectTarget_Structures : md_str_eq_cstr(v, "frames") ? SelectTarget_Frames : SelectTarget_Both;
                }
            } else if (md_str_eq_cstr(ident, "Layers")) {
                layers_declared = true;
            }
            // From before there were layers
            else if (md_str_eq_cstr(ident, "Filter")) {
                legacy.present = true;
                viamd::extract_to_char_buf(legacy.filter, sizeof(legacy.filter), arg);
            } else if (md_str_eq_cstr(ident, "MarkerSize")) {
                legacy.present = true;
                viamd::extract_flt(legacy.marker_size, arg);
            } else if (md_str_eq_cstr(ident, "UseMass")) {
                legacy.present = true;
                viamd::extract_bool(legacy.use_mass, arg);
            } else if (md_str_eq_cstr(ident, "Period")) {
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
        md_str_t ident, arg;
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
        md_str_t ident, arg;
        if (!l) {
            while (viamd::next_entry(ident, arg, state)) {}
            return;
        }
        bool has_colormap = false, has_density_colormap = false;
        while (viamd::next_entry(ident, arg, state)) {
            md_str_t s = {};
            if (md_str_eq_cstr(ident, "Name")) {
                viamd::extract_to_char_buf(l->name, sizeof(l->name), arg);
            } else if (md_str_eq_cstr(ident, "Group")) {
                int file_id = 0;
                if (viamd::extract_int(file_id, arg)) {
                    for (int k = 0; k < num_read_groups; ++k) {
                        if (read_group_file_id[k] == file_id) {
                            l->node.parent = read_group_id[k];
                            tree_changed();     // Last in it
                        }
                    }
                }
            } else if (md_str_eq_cstr(ident, "Visible")) {
                viamd::extract_bool(l->visible, arg);
            } else if (md_str_eq_cstr(ident, "Filter")) {
                viamd::extract_to_char_buf(l->filter.text, sizeof(l->filter.text), arg);
            } else if (md_str_eq_cstr(ident, "UseMass")) {
                viamd::extract_bool(l->use_mass, arg);
            } else if (md_str_eq_cstr(ident, "Frames") && viamd::extract_str(s, arg)) {
                // (A layer of the frames selected, "selection", as there was: of all of them)
                l->frames = md_str_eq_cstr(s, "current") ? LayerFrames_Current : md_str_eq_cstr(s, "period") ? LayerFrames_Period : LayerFrames_All;
            } else if (md_str_eq_cstr(ident, "Period")) {
                int id = 0;
                if (viamd::extract_int(id, arg)) l->period = (uint32_t)MD_MAX(id, 0);
            } else if (md_str_eq_cstr(ident, "Style") && viamd::extract_str(s, arg)) {
                l->style = md_str_eq_cstr(s, "contours") ? LayerStyle_Contours : md_str_eq_cstr(s, "heatmap") ? LayerStyle_Heatmap : LayerStyle_Points;
            } else if (md_str_eq_cstr(ident, "Color")) {
                md_float4_t c;
                if (viamd::extract_vec4(c, arg)) l->color = ImVec4(c.x, c.y, c.z, c.w);
            } else if (md_str_eq_cstr(ident, "Opacity")) {
                viamd::extract_flt(l->opacity, arg);
            } else if (md_str_eq_cstr(ident, "Fill")) {
                viamd::extract_bool(l->fill, arg);
            } else if (md_str_eq_cstr(ident, "Lines")) {
                viamd::extract_bool(l->lines, arg);
            } else if (md_str_eq_cstr(ident, "LineWidth")) {
                viamd::extract_flt(l->line_width, arg);
            } else if (md_str_eq_cstr(ident, "Levels")) {
                float v[MAX_LEVELS];
                for (int k = MAX_LEVELS; k > 0; --k) {
                    if (viamd::extract_flt_vec(v, (size_t)k, arg)) {
                        l->num_levels = k;
                        md_memcpy(l->levels, v, sizeof(float) * k);
                        break;
                    }
                }
            } else if (md_str_eq_cstr(ident, "Bandwidth")) {
                if (!viamd::extract_flt(l->bandwidth, arg)) l->bandwidth = 0.0f;  // auto
            } else if (md_str_eq_cstr(ident, "Colormap") && viamd::extract_str(s, arg)) {
                char name[64];
                md_str_copy_to_char_buf(name, sizeof(name), s);
                const ImPlotColormap cm = ImPlot::GetColormapIndex(name);
                if (cm >= 0) l->colormap = cm;
                has_colormap = cm >= 0;
            } else if (md_str_eq_cstr(ident, "DensityColormap") && viamd::extract_str(s, arg)) {
                char name[64];
                md_str_copy_to_char_buf(name, sizeof(name), s);
                const ImPlotColormap cm = ImPlot::GetColormapIndex(name);
                if (cm >= 0) l->density_colormap = cm;
                has_density_colormap = true;
            } else if (md_str_eq_cstr(ident, "StructureColormap") && viamd::extract_str(s, arg)) {
                char name[64];
                md_str_copy_to_char_buf(name, sizeof(name), s);
                const ImPlotColormap cm = ImPlot::GetColormapIndex(name);
                if (cm >= 0) l->structure_colormap = cm;
            } else if (md_str_eq_cstr(ident, "LogScale")) {
                viamd::extract_bool(l->log_scale, arg);
            } else if (md_str_eq_cstr(ident, "PointColoring") && viamd::extract_str(s, arg)) {
                l->point_coloring = md_str_eq_cstr(s, "structure") ? PointColoring_Structure : md_str_eq_cstr(s, "time") ? PointColoring_Time : PointColoring_Layer;
            } else if (md_str_eq_cstr(ident, "PointSize")) {
                viamd::extract_flt(l->point_size, arg);
            } else if (md_str_eq_cstr(ident, "Trajectory") && viamd::extract_str(s, arg)) {
                l->trajectory = md_str_eq_cstr(s, "points") ? Trajectory_Points : md_str_eq_cstr(s, "both") ? Trajectory_Both : Trajectory_Lines;
            } else if (md_str_eq_cstr(ident, "TrajectoryWidth")) {
                viamd::extract_flt(l->trajectory_width, arg);
            }
        }
        // Of a workspace from when a heatmap was in the colours of time: as it was
        if (has_colormap && !has_density_colormap) l->density_colormap = l->colormap;
        sanitize_layer(*l);
        if (!selected_layer) selected_layer = l->id;
    }

    static void sanitize_layer(Layer& l) {
        l.frames = (LayerFrames)MD_CLAMP((int)l.frames, 0, LayerFrames_Count - 1);
        l.style = (LayerStyle)MD_CLAMP((int)l.style, 0, LayerStyle_Count - 1);
        l.opacity = MD_CLAMP(l.opacity, 0.0f, 1.0f);
        l.line_width = MD_CLAMP(l.line_width, 0.5f, 8.0f);
        l.bandwidth = l.bandwidth > 0.0f ? MD_CLAMP(l.bandwidth, 0.005f, MAX_BANDWIDTH) : 0.0f;
        l.point_size = MD_CLAMP(l.point_size, 1.0f, 16.0f);
        l.trajectory = (Trajectory)MD_CLAMP((int)l.trajectory, 0, Trajectory_Count - 1);
        l.trajectory_width = MD_CLAMP(l.trajectory_width, 0.5f, 6.0f);
        l.num_levels = MD_CLAMP(l.num_levels, 1, MAX_LEVELS);
        for (int i = 0; i < l.num_levels; ++i) l.levels[i] = MD_CLAMP(l.levels[i], 0.01f, 0.999f);
        if (l.colormap < 0 || l.colormap >= ImPlot::GetColormapCount()) l.colormap = ImPlotColormap_Spectral;
        if (l.density_colormap < 0 || l.density_colormap >= ImPlot::GetColormapCount()) l.density_colormap = ImPlotColormap_Viridis;
        if (l.structure_colormap < 0 || l.structure_colormap >= ImPlot::GetColormapCount()) l.structure_colormap = ImPlotColormap_Deep;
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
        pts->point_size = MD_CLAMP(legacy.marker_size * 2.0f, 1.0f, 16.0f);
        uint32_t period = (uint32_t)MD_MAX(legacy.period, 0);
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
                snprintf(buf, cap, "%zu %s \xc2\xb7 %u%% evaluated", n, what, (uint32_t)(100.0 * l.eval->frames_done.load() / MD_MAX((size_t)1, l.eval->num_frames)));
                return;
            }
            md_bitset_t frames = {};
            if (layer_frames(l, &frames)) {
                snprintf(buf, cap, "%zu %s \xc2\xb7 %u frames", n, what, (uint32_t)md_bitset_count(frames));
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
        if (l.style == LayerStyle_Heatmap) density_map_draw_colormap(d, l.density_colormap, l.log_scale);
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
        auto u = [&](md_float2_t p) { return (p.x - D.x_min) / (D.x_max - D.x_min); };
        auto v = [&](md_float2_t p) { return (p.y - D.y_min) / (D.y_max - D.y_min); };
        const float clip[6] = { u(tri_lin), v(tri_lin), u(tri_pla), v(tri_pla), u(tri_iso), v(tri_iso) };
        return density_map_composite(&canvas, w, h, view_uv, draws, num_draws, clip);
    }

    // ### Points

    // The points of a layer: the shapes of its structures in the frames drawn, which are the layer's frames with some
    // skipped to keep to a number of points. A point drawn stands for the layer's frames up to the next one drawn.
    struct PointSource {
        const md_float2_t* coords = nullptr;     // [structure * stride + frame]
        size_t   stride = 0;
        size_t   num_structures = 0;
        const uint32_t* frames = nullptr;        // Drawn, ascending, in the frame allocator; of the frame shown {0}
        size_t   num_frames = 0;
        uint32_t frames_end = 0;                 // One past the last frame of the layer's
        md_bitset_t set = {};                    // The layer's frames (for this frame); empty for the frame shown
        uint32_t step = 1;
        bool     current = false;
        bool     smoothed = false;               // The shapes are smoothed over time (smoothing)
        const md_bitset_t* structures = nullptr; // Of the points, when not the layer's (smoothed, of the frame before)
        // Of a window narrower than the layer's frames: the frames the paths would join next past its ends (the ones
        // drawn were the window wider), UINT32_MAX when none. Their points are not of the source: the paths are cut at
        // the window's ends towards them.
        uint32_t before = UINT32_MAX;
        uint32_t after = UINT32_MAX;
        // The first and the last of all of the layer's frames, not only the window's: what the colours of time span, so
        // that a point is in the same colour however much of time is shown, and in 2D and 3D alike
        uint32_t time_beg = 0;
        uint32_t time_last = 0;

        uint32_t frame(size_t k) const { return frames[k]; }
        uint32_t next(size_t k) const { return k + 1 < num_frames ? frames[k + 1] : frames_end; }
    };

    // The points of the layer, at most about max_points of them, of its shapes as they are shown (smoothed or not:
    // eval_shapes). Of the layer's frames only those in [frame_beg, frame_end) (the time window, in 3D), when given.
    bool point_source(PointSource* ps, const Layer& l, size_t max_points, uint32_t frame_beg = 0, uint32_t frame_end = UINT32_MAX) const {
        if (!l.structs_valid) return false;
        static const uint32_t zero = 0;
        // Of the frame shown, smoothed: the smoothed shapes there, between the frames on either side of it as the paths
        // are drawn, so that the points are on them
        if (window_smoothed(l)) {
            if (md_array_size(l.window_coords) == 0) return false;
            *ps = PointSource{};
            ps->coords = l.window_coords;
            ps->stride = 1;
            ps->num_structures = md_array_size(l.window_coords);
            ps->frames = &zero;
            ps->num_frames = 1;
            ps->frames_end = 1;
            ps->current = true;
            ps->smoothed = true;
            ps->structures = l.window_structs;
            return true;
        }
        if (current_smoothed(l) && l.eval && eval_shapes(*l.eval) != l.eval->coords) {
            const Evaluation& e = *l.eval;
            const md_float2_t* c = eval_shapes(e);
            const size_t nf = e.num_frames;
            const double f = MD_CLAMP(app->animation.frame, 0.0, (double)(nf - 1));
            const size_t f0 = (size_t)f, f1 = MD_MIN(f0 + 1, nf - 1);
            const float a = (float)(f - (double)f0);
            md_float2_t* at = (md_float2_t*)md_alloc(app->allocator.frame, sizeof(md_float2_t) * MD_MAX(e.num_structures, (size_t)1));
            for (size_t s = 0; s < e.num_structures; ++s) {
                const md_float2_t p0 = c[s * nf + f0], p1 = c[s * nf + f1];
                at[s] = !coord_valid(p1) ? p0 : !coord_valid(p0) ? p1 : md_float2_t{ p0.x + (p1.x - p0.x) * a, p0.y + (p1.y - p0.y) * a };
            }
            *ps = PointSource{};
            ps->coords = at;
            ps->stride = 1;
            ps->num_structures = e.num_structures;
            ps->frames = &zero;
            ps->num_frames = 1;
            ps->frames_end = 1;
            ps->current = true;
            ps->smoothed = true;
            return true;
        }
        if (l.frames == LayerFrames_Current || run_num_frames(app) < 2) {
            if (md_array_size(l.current_coords) == 0) return false;
            *ps = PointSource{};
            ps->coords = l.current_coords;
            ps->stride = 1;
            ps->num_structures = md_array_size(l.current_coords);
            ps->frames = &zero;
            ps->num_frames = 1;
            ps->frames_end = 1;
            ps->current = true;
            return true;
        }
        if (!l.eval) return false;
        md_bitset_t set = {};
        if (!layer_frames(l, &set)) return false;
        const uint32_t nf = (uint32_t)l.eval->num_frames;
        // Frames skipped are counted from the first of the layer's, not of the window's: as the window moves (following
        // the frame shown) the frames drawn stay the same ones
        uint64_t rank = 0;
        const md_bitset_t all = set;
        const bool windowed = frame_beg > 0 || frame_end < nf;
        if (windowed) {
            const uint32_t e = MD_MIN(frame_end, nf);
            if (frame_beg >= e) return false;
            rank = frame_beg > 0 ? md_bitset_count_range(set, 0, frame_beg) : 0;
            set = md_bitset_and(set, md_bitset_range(frame_beg, e), app->allocator.frame);
            if (md_bitset_empty(set)) return false;
        }
        const uint64_t count = md_bitset_count_range(set, 0, nf);
        if (count == 0) return false;
        const uint64_t total = (uint64_t)l.eval->num_structures * count;
        const uint32_t step = (uint32_t)MD_MAX((uint64_t)1, (total + max_points - 1) / MD_MAX((size_t)1, max_points));
        const size_t cap = (size_t)((count + step - 1) / step);
        uint32_t* list = (uint32_t*)md_alloc(app->allocator.frame, sizeof(uint32_t) * cap);
        size_t n = 0;
        uint64_t i = rank;
        uint32_t last = 0;
        md_bitset_iter_t it = md_bitset_iter_range(set, 0, nf);
        while (md_bitset_iter_next(&it)) {
            last = it.idx;
            if (i++ % step == 0 && n < cap) list[n++] = it.idx;
        }
        *ps = PointSource{};
        ps->coords = eval_shapes(*l.eval);
        ps->smoothed = ps->coords != l.eval->coords;
        ps->stride = nf;
        ps->num_structures = l.eval->num_structures;
        ps->frames = list;
        ps->num_frames = n;
        ps->frames_end = last + 1;
        ps->set = set;
        ps->step = step;
        ps->current = false;
        ps->time_beg = md_bitset_first(all);
        ps->time_last = MD_MAX(md_bitset_last_end(all), 1u) - 1;
        if (windowed && n > 0) {
            // Drawn are the frames of a rank a multiple of step, in all of the layer's: the one before the first drawn,
            // and the one after the last
            const uint64_t r_first = rank + (step - rank % step) % step;
            const uint64_t r_before = r_first >= step ? r_first - step : UINT64_MAX;
            const uint64_t r_after = r_first + (uint64_t)(n - 1) * step + step;
            uint64_t r = 0;
            md_bitset_iter_t ia = md_bitset_iter_range(all, 0, nf);
            while (md_bitset_iter_next(&ia) && r <= r_after) {
                if (r == r_before) ps->before = ia.idx;
                if (r == r_after) ps->after = ia.idx;
                ++r;
            }
        }
        return n > 0;
    }

    // f(structure, k, shape) for each point there is a shape of: in the frame ps.frames[k]
    template <typename F>
    static void for_each_point(const PointSource& ps, F&& f) {
        for (size_t s = 0; s < ps.num_structures; ++s) {
            const md_float2_t* c = ps.coords + s * ps.stride;
            for (size_t k = 0; k < ps.num_frames; ++k) {
                const md_float2_t p = ps.current ? c[0] : c[ps.frames[k]];
                if (coord_valid(p)) f((uint32_t)s, k, p);
            }
        }
    }

    // The frames point k stands for, into m: of the layer's, from its frame to the next drawn. None of the frame shown.
    static void add_point_frames(md_bitmask_t* m, const PointSource& ps, size_t k) {
        if (ps.current || md_bitset_empty(ps.set)) return;
        const uint32_t f = ps.frames[k], e = ps.next(k);
        if (!ps.set.words) {
            const uint32_t b = MD_MAX(f, ps.set.beg), x = MD_MIN(e, ps.set.end);
            if (b < x) md_bitmask_set_range(m, b, x);
            return;
        }
        md_bitset_iter_t it = md_bitset_iter_range(ps.set, f, e);
        while (md_bitset_iter_next(&it)) md_bitmask_set_bit(m, it.idx);
    }

    // Whether any of the frames point k stands for is selected
    bool point_frames_selected(const PointSource& ps, size_t k) const {
        if (ps.current) return false;
        return md_bitmask_test_any_range(&app->selection.frame_mask, ps.frames[k], ps.next(k));
    }

    ImVec4 point_color(const Layer& l, const PointSource& ps, uint32_t s, size_t k) const {
        switch (l.point_coloring) {
        case PointColoring_Structure: {
            // A qualitative map's colours in turn, over again past its last; a continuous one across the structures, in
            // their order
            const ImPlotColormap cm = l.structure_colormap;
            if (ImPlot::ColormapQualitative(cm)) return ImPlot::GetColormapColor((int)(s % (uint32_t)ImPlot::GetColormapSize(cm)), cm);
            return ImPlot::SampleColormap(ps.num_structures > 1 ? (float)s / (float)(ps.num_structures - 1) : 0.5f, cm);
        }
        case PointColoring_Time: {
            if (ps.current) return l.color;
            // Over all of the layer's frames: not of the window of time shown (PointSource::time_beg)
            const uint32_t f0 = ps.time_beg, f1 = ps.time_last;
            if (f1 <= f0) return l.color;
            return ImPlot::SampleColormap((float)(ps.frames[k] - f0) / (float)(f1 - f0), l.colormap);
        }
        default:
            return l.color;
        }
    }

    // The atoms of the structure s of layer l
    md_bitset_t structure_atoms(const Layer& l, const PointSource& ps, uint32_t s) const {
        if (ps.structures) return s < ps.num_structures ? ps.structures[s] : md_bitset_empty_set();
        if (!ps.current && l.eval && s < l.eval->num_structures) return l.eval->structures[s];
        if (s < md_array_size(l.structures)) return l.structures[s];
        return md_bitset_empty_set();
    }

    // Whether the structure is of the atoms selected: with none selected, every structure is, as far as which points are
    // of the selection goes (the frames selected, of every structure)
    bool structure_selected(const Layer& l, const PointSource& ps, uint32_t s, bool any_atoms) const {
        if (!any_atoms) return true;
        const md_bitset_t atoms = structure_atoms(l, ps, s);
        return !md_bitset_empty(atoms) && md_bitmask_test_bit(&app->selection.selection_mask, md_bitset_first(atoms));
    }

    bool selects_structures() const { return select_target != SelectTarget_Frames; }
    bool selects_frames() const { return select_target != SelectTarget_Structures; }

    // What is highlighted from now on in this frame is ours: the atoms', and the frames'
    void begin_highlight() {
        md_bitmask_clear(&app->selection.highlight_mask);
        md_bitmask_clear(&app->selection.frame_highlight_mask);
        highlights_now = true;
        frame_highlights_now = true;
    }

    // Point k of structure s of a layer into the highlight: its structure, and (when frames are selected) its frames
    void highlight_point(const Layer& l, const PointSource& ps, uint32_t s, size_t k) {
        md_bitmask_or(&app->selection.highlight_mask, structure_atoms(l, ps, s));
        if (selects_frames()) add_point_frames(&app->selection.frame_highlight_mask, ps, k);
    }

    // What is highlighted, selected (op): its structures and its frames, as far as the selection is of them
    void commit_highlight(SelectionOperator op) {
        if (selects_structures()) {
            grow_mask_by_selection_granularity(&app->selection.highlight_mask, app->selection.granularity, app->mold.sys);
            modify_field(&app->selection.selection_mask, md_bitmask_view(&app->selection.highlight_mask), op);
        }
        if (selects_frames() && !md_bitmask_empty(&app->selection.frame_highlight_mask)) {
            frame_selection_modify(app, md_bitmask_view(&app->selection.frame_highlight_mask), op);
        }
    }

    // A right click on nothing: the selection cleared, as far as it is of what is selected here
    void clear_selection() {
        if (selects_structures()) md_bitmask_clear(&app->selection.selection_mask);
        if (selects_frames()) frame_selection_clear(app);
    }

    struct Hover {
        int      layer = -1;
        uint32_t structure = 0;
        size_t   k = 0;             // Of the source's frames
        uint32_t frame = 0;
        uint32_t next = 0;
        bool     current = true;
        md_float2_t   coord = {};
    };

    Hover make_hover(int layer, const PointSource& ps, uint32_t s, size_t k, md_float2_t c) const {
        Hover h;
        h.layer = layer;
        h.structure = s;
        h.k = k;
        h.frame = ps.current ? 0 : ps.frames[k];
        h.next = ps.current ? 0 : ps.next(k);
        h.current = ps.current;
        h.coord = c;
        return h;
    }

    void draw_tooltip(const Hover& h, const PointSource& ps) const {
        const Layer& l = layers[h.layer];
        const size_t ns = ps.num_structures;
        const md_float3_t w = shape_weights_of(h.coord);
        const md_bitset_t atoms = structure_atoms(l, ps, h.structure);
        ImGui::BeginTooltip();
        if (ns > 1) ImGui::Text("Structure %u of %zu  \xc2\xb7  %u atoms", h.structure + 1, ns, (uint32_t)md_bitset_count(atoms));
        else        ImGui::Text("%u atoms", (uint32_t)md_bitset_count(atoms));
        ImGui::Text("linear %.2f   planar %.2f   spherical %.2f", w.x, w.y, w.z);
        if (h.current) {
            ImGui::TextDisabled("frame shown");
        } else {
            // The run's time in the display unit, as the timeline shows it
            const double* t = run_frame_times(app);
            char unit[32] = "";
            const double scl = display_units::factor_print(unit, sizeof(unit), run_time_unit(app));
            const bool has_time = t && run_num_frames(app) > h.frame && unit[0] != '\0';
            if (has_time) ImGui::TextDisabled("frame %u (%.6g %s)  \xc2\xb7  click to go there", h.frame, t[h.frame] * scl, unit);
            else          ImGui::TextDisabled("frame %u  \xc2\xb7  click to go there", h.frame);
            if (h.next > h.frame + 1) ImGui::TextDisabled("standing for the frames to %u, which are not drawn", h.next - 1);
        }
        if (ps.smoothed) ImGui::TextDisabled(ICON_FA_WAVE_SQUARE " smoothed over %d frames", smoothing);
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
        auto to_pixels = [&](md_float2_t c) {
            return ImVec2((float)(pos.x + (c.x - lim.X.Min) * sx), (float)(pos.y + (lim.Y.Max - c.y) * sy));
        };

        PointSource sources[MAX_LAYERS];
        bool has_source[MAX_LAYERS] = {};
        for (int i = 0; i < num_layers; ++i) {
            layers[i].point_count = 0;
            if (layer_shown(layers[i]) && layers[i].style == LayerStyle_Points) {
                has_source[i] = point_source(&sources[i], layers[i], MAX_POINTS);
                if (has_source[i]) layers[i].point_stride = sources[i].step;
            }
        }

        // What the pointer is on: the nearest point within reach, of the layer on top
        Hover hover;
        if (hovered && !box.active) {
            const ImVec2 mouse = ImGui::GetMousePos();
            for (int i = 0; i < num_layers && hover.layer < 0; ++i) {
                if (!has_source[i]) continue;
                const float r = MD_MAX(8.0f, layers[i].point_size);
                float best = r * r;
                for_each_point(sources[i], [&](uint32_t s, size_t k, md_float2_t c) {
                    const ImVec2 p = to_pixels(c);
                    const float d2 = (p.x - mouse.x) * (p.x - mouse.x) + (p.y - mouse.y) * (p.y - mouse.y);
                    if (d2 < best) {
                        best = d2;
                        hover = make_hover(i, sources[i], s, k, c);
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

        // The highlight: what the box covers, or the point under the pointer
        if (hovered || box.active) begin_highlight();
        if (box.active && box_is_click) {
            if (box.click_layer >= 0 && box.click_layer < num_layers && box.click_structure >= 0 && has_source[box.click_layer]) {
                const PointSource& ps = sources[box.click_layer];
                md_bitmask_or(highlight_mask, structure_atoms(layers[box.click_layer], ps, (uint32_t)box.click_structure));
                if (selects_frames() && !ps.current && box.click_next > box.click_frame) {
                    md_bitmask_set_range(&app->selection.frame_highlight_mask, box.click_frame, box.click_next);
                }
            }
        } else if (box.active) {
            const double x0 = MD_MIN(box.p0.x, box.p1.x), x1 = MD_MAX(box.p0.x, box.p1.x);
            const double y0 = MD_MIN(box.p0.y, box.p1.y), y1 = MD_MAX(box.p0.y, box.p1.y);
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                // A structure once, however many of its points are in
                md_allocator_i* temp = md_get_heap_allocator();
                md_bitmask_t taken = md_bitmask_create(0, (uint32_t)sources[i].num_structures, temp);
                for_each_point(sources[i], [&](uint32_t s, size_t k, md_float2_t c) {
                    if (x0 <= c.x && c.x <= x1 && y0 <= c.y && c.y <= y1) {
                        if (!md_bitmask_test_bit(&taken, s)) {
                            md_bitmask_set_bit(&taken, s);
                            md_bitmask_or(highlight_mask, structure_atoms(layers[i], sources[i], s));
                        }
                        if (selects_frames()) add_point_frames(&app->selection.frame_highlight_mask, sources[i], k);
                    }
                });
                md_bitmask_free(&taken);
            }
        } else if (hover.layer >= 0) {
            highlight_point(layers[hover.layer], sources[hover.layer], hover.structure, hover.k);
        }

        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        density_map_marks_begin(dl);
        for (int i = num_layers - 1; i >= 0; --i) {
            if (!has_source[i]) continue;
            Layer& l = layers[i];
            const PointSource& ps = sources[i];
            const float r = l.point_size * 0.5f;
            for_each_point(ps, [&](uint32_t s, size_t k, md_float2_t c) {
                ImVec4 col = point_color(l, ps, s, k);
                col.w *= l.opacity;
                const ImU32 u = ImGui::ColorConvertFloat4ToU32(col);
                l.point_count += 1;
                if (ps.current) density_map_mark(dl, to_pixels(c), r + 1.0f, DensityMapMark_Disc, u);
                else            density_map_mark(dl, to_pixels(c), r, DensityMapMark_Square, u);
            });
        }
        const md_float4_t sel_col = app->selection.color.selection.outline;
        const md_float4_t hl_col = app->selection.color.highlight.outline;
        const bool any_sel = !md_bitmask_empty(selection_mask);
        const bool any_hl = !md_bitmask_empty(highlight_mask);
        // The points of the selection over frames, in its colour over the others: of the frames selected, of the
        // structures selected (or of any)
        if (!md_bitmask_empty(&app->selection.frame_mask)) {
            const ImU32 u = ImGui::ColorConvertFloat4ToU32(vec_cast(sel_col));
            for (int i = num_layers - 1; i >= 0; --i) {
                if (!has_source[i] || sources[i].current) continue;
                const Layer& l = layers[i];
                const PointSource& ps = sources[i];
                const float r = l.point_size * 0.5f;
                for (uint32_t s = 0; s < (uint32_t)ps.num_structures; ++s) {
                    if (!structure_selected(l, ps, s, any_sel)) continue;
                    const md_float2_t* c = ps.coords + s * ps.stride;
                    for (size_t k = 0; k < ps.num_frames; ++k) {
                        const md_float2_t p = c[ps.frames[k]];
                        if (coord_valid(p) && point_frames_selected(ps, k)) density_map_mark(dl, to_pixels(p), r, DensityMapMark_Square, u);
                    }
                }
            }
        }
        // Rings around the points of the frame shown of what is selected and highlighted
        for (int i = num_layers - 1; i >= 0 && (any_sel || any_hl); --i) {
            if (!has_source[i] || !sources[i].current) continue;
            const Layer& l = layers[i];
            const float r = MD_MAX(l.point_size * 0.5f + 3.0f, 5.0f);
            for_each_point(sources[i], [&](uint32_t s, size_t, md_float2_t c) {
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
            density_map_mark(dl, to_pixels(hover.coord), MD_MAX(layers[hover.layer].point_size * 0.5f + 3.0f, 5.0f), DensityMapMark_Ring, ImGui::ColorConvertFloat4ToU32(vec_cast(hl_col)));
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
            box.click_frame = hover.frame;
            box.click_next = hover.next;
        }
        if (box.active) {
            box.p1 = mouse;
            if (!box_is_click) {
                const ImVec2 a = ImPlot::PlotToPixels(box.p0);
                const ImVec2 b = ImPlot::PlotToPixels(box.p1);
                const md_float4_t c = box.op == SelectionOperator::Or ? sel_col : md_float4_t{1.0f, 0.4f, 0.4f, 0.9f};
                ImPlot::PushPlotClipRect();
                dl->AddRectFilled(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 0.12f)));
                dl->AddRect(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(vec_cast(c)), 0.0f, 0, 1.0f);
                ImPlot::PopPlotClipRect();
            }
            const ImGuiMouseButton button = box.op == SelectionOperator::Or ? ImGuiMouseButton_Left : ImGuiMouseButton_Right;
            if (ImGui::IsMouseReleased(button) || !ImGui::IsMouseDown(button)) {
                if (!box_is_click || box.click_structure >= 0) {
                    commit_highlight(box.op);
                } else if (box.op == SelectionOperator::AndNot) {
                    clear_selection();
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
        auto px = [](md_float2_t p) { return ImPlot::PlotToPixels(ImPlotPoint(p.x, p.y)); };
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
            if (group == 0 && like->node.parent == 0) md_memcpy(filter, like->filter.text, sizeof(filter));
            use_mass = like->use_mass;
        }
        Layer* l = add_layer(frames, style, 0, group);
        if (!l) return nullptr;
        md_memcpy(l->filter.text, filter, sizeof(filter));
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
        const double s = MD_MAX(1.0, MD_MIN((size.x - 2.0 * mx) / 1.0, (size.y - mt - mb) / h_tri));     // Pixels per unit
        const double w = size.x / s;
        const double h = size.y / s;
        const double extra = (size.y - s * h_tri - mt - mb) * 0.5;
        const double y0 = -(mb + MD_MAX(0.0, extra)) / s;
        return ImPlotRect(0.5 - w * 0.5, 0.5 + w * 0.5, y0, y0 + h);
    }

    void draw_plot() {
        constexpr ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
        constexpr ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoDecorations;

        ImPlotInputMap& input = ImPlot::GetInputMap();
        const ImPlotInputMap saved_input = input;
        // (Shift is for selecting, so the plot does not pan or zoom meanwhile: plot_input.h)
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
            // Nothing hovered while it is panned (plot_input.h)
            draw_points(plot_input_hovered());
            {
                ImPlot::PushPlotClipRect();
                draw_smoothing_badge(ImPlot::GetPlotDrawList(), ImPlot::GetPlotPos());
                ImPlot::PopPlotClipRect();
            }

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
                const double eps = 1e-6 * MD_MAX(1e-6, fitted_rect.X.Size());
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

    // ### 3D: the shapes over time
    //
    // In the plot's coordinates x is time (the timeline's), y and z the triangle's plane: the triangle stood up across
    // the time axis.

    // The time of a frame in the display unit (time_scale, set where the 3D plot is drawn), as the timeline shows it;
    // the frame itself without times
    double time_scale = 1.0;
    double frame_time(double frame) const {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        if (!t || n == 0) return frame;
        return time_at_frame(t, n, frame) * time_scale;
    }

    // A point drawn in 3D: of structure s in the frame ps.frames[k] (or the frame shown), where it is on the screen
    struct Point3D {
        uint32_t s;
        uint32_t k;
        float    x, y, z;
        ImVec2   px;
    };

    // Whether a line goes from point a to b, the next of the points of the layer: of the same structure, in the next frame
    // drawn (no frame without a shape between them), not across frames the layer leaves out (a gap in the frames selected)
    static bool joined(const PointSource& ps, const Point3D& a, const Point3D& b) {
        if (a.s != b.s || b.k != a.k + 1) return false;
        return ps.frames[b.k] - ps.frames[a.k] <= ps.step;
    }

    static ImU32 color_u32(ImVec4 c, float alpha) {
        c.w *= alpha;
        return ImGui::ColorConvertFloat4ToU32(c);
    }

    // ### The time window of the 3D plot
    //
    // The time the 3D plot shows, and cuts the paths to, is a window of the run: the timeline's visible range (linked,
    // both ways), one of its own, or one that follows the frame shown. All of the run is the window when nothing is cut.

    // The window, in the run's own unit
    void time_window(double* tb, double* te) const {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        const double r0 = (t && n > 0) ? t[0] : 0.0;
        const double r1 = (t && n > 0) ? t[n - 1] : 0.0;
        switch (time_win.mode) {
        case TimeWindow_Timeline:
            *tb = MD_CLAMP(app->timeline.view_range.beg_x, r0, r1);
            *te = MD_CLAMP(app->timeline.view_range.end_x, *tb, r1);
            break;
        case TimeWindow_Follow: {
            const double tc = run_time_at_frame(app, app->animation.frame);
            *tb = tc + time_win.beg;
            *te = tc + time_win.end;
            break;
        }
        default:
            if (time_win.end > time_win.beg) {
                *tb = time_win.beg;
                *te = time_win.end;
            } else {
                *tb = r0;
                *te = r1;
            }
            break;
        }
    }

    // The frames whose times are in the window, [*fb, *fe): empty when it is beside the run
    void time_window_frames(uint32_t* fb, uint32_t* fe) const {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        *fb = *fe = 0;
        if (!t || n == 0) return;
        double tb, te;
        time_window(&tb, &te);
        if (te < t[0] || tb > t[n - 1] || te < tb) return;
        // The first frame at or after tb, and one past the last at or before te
        uint32_t b = (uint32_t)ceil(frame_at_time(t, n, tb) - 1e-9);
        uint32_t e = (uint32_t)floor(frame_at_time(t, n, te) + 1e-9);
        while (e + 1 < n && t[e + 1] <= te) ++e;
        *fb = b;
        *fe = MD_MAX(b, e + 1);
    }

    // Another mode, from the window as it is: what is shown stays shown. Following, as long a tail behind the frame
    // shown; else the window as times (linked, into the timeline's range).
    void time_window_set_mode(TimeWindow mode) {
        if (mode == time_win.mode) return;
        double tb, te;
        time_window(&tb, &te);
        time_win.mode = mode;
        if (mode == TimeWindow_Follow) {
            time_win.beg = -(te - tb);
            time_win.end = 0.0;
        } else {
            time_window_set(tb, te);
        }
    }

    // The window set to the times [b, e], in the run's own unit, as zooming, panning and the range in the toolbar move
    // it: at least two frames long. Linked and of its own, within the run (all of it at most); following, as offsets.
    void time_window_set(double b, double e) {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        if (!t || n < 2) return;
        const double span = t[n - 1] - t[0];
        const double min_len = 2.0 * span / (double)(n - 1);
        if (e - b < min_len) {
            const double c = 0.5 * (b + e);
            b = c - 0.5 * min_len;
            e = c + 0.5 * min_len;
        }
        if (time_win.mode == TimeWindow_Follow) {
            const double tc = run_time_at_frame(app, app->animation.frame);
            time_win.beg = MD_CLAMP(b - tc, -span, span);
            time_win.end = MD_CLAMP(e - tc, time_win.beg, span);
            return;
        }
        if (e - b >= span) {
            b = t[0];
            e = t[n - 1];
        }
        // Moved back into the run, as long as it is
        if (b < t[0])     { e += t[0] - b; b = t[0]; }
        if (e > t[n - 1]) { b -= e - t[n - 1]; e = t[n - 1]; }
        if (time_win.mode == TimeWindow_Timeline) {
            app->timeline.view_range.beg_x = b;
            app->timeline.view_range.end_x = e;
        } else {
            time_win.beg = b;
            time_win.end = e;
        }
    }

    // All of it: the run (linked, the timeline's too), or following, all of the run behind the frame shown
    void time_window_fit() {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        if (!t || n < 2) return;
        if (time_win.mode == TimeWindow_Follow) {
            time_win.beg = -(t[n - 1] - t[0]);
            time_win.end = 0.0;
        } else {
            time_window_set(t[0], t[n - 1]);
        }
    }

    // The toolbar's controls of the window: its mode, and its ends in the display unit, as the timeline shows time
    void draw_time_window_controls() {
        const double* t = run_frame_times(app);
        const size_t n = run_num_frames(app);
        if (!t || n < 2) return;
        char unit[32] = "";
        const double scl = display_units::factor_print(unit, sizeof(unit), run_time_unit(app));

        static const char* const icons[TimeWindow_Count] = { ICON_FA_LINK, ICON_FA_THUMBTACK, ICON_FA_PERSON_RUNNING };
        static const char* const tips[TimeWindow_Count] = {
            "Time in 3D: the timeline's visible range, both ways (zoom and pan either)",
            "Time in 3D: a window of its own",
            "Time in 3D: following the frame shown (a tail behind it when it ends at 0)",
        };
        for (int k = 0; k < TimeWindow_Count; ++k) {
            if (k > 0) ImGui::SameLine(0.0f, 1.0f);
            const bool on = time_win.mode == k;
            if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
            if (ImGui::Button(icons[k])) time_window_set_mode((TimeWindow)k);
            if (on) ImGui::PopStyleColor();
            ImGui::SetItemTooltip("%s", tips[k]);
        }

        const bool follow = time_win.mode == TimeWindow_Follow;
        const double span = t[n - 1] - t[0];
        double tb, te;
        time_window(&tb, &te);
        if (follow) {
            const double tc = run_time_at_frame(app, app->animation.frame);
            tb -= tc;
            te -= tc;
        }
        const double lo_bound = follow ? -span : t[0];
        const double hi_bound = follow ?  span : t[n - 1];
        float v[2] = { (float)(tb * scl), (float)(te * scl) };
        char fmt[48];
        snprintf(fmt, sizeof(fmt), follow ? "%%+.6g%s%s" : "%%.6g%s%s", unit[0] ? " " : "", unit);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(MD_CLAMP(ImGui::GetContentRegionAvail().x * 0.35f, 120.0f, 260.0f));
        if (ImGui::DragFloatRange2("##time_window", &v[0], &v[1], (float)(span * scl / 500.0), (float)(lo_bound * scl), (float)(hi_bound * scl), fmt, fmt, ImGuiSliderFlags_AlwaysClamp)) {
            double b = v[0] / scl, e = MD_MAX(v[1], v[0]) / scl;
            if (follow) {
                const double tc = run_time_at_frame(app, app->animation.frame);
                b += tc;
                e += tc;
            }
            time_window_set(b, e);
        }
        if (ImGui::BeginPopupContextItem("##time_window_menu")) {
            if (ImGui::MenuItem(follow ? "All of the run behind the frame shown" : "The whole run")) time_window_fit();
            ImGui::EndPopup();
        }
        uint32_t fb, fe;
        time_window_frames(&fb, &fe);
        if (fe > fb) ImGui::SetItemTooltip("Frames %u to %u, %u frames. Drag the ends or ctrl-click to type; right-click for all of it.", fb, fe - 1, fe - fb);
        else         ImGui::SetItemTooltip("No frames in the window. Drag the ends or ctrl-click to type; right-click for all of it.");
    }

    // In a corner of a plot, while what it shows is smoothed over time: said with the data, so that no view of it (a
    // screenshot of the plot alone) passes for the shapes of the frames
    void draw_smoothing_badge(ImDrawList* dl, ImVec2 corner) const {
        if (smoothing <= 1 || run_num_frames(app) < 2) return;
        char text[64];
        snprintf(text, sizeof(text), ICON_FA_WAVE_SQUARE " smoothed over %d frames", smoothing);
        const ImVec2 pad(5.0f, 2.0f);
        const ImVec2 p0 = corner + ImVec2(6.0f, 6.0f);
        const ImVec2 p1 = p0 + ImGui::CalcTextSize(text) + pad * 2.0f;
        dl->AddRectFilled(p0, p1, IM_COL32(0, 0, 0, 150), 3.0f);
        dl->AddRect(p0, p1, ImGui::GetColorU32(ImGuiCol_ButtonActive), 3.0f);
        dl->AddText(p0 + pad, ImGui::GetColorU32(ImGuiCol_Text), text);
    }

    // The smoothing over time, in the toolbar: in sight, as it changes what every point is. Marked while it is on.
    void draw_smoothing_control() {
        const bool on = smoothing > 1;
        if (on) {
            ImGui::PushStyleColor(ImGuiCol_FrameBg, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
            ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, ImGui::GetStyleColorVec4(ImGuiCol_ButtonHovered));
        }
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6.0f);
        if (ImGui::SliderInt("##smoothing", &smoothing, 1, MAX_SMOOTHING, on ? ICON_FA_WAVE_SQUARE " %d" : ICON_FA_WAVE_SQUARE " off",
                             ImGuiSliderFlags_AlwaysClamp | ImGuiSliderFlags_Logarithmic)) {
            smoothing |= 1;     // Centred: odd
            // Made now, before the points are drawn this frame: not at the next update
            for (int i = 0; i < num_layers; ++i) {
                if (layers[i].eval) update_smoothed(*layers[i].eval);
            }
        }
        if (on) ImGui::PopStyleColor(2);
        if (ImGui::BeginPopupContextItem("##smoothing_menu")) {
            if (ImGui::MenuItem("Off", nullptr, false, on)) smoothing = 1;
            ImGui::EndPopup();
        }
        ImGui::SetItemTooltip("Smoothing over time: a moving average of the shapes over this many frames, centred.\n"
                              "Of every layer and view alike: the points and densities in 2D, the paths in 3D, the points\n"
                              "of the frame shown, what is selected and what is exported. As smooth(shape_weights(...), n)\n"
                              "in a script. Ctrl-click to type; right-click to turn it off.");
    }

    void draw_plot3d() {
        v3.drawn = true;
        if (ImGui::GetContentRegionAvail().x < 32.0f || ImGui::GetContentRegionAvail().y < 32.0f) return;
        md_allocator_i* fa = app->allocator.frame;
        const uint32_t nf = (uint32_t)run_num_frames(app);
        char unit[32] = "";
        time_scale = display_units::factor_print(unit, sizeof(unit), run_time_unit(app));
        // The run, and the time axis: the time window, which the paths are cut to
        const double r0 = frame_time(0.0);
        double r1 = nf > 1 ? frame_time((double)(nf - 1)) : r0 + 1.0;
        if (!(r1 > r0)) r1 = r0 + 1.0;
        const double tc = nf > 0 ? frame_time(app->animation.frame) : r0;
        uint32_t filter_beg = 0, filter_end = UINT32_MAX;
        double t0 = r0, t1 = r1;
        if (nf > 1) {
            time_window_frames(&filter_beg, &filter_end);
            double tb, te;
            time_window(&tb, &te);
            t0 = tb * time_scale;
            t1 = te * time_scale;
            if (!(t1 > t0)) t1 = t0 + (r1 - r0) * 1e-3;
        }
        // The prism of the run, within the axis
        const double p0 = MD_MAX(r0, t0), p1 = MD_MIN(r1, t1);

        char x_label[64] = "Frame";
        if (unit[0] != '\0') snprintf(x_label, sizeof(x_label), "Time (%s)", unit);

        const ImGuiIO& io = ImGui::GetIO();
        const bool shift = io.KeyShift;
        const bool ctrl = io.KeyCtrl && !shift;
        // ImPlot3D only turns the box, with the left button as the viewport does; zooming and panning (in time, of the time
        // filter's window) are ours. Shift is for selecting: the box does not turn meanwhile.
        ImPlot3DFlags flags = ImPlot3DFlags_CanvasOnly | ImPlot3DFlags_NoMenus | ImPlot3DFlags_NoPan | ImPlot3DFlags_NoZoom;
        if (shift || ctrl || v3.active || v3.panning || v3.dollying || v3.scrubbing) flags |= ImPlot3DFlags_NoRotate;
        // Turning is ImPlot3D's, as every 3D view's (plot_input.h); panning, fitting and zooming in time are ours, of the
        // window, and there is no menu
        ImPlot3DInputMap& input_map = ImPlot3D::GetInputMap();
        const ImPlot3DInputMap saved_map = input_map;
        input_map.Pan = -1;
        input_map.Fit = -1;
        input_map.Menu = -1;
        // A double-click on a plane of the box turns to face it (ImPlot3D); off them it is ours alone, of time: the box
        // is not turned back. Of the last frame, as ImPlot3D takes it before the box is known: a double-click is still.
        if (!v3.over_box) input_map.ResetRotation = -1;

        ImPlot3D::PushStyleColor(ImPlot3DCol_PlotBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));
        ImPlot3D::PushStyleColor(ImPlot3DCol_FrameBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));
        ImPlot3D::PushStyleVar(ImPlot3DStyleVar_ViewScaleFactor, v3.view_scale);
        ImPlot3D::PushStyleVar(ImPlot3DStyleVar_Perspective, perspective ? PERSPECTIVE_FOV : 0.0f);
        MD_DEFER {
            ImPlot3D::PopStyleVar(2);
            ImPlot3D::PopStyleColor(2);
            ImPlot3D::GetInputMap() = saved_map;
        };

        if (!ImPlot3D::BeginPlot("##Shapespace3D", ImVec2(-1, -1), flags)) return;

        // The triangle all of it, always: the shape axes are fixed. The time axis is the time window.
        const ImPlot3DAxisFlags shape_axis = ImPlot3DAxisFlags_NoLabel | ImPlot3DAxisFlags_NoTickLabels | ImPlot3DAxisFlags_NoTickMarks |
                                             ImPlot3DAxisFlags_NoGridLines | ImPlot3DAxisFlags_Lock;
        ImPlot3D::SetupAxes(x_label, nullptr, nullptr, ImPlot3DAxisFlags_None, shape_axis, shape_axis);
        ImPlot3D::SetupAxesLimits(t0, t1, -0.04, 1.04, -0.04, tri_iso.y + 0.04, ImPlot3DCond_Always);
        // The triangle equilateral, time longer
        ImPlot3D::SetupBoxScale(1.3, 0.7, 0.7 * tri_iso.y);
        ImPlot3D::SetupBoxInitialRotation(20.0, -35.0);
        if (v3.preset >= 0) {
            // Along the time axis (the triangle as in 2D), from the side, at an angle
            static const double el[3] = { 0.0, 0.0, 20.0 };
            static const double az[3] = { -90.0, 0.0, -35.0 };
            ImPlot3D::SetupBoxRotation(el[v3.preset], az[v3.preset], true, ImPlot3DCond_Always);
            v3.preset = -1;
        } else {
            ImPlot3D::SetupBoxRotation(20.0, -35.0, false, ImPlot3DCond_Once);
        }
        // (The setup is finished by the first call that is not of it: GetPlotRectPos)

        const ImVec2 rect_pos = ImPlot3D::GetPlotRectPos();
        const ImVec2 rect_size = ImPlot3D::GetPlotRectSize();
        const ImRect plot_rect(rect_pos, rect_pos + rect_size);
        const ImVec2 mouse = ImGui::GetMousePos();
        // The plot is an item that takes the press (to turn and pan): hovered while it holds the mouse too
        const ImGuiContext& g = *ImGui::GetCurrentContext();
        const bool held_elsewhere = g.ActiveId != 0 && g.ActiveIdWindow != ImGui::GetCurrentWindow();
        const bool hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem) && !held_elsewhere && plot_rect.Contains(mouse);

        const ImPlot3DItemFlags item_flags = ImPlot3DItemFlags_NoLegend | ImPlot3DItemFlags_NoFit;
        const md_float2_t corners[3] = { tri_lin, tri_pla, tri_iso };

        // The prism of the run: the triangle at its first and last frames (of those in the window), its corners along time
        {
            const double pa = p1 > p0 ? p0 : t0, pb = p1 > p0 ? p1 : t1;
            double xs[18], ys[18], zs[18];
            int n = 0;
            auto seg = [&](double xa, md_float2_t a, double xb, md_float2_t b) {
                xs[n] = xa; ys[n] = a.x; zs[n] = a.y; ++n;
                xs[n] = xb; ys[n] = b.x; zs[n] = b.y; ++n;
            };
            for (int c = 0; c < 3; ++c) seg(pa, corners[c], pb, corners[c]);
            ImPlot3DSpec spec;
            spec.LineColor = ImVec4(1, 1, 1, 0.22f);
            spec.LineWeight = 1.0f;
            spec.Flags = item_flags | ImPlot3DLineFlags_Segments;
            ImPlot3D::PlotLine("##prism", xs, ys, zs, n, spec);

            // The reference triangle at either end, with its grid of equal weights (as in 2D): where the shapes are
            double gx[2 * 15], gy[2 * 15], gz[2 * 15];
            double ex[2 * 3], ey[2 * 3], ez[2 * 3];
            for (int end = 0; end < 2; ++end) {
                const double x = end == 0 ? pa : pb;
                int ng = 0, ne = 0;
                auto line = [&](double* lx, double* ly, double* lz, int& k, md_float2_t a, md_float2_t b) {
                    lx[k] = x; ly[k] = a.x; lz[k] = a.y; ++k;
                    lx[k] = x; ly[k] = b.x; lz[k] = b.y; ++k;
                };
                for (int c = 0; c < 3; ++c) line(ex, ey, ez, ne, corners[c], corners[(c + 1) % 3]);
                if (show_grid) {
                    auto lerp2 = [](md_float2_t a, md_float2_t b, float u) { return md_float2_t{ a.x + (b.x - a.x) * u, a.y + (b.y - a.y) * u }; };
                    const md_float2_t a = corners[0], b = corners[1], c = corners[2];
                    for (int k = 1; k < 4; ++k) {
                        const float u = k / 4.0f;
                        line(gx, gy, gz, ng, lerp2(a, b, u), lerp2(a, c, u));
                        line(gx, gy, gz, ng, lerp2(b, a, u), lerp2(b, c, u));
                        line(gx, gy, gz, ng, lerp2(c, a, u), lerp2(c, b, u));
                    }
                }
                char id[24];
                ImPlot3DSpec es;
                es.LineColor = ImVec4(1, 1, 1, 0.45f);
                es.LineWeight = 1.0f;
                es.Flags = item_flags | ImPlot3DLineFlags_Segments;
                snprintf(id, sizeof(id), "##tri%d", end);
                ImPlot3D::PlotLine(id, ex, ey, ez, ne, es);
                if (ng > 0) {
                    es.LineColor = ImVec4(1, 1, 1, 0.12f);
                    snprintf(id, sizeof(id), "##grid%d", end);
                    ImPlot3D::PlotLine(id, gx, gy, gz, ng, es);
                }
            }
            ImPlot3D::PlotText("Linear", pa, tri_lin.x, tri_lin.y, 0.0, ImVec2(0, 12));
            ImPlot3D::PlotText("Planar", pa, tri_pla.x, tri_pla.y, 0.0, ImVec2(0, 12));
            ImPlot3D::PlotText("Spherical", pa, tri_iso.x, tri_iso.y, 0.0, ImVec2(0, -12));
        }

        // The frame shown, a slice across; the period under the pointer (in the timeline or a layer's row), two faces
        auto slice = [&](const char* id, double t, ImVec4 line, ImVec4 fill) {
            double xs[3] = { t, t, t }, ys[3], zs[3];
            for (int c = 0; c < 3; ++c) { ys[c] = corners[c].x; zs[c] = corners[c].y; }
            ImPlot3DSpec spec;
            spec.LineColor = line;
            spec.FillColor = fill;
            spec.FillAlpha = 1.0f;
            spec.LineWeight = 1.5f;
            spec.Flags = item_flags | ImPlot3DTriangleFlags_NoMarkers;
            ImPlot3D::PlotTriangle(id, xs, ys, zs, 3, spec);
        };
        if (nf > 0) slice("##current", tc, ImVec4(1, 1, 0, 0.75f), ImVec4(1, 1, 0, 0.05f));
        // The timeline's periods along time, as it shades them: a band each across a side of the box along the edge the
        // time axis is drawn on (as the box is turned), out to that axis, a line for a marker; named along the axis' edge,
        // by its ticks. The one hovered (in the timeline, or a layer's) brighter.
        const uint32_t hot_period = timeline_periods_hovered();
        if (time_win.show_periods && nf > 1) {
            const double y_lo = -0.04, y_hi = 1.04, z_lo = -0.04, z_hi = tri_iso.y + 0.04;
            // The edge the time axis is drawn along, as the box is turned; looking down the time axis there is none: the
            // floor's near edge
            ImPlot3DPoint ea, eb;
            double y_edge = y_lo, z_edge = z_lo;
            bool names = false;     // Not where time runs (nearly) at the viewer: they would be on top of each other
            if (ImPlot3D::GetAxisEdge(ImAxis3D_X, &ea, &eb)) {
                y_edge = ea.y < 0.5 ? y_lo : y_hi;
                z_edge = ea.z < 0.5 * (z_lo + z_hi) ? z_lo : z_hi;
                names = ImLengthSqr(ImPlot3D::PlotToPixels(ea) - ImPlot3D::PlotToPixels(eb)) > 120.0f * 120.0f;
            }
            // Of the two sides of the box along that edge (across the triangle's plane, z = z_edge, or upright, y = y_edge),
            // one behind what is in the box (its outside facing away, as ImPlot3D's grids are), and of two such the one seen
            // more of: the floor from above, the far wall from the side
            const ImVec2 centre = ImPlot3D::GetPlotRectPos() + ImPlot3D::GetPlotRectSize() * 0.5f;
            const ImPlot3DPoint view = ImPlot3D::PixelsToPlotRay(centre).Direction;    // Into the box
            const double floor_out = z_edge < 0.5 * (z_lo + z_hi) ? -view.z : view.z;  // Its outside's way, along the view
            const double wall_out = y_edge < 0.5 ? -view.y : view.y;
            auto face_area = [&](bool floor) {
                const ImVec2 a = ImPlot3D::PlotToPixels(t0, floor ? y_lo : y_edge, floor ? z_edge : z_lo);
                const ImVec2 b = ImPlot3D::PlotToPixels(t1, floor ? y_lo : y_edge, floor ? z_edge : z_lo);
                const ImVec2 c = ImPlot3D::PlotToPixels(t0, floor ? y_hi : y_edge, floor ? z_edge : z_hi);
                const ImVec2 u = b - a, v = c - a;
                return fabsf(u.x * v.y - u.y * v.x);
            };
            const bool floor_behind = floor_out >= -1e-6, wall_behind = wall_out >= -1e-6;
            const bool on_floor = floor_behind != wall_behind ? floor_behind : face_area(true) >= face_area(false);
            // A point of the side at time x, a fraction f across it from the far edge to the axis' edge
            auto side_point = [&](double f, double* y, double* z) {
                if (on_floor) {
                    *y = y_edge < 0.5 ? y_hi + (y_lo - y_hi) * f : y_lo + (y_hi - y_lo) * f;
                    *z = z_edge;
                } else {
                    *y = y_edge;
                    *z = z_edge < 0.5 * (z_lo + z_hi) ? z_hi + (z_lo - z_hi) * f : z_lo + (z_hi - z_lo) * f;
                }
            };
            double sy0, sz0, sy1, sz1, ny, nz;
            side_point(0.0, &sy0, &sz0);
            side_point(1.0, &sy1, &sz1);
            side_point(0.94, &ny, &nz);     // The names a little within the side from the axis' edge, clear of the ticks
            const TimelinePeriods& rs = app->timeline.periods;
            for (int i = 0; i < rs.count; ++i) {
                const TimelinePeriod& p = rs.items[i];
                int b = 0, e = 0;
                if (!timeline_period_frames(&b, &e, p, app->animation.frame, (int)nf)) continue;
                const double xb = frame_time(b), xe = frame_time(e);
                if (xe < t0 || xb > t1) continue;
                const bool hot = p.id == hot_period;
                const ImVec4 col(p.color.x, p.color.y, p.color.z, hot ? 0.95f : 0.7f);
                char id[32];
                snprintf(id, sizeof(id), "##period%u", p.id);
                ImPlot3DSpec spec;
                spec.Flags = item_flags;
                const double xa = MD_MAX(xb, t0), xz = MD_MIN(xe, t1);
                if (b == e) {
                    const double xs[2] = { xb, xb }, ys[2] = { sy0, sy1 }, zs[2] = { sz0, sz1 };
                    spec.LineColor = col;
                    spec.LineWeight = hot ? 2.5f : 1.5f;
                    ImPlot3D::PlotLine(id, xs, ys, zs, 2, spec);
                } else {
                    const double xs[4] = { xa, xz, xz, xa }, ys[4] = { sy0, sy0, sy1, sy1 }, zs[4] = { sz0, sz0, sz1, sz1 };
                    spec.FillColor = ImVec4(p.color.x, p.color.y, p.color.z, hot ? 0.35f : 0.16f);
                    spec.FillAlpha = 1.0f;
                    spec.LineColor = col;
                    spec.LineWeight = 1.0f;
                    ImPlot3D::PlotQuad(id, xs, ys, zs, 4, spec);
                }
                if (names) ImPlot3D::PlotText(p.name, 0.5 * (xa + xz), ny, nz, 0.0, ImVec2(0, 0));
            }
        }
        if (const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, hot_period)) {
            int b = 0, e = 0;
            if (timeline_period_frames(&b, &e, *p, app->animation.frame, (int)nf)) {
                const ImVec4 col(p->color.x, p->color.y, p->color.z, 0.8f);
                const ImVec4 fill(p->color.x, p->color.y, p->color.z, 0.10f);
                slice("##period_beg", frame_time(b), col, fill);
                if (e != b) slice("##period_end", frame_time(e), col, fill);
            }
        }

        // f(structure, k, a, b) for the paths out of the window, to the points past its ends: from the point k inside, a,
        // to where the path to the point outside crosses the window's end, b ({ time, x, y }), as joined() joins points
        auto for_each_edge = [&](const PointSource& ps, auto&& f) {
            if (ps.current || ps.num_frames == 0) return;
            for (int end = 0; end < 2; ++end) {
                const uint32_t out = end == 0 ? ps.before : ps.after;
                if (out == UINT32_MAX) continue;
                const size_t k = end == 0 ? 0 : ps.num_frames - 1;
                const uint32_t in = ps.frames[k];
                if ((end == 0 ? in - out : out - in) > ps.step) continue;
                const double ti = frame_time(in), to = frame_time(out);
                const double te = end == 0 ? t0 : t1;
                if (!(to != ti)) continue;
                const double u = MD_CLAMP((te - ti) / (to - ti), 0.0, 1.0);
                for (size_t st = 0; st < ps.num_structures; ++st) {
                    const md_float2_t pi = ps.coords[st * ps.stride + in], po = ps.coords[st * ps.stride + out];
                    if (!coord_valid(pi) || !coord_valid(po)) continue;
                    const float a[3] = { (float)ti, pi.x, pi.y };
                    const float b[3] = { (float)(ti + (to - ti) * u), pi.x + (po.x - pi.x) * (float)u, pi.y + (po.y - pi.y) * (float)u };
                    f((uint32_t)st, k, a, b);
                }
            }
        };

        // The layers of points: over frames a line (or points) for each structure, of the frame shown points on its slice
        PointSource sources[MAX_LAYERS];
        bool has_source[MAX_LAYERS] = {};
        Point3D* points[MAX_LAYERS] = {};
        size_t num_points[MAX_LAYERS] = {};
        const bool any_sel = !md_bitmask_empty(&app->selection.selection_mask);
        const bool any_frames_sel = !md_bitmask_empty(&app->selection.frame_mask);
        for (int i = num_layers - 1; i >= 0; --i) {
            Layer& l = layers[i];
            if (!layer_shown(l) || l.style != LayerStyle_Points || l.opacity <= 0.0f) continue;
            if (!point_source(&sources[i], l, MAX_LINE_POINTS, filter_beg, filter_end)) continue;
            has_source[i] = true;
            const PointSource& ps = sources[i];
            if (!ps.current) l.line_stride = ps.step;

            // Every point that is, in the plot's coordinates
            const size_t cap = ps.num_structures * ps.num_frames;
            Point3D* pts = (Point3D*)md_alloc(fa, sizeof(Point3D) * MD_MAX(cap, (size_t)1));
            size_t np = 0;
            for_each_point(ps, [&](uint32_t s, size_t k, md_float2_t c) {
                pts[np++] = Point3D{ s, (uint32_t)k, (float)(ps.current ? tc : frame_time(ps.frames[k])), c.x, c.y, ImVec2() };
            });
            points[i] = pts;
            num_points[i] = np;
            if (np == 0) continue;

            // Of what is not selected, while something is, a fainter colour: as in the viewport
            auto alpha_of = [&](uint32_t s) { return l.opacity * (structure_selected(l, ps, s, any_sel) ? 1.0f : 0.3f); };

            char id[32];
            const bool lines = !ps.current && l.trajectory != Trajectory_Points;
            const bool marks = ps.current || l.trajectory != Trajectory_Lines;
            if (lines) {
                // Segments between the points of a structure one after the other, in the colour of the first
                const size_t max_seg = np;
                float* xs = (float*)md_alloc(fa, sizeof(float) * 2 * max_seg);
                float* ys = (float*)md_alloc(fa, sizeof(float) * 2 * max_seg);
                float* zs = (float*)md_alloc(fa, sizeof(float) * 2 * max_seg);
                ImU32* cols = (ImU32*)md_alloc(fa, sizeof(ImU32) * max_seg);
                int ns = 0;
                for (size_t j = 0; j + 1 < np; ++j) {
                    const Point3D& a = pts[j];
                    const Point3D& b = pts[j + 1];
                    if (!joined(ps, a, b)) continue;
                    xs[2 * ns] = a.x; ys[2 * ns] = a.y; zs[2 * ns] = a.z;
                    xs[2 * ns + 1] = b.x; ys[2 * ns + 1] = b.y; zs[2 * ns + 1] = b.z;
                    cols[ns] = color_u32(point_color(l, ps, a.s, a.k), alpha_of(a.s));
                    ++ns;
                }
                // Out to the points past the window, cut at its ends: in the colour of the point inside
                float* exs = (float*)md_alloc(fa, sizeof(float) * 4 * ps.num_structures);
                float* eys = (float*)md_alloc(fa, sizeof(float) * 4 * ps.num_structures);
                float* ezs = (float*)md_alloc(fa, sizeof(float) * 4 * ps.num_structures);
                ImU32* ecols = (ImU32*)md_alloc(fa, sizeof(ImU32) * 2 * ps.num_structures);
                int ne = 0;
                for_each_edge(ps, [&](uint32_t s, size_t k, const float* a, const float* b) {
                    exs[2 * ne] = a[0]; eys[2 * ne] = a[1]; ezs[2 * ne] = a[2];
                    exs[2 * ne + 1] = b[0]; eys[2 * ne + 1] = b[1]; ezs[2 * ne + 1] = b[2];
                    ecols[ne] = color_u32(point_color(l, ps, s, k), alpha_of(s));
                    ++ne;
                });
                if (ne > 0) {
                    ImPlot3DSpec spec;
                    spec.LineColors = ecols;
                    spec.LineWeight = l.trajectory_width;
                    spec.Flags = item_flags | ImPlot3DLineFlags_Segments;
                    snprintf(id, sizeof(id), "##edges%u", l.id);
                    ImPlot3D::PlotLine(id, exs, eys, ezs, 2 * ne, spec);
                }
                if (ns > 0) {
                    ImPlot3DSpec spec;
                    spec.LineColors = cols;
                    spec.LineWeight = l.trajectory_width;
                    spec.Flags = item_flags | ImPlot3DLineFlags_Segments;
                    snprintf(id, sizeof(id), "##lines%u", l.id);
                    ImPlot3D::PlotLine(id, xs, ys, zs, 2 * ns, spec);
                }
            }
            if (marks) {
                float* xs = (float*)md_alloc(fa, sizeof(float) * np);
                float* ys = (float*)md_alloc(fa, sizeof(float) * np);
                float* zs = (float*)md_alloc(fa, sizeof(float) * np);
                ImU32* cols = (ImU32*)md_alloc(fa, sizeof(ImU32) * np);
                for (size_t j = 0; j < np; ++j) {
                    xs[j] = pts[j].x; ys[j] = pts[j].y; zs[j] = pts[j].z;
                    cols[j] = color_u32(point_color(l, ps, pts[j].s, pts[j].k), alpha_of(pts[j].s));
                }
                ImPlot3DSpec spec;
                spec.Marker = ps.current ? ImPlot3DMarker_Circle : ImPlot3DMarker_Square;
                spec.MarkerSize = ps.current ? l.point_size * 0.5f + 1.0f : l.point_size * 0.5f;
                spec.MarkerFillColors = cols;
                spec.MarkerLineColors = cols;
                spec.FillAlpha = 1.0f;
                spec.Flags = item_flags;
                snprintf(id, sizeof(id), "##marks%u", l.id);
                ImPlot3D::PlotScatter(id, xs, ys, zs, (int)np, spec);
            }
        }

        // Where the points are on the screen, for what the pointer is on and what a box covers
        Hover hover;
        // Nothing hovered while it is turned, panned or zoomed by a drag (plot_input.h); the box still covers what it covers
        const bool hover_points = hovered && !plot_input_navigating() && !v3.scrubbing;
        if (hover_points || v3.active) {
            float best = 0.0f;
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                const float r = MD_MAX(8.0f, layers[i].point_size);
                if (hover.layer < 0) best = r * r;
                for (size_t j = 0; j < num_points[i]; ++j) {
                    Point3D& p = points[i][j];
                    p.px = ImPlot3D::PlotToPixels(p.x, p.y, p.z);
                    if (!hover_points || v3.active || (hover.layer >= 0 && hover.layer < i)) continue;
                    const float d2 = ImLengthSqr(p.px - mouse);
                    if (d2 < best) {
                        best = d2;
                        hover = make_hover(i, sources[i], p.s, p.k, { p.y, p.z });
                    }
                }
            }
        }

        // The highlight: what the box covers, or the point under the pointer
        const bool box_is_click = v3.active && ImLengthSqr(mouse - v3.p0) < 9.0f;
        if (hover_points || v3.active) begin_highlight();
        if (v3.active && box_is_click) {
            if (v3.click_layer >= 0 && v3.click_layer < num_layers && v3.click_structure >= 0 && has_source[v3.click_layer]) {
                const PointSource& ps = sources[v3.click_layer];
                md_bitmask_or(&app->selection.highlight_mask, structure_atoms(layers[v3.click_layer], ps, (uint32_t)v3.click_structure));
                if (selects_frames() && !ps.current && v3.click_next > v3.click_frame) {
                    md_bitmask_set_range(&app->selection.frame_highlight_mask, v3.click_frame, v3.click_next);
                }
            }
        } else if (v3.active) {
            const ImRect r(ImMin(v3.p0, mouse), ImMax(v3.p0, mouse));
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                md_bitmask_t taken = md_bitmask_create(0, (uint32_t)sources[i].num_structures, md_get_heap_allocator());
                for (size_t j = 0; j < num_points[i]; ++j) {
                    const Point3D& p = points[i][j];
                    if (!r.Contains(p.px)) continue;
                    if (!md_bitmask_test_bit(&taken, p.s)) {
                        md_bitmask_set_bit(&taken, p.s);
                        md_bitmask_or(&app->selection.highlight_mask, structure_atoms(layers[i], sources[i], p.s));
                    }
                    if (selects_frames()) add_point_frames(&app->selection.frame_highlight_mask, sources[i], p.k);
                }
                md_bitmask_free(&taken);
            }
        } else if (hover.layer >= 0) {
            highlight_point(layers[hover.layer], sources[hover.layer], hover.structure, hover.k);
        }

        // What is highlighted, and the points of the selection over frames, over the layers: a line in the highlight's
        // colour through the points of each structure highlighted, and marks on the points of the frames selected
        {
            const bool any_hl = !md_bitmask_empty(&app->selection.highlight_mask);
            const ImVec4 hl_col = vec_cast(app->selection.color.highlight.outline);
            const ImVec4 sel_col = vec_cast(app->selection.color.selection.outline);
            for (int i = 0; i < num_layers && (any_hl || any_frames_sel); ++i) {
                if (!has_source[i] || num_points[i] == 0) continue;
                const Layer& l = layers[i];
                const PointSource& ps = sources[i];
                const Point3D* pts = points[i];
                const size_t np = num_points[i];
                auto structure_highlighted = [&](uint32_t s) {
                    const md_bitset_t atoms = structure_atoms(l, ps, s);
                    return !md_bitset_empty(atoms) && md_bitmask_test_bit(&app->selection.highlight_mask, md_bitset_first(atoms));
                };
                char id[32];
                if (any_hl && !ps.current) {
                    const size_t cap = 2 * (np + 2 * ps.num_structures);
                    float* xs = (float*)md_alloc(fa, sizeof(float) * cap);
                    float* ys = (float*)md_alloc(fa, sizeof(float) * cap);
                    float* zs = (float*)md_alloc(fa, sizeof(float) * cap);
                    int ns = 0;
                    if (l.trajectory != Trajectory_Points) {
                        for_each_edge(ps, [&](uint32_t s, size_t, const float* a, const float* b) {
                            if (!structure_highlighted(s)) return;
                            xs[2 * ns] = a[0]; ys[2 * ns] = a[1]; zs[2 * ns] = a[2];
                            xs[2 * ns + 1] = b[0]; ys[2 * ns + 1] = b[1]; zs[2 * ns + 1] = b[2];
                            ++ns;
                        });
                    }
                    uint32_t last_s = UINT32_MAX;
                    bool last_hl = false;
                    for (size_t j = 0; j + 1 < np; ++j) {
                        const Point3D& a = pts[j];
                        const Point3D& b = pts[j + 1];
                        if (a.s != last_s) {
                            last_s = a.s;
                            last_hl = structure_highlighted(a.s);
                        }
                        if (!last_hl || !joined(ps, a, b)) continue;
                        xs[2 * ns] = a.x; ys[2 * ns] = a.y; zs[2 * ns] = a.z;
                        xs[2 * ns + 1] = b.x; ys[2 * ns + 1] = b.y; zs[2 * ns + 1] = b.z;
                        ++ns;
                    }
                    if (ns > 0) {
                        ImPlot3DSpec spec;
                        spec.LineColor = hl_col;
                        spec.LineWeight = l.trajectory_width + 1.5f;
                        spec.Flags = item_flags | ImPlot3DLineFlags_Segments;
                        snprintf(id, sizeof(id), "##hl%u", l.id);
                        ImPlot3D::PlotLine(id, xs, ys, zs, 2 * ns, spec);
                    }
                }
                // Of the selection over frames, in its colour over the rest: the stretches of the lines in the frames
                // selected (of the structures selected, or of any), and the points where the layer draws points
                if (any_frames_sel && !ps.current) {
                    auto selected = [&](const Point3D& p) { return point_frames_selected(ps, p.k) && structure_selected(l, ps, p.s, any_sel); };
                    if (l.trajectory != Trajectory_Points) {
                        float* xs = (float*)md_alloc(fa, sizeof(float) * 2 * np);
                        float* ys = (float*)md_alloc(fa, sizeof(float) * 2 * np);
                        float* zs = (float*)md_alloc(fa, sizeof(float) * 2 * np);
                        int ns = 0;
                        for (size_t j = 0; j + 1 < np; ++j) {
                            const Point3D& a = pts[j];
                            const Point3D& b = pts[j + 1];
                            if (!joined(ps, a, b) || !selected(a)) continue;
                            xs[2 * ns] = a.x; ys[2 * ns] = a.y; zs[2 * ns] = a.z;
                            xs[2 * ns + 1] = b.x; ys[2 * ns + 1] = b.y; zs[2 * ns + 1] = b.z;
                            ++ns;
                        }
                        if (ns > 0) {
                            ImPlot3DSpec spec;
                            spec.LineColor = sel_col;
                            spec.LineWeight = l.trajectory_width + 1.0f;
                            spec.Flags = item_flags | ImPlot3DLineFlags_Segments;
                            snprintf(id, sizeof(id), "##sel_lines%u", l.id);
                            ImPlot3D::PlotLine(id, xs, ys, zs, 2 * ns, spec);
                        }
                    }
                    if (l.trajectory != Trajectory_Lines) {
                        float* xs = (float*)md_alloc(fa, sizeof(float) * np);
                        float* ys = (float*)md_alloc(fa, sizeof(float) * np);
                        float* zs = (float*)md_alloc(fa, sizeof(float) * np);
                        int n = 0;
                        for (size_t j = 0; j < np; ++j) {
                            const Point3D& p = pts[j];
                            if (!selected(p)) continue;
                            xs[n] = p.x; ys[n] = p.y; zs[n] = p.z;
                            ++n;
                        }
                        if (n > 0) {
                            ImPlot3DSpec spec;
                            spec.Marker = ImPlot3DMarker_Square;
                            spec.MarkerSize = l.point_size * 0.5f;
                            spec.MarkerLineColor = sel_col;
                            spec.MarkerFillColor = sel_col;
                            spec.FillAlpha = 1.0f;
                            spec.Flags = item_flags;
                            snprintf(id, sizeof(id), "##sel_marks%u", l.id);
                            ImPlot3D::PlotScatter(id, xs, ys, zs, n, spec);
                        }
                    }
                }
            }
        }

        const bool any_points = [&] { for (int i = 0; i < num_layers; ++i) if (num_points[i] > 0) return true; return false; }();
        // Where time runs on the screen (along the middle of the prism), and how much of the plot the box takes: of this
        // frame, for zooming, panning and the box's size
        const md_float2_t mid = { 0.5f, tri_iso.y / 3.0f };
        const ImVec2 time_a = ImPlot3D::PlotToPixels(t0, mid.x, mid.y);
        const ImVec2 time_b = ImPlot3D::PlotToPixels(t1, mid.x, mid.y);
        const ImVec2 time_m = ImPlot3D::PlotToPixels(0.5 * (t0 + t1), mid.x, mid.y);
        ImRect box_px(ImVec2(FLT_MAX, FLT_MAX), ImVec2(-FLT_MAX, -FLT_MAX));
        for (int k = 0; k < 8; ++k) {
            const ImVec2 c = ImPlot3D::PlotToPixels(k & 1 ? t1 : t0, k & 2 ? 1.04 : -0.04, k & 4 ? tri_iso.y + 0.04 : -0.04);
            box_px.Add(c);
        }

        // Ctrl + left-click or drag: to the frame at the pointer's time, as in the timeline (plot_input.h). The time is
        // where the pointer is on a wall of the box along time, behind what is in it (where the ray through the pointer
        // leaves the box: the walls ImPlot3D grids), or else on the time axis, by its ticks. None through an end of the
        // box (looking along time), or off the box and its axis.
        // Where the ray through the pointer leaves the box (its parameter), and through which of its faces (the axis
        // across it, -1 when it misses): one of the planes ImPlot3D draws behind what is in the box, as it is turned
        const ImPlot3DRay ray = ImPlot3D::PixelsToPlotRay(mouse);
        int exit_face = -1;
        const double exit_t = [&] {
            const double lo[3] = { t0, -0.04, -0.04 }, hi[3] = { t1, 1.04, tri_iso.y + 0.04 };
            const double o[3] = { ray.Origin.x, ray.Origin.y, ray.Origin.z };
            const double d[3] = { ray.Direction.x, ray.Direction.y, ray.Direction.z };
            double t_in = -DBL_MAX, t_out = DBL_MAX;
            int out = -1;
            for (int k = 0; k < 3; ++k) {
                if (fabs(d[k]) < 1e-12) {
                    if (o[k] < lo[k] || o[k] > hi[k]) t_in = DBL_MAX;    // Beside the box, along it: a miss
                    continue;
                }
                double a = (lo[k] - o[k]) / d[k], b = (hi[k] - o[k]) / d[k];
                if (a > b) { const double tmp = a; a = b; b = tmp; }
                t_in = MD_MAX(t_in, a);
                if (b < t_out) { t_out = b; out = k; }
            }
            if (!(t_in <= t_out)) return (double)NAN;
            exit_face = out;
            return t_out;
        }();
        v3.over_box = hovered && exit_face >= 0;
        // Off the box, a double-click zooms out fully: the whole run. On one of its planes ImPlot3D turns to face it, and
        // the time shown stays.
        if (hovered && !shift && !ctrl && !v3.active && exit_face < 0 && nf > 1 && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) time_window_fit();
        auto pointer_time = [&]() -> double {
            if (exit_face > 0) return ray.Origin.x + ray.Direction.x * exit_t;
            ImPlot3DPoint ea, eb;
            if (!ImPlot3D::GetAxisEdge(ImAxis3D_X, &ea, &eb)) return NAN;
            const ImVec2 a = ImPlot3D::PlotToPixels(ea), b = ImPlot3D::PlotToPixels(eb);
            const ImVec2 m = ImPlot3D::PlotToPixels((ea + eb) * 0.5);
            const ImVec2 ab = b - a;
            const float len2 = ImLengthSqr(ab);
            if (len2 < 400.0f) return NAN;
            const float l = ImDot(mouse - a, ab) / len2;
            const float reach = 3.0f * ImGui::GetFontSize();    // Over the ticks and their labels
            if (l < 0.0f || l > 1.0f || ImLengthSqr(mouse - (a + ab * l)) > reach * reach) return NAN;
            // Foreshortened in perspective, as the time axis is for zooming
            const float lm = ImDot(m - a, ab) / len2;
            const float r = lm > 0.01f && lm < 0.99f ? 1.0f / lm - 1.0f : 1.0f;
            const double u = l * r / (1.0f - l + l * r);
            return ea.x + (eb.x - ea.x) * u;
        };
        if (nf > 1 && ctrl && hovered && !v3.active && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && isfinite(pointer_time())) {
            v3.scrubbing = true;
            v3.pressed = false;
        }
        if (v3.scrubbing && (!ImGui::IsMouseDown(ImGuiMouseButton_Left) || !io.KeyCtrl)) v3.scrubbing = false;
        if (v3.scrubbing) {
            const double t = pointer_time();
            if (isfinite(t)) {
                const double* times = run_frame_times(app);
                const double rt = MD_CLAMP(t, r0, r1) / time_scale;
                const double f = times ? frame_at_time(times, nf, rt) : rt;
                app->animation.frame = MD_CLAMP(f, 0.0, (double)(nf - 1));
            }
        }
        ImPlot3D::EndPlot();

        // The box fills the plot, with room for the time axis' labels. Not while it is turned by hand: it settles after.
        if (box_px.GetWidth() > 1.0f && box_px.GetHeight() > 1.0f && !(ImGui::IsMouseDown(ImGuiMouseButton_Left) && ImGui::IsItemActive())) {
            const float fit = MD_MIN(plot_rect.GetWidth() * 0.86f / box_px.GetWidth(), plot_rect.GetHeight() * 0.80f / box_px.GetHeight());
            const float target = MD_CLAMP(v3.view_scale * fit, 0.2f, 8.0f);
            v3.view_scale += (target - v3.view_scale) * 0.5f;
        }

        // Zooming and panning in time: of the time window (linked, the timeline's too), about the time under the pointer
        {
            const ImVec2 axis = time_b - time_a;
            const float len2 = ImLengthSqr(axis);
            const bool along = len2 > 400.0f;     // Not looking (nearly) down the time axis
            const double w0 = t0, w1 = t1;
            // Of a point on the screen, how far along the window the time axis is under it, in [0, 1]. The axis drawn is
            // foreshortened in perspective: its middle (time_m) is not halfway along it on the screen, and gives by how
            // much (the ratio of the depths of its ends). Seen from afar, that is where the point is along it on the screen.
            const float lm = len2 > 0.0f ? ImDot(time_m - time_a, axis) / len2 : 0.5f;
            const float r = lm > 0.01f && lm < 0.99f ? 1.0f / lm - 1.0f : 1.0f;
            auto time_u = [&](ImVec2 p) {
                const float l = MD_CLAMP(ImDot(p - time_a, axis) / len2, 0.0f, 1.0f);
                return l * r / (1.0f - l + l * r);
            };
            if (hovered && !shift && !v3.active && nf > 1) {
                ImGui::SetKeyOwner(ImGuiKey_MouseWheelY, ImGui::GetCurrentWindow()->ID);
                if (io.MouseWheel != 0.0f) {
                    const float u = along ? time_u(mouse) : 0.5f;
                    const double tm = w0 + u * (w1 - w0);
                    const double f = pow(0.85, (double)io.MouseWheel);
                    time_window_set((tm - (tm - w0) * f) / time_scale, (tm + (w1 - tm) * f) / time_scale);
                }
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) v3.panning = true;
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    v3.dollying = true;
                    v3.dolly_u = along ? time_u(mouse) : 0.5f;
                }
            }
            // Middle-drag, as the viewport dollies: up zooms in
            if (v3.dollying) {
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
                    v3.dollying = false;
                } else if (io.MouseDelta.y != 0.0f && nf > 1) {
                    const double tm = w0 + v3.dolly_u * (w1 - w0);
                    const double f = exp(0.01 * (double)io.MouseDelta.y);
                    time_window_set((tm - (tm - w0) * f) / time_scale, (tm + (w1 - tm) * f) / time_scale);
                }
            }
            if (v3.panning) {
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
                    v3.panning = false;
                } else if (along && (io.MouseDelta.x != 0.0f || io.MouseDelta.y != 0.0f)) {
                    // The paths move with the pointer: the window the other way. Past the ends of the axis too, as the
                    // pointer leaves it, but not as far as where the axis would vanish in perspective
                    auto u_of = [&](ImVec2 p) {
                        const float l = MD_CLAMP(ImDot(p - time_a, axis) / len2, -0.5f, 1.5f);
                        return l * r / MD_MAX(1.0f - l + l * r, 0.05f);
                    };
                    const double dt = -(double)(u_of(mouse) - u_of(mouse - io.MouseDelta)) * (w1 - w0);
                    time_window_set((w0 + dt) / time_scale, (w1 + dt) / time_scale);
                }
            }
        }

        // Over the plot: the point under the pointer, the box
        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->PushClipRect(plot_rect.Min, plot_rect.Max, true);
        draw_smoothing_badge(dl, plot_rect.Min);
        const md_float4_t sel_c = app->selection.color.selection.outline;
        const md_float4_t hl_c = app->selection.color.highlight.outline;
        if (hover.layer >= 0) {
            // Found again for where it is: the point's pixel position is of this frame
            for (size_t j = 0; j < num_points[hover.layer]; ++j) {
                const Point3D& p = points[hover.layer][j];
                if (p.s == hover.structure && p.k == hover.k) {
                    dl->AddCircle(p.px, MD_MAX(layers[hover.layer].point_size * 0.5f + 3.0f, 5.0f), ImGui::ColorConvertFloat4ToU32(vec_cast(hl_c)), 0, 1.5f);
                    break;
                }
            }
        }
        if (v3.active && !box_is_click) {
            const md_float4_t c = v3.op == SelectionOperator::Or ? sel_c : md_float4_t{1.0f, 0.4f, 0.4f, 0.9f};
            const ImVec2 a = ImMin(v3.p0, mouse), b = ImMax(v3.p0, mouse);
            dl->AddRectFilled(a, b, ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 0.12f)));
            dl->AddRect(a, b, ImGui::ColorConvertFloat4ToU32(vec_cast(c)), 0.0f, 0, 1.0f);
        }
        if (!any_points && num_layers > 0) {
            const char* msg = nf > 1 && filter_end <= filter_beg ? "No frames in the time shown" : "Layers of points show here, over time";
            const ImVec2 ts = ImGui::CalcTextSize(msg);
            dl->AddText(ImVec2(plot_rect.Min.x + 8.0f, plot_rect.Max.y - ts.y - 6.0f), ImGui::GetColorU32(ImGuiCol_TextDisabled), msg);
        }
        dl->PopClipRect();

        if (hover.layer >= 0 && !v3.active) {
            draw_tooltip(hover, sources[hover.layer]);
            // Pressed on a point of a frame and let go there: that frame (a drag pans)
            if (!hover.current && !shift && !ctrl && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                v3.pressed = true;
                v3.press_pos = mouse;
                v3.click_frame = hover.frame;
            }
        }
        if (v3.pressed && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            if (ImLengthSqr(mouse - v3.press_pos) < 9.0f) app->animation.frame = v3.click_frame;
            v3.pressed = false;
        }
        if (v3.pressed && !ImGui::IsMouseDown(ImGuiMouseButton_Left)) v3.pressed = false;

        // Selecting: shift-drag a box, shift-click a point; the right button takes away; a right click on nothing clears
        if (!v3.active && hovered && shift && (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
            v3.active = true;
            v3.pressed = false;
            v3.op = ImGui::IsMouseClicked(ImGuiMouseButton_Left) ? SelectionOperator::Or : SelectionOperator::AndNot;
            v3.p0 = v3.p1 = mouse;
            v3.click_layer = hover.layer;
            v3.click_structure = hover.layer >= 0 ? (int64_t)hover.structure : -1;
            v3.click_frame = hover.frame;
            v3.click_next = hover.next;
        }
        if (v3.active) {
            v3.p1 = mouse;
            const ImGuiMouseButton button = v3.op == SelectionOperator::Or ? ImGuiMouseButton_Left : ImGuiMouseButton_Right;
            if (ImGui::IsMouseReleased(button) || !ImGui::IsMouseDown(button)) {
                if (!box_is_click || v3.click_structure >= 0) {
                    commit_highlight(v3.op);
                } else if (v3.op == SelectionOperator::AndNot) {
                    clear_selection();
                }
                v3.active = false;
                v3.click_layer = -1;
                v3.click_structure = -1;
            }
        }
    }

    // ### The window

    void export_layer(const Layer& l, md_str_t ext) {
        const EvaluationRef& eval = l.eval;
        if (!eval || !eval->complete()) return;
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) return;
        const md_str_t path = { path_buf, strnlen(path_buf, sizeof(path_buf)) };

        md_temp_scope_t temp_scope = md_temp_begin();
        md_allocator_i* temp = md_temp_allocator(temp_scope);
        MD_DEFER { md_temp_end(temp_scope); };

        // The time column is the run's time in the display unit, as the timeline shows it
        char unit_buf[64];
        const double time_scl = display_units::factor_print(unit_buf, sizeof(unit_buf), run_time_unit(app));
        md_str_t x_label = MD_STR_INIT("Frame");
        if (unit_buf[0] != '\0') {
            x_label = md_alloc_printf(temp, "Time (%s)", unit_buf);
        }
        md_array(md_str_t)  labels = 0;
        md_array(float*) values = 0;
        md_array_push(labels, x_label, temp);
        md_array_push(values, NULL, temp);
        const size_t ns = eval->num_structures;
        const size_t nf = eval->num_frames;
        const md_float2_t* shapes = eval_shapes(*eval);     // As shown: smoothed, when they are
        // Said in every column's label, so that a file of smoothed shapes cannot be taken for the shapes of the frames
        char sm[48] = "";
        if (shapes != eval->coords) snprintf(sm, sizeof(sm), " (smoothed over %d frames)", smoothing);
        for (size_t s = 0; s < ns; ++s) {
            static const char* const kinds[3] = { "linear", "planar", "isotropic" };
            for (int k = 0; k < 3; ++k) {
                md_array_push(labels, ns > 1 ? md_alloc_printf(temp, "%s %zu%s", kinds[k], s + 1, sm) : md_alloc_printf(temp, "%s%s", kinds[k], sm), temp);
            }
            for (int k = 0; k < 3; ++k) md_array_push(values, NULL, temp);
        }
        const double* time = run_frame_times(app);
        const bool has_time = time && run_num_frames(app) >= nf;
        for (size_t f = 0; f < nf; ++f) {
            md_array_push(values[0], has_time ? (float)(time[f] * time_scl) : (float)f, temp);
            for (size_t s = 0; s < ns; ++s) {
                const md_float3_t w = shape_weights_of(shapes[s * nf + f]);
                md_array_push(values[1 + s * 3 + 0], w.x, temp);
                md_array_push(values[1 + s * 3 + 1], w.y, temp);
                md_array_push(values[1 + s * 3 + 2], w.z, temp);
            }
        }
        const size_t num_cols = ns * 3 + 1;
        if (md_str_eq(ext, MD_STR_LIT("csv"))) {
            md_csv_write_to_file(values, labels, num_cols, nf, path);
        } else {
            const md_str_t header = md_xvg_format_header(MD_STR_LIT("Shape Space"), x_label, MD_STR_LIT("Shape Weight"), num_cols - 1, labels + 1, temp);
            const md_str_t xvg = md_xvg_format(header, num_cols, nf, values, temp);
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
            if (ImGui::MenuItem("CSV", nullptr, false, can_export)) export_layer(*sel, MD_STR_LIT("csv"));
            if (ImGui::MenuItem("XVG", nullptr, false, can_export)) export_layer(*sel, MD_STR_LIT("xvg"));
            if (!can_export) ImGui::TextDisabled("Select a layer over frames, once evaluated");
            ImGui::EndPopup();
        }

        // The plot in 2D, the shapes over time in 3D, or both
        ImGui::SameLine();
        {
            int vm = (int)view_mode;
            static const int values[] = { ViewMode_2D, ViewMode_3D, ViewMode_Split };
            static const char* const tips[] = { "The triangle", "The triangle over time: a line through the frames of each structure", "Both, side by side" };
            density_map_segmented("##view", view_mode_names, values, ViewMode_Count, &vm, tips, false);
            view_mode = (ViewMode)vm;
        }
        if (view_mode != ViewMode_2D) {
            ImGui::SameLine();
            if (ImGui::Button(ICON_FA_CUBE " " ICON_FA_CARET_DOWN)) ImGui::OpenPopup("##view3d");
            ImGui::SetItemTooltip("Turn the 3D view");
            if (ImGui::BeginPopup("##view3d")) {
                if (ImGui::MenuItem("Along the time axis")) v3.preset = 0;
                ImGui::SetItemTooltip("The triangle as in 2D, every frame behind the first");
                if (ImGui::MenuItem("From the side")) v3.preset = 1;
                ImGui::SetItemTooltip("Time across, the corners up and down");
                if (ImGui::MenuItem("At an angle")) v3.preset = 2;
                ImGui::Separator();
                ImGui::MenuItem("Perspective", nullptr, &perspective);
                ImGui::SetItemTooltip("Seen from nearby: what is nearer is larger. Off, seen from afar: sizes and parallels\n"
                                      "are kept, and the lengths along an axis compare anywhere in the box.");
                ImGui::EndPopup();
            }
            ImGui::SameLine();
            draw_time_window_controls();
        }
        if (run_num_frames(app) >= 2) {
            ImGui::SameLine();
            draw_smoothing_control();
        }

        ImGui::SameLine();
        const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + ImGui::CalcTextSize(ICON_FA_GEAR).x + ImGui::GetStyle().FramePadding.x * 4 + ImGui::GetStyle().ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MD_MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
        if (ImGui::Button(ICON_FA_EXPAND)) {
            fit_view = true;
            if (view_mode != ViewMode_2D) time_window_fit();
        }
        ImGui::SetItemTooltip(view_mode == ViewMode_2D ? "Fit the view (or double-click the plot)" : "Fit the view, in 3D all of the run (or double-click off the box)");
        ImGui::SameLine();
        if (ImGui::Button(ICON_FA_GEAR)) ImGui::OpenPopup("##settings");
        if (ImGui::BeginPopup("##settings")) {
            ImGui::Checkbox("Grid of equal weights", &show_grid);
            ImGui::SeparatorText("A selection selects");
            for (int k = 0; k < SelectTarget_Count; ++k) {
                if (ImGui::RadioButton(select_target_names[k], select_target == k)) select_target = (SelectTarget)k;
            }
            ImGui::Separator();
            ImGui::TextDisabled("Click a point of a frame to go there.");
            ImGui::TextDisabled("Shift-drag to select the structures and frames of the points, with the right button to deselect them.");
            ImGui::TextDisabled("Drop a period from the timeline on the plot to add its density.");
            ImGui::TextDisabled("In 3D: drag turns; in time, right-drag pans and the wheel zooms; double-click off the box shows the\n"
                               "whole run, on a plane turns to face it; ctrl-click or drag on the walls\n"
                               "or the time axis goes to the frame there.");
            ImGui::TextDisabled("The time shown in 3D is the timeline's (" ICON_FA_LINK "), its own (" ICON_FA_THUMBTACK ") or a tail following the frame shown (" ICON_FA_PERSON_RUNNING ").");
            ImGui::Checkbox("Periods along time in 3D", &time_win.show_periods);
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
        const bool by_map = l.style == LayerStyle_Points && l.point_coloring != PointColoring_Layer;
        const ImPlotColormap sw_map = l.style != LayerStyle_Points ? l.density_colormap
                                    : l.point_coloring == PointColoring_Structure ? l.structure_colormap : l.colormap;
        density_map_swatch(dl, sw_min, sw_max, by_map ? DensityMapSwatch_Heatmap : kind, l.color, l.fill, sw_map, l.visible ? 1.0f : 0.35f);

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
        const size_t nk = representation_tree_children(kids, MD_COUNTOF(kids), nodes, num_nodes, g.node.id);
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
            md_memcpy(saved, l.name, sizeof(saved));
            l.name[0] = '\0';
            layer_label(auto_name, sizeof(auto_name), l);
            md_memcpy(l.name, saved, sizeof(saved));
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
            filter_input("##filter", &l.filter, app->script.language, layer_group_filter_context(app, in_context ? g : nullptr), 0, nullptr, none);
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
        if (l.style == LayerStyle_Heatmap) {
            prop::row("colormap", "Of the density");
            ImPlot::ColormapSelection("##dcmap", &l.density_colormap);
        } else if (l.style == LayerStyle_Points && l.point_coloring == PointColoring_Time) {
            prop::row("colormap", "Of time: over all of the layer's frames");
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
            if (l.point_coloring == PointColoring_Structure) {
                prop::row("colormap", "Of the structures. A qualitative map gives them its colours in turn, again from the first\n"
                                      "past its last; a continuous one spreads them across it, in their order.");
                ImPlot::ColormapSelection("##scmap", &l.structure_colormap);
            }
            if (l.frames != LayerFrames_Current && l.point_stride > 1) {
                prop::hint("Every %u frames, to keep to %u points", l.point_stride, (uint32_t)MAX_POINTS);
            }
            if (l.frames != LayerFrames_Current) {
                prop::row("in 3D", "Over time: a line through the frames of each structure, its points, or both");
                int t = (int)l.trajectory;
                static const int values[] = { Trajectory_Lines, Trajectory_Points, Trajectory_Both };
                density_map_segmented("##trajectory", trajectory_names, values, Trajectory_Count, &t, nullptr, true);
                l.trajectory = (Trajectory)t;
                if (l.trajectory != Trajectory_Points) {
                    prop::row("line width");
                    ImGui::SliderFloat("##traj_width", &l.trajectory_width, 0.5f, 6.0f, "%.1f px");
                }
                if (l.line_stride > 1) prop::hint("In 3D every %u frames, to keep to %u points", l.line_stride, (uint32_t)MAX_LINE_POINTS);
            }
        }

        ImGui::PopID();
        prop::end();
    }

    void draw_sidebar() {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("LAYERS");
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MD_MAX(0.0f, ImGui::GetContentRegionAvail().x - ImGui::GetFrameHeight()));
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
            const size_t nt = representation_tree_children(top, MD_COUNTOF(top), nodes, n, 0);
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
        frame_highlights_now = false;
        plot_drawn = false;
        v3.drawn = false;
        draw_window();
        // Collapsed, hidden or closed in the middle of a box: no selection from it once it is back
        if (!plot_drawn && box.active) {
            box.active = false;
            box.click_layer = -1;
            box.click_structure = -1;
        }
        if (!v3.drawn) {
            v3.active = false;
            v3.pressed = false;
        }
        if (app && owns_highlight && !highlights_now) md_bitmask_clear(&app->selection.highlight_mask);
        owns_highlight = highlights_now;
        if (app && owns_frame_highlight && !frame_highlights_now) md_bitmask_clear(&app->selection.frame_highlight_mask);
        owns_frame_highlight = frame_highlights_now;
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
            if (view_mode == ViewMode_2D || view_mode == ViewMode_Split) {
                const float w = view_mode == ViewMode_Split ? MD_MAX(120.0f, ImGui::GetContentRegionAvail().x * 0.42f) : 0.0f;
                // Split, a child of its own: its width is the one resized to, not that of the 2D plot alone
                ImGui::BeginChild(view_mode == ViewMode_Split ? "##plot_split" : "##plot", ImVec2(w, 0), view_mode == ViewMode_Split ? ImGuiChildFlags_ResizeX : ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
                draw_plot();
                ImGui::EndChild();
                if (view_mode == ViewMode_Split) ImGui::SameLine();
            }
            if (view_mode == ViewMode_3D || view_mode == ViewMode_Split) {
                ImGui::BeginChild("##plot3d", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
                draw_plot3d();
                ImGui::EndChild();
            }
        }
        ImGui::End();
    }
};

static Shapespace instance = {};

}  // namespace shapespace
