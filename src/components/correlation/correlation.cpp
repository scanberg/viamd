// CORRELATION
//
// One series plotted against another over the frames of the trajectory, in one scatter plot.
//
// What is plotted is a stack of LAYERS the user defines, as in the Ramachandran window (without its reference). A layer
// pairs series dropped from the timeline (or the dataset) on its X and Y axes: one against one, one against several
// (1 x n, n x 1), or several against as many, pairwise (x_i against y_i). A series of a population (one value per
// member, as a distance per pair of atoms) is taken member by member, paired by index with the other side's members, or
// against the other side's single value. A layer is over some of the frames: the current one, all of them, or those of
// a timeline period; and drawn as points, as contours of their density (at levels of the fraction of the samples they
// enclose), or as a heat map of it. The layers are drawn bottom up, the first one on top.
//
// The samples of a layer are gathered on the main thread when what they are of changes: the X series is sampled at its
// own samples within the layer's frames, and the Y series at the same time (interpolated when it is sampled elsewhere).
// Their Pearson correlation, and the least squares line through them, are found pair by pair, and their Spearman rank
// correlation (of an evenly spaced subsample past a million). A periodic series (an angle) is cut at its seam, and its
// correlation is the circular one: Jammalamadaka-SenGupta's against another angle, Mardia's (unsigned) against a linear
// value. Distributions (script properties over the frames) are correlated with each other bin by bin, over their common
// range, of the layer's frames. Their density is
// estimated off the main thread (density_map.h), over a bounded domain of their own extent, smoothed with a gaussian of
// a fraction of it.
//
//   drag a series          from the timeline or the distributions onto an axis of the plot (or a list of the
//                          inspector): it is added to the selected layer
//   hover                  says what the point under the pointer is, and shows it in the viewport
//   click                  on a point over time, goes to its frame
//   double-click           fits the view
//   a period dropped       from the timeline on the plot, a new layer of the selected one's series over it; on a
//                          layer, its frames

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <md_attributes.h>

#include <viamd_event.h>
#include <viamd.h>
#include <serialization_utils.h>
#include <timeline_periods_ui.h>
#include <property_grid.h>
#include <gfx/gl.h>
#include <task_system.h>
#include <density_map.h>
#include <layer_groups.h>
#include <plot_series.h>
#include <plot_periodic_core.h>

#include "correlation_core.h"

#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <imgui_internal.h>
#include <implot_internal.h>

#include <algorithm>
#include <atomic>
#include <memory>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdlib.h>

namespace correlation {

// ## Constants

constexpr int MAX_LAYERS        = 16;
constexpr int MAX_AXIS_SERIES   = 16;       // On each axis of a layer
constexpr int MAX_PAIRS         = MAX_AXIS_SERIES;
constexpr int MAX_LEVELS        = DENSITY_MAP_MAX_LEVELS;
constexpr size_t MAX_POINTS     = 100000;   // Of a layer of points: samples are skipped beyond
constexpr size_t MAX_RANKED     = 1000000;  // Of a pair, ranked for its Spearman correlation: samples are skipped beyond
constexpr int    NUM_BINS       = 256;      // Asked of a distribution: it may give fewer
constexpr uint32_t NO_FRAME     = UINT32_MAX;   // Of a sample of a distribution: a bin, not a frame
constexpr int    MAX_LAG        = 1000;     // Frames, either way, of the correlation against the lag
constexpr double LAG_BUDGET     = 2e7;      // Samples times lags, of a pair: members are skipped beyond

// Colours for new layers, the least used first
static const vec4_t layer_palette[] = {
    {0.95f, 0.95f, 0.95f, 1.0f},
    {1.00f, 0.62f, 0.20f, 1.0f},
    {0.35f, 0.80f, 1.00f, 1.0f},
    {0.95f, 0.42f, 0.80f, 1.0f},
    {0.55f, 0.90f, 0.35f, 1.0f},
    {1.00f, 0.88f, 0.30f, 1.0f},
    {0.68f, 0.58f, 1.00f, 1.0f},
    {1.00f, 0.42f, 0.42f, 1.0f},
};

static const ImVec4 error_color = {1.00f, 0.43f, 0.43f, 1.0f};

// ## Layers

enum LayerFrames : int {
    LayerFrames_Current = 0,    // The frame shown
    LayerFrames_All,            // The whole trajectory
    LayerFrames_Period,         // The frames of a timeline period, where they are now
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
    PointColoring_Pair,         // Each pair (or member) its own
    PointColoring_Time,         // The time, through the layer's colour map
    PointColoring_Count,
};

static const char* style_names[LayerStyle_Count] = { "Points", "Contours", "Heatmap" };
static const char* point_coloring_names[PointColoring_Count] = { "Layer color", "Pair", "Time" };

// What a series on an axis of a layer resolved to, the last time its samples were gathered
struct SlotInfo {
    bool  ok = false;
    bool  binned = false;           // A distribution as evaluated: its bins are its samples
    int   dim = 1;
    SeriesPeriodic periodic = {};   // Of periodic values, in the display unit
    char  label[96] = "";
    char  unit[32] = "";
    char  bin_unit[32] = "";        // Of the bins of a distribution
    char  script_ident[64] = "";
    md_script_vis_ref_t vis_ref = {};
    ImVec4 color = {1, 1, 1, 1};
};

// The samples of one pair, pooled over its members: what its correlation and line are found from
struct PairStats {
    uint64_t n = 0;
    double   mx = 0, my = 0;
    double   sxx = 0, syy = 0, sxy = 0;    // Sums of the products of the deviations
    float    x_min = 0, x_max = 0;
    int      members = 0;

    inline void add(double x, double y) {
        if (n == 0) {
            x_min = x_max = (float)x;
        } else {
            x_min = MIN(x_min, (float)x);
            x_max = MAX(x_max, (float)x);
        }
        n += 1;
        const double dx = x - mx;
        const double dy = y - my;
        mx += dx / (double)n;
        my += dy / (double)n;
        sxx += dx * (x - mx);
        syy += dy * (y - my);
        sxy += dx * (y - my);
    }
    // Found once the samples are in (pair_statistics)
    size_t beg = 0, end = 0;        // Its samples, in the layer's arrays
    bool   rho_valid = false;       // Spearman's, of the ranks
    double rho = 0;
    uint64_t rho_n = 0;             // The samples ranked: at most MAX_RANKED
    int    circular = 0;            // Axes of periodic values: then rc in place of Pearson's
    bool   rc_valid = false;        // One: Mardia's circular-linear (0..1). Two: Jammalamadaka-SenGupta's circular (-1..1)
    double rc = 0;

    bool   has_r()     const { return n >= 3 && sxx > 0.0 && syy > 0.0; }
    double r()         const { return has_r() ? sxy / sqrt(sxx * syy) : 0.0; }
    bool   has_corr()  const { return circular ? rc_valid : has_r(); }
    double corr()      const { return circular ? rc : r(); }
    // No line through values around a circle
    bool   has_fit()   const { return n >= 2 && sxx > 0.0 && !circular; }
    double slope()     const { return has_fit() ? sxy / sxx : 0.0; }
    double intercept() const { return my - slope() * mx; }
};

struct Layer {
    uint32_t id = 0;                // Unique while it lives, never 0
    RepresentationNode node = {};   // Where it is in the list (layer_groups.h): at the top level, there are no groups
    char     name[64] = "";         // Empty: named after what it shows (layer_label)
    bool     visible = true;

    SeriesKey xs[MAX_AXIS_SERIES] = {};
    SeriesKey ys[MAX_AXIS_SERIES] = {};
    int      num_x = 0;
    int      num_y = 0;

    LayerFrames frames = LayerFrames_All;
    uint32_t    period = 0;         // Of LayerFrames_Period, by id

    LayerStyle style = LayerStyle_Points;
    vec4_t   color = {1, 1, 1, 1};
    float    opacity = 1.0f;

    // Contours
    bool     fill = false;
    bool     lines = true;
    float    line_width = 1.5f;            // Pixels
    int      num_levels = 3;
    float    levels[MAX_LEVELS] = { 0.5f, 0.8f, 0.95f };   // Fractions of the samples enclosed, increasing
    float    bandwidth = 3.0f;             // % of the extent of the samples along each axis: the sigma of the smoothing

    // Heatmap (and points by time)
    ImPlotColormap colormap = ImPlotColormap_Viridis;
    bool     log_scale = false;

    // Points
    PointColoring point_coloring = PointColoring_Pair;
    float    point_size = 4.0f;            // Pixels across
    bool     show_fit = false;             // The least squares line of each pair

    // ### Derived, not saved

    // The samples, gathered when what they are of changes
    md_array(float)    sx = 0;
    md_array(float)    sy = 0;
    md_array(float)    st = 0;      // Time
    md_array(uint32_t) sframe = 0;  // The frame nearest
    md_array(uint16_t) spair = 0;
    md_array(uint16_t) smember = 0;
    uint64_t data_key = 0;          // Of what they were gathered from
    uint64_t data_hash = 0;         // Of the samples themselves
    bool     data_valid = false;
    const char* data_error = "";    // Why there are none
    int      num_pairs = 0;
    PairStats pairs[MAX_PAIRS] = {};
    SlotInfo x_info[MAX_AXIS_SERIES] = {};
    SlotInfo y_info[MAX_AXIS_SERIES] = {};
    float    x_min = 0, x_max = 0, y_min = 0, y_max = 0;
    float    t_min = 0, t_max = 0;
    bool     members_truncated = false;     // Populations of different sizes paired: the fewer members of the two
    bool     binned = false;                // Of distributions: bin by bin, the samples have no frame (NO_FRAME)

    // The density, on the GPU, over its domain, and what it is of
    uint32_t tex = 0;
    DensityMapResult result = {};
    DensityMapDomain result_domain = {};
    uint64_t result_samples = 0;
    bool     has_result = false;
    uint64_t shown_hash = 0;        // Of the inputs of the density in tex
    uint64_t want_hash = 0;         // Of the inputs it should be of, 0 for none
    DensityMapJob* job = nullptr;   // Under way
    uint64_t job_serial = 0;        // Of the job launched last: the callbacks of cancelled ones are not its

    uint32_t point_stride = 1;      // Samples skipped to keep the points within MAX_POINTS

    // The correlation of each pair against the lag of Y behind X (update_lag): num_pairs rows of 2 lag_max + 1, from
    // -lag_max, NaN where there is none
    std::vector<float> lag_r;
    int      lag_max = 0;
    uint64_t lag_key = 0;
    bool     lag_skipped = false;   // Some members were left out, to keep within LAG_BUDGET
};

// The shape of a layer's axes
enum LayerShape : int {
    LayerShape_Empty = 0,       // An axis without series
    LayerShape_Mismatch,        // Several against a different number of several
    LayerShape_Valid,
};

static LayerShape layer_shape(const Layer& l) {
    if (l.num_x == 0 || l.num_y == 0) return LayerShape_Empty;
    if (l.num_x > 1 && l.num_y > 1 && l.num_x != l.num_y) return LayerShape_Mismatch;
    return LayerShape_Valid;
}

static inline int layer_num_pairs(const Layer& l) { return layer_shape(l) == LayerShape_Valid ? MAX(l.num_x, l.num_y) : 0; }
static inline int pair_x(const Layer& l, int p) { return l.num_x == 1 ? 0 : p; }
static inline int pair_y(const Layer& l, int p) { return l.num_y == 1 ? 0 : p; }

// ## Series keys in the workspace: "source;variant;path"

static const char* variant_names[SeriesVariant_Count] = { "values", "mean", "sigma", "extent", "aggregate" };

static void key_to_str(char* buf, size_t cap, const SeriesKey& k) {
    snprintf(buf, cap, "%s;%s;%s", series_source_name(k.source), variant_names[k.variant < SeriesVariant_Count ? k.variant : 0], k.path);
}

static bool key_from_str(SeriesKey* out, str_t s) {
    SeriesKey k = {};
    str_t parts[2];
    for (int i = 0; i < 2; ++i) {
        size_t sep = 0;
        if (!str_find_char(&sep, s, ';')) return false;
        parts[i] = str_substr(s, 0, sep);
        s = str_substr(s, sep + 1);
    }
    if (!series_source_from_name(&k.source, parts[0], nullptr)) return false;
    bool found = false;
    for (int v = 0; v < SeriesVariant_Count; ++v) {
        if (str_eq_cstr(parts[1], variant_names[v])) {
            k.variant = (SeriesVariant)v;
            found = true;
        }
    }
    if (!found || str_empty(s)) return false;
    str_copy_to_char_buf(k.path, sizeof(k.path), s);
    *out = k;
    return true;
}

// The value of member k of a series at fractional sample index fi, in its display unit. Of periodic values, the short way
// around the circle between the samples, and possibly outside the period they are given in (cut_periodic).
static inline double sample_at(const SeriesTemporalView& v, double fi, int k) {
    const int n = v.num_samples;
    const int i0 = CLAMP((int)floor(fi), 0, n - 1);
    const int i1 = MIN(i0 + 1, n - 1);
    const double t = fi - (double)i0;
    const double a = (double)v.y[(size_t)i0 * v.stride + k] * v.y_scale;
    double b = (double)v.y[(size_t)i1 * v.stride + k] * v.y_scale;
    if (t == 0.0) return a;
    if (series_periodic_valid(v.periodic)) b = periodic_near(b, a, series_periodic_length(v.periodic));
    return a + (b - a) * t;
}

// A periodic value laid out where the circle is cut: at the seam of its series, where its values are least likely to be
static inline double cut_periodic(double x, const SeriesPeriodic& p) {
    return series_periodic_valid(p) ? periodic_wrap(x, p.seam, series_periodic_length(p)) : x;
}

// A value as it is correlated: of periodic values, its angle in radians from the beginning of the period
static inline double to_corr_value(double v, const SeriesPeriodic& p) {
    return series_periodic_valid(p) ? 2.0 * M_PI * (v - p.beg) / series_periodic_length(p) : v;
}

// The value of member k of a distribution at x, in its display unit: between the centres of its bins
static inline double bin_value_at(const SeriesHistogramView& v, double x, int k) {
    const double fi = series_histogram_index_at(v, x);
    const int n = v.num_bins;
    const int i0 = CLAMP((int)floor(fi), 0, n - 1);
    const int i1 = CLAMP(i0 + 1, 0, n - 1);
    const double t = CLAMP(fi - (double)i0, 0.0, 1.0);
    const double a = v.bins[(size_t)k * n + i0];
    const double b = v.bins[(size_t)k * n + i1];
    return (a + (b - a) * t) * v.y_scale;
}

// A key as it is plotted here: the values of a population pooled over its members are its values (each member is a
// sample anyway), and the frames are the layer's (the period of a key goes to its layer, see drop_series)
static SeriesKey plotted_key(const SeriesKey& k) {
    SeriesKey out = k;
    if (out.variant == SeriesVariant_Aggregate) out.variant = SeriesVariant_Values;
    out.period = 0;
    return out;
}

// Whether key is of a script distribution, evaluated as such (not a property over time binned)
static bool is_script_distribution(const ApplicationState* app, const SeriesKey& k) {
    if (k.source != SeriesSource_Script || !app->script.eval_ir) return false;
    const md_script_property_flags_t f = md_script_ir_property_flags(app->script.eval_ir, series_script_ident(k));
    return (f & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) && !(f & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL);
}

// The time of a fractional frame, along the timeline
static double frame_time(const ApplicationState* app, double frame) {
    const int64_t n = (int64_t)md_array_size(app->timeline.x_values);
    if (n <= 0) return frame;
    const int64_t f0 = CLAMP((int64_t)frame, (int64_t)0, n - 1);
    const int64_t f1 = CLAMP((int64_t)frame + 1, (int64_t)0, n - 1);
    const double t = frame - floor(frame);
    return app->timeline.x_values[f0] + (app->timeline.x_values[f1] - app->timeline.x_values[f0]) * t;
}

// The frame whose time is nearest t
static uint32_t nearest_frame(const ApplicationState* app, double t) {
    const float* xs = app->timeline.x_values;
    const int64_t n = (int64_t)md_array_size(xs);
    if (n <= 0) return 0;
    if (t <= xs[0]) return 0;
    if (t >= xs[n - 1]) return (uint32_t)(n - 1);
    int64_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        const int64_t mid = (lo + hi) / 2;
        if (xs[mid] <= t) lo = mid;
        else hi = mid;
    }
    return (uint32_t)(t - xs[lo] <= xs[hi] - t ? lo : hi);
}

// ## The component

struct Correlation : viamd::EventHandler {
    ApplicationState* app = nullptr;
    bool show_window = false;

    // In the order of the list (apply_tree): what is listed first is drawn on top
    Layer layers[MAX_LAYERS];
    int   num_layers = 0;
    uint32_t next_layer_id = 1;
    uint32_t selected_layer = 0;    // By id: the one the inspector shows, and that series dropped on the plot go to
    bool  tree_dirty = false;       // The layers are to be put in the order of the list again
    static constexpr uint32_t LIST_ID = 2;  // What its rows are dragged as (LayerNodeRef)

    uint64_t system_generation = 1;  // Moves with every system loaded or freed: samples are of one of them

    // How it is shown
    bool show_sidebar = true;
    bool show_legend = true;
    bool show_lag = false;          // The correlation against the lag, under the plot
    int  max_lag = 50;              // Frames, either way
    float lag_frac = 0.3f;          // Of the height, the lag's plot
    int  lag_shown = 0;             // The max_lag its axis was fitted to
    ImPlotRect viewrect = ImPlotRect(0, 1, 0, 1);
    bool  fit_view = true;          // To the samples of the layers shown, next frame
    uint64_t fit_key = 0;           // Of the series of the layers shown: fitted again when they change

    // The highlight is whoever sets it to clear: what this window highlighted and no longer does is cleared once
    bool  owns_highlight = false;
    bool  highlights_now = false;

    DensityMapCanvas canvas = {};

    bool read_layers = false;       // While a workspace is read: whether it had any
    bool layers_declared = false;   // and whether it said how many (possibly none)

    Correlation() { viamd::event_system_register_handler(*this); }

    // ## Events

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app = (ApplicationState*)e.payload;
                create_default_layers();
                workspace_register_window("Correlation", &show_window);
                break;
            }
            case viamd::EventType_ViamdShutdown:
                shutdown();
                break;
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                clear_layers();
                reset_settings();
                read_layers = false;
                layers_declared = false;
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                const str_t section = viamd::section_header(state);
                if (str_eq(section, STR_LIT("Correlation"))) {
                    deserialize_settings(state);
                } else if (str_eq(section, STR_LIT("CorrelationLayer"))) {
                    deserialize_layer(state);
                }
                break;
            }
            case viamd::EventType_ViamdDeserializeEnd:
                if (!read_layers && !layers_declared) create_default_layers();
                fit_view = true;
                break;
            case viamd::EventType_ViamdSystemInit:
                system_generation += 1;
                break;
            case viamd::EventType_ViamdSystemFree:
            case viamd::EventType_ViamdTrajectoryFree:
                // What the samples were of is gone: nothing of it is shown again
                cancel_all_jobs(true);
                for (int i = 0; i < num_layers; ++i) {
                    layers[i].has_result = false;
                    layers[i].shown_hash = 0;
                    layers[i].want_hash = 0;
                    layers[i].data_key = 0;
                }
                system_generation += 1;
                density_map_canvas_invalidate(&canvas);
                break;
            case viamd::EventType_ViamdFrameTick:
                update();
                draw();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Correlation", &show_window);
                break;
            default:
                break;
            }
        }
    }

    void shutdown() {
        cancel_all_jobs(true);
        clear_layers();
        density_map_canvas_free(&canvas);
        density_map_shutdown();
    }

    // ## Layers

    Layer* find_layer(uint32_t id) {
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].id == id) return &layers[i];
        }
        return nullptr;
    }

    const Layer* find_layer(uint32_t id) const {
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

    vec4_t next_layer_color() const {
        int best = 0;
        int best_uses = INT32_MAX;
        for (int p = 0; p < (int)ARRAY_SIZE(layer_palette); ++p) {
            int uses = 0;
            for (int i = 0; i < num_layers; ++i) {
                if (vec4_equal(layers[i].color, layer_palette[p])) ++uses;
            }
            if (uses < best_uses) {
                best = p;
                best_uses = uses;
            }
        }
        return layer_palette[best];
    }

    // ### The list

    size_t tree_nodes(RepresentationNode* out[MAX_LAYERS]) {
        for (int i = 0; i < num_layers; ++i) out[i] = &layers[i].node;
        return (size_t)num_layers;
    }

    // The tree made canonical, and the layers put in its order (moved: what points at them is no longer valid)
    void apply_tree() {
        RepresentationNode* nodes[MAX_LAYERS];
        representation_tree_normalize(nodes, tree_nodes(nodes));
        bool moved = false;
        for (int i = 1; i < num_layers; ++i) {
            for (int j = i; j > 0 && layers[j - 1].node.order > layers[j].node.order; --j) {
                const Layer tmp = layers[j];
                layers[j] = layers[j - 1];
                layers[j - 1] = tmp;
                moved = true;
            }
        }
        if (moved) density_map_canvas_invalidate(&canvas);
        tree_dirty = false;
    }

    void tree_changed() {
        RepresentationNode* nodes[MAX_LAYERS];
        representation_tree_normalize(nodes, tree_nodes(nodes));
        tree_dirty = true;
    }

    // Adds a layer: first for 'at' 0, last for -1, or before what is at position 'at' of the list. NULL when there is no
    // room. In place in the list once the tree is applied (apply_tree).
    Layer* add_layer(LayerFrames frames, LayerStyle style, int at = -1) {
        if (num_layers >= MAX_LAYERS) return nullptr;
        Layer& l = layers[num_layers++];
        l = Layer{};
        l.id = next_layer_id++;
        l.node.id = l.id;
        if (at < 0) {
            l.node.order = LAYER_TREE_LAST;
        } else if (at == 0) {
            RepresentationNode* nodes[MAX_LAYERS];
            l.node.order = layer_tree_first_order(nodes, tree_nodes(nodes), 0);
        } else {
            l.node.order = 2 * at - 1;  // Between those at at - 1 and at, the orders being two apart
        }
        l.frames = frames;
        l.style = style;
        l.color = next_layer_color();
        if (style == LayerStyle_Contours) {
            l.fill = false;
            l.lines = true;
        }
        tree_changed();
        return &l;
    }

    // The series of another layer
    static void copy_series(Layer& dst, const Layer& src) {
        MEMCPY(dst.xs, src.xs, sizeof(dst.xs));
        MEMCPY(dst.ys, src.ys, sizeof(dst.ys));
        dst.num_x = src.num_x;
        dst.num_y = src.num_y;
        dst.data_key = 0;
    }

    void free_derived(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        cancel_job(l, true);
        if (l.tex) density_map_texture_free(&l.tex);
        l.tex = 0;
        md_array_free(l.sx, heap);
        md_array_free(l.sy, heap);
        md_array_free(l.st, heap);
        md_array_free(l.sframe, heap);
        md_array_free(l.spair, heap);
        md_array_free(l.smember, heap);
        l.lag_r.clear();
        l.lag_r.shrink_to_fit();
        l.lag_key = 0;
    }

    void remove_layer(uint32_t id) {
        const int idx = layer_index(id);
        if (idx < 0) return;
        free_derived(layers[idx]);
        for (int i = idx; i + 1 < num_layers; ++i) layers[i] = layers[i + 1];
        num_layers -= 1;
        layers[num_layers] = Layer{};
        if (selected_layer == id) {
            selected_layer = num_layers > 0 ? layers[MIN(idx, num_layers - 1)].id : 0;
        }
        tree_changed();
        density_map_canvas_invalidate(&canvas);
    }

    // A copy of what is set, without what is derived from it (found again), just above it
    void duplicate_layer(uint32_t id) {
        const Layer* src = find_layer(id);
        if (!src || num_layers >= MAX_LAYERS) return;
        const Layer copy = *src;
        Layer& l = layers[num_layers++];
        l = Layer{};
        l.id = next_layer_id++;
        l.node = {l.id, 0, copy.node.order - 1, false};
        MEMCPY(l.name, copy.name, sizeof(l.name));
        l.visible = copy.visible;
        copy_series(l, copy);
        l.frames = copy.frames;
        l.period = copy.period;
        l.style = copy.style;
        l.color = copy.color;
        l.opacity = copy.opacity;
        l.fill = copy.fill;
        l.lines = copy.lines;
        l.line_width = copy.line_width;
        l.num_levels = copy.num_levels;
        MEMCPY(l.levels, copy.levels, sizeof(l.levels));
        l.bandwidth = copy.bandwidth;
        l.colormap = copy.colormap;
        l.log_scale = copy.log_scale;
        l.point_coloring = copy.point_coloring;
        l.point_size = copy.point_size;
        l.show_fit = copy.show_fit;
        selected_layer = l.id;
        tree_changed();
    }

    void clear_layers() {
        for (int i = 0; i < num_layers; ++i) {
            free_derived(layers[i]);
            layers[i] = Layer{};
        }
        num_layers = 0;
        selected_layer = 0;
        tree_dirty = false;
        density_map_canvas_invalidate(&canvas);
    }

    void move_node(LayerNodeRef src, LayerNodeRef target, RepresentationDrop drop) {
        RepresentationNode* nodes[MAX_LAYERS];
        const size_t n = tree_nodes(nodes);
        const int si = representation_tree_find(nodes, n, false, src.id);
        const int ti = target.id ? representation_tree_find(nodes, n, false, target.id) : -1;    // None: the end of the list
        if (si < 0 || (target.id && ti < 0) || src.group) return;
        if (layer_tree_move(nodes, n, (size_t)si, ti, drop)) tree_dirty = true;
    }

    void reset_settings() {
        show_sidebar = true;
        show_legend = true;
        show_lag = false;
        max_lag = 50;
        lag_frac = 0.3f;
        viewrect = ImPlotRect(0, 1, 0, 1);
        fit_view = true;
    }

    // One layer of points over all frames, waiting for its series
    void create_default_layers() {
        clear_layers();
        if (Layer* l = add_layer(LayerFrames_All, LayerStyle_Points)) {
            l->color = {1, 1, 1, 1};
            selected_layer = l->id;
        }
    }

    // ### Series on the axes

    // Added to an axis of the layer, at its end. Allowed to be a shape that is not valid (yet): several on each axis,
    // different in number, while they are being put together.
    bool add_series(Layer& l, int axis, const SeriesKey& key) {
        SeriesKey* keys = axis == 0 ? l.xs : l.ys;
        int* count = axis == 0 ? &l.num_x : &l.num_y;
        if (*count >= MAX_AXIS_SERIES) return false;
        keys[(*count)++] = key;
        l.data_key = 0;
        fit_view = true;
        return true;
    }

    void remove_series(Layer& l, int axis, int idx) {
        SeriesKey* keys = axis == 0 ? l.xs : l.ys;
        int* count = axis == 0 ? &l.num_x : &l.num_y;
        if (idx < 0 || idx >= *count) return;
        for (int i = idx; i + 1 < *count; ++i) keys[i] = keys[i + 1];
        *count -= 1;
        keys[*count] = SeriesKey{};
        l.data_key = 0;
        fit_view = true;
    }

    void swap_axes(Layer& l) {
        for (int i = 0; i < MAX_AXIS_SERIES; ++i) {
            const SeriesKey tmp = l.xs[i];
            l.xs[i] = l.ys[i];
            l.ys[i] = tmp;
        }
        const int n = l.num_x;
        l.num_x = l.num_y;
        l.num_y = n;
        l.data_key = 0;
        fit_view = true;
    }

    // A series dropped on an axis of the plot: to the selected layer, or a new one when there is none. Taken over a period,
    // it brings the period to a layer that has no series yet and is over every frame.
    void drop_series(int axis, const SeriesKey& key) {
        Layer* l = find_layer(selected_layer);
        if (!l) {
            l = add_layer(LayerFrames_All, LayerStyle_Points, 0);
            if (!l) return;
            selected_layer = l->id;
        }
        drop_series_on(*l, axis, key);
    }

    void drop_series_on(Layer& l, int axis, const SeriesKey& key) {
        if (key.period != 0 && l.num_x + l.num_y == 0 && l.frames == LayerFrames_All && timeline_periods_find(&app->timeline.periods, key.period)) {
            l.frames = LayerFrames_Period;
            l.period = key.period;
        }
        add_series(l, axis, plotted_key(key));
    }

    // A series can be plotted here when it can be over time, or is a script distribution (correlated bin by bin with
    // another), and is not a range of values
    bool series_droppable(const SeriesKey& dropped, const char** reason) const {
        *reason = nullptr;
        if (dropped.variant == SeriesVariant_Extent) {
            *reason = "Its range: drop its values or mean instead";
            return false;
        }
        const SeriesKey key = plotted_key(dropped);
        if (is_script_distribution(app, key)) return true;
        if (!(series_views(app, key, reason) & SeriesView_Timeline)) {
            if (!*reason) *reason = "Not over time, and not a distribution";
            return false;
        }
        return true;
    }

    // A payload of a series, from the timeline or the distributions
    static const SeriesDragPayload* series_payload(const ImGuiPayload* p) {
        if (!p || p->DataSize != (int)sizeof(SeriesDragPayload)) return nullptr;
        if (!p->IsDataType(TIMELINE_SERIES_DND) && !p->IsDataType(DISTRIBUTION_SERIES_DND)) return nullptr;
        return (const SeriesDragPayload*)p->Data;
    }

    // Accepts a series of either kind on the drop target open: the payload once it is delivered
    static const SeriesDragPayload* accept_series(ImGuiDragDropFlags flags = 0) {
        const ImGuiPayload* p = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND, flags);
        if (!p) p = ImGui::AcceptDragDropPayload(DISTRIBUTION_SERIES_DND, flags);
        return series_payload(p);
    }

    // ## Jobs

    void cancel_job(Layer& l, bool wait) {
        if (!l.job) return;
        // Its callback still comes, and finds the layer gone or the job not its own
        density_map_cancel(l.job, wait);
        l.job = nullptr;
        l.want_hash = 0;
    }

    void cancel_all_jobs(bool wait) {
        for (int i = 0; i < num_layers; ++i) {
            cancel_job(layers[i], wait);
        }
    }

    // The domain of a layer's density: its samples' extent, with room for the smoothing around them
    static DensityMapDomain layer_domain(const Layer& l) {
        DensityMapDomain d = {};
        d.periodic = false;
        const float bw = l.bandwidth * 0.01f;
        auto range = [bw](float lo, float hi, float* out_lo, float* out_hi) {
            float r = hi - lo;
            if (!(r > 0.0f)) r = MAX(fabsf(lo), 1.0f) * 0.1f;
            const float pad = r * (4.0f * bw + 0.02f);
            *out_lo = lo - pad;
            *out_hi = hi + pad;
        };
        range(l.x_min, l.x_max, &d.x_min, &d.x_max);
        range(l.y_min, l.y_max, &d.y_min, &d.y_max);
        return d;
    }

    // The density of the samples of the layer, pooled over its pairs, in channel 0
    void launch_job(Layer& l, uint64_t hash) {
        ASSERT(!l.job);
        struct Input {
            std::vector<float> x, y;
        };
        std::shared_ptr<Input> in = std::make_shared<Input>();
        const size_t n = md_array_size(l.sx);
        in->x.assign(l.sx, l.sx + n);
        in->y.assign(l.sy, l.sy + n);
        const DensityMapDomain domain = layer_domain(l);
        const float sigma = (l.x_max > l.x_min ? l.x_max - l.x_min : (domain.x_max - domain.x_min) * 0.5f) * l.bandwidth * 0.01f;
        const uint32_t num_chunks = (uint32_t)MIN((uint64_t)8, MAX((uint64_t)1, (uint64_t)(n >> 18)));

        auto source = [in](uint32_t chunk, uint32_t chunks, DensityMapBins& bins, const std::atomic<int32_t>& cancel) {
            const size_t count = in->x.size();
            const size_t beg = count * chunk / chunks;
            const size_t end = count * (chunk + 1) / chunks;
            for (size_t i = beg; i < end; ++i) {
                if ((i & 4095) == 0 && cancel.load(std::memory_order_relaxed)) return;
                bins.add(0, in->x[i], in->y[i]);
            }
        };
        const uint32_t id = l.id;
        const uint64_t serial = ++l.job_serial;
        auto done = [this, id, hash, serial, domain](const DensityMapResult* result, const float* grid) {
            Layer* layer = find_layer(id);
            if (!layer || !layer->job || layer->job_serial != serial || layer->want_hash != hash) return;
            layer->job = nullptr;
            if (!result) {
                // Interrupted: tried again
                layer->want_hash = 0;
                return;
            }
            if (!layer->tex) layer->tex = density_map_texture_create(false);
            density_map_texture_upload(layer->tex, grid);
            layer->result = *result;
            layer->result_domain = domain;
            layer->result_samples = result->samples[0];
            layer->has_result = true;
            layer->shown_hash = hash;
        };
        l.want_hash = hash;
        l.job = density_map_launch(STR_LIT("Correlation density"), domain, sigma, num_chunks, source, done);
    }

    // ## Update

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

    // Drawn as a density: over frames, of a trajectory; or the bins of distributions
    bool layer_draws_density(const Layer& l) const {
        if (l.style == LayerStyle_Points) return false;
        if (l.binned) return true;
        if (l.frames == LayerFrames_Current) return false;
        return run_num_frames(app) >= 2;
    }

    // Drawn as points: those of the frame shown are, whatever the style (a density of one point per pair says nothing)
    bool layer_draws_points(const Layer& l) const {
        if (l.binned) return l.style == LayerStyle_Points;
        return l.style == LayerStyle_Points || l.frames == LayerFrames_Current || run_num_frames(app) < 2;
    }

    static uint64_t series_version(const ApplicationState* app, const SeriesKey& k) {
        const str_t path = str_from_cstr(k.path);
        const md_attributes_t* t = series_table(app, k.source, path);
        const md_attribute_t* a = t ? md_attributes_find(t, path) : nullptr;
        return a ? a->version : 0;
    }

    // What a series on an axis is, and what of it its samples depend on: over time into v, or a distribution (over the
    // layer's period, 0 for every frame) into hv
    uint64_t resolve_slot(SeriesTemporalView* v, SeriesHistogramView* hv, SlotInfo* info, const SeriesKey& key, uint32_t period) {
        *info = SlotInfo{};
        const SeriesKey k = plotted_key(key);
        uint64_t h = md_hash64(&k.source, sizeof(k.source), 0);
        h = md_hash64(&k.variant, sizeof(k.variant), h);
        h = md_hash64(k.path, strlen(k.path), h);
        h = md_hash64_combine(h, series_version(app, k));
        const SeriesState st = series_state(app, k);
        h = md_hash64_combine(h, ((uint64_t)st.num_frames << 8) | st.status);
        info->color = series_default_color(app, k);

        if (is_script_distribution(app, k)) {
            info->binned = true;
            info->ok = series_resolve_histogram(hv, app, k, NUM_BINS, period) && hv->bins && hv->num_bins >= 2 && hv->x_max > hv->x_min;
            info->dim = MAX(1, hv->dim);
            snprintf(info->label, sizeof(info->label), "%s", hv->label);
            snprintf(info->unit, sizeof(info->unit), "%s", hv->y_unit_str);
            snprintf(info->bin_unit, sizeof(info->bin_unit), "%s", hv->x_unit_str);
            snprintf(info->script_ident, sizeof(info->script_ident), "%s", hv->script_ident);
            info->vis_ref = hv->vis_ref;
            const uintptr_t ptr = (uintptr_t)hv->bins;
            h = md_hash64(&ptr, sizeof(ptr), h);
            const double dims[6] = { (double)hv->num_bins, (double)hv->dim, hv->x_min, hv->x_max, hv->y_scale, (double)info->ok };
            h = md_hash64(dims, sizeof(dims), h);
            h = md_hash64(&period, sizeof(period), h);
            return h;
        }

        // A band's centre (Sigma: its mean) is its value here; a range has none
        info->ok = series_resolve_temporal(v, app, k) && k.variant != SeriesVariant_Extent && v->x && v->y && v->stride >= 1;
        info->dim = MAX(1, v->dim);
        info->periodic = v->periodic;
        snprintf(info->label, sizeof(info->label), "%s", v->label);
        snprintf(info->unit, sizeof(info->unit), "%s", v->unit_str);
        snprintf(info->script_ident, sizeof(info->script_ident), "%s", v->script_ident);
        info->vis_ref = v->vis_ref;
        const uintptr_t ptrs[2] = { (uintptr_t)v->x, (uintptr_t)v->y };
        h = md_hash64(ptrs, sizeof(ptrs), h);
        const int64_t dims[3] = { v->num_samples, v->stride, info->ok };
        h = md_hash64(dims, sizeof(dims), h);
        h = md_hash64(&v->y_scale, sizeof(v->y_scale), h);
        return h;
    }

    // The samples of the layer, gathered again when what they are of changes
    void gather(Layer& l) {
        static SeriesTemporalView vx[MAX_AXIS_SERIES];
        static SeriesTemporalView vy[MAX_AXIS_SERIES];
        static SeriesHistogramView hx[MAX_AXIS_SERIES];
        static SeriesHistogramView hy[MAX_AXIS_SERIES];

        const uint32_t dist_period = l.frames == LayerFrames_Period ? l.period : 0;
        uint64_t key = md_hash64_combine(system_generation, app->timeline.units_version);
        for (int i = 0; i < l.num_x; ++i) key = md_hash64_combine(key, resolve_slot(&vx[i], &hx[i], &l.x_info[i], l.xs[i], dist_period) ^ 0x1111);
        for (int i = 0; i < l.num_y; ++i) key = md_hash64_combine(key, resolve_slot(&vy[i], &hy[i], &l.y_info[i], l.ys[i], dist_period) ^ 0x2222);
        const int num_frames = (int)run_num_frames(app);
        uint32_t beg = 0, end = 0;
        const bool has_range = layer_frame_range(l, &beg, &end);
        // Distributions are of every frame unless of a period: the frame shown does not change them
        bool all_binned = l.num_x + l.num_y > 0;
        for (int i = 0; i < l.num_x; ++i) all_binned &= l.x_info[i].binned;
        for (int i = 0; i < l.num_y; ++i) all_binned &= l.y_info[i].binned;
        const bool by_frame = l.frames == LayerFrames_Current && !all_binned;
        const uint64_t frames[5] = { (uint64_t)l.frames, by_frame || l.frames != LayerFrames_Current ? beg : 0, by_frame || l.frames != LayerFrames_Current ? end : 0, (uint64_t)has_range, (uint64_t)num_frames };
        key = md_hash64(frames, sizeof(frames), key);
        if (by_frame) key = md_hash64(&app->animation.frame, sizeof(app->animation.frame), key);
        key |= 1;
        if (key == l.data_key) return;
        l.data_key = key;

        md_array_shrink(l.sx, 0);
        md_array_shrink(l.sy, 0);
        md_array_shrink(l.st, 0);
        md_array_shrink(l.sframe, 0);
        md_array_shrink(l.spair, 0);
        md_array_shrink(l.smember, 0);
        l.data_valid = false;
        l.data_hash = 0;
        l.num_pairs = 0;
        l.members_truncated = false;
        l.binned = false;
        l.point_stride = 1;
        for (int p = 0; p < MAX_PAIRS; ++p) l.pairs[p] = PairStats{};

        const LayerShape shape = layer_shape(l);
        if (shape == LayerShape_Empty) {
            l.data_error = "drop series on the X and Y axes";
            return;
        }
        if (shape == LayerShape_Mismatch) {
            l.data_error = "pair as many X as Y, or one against several";
            return;
        }
        int num_binned = 0;
        for (int i = 0; i < l.num_x; ++i) num_binned += l.x_info[i].binned;
        for (int i = 0; i < l.num_y; ++i) num_binned += l.y_info[i].binned;
        if (num_binned > 0 && num_binned < l.num_x + l.num_y) {
            l.data_error = "a distribution pairs only with distributions";
            return;
        }
        l.binned = num_binned > 0;
        if (!l.binned && (!has_range || num_frames <= 0)) {
            l.data_error = l.frames == LayerFrames_Period ? "the period is missing" : "no frames";
            return;
        }
        l.num_pairs = layer_num_pairs(l);
        if (l.binned) {
            gather_binned(l, hx, hy);
            finish_gather(l, nullptr);
            return;
        }

        const float* run_x = app->timeline.x_values;
        const bool bounded = l.frames == LayerFrames_Period;
        const double t_beg = bounded ? run_x[beg] : -HUGE_VAL;
        const double t_end = bounded ? run_x[end - 1] : HUGE_VAL;
        const double t_eps = bounded ? 1e-6 * MAX(fabs(t_end - t_beg), 1e-6) : 0.0;
        const bool current = l.frames == LayerFrames_Current;
        const double t_cur = frame_time(app, app->animation.frame);
        const uint32_t frame_cur = (uint32_t)CLAMP((int)(app->animation.frame + 0.5), 0, num_frames - 1);

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_allocator_i* alloc = md_temp_allocator(temp);

        bool any_resolved = false;
        for (int p = 0; p < l.num_pairs; ++p) {
            const int xi = pair_x(l, p);
            const int yi = pair_y(l, p);
            if (!l.x_info[xi].ok || !l.y_info[yi].ok) continue;
            any_resolved = true;
            const SeriesTemporalView& X = vx[xi];
            const SeriesTemporalView& Y = vy[yi];
            if (X.num_samples <= 0 || Y.num_samples <= 0) continue;
            const int dx = l.x_info[xi].dim;
            const int dy = l.y_info[yi].dim;
            const int members = dx == 1 ? dy : dy == 1 ? dx : MIN(dx, dy);
            if (dx > 1 && dy > 1 && dx != dy) l.members_truncated = true;
            l.pairs[p].members = members;

            md_bitset_t mask_x = {}, mask_y = {};
            const bool has_mask_x = series_frame_mask(&mask_x, app, l.xs[xi], alloc);
            const bool has_mask_y = series_frame_mask(&mask_y, app, l.ys[yi], alloc);
            const bool same_axis = X.x == Y.x && X.num_samples == Y.num_samples;
            const bool x_on_frames = X.x == run_x && (int)md_array_size(run_x) == X.num_samples;
            const double y_lo = Y.x[0], y_hi = Y.x[Y.num_samples - 1];

            // Whether Y has a value at fractional sample index fy: both samples around it evaluated
            auto y_held = [&](double fy) {
                if (!has_mask_y) return true;
                const uint32_t i0 = (uint32_t)floor(fy);
                const uint32_t i1 = MIN(i0 + 1, (uint32_t)Y.num_samples - 1);
                return md_bitset_test(mask_y, i0) && (fy == (double)i0 || md_bitset_test(mask_y, i1));
            };
            auto x_held = [&](double fx) {
                if (!has_mask_x) return true;
                const uint32_t i0 = (uint32_t)floor(fx);
                const uint32_t i1 = MIN(i0 + 1, (uint32_t)X.num_samples - 1);
                return md_bitset_test(mask_x, i0) && (fx == (double)i0 || md_bitset_test(mask_x, i1));
            };
            auto push = [&](double x, double y, double t, uint32_t frame, int m) {
                push_sample(l, p, x, y, t, frame, m);
            };

            if (current) {
                if (t_cur < X.x[0] || t_cur > X.x[X.num_samples - 1] || t_cur < y_lo || t_cur > y_hi) continue;
                const double fx = series_temporal_index_at(X, t_cur);
                const double fy = series_temporal_index_at(Y, t_cur);
                if (!x_held(fx) || !y_held(fy)) continue;
                for (int m = 0; m < members; ++m) {
                    push(sample_at(X, fx, dx == 1 ? 0 : m), sample_at(Y, fy, dy == 1 ? 0 : m), t_cur, frame_cur, m);
                }
                continue;
            }

            for (int i = 0; i < X.num_samples; ++i) {
                const double t = X.x[i];
                if (t < t_beg - t_eps || t > t_end + t_eps) continue;
                if (has_mask_x && !md_bitset_test(mask_x, (uint32_t)i)) continue;
                double fy = (double)i;
                if (!same_axis) {
                    if (t < y_lo || t > y_hi) continue;
                    fy = series_temporal_index_at(Y, t);
                }
                if (!y_held(fy)) continue;
                const uint32_t frame = x_on_frames ? (uint32_t)i : nearest_frame(app, t);
                for (int m = 0; m < members; ++m) {
                    const double x = (double)X.y[(size_t)i * X.stride + (dx == 1 ? 0 : m)] * X.y_scale;
                    push(x, sample_at(Y, fy, dy == 1 ? 0 : m), t, frame, m);
                }
            }
        }

        finish_gather(l, any_resolved ? "no samples over its frames" : "the series are not available");
    }

    // A sample of pair p, its periodic values laid out where their circle is cut
    static void push_sample(Layer& l, int p, double x, double y, double t, uint32_t frame, int m) {
        if (!isfinite(x) || !isfinite(y)) return;
        x = cut_periodic(x, l.x_info[pair_x(l, p)].periodic);
        y = cut_periodic(y, l.y_info[pair_y(l, p)].periodic);
        md_allocator_i* heap = md_get_heap_allocator();
        md_array_push(l.sx, (float)x, heap);
        md_array_push(l.sy, (float)y, heap);
        md_array_push(l.st, (float)t, heap);
        md_array_push(l.sframe, frame, heap);
        md_array_push(l.spair, (uint16_t)p, heap);
        md_array_push(l.smember, (uint16_t)m, heap);
        l.pairs[p].add(x, y);
    }

    // Distributions, bin by bin over the range they share: each pair at as many places as the finer of the two has bins,
    // its values there interpolated between the centres of the bins. What a sample is at is its place (in st).
    void gather_binned(Layer& l, const SeriesHistogramView* hx, const SeriesHistogramView* hy) {
        for (int p = 0; p < l.num_pairs; ++p) {
            const int xi = pair_x(l, p);
            const int yi = pair_y(l, p);
            if (!l.x_info[xi].ok || !l.y_info[yi].ok) continue;
            const SeriesHistogramView& X = hx[xi];
            const SeriesHistogramView& Y = hy[yi];
            const double lo = MAX(X.x_min, Y.x_min);
            const double hi = MIN(X.x_max, Y.x_max);
            if (!(hi > lo)) continue;
            const int dx = l.x_info[xi].dim;
            const int dy = l.y_info[yi].dim;
            const int members = dx == 1 ? dy : dy == 1 ? dx : MIN(dx, dy);
            if (dx > 1 && dy > 1 && dx != dy) l.members_truncated = true;
            l.pairs[p].members = members;
            const int n = MAX(X.num_bins, Y.num_bins);
            for (int j = 0; j < n; ++j) {
                const double at = lo + (j + 0.5) * (hi - lo) / n;
                for (int m = 0; m < members; ++m) {
                    push_sample(l, p, bin_value_at(X, at, dx == 1 ? 0 : m), bin_value_at(Y, at, dy == 1 ? 0 : m), at, NO_FRAME, m);
                }
            }
        }
        bool any = false;
        for (int p = 0; p < l.num_pairs; ++p) any |= l.x_info[pair_x(l, p)].ok && l.y_info[pair_y(l, p)].ok;
        l.data_error = any ? "the distributions share no range of values" : "the distributions are not available";
    }

    // Once the samples are in: their extent, their statistics, and how many of them are drawn
    void finish_gather(Layer& l, const char* empty_error) {
        const size_t n = md_array_size(l.sx);
        if (n == 0) {
            if (!l.binned) l.data_error = empty_error;
            return;
        }
        l.x_min = l.x_max = l.sx[0];
        l.y_min = l.y_max = l.sy[0];
        l.t_min = l.t_max = l.st[0];
        for (size_t i = 1; i < n; ++i) {
            l.x_min = MIN(l.x_min, l.sx[i]);
            l.x_max = MAX(l.x_max, l.sx[i]);
            l.y_min = MIN(l.y_min, l.sy[i]);
            l.y_max = MAX(l.y_max, l.sy[i]);
            l.t_min = MIN(l.t_min, l.st[i]);
            l.t_max = MAX(l.t_max, l.st[i]);
        }
        // The samples of a pair are together, in the order of the pairs
        for (int p = 0; p < l.num_pairs; ++p) l.pairs[p].beg = l.pairs[p].end = 0;
        for (size_t i = 0; i < n;) {
            const int p = l.spair[i];
            size_t j = i;
            while (j < n && l.spair[j] == p) ++j;
            l.pairs[p].beg = i;
            l.pairs[p].end = j;
            i = j;
        }
        for (int p = 0; p < l.num_pairs; ++p) pair_statistics(l, p);

        uint64_t h = md_hash64(l.sx, md_array_bytes(l.sx), n);
        h = md_hash64(l.sy, md_array_bytes(l.sy), h);
        h = md_hash64(l.sframe, md_array_bytes(l.sframe), h);
        l.data_hash = h | 1;
        l.data_valid = true;
        l.data_error = "";
        l.point_stride = (uint32_t)MAX((size_t)1, (n + MAX_POINTS - 1) / MAX_POINTS);
    }

    // Spearman's correlation of a pair, of at most MAX_RANKED of its samples; of periodic values, the circular one
    static void pair_statistics(Layer& l, int p) {
        PairStats& s = l.pairs[p];
        const size_t count = s.end - s.beg;
        if (count < 3) return;
        const SeriesPeriodic& px = l.x_info[pair_x(l, p)].periodic;
        const SeriesPeriodic& py = l.y_info[pair_y(l, p)].periodic;
        const bool cx = series_periodic_valid(px);
        const bool cy = series_periodic_valid(py);
        s.circular = (int)cx + (int)cy;

        if (!s.circular) {
            const size_t step = MAX((size_t)1, (count + MAX_RANKED - 1) / MAX_RANKED);
            const size_t m = (count + step - 1) / step;
            std::vector<double> x(m), y(m);
            for (size_t i = 0; i < m; ++i) {
                x[i] = l.sx[s.beg + i * step];
                y[i] = l.sy[s.beg + i * step];
            }
            s.rho_valid = corr_spearman(&s.rho, x.data(), y.data(), m);
            s.rho_n = m;
            return;
        }

        std::vector<double> x(count), y(count);
        for (size_t i = 0; i < count; ++i) {
            x[i] = to_corr_value(l.sx[s.beg + i], px);
            y[i] = to_corr_value(l.sy[s.beg + i], py);
        }
        s.rc_valid = corr_of(&s.rc, corr_kind(cx, cy), x.data(), y.data(), count);
    }

    // Has a correlation against the lag: samples at frames, of more than one
    bool layer_has_lag(const Layer& l) const {
        return l.data_valid && !l.binned && l.frames != LayerFrames_Current && run_num_frames(app) >= 2;
    }

    // The correlation of each pair of X at frame f against Y at frame f + lag, pooled over members (paired with
    // themselves), for lags of up to max_lag frames either way
    void update_lag(Layer& l) {
        if (!show_lag || !l.visible || !layer_has_lag(l)) return;
        const int L = CLAMP(max_lag, 1, MAX_LAG);
        const uint64_t key = md_hash64(&L, sizeof(L), l.data_hash) | 1;
        if (key == l.lag_key) return;
        l.lag_key = key;
        l.lag_max = L;
        l.lag_skipped = false;
        const int nl = 2 * L + 1;
        l.lag_r.assign((size_t)l.num_pairs * nl, NAN);

        std::vector<uint32_t> ord, seg, ia, ib;
        std::vector<uint32_t> fr;
        std::vector<double> xv, yv, a, b;
        for (int p = 0; p < l.num_pairs; ++p) {
            const PairStats& s = l.pairs[p];
            const size_t count = s.end - s.beg;
            if (count < 3) continue;
            const SeriesPeriodic& px = l.x_info[pair_x(l, p)].periodic;
            const SeriesPeriodic& py = l.y_info[pair_y(l, p)].periodic;
            const CorrKind kind = corr_kind(series_periodic_valid(px), series_periodic_valid(py));

            // The samples of each member together, in the order of their frames (as they were gathered): every
            // m_step'th member, to keep within the budget
            const int members = MAX(1, s.members);
            const int m_step = (int)MIN((double)members, MAX(1.0, ceil((double)count * nl / LAG_BUDGET)));
            if (m_step > 1) l.lag_skipped = true;
            seg.assign((size_t)members + 1, 0);
            for (size_t k = s.beg; k < s.end; ++k) {
                const int m = l.smember[k];
                if (m < members && m % m_step == 0) seg[m + 1] += 1;
            }
            for (int m = 0; m < members; ++m) seg[m + 1] += seg[m];
            const size_t total = seg[members];
            ord.assign(total, 0);
            {
                std::vector<uint32_t> at(seg.begin(), seg.end() - 1);
                for (size_t k = s.beg; k < s.end; ++k) {
                    const int m = l.smember[k];
                    if (m < members && m % m_step == 0) ord[at[m]++] = (uint32_t)k;
                }
            }
            fr.resize(total);
            xv.resize(total);
            yv.resize(total);
            for (size_t i = 0; i < total; ++i) {
                fr[i] = l.sframe[ord[i]];
                xv[i] = to_corr_value(l.sx[ord[i]], px);
                yv[i] = to_corr_value(l.sy[ord[i]], py);
            }
            ia.resize(total);
            ib.resize(total);
            a.resize(total);
            b.resize(total);

            float* row = l.lag_r.data() + (size_t)p * nl;
            for (int lag = -L; lag <= L; ++lag) {
                size_t n = 0;
                for (int m = 0; m < members; ++m) {
                    const size_t mb = seg[m], len = seg[m + 1] - seg[m];
                    if (len == 0) continue;
                    const size_t c = corr_lag_match(ia.data(), ib.data(), fr.data() + mb, len, lag);
                    for (size_t k = 0; k < c; ++k) {
                        a[n] = xv[mb + ia[k]];
                        b[n] = yv[mb + ib[k]];
                        ++n;
                    }
                }
                double r = 0;
                if (corr_of(&r, kind, a.data(), b.data(), n)) row[lag + L] = (float)r;
            }
        }
    }

    void update() {
        if (!app || !show_window) return;
        if (tree_dirty) apply_tree();

        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            gather(l);
            update_lag(l);
            if (!l.visible || !layer_draws_density(l) || !l.data_valid || md_array_size(l.sx) < 3) continue;

            uint64_t h = md_hash64_combine(l.data_hash, system_generation);
            h = md_hash64(&l.bandwidth, sizeof(l.bandwidth), h);
            h |= 1;
            if (h == l.shown_hash || h == l.want_hash) continue;
            // One at a time: a job under way is let finish (its result is shown meanwhile), the latest is asked for after
            if (l.job) continue;
            launch_job(l, h);
        }

        // Fitted again when the series shown change
        uint64_t fk = 0;
        for (int i = 0; i < num_layers; ++i) {
            const Layer& l = layers[i];
            if (!l.visible) continue;
            fk = md_hash64(l.xs, sizeof(SeriesKey) * l.num_x, fk ^ (uint64_t)l.num_x);
            fk = md_hash64(l.ys, sizeof(SeriesKey) * l.num_y, fk ^ ((uint64_t)l.num_y << 8));
            fk = md_hash64_combine(fk, l.data_valid);
        }
        if (fk != fit_key) {
            fit_key = fk;
            fit_view = true;
        }
    }

    // ## Workspace

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("Correlation"));
        viamd::write_bool(state, STR_LIT("ShowLayers"), show_sidebar);
        viamd::write_bool(state, STR_LIT("ShowLegend"), show_legend);
        viamd::write_bool(state, STR_LIT("ShowLag"), show_lag);
        viamd::write_int(state, STR_LIT("MaxLag"), max_lag);
        viamd::write_flt(state, STR_LIT("LagHeight"), lag_frac);
        // Said even when there are none, so that none stays none (rather than the defaults) once read again
        viamd::write_int(state, STR_LIT("Layers"), num_layers);
        for (int i = 0; i < num_layers; ++i) serialize_layer(state, layers[i]);
    }

    void serialize_layer(viamd::serialization_state_t& state, const Layer& l) {
        viamd::write_section_header(state, STR_LIT("CorrelationLayer"));
        if (l.name[0]) viamd::write_str(state, STR_LIT("Name"), str_from_cstr(l.name));
        viamd::write_bool(state, STR_LIT("Visible"), l.visible);
        char buf[SERIES_PATH_CAP + 64];
        for (int i = 0; i < l.num_x; ++i) {
            key_to_str(buf, sizeof(buf), l.xs[i]);
            viamd::write_str(state, STR_LIT("X"), str_from_cstr(buf));
        }
        for (int i = 0; i < l.num_y; ++i) {
            key_to_str(buf, sizeof(buf), l.ys[i]);
            viamd::write_str(state, STR_LIT("Y"), str_from_cstr(buf));
        }
        static const str_t frames[LayerFrames_Count] = { STR_LIT("current"), STR_LIT("all"), STR_LIT("period") };
        viamd::write_str(state, STR_LIT("Frames"), frames[l.frames]);
        if (l.frames == LayerFrames_Period) viamd::write_int(state, STR_LIT("Period"), (int)l.period);
        static const str_t styles[LayerStyle_Count] = { STR_LIT("points"), STR_LIT("contours"), STR_LIT("heatmap") };
        viamd::write_str(state, STR_LIT("Style"), styles[l.style]);
        viamd::write_vec4(state, STR_LIT("Color"), l.color);
        viamd::write_flt(state, STR_LIT("Opacity"), l.opacity);
        viamd::write_bool(state, STR_LIT("Fill"), l.fill);
        viamd::write_bool(state, STR_LIT("Lines"), l.lines);
        viamd::write_flt(state, STR_LIT("LineWidth"), l.line_width);
        viamd::write_flt_vec(state, STR_LIT("Levels"), l.levels, (size_t)l.num_levels);
        viamd::write_flt(state, STR_LIT("Bandwidth"), l.bandwidth);
        viamd::write_str(state, STR_LIT("Colormap"), str_from_cstr(ImPlot::GetColormapName(l.colormap)));
        viamd::write_bool(state, STR_LIT("LogScale"), l.log_scale);
        static const str_t colorings[PointColoring_Count] = { STR_LIT("layer"), STR_LIT("pair"), STR_LIT("time") };
        viamd::write_str(state, STR_LIT("PointColoring"), colorings[l.point_coloring]);
        viamd::write_flt(state, STR_LIT("PointSize"), l.point_size);
        viamd::write_bool(state, STR_LIT("ShowFit"), l.show_fit);
    }

    void deserialize_settings(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "ShowLayers")) {
                viamd::extract_bool(show_sidebar, arg);
            } else if (str_eq_cstr(ident, "ShowLegend")) {
                viamd::extract_bool(show_legend, arg);
            } else if (str_eq_cstr(ident, "ShowLag")) {
                viamd::extract_bool(show_lag, arg);
            } else if (str_eq_cstr(ident, "MaxLag")) {
                viamd::extract_int(max_lag, arg);
                max_lag = CLAMP(max_lag, 1, MAX_LAG);
            } else if (str_eq_cstr(ident, "LagHeight")) {
                viamd::extract_flt(lag_frac, arg);
                lag_frac = CLAMP(lag_frac, 0.1f, 0.8f);
            } else if (str_eq_cstr(ident, "Layers")) {
                layers_declared = true;
            }
        }
    }

    // As the list is read: each last
    void deserialize_layer(viamd::deserialization_state_t& state) {
        if (!read_layers) {
            // The first of the workspace: the ones there are go
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
            } else if (str_eq_cstr(ident, "Visible")) {
                viamd::extract_bool(l->visible, arg);
            } else if ((str_eq_cstr(ident, "X") || str_eq_cstr(ident, "Y")) && viamd::extract_str(s, arg)) {
                SeriesKey k = {};
                if (key_from_str(&k, s)) {
                    const bool x = str_eq_cstr(ident, "X");
                    int* count = x ? &l->num_x : &l->num_y;
                    if (*count < MAX_AXIS_SERIES) (x ? l->xs : l->ys)[(*count)++] = k;
                }
            } else if (str_eq_cstr(ident, "Frames") && viamd::extract_str(s, arg)) {
                l->frames = str_eq_cstr(s, "current") ? LayerFrames_Current : str_eq_cstr(s, "period") ? LayerFrames_Period : LayerFrames_All;
            } else if (str_eq_cstr(ident, "Period")) {
                int id = 0;
                if (viamd::extract_int(id, arg)) l->period = (uint32_t)MAX(id, 0);
            } else if (str_eq_cstr(ident, "Style") && viamd::extract_str(s, arg)) {
                l->style = str_eq_cstr(s, "contours") ? LayerStyle_Contours : str_eq_cstr(s, "heatmap") ? LayerStyle_Heatmap : LayerStyle_Points;
            } else if (str_eq_cstr(ident, "Color")) {
                viamd::extract_vec4(l->color, arg);
            } else if (str_eq_cstr(ident, "Opacity")) {
                viamd::extract_flt(l->opacity, arg);
            } else if (str_eq_cstr(ident, "Fill")) {
                viamd::extract_bool(l->fill, arg);
            } else if (str_eq_cstr(ident, "Lines")) {
                viamd::extract_bool(l->lines, arg);
            } else if (str_eq_cstr(ident, "LineWidth")) {
                viamd::extract_flt(l->line_width, arg);
            } else if (str_eq_cstr(ident, "Levels")) {
                // As many as there are, up to MAX_LEVELS
                float v[MAX_LEVELS];
                int n = 0;
                for (int k = MAX_LEVELS; k > 0; --k) {
                    if (viamd::extract_flt_vec(v, (size_t)k, arg)) {
                        n = k;
                        break;
                    }
                }
                if (n > 0) {
                    l->num_levels = n;
                    MEMCPY(l->levels, v, sizeof(float) * n);
                }
            } else if (str_eq_cstr(ident, "Bandwidth")) {
                viamd::extract_flt(l->bandwidth, arg);
            } else if (str_eq_cstr(ident, "Colormap") && viamd::extract_str(s, arg)) {
                char name[64];
                str_copy_to_char_buf(name, sizeof(name), s);
                const ImPlotColormap cm = ImPlot::GetColormapIndex(name);
                if (cm >= 0) l->colormap = cm;
            } else if (str_eq_cstr(ident, "LogScale")) {
                viamd::extract_bool(l->log_scale, arg);
            } else if (str_eq_cstr(ident, "PointColoring") && viamd::extract_str(s, arg)) {
                l->point_coloring = str_eq_cstr(s, "layer") ? PointColoring_Layer : str_eq_cstr(s, "time") ? PointColoring_Time : PointColoring_Pair;
            } else if (str_eq_cstr(ident, "PointSize")) {
                viamd::extract_flt(l->point_size, arg);
            } else if (str_eq_cstr(ident, "ShowFit")) {
                viamd::extract_bool(l->show_fit, arg);
            }
        }
        sanitize_layer(*l);
        if (!selected_layer) selected_layer = l->id;
    }

    static void sort_levels(Layer& l) {
        for (int i = 1; i < l.num_levels; ++i) {
            for (int j = i; j > 0 && l.levels[j - 1] > l.levels[j]; --j) {
                const float tmp = l.levels[j];
                l.levels[j] = l.levels[j - 1];
                l.levels[j - 1] = tmp;
            }
        }
    }

    static void sanitize_layer(Layer& l) {
        l.frames = (LayerFrames)CLAMP((int)l.frames, 0, LayerFrames_Count - 1);
        l.style = (LayerStyle)CLAMP((int)l.style, 0, LayerStyle_Count - 1);
        l.opacity = CLAMP(l.opacity, 0.0f, 1.0f);
        l.line_width = CLAMP(l.line_width, 0.5f, 8.0f);
        l.bandwidth = CLAMP(l.bandwidth, 0.5f, 20.0f);
        l.point_size = CLAMP(l.point_size, 1.0f, 16.0f);
        l.num_levels = CLAMP(l.num_levels, 1, MAX_LEVELS);
        for (int i = 0; i < l.num_levels; ++i) l.levels[i] = CLAMP(l.levels[i], 0.01f, 0.999f);
        sort_levels(l);
        if (l.colormap < 0 || l.colormap >= ImPlot::GetColormapCount()) l.colormap = ImPlotColormap_Viridis;
    }

    // ## Drawing

    // What a layer is over, as its name when it has none of its own
    void layer_label(char* buf, size_t cap, const Layer& l) const {
        if (l.name[0]) {
            snprintf(buf, cap, "%s", l.name);
            return;
        }
        if (l.frames == LayerFrames_Current) {
            snprintf(buf, cap, "Current frame");
        } else if (l.frames == LayerFrames_All) {
            snprintf(buf, cap, "All frames");
        } else {
            const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, l.period);
            snprintf(buf, cap, "%s", p ? p->name : "No period");
        }
    }

    // The series of a slot as it is called
    void slot_label(char* buf, size_t cap, const Layer& l, int axis, int i) const {
        const SlotInfo& info = axis == 0 ? l.x_info[i] : l.y_info[i];
        if (info.label[0]) snprintf(buf, cap, "%s", info.label);
        else series_label(buf, cap, app, axis == 0 ? l.xs[i] : l.ys[i]);
    }

    // An axis of a layer as it is labelled: its series, or how many, and their unit when they share it
    void axis_label(char* buf, size_t cap, const Layer& l, int axis) const {
        const int n = axis == 0 ? l.num_x : l.num_y;
        const SlotInfo* info = axis == 0 ? l.x_info : l.y_info;
        if (n == 0) {
            snprintf(buf, cap, "%s", axis == 0 ? "X" : "Y");
            return;
        }
        const char* unit = info[0].unit;
        for (int i = 1; i < n; ++i) {
            if (strcmp(info[i].unit, unit) != 0) unit = "";
        }
        char name[128];
        if (n == 1) slot_label(name, sizeof(name), l, axis, 0);
        else snprintf(name, sizeof(name), "%d series", n);
        if (unit[0]) snprintf(buf, cap, "%s (%s)", name, unit);
        else snprintf(buf, cap, "%s", name);
    }

    // A pair as it is called: what varies along it, or both of its series
    void pair_label(char* buf, size_t cap, const Layer& l, int p) const {
        char x[96], y[96];
        slot_label(x, sizeof(x), l, 0, pair_x(l, p));
        slot_label(y, sizeof(y), l, 1, pair_y(l, p));
        snprintf(buf, cap, "%s vs %s", y, x);
    }

    // The colour of a pair (and member): of the series that varies between pairs, or of its members by the colour map
    ImVec4 pair_color(const Layer& l, int p, int m) const {
        if (l.num_pairs > 1) {
            return l.num_x > 1 && l.num_y == 1 ? l.x_info[p].color : l.y_info[pair_y(l, p)].color;
        }
        const int members = l.pairs[0].members;
        if (members > 1) return ImPlot::SampleColormap((float)m / (float)(members - 1), l.colormap);
        return vec_cast(l.color);
    }

    ImVec4 point_color(const Layer& l, size_t idx) const {
        switch (l.point_coloring) {
        case PointColoring_Pair:
            return pair_color(l, l.spair[idx], l.smember[idx]);
        case PointColoring_Time: {
            if (!(l.t_max > l.t_min)) return vec_cast(l.color);
            return ImPlot::SampleColormap((l.st[idx] - l.t_min) / (l.t_max - l.t_min), l.colormap);
        }
        default:
            return vec_cast(l.color);
        }
    }

    // Fills what layer l is drawn with, false when there is nothing of it to draw
    bool layer_draw(DensityMapDraw* d, const Layer& l) const {
        density_map_draw_init(d);
        if (!layer_draws_density(l) || !l.has_result || !l.tex) return false;
        const DensityMapResult& r = l.result;
        if (r.samples[0] == 0 || r.max[0] <= 0.0f) return false;
        d->tex = l.tex;
        d->opacity = l.opacity;
        d->line_width = l.line_width;
        MEMCPY(d->line_color, l.color.elem, sizeof(d->line_color));
        d->mode = l.style == LayerStyle_Heatmap ? 1 : 0;
        d->fill = l.fill;
        d->lines = l.lines;
        if (l.style == LayerStyle_Heatmap) density_map_draw_colormap(d, l.colormap, l.log_scale);
        density_map_draw_levels(d, r, 0, l.levels, l.num_levels);
        density_map_draw_fills_from(d, vec_cast(l.color));
        // Over a domain of its own: each layer's is its samples' extent
        d->own_view_uv = 1;
        density_map_plot_view_uv(d->view_uv, l.result_domain);
        d->content = l.shown_hash;
        return true;
    }

    // The densities of the layers, composited (density_map_composite): 0 when there are none
    uint32_t render_canvas() {
        DensityMapDraw draws[MAX_LAYERS];
        int num_draws = 0;
        // Bottom up: the first layer is on top
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!l.visible || l.opacity <= 0.0f) continue;
            if (layer_draw(&draws[num_draws], l)) ++num_draws;
        }
        int w = 0, h = 0;
        density_map_plot_canvas_size(&w, &h);
        if (num_draws == 0 || w <= 0 || h <= 0) return 0;
        return density_map_composite(&canvas, w, h, draws[0].view_uv, draws, num_draws);
    }

    // ### Points

    // What the pointer is on
    struct Hover {
        int    layer = -1;
        size_t idx = 0;
    };

    void draw_tooltip(const Hover& h) const {
        const Layer& l = layers[h.layer];
        const int p = l.spair[h.idx];
        const int m = l.smember[h.idx];
        const int xi = pair_x(l, p);
        const int yi = pair_y(l, p);
        char x[96], y[96];
        slot_label(x, sizeof(x), l, 0, xi);
        slot_label(y, sizeof(y), l, 1, yi);
        ImGui::BeginTooltip();
        ImGui::Text("%s  %.4g %s", x, l.sx[h.idx], l.x_info[xi].unit);
        ImGui::Text("%s  %.4g %s", y, l.sy[h.idx], l.y_info[yi].unit);
        if (l.pairs[p].members > 1) {
            const SeriesKey& key = l.y_info[yi].dim > 1 ? l.ys[yi] : l.xs[xi];
            const str_t member = series_member_label(app, key, m);
            if (!str_empty(member)) ImGui::TextDisabled("member %d \xc2\xb7 " STR_FMT, m + 1, STR_ARG(member));
            else ImGui::TextDisabled("member %d", m + 1);
        }
        const bool binned = l.sframe[h.idx] == NO_FRAME;
        if (binned) ImGui::TextDisabled("at %.4g %s", l.st[h.idx], l.x_info[xi].bin_unit);
        else ImGui::TextDisabled("frame %u \xc2\xb7 t = %.4g", l.sframe[h.idx], l.st[h.idx]);
        char label[128];
        layer_label(label, sizeof(label), l);
        if (binned) ImGui::TextDisabled("%s", label);
        else ImGui::TextDisabled("%s \xc2\xb7 click to go to the frame", label);
        ImGui::EndTooltip();
    }

    // A pair's correlation as text: Pearson's and Spearman's, or of periodic values the circular one
    static void stats_text(char* buf, size_t cap, const PairStats& s) {
        if (s.circular) {
            if (s.has_corr()) snprintf(buf, cap, "circular r = %.3f", s.corr());
            else snprintf(buf, cap, "circular r = -");
            return;
        }
        char r[24] = "r = -", rho[24] = "";
        if (s.has_r()) snprintf(r, sizeof(r), "r = %.3f", s.r());
        if (s.rho_valid) snprintf(rho, sizeof(rho), " \xc2\xb7 \xcf\x81 = %.3f", s.rho);
        snprintf(buf, cap, "%s%s", r, rho);
    }

    // What the series of the hovered point show in the viewport: the script's visualization of Y, else of X
    void visualize(const Hover& h) {
        const Layer& l = layers[h.layer];
        const int p = l.spair[h.idx];
        const int m = l.smember[h.idx];
        const SlotInfo* infos[2] = { &l.y_info[pair_y(l, p)], &l.x_info[pair_x(l, p)] };
        for (const SlotInfo* info : infos) {
            if (!info->ok || !info->script_ident[0]) continue;
            const int k = info->dim > 1 ? m : -1;
            script_visualize_ref(app, info->vis_ref, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(app, str_from_cstr(info->script_ident), k);
            highlights_now = true;
            return;
        }
    }

    // The points of the visible layers: hovered, then drawn; with the lines fitted through them
    void draw_points(bool hovered) {
        const ImVec2 mouse = ImGui::GetMousePos();
        Hover hover;
        if (hovered) {
            float best = FLT_MAX;
            // The layer on top has it, where it has a point at all
            for (int i = 0; i < num_layers && hover.layer < 0; ++i) {
                const Layer& l = layers[i];
                if (!l.visible || !l.data_valid || !layer_draws_points(l)) continue;
                const float r = MAX(6.0f, l.point_size);
                const float r2 = r * r;
                const size_t n = md_array_size(l.sx);
                for (size_t k = 0; k < n; k += l.point_stride) {
                    const ImVec2 p = ImPlot::PlotToPixels(l.sx[k], l.sy[k]);
                    const float dx = p.x - mouse.x;
                    const float dy = p.y - mouse.y;
                    const float d2 = dx * dx + dy * dy;
                    if (d2 < best && d2 < r2) {
                        best = d2;
                        hover = Hover{ i, k };
                    }
                }
            }
        }

        const uint32_t frame_cur = (uint32_t)MAX(0.0, app->animation.frame + 0.5);
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        density_map_marks_begin(dl);
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!l.visible || !l.data_valid || !layer_draws_points(l)) continue;
            const float r = l.point_size * 0.5f;
            const bool many = l.frames != LayerFrames_Current;   // Over frames: small squares, no outlines
            const size_t n = md_array_size(l.sx);
            for (size_t k = 0; k < n; k += l.point_stride) {
                ImVec4 c = point_color(l, k);
                c.w *= l.opacity;
                const ImVec2 p = ImPlot::PlotToPixels(l.sx[k], l.sy[k]);
                if (many) density_map_mark(dl, p, r, DensityMapMark_Square, ImGui::ColorConvertFloat4ToU32(c));
                else density_map_mark(dl, p, r + 1.0f, DensityMapMark_Disc, ImGui::ColorConvertFloat4ToU32(c));
            }
        }
        // Where the frame shown is among the points over frames, and what is hovered
        const ImU32 ring = ImGui::ColorConvertFloat4ToU32(vec_cast(app->selection.color.highlight.outline));
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!l.visible || !l.data_valid || !layer_draws_points(l) || l.frames == LayerFrames_Current) continue;
            const float r = MAX(l.point_size * 0.5f + 2.5f, 4.0f);
            const size_t n = md_array_size(l.sx);
            for (size_t k = 0; k < n; ++k) {
                if (l.sframe[k] != frame_cur) continue;
                density_map_mark(dl, ImPlot::PlotToPixels(l.sx[k], l.sy[k]), r + 1.0f, DensityMapMark_Ring, IM_COL32(255, 255, 255, 220));
            }
        }
        if (hover.layer >= 0) {
            const Layer& l = layers[hover.layer];
            const float r = MAX(l.point_size * 0.5f + 2.5f, 4.0f);
            density_map_mark(dl, ImPlot::PlotToPixels(l.sx[hover.idx], l.sy[hover.idx]), r + 1.0f, DensityMapMark_Ring, ring);
        }
        density_map_marks_end(dl);

        // The least squares line of each pair, over the range of its X
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!l.visible || !l.data_valid || !l.show_fit) continue;
            for (int p = 0; p < l.num_pairs; ++p) {
                const PairStats& s = l.pairs[p];
                if (!s.has_fit()) continue;
                const double a = s.intercept(), b = s.slope();
                const ImVec2 p0 = ImPlot::PlotToPixels(s.x_min, a + b * s.x_min);
                const ImVec2 p1 = ImPlot::PlotToPixels(s.x_max, a + b * s.x_max);
                ImVec4 c = l.point_coloring == PointColoring_Pair && l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color);
                c.w *= MAX(l.opacity, 0.5f);
                dl->AddLine(p0, p1, IM_COL32(0, 0, 0, 160), 3.5f);
                dl->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32(c), 1.75f);
            }
        }
        ImPlot::PopPlotClipRect();

        if (hover.layer >= 0) {
            draw_tooltip(hover);
            visualize(hover);
            // A click (not a pan) on a point: its frame
            if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && ImGui::GetIO().MouseDragMaxDistanceSqr[ImGuiMouseButton_Left] < 9.0f &&
                layers[hover.layer].sframe[hover.idx] != NO_FRAME) {
                app->animation.frame = (double)layers[hover.layer].sframe[hover.idx];
            }
        }
    }

    // The layer whose series name the axes: the selected one if it is shown, else the first shown with series
    const Layer* axes_layer() const {
        const Layer* sel = find_layer(selected_layer);
        if (sel && sel->visible && (sel->num_x > 0 || sel->num_y > 0)) return sel;
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].visible && (layers[i].num_x > 0 || layers[i].num_y > 0)) return &layers[i];
        }
        return sel;
    }

    // The view around the samples of the layers shown
    void fit_to_layers() {
        bool any = false;
        double x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        for (int i = 0; i < num_layers; ++i) {
            const Layer& l = layers[i];
            if (!l.visible || !l.data_valid) continue;
            if (!any) {
                x0 = l.x_min; x1 = l.x_max; y0 = l.y_min; y1 = l.y_max;
                any = true;
            } else {
                x0 = MIN(x0, (double)l.x_min); x1 = MAX(x1, (double)l.x_max);
                y0 = MIN(y0, (double)l.y_min); y1 = MAX(y1, (double)l.y_max);
            }
        }
        if (!any) return;
        auto pad = [](double* lo, double* hi) {
            double r = *hi - *lo;
            if (!(r > 0.0)) r = MAX(fabs(*lo), 1.0) * 0.1;
            *lo -= r * 0.05;
            *hi += r * 0.05;
        };
        pad(&x0, &x1);
        pad(&y0, &y1);
        viewrect = ImPlotRect(x0, x1, y0, y1);
        fit_view = false;
    }

    // While a series is dragged: the strips along the axes it is dropped on (and the axes themselves)
    void draw_series_drop_targets(ImVec2 pos, ImVec2 size) {
        const ImGuiPayload* payload = ImGui::GetDragDropPayload();
        const SeriesDragPayload* sp = series_payload(payload);
        if (!sp) return;
        const char* reason = nullptr;
        const bool ok = series_droppable(sp->key, &reason);

        const float band_h = CLAMP(size.y * 0.25f, 36.0f, 120.0f);
        const float band_w = CLAMP(size.x * 0.25f, 36.0f, 140.0f);
        const ImRect rx(ImVec2(pos.x + band_w, pos.y + size.y - band_h), ImVec2(pos.x + size.x, pos.y + size.y));
        const ImRect ry(ImVec2(pos.x, pos.y), ImVec2(pos.x + band_w, pos.y + size.y - band_h));
        const ImRect rects[2] = { rx, ry };
        const char* names[2] = { "X", "Y" };
        const ImAxis axes[2] = { ImAxis_X1, ImAxis_Y1 };

        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        for (int a = 0; a < 2; ++a) {
            const bool inside = rects[a].Contains(ImGui::GetMousePos());
            const ImVec4 c = ok ? ImVec4(0.35f, 0.75f, 1.0f, inside ? 0.30f : 0.12f) : ImVec4(1.0f, 0.4f, 0.4f, 0.12f);
            dl->AddRectFilled(rects[a].Min, rects[a].Max, ImGui::ColorConvertFloat4ToU32(c), 4.0f);
            dl->AddRect(rects[a].Min, rects[a].Max, ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 0.8f)), 4.0f, 0, 1.0f);
            char text[64];
            snprintf(text, sizeof(text), ok ? "Drop on %s" : "Cannot plot on %s", names[a]);
            const ImVec2 ts = ImGui::CalcTextSize(text);
            dl->AddText(ImVec2(rects[a].GetCenter().x - ts.x * 0.5f, rects[a].GetCenter().y - ts.y * 0.5f), IM_COL32(255, 255, 255, 230), text);
        }
        ImPlot::PopPlotClipRect();
        if (!ok && reason) ImGui::SetTooltip("%s", reason);

        for (int a = 0; a < 2; ++a) {
            ImGui::PushID(a);
            bool open = ImGui::BeginDragDropTargetCustom(rects[a], ImGui::GetID("##series_band"));
            if (open) {
                if (const SeriesDragPayload* dp = accept_series()) {
                    if (ok) drop_series(a, dp->key);
                }
                ImGui::EndDragDropTarget();
            }
            if (ImPlot::BeginDragDropTargetAxis(axes[a])) {
                if (const SeriesDragPayload* dp = accept_series()) {
                    if (ok) drop_series(a, dp->key);
                }
                ImPlot::EndDragDropTarget();
            }
            ImGui::PopID();
        }
    }

    // The pairs of the layer the axes are of, with their correlations, in the corner
    void draw_legend(const Layer& l, ImVec2 pos) {
        if (!l.data_valid || l.num_pairs == 0) return;
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        const float lh = ImGui::GetTextLineHeight();
        const float pad = 6.0f;
        char lines[MAX_PAIRS][192];
        float w = 0.0f;
        for (int p = 0; p < l.num_pairs; ++p) {
            char name[160];
            pair_label(name, sizeof(name), l, p);
            const PairStats& s = l.pairs[p];
            char st[64];
            stats_text(st, sizeof(st), s);
            snprintf(lines[p], sizeof(lines[p]), "%s   %s", name, st);
            w = MAX(w, ImGui::CalcTextSize(lines[p]).x);
        }
        const ImVec2 p0 = ImVec2(pos.x + 8, pos.y + 8);
        const ImVec2 p1 = ImVec2(p0.x + w + lh + pad * 3, p0.y + l.num_pairs * (lh + 2) + pad * 2 - 2);
        ImPlot::PushPlotClipRect();
        dl->AddRectFilled(p0, p1, IM_COL32(0, 0, 0, 140), 4.0f);
        for (int p = 0; p < l.num_pairs; ++p) {
            const float y = p0.y + pad + p * (lh + 2);
            const ImVec4 c = l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color);
            dl->AddRectFilled(ImVec2(p0.x + pad, y + lh * 0.2f), ImVec2(p0.x + pad + lh * 0.6f, y + lh * 0.8f), ImGui::ColorConvertFloat4ToU32(c), 2.0f);
            dl->AddText(ImVec2(p0.x + pad * 2 + lh * 0.6f, y), IM_COL32(255, 255, 255, 230), lines[p]);
        }
        ImPlot::PopPlotClipRect();
    }

    void draw_plot() {
        if (fit_view) fit_to_layers();

        constexpr ImPlotFlags flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
        // Opaque, whatever the windows are: what is behind the window is not to show through the data
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));
        if (!ImPlot::BeginPlot("##correlation", ImVec2(-1, -1), flags)) {
            ImPlot::PopStyleColor();
            return;
        }
        const Layer* al = axes_layer();
        char xl[160] = "X", yl[160] = "Y";
        if (al) {
            axis_label(xl, sizeof(xl), *al, 0);
            axis_label(yl, sizeof(yl), *al, 1);
        }
        ImPlot::SetupAxes(xl, yl, 0, 0);
        ImPlot::SetupAxisLinks(ImAxis_X1, &viewrect.X.Min, &viewrect.X.Max);
        ImPlot::SetupAxisLinks(ImAxis_Y1, &viewrect.Y.Min, &viewrect.Y.Max);
        ImPlot::SetupFinish();

        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 psize = ImPlot::GetPlotSize();

        // The densities, composited, under everything else
        density_map_plot_canvas(render_canvas());

        draw_points(ImPlot::IsPlotHovered());

        if (show_legend && al) draw_legend(*al, pos);

        // Why there is nothing to see
        bool any = false;
        for (int i = 0; i < num_layers; ++i) any |= layers[i].visible && layers[i].data_valid;
        if (!any) {
            const char* text = num_layers == 0 ? "Add a layer with +" : al && al->data_error[0] ? al->data_error : "drop series on the X and Y axes";
            char buf[128];
            snprintf(buf, sizeof(buf), "%c%s", toupper(text[0]), text + 1);
            const ImVec2 ts = ImGui::CalcTextSize(buf);
            ImPlot::GetPlotDrawList()->AddText(ImVec2(pos.x + (psize.x - ts.x) * 0.5f, pos.y + (psize.y - ts.y) * 0.5f), ImGui::GetColorU32(ImGuiCol_TextDisabled), buf);
        }

        draw_series_drop_targets(pos, psize);

        // A period from the timeline: a new layer of the selected one's series over it
        if (ImPlot::BeginDragDropTargetPlot()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                ASSERT(payload->DataSize == sizeof(uint32_t));
                add_period_layer(*(const uint32_t*)payload->Data);
            }
            ImPlot::EndDragDropTarget();
        }

        // A double-click fits the view, as the plots have it
        if (ImPlot::FitThisFrame()) fit_view = true;

        ImPlot::EndPlot();
        ImPlot::PopStyleColor();
    }

    // Of the selected layer's series, first in the list
    Layer* add_period_layer(uint32_t period_id) {
        const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, period_id);
        if (!p) return nullptr;
        const Layer* src = find_layer(selected_layer);
        const uint32_t src_id = src ? src->id : 0;
        Layer* l = add_layer(LayerFrames_Period, LayerStyle_Contours, 0);
        if (!l) return nullptr;
        if (const Layer* s = find_layer(src_id)) copy_series(*l, *s);
        l->period = period_id;
        l->color = p->color;
        l->color.w = 1.0f;
        selected_layer = l->id;
        return l;
    }

    // The correlation of each pair against the lag of Y behind X, under the plot: of the layers shown over frames
    void draw_lag_plot() {
        constexpr ImPlotFlags flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));
        if (!ImPlot::BeginPlot("##lag", ImVec2(-1, -1), flags)) {
            ImPlot::PopStyleColor();
            return;
        }
        const int L = CLAMP(max_lag, 1, MAX_LAG);
        bool any_circular = false, any_linear = false;
        for (int i = 0; i < num_layers; ++i) {
            const Layer& l = layers[i];
            if (!l.visible || !layer_has_lag(l)) continue;
            for (int p = 0; p < l.num_pairs; ++p) (l.pairs[p].circular ? any_circular : any_linear) = true;
        }
        ImPlot::SetupAxis(ImAxis_X1, "lag of Y behind X (frames)");
        ImPlot::SetupAxis(ImAxis_Y1, any_circular && !any_linear ? "circular r" : "r", ImPlotAxisFlags_Lock);
        ImPlot::SetupAxisLimits(ImAxis_X1, -L - 0.5, L + 0.5, L != lag_shown ? ImPlotCond_Always : ImPlotCond_Once);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.05, 1.05, ImPlotCond_Always);
        ImPlot::SetupFinish();
        lag_shown = L;

        static std::vector<float> xs;
        xs.resize(2 * L + 1);
        for (int k = 0; k <= 2 * L; ++k) xs[k] = (float)(k - L);

        // Zero lag, and no correlation
        const double zero = 0.0;
        ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.25f), 1.0f);
        ImPlot::PlotInfLines("##zero_lag", &zero, 1);
        ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.25f), 1.0f);
        ImPlot::PlotInfLines("##zero_r", &zero, 1, ImPlotInfLinesFlags_Horizontal);

        const bool hovered = ImPlot::IsPlotHovered();
        const int tau = CLAMP((int)lround(ImPlot::GetPlotMousePos().x), -L, L);
        int curves = 0;
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!l.visible || !layer_has_lag(l) || l.lag_max != L || l.lag_r.size() != (size_t)l.num_pairs * (2 * L + 1)) continue;
            for (int p = 0; p < l.num_pairs; ++p) {
                const float* row = l.lag_r.data() + (size_t)p * (2 * L + 1);
                ImVec4 c = l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color);
                c.w = MAX(l.opacity, 0.5f);
                char id[32];
                snprintf(id, sizeof(id), "##lag%u_%d", l.id, p);
                ImPlot::SetNextLineStyle(c, 1.75f);
                ImPlot::PlotLine(id, xs.data(), row, 2 * L + 1);
                // Where it is strongest
                int best = -1;
                for (int k = 0; k <= 2 * L; ++k) {
                    if (isfinite(row[k]) && (best < 0 || fabsf(row[k]) > fabsf(row[best]))) best = k;
                }
                if (best >= 0) {
                    snprintf(id, sizeof(id), "##peak%u_%d", l.id, p);
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3.5f, c, 1.0f, ImVec4(0, 0, 0, 0.8f));
                    ImPlot::PlotScatter(id, &xs[best], &row[best], 1);
                }
                ++curves;
            }
        }

        if (curves == 0) {
            const char* text = "Of layers over the frames of a trajectory, not of distributions";
            const ImVec2 ts = ImGui::CalcTextSize(text);
            const ImVec2 pos = ImPlot::GetPlotPos(), size = ImPlot::GetPlotSize();
            ImPlot::GetPlotDrawList()->AddText(ImVec2(pos.x + (size.x - ts.x) * 0.5f, pos.y + (size.y - ts.y) * 0.5f), ImGui::GetColorU32(ImGuiCol_TextDisabled), text);
        } else if (hovered) {
            const double t = (double)tau;
            ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.6f), 1.0f);
            ImPlot::PlotInfLines("##cursor", &t, 1);
            draw_lag_tooltip(tau, L);
        }

        ImPlot::EndPlot();
        ImPlot::PopStyleColor();
    }

    void draw_lag_tooltip(int tau, int L) const {
        ImGui::BeginTooltip();
        const float* xv = app->timeline.x_values;
        const int64_t nf = (int64_t)md_array_size(xv);
        char unit[32] = "";
        md_unit_print(unit, sizeof(unit), app->timeline.time_unit);
        if (nf >= 2) {
            const double dt = (xv[nf - 1] - xv[0]) / (double)(nf - 1);
            ImGui::Text("Y %d frames %s X  (%.4g %s)", abs(tau), tau >= 0 ? "after" : "before", fabs(tau * dt), unit);
        } else {
            ImGui::Text("Y %d frames %s X", abs(tau), tau >= 0 ? "after" : "before");
        }
        bool skipped = false;
        for (int i = 0; i < num_layers; ++i) {
            const Layer& l = layers[i];
            if (!l.visible || !layer_has_lag(l) || l.lag_max != L || l.lag_r.size() != (size_t)l.num_pairs * (2 * L + 1)) continue;
            skipped |= l.lag_skipped;
            for (int p = 0; p < l.num_pairs; ++p) {
                const float* row = l.lag_r.data() + (size_t)p * (2 * L + 1);
                int best = -1;
                for (int k = 0; k <= 2 * L; ++k) {
                    if (isfinite(row[k]) && (best < 0 || fabsf(row[k]) > fabsf(row[best]))) best = k;
                }
                const ImVec4 c = l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color);
                ImGui::ColorButton("##c", c, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(ImGui::GetTextLineHeight() * 0.7f, ImGui::GetTextLineHeight() * 0.7f));
                ImGui::SameLine();
                char name[192];
                pair_label(name, sizeof(name), l, p);
                const float r = row[tau + L];
                if (isfinite(r)) ImGui::Text("%s   %s = %.3f", name, l.pairs[p].circular ? "circular r" : "r", r);
                else ImGui::Text("%s   -", name);
                if (best >= 0) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("(strongest at %+d: %.3f)", best - L, row[best]);
                }
            }
        }
        if (skipped) ImGui::TextDisabled("Some members of populations left out: too many samples");
        ImGui::EndTooltip();
    }

    void draw_toolbar() {
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6, 2));
        if (ImGui::Button(show_sidebar ? ICON_FA_LAYER_GROUP "##layers" : ICON_FA_LAYER_GROUP "##layers_off")) {
            show_sidebar = !show_sidebar;
        }
        ImGui::SetItemTooltip(show_sidebar ? "Hide the layers" : "Show the layers");
        ImGui::SameLine();
        if (show_lag) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
        if (ImGui::Button(ICON_FA_CHART_LINE "##lag")) show_lag = !show_lag;
        if (show_lag) ImGui::PopStyleColor();
        ImGui::SetItemTooltip(show_lag ? "Hide the correlation against the lag" : "Show the correlation against the lag: of X at a frame with Y frames later (or earlier)");

        // On the right: fitting the view and the settings
        ImGui::SameLine();
        const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + ImGui::CalcTextSize(ICON_FA_GEAR).x + ImGui::GetStyle().FramePadding.x * 4 + ImGui::GetStyle().ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
        if (ImGui::Button(ICON_FA_EXPAND)) {
            fit_view = true;
        }
        ImGui::SetItemTooltip("Fit the view (or double-click the plot)");
        ImGui::SameLine();
        if (ImGui::Button(ICON_FA_GEAR)) ImGui::OpenPopup("##settings");
        if (ImGui::BeginPopup("##settings")) {
            ImGui::Checkbox("Legend", &show_legend);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
            ImGui::SliderInt("Lag up to", &max_lag, 1, MAX_LAG, "%d frames", ImGuiSliderFlags_Logarithmic | ImGuiSliderFlags_AlwaysClamp);
            ImGui::SetItemTooltip("How far either way the correlation against the lag goes");
            ImGui::Separator();
            ImGui::TextDisabled("Drag series from the timeline onto the X or Y axis of the plot,");
            ImGui::TextDisabled("or onto the lists of the selected layer.");
            ImGui::TextDisabled("Drop a period from the timeline on the plot to add a layer over it.");
            ImGui::TextDisabled("Click a point to go to its frame.");
            ImGui::EndPopup();
        }
        ImGui::PopStyleVar();
    }

    // Statistics of a layer as one line: its correlation, or how many pairs
    void layer_stats(char* buf, size_t cap, const Layer& l) const {
        if (!l.data_valid) {
            snprintf(buf, cap, "%s", l.data_error[0] ? l.data_error : "no samples");
            return;
        }
        const size_t n = md_array_size(l.sx);
        const char* what = l.binned ? "bins" : "samples";
        if (l.num_pairs == 1) {
            char st[64];
            stats_text(st, sizeof(st), l.pairs[0]);
            snprintf(buf, cap, "%s \xc2\xb7 %zu %s", st, n, what);
        } else {
            snprintf(buf, cap, "%d pairs \xc2\xb7 %zu %s", l.num_pairs, n, what);
        }
    }

    // The pairs of a layer, with their statistics
    void draw_pair_table(const Layer& l) const {
        if (!ImGui::BeginTable("##pairs", 5, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg)) return;
        ImGui::TableSetupColumn("pair", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("n");
        ImGui::TableSetupColumn("r");
        ImGui::TableSetupColumn("\xcf\x81");
        ImGui::TableSetupColumn("slope");
        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        static const char* tips[5] = {
            nullptr,
            "Samples: frames (times members), or places along the bins of distributions",
            "Pearson's correlation. Of periodic values (angles), the circular one: Jammalamadaka-SenGupta's of two (-1..1), Mardia's of one against a linear value (0..1, no sign)",
            "Spearman's rank correlation: how monotonic, whatever the shape",
            "Of the least squares line, Y on X",
        };
        for (int c = 0; c < 5; ++c) {
            ImGui::TableSetColumnIndex(c);
            ImGui::TableHeader(ImGui::TableGetColumnName(c));
            if (tips[c]) ImGui::SetItemTooltip("%s", tips[c]);
        }
        for (int p = 0; p < l.num_pairs; ++p) {
            const PairStats& s = l.pairs[p];
            char name[192];
            pair_label(name, sizeof(name), l, p);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            const ImVec4 c = l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color);
            ImGui::ColorButton("##c", c, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(ImGui::GetTextLineHeight() * 0.7f, ImGui::GetTextLineHeight() * 0.7f));
            ImGui::SameLine();
            ImGui::TextUnformatted(name);
            if (s.members > 1) {
                ImGui::SameLine();
                ImGui::TextDisabled("(%d)", s.members);
            }
            ImGui::TableNextColumn();
            ImGui::Text("%llu", (unsigned long long)s.n);
            ImGui::TableNextColumn();
            if (s.has_corr()) ImGui::Text("%.3f", s.corr());
            else ImGui::TextDisabled("-");
            if (s.circular) ImGui::SetItemTooltip(s.circular == 2 ? "Circular (both periodic)" : "Circular-linear: how much the one predicts the other, no sign");
            ImGui::TableNextColumn();
            if (s.rho_valid) ImGui::Text("%.3f", s.rho);
            else ImGui::TextDisabled("-");
            if (s.rho_valid && s.rho_n < s.n) ImGui::SetItemTooltip("Of %llu samples, evenly spaced", (unsigned long long)s.rho_n);
            ImGui::TableNextColumn();
            if (s.has_fit()) ImGui::Text("%.4g", s.slope());
            else ImGui::TextDisabled("-");
        }
        ImGui::EndTable();
    }

    void draw_layer_tooltip(const Layer& l) const {
        char label[128];
        layer_label(label, sizeof(label), l);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(label);
        if (!l.data_valid) {
            ImGui::TextColored(error_color, "%s", l.data_error[0] ? l.data_error : "no samples");
            ImGui::EndTooltip();
            return;
        }
        draw_pair_table(l);
        if (layer_draws_points(l) && l.point_stride > 1) {
            ImGui::TextDisabled("Every %u samples are shown", l.point_stride);
        }
        ImGui::EndTooltip();
    }

    // What the rows asked for, done once the list is drawn
    struct SidebarActions {
        uint32_t remove = 0;
        uint32_t duplicate = 0;
        bool     move = false;
        LayerNodeRef move_src, move_target;
        RepresentationDrop move_drop = RepresentationDrop_Before;
    };

    void draw_layer_row(Layer& l, RepresentationNode* const nodes[], size_t num_nodes, SidebarActions* act) {
        ImGui::PushID("layer");
        ImGui::PushID((int)l.id);

        const float line_h = ImGui::GetTextLineHeight();
        const float row_h = line_h * 2.0f + 8.0f;
        const ImVec2 row_min = ImGui::GetCursorScreenPos();
        const float row_w = ImGui::GetContentRegionAvail().x;
        const bool selected = selected_layer == l.id;
        const ImRect row_rect(row_min, ImVec2(row_min.x + row_w, row_min.y + row_h));

        char label[128];
        layer_label(label, sizeof(label), l);

        // The row, a selectable under the eye, the swatch and the text
        ImGui::SetNextItemAllowOverlap();
        if (ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h))) {
            selected_layer = l.id;
        }
        const bool row_hovered = ImGui::IsItemHovered();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal | ImGuiHoveredFlags_NoSharedDelay) && !ImGui::IsMouseDragging(0)) {
            draw_layer_tooltip(l);
        }
        layer_tree_drag_source({LIST_ID, l.id, false}, label);
        if (ImGui::BeginPopupContextItem("##menu")) {
            selected_layer = l.id;
            if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false, num_layers < MAX_LAYERS)) act->duplicate = l.id;
            if (ImGui::MenuItem(ICON_FA_RIGHT_LEFT " Swap X and Y", nullptr, false, l.num_x + l.num_y > 0)) swap_axes(l);
            ImGui::Separator();
            if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) act->remove = l.id;
            ImGui::EndPopup();
        }
        // A period from the timeline: this layer over its frames; a series: onto the axis it fits
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                l.frames = LayerFrames_Period;
                l.period = *(const uint32_t*)payload->Data;
                selected_layer = l.id;
            }
            ImGui::EndDragDropTarget();
        }
        // A layer of the list: moved here
        LayerNodeRef src;
        RepresentationDrop drop;
        if (layer_tree_drop_target(row_rect, row_min.x, {LIST_ID, l.id, false}, false, nodes, num_nodes, &src, &drop)) {
            act->move = true;
            act->move_src = src;
            act->move_target = {LIST_ID, l.id, false};
            act->move_drop = drop;
        }
        if (row_hovered && l.frames == LayerFrames_Period) {
            timeline_periods_highlight(l.period);
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float pad = 4.0f;
        const ImVec2 text_col_min = ImVec2(row_min.x + pad + line_h + 6 + line_h + 6, row_min.y + 4);

        // The eye
        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, row_min.y + (row_h - ImGui::GetFrameHeight()) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        if (ImGui::Button(l.visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) {
            l.visible = !l.visible;
        }
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip(l.visible ? "Hide" : "Show");

        // The swatch: what it looks like
        const ImVec2 sw_min = ImVec2(row_min.x + pad + line_h + 8, row_min.y + (row_h - line_h) * 0.5f);
        const ImVec2 sw_max = ImVec2(sw_min.x + line_h, sw_min.y + line_h);
        const float vis_alpha = l.visible ? 1.0f : 0.35f;
        const DensityMapSwatch kind = l.style == LayerStyle_Points ? DensityMapSwatch_Points : l.style == LayerStyle_Contours ? DensityMapSwatch_Contours : DensityMapSwatch_Heatmap;
        density_map_swatch(dl, sw_min, sw_max, kind, vec_cast(l.color), l.fill, l.colormap, vis_alpha);

        // The name, and below it what it holds
        char stats[160];
        layer_stats(stats, sizeof(stats), l);
        const ImU32 text_col = ImGui::GetColorU32(l.visible ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        // On the selection's colour the disabled one is too faint to read
        const ImU32 dim_col = selected ? ImGui::GetColorU32(ImGuiCol_Text, 0.7f) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        const ImVec4 clip = ImVec4(text_col_min.x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), text_col_min, text_col, label, nullptr, 0.0f, &clip);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(text_col_min.x, text_col_min.y + line_h + 1),
                    !l.data_valid ? ImGui::ColorConvertFloat4ToU32(error_color) : dim_col, stats, nullptr, 0.0f, &clip);

        // Working on its density
        if (l.job) {
            const float t = (float)ImGui::GetTime();
            const ImVec2 c = ImVec2(row_min.x + row_w - pad - 6, row_min.y + 9);
            dl->PathArcTo(c, 4.5f, t * 6.0f, t * 6.0f + 4.0f, 10);
            dl->PathStroke(dim_col, 0, 1.5f);
        }

        ImGui::SetCursorScreenPos(ImVec2(row_min.x, row_min.y + row_h + 2));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
        ImGui::PopID();
    }

    // A new layer, first in the list, of the selected layer's series
    Layer* add_layer_like_selected(LayerFrames frames, LayerStyle style) {
        const Layer* src = find_layer(selected_layer);
        const uint32_t src_id = src ? src->id : 0;
        Layer* l = add_layer(frames, style, 0);
        if (!l) return nullptr;
        if (const Layer* s = find_layer(src_id)) copy_series(*l, *s);
        selected_layer = l->id;
        return l;
    }

    void draw_add_menu() {
        const bool room = num_layers < MAX_LAYERS;
        if (ImGui::MenuItem("Current frame", nullptr, false, room)) {
            if (Layer* l = add_layer_like_selected(LayerFrames_Current, LayerStyle_Points)) l->point_size = 8.0f;
        }
        if (ImGui::MenuItem("All frames", nullptr, false, room)) {
            add_layer_like_selected(LayerFrames_All, LayerStyle_Points);
        }
        const TimelinePeriods& ps = app->timeline.periods;
        if (ImGui::BeginMenu("Period", ps.count > 0 && room)) {
            for (int i = 0; i < ps.count; ++i) {
                char label[TIMELINE_PERIOD_NAME_CAP + 16];
                snprintf(label, sizeof(label), "%s###p%u", ps.items[i].name, ps.items[i].id);
                if (ImGui::MenuItem(label)) add_period_layer(ps.items[i].id);
                if (ImGui::IsItemHovered()) timeline_periods_highlight(ps.items[i].id);
            }
            ImGui::EndMenu();
        }
        if (ps.count == 0) ImGui::SetItemTooltip("Mark periods in the timeline (shift-drag) to correlate the series over them");
        ImGui::TextDisabled("New layers start with the selected layer's series");
    }

    void draw_levels(Layer& l) {
        prop::row("levels", "The contours enclose these fractions of the samples (highest density regions)");
        density_map_levels_editor(l.levels, &l.num_levels);
    }

    // The series of an axis of a layer: listed, removed, and dropped on from the timeline
    void draw_series_list(Layer& l, int axis) {
        ImGui::PushID(axis);
        const int n = axis == 0 ? l.num_x : l.num_y;
        const SlotInfo* infos = axis == 0 ? l.x_info : l.y_info;
        const float line_h = ImGui::GetFrameHeight();
        const float w = prop::control_width();
        const float h = MAX(1, n) * (line_h + 2.0f) + 4.0f;
        const ImVec2 p0 = ImGui::GetCursorScreenPos();

        ImGui::BeginGroup();
        int remove = -1;
        if (n == 0) {
            ImGui::Dummy(ImVec2(w, 2));
            ImGui::TextDisabled("  drop series here");
        }
        for (int i = 0; i < n; ++i) {
            ImGui::PushID(i);
            const SlotInfo& info = infos[i];
            char name[128];
            slot_label(name, sizeof(name), l, axis, i);
            ImGui::ColorButton("##c", info.color, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(line_h * 0.5f, line_h * 0.5f));
            ImGui::SameLine();
            const float bw = line_h;
            const ImVec2 tp = ImGui::GetCursorScreenPos();
            ImGui::PushClipRect(tp, ImVec2(p0.x + w - bw - 4, tp.y + line_h), true);
            if (!info.ok) ImGui::TextColored(error_color, "%s", name);
            else if (info.unit[0]) ImGui::Text("%s  \xc2\xb7  %s", name, info.unit);
            else ImGui::TextUnformatted(name);
            ImGui::PopClipRect();
            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                ImGui::TextUnformatted(name);
                if (!info.ok) ImGui::TextColored(error_color, "Not available");
                else if (info.dim > 1) ImGui::TextDisabled("%d members", info.dim);
                ImGui::EndTooltip();
            }
            ImGui::SameLine();
            ImGui::SetCursorScreenPos(ImVec2(p0.x + w - bw, ImGui::GetCursorScreenPos().y));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
            if (ImGui::Button(ICON_FA_XMARK "##rm", ImVec2(bw, 0))) remove = i;
            ImGui::PopStyleColor();
            ImGui::SetItemTooltip("Remove");
            ImGui::PopID();
        }
        ImGui::EndGroup();

        const ImRect r(p0, ImVec2(p0.x + w, MAX(p0.y + h, ImGui::GetItemRectMax().y)));
        ImGui::GetWindowDrawList()->AddRect(r.Min, r.Max, ImGui::GetColorU32(ImGuiCol_Border), 3.0f);
        if (ImGui::BeginDragDropTargetCustom(r, ImGui::GetID("##drop"))) {
            const ImGuiPayload* p = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND, ImGuiDragDropFlags_AcceptBeforeDelivery);
            if (!p) p = ImGui::AcceptDragDropPayload(DISTRIBUTION_SERIES_DND, ImGuiDragDropFlags_AcceptBeforeDelivery);
            if (const SeriesDragPayload* sp = series_payload(p)) {
                const SeriesKey key = sp->key;
                const char* reason = nullptr;
                const bool ok = series_droppable(key, &reason);
                if (!ok && reason) ImGui::SetTooltip("%s", reason);
                if (p->IsDelivery() && ok) drop_series_on(l, axis, key);
            }
            ImGui::EndDragDropTarget();
        }
        if (remove >= 0) remove_series(l, axis, remove);
        ImGui::PopID();
    }

    void draw_inspector(Layer& l) {
        prop::begin("correlation");
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

        prop::row("X", "Series along X: one, or as many as Y to pair them in order");
        draw_series_list(l, 0);
        prop::row("Y", "Series along Y: one, or as many as X to pair them in order");
        draw_series_list(l, 1);
        prop::row("");
        if (ImGui::Button(ICON_FA_RIGHT_LEFT " Swap X and Y")) swap_axes(l);
        if (layer_shape(l) == LayerShape_Mismatch) {
            prop::hint("One series on an axis against several on the other, or as many on both (paired in order)");
        }
        if (l.members_truncated) {
            prop::hint("Populations of different sizes: paired up to the smaller");
        }
        if (l.binned) {
            prop::hint("Distributions: correlated bin by bin over their common range, of the frames below");
        } else {
            bool periodic = false;
            for (int i = 0; i < l.num_x; ++i) periodic |= series_periodic_valid(l.x_info[i].periodic);
            for (int i = 0; i < l.num_y; ++i) periodic |= series_periodic_valid(l.y_info[i].periodic);
            if (periodic) prop::hint("Periodic values are cut at their seam; their correlation is the circular one");
        }

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
                if (ps.count == 0) ImGui::TextDisabled("Mark periods in the timeline to use them here");
                ImGui::EndCombo();
            }
            if (l.frames != LayerFrames_Current && run_num_frames(app) < 2 && l.style != LayerStyle_Points) {
                prop::hint("No trajectory to find a density over: drawn as points");
            }
            if (l.frames == LayerFrames_Current && l.style != LayerStyle_Points) {
                prop::hint("One frame: drawn as points");
            }
        }

        prop::row("style");
        {
            int s = (int)l.style;
            static const int values[] = { LayerStyle_Points, LayerStyle_Contours, LayerStyle_Heatmap };
            density_map_segmented("##style", style_names, values, LayerStyle_Count, &s, nullptr, true);
            l.style = (LayerStyle)s;
        }

        // A heat map is coloured by its colour map, points by time too
        if (l.style != LayerStyle_Heatmap && !(layer_draws_points(l) && l.point_coloring == PointColoring_Time)) {
            prop::row("color");
            prop::color("##color", l.color.elem, true);
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
            draw_levels(l);
        }
        if (l.style == LayerStyle_Heatmap) {
            prop::row("colormap");
            ImPlot::ColormapSelection("##cmap", &l.colormap);
            prop::row("log scale", "Four decades below the maximum");
            ImGui::Checkbox("##log", &l.log_scale);
        }
        if (l.style != LayerStyle_Points) {
            prop::row("bandwidth", "The standard deviation of the gaussian the density is smoothed with, in % of the extent of the samples");
            ImGui::SliderFloat("##bandwidth", &l.bandwidth, 0.5f, 20.0f, "%.1f %%", ImGuiSliderFlags_AlwaysClamp);
        }
        if (layer_draws_points(l)) {
            prop::row("point size");
            ImGui::SliderFloat("##point_size", &l.point_size, 1.0f, 12.0f, "%.0f px");
            prop::row("color by", "Layer: its colour. Pair: the series that varies between pairs, or the members of a population. Time: by the colour map");
            int c = (int)l.point_coloring;
            if (ImGui::Combo("##coloring", &c, point_coloring_names, PointColoring_Count)) l.point_coloring = (PointColoring)c;
            const bool members = l.num_pairs == 1 && l.pairs[0].members > 1;
            if (l.point_coloring == PointColoring_Time || (l.point_coloring == PointColoring_Pair && members)) {
                prop::row("colormap");
                ImPlot::ColormapSelection("##cmap_points", &l.colormap);
            }
            if (l.point_stride > 1) {
                prop::hint("Every %u samples, to keep to %u points", l.point_stride, (uint32_t)MAX_POINTS);
            }
        }
        prop::row("fit", "The least squares line of each pair");
        ImGui::Checkbox("##fit", &l.show_fit);

        ImGui::PopID();
        prop::end();

        if (l.data_valid && l.num_pairs > 0) {
            ImGui::Spacing();
            draw_pair_table(l);
        } else if (l.data_error[0]) {
            ImGui::Spacing();
            ImGui::TextColored(error_color, "%s", l.data_error);
        }
    }

    void draw_sidebar() {
        // The header: the title and adding
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("LAYERS");
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - ImGui::GetFrameHeight()));
        if (ImGui::Button(ICON_FA_PLUS "##add", ImVec2(ImGui::GetFrameHeight(), 0))) ImGui::OpenPopup("##add_menu");
        ImGui::SetItemTooltip("Add a layer");
        if (ImGui::BeginPopup("##add_menu")) {
            draw_add_menu();
            ImGui::EndPopup();
        }

        SidebarActions act;
        {
            RepresentationNode* nodes[MAX_LAYERS];
            const size_t n = tree_nodes(nodes);
            for (int i = 0; i < num_layers; ++i) draw_layer_row(layers[i], nodes, n, &act);
        }
        if (num_layers == 0) {
            ImGui::TextDisabled("No layers: add one with +");
        }
        // Below the list: the end of it
        {
            RepresentationNode* nodes[MAX_LAYERS];
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
        if (act.remove) remove_layer(act.remove);

        if (Layer* l = find_layer(selected_layer)) {
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();
            draw_inspector(*l);
        }

        // Delete takes the selected layer away
        if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Delete)) {
            if (selected_layer) remove_layer(selected_layer);
        }
        if (tree_dirty) apply_tree();
    }

    void draw() {
        highlights_now = false;
        draw_window();
        if (app && owns_highlight && !highlights_now) md_bitmask_clear(&app->selection.highlight_mask);
        owns_highlight = highlights_now;
    }

    void draw_window() {
        if (!show_window || !app) return;

        ImGui::SetNextWindowSize(ImVec2(760, 480), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Correlation", &show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar)) {
            draw_toolbar();
            if (show_sidebar) {
                ImGui::BeginChild("##sidebar", ImVec2(320, 0), ImGuiChildFlags_ResizeX | ImGuiChildFlags_Borders);
                draw_sidebar();
                ImGui::EndChild();
                ImGui::SameLine();
            }
            ImGui::BeginGroup();
            const float avail = ImGui::GetContentRegionAvail().y;
            const float split = 6.0f;
            const float lag_h = show_lag ? MAX(80.0f, avail * lag_frac) : 0.0f;
            ImGui::BeginChild("##plot", ImVec2(0, show_lag ? MAX(80.0f, avail - lag_h - split) : 0.0f), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
            draw_plot();
            ImGui::EndChild();
            if (show_lag) {
                // Dragged to share the height between the two
                ImGui::InvisibleButton("##split", ImVec2(-1, split));
                if (ImGui::IsItemHovered() || ImGui::IsItemActive()) ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeNS);
                if (ImGui::IsItemActive() && avail > 0.0f) lag_frac = CLAMP(lag_frac - ImGui::GetIO().MouseDelta.y / avail, 0.1f, 0.8f);
                ImGui::BeginChild("##lag", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
                draw_lag_plot();
                ImGui::EndChild();
            }
            ImGui::EndGroup();
        }
        ImGui::End();
    }
};

static Correlation instance = {};

}  // namespace correlation
