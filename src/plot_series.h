#pragma once

// What the plotting windows (Timelines, Distributions), the property export and the density
// volume show: series of values, named by the PATH of an attribute.
//
// A series is a script property ("script/<ident>" in the evaluation's table), over every frame or
// over the frames of a timeline period, or a quantity loaded along the run ("<run>/edr/potential",
// "<run>/xvg/<file>/<column>") in the system's table. It is held by path, never by pointer or index: both tables are rebuilt under it -
// a script recompiles, an energy file is loaded again - and a path survives that where a pointer
// does not. So the views own paths, resolve them when they draw, and a path with nothing behind it
// is shown as an empty legend entry that can still be removed.
//
// How a series is drawn is a view of the same attribute:
//   - over time (Timelines): on the attribute's own frame axis (md_attributes_axis), converted to
//     the timeline's time unit. An energy file written every 10 steps against coordinates every
//     1000 is drawn at its own resolution, not resampled onto the trajectory's frames.
//   - as a distribution (Distributions): a temporal attribute binned over its frames (the frames the
//     evaluation covered, for a script property), or a script distribution as it was evaluated.
// Derived data - converted axes, histograms - is cached by path and attribute version.

#include <core/md_str.h>
#include <core/md_unit.h>
#include <md_script.h>

// The rest of the application uses ImGui's vector operators, which have to be asked for before
// imgui.h is first included; this header may be that first include
#ifndef IMGUI_DEFINE_MATH_OPERATORS
#define IMGUI_DEFINE_MATH_OPERATORS
#endif
#include <imgui.h>
#include <implot.h>

#include <bitset>
#include <stdint.h>
#include <stddef.h>

#define PLOT_MAX_SUBPLOTS            10
#define PLOT_MAX_SERIES_PER_SUBPLOT  32
#define PLOT_MAX_POPULATION          256
#define PLOT_MAX_PERIODS             8
#define SERIES_PATH_CAP              192

// ImGui drag and drop payload types, all carrying a SeriesDragPayload
#define TIMELINE_SERIES_DND     "TIMELINE_SERIES"
#define DISTRIBUTION_SERIES_DND "DISTRIBUTION_SERIES"
#define VOLUME_SERIES_DND       "VOLUME_SERIES"     // a script volume, for the density volume

struct ApplicationState;
struct md_attribute_t;
struct md_attributes_t;

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

enum SeriesSource : uint8_t {
    SeriesSource_System = 0,        // mold.sys.attributes: loaded along the run
    SeriesSource_Script,            // the script evaluation
    SeriesSource_Count
};

// What of the attribute is drawn. The summaries read the siblings a population carries
// ("<path>/mean", "<path>/variance", "<path>/extent"), so they apply to whatever has them.
enum SeriesVariant : uint8_t {
    SeriesVariant_Values = 0,
    SeriesVariant_Mean,
    SeriesVariant_Sigma,        // mean +- one standard deviation, as a band
    SeriesVariant_Extent,       // min .. max, as a band
    SeriesVariant_Aggregate,    // a distribution over every member of a population at once
    SeriesVariant_Count
};

// What the timeline filter was shown over, in a workspace from before there were periods: the period the filter followed,
// once the workspace is read (series_legacy_period)
#define SERIES_PERIOD_LEGACY_FILTER UINT32_MAX

struct SeriesKey {
    SeriesSource  source  = SeriesSource_System;
    SeriesVariant variant = SeriesVariant_Values;
    char path[SERIES_PATH_CAP] = "";
    // The timeline period a script property is taken over, by id, 0 for every frame. A temporal property is binned over
    // the period's frames, a distribution or a volume read from the evaluation's slot of it (SCRIPT_SLOT_PERIOD).
    uint32_t period = 0;
};

bool      series_key_equal(const SeriesKey& a, const SeriesKey& b);

// The source's name in a workspace ("system", "script"), and back. "script_filtered", of a workspace from before there
// were periods, is the script over SERIES_PERIOD_LEGACY_FILTER, which is written to period when it is given.
const char* series_source_name(SeriesSource source);
bool        series_source_from_name(SeriesSource* out, str_t name, uint32_t* period = nullptr);
SeriesKey series_key(SeriesSource source, str_t path, SeriesVariant variant = SeriesVariant_Values, uint32_t period = 0);

// The period a workspace's SERIES_PERIOD_LEGACY_FILTER stands for, once its periods are read: the period the filter
// followed, or every frame. Any other period as it is.
uint32_t series_legacy_period(const ApplicationState* app, uint32_t period);

// The table a source reads from, NULL when it has none right now. The script's is the evaluation's: over a period, a
// temporal property is taken over the period's frames (series_frame_mask), and a distribution or a volume read from the
// evaluation's slot of it over them, not from the table.
const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source);
// The table which holds path for source
const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source, str_t path);
// The frames the temporal script property of key is taken over: those it has evaluated, of its period's, written to
// out, allocated with alloc. Its own: a property taken over from a previous evaluation holds its frames while one added
// since is still being evaluated. A copy taken under the evaluation's lock, so it can be read while the evaluation
// runs. False for any other source, which is taken over every frame.
bool series_frame_mask(struct md_bitset_t* out, const ApplicationState* app, const SeriesKey& key, struct md_allocator_i* alloc);

// The name of a period, "missing period" when it is gone
const char* series_period_name(const ApplicationState* app, uint32_t period);

// "d1", "d1 (mean)", "r (Production)", "Potential"
void   series_label(char* buf, size_t cap, const ApplicationState* app, const SeriesKey& key);
ImVec4 series_default_color(const ApplicationState* app, const SeriesKey& key);

// "script/<ident>" -> "<ident>" for a script series, empty otherwise
str_t series_script_ident(const SeriesKey& key);

// How what a script series shows stands against the script as it is now in the editor, and how far its evaluation has
// come. Anything else is always current.
enum SeriesStatus : uint8_t {
    SeriesStatus_Current = 0,   // what the script says, over every frame it is shown over
    SeriesStatus_Partial,       // what the script says, over some of those frames so far
    SeriesStatus_Outdated,      // as evaluated from an earlier version of the script: what it is computed from changed
    SeriesStatus_Removed,       // as evaluated, but no longer in the script
};

struct SeriesState {
    SeriesStatus status = SeriesStatus_Current;
    uint32_t num_frames = 0;    // the frames it holds of those it is shown over
    uint32_t num_wanted = 0;    // the frames it is shown over
    uint32_t num_extra = 0;     // frames it still holds and is no longer shown over: a slot moved, still taking them out
};

SeriesState series_state(const ApplicationState* app, const SeriesKey& key);
// " · outdated", " · 380/500": what is appended to the label of a series in its state, empty when it is current
void series_state_suffix(char* buf, size_t cap, const SeriesState& state);
// A sentence on it for a tooltip, empty when it is current
void series_state_describe(char* buf, size_t cap, const SeriesState& state);

// Call between a BeginDragDropSource and its EndDragDropSource: sets the payload and draws the
// preview. src_subplot is -1 unless the series is dragged out of a plot.
struct SeriesDragPayload {
    SeriesKey key;
    int src_subplot = -1;
};
void series_set_drag_payload(const char* dnd_type, const SeriesKey& key, int src_subplot, const char* label, ImVec4 color);

// ## Subplots

enum PlotType : uint8_t {
    PlotType_Line = 0,
    PlotType_Scatter,
    PlotType_Area,          // a distribution filled down to zero
    PlotType_Bars,
    PlotType_Count
};

// One entry of a subplot: which series, and how it is drawn there. The same series in two
// subplots is two entries, each with its own style.
struct PlotSeries {
    SeriesKey key;

    ImVec4          color = {1,1,1,1};
    bool            use_colormap = false;       // one colour per population member
    ImPlotColormap  colormap = ImPlotColormap_Plasma;
    float           colormap_alpha = 1.0f;

    PlotType        plot_type = PlotType_Line;
    ImPlotMarker    marker = ImPlotMarker_Square;
    float           marker_size = 1.0f;
    float           bar_width = 1.0f;           // of a bin, for bars
    int             num_bins = 128;             // for a distribution

    // Which members of a population are drawn
    std::bitset<PLOT_MAX_POPULATION> population_mask = {};

    // Distributions: the periods it is shown over, one curve each, to compare the same property between them. A period
    // is a timeline period, by id, or 0 for the frames of the key: every frame, or its period's. A script property only;
    // anything else is shown over every frame.
    uint32_t periods[PLOT_MAX_PERIODS] = {0};
    int      num_periods = 1;
};

// Whether the series is shown over anything but the frames of its key
static inline bool plot_series_has_periods(const PlotSeries& s) {
    return s.num_periods > 1 || (s.num_periods == 1 && s.periods[0] != 0);
}
bool plot_series_has_period(const PlotSeries& s, uint32_t id);
// False when it has it already, or has no room for another
bool plot_series_add_period(PlotSeries& s, uint32_t id);
// One is always left: the last is not removed
void plot_series_remove_period(PlotSeries& s, uint32_t id);

struct PlotSubplot {
    PlotSeries series[PLOT_MAX_SERIES_PER_SUBPLOT];
    int count = 0;
};

int  plot_find_series(const PlotSubplot& sp, const SeriesKey& key);
// Adds the series with its default style and returns its index; the index it already has if it
// is there, -1 when the subplot is full.
int  plot_add_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key);
void plot_remove_series(PlotSubplot& sp, int index);
// Out of one subplot and into another, keeping its style
void plot_move_series(const ApplicationState* app, PlotSubplot& src, PlotSubplot& dst, const SeriesKey& key);
void plot_clear(PlotSubplot* subplots, int count);
// What a workspace from before there were periods showed over the timeline filter (SERIES_PERIOD_LEGACY_FILTER) over
// the period the filter followed, once its periods are read. A series that then is one shown already is left out.
void plot_resolve_legacy_periods(const ApplicationState* app, PlotSubplot* subplots, int count);

// Workspace: the subplot count and every entry with its style, as one section plus one section
// per entry. The layout is cleared by the caller before a workspace is read.
void plot_layout_serialize(viamd::serialization_state_t& state, str_t section, str_t series_section, const PlotSubplot* subplots, int num_subplots);
// Reads the section the state is at if it is one of the two. False otherwise.
bool plot_layout_deserialize(viamd::deserialization_state_t& state, str_t section, str_t series_section, const ApplicationState* app, PlotSubplot* subplots, int* num_subplots);

// ## Series loaded along the run
//
// Grouped by the file they came from: a group is any attribute group below the run with a
// "source" string beside its members ("<run>/edr", "<run>/xvg/energy"). Paths are views into the
// system's table.

size_t system_series_groups(str_t out_groups[], size_t cap, const ApplicationState* app);
str_t  system_series_group_source(const ApplicationState* app, str_t group);
// "edr  (md.edr)": the group below the run and the name of the file
void   system_series_group_label(char* buf, size_t cap, const ApplicationState* app, str_t group);
// The plottable members of a group: temporal, numeric, resident, not an axis. Returns the total.
size_t system_series_members(str_t out_paths[], size_t cap, const ApplicationState* app, str_t group);

// ## What a series can be shown as

enum : uint32_t {
    SeriesView_Timeline     = 1,    // over time: sampled along a frame axis the timeline can place
    SeriesView_Distribution = 2,    // binned over its frames, or a script distribution as evaluated
};

// Which views can draw the attribute at key.path (the Values variant of it). When none can, why
// not, as a static string for a tooltip. Cheap enough to ask for every attribute in a table.
uint32_t series_views(const ApplicationState* app, const SeriesKey& key, const char** out_reason = nullptr);

// ## Over time

// Resolved for the frame it was resolved in: what a plot reads.
struct SeriesTemporalView {
    const float* x = nullptr;       // num_samples, in the timeline's time unit
    const float* y = nullptr;       // num_samples * stride, in the attribute's stored unit
    const float* y2 = nullptr;      // Sigma: the variance beside the mean in y
    int   num_samples = 0;
    int   stride = 1;               // values per sample in y
    int   dim = 1;                  // population members drawn (1 for the summaries)
    bool  band = false;             // drawn as an area between series_temporal_band's lo and hi
    float y_scale = 1.0f;           // stored unit -> display unit

    char  label[96] = "";                       // with the state appended (series_state_suffix)
    char  plot_id[SERIES_PATH_CAP + 104] = "";  // label###identity, what ImPlot is handed: the label may change
    char  unit_str[32] = "";                    // the display unit of the values
    char  script_ident[64] = "";                // a script property's identifier, empty otherwise
    md_script_vis_ref_t vis_ref = {};
    SeriesState state = {};
};

// False when nothing is published at the path (yet), or it cannot be placed on the timeline: not
// temporal, not numeric, virtual, or on an axis whose unit does not relate to the timeline's. The
// labels and the script identity are filled in either way.
bool series_resolve_temporal(SeriesTemporalView* out, ApplicationState* app, const SeriesKey& key);

static inline ImPlotPoint series_temporal_point(const SeriesTemporalView& v, int i, int k) {
    return ImPlotPoint(v.x[i], v.y[(size_t)i * v.stride + k] * v.y_scale);
}
void   series_temporal_band(const SeriesTemporalView& v, int i, double* lo, double* hi);
// Fractional sample index of x along the series' axis, clamped to its range
double series_temporal_index_at(const SeriesTemporalView& v, double x);
// The value(s) at sample i as text, in the display unit (the unit itself is not appended)
int    series_temporal_print_value(char* buf, size_t cap, const SeriesTemporalView& v, int i, int k);

// ## As a distribution

// One period of a series over periods (series_resolve_histogram_periods)
struct SeriesPeriodInfo {
    uint32_t id = 0;            // the timeline period, 0 for the frames of the series' key
    bool     valid = false;     // there is something to show over it
    bool     own_color = false; // in the colour of its period; otherwise in the series' colour
    ImVec4   color = {};
    char     name[64] = "";
    uint32_t num_frames = 0;    // the frames it holds so far, of
    uint32_t num_wanted = 0;    // the frames it covers
};

struct SeriesHistogramView {
    const float* bins = nullptr;    // dim * num_bins, member k's bins at k * num_bins, stored units
    int    num_bins = 0;
    int    dim = 1;
    double x_min = 0, x_max = 0;    // in the display unit
    double y_scale = 1.0;           // bins -> display unit
    double y_min = 0, y_max = 0;    // in the display unit

    char  label[96] = "";                       // with the state appended (series_state_suffix)
    char  plot_id[SERIES_PATH_CAP + 104] = "";
    char  x_unit_str[32] = "";
    char  y_unit_str[32] = "";
    char  script_ident[64] = "";
    md_script_vis_ref_t vis_ref = {};
    SeriesState state = {};

    // Over periods: member k is period k. 0 when the members are those of a population.
    int num_periods = 0;
    SeriesPeriodInfo periods[PLOT_MAX_PERIODS];
};

// num_bins is what is asked for; a script distribution is evaluated at a resolution of its own and
// gives at most that. False, with labels filled in, when there is nothing to bin.
// period: a timeline period (see PlotSeries::periods), over whose frames it is binned, or 0 for the frames of the key.
// A script distribution over a period is read from the evaluation's slot of it (SCRIPT_SLOT_PERIOD), which is asked for
// here: nothing is evaluated for it unless it is not kept per frame, and then only when Evaluate is pressed.
bool series_resolve_histogram(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins, uint32_t period = 0);

// Over several periods, one member each, in their order: all binned alike, so that they compare. A population is
// pooled (the aggregate) over periods. Without periods, the same as series_resolve_histogram.
bool series_resolve_histogram_periods(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins,
                                      const uint32_t* periods, int num_periods);

// The frames of a period now, [beg, end[: of the period, or for 0 of the key's own period. False when it covers no
// frames, or the period is gone.
bool series_period_frames(uint32_t* beg, uint32_t* end, const ApplicationState* app, const SeriesKey& key, uint32_t period);

static inline double series_histogram_bin_width(const SeriesHistogramView& v) {
    return v.num_bins > 0 ? (v.x_max - v.x_min) / v.num_bins : 0.0;
}
static inline ImPlotPoint series_histogram_point(const SeriesHistogramView& v, int i, int k) {
    const double w = series_histogram_bin_width(v);
    return ImPlotPoint(v.x_min + (i + 0.5) * w, v.bins[(size_t)k * v.num_bins + i] * v.y_scale);
}
// Fractional bin index of x, unclamped
static inline double series_histogram_index_at(const SeriesHistogramView& v, double x) {
    const double w = series_histogram_bin_width(v);
    return w > 0.0 ? (x - v.x_min) / w - 0.5 : 0.0;
}

// ## Volumes

// A script property evaluated as a volume (rank 3), for the density volume and the export
struct SeriesVolumeView {
    const float* data = nullptr;    // dim[0] * dim[1] * dim[2] values, as the script publishes them. Valid this frame.
    uint32_t dim[3] = {};
    uint64_t version = 0;           // of the data: changes whenever it does
    char  label[64] = "";
    char  script_ident[64] = "";
    md_script_vis_ref_t vis_ref = {};
    SeriesState state = {};
};
bool series_resolve_volume(SeriesVolumeView* out, ApplicationState* app, const SeriesKey& key);

// ## Listing script properties

// Calls fn(key) for each property of the compiled script whose kind matches any of the flags
// (MD_SCRIPT_PROPERTY_FLAG_*) and that the source's evaluation has published, in script order. Over every frame: a
// period is set on the key by whoever wants one.
template <typename Fn>
void series_for_each_script_property(const ApplicationState* app, SeriesSource source, uint32_t kind_flags, Fn fn);

// ## Caches

// Drops cached conversions and histograms nothing has asked for in a while, and all on shutdown
void series_cache_gc(ApplicationState* app);
void series_cache_free(ApplicationState* app);

// Implementation of the template above; needs the script and the application
size_t series_script_property_count(const ApplicationState* app);
bool   series_script_property(SeriesKey* out, const ApplicationState* app, SeriesSource source, size_t index, uint32_t kind_flags);

template <typename Fn>
void series_for_each_script_property(const ApplicationState* app, SeriesSource source, uint32_t kind_flags, Fn fn) {
    const size_t num = series_script_property_count(app);
    for (size_t i = 0; i < num; ++i) {
        SeriesKey key;
        if (series_script_property(&key, app, source, i, kind_flags)) {
            fn(key);
        }
    }
}
