#pragma once

// What the plotting windows (Timelines, Distributions), the property export and the density
// volume show: series of values, named by the PATH of an attribute.
//
// A series is a script property ("script/<ident>" in the full or the filtered evaluation's table)
// or a quantity loaded along the run ("<run>/edr/potential", "<run>/xvg/<file>/<column>") in the
// system's table. It is held by path, never by pointer or index: both tables are rebuilt under it -
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
#define SERIES_PATH_CAP              192

// ImGui drag and drop payload types, both carrying a SeriesDragPayload
#define TIMELINE_SERIES_DND     "TIMELINE_SERIES"
#define DISTRIBUTION_SERIES_DND "DISTRIBUTION_SERIES"
#define VOLUME_SERIES_DND       "VOLUME_SERIES"
#define VECTOR_SERIES_DND       "VECTOR_SERIES"

struct ApplicationState;
struct md_attribute_t;
struct md_attributes_t;
struct md_script_vis_payload_o;

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

enum SeriesSource : uint8_t {
    SeriesSource_System = 0,        // mold.sys.attributes: loaded along the run
    SeriesSource_Script,            // the full script evaluation
    SeriesSource_ScriptFiltered,    // the script evaluated over the timeline filter's frames
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

struct SeriesKey {
    SeriesSource  source  = SeriesSource_System;
    SeriesVariant variant = SeriesVariant_Values;
    char path[SERIES_PATH_CAP] = "";
};

bool      series_key_equal(const SeriesKey& a, const SeriesKey& b);

// The source's name in a workspace ("system", "script", "script_filtered"), and back
const char* series_source_name(SeriesSource source);
bool        series_source_from_name(SeriesSource* out, str_t name);
SeriesKey series_key(SeriesSource source, str_t path, SeriesVariant variant = SeriesVariant_Values);

// The table a source reads from, NULL when it has none right now
const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source);

// "d1", "d1 (mean)", "d1 (filtered)", "Potential"
void   series_label(char* buf, size_t cap, const ApplicationState* app, const SeriesKey& key);
ImVec4 series_default_color(const ApplicationState* app, const SeriesKey& key);

// "script/<ident>" -> "<ident>" for a script series, empty otherwise
str_t series_script_ident(const SeriesKey& key);

// Call between a BeginDragDropSource and its EndDragDropSource: sets the payload and draws the
// preview. src_subplot is -1 unless the series is dragged out of a plot.
struct SeriesDragPayload {
    SeriesKey key;
    int src_subplot = -1;
    uint32_t src_uid = 0;       // The series dragged out of a plot, when it is one of several of the same key
};
void series_set_drag_payload(const char* dnd_type, const SeriesKey& key, int src_subplot, const char* label, ImVec4 color, uint32_t src_uid = 0);

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
    uint32_t uid = 0;           // Names this entry in its plot for as long as it lives: the same key can be there twice

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

    // Vectors: the members at the displayed frame, and layers of the density of every member over the frames
    bool            show_points = true;
    bool            points_all = false;         // every frame, not only the displayed one
    bool            iso_lines = true;           // contours, in the colour of the series
    bool            iso_levels = false;         // the bands between them, filled
    bool            density_map = false;        // a colormap of the density
    int             iso_count = 4;              // contours, at evenly spaced fractions of the peak density
    ImPlotColormap  density_colormap = ImPlotColormap_Viridis;
    float           density_alpha = 0.8f;
    float           smoothing = 2.0f;
};

struct PlotSubplot {
    PlotSeries series[PLOT_MAX_SERIES_PER_SUBPLOT];
    int count = 0;
};

int  plot_find_series(const PlotSubplot& sp, const SeriesKey& key);
int  plot_find_uid(const PlotSubplot& sp, uint32_t uid);
// Adds the series even when its key is already in the subplot, which is how a series is shown twice
int  plot_add_series_copy(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key);
// Adds the series with its default style and returns its index; the index it already has if it
// is there, -1 when the subplot is full.
int  plot_add_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key);
void plot_remove_series(PlotSubplot& sp, int index);
// Out of one subplot and into another, keeping its style
void plot_move_series(const ApplicationState* app, PlotSubplot& src, PlotSubplot& dst, const SeriesKey& key);
void plot_clear(PlotSubplot* subplots, int count);

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
    SeriesView_Volume       = 4,    // resident scalar script volume with reference-space metadata
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

    char  label[64] = "";
    char  plot_id[SERIES_PATH_CAP + 72] = "";   // label##identity, what ImPlot is handed
    char  unit_str[32] = "";                    // the display unit of the values
    char  script_ident[64] = "";                // a script property's identifier, empty otherwise
    const md_script_vis_payload_o* vis_payload = nullptr;
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

struct SeriesHistogramView {
    const float* bins = nullptr;    // dim * num_bins, member k's bins at k * num_bins, stored units
    int    num_bins = 0;
    int    dim = 1;
    double x_min = 0, x_max = 0;    // in the display unit
    double y_scale = 1.0;           // bins -> display unit
    double y_min = 0, y_max = 0;    // in the display unit

    char  label[64] = "";
    char  plot_id[SERIES_PATH_CAP + 72] = "";
    char  x_unit_str[32] = "";
    char  y_unit_str[32] = "";
    char  script_ident[64] = "";
    const md_script_vis_payload_o* vis_payload = nullptr;
};

// num_bins is what is asked for; a script distribution is evaluated at a resolution of its own and
// gives at most that. False, with labels filled in, when there is nothing to bin.
bool series_resolve_histogram(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins);

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
    const md_attribute_t* attr = nullptr;
    uint64_t version = 0;
    char  label[64] = "";
    char  script_ident[64] = "";
    const md_script_vis_payload_o* vis_payload = nullptr;
};
bool series_resolve_volume(SeriesVolumeView* out, ApplicationState* app, const SeriesKey& key);

// ## Editing series

// The style of one series: colour or colormap, and which of its 'dim' population members are drawn.
// Hovering a member sets *hovered_pop_idx and highlights it in the 3D view.
void plot_series_style_popup(ApplicationState* app, PlotSeries& s, int dim, const char* script_ident, const md_script_vis_payload_o* vis_payload, int* hovered_pop_idx);

// One entry of a list of series: a click toggles it in 'first', a drag carries it (dnd_type)
void plot_series_list_item(ApplicationState* app, const char* dnd_type, PlotSubplot& first, const SeriesKey& key);

// ## Vectors

// The style of a series of vectors: points, the layers of the density, which frames. Call inside a popup.
void plot_vector_style_popup(PlotSeries& s, bool point_color = true, bool all_frames = false);

// A script property with a vector per member (ramachandran(), shapespace()): num_frames * count
// vectors of 'components' floats, the frame first
struct SeriesVectorView {
    const float* data = nullptr;
    uint32_t num_frames = 0;
    uint32_t count = 0;         // members per frame
    uint32_t components = 0;
    uint64_t version = 0;
    char  label[64] = "";
    char  script_ident[64] = "";
    const md_script_vis_payload_o* vis_payload = nullptr;
};
// False unless the key is a resident script property of that many components
bool series_resolve_vector(SeriesVectorView* out, ApplicationState* app, const SeriesKey& key, uint32_t components);

// A window of series of vectors of that many components shares these: the menus, the bar of series under
// its toolbar, which also takes the properties dropped on the window, and the removal of a series dragged out.

// Adds a series with the style of a vector, the same property again when allow_duplicate. Returns its index.
int  plot_add_vector_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key, bool allow_duplicate);
// The properties which can be shown: a click toggles one, a drag carries it
void series_vector_properties_menu(ApplicationState* app, PlotSubplot& sp, uint32_t components);
// The layers of every series. Call inside a menu.
void series_vector_layers_menu(ApplicationState* app, PlotSubplot& sp, uint32_t components);
// A chip for each series: right click edits it, drag moves it out of the window to remove it. Changes are applied before it returns.
// members: the vectors of a series are the members of a population, each with a colour of its own.
void series_vector_bar(ApplicationState* app, PlotSubplot& sp, uint32_t components, bool members = true);
// Call at the end of the window: removes the series dragged out of it and dropped where nothing takes it
void series_vector_drag_out(PlotSubplot& sp);

// ## Listing script properties

// Calls fn(key) for each property of the compiled script whose kind matches any of the flags
// (MD_SCRIPT_PROPERTY_FLAG_*) and that the source's evaluation has published, in script order.
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
