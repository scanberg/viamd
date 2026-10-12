#include <plot_series.h>
#include <imgui_internal.h>
#include <implot_widgets.h>
#include <imgui_widgets.h>

#include <viamd.h>
#include <display_units.h>
#include <serialization_utils.h>

#include <md_system.h>
#include <md_script.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_str.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// ## Keys and names

bool series_key_equal(const SeriesKey& a, const SeriesKey& b) {
    return a.source == b.source && a.variant == b.variant && strcmp(a.path, b.path) == 0;
}

SeriesKey series_key(SeriesSource source, str_t path, SeriesVariant variant) {
    SeriesKey key;
    key.source  = source;
    key.variant = variant;
    // A path that does not fit is left empty rather than cut: a truncated path names something else
    if (path.len < sizeof(key.path)) {
        str_copy_to_char_buf(key.path, sizeof(key.path), path);
    }
    return key;
}

// The child of the series' path each variant is read from; empty for the series itself.
static const char* variant_leaf[SeriesVariant_Count] = {
    "",
    "mean",
    "variance",
    "extent",
    "",
};

static const char* variant_label[SeriesVariant_Count] = {
    "",
    " (mean)",
    " (mean+-sd)",
    " (min/max)",
    " (agg)",
};

static const char* variant_name[SeriesVariant_Count] = {
    "values",
    "mean",
    "sigma",
    "extent",
    "aggregate",
};

static const char* source_name[SeriesSource_Count] = {
    "system",
    "script",
    "script_filtered",
};

const char* series_source_name(SeriesSource source) {
    return source_name[source < SeriesSource_Count ? source : 0];
}

bool series_source_from_name(SeriesSource* out, str_t name) {
    for (int k = 0; k < SeriesSource_Count; ++k) {
        if (str_eq_cstr(name, source_name[k])) {
            *out = (SeriesSource)k;
            return true;
        }
    }
    return false;
}

static const char* plot_type_name[PlotType_Count] = {
    "line",
    "scatter",
    "area",
    "bars",
};

static const str_t script_prefix = STR_INIT("script/");

static bool is_script(SeriesSource source) {
    return source == SeriesSource_Script || source == SeriesSource_ScriptFiltered;
}

static const md_script_eval_t* series_eval(const ApplicationState* app, SeriesSource source) {
    switch (source) {
    case SeriesSource_Script:         return app->script.full_eval;
    case SeriesSource_ScriptFiltered: return app->script.filt_eval;
    default: return nullptr;
    }
}

const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source) {
    ASSERT(app);
    if (source == SeriesSource_System) {
        return &app->mold.sys.attributes;
    }
    const md_script_eval_t* eval = series_eval(app, source);
    return eval ? md_script_eval_attributes(eval) : nullptr;
}

str_t series_script_ident(const SeriesKey& key) {
    str_t path = str_from_cstr(key.path);
    if (is_script(key.source) && str_begins_with(path, script_prefix)) {
        return str_substr(path, script_prefix.len);
    }
    return {};
}

static str_t path_leaf(str_t path) {
    size_t loc;
    if (str_rfind_char(&loc, path, '/')) {
        return str_substr(path, loc + 1);
    }
    return path;
}

void series_label(char* buf, size_t cap, const ApplicationState* app, const SeriesKey& key) {
    ASSERT(buf && cap > 0);
    str_t base = {};
    if (is_script(key.source)) {
        base = series_script_ident(key);
    } else {
        const md_attributes_t* table = series_table(app, key.source);
        const md_attribute_t*  attr  = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
        base = (attr && !str_empty(attr->label)) ? attr->label : path_leaf(str_from_cstr(key.path));
    }
    const int v = key.variant < SeriesVariant_Count ? key.variant : 0;
    const char* filtered = key.source == SeriesSource_ScriptFiltered ? " (filtered)" : "";
    snprintf(buf, cap, STR_FMT "%s%s", STR_ARG(base), variant_label[v], filtered);
}

ImVec4 series_default_color(const ApplicationState* app, const SeriesKey& key) {
    size_t idx = 0;
    if (is_script(key.source)) {
        // The property's place in the script, the same in every view of it
        const str_t ident = series_script_ident(key);
        const md_script_ir_t* ir = app->script.eval_ir;
        const size_t num = ir ? md_script_ir_property_count(ir) : 0;
        const str_t* names = ir ? md_script_ir_property_names(ir) : nullptr;
        for (size_t i = 0; i < num; ++i) {
            if (str_eq(names[i], ident)) { idx = i; break; }
        }
    } else {
        idx = (size_t)md_hash64_str(str_from_cstr(key.path), 0);
    }
    ImVec4 color = ImGui::ColorConvertU32ToFloat4(PROPERTY_COLORS[idx % ARRAY_SIZE(PROPERTY_COLORS)]);
    if (key.variant == SeriesVariant_Sigma)  color.w *= 0.4f;
    if (key.variant == SeriesVariant_Extent) color.w *= 0.2f;
    return color;
}

void series_set_drag_payload(const char* dnd_type, const SeriesKey& key, int src_subplot, const char* label, ImVec4 color, uint32_t src_uid) {
    SeriesDragPayload payload;
    payload.key = key;
    payload.src_subplot = src_subplot;
    payload.src_uid = src_uid;
    ImGui::SetDragDropPayload(dnd_type, &payload, sizeof(payload));
    ImPlot::ItemIcon(color);
    ImGui::SameLine();
    ImGui::TextUnformatted(label);
}

// ## Attributes

// Something a view can be drawn from: numbers, and resident. A virtual attribute (a trajectory's
// coordinates) would be produced whole to draw it, which is not something anyone means by
// dragging it into a plot.
static bool attribute_numeric_resident(const md_attribute_t* attr) {
    if (!attr) return false;
    if (attr->format.rank < 1 || attr->format.shape[0] == 0) return false;
    return md_attribute_type_is_numeric(attr->format.type) && !md_attribute_is_virtual(attr);
}

static bool attribute_temporal(const md_attribute_t* attr) {
    return attribute_numeric_resident(attr) && (attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
}

static md_script_property_flags_t script_property_kind(const ApplicationState* app, const SeriesKey& key) {
    const md_script_ir_t* ir = app->script.eval_ir;
    return (ir && is_script(key.source)) ? md_script_ir_property_flags(ir, series_script_ident(key)) : MD_SCRIPT_PROPERTY_FLAG_NONE;
}

// Labels and script identity, which a view has whether or not the series resolves
template <typename View>
static void fill_names(View* out, const ApplicationState* app, const SeriesKey& key) {
    series_label(out->label, sizeof(out->label), app, key);
    snprintf(out->plot_id, sizeof(out->plot_id), "%s##%d:%d:%s", out->label, (int)key.source, (int)key.variant, key.path);
    if (is_script(key.source)) {
        const str_t ident = series_script_ident(key);
        str_copy_to_char_buf(out->script_ident, sizeof(out->script_ident), ident);
        out->vis_payload = app->script.eval_ir ? md_script_ir_property_vis_payload(app->script.eval_ir, ident) : nullptr;
    }
}

// ## Caches
//
// What is derived from an attribute is kept until the attribute changes. An attribute's version
// only orders the changes within one table, and a table can be freed and another allocated in its
// place (a new dataset, a recompiled script), so the caches are dropped whole when that happens
// (series_cache_free) rather than trusted to tell.

struct ConvertEntry {
    uint64_t               key = 0;
    const md_attributes_t* table = nullptr;
    uint64_t               version = 0;
    bool                   valid = false;
    md_array(float)        data = nullptr;
    int                    last_used = 0;
};

struct HistogramEntry {
    uint64_t               key = 0;
    const md_attributes_t* table = nullptr;
    uint64_t               version = 0;     // of the values, and of the weights where there are any
    bool                   valid = false;
    int                    num_bins = 0;
    int                    dim = 0;
    double                 x_min = 0, x_max = 0;    // stored unit
    float                  y_min = 0, y_max = 0;
    md_array(float)        bins = nullptr;
    int                    last_used = 0;
};

static md_array(ConvertEntry)   convert_cache = nullptr;
static md_array(HistogramEntry) histogram_cache = nullptr;

// A converted copy of the whole attribute as floats in dst_unit (as stored for md_unit_none())
static const float* cache_convert(ApplicationState* app, const md_attributes_t* table, const md_attribute_t* attr, md_unit_t dst_unit, size_t* out_count) {
    uint64_t key = md_hash64(&dst_unit.base.raw_bits, sizeof(dst_unit.base.raw_bits), 0);
    key = md_hash64(&dst_unit.mult, sizeof(dst_unit.mult), key);
    key = md_hash64_str(attr->path, key);

    const int frame = ImGui::GetFrameCount();
    ConvertEntry* entry = nullptr;
    for (size_t i = 0; i < md_array_size(convert_cache); ++i) {
        if (convert_cache[i].key == key && convert_cache[i].table == table) {
            entry = &convert_cache[i];
            break;
        }
    }
    if (!entry) {
        md_array_push(convert_cache, ConvertEntry{}, app->allocator.persistent);
        entry = &convert_cache[md_array_size(convert_cache) - 1];
        entry->key   = key;
        entry->table = table;
    } else if (entry->version == attr->version) {
        entry->last_used = frame;
        *out_count = md_array_size(entry->data);
        return entry->valid ? entry->data : nullptr;
    }

    entry->version   = attr->version;
    entry->last_used = frame;
    const size_t count = md_attribute_element_count(&attr->format);
    md_array_resize(entry->data, count, app->allocator.persistent);
    entry->valid = md_attribute_extract_f32(entry->data, count, attr, md_attribute_slice_all(), dst_unit) == count;
    if (!entry->valid) {
        md_array_free(entry->data, app->allocator.persistent);
        entry->data = nullptr;
    }
    *out_count = md_array_size(entry->data);
    return entry->valid ? entry->data : nullptr;
}

// The attribute's values as floats: in place when they are stored so (a script property, which
// is read while the evaluation fills it in), converted and cached otherwise
static const float* attribute_floats(ApplicationState* app, const md_attributes_t* table, const md_attribute_t* attr) {
    if (attr->format.type == MD_ATTRIBUTE_TYPE_F32) {
        return (const float*)attr->data;
    }
    size_t count = 0;
    const float* data = cache_convert(app, table, attr, md_unit_none(), &count);
    return (data && count == md_attribute_element_count(&attr->format)) ? data : nullptr;
}

void series_cache_gc(ApplicationState* app) {
    // Several seconds unused: nothing draws it any more
    const int frame = ImGui::GetFrameCount();
    for (size_t i = 0; i < md_array_size(convert_cache);) {
        if (frame - convert_cache[i].last_used > 600) {
            md_array_free(convert_cache[i].data, app->allocator.persistent);
            md_array_swap_back_and_pop(convert_cache, i);
        } else {
            ++i;
        }
    }
    for (size_t i = 0; i < md_array_size(histogram_cache);) {
        if (frame - histogram_cache[i].last_used > 600) {
            md_array_free(histogram_cache[i].bins, app->allocator.persistent);
            md_array_swap_back_and_pop(histogram_cache, i);
        } else {
            ++i;
        }
    }
}

void series_cache_free(ApplicationState* app) {
    for (size_t i = 0; i < md_array_size(convert_cache); ++i) {
        md_array_free(convert_cache[i].data, app->allocator.persistent);
    }
    md_array_free(convert_cache, app->allocator.persistent);
    convert_cache = nullptr;

    for (size_t i = 0; i < md_array_size(histogram_cache); ++i) {
        md_array_free(histogram_cache[i].bins, app->allocator.persistent);
    }
    md_array_free(histogram_cache, app->allocator.persistent);
    histogram_cache = nullptr;
}

// ## What a series can be shown as

uint32_t series_views(const ApplicationState* app, const SeriesKey& key, const char** out_reason) {
    const char* reason_buf = nullptr;
    const char** reason = out_reason ? out_reason : &reason_buf;
    *reason = nullptr;

    const md_attributes_t* table = series_table(app, key.source);
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr) {
        *reason = "Nothing is published at this path";
        return 0;
    }
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {
        *reason = "Text, not numbers";
        return 0;
    }
    if (is_script(key.source) && (script_property_kind(app, key) & MD_SCRIPT_PROPERTY_FLAG_VOLUME)) {
        if (attr->format.type == MD_ATTRIBUTE_TYPE_F32 && attr->format.rank == 3 &&
            attr->format.components == 1 && attribute_numeric_resident(attr) &&
            attr->format.shape[0] && attr->format.shape[1] && attr->format.shape[2]) {
            return SeriesView_Volume;
        }
        *reason = "A volume must be a resident, nonempty scalar f32 grid";
        return 0;
    }

    const bool temporal = attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL;
    if (temporal) {
        if (!attribute_numeric_resident(attr)) {
            *reason = "Produced on demand, a frame at a time: showing it would produce every frame";
            return 0;
        }
        if (md_attributes_axis(table, attr) == attr) {
            *reason = "A frame axis: the time the others are sampled at";
            return 0;
        }

        uint32_t views = SeriesView_Distribution;
        const size_t n = attr->format.shape[0];
        if (is_script(key.source)) {
            if (n == md_array_size(app->timeline.x_values)) views |= SeriesView_Timeline;
        } else {
            // What series_resolve_temporal needs of the axis, short of converting it
            const md_attribute_t* axis = md_attributes_axis(table, attr);
            const md_attribute_t* run_axis = run_time_axis(app);
            if (axis && axis->format.shape[0] == n) {
                if ((run_axis && md_attribute_same_data(axis, run_axis)) ||
                    md_unit_is_none(app->timeline.time_unit) == md_unit_is_none(axis->unit)) {
                    views |= SeriesView_Timeline;
                }
            }
        }
        if (!(views & SeriesView_Timeline)) {
            *reason = "Its frame axis cannot be placed on the timeline (no trajectory, or time without a unit on one side)";
        }
        return views;
    }

    if (is_script(key.source) && (script_property_kind(app, key) & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) && attribute_numeric_resident(attr)) {
        return SeriesView_Distribution;
    }

    *reason = "Not sampled over time, and not a distribution";
    return 0;
}

// ## Over time

bool series_resolve_temporal(SeriesTemporalView* out, ApplicationState* app, const SeriesKey& key) {
    ASSERT(out && app);
    *out = SeriesTemporalView{};
    fill_names(out, app, key);

    if (key.variant >= SeriesVariant_Count || key.variant == SeriesVariant_Aggregate || key.path[0] == '\0') return false;

    const md_attributes_t* table = series_table(app, key.source);
    if (!table) return false;

    const str_t path = str_from_cstr(key.path);
    const md_attribute_t* attr = md_attributes_find(table, path);
    if (!attribute_temporal(attr)) return false;

    // What is drawn, and for a band around a centre, what spreads it
    const md_attribute_t* src  = attr;
    const md_attribute_t* src2 = nullptr;
    switch (key.variant) {
    case SeriesVariant_Mean:
        src = md_attributes_find_in(table, path, str_from_cstr(variant_leaf[SeriesVariant_Mean]));
        break;
    case SeriesVariant_Sigma:
        src  = md_attributes_find_in(table, path, str_from_cstr(variant_leaf[SeriesVariant_Mean]));
        src2 = md_attributes_find_in(table, path, str_from_cstr(variant_leaf[SeriesVariant_Sigma]));
        if (!attribute_temporal(src2)) return false;
        break;
    case SeriesVariant_Extent:
        src = md_attributes_find_in(table, path, str_from_cstr(variant_leaf[SeriesVariant_Extent]));
        break;
    default:
        break;
    }
    if (!attribute_temporal(src)) return false;

    const size_t n = src->format.shape[0];
    const size_t stride = md_attribute_element_count(&src->format) / n;
    if (src2 && md_attribute_element_count(&src2->format) != n) return false;
    if (key.variant == SeriesVariant_Extent && stride != 2) return false;
    if ((key.variant == SeriesVariant_Mean || key.variant == SeriesVariant_Sigma) && stride != 1) return false;

    // Shown in the display unit of the values themselves. The variance is in their square and is
    // read through its square root, so the same factor applies to the band.
    md_unit_t shown;
    out->y_scale = (float)display_units::factor(&shown, attr->unit);
    md_unit_print(out->unit_str, sizeof(out->unit_str), shown);

    // ### x
    const size_t num_frames = md_array_size(app->timeline.x_values);
    const float* x = nullptr;
    if (is_script(key.source)) {
        // Evaluated over the run's frames, so its axis is the timeline's
        if (n != num_frames) return false;
        x = app->timeline.x_values;
    } else {
        const md_attribute_t* axis = md_attributes_axis(table, src);
        if (!axis || axis->format.shape[0] != n) return false;

        const md_attribute_t* run_axis = run_time_axis(app);
        if (run_axis && md_attribute_same_data(axis, run_axis) && num_frames == n) {
            x = app->timeline.x_values;
        } else {
            // Ordinals place only against ordinals: a run without time has nothing to put a
            // picosecond on, and a file without time nothing to put on a run with it
            const md_unit_t time_unit = app->timeline.time_unit;
            if (md_unit_is_none(time_unit) != md_unit_is_none(axis->unit)) return false;
            size_t count = 0;
            x = cache_convert(app, table, axis, time_unit, &count);
            if (!x || count != n) return false;
        }
    }

    // ### y
    const float* y  = attribute_floats(app, table, src);
    const float* y2 = src2 ? attribute_floats(app, table, src2) : nullptr;
    if (!y || (src2 && !y2)) return false;

    out->x = x;
    out->y = y;
    out->y2 = y2;
    out->num_samples = (int)n;
    out->stride = (int)stride;

    switch (key.variant) {
    case SeriesVariant_Values:
        out->dim = (int)MIN(stride, (size_t)PLOT_MAX_POPULATION);
        break;
    case SeriesVariant_Mean:
        out->dim = 1;
        break;
    case SeriesVariant_Sigma:
    case SeriesVariant_Extent:
        out->dim  = 1;
        out->band = true;
        break;
    default:
        return false;
    }
    return true;
}

void series_temporal_band(const SeriesTemporalView& v, int i, double* lo, double* hi) {
    if (v.y2) {
        const double mean = v.y[i];
        const double sd   = sqrt(MAX(0.0, (double)v.y2[i]));
        *lo = (mean - sd) * v.y_scale;
        *hi = (mean + sd) * v.y_scale;
    } else {
        *lo = v.y[(size_t)i * v.stride + 0] * v.y_scale;
        *hi = v.y[(size_t)i * v.stride + 1] * v.y_scale;
    }
}

double series_temporal_index_at(const SeriesTemporalView& v, double x) {
    const int n = v.num_samples;
    if (n <= 0 || !v.x) return 0.0;
    if (x <= v.x[0])     return 0.0;
    if (x >= v.x[n - 1]) return (double)(n - 1);
    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (v.x[mid] <= x) lo = mid;
        else hi = mid;
    }
    const double dx = (double)v.x[hi] - (double)v.x[lo];
    return lo + (dx > 0.0 ? (x - v.x[lo]) / dx : 0.0);
}

int series_temporal_print_value(char* buf, size_t cap, const SeriesTemporalView& v, int i, int k) {
    if (i < 0 || i >= v.num_samples) return snprintf(buf, cap, "-");
    if (v.band) {
        double lo, hi;
        series_temporal_band(v, i, &lo, &hi);
        if (v.y2) {
            return snprintf(buf, cap, "%.4g +- %.4g", 0.5 * (lo + hi), 0.5 * (hi - lo));
        }
        return snprintf(buf, cap, "%.4g .. %.4g", lo, hi);
    }
    return snprintf(buf, cap, "%.4g", series_temporal_point(v, i, k).y);
}

// ## As a distribution

// Density of values[f * stride + k] over the frames f in mask (all num_frames without one), per
// member k or over all of them at once. Values outside [lo, hi] are not counted.
static void bin_temporal(float* bins, int num_bins, int dim_out, double lo, double hi, const float* values, size_t num_frames, size_t stride, const md_bitfield_t* mask) {
    const bool aggregate = dim_out == 1 && stride > 1;
    MEMSET(bins, 0, sizeof(float) * (size_t)num_bins * dim_out);

    md_temp_scope_t temp = md_temp_begin();
    size_t* count = md_temp_alloc_array(temp, size_t, dim_out);
    MEMSET(count, 0, sizeof(size_t) * dim_out);

    const double inv_range = hi > lo ? 1.0 / (hi - lo) : 0.0;
    auto add_frame = [&](size_t f) {
        const float* row = values + f * stride;
        const size_t members = aggregate ? stride : (size_t)dim_out;
        for (size_t k = 0; k < members; ++k) {
            const double val = row[k];
            if (!(lo <= val && val <= hi)) continue;
            const int b = CLAMP((int)((val - lo) * inv_range * num_bins), 0, num_bins - 1);
            const size_t d = aggregate ? 0 : k;
            bins[d * num_bins + b] += 1.0f;
            count[d] += 1;
        }
    };

    if (mask) {
        md_bitfield_iter_t it = md_bitfield_iter_create(mask);
        while (md_bitfield_iter_next(&it)) {
            const size_t f = md_bitfield_iter_idx(&it);
            if (f >= num_frames) break;
            add_frame(f);
        }
    } else {
        for (size_t f = 0; f < num_frames; ++f) {
            add_frame(f);
        }
    }

    const double width = (hi - lo) / num_bins;
    for (int d = 0; d < dim_out; ++d) {
        if (count[d] == 0 || width <= 0.0) continue;
        const float scl = (float)(1.0 / (width * count[d]));
        for (int b = 0; b < num_bins; ++b) {
            bins[d * num_bins + b] *= scl;
        }
    }
    md_temp_end(temp);
}

// A distribution evaluated at num_src bins, averaged down to num_dst (a divisor of it): each
// destination bin is the weighted mean of the source bins it covers
static void downsample_bins(float* dst, int num_dst, const float* src, const float* weights, int num_src) {
    const int factor = MAX(1, num_src / num_dst);
    for (int d = 0; d < num_dst; ++d) {
        double sum = 0.0;
        double weight = 0.0;
        for (int i = 0; i < factor; ++i) {
            const int s = d * factor + i;
            if (s >= num_src) break;
            sum    += src[s];
            weight += weights ? weights[s] : 1.0;
        }
        dst[d] = weight > 0.0 ? (float)(sum / weight) : 0.0f;
    }
}

// A (min, max) pair: rank 0 with two components, as the script publishes it, or two values of any shape
static bool read_range(double out[2], const md_attribute_t* range) {
    if (!range || !range->data || range->format.type == MD_ATTRIBUTE_TYPE_NONE || range->format.type == MD_ATTRIBUTE_TYPE_STR) return false;
    if (md_attribute_element_count(&range->format) != 2) return false;
    return md_attribute_extract_f64(out, 2, range, md_attribute_slice_all(), md_unit_none()) == 2 && out[1] > out[0];
}

bool series_resolve_histogram(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins) {
    ASSERT(out && app);
    *out = SeriesHistogramView{};
    fill_names(out, app, key);

    if (key.path[0] == '\0') return false;
    if (key.variant != SeriesVariant_Values && key.variant != SeriesVariant_Aggregate) return false;

    const md_attributes_t* table = series_table(app, key.source);
    if (!table) return false;

    const str_t path = str_from_cstr(key.path);
    const md_attribute_t* attr = md_attributes_find(table, path);
    if (!attribute_numeric_resident(attr)) return false;

    const bool temporal     = attribute_temporal(attr);
    const bool distribution = !temporal && (script_property_kind(app, key) & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION);
    if (!temporal && !distribution) return false;

    const md_attribute_t* weight = distribution ? md_attributes_find_in(table, path, STR_LIT("weight")) : nullptr;
    if (weight && (!attribute_numeric_resident(weight) || md_attribute_element_count(&weight->format) != md_attribute_element_count(&attr->format))) {
        weight = nullptr;
    }

    const size_t n = attr->format.shape[0];
    const size_t stride = md_attribute_element_count(&attr->format) / n;

    int dim = 1;
    if (temporal) {
        dim = (key.variant == SeriesVariant_Aggregate) ? 1 : (int)MIN(stride, (size_t)PLOT_MAX_POPULATION);
        num_bins = CLAMP(num_bins, 2, 4096);
    } else {
        // Evaluated at a resolution of its own; only coarser is on offer
        if (stride != 1 || key.variant != SeriesVariant_Values) return false;
        num_bins = CLAMP(num_bins, 1, (int)n);
        while (n % (size_t)num_bins) --num_bins;
    }

    // ### Cached bins
    uint64_t hkey = md_hash64_str(path, (uint64_t)key.source * 31 + key.variant);
    hkey = md_hash64(&num_bins, sizeof(num_bins), hkey);
    const uint64_t version = attr->version + (weight ? weight->version << 1 : 0);

    const int frame = ImGui::GetFrameCount();
    HistogramEntry* entry = nullptr;
    for (size_t i = 0; i < md_array_size(histogram_cache); ++i) {
        if (histogram_cache[i].key == hkey && histogram_cache[i].table == table) {
            entry = &histogram_cache[i];
            break;
        }
    }
    if (!entry) {
        md_array_push(histogram_cache, HistogramEntry{}, app->allocator.persistent);
        entry = &histogram_cache[md_array_size(histogram_cache) - 1];
        entry->key = hkey;
        entry->table = table;
        entry->version = ~version;
    }
    entry->last_used = frame;

    if (entry->version != version) {
        entry->version = version;
        entry->valid = false;
        entry->num_bins = num_bins;
        entry->dim = dim;
        md_array_resize(entry->bins, (size_t)num_bins * dim, app->allocator.persistent);

        const float* values = attribute_floats(app, table, attr);
        if (values) {
            double range[2] = {0, 0};
            const bool has_range = read_range(range, md_attributes_find_in(table, path, STR_LIT("range")));
            if (temporal) {
                const md_bitfield_t* mask = nullptr;
                if (is_script(key.source)) {
                    const md_script_eval_t* eval = series_eval(app, key.source);
                    mask = eval ? md_script_eval_frame_mask(eval) : nullptr;
                }
                if (!has_range) {
                    // The range the values cover, over the frames that count
                    double lo = DBL_MAX, hi = -DBL_MAX;
                    for (size_t i = 0; i < n * stride; ++i) {
                        const double v = values[i];
                        if (v == v) { lo = MIN(lo, v); hi = MAX(hi, v); }
                    }
                    if (lo > hi) { lo = 0; hi = 1; }
                    if (hi <= lo) { const double pad = MAX(fabs(lo) * 1e-3, 1e-6); lo -= pad; hi += pad; }
                    range[0] = lo;
                    range[1] = hi;
                }
                bin_temporal(entry->bins, num_bins, dim, range[0], range[1], values, n, stride, mask);
                entry->valid = true;
            } else if (has_range) {
                const float* w = weight ? attribute_floats(app, table, weight) : nullptr;
                downsample_bins(entry->bins, num_bins, values, w, (int)n);
                entry->valid = true;
            }
            entry->x_min = range[0];
            entry->x_max = range[1];

            float y_min = FLT_MAX, y_max = -FLT_MAX;
            for (size_t i = 0; i < (size_t)num_bins * dim; ++i) {
                y_min = MIN(y_min, entry->bins[i]);
                y_max = MAX(y_max, entry->bins[i]);
            }
            entry->y_min = entry->valid ? y_min : 0.0f;
            entry->y_max = entry->valid ? y_max : 0.0f;
        }
    }
    if (!entry->valid) return false;

    // ### Units
    // A temporal attribute binned: its values along x, and along y a probability density, which is
    // per unit of x and so rescales against it - it integrates to one in the unit it is shown in.
    // A script distribution: its bin axis along x, and its own values along y.
    double x_scl = 1.0;
    double y_scl = 1.0;
    if (distribution) {
        const md_attribute_t* bin = md_attributes_find_in(table, path, STR_LIT("bin"));
        x_scl = display_units::factor_print(out->x_unit_str, sizeof(out->x_unit_str), bin ? bin->unit : md_unit_none());
        y_scl = display_units::factor_print(out->y_unit_str, sizeof(out->y_unit_str), attr->unit);
    } else {
        x_scl = display_units::factor_print(out->x_unit_str, sizeof(out->x_unit_str), attr->unit);
        y_scl = x_scl != 0.0 ? 1.0 / x_scl : 1.0;
        if (out->x_unit_str[0] != '\0') {
            // "1/nm", but "1/(kJ/mol)": a compound unit is inverted whole
            const bool compound = strpbrk(out->x_unit_str, "/ *^") != nullptr;
            snprintf(out->y_unit_str, sizeof(out->y_unit_str), compound ? "1/(%s)" : "1/%s", out->x_unit_str);
        }
    }

    out->bins = entry->bins;
    out->num_bins = entry->num_bins;
    out->dim = entry->dim;
    out->x_min = entry->x_min * x_scl;
    out->x_max = entry->x_max * x_scl;
    out->y_scale = y_scl;
    out->y_min = entry->y_min * y_scl;
    out->y_max = entry->y_max * y_scl;
    return true;
}

// ## Volumes

bool series_resolve_volume(SeriesVolumeView* out, ApplicationState* app, const SeriesKey& key) {
    ASSERT(out && app);
    *out = SeriesVolumeView{};
    series_label(out->label, sizeof(out->label), app, key);
    if (!is_script(key.source) || key.variant != SeriesVariant_Values ||
        !(series_views(app, key) & SeriesView_Volume)) return false;
    const md_script_eval_t* eval = series_eval(app, key.source);
    if (!eval || md_script_eval_ir_fingerprint(eval) != md_script_ir_fingerprint(app->script.eval_ir)) return false;
    if (key.source == SeriesSource_ScriptFiltered && !app->timeline.filter.enabled) return false;
    const str_t ident = series_script_ident(key);
    str_copy_to_char_buf(out->script_ident, sizeof(out->script_ident), ident);
    out->vis_payload = app->script.eval_ir ? md_script_ir_property_vis_payload(app->script.eval_ir, ident) : nullptr;

    const md_attributes_t* table = series_table(app, key.source);
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr || !attr->data || !out->vis_payload || !md_script_ir_valid(app->script.eval_ir)) return false;
    out->attr = attr;
    out->version = attr->version;
    return true;
}

bool series_resolve_vector(SeriesVectorView* out, ApplicationState* app, const SeriesKey& key, uint32_t components) {
    ASSERT(out && app);
    *out = SeriesVectorView{};
    series_label(out->label, sizeof(out->label), app, key);
    if (key.source != SeriesSource_Script || key.variant != SeriesVariant_Values || !app->script.eval_ir ||
        !md_script_ir_valid(app->script.eval_ir)) return false;
    const md_script_eval_t* eval = series_eval(app, key.source);
    if (!eval || md_script_eval_ir_fingerprint(eval) != md_script_ir_fingerprint(app->script.eval_ir)) return false;
    const str_t ident = series_script_ident(key);
    if (!(md_script_ir_property_flags(app->script.eval_ir, ident) & MD_SCRIPT_PROPERTY_FLAG_VECTOR)) return false;
    str_copy_to_char_buf(out->script_ident, sizeof(out->script_ident), ident);
    out->vis_payload = md_script_ir_property_vis_payload(app->script.eval_ir, ident);

    const md_attributes_t* table = series_table(app, key.source);
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr || !attr->data || attr->format.type != MD_ATTRIBUTE_TYPE_F32 || attr->format.rank != 2 ||
        attr->format.components != components || md_attribute_is_virtual(attr)) return false;
    out->data       = (const float*)attr->data;
    out->num_frames = attr->format.shape[0];
    out->count      = attr->format.shape[1];
    out->components = components;
    out->version    = attr->version;
    return true;
}

void plot_series_style_popup(ApplicationState* data, PlotSeries& s, int dim, const char* script_ident, const md_script_vis_payload_o* vis_payload, int* hovered_pop_idx) {
    if (dim > 1) {
        const char* color_type_labels[] = {"Solid", "Colormap"};
        int color_type = s.use_colormap ? 1 : 0;
        if (ImGui::Combo("Color Type", &color_type, color_type_labels, IM_ARRAYSIZE(color_type_labels))) {
            s.use_colormap = color_type == 1;
        }
    }
    if (s.use_colormap && dim > 1) {
        ImPlot::ColormapSelection("##Colormap", &s.colormap);
        ImGui::SliderFloat("Alpha", &s.colormap_alpha, 0.0f, 1.0f);
    } else {
        ImGui::ColorEdit4("Color", &s.color.x);
    }

    if (dim > 1) {
        ImGui::Separator();
        if (ImGui::Button("Set All")) {
            s.population_mask.set();
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear All")) {
            s.population_mask.reset();
        }

        const float sz = ImGui::GetFontSize() * 1.5f;
        ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
        for (int k = 0; k < MIN(dim, PLOT_MAX_POPULATION); ++k) {
            char lbl[32];
            snprintf(lbl, sizeof(lbl), "%d", k+1);
            if (ImGui::Selectable(lbl, s.population_mask.test(k), ImGuiSelectableFlags_DontClosePopups, ImVec2(sz, sz))) {
                s.population_mask.flip(k);
            }
            if (ImGui::IsItemHovered()) {
                if (script_ident[0] != '\0') {
                    script_visualize_payload(data, vis_payload, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                    script_set_hovered_property(data, str_from_cstr(script_ident), k);
                }
                *hovered_pop_idx = k;
            }
            if (!k || ((k+1) % 10)) {
                ImGui::SameLine();
            }
        }
        ImGui::PopStyleVar();
    }
}

// One entry in a list of series: a click toggles it in the first subplot, a drag puts it in any
void plot_series_list_item(ApplicationState* data, const char* dnd_type, PlotSubplot& first, const SeriesKey& key) {
    char label[96];
    series_label(label, sizeof(label), data, key);
    const ImVec4 color = series_default_color(data, key);
    const int idx = plot_find_series(first, key);

    ImGui::PushID(key.path);
    ImGui::PushID((int)key.variant * SeriesSource_Count + (int)key.source);
    ImPlot::ItemIcon(color);
    ImGui::SameLine();
    if (ImGui::Selectable(label, idx != -1, ImGuiSelectableFlags_DontClosePopups)) {
        if (idx != -1) {
            plot_remove_series(first, idx);
        } else {
            plot_add_series(data, first, key);
        }
    }
    if (ImGui::IsItemHovered()) {
        const str_t ident = series_script_ident(key);
        if (!str_empty(ident)) {
            const md_script_vis_payload_o* vis = data->script.eval_ir ? md_script_ir_property_vis_payload(data->script.eval_ir, ident) : nullptr;
            script_visualize_payload(data, vis, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(data, ident);
        } else {
            ImGui::SetTooltip("%s\nClick to show in the first subplot, drag into any", key.path);
        }
    }
    if (ImGui::BeginDragDropSource()) {
        series_set_drag_payload(dnd_type, key, -1, label, color);
        ImGui::EndDragDropSource();
    }
    ImGui::PopID();
    ImGui::PopID();
}

void plot_vector_style_popup(PlotSeries& s, bool point_color) {
    ImGui::SeparatorText("Points");
    ImGui::Checkbox("Displayed frame", &s.show_points);
    if (s.show_points) {
        if (ImGui::BeginCombo("Marker", ImPlot::GetMarkerName(s.marker))) {
            for (int m = 0; m < ImPlotMarker_COUNT; ++m) {
                if (ImGui::Selectable(ImPlot::GetMarkerName(m), s.marker == m)) s.marker = (ImPlotMarker)m;
            }
            ImGui::EndCombo();
        }
        ImGui::SliderFloat("Marker Size", &s.marker_size, 0.1f, 10.0f, "%.2f");
    }
    if (point_color) ImGui::ColorEdit4("Color", &s.color.x);

    ImGui::SeparatorText("Density");
    ImGui::Checkbox("Isolines", &s.iso_lines);
    ImGui::Checkbox("Filled levels", &s.iso_levels);
    ImGui::Checkbox("Colormap", &s.density_map);
    if (s.iso_lines || s.iso_levels || s.density_map) {
        if (s.density_map) ImPlot::ColormapSelection("Density Color Map", &s.density_colormap);
        ImGui::SliderFloat("Opacity", &s.density_alpha, 0.0f, 1.0f);
        ImGui::SliderFloat("Smoothing", &s.smoothing, 0.25f, 6.0f, "%.2f");
        const char* range_names[VectorRange_Count] = {"Full trajectory", "Filtered trajectory"};
        int range = (int)s.range;
        if (ImGui::Combo("Frames", &range, range_names, VectorRange_Count)) s.range = (VectorRange)range;
    }
}

void series_vector_frame_range(const ApplicationState* app, const PlotSeries& s, uint32_t num_frames, uint32_t* beg, uint32_t* end) {
    *beg = 0;
    *end = num_frames;
    if (num_frames == 0) return;
    if (s.range == VectorRange_Filtered && app->timeline.filter.enabled) {
        *beg = MIN((uint32_t)MAX(0.0, app->timeline.filter.beg_frame), num_frames - 1);
        *end = MIN((uint32_t)MAX(0.0, app->timeline.filter.end_frame) + 1, num_frames);
        if (*end < *beg) *end = *beg;
    }
}

int plot_add_vector_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key, bool allow_duplicate) {
    if (!allow_duplicate) {
        const int existing = plot_find_series(sp, key);
        if (existing != -1) return existing;
    }
    const int idx = plot_add_series_copy(app, sp, key);
    if (idx != -1) {
        sp.series[idx].marker = ImPlotMarker_Circle;
        sp.series[idx].marker_size = 5.0f;
    }
    return idx;
}

void series_vector_properties_menu(ApplicationState* app, PlotSubplot& sp, uint32_t components) {
    int num_listed = 0;
    series_for_each_script_property(app, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_VECTOR, [&](const SeriesKey& key) {
        SeriesVectorView v;
        if (!series_resolve_vector(&v, app, key, components)) return;
        num_listed += 1;

        char label[96];
        series_label(label, sizeof(label), app, key);
        const ImVec4 color = series_default_color(app, key);
        const int idx = plot_find_series(sp, key);

        ImGui::PushID(key.path);
        ImPlot::ItemIcon(color);
        ImGui::SameLine();
        if (ImGui::Selectable(label, idx != -1, ImGuiSelectableFlags_DontClosePopups)) {
            if (idx != -1) plot_remove_series(sp, idx);
            else           plot_add_vector_series(app, sp, key, false);
        }
        if (ImGui::IsItemHovered()) {
            script_visualize_payload(app, v.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(app, str_from_cstr(v.script_ident));
        }
        if (ImGui::BeginDragDropSource()) {
            series_set_drag_payload(VECTOR_SERIES_DND, key, -1, label, color);
            ImGui::EndDragDropSource();
        }
        ImGui::PopID();
    });
    if (num_listed == 0) {
        ImGui::TextDisabled("No properties with %u components, create one below or in the script editor", components);
    }
}

void series_vector_layers_menu(ApplicationState* app, PlotSubplot& sp, uint32_t components) {
    (void)components;
    if (sp.count == 0) {
        ImGui::TextDisabled("No series");
        return;
    }
    for (int j = 0; j < sp.count; ++j) {
        PlotSeries& s = sp.series[j];
        char label[96];
        series_label(label, sizeof(label), app, s.key);
        ImGui::PushID((int)s.uid);
        ImPlot::ItemIcon(s.color);
        ImGui::SameLine();
        if (ImGui::TreeNode(label)) {
            plot_vector_style_popup(s);
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
}

void series_vector_bar(ApplicationState* app, PlotSubplot& sp, uint32_t components, bool members) {
    // A property dropped anywhere on the window is added. One dragged out of a plot is already in a window.
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (ImGui::BeginDragDropTargetCustom(window->InnerRect, ImGui::GetID("##vector-drop"))) {
        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(VECTOR_SERIES_DND)) {
            if (payload->DataSize == sizeof(SeriesDragPayload)) {
                SeriesVectorView view;
                const SeriesDragPayload* dnd = (const SeriesDragPayload*)payload->Data;
                if (dnd->src_subplot < 0 && series_resolve_vector(&view, app, dnd->key, components)) {
                    plot_add_vector_series(app, sp, dnd->key, true);
                }
            }
        }
        ImGui::EndDragDropTarget();
    }

    int remove_idx = -1;
    int duplicate_idx = -1;
    const float spacing = ImGui::GetStyle().ItemSpacing.x;
    const float icon = ImGui::GetFontSize();

    for (int j = 0; j < sp.count; ++j) {
        PlotSeries& s = sp.series[j];
        SeriesVectorView v;
        const bool resolved = series_resolve_vector(&v, app, s.key, components);
        char label[96];
        series_label(label, sizeof(label), app, s.key);

        const float text_w = ImGui::CalcTextSize(label).x;
        if (j > 0) {
            const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
            if (ImGui::GetItemRectMax().x + spacing + icon + spacing + text_w < right) ImGui::SameLine();
        }

        ImGui::PushID((int)s.uid);
        ImGui::ColorButton("##color", s.color, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoDragDrop, ImVec2(icon, icon));
        ImGui::SameLine();
        ImGui::Selectable(label, false, ImGuiSelectableFlags_None, ImVec2(text_w, 0));
        if (ImGui::IsItemHovered() && resolved) {
            script_visualize_payload(app, v.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(app, str_from_cstr(v.script_ident));
        }
        if (ImGui::BeginDragDropSource()) {
            series_set_drag_payload(VECTOR_SERIES_DND, s.key, 0, label, s.color, s.uid);
            ImGui::EndDragDropSource();
        }
        if (ImGui::BeginPopupContextItem("##series-style")) {
            if (ImGui::DeleteButton("Remove")) {
                remove_idx = j;
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Duplicate")) {
                duplicate_idx = j;
                ImGui::CloseCurrentPopup();
            }
            // With several members the colour is one per member, from a colormap or the series'
            const bool population = members && resolved && v.count > 1;
            plot_vector_style_popup(s, !population);
            if (population) {
                ImGui::SeparatorText("Members");
                int hovered_pop_idx = -1;
                plot_series_style_popup(app, s, (int)v.count, v.script_ident, v.vis_payload, &hovered_pop_idx);
            }
            ImGui::EndPopup();
        }
        ImGui::PopID();
    }

    if (duplicate_idx != -1 && duplicate_idx < sp.count) {
        const PlotSeries original = sp.series[duplicate_idx];
        const int idx = plot_add_series_copy(app, sp, original.key);
        if (idx != -1) {
            const uint32_t uid = sp.series[idx].uid;
            sp.series[idx] = original;
            sp.series[idx].uid = uid;
        }
    }
    if (remove_idx != -1) plot_remove_series(sp, remove_idx);
}

void series_vector_drag_out(PlotSubplot& sp) {
    if (!ImGui::IsMouseReleased(ImGuiMouseButton_Left)) return;
    const ImGuiPayload* payload = ImGui::GetDragDropPayload();
    if (payload && payload->IsDataType(VECTOR_SERIES_DND) && !ImGui::IsDragDropPayloadBeingAccepted()) {
        const SeriesDragPayload* dnd = (const SeriesDragPayload*)payload->Data;
        const int idx = dnd && dnd->src_subplot == 0 ? plot_find_uid(sp, dnd->src_uid) : -1;
        if (idx != -1) plot_remove_series(sp, idx);
    }
}

// ## Script properties

size_t series_script_property_count(const ApplicationState* app) {
    return app->script.eval_ir ? md_script_ir_property_count(app->script.eval_ir) : 0;
}

bool series_script_property(SeriesKey* out, const ApplicationState* app, SeriesSource source, size_t index, uint32_t kind_flags) {
    const md_script_ir_t* ir = app->script.eval_ir;
    if (!ir || !is_script(source) || index >= md_script_ir_property_count(ir)) return false;
    const str_t name = md_script_ir_property_names(ir)[index];
    if (!(md_script_ir_property_flags(ir, name) & kind_flags)) return false;

    const md_attributes_t* table = series_table(app, source);
    if (!table) return false;

    char path[SERIES_PATH_CAP];
    const int len = snprintf(path, sizeof(path), STR_FMT STR_FMT, STR_ARG(script_prefix), STR_ARG(name));
    if (len <= 0 || (size_t)len >= sizeof(path)) return false;
    if (!md_attributes_find(table, str_t{path, (size_t)len})) return false;

    *out = series_key(source, str_t{path, (size_t)len});
    return true;
}

// ## Subplots

int plot_find_series(const PlotSubplot& sp, const SeriesKey& key) {
    for (int i = 0; i < sp.count; ++i) {
        if (series_key_equal(sp.series[i].key, key)) return i;
    }
    return -1;
}

static uint32_t next_series_uid() {
    static uint32_t counter = 0;
    return ++counter;
}

int plot_find_uid(const PlotSubplot& sp, uint32_t uid) {
    for (int i = 0; uid != 0 && i < sp.count; ++i) {
        if (sp.series[i].uid == uid) return i;
    }
    return -1;
}

int plot_add_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key) {
    const int existing = plot_find_series(sp, key);
    if (existing != -1) return existing;
    return plot_add_series_copy(app, sp, key);
}

int plot_add_series_copy(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key) {
    if (key.path[0] == '\0') return -1;
    if (sp.count >= PLOT_MAX_SERIES_PER_SUBPLOT) {
        MD_LOG_INFO("A subplot holds at most %d series", PLOT_MAX_SERIES_PER_SUBPLOT);
        return -1;
    }
    PlotSeries& s = sp.series[sp.count];
    s = PlotSeries{};
    s.key = key;
    s.uid = next_series_uid();
    s.color = series_default_color(app, key);
    s.population_mask.set();
    return sp.count++;
}

void plot_remove_series(PlotSubplot& sp, int index) {
    if (index < 0 || index >= sp.count) return;
    for (int i = index; i < sp.count - 1; ++i) {
        sp.series[i] = sp.series[i + 1];
    }
    sp.count -= 1;
}

void plot_move_series(const ApplicationState* app, PlotSubplot& src, PlotSubplot& dst, const SeriesKey& key) {
    if (&src == &dst) return;
    const int src_idx = plot_find_series(src, key);
    if (src_idx == -1) {
        plot_add_series(app, dst, key);
        return;
    }
    if (plot_find_series(dst, key) == -1) {
        if (dst.count >= PLOT_MAX_SERIES_PER_SUBPLOT) return;
        dst.series[dst.count++] = src.series[src_idx];
    }
    plot_remove_series(src, src_idx);
}

void plot_clear(PlotSubplot* subplots, int count) {
    for (int i = 0; i < count; ++i) {
        subplots[i].count = 0;
    }
}

// ## Series loaded along the run

size_t system_series_groups(str_t out_groups[], size_t cap, const ApplicationState* app) {
    ASSERT(app);
    if (app->mold.run[0] == '\0') return 0;
    const md_attributes_t* table = &app->mold.sys.attributes;
    const str_t run = str_from_cstr(app->mold.run);

    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter(table, run); md_attributes_next(&it);) {
        const md_attribute_t* attr = it.attr;
        if (attr->format.type != MD_ATTRIBUTE_TYPE_STR) continue;
        if (!str_eq_cstr(md_attribute_leaf(attr), "source")) continue;
        const str_t group = md_attribute_group(attr);
        if (str_empty(group) || str_eq(group, run)) continue;
        if (out_groups && count < cap) out_groups[count] = group;
        count += 1;
    }
    return count;
}

str_t system_series_group_source(const ApplicationState* app, str_t group) {
    const md_attributes_t* table = &app->mold.sys.attributes;
    const md_attribute_t* src = md_attributes_find_in(table, group, STR_LIT("source"));
    if (!src || src->format.type != MD_ATTRIBUTE_TYPE_STR) return {};
    return md_attribute_str(table, src, 0);
}

void system_series_group_label(char* buf, size_t cap, const ApplicationState* app, str_t group) {
    ASSERT(buf && cap > 0);
    str_t rel = group;
    const str_t run = str_from_cstr(app->mold.run);
    if (!str_empty(run) && str_begins_with(group, run) && group.len > run.len + 1) {
        rel = str_substr(group, run.len + 1);
    }
    str_t file = system_series_group_source(app, group);
    size_t loc;
    if (str_rfind_char(&loc, file, '/') || str_rfind_char(&loc, file, '\\')) {
        file = str_substr(file, loc + 1);
    }
    if (str_empty(file)) {
        snprintf(buf, cap, STR_FMT, STR_ARG(rel));
    } else {
        snprintf(buf, cap, STR_FMT "  (" STR_FMT ")", STR_ARG(rel), STR_ARG(file));
    }
}

size_t system_series_members(str_t out_paths[], size_t cap, const ApplicationState* app, str_t group) {
    ASSERT(app);
    const md_attributes_t* table = &app->mold.sys.attributes;

    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter(table, group); md_attributes_next(&it);) {
        const md_attribute_t* attr = it.attr;
        if (!attribute_temporal(attr)) continue;
        if (md_attributes_axis(table, attr) == attr) continue;  // The group's own time
        if (out_paths && count < cap) out_paths[count] = attr->path;
        count += 1;
    }
    return count;
}

// ## Workspace

void plot_layout_serialize(viamd::serialization_state_t& state, str_t section, str_t series_section, const PlotSubplot* subplots, int num_subplots) {
    viamd::write_section_header(state, section);
    viamd::write_int(state, STR_LIT("NumSubplots"), num_subplots);

    for (int i = 0; i < PLOT_MAX_SUBPLOTS; ++i) {
        const PlotSubplot& sp = subplots[i];
        for (int j = 0; j < sp.count; ++j) {
            const PlotSeries& s = sp.series[j];
            viamd::write_section_header(state, series_section);
            viamd::write_int (state, STR_LIT("Subplot"), i);
            viamd::write_str (state, STR_LIT("Source"),  str_from_cstr(source_name[s.key.source < SeriesSource_Count ? s.key.source : 0]));
            viamd::write_str (state, STR_LIT("Variant"), str_from_cstr(variant_name[s.key.variant < SeriesVariant_Count ? s.key.variant : 0]));
            viamd::write_str (state, STR_LIT("Path"),    str_from_cstr(s.key.path));
            viamd::write_vec4(state, STR_LIT("Color"),   vec4_set(s.color.x, s.color.y, s.color.z, s.color.w));
            viamd::write_str (state, STR_LIT("PlotType"), str_from_cstr(plot_type_name[s.plot_type < PlotType_Count ? s.plot_type : 0]));
            viamd::write_int (state, STR_LIT("Marker"),  (int)s.marker);
            viamd::write_flt (state, STR_LIT("MarkerSize"), s.marker_size);
            viamd::write_flt (state, STR_LIT("BarWidth"), s.bar_width);
            viamd::write_int (state, STR_LIT("NumBins"), s.num_bins);
            viamd::write_bool(state, STR_LIT("UseColormap"), s.use_colormap);
            viamd::write_int (state, STR_LIT("Colormap"), (int)s.colormap);
            viamd::write_flt (state, STR_LIT("ColormapAlpha"), s.colormap_alpha);
            viamd::write_bool(state, STR_LIT("ShowPoints"), s.show_points);
            viamd::write_bool(state, STR_LIT("IsoLines"), s.iso_lines);
            viamd::write_bool(state, STR_LIT("IsoLevels"), s.iso_levels);
            viamd::write_bool(state, STR_LIT("DensityMap"), s.density_map);
            viamd::write_int (state, STR_LIT("Range"), (int)s.range);
            viamd::write_int (state, STR_LIT("DensityColormap"), (int)s.density_colormap);
            viamd::write_flt (state, STR_LIT("DensityAlpha"), s.density_alpha);
            viamd::write_flt (state, STR_LIT("Smoothing"), s.smoothing);
        }
    }
}

bool plot_layout_deserialize(viamd::deserialization_state_t& state, str_t section, str_t series_section, const ApplicationState* app, PlotSubplot* subplots, int* num_subplots) {
    const str_t cur = viamd::section_header(state);
    str_t ident, arg;

    if (str_eq(cur, section)) {
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "NumSubplots")) {
                int n;
                if (viamd::extract_int(n, arg)) {
                    *num_subplots = CLAMP(n, 1, PLOT_MAX_SUBPLOTS);
                }
            }
        }
        return true;
    }

    if (str_eq(cur, series_section)) {
        PlotSeries s = {};
        s.population_mask.set();
        int  subplot = -1;
        bool has_color = false;
        while (viamd::next_entry(ident, arg, state)) {
            str_t str;
            int   i;
            if (str_eq_cstr(ident, "Subplot")) {
                viamd::extract_int(subplot, arg);
            } else if (str_eq_cstr(ident, "Source") && viamd::extract_str(str, arg)) {
                for (int k = 0; k < SeriesSource_Count; ++k) {
                    if (str_eq_cstr(str, source_name[k])) s.key.source = (SeriesSource)k;
                }
            } else if (str_eq_cstr(ident, "Variant") && viamd::extract_str(str, arg)) {
                for (int k = 0; k < SeriesVariant_Count; ++k) {
                    if (str_eq_cstr(str, variant_name[k])) s.key.variant = (SeriesVariant)k;
                }
            } else if (str_eq_cstr(ident, "Path")) {
                viamd::extract_to_char_buf(s.key.path, sizeof(s.key.path), arg);
            } else if (str_eq_cstr(ident, "Color")) {
                vec4_t c = {};
                viamd::extract_vec4(c, arg);
                s.color = ImVec4(c.x, c.y, c.z, c.w);
                has_color = true;
            } else if (str_eq_cstr(ident, "PlotType") && viamd::extract_str(str, arg)) {
                for (int k = 0; k < PlotType_Count; ++k) {
                    if (str_eq_cstr(str, plot_type_name[k])) s.plot_type = (PlotType)k;
                }
            } else if (str_eq_cstr(ident, "Marker") && viamd::extract_int(i, arg)) {
                s.marker = (ImPlotMarker)CLAMP(i, 0, ImPlotMarker_COUNT - 1);
            } else if (str_eq_cstr(ident, "MarkerSize")) {
                viamd::extract_flt(s.marker_size, arg);
            } else if (str_eq_cstr(ident, "BarWidth")) {
                viamd::extract_flt(s.bar_width, arg);
            } else if (str_eq_cstr(ident, "NumBins") && viamd::extract_int(i, arg)) {
                s.num_bins = CLAMP(i, 2, 4096);
            } else if (str_eq_cstr(ident, "UseColormap")) {
                viamd::extract_bool(s.use_colormap, arg);
            } else if (str_eq_cstr(ident, "Colormap") && viamd::extract_int(i, arg)) {
                s.colormap = (ImPlotColormap)i;
            } else if (str_eq_cstr(ident, "ColormapAlpha")) {
                viamd::extract_flt(s.colormap_alpha, arg);
            } else if (str_eq_cstr(ident, "ShowPoints")) {
                viamd::extract_bool(s.show_points, arg);
            } else if (str_eq_cstr(ident, "IsoLines")) {
                viamd::extract_bool(s.iso_lines, arg);
            } else if (str_eq_cstr(ident, "IsoLevels")) {
                viamd::extract_bool(s.iso_levels, arg);
            } else if (str_eq_cstr(ident, "DensityMap")) {
                viamd::extract_bool(s.density_map, arg);
            } else if (str_eq_cstr(ident, "Range") && viamd::extract_int(i, arg)) {
                s.range = i == 1 ? VectorRange_Filtered : VectorRange_Full;
            } else if (str_eq_cstr(ident, "DensityColormap") && viamd::extract_int(i, arg)) {
                s.density_colormap = (ImPlotColormap)i;
            } else if (str_eq_cstr(ident, "DensityAlpha")) {
                viamd::extract_flt(s.density_alpha, arg);
            } else if (str_eq_cstr(ident, "Smoothing")) {
                viamd::extract_flt(s.smoothing, arg);
            }
        }

        if (0 <= subplot && subplot < PLOT_MAX_SUBPLOTS) {
            const int idx = plot_add_series_copy(app, subplots[subplot], s.key);
            if (idx != -1) {
                PlotSeries& dst = subplots[subplot].series[idx];
                const ImVec4 default_color = dst.color;
                const uint32_t uid = dst.uid;
                dst = s;
                dst.uid = uid;
                if (!has_color) dst.color = default_color;
            }
        }
        return true;
    }

    return false;
}
