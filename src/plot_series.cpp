#include <plot_series.h>

#include <viamd.h>
#include <display_units.h>
#include <serialization_utils.h>

#include <md_system.h>
#include <md_script.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_str.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// ## Keys and names

bool series_key_equal(const SeriesKey& a, const SeriesKey& b) {
    return a.source == b.source && a.variant == b.variant && a.period == b.period && strcmp(a.path, b.path) == 0;
}

SeriesKey series_key(SeriesSource source, str_t path, SeriesVariant variant, uint32_t period) {
    SeriesKey key;
    key.source  = source;
    key.variant = variant;
    key.period  = period;
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
};

const char* series_source_name(SeriesSource source) {
    return source_name[source < SeriesSource_Count ? source : 0];
}

bool series_source_from_name(SeriesSource* out, str_t name, uint32_t* period) {
    for (int k = 0; k < SeriesSource_Count; ++k) {
        if (str_eq_cstr(name, source_name[k])) {
            *out = (SeriesSource)k;
            return true;
        }
    }
    if (str_eq_cstr(name, "script_filtered")) {
        *out = SeriesSource_Script;
        if (period) *period = SERIES_PERIOD_LEGACY_FILTER;
        return true;
    }
    return false;
}

uint32_t series_legacy_period(const ApplicationState* app, uint32_t period) {
    if (period != SERIES_PERIOD_LEGACY_FILTER) return period;
    const uint32_t id = app->timeline.periods.legacy_filter_id;
    return timeline_periods_find(&app->timeline.periods, id) ? id : 0;
}

static const char* plot_type_name[PlotType_Count] = {
    "line",
    "scatter",
    "area",
    "bars",
};

static const str_t script_prefix = STR_INIT("script/");

static bool is_script(SeriesSource source) {
    return source == SeriesSource_Script;
}

// The script's is the evaluation's table. Over a period, a temporal property is taken over the period's frames, and a
// distribution or a volume from the slot over them (series_slot).
const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source) {
    ASSERT(app);
    if (source == SeriesSource_System) {
        return &app->mold.sys.attributes;
    }
    const md_script_eval_t* eval = is_script(source) ? app->script.eval : nullptr;
    return eval ? md_script_eval_attributes(eval) : nullptr;
}

const md_attributes_t* series_table(const ApplicationState* app, SeriesSource source, str_t path) {
    (void)path;
    return series_table(app, source);
}

// The frames of the key's period, [beg, end[: every frame without one. False when it covers no frames.
static bool key_period_frames(uint32_t* beg, uint32_t* end, const ApplicationState* app, const SeriesKey& key) {
    return series_period_frames(beg, end, app, key, 0);
}

bool series_frame_mask(md_bitset_t* out, const ApplicationState* app, const SeriesKey& key, md_allocator_i* alloc) {
    ASSERT(out);
    const md_script_eval_t* eval = is_script(key.source) ? app->script.eval : nullptr;
    if (!eval) return false;
    md_bitset_t held = md_script_eval_property_frames(eval, series_script_ident(key), alloc);
    if (key.period != 0) {
        uint32_t beg, end;
        held = md_bitset_and(held, key_period_frames(&beg, &end, app, key) ? md_bitset_range(beg, end) : md_bitset_empty_set(), alloc);
    }
    *out = held;
    return true;
}

const char* series_period_name(const ApplicationState* app, uint32_t period) {
    const TimelinePeriod* r = timeline_periods_find(&app->timeline.periods, period);
    return r ? r->name : "missing period";
}

str_t series_script_ident(const SeriesKey& key) {
    str_t path = str_from_cstr(key.path);
    if (is_script(key.source) && str_begins_with(path, script_prefix)) {
        return str_substr(path, script_prefix.len);
    }
    return {};
}

SeriesState series_state(const ApplicationState* app, const SeriesKey& key) {
    SeriesState st;
    const md_script_eval_t* eval = is_script(key.source) ? app->script.eval : nullptr;
    if (!eval) return st;
    const str_t ident = series_script_ident(key);

    // Against the script in the editor, when it compiles and is not the one evaluated: computed from the same thing,
    // or not (md_script_ir_property_fingerprint). A property neither has is nothing to speak of: nothing is shown.
    const md_script_ir_t* ir = md_script_ir_valid(app->script.ir) ? app->script.ir : app->script.eval_ir;
    const uint64_t evaluated = md_script_eval_property_fingerprint(eval, ident);
    const uint64_t now = ir ? md_script_ir_property_fingerprint(ir, ident) : evaluated;
    if (!evaluated && !now) return st;
    if (evaluated && !now) {
        st.status = SeriesStatus_Removed;
        return st;
    }
    if (evaluated && now != evaluated) {
        st.status = SeriesStatus_Outdated;
        return st;
    }

    // How far it has come: over every frame, or over its period's
    uint32_t beg = 0, end = 0;
    const bool over_period = key.period != 0;
    st.num_wanted = key_period_frames(&beg, &end, app, key) ? end - beg : 0;
    const uint32_t kinds = (uint32_t)md_script_ir_property_flags(app->script.eval_ir, ident);
    if (over_period && (kinds & (MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION | MD_SCRIPT_PROPERTY_FLAG_VOLUME))) {
        // Its slot over the period's frames. What it has yet to do is the frames it wants and does not hold, and those
        // it holds and does not want: of what it holds, (wanted + held - pending) / 2 are wanted.
        md_script_eval_slot_info_t info;
        if (md_script_eval_slot_info(&info, eval, SCRIPT_SLOT_PERIOD(key.period), ident)) {
            const uint32_t both = (st.num_wanted + info.num_frames - MIN(info.num_pending, st.num_wanted + info.num_frames)) / 2;
            st.num_frames = MIN(both, st.num_wanted);
            st.num_extra  = info.num_frames - MIN(both, info.num_frames);
            if (info.num_pending > 0) st.status = SeriesStatus_Partial;
        }
        return st;
    }
    md_bitset_t held = md_script_eval_property_frames(eval, ident, app->allocator.frame);
    if (over_period) {
        held = md_bitset_and(held, st.num_wanted ? md_bitset_range(beg, end) : md_bitset_empty_set(), app->allocator.frame);
    }
    st.num_frames = (uint32_t)md_bitset_count(held);
    if (st.num_frames < st.num_wanted) st.status = SeriesStatus_Partial;
    return st;
}

void series_state_suffix(char* buf, size_t cap, const SeriesState& state) {
    switch (state.status) {
    case SeriesStatus_Partial:
        if (state.num_frames < state.num_wanted) snprintf(buf, cap, " \xC2\xB7 %u/%u", state.num_frames, state.num_wanted);
        else snprintf(buf, cap, " \xC2\xB7 %u to take out", state.num_extra);
        break;
    case SeriesStatus_Outdated: snprintf(buf, cap, " \xC2\xB7 outdated"); break;
    case SeriesStatus_Removed:  snprintf(buf, cap, " \xC2\xB7 removed"); break;
    default: if (cap) buf[0] = '\0'; break;
    }
}

void series_state_describe(char* buf, size_t cap, const SeriesState& state) {
    switch (state.status) {
    case SeriesStatus_Partial:
        if (state.num_extra > 0) {
            snprintf(buf, cap, "Holds %u of its %u frames so far, and %u frames it no longer covers,\nwhich are taken out as they are evaluated: Evaluate to finish it.",
                     state.num_frames, state.num_wanted, state.num_extra);
        } else {
            snprintf(buf, cap, "Evaluated over %u of its %u frames so far", state.num_frames, state.num_wanted);
        }
        break;
    case SeriesStatus_Outdated:
        snprintf(buf, cap, "Shown as evaluated from an earlier version of the script.\nWhat it is computed from has changed: Evaluate to update it.");
        break;
    case SeriesStatus_Removed:
        snprintf(buf, cap, "Shown as last evaluated: it is no longer in the script.");
        break;
    default:
        if (cap) buf[0] = '\0';
        break;
    }
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
    if (key.period != 0) {
        snprintf(buf, cap, STR_FMT "%s (%s)", STR_ARG(base), variant_label[v], series_period_name(app, key.period));
    } else {
        snprintf(buf, cap, STR_FMT "%s", STR_ARG(base), variant_label[v]);
    }
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

void series_set_drag_payload(const char* dnd_type, const SeriesKey& key, int src_subplot, const char* label, ImVec4 color) {
    SeriesDragPayload payload;
    payload.key = key;
    payload.src_subplot = src_subplot;
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
    char base[64], suffix[32];
    series_label(base, sizeof(base), app, key);
    out->state = series_state(app, key);
    series_state_suffix(suffix, sizeof(suffix), out->state);
    snprintf(out->label, sizeof(out->label), "%s%s", base, suffix);
    snprintf(out->plot_id, sizeof(out->plot_id), "%s###%d:%d:%u:%s", out->label, (int)key.source, (int)key.variant, key.period, key.path);
    if (is_script(key.source)) {
        const str_t ident = series_script_ident(key);
        str_copy_to_char_buf(out->script_ident, sizeof(out->script_ident), ident);
        out->vis_ref = md_script_ir_property_vis_ref(app->script.eval_ir, ident);
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

// What the slot of a distribution or volume holds, copied out of it whenever it holds something new
// (md_script_eval_slot_extract): over a period (SCRIPT_SLOT_PERIOD)
struct SlotEntry {
    uint64_t        key = 0;
    char            ident[64] = "";
    uint64_t        version = 0;
    md_script_eval_slot_info_t info = {};
    md_array(float) values = nullptr;
    md_array(float) weights = nullptr;
    int             last_used = 0;
};

static md_array(ConvertEntry)   convert_cache = nullptr;
static md_array(HistogramEntry) histogram_cache = nullptr;
static md_array(SlotEntry)      slot_cache = nullptr;

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

// The slot of the script property under key, as it is now. NULL while there is none: not asked for, or not a
// distribution or a volume. The values stay put until the slot is next read. A period's slot is asked for by what shows
// it (series_period_slot), and let go when nothing has read it for a while (series_cache_gc).
static const SlotEntry* series_slot(ApplicationState* app, uint64_t key, str_t ident) {
    const md_script_eval_t* eval = app->script.eval;
    md_script_eval_slot_info_t info;
    if (!eval || str_empty(ident) || !md_script_eval_slot_info(&info, eval, key, ident)) return nullptr;

    SlotEntry* entry = nullptr;
    for (size_t i = 0; i < md_array_size(slot_cache); ++i) {
        if (slot_cache[i].key == key && str_eq(str_from_cstr(slot_cache[i].ident), ident)) {
            entry = &slot_cache[i];
            break;
        }
    }
    if (!entry) {
        md_array_push(slot_cache, SlotEntry{}, app->allocator.persistent);
        entry = &slot_cache[md_array_size(slot_cache) - 1];
        entry->key = key;
        str_copy_to_char_buf(entry->ident, sizeof(entry->ident), ident);
    }
    entry->last_used = ImGui::GetFrameCount();

    if (entry->version != info.version || md_array_size(entry->values) != info.num_values) {
        md_array_resize(entry->values,  info.num_values, app->allocator.persistent);
        md_array_resize(entry->weights, info.num_values, app->allocator.persistent);
        if (!md_script_eval_slot_extract(&entry->info, entry->values, entry->weights, info.num_values, eval, key, ident)) {
            entry->version = 0;
            return nullptr;
        }
        entry->version = entry->info.version;
    }
    return entry;
}

// The slot of a script property over a period, asked for over the frames [beg, end[ it covers now: moved, the period
// takes its slot along, and the slot keeps what it has in common with its new frames
static const SlotEntry* series_period_slot(ApplicationState* app, str_t ident, uint32_t period, uint32_t beg, uint32_t end) {
    const md_bitset_t frames = md_bitset_range(beg, end);
    if (!app->script.eval || !md_script_eval_slot_request(app->script.eval, SCRIPT_SLOT_PERIOD(period), ident, &frames)) return nullptr;
    return series_slot(app, SCRIPT_SLOT_PERIOD(period), ident);
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
    for (size_t i = 0; i < md_array_size(slot_cache);) {
        if (frame - slot_cache[i].last_used > 600) {
            if (app->script.eval) {
                // Nothing shows it any more
                md_script_eval_slot_release(app->script.eval, slot_cache[i].key, str_from_cstr(slot_cache[i].ident));
            }
            md_array_free(slot_cache[i].values,  app->allocator.persistent);
            md_array_free(slot_cache[i].weights, app->allocator.persistent);
            md_array_swap_back_and_pop(slot_cache, i);
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

    for (size_t i = 0; i < md_array_size(slot_cache); ++i) {
        md_array_free(slot_cache[i].values,  app->allocator.persistent);
        md_array_free(slot_cache[i].weights, app->allocator.persistent);
    }
    md_array_free(slot_cache, app->allocator.persistent);
    slot_cache = nullptr;
}

// ## What a series can be shown as

uint32_t series_views(const ApplicationState* app, const SeriesKey& key, const char** out_reason) {
    const char* reason_buf = nullptr;
    const char** reason = out_reason ? out_reason : &reason_buf;
    *reason = nullptr;

    const md_attributes_t* table = series_table(app, key.source, str_from_cstr(key.path));
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr) {
        *reason = "Nothing is published at this path";
        return 0;
    }
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {
        *reason = "Text, not numbers";
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
            // Over a period it is the same over time, so only the one over every frame is offered for that
            if (key.period == 0 && n == md_array_size(app->timeline.x_values)) views |= SeriesView_Timeline;
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
        if (!(views & SeriesView_Timeline) && key.period == 0) {
            *reason = "Its frame axis cannot be placed on the timeline (no trajectory, or time without a unit on one side)";
        }
        return views;
    }

    if (is_script(key.source) && (script_property_kind(app, key) & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) && attribute_numeric_resident(attr)) {
        if (attr->format.rank != 1) {
            *reason = "An array of distributions, [n][bins]: not plotted yet";
            return 0;
        }
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
    // Over a period it is the same over time: shown over every frame
    if (key.period != 0) return false;

    const str_t path = str_from_cstr(key.path);
    const md_attributes_t* table = series_table(app, key.source, path);
    if (!table) return false;

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
static void bin_temporal(float* bins, int num_bins, int dim_out, double lo, double hi, const float* values, size_t num_frames, size_t stride, const md_bitset_t* mask) {
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
        md_bitset_iter_t it = md_bitset_iter(*mask);
        while (md_bitset_iter_next(&it)) {
            const size_t f = it.idx;
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

bool series_period_frames(uint32_t* beg, uint32_t* end, const ApplicationState* app, const SeriesKey& key, uint32_t period) {
    const uint32_t num_frames = app->script.eval ? (uint32_t)md_script_eval_frame_count(app->script.eval) : (uint32_t)run_num_frames(app);
    if (num_frames == 0) return false;
    // The key's own: every frame, or its period's
    if (period == 0) period = key.period;
    if (period == 0) {
        *beg = 0;
        *end = num_frames;
        return true;
    }
    const TimelinePeriod* r = timeline_periods_find(&app->timeline.periods, period);
    int b = 0, e = 0;
    if (!r || !timeline_period_frames(&b, &e, *r, app->animation.frame, (int)num_frames)) return false;
    *beg = (uint32_t)b;
    *end = (uint32_t)e + 1;
    return true;
}

bool series_resolve_histogram(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins, uint32_t period) {
    ASSERT(out && app);
    *out = SeriesHistogramView{};
    fill_names(out, app, key);

    if (key.path[0] == '\0') return false;
    if (key.variant != SeriesVariant_Values && key.variant != SeriesVariant_Aggregate) return false;
    // Over a period, the one asked for or the key's own: a script property's frames, which are the trajectory's
    if (period == 0) period = key.period;
    if (period != 0 && !is_script(key.source)) return false;
    uint32_t period_beg = 0, period_end = 0;
    if (period != 0 && !series_period_frames(&period_beg, &period_end, app, key, period)) return false;

    const str_t path = str_from_cstr(key.path);
    const md_attributes_t* table = series_table(app, key.source, path);
    if (!table) return false;

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

    // What is binned: a temporal attribute over the frames that count (the frames the evaluation holds, of the period's),
    // a distribution as published, or over the period's frames as its slot holds it. A period's slot is asked for here,
    // for the frames the period covers now: moved, the period takes its slot along.
    md_bitset_t frames = md_bitset_empty_set();
    const md_bitset_t* mask = nullptr;
    if (temporal) {
        SeriesKey over = key;
        over.period = 0;
        if (series_frame_mask(&frames, app, over, app->allocator.frame)) {
            if (period != 0) frames = md_bitset_and(frames, md_bitset_range(period_beg, period_end), app->allocator.frame);
            mask = &frames;
        }
    }
    const SlotEntry* slot = nullptr;
    if (distribution && period != 0) {
        slot = series_period_slot(app, series_script_ident(key), period, period_beg, period_end);
        if (!slot || md_array_size(slot->values) != md_attribute_element_count(&attr->format)) return false;
    }

    // ### Cached bins
    uint64_t hkey = md_hash64_str(path, (uint64_t)key.source * 31 + key.variant);
    hkey = md_hash64(&num_bins, sizeof(num_bins), hkey);
    hkey = md_hash64(&period, sizeof(period), hkey);
    // Of what is binned: the values and weights, and which frames count
    uint64_t version = slot ? slot->version : attr->version + (weight ? weight->version << 1 : 0);
    if (mask) version = md_bitset_hash64(*mask, version);

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

        const float* values = slot ? slot->values : attribute_floats(app, table, attr);
        if (values) {
            double range[2] = {0, 0};
            const bool has_range = read_range(range, md_attributes_find_in(table, path, STR_LIT("range")));
            if (temporal) {
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
                const float* w = slot ? slot->weights : weight ? attribute_floats(app, table, weight) : nullptr;
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

bool series_resolve_histogram_periods(SeriesHistogramView* out, ApplicationState* app, const SeriesKey& key, int num_bins,
                                      const uint32_t* periods, int num_periods) {
    ASSERT(out && app);
    const bool over_periods = is_script(key.source) && (num_periods > 1 || (num_periods == 1 && periods[0] != 0));
    if (!over_periods) {
        return series_resolve_histogram(out, app, key, num_bins);
    }

    // A population is pooled: one curve per period, not one per member and period
    SeriesKey k = key;
    {
        SeriesHistogramView probe;
        if (series_resolve_histogram(&probe, app, key, num_bins) && probe.dim > 1) {
            k.variant = SeriesVariant_Aggregate;
        }
    }

    num_periods = MIN(num_periods, PLOT_MAX_PERIODS);
    SeriesHistogramView views[PLOT_MAX_PERIODS];
    bool any = false;
    int first = -1;
    for (int p = 0; p < num_periods; ++p) {
        const bool ok = series_resolve_histogram(&views[p], app, k, num_bins, periods[p]);
        any |= ok;
        if (ok && first < 0) first = p;
    }

    // Named as the series is, whichever resolved
    *out = views[first >= 0 ? first : 0];
    out->num_periods = num_periods;
    out->dim = num_periods;
    if (!any) {
        out->bins = nullptr;
        out->num_bins = 0;
    }

    // Every period binned alike: the bins of each, one after the other
    const SeriesHistogramView& ref = views[first >= 0 ? first : 0];
    float* bins = any ? (float*)md_alloc(app->allocator.frame, (size_t)ref.num_bins * num_periods * sizeof(float)) : nullptr;
    double y_min = DBL_MAX, y_max = -DBL_MAX;
    const md_script_eval_t* eval = app->script.eval;
    const str_t ident = series_script_ident(key);
    const bool distribution = script_property_kind(app, key) & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION;
    for (int p = 0; p < num_periods; ++p) {
        SeriesPeriodInfo& info = out->periods[p];
        info = SeriesPeriodInfo{};
        info.id = periods[p];
        const SeriesHistogramView& v = views[p];
        info.valid = any && v.bins && v.num_bins == ref.num_bins && v.x_min == ref.x_min && v.x_max == ref.x_max;
        if (bins) {
            float* dst = bins + (size_t)p * ref.num_bins;
            if (info.valid) {
                MEMCPY(dst, v.bins, (size_t)ref.num_bins * sizeof(float));
                y_min = MIN(y_min, v.y_min);
                y_max = MAX(y_max, v.y_max);
            } else {
                MEMSET(dst, 0, (size_t)ref.num_bins * sizeof(float));
            }
        }

        // What it is called and coloured, and how much of it there is. The key's own period in the series' colour.
        const uint32_t period = periods[p] != 0 ? periods[p] : key.period;
        if (period == 0) {
            snprintf(info.name, sizeof(info.name), "All frames");
        } else if (const TimelinePeriod* r = timeline_periods_find(&app->timeline.periods, period)) {
            snprintf(info.name, sizeof(info.name), "%s", r->name);
            if (periods[p] != 0) {
                info.color = ImVec4(r->color.x, r->color.y, r->color.z, r->color.w);
                info.own_color = true;
            }
        }
        uint32_t beg = 0, end = 0;
        if (eval && series_period_frames(&beg, &end, app, key, period)) {
            info.num_wanted = end - beg;
            if (distribution && period != 0) {
                md_script_eval_slot_info_t si;
                if (md_script_eval_slot_info(&si, eval, SCRIPT_SLOT_PERIOD(period), ident)) {
                    // Of what it holds, what it wants (see series_state)
                    const uint32_t both = (info.num_wanted + si.num_frames - MIN(si.num_pending, info.num_wanted + si.num_frames)) / 2;
                    info.num_frames = MIN(both, info.num_wanted);
                }
            } else {
                const md_bitset_t held = md_script_eval_property_frames(eval, ident, app->allocator.frame);
                info.num_frames = (uint32_t)md_bitset_count(md_bitset_and(held, md_bitset_range(beg, end), app->allocator.frame));
            }
        }
    }
    out->bins = bins;
    if (any) {
        out->y_min = y_min;
        out->y_max = y_max;
    }
    return any;
}

// ## Volumes

// A volume is {x, y, z}, an array of n of them {n, x, y, z}: the member is a slice of the first axis. 0 for anything else.
static uint32_t volume_shape(uint32_t dim[3], const md_attribute_t* attr) {
    if (!attr || attr->format.type != MD_ATTRIBUTE_TYPE_F32 || attr->format.components != 1) return 0;
    const uint32_t rank = attr->format.rank;
    if (rank != 3 && rank != 4) return 0;
    const uint32_t* shape = attr->format.shape + (rank - 3);
    for (int i = 0; i < 3; ++i) {
        if (dim) dim[i] = shape[i];
    }
    return rank == 4 ? attr->format.shape[0] : 1;
}

uint32_t series_volume_count(const ApplicationState* app, const SeriesKey& key) {
    if (!is_script(key.source)) return 0;
    const md_attributes_t* table = series_table(app, key.source, str_from_cstr(key.path));
    return volume_shape(nullptr, table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr);
}

bool series_resolve_volume(SeriesVolumeView* out, ApplicationState* app, const SeriesKey& key, uint32_t member) {
    ASSERT(out && app);
    *out = SeriesVolumeView{};
    series_label(out->label, sizeof(out->label), app, key);
    if (!is_script(key.source)) return false;
    const str_t ident = series_script_ident(key);
    str_copy_to_char_buf(out->script_ident, sizeof(out->script_ident), ident);
    out->vis_ref = md_script_ir_property_vis_ref(app->script.eval_ir, ident);
    out->state = series_state(app, key);

    const md_attributes_t* table = series_table(app, key.source, str_from_cstr(key.path));
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr || !attr->data) return false;
    out->count = volume_shape(out->dim, attr);
    if (out->count == 0) return false;
    out->member = MIN(member, out->count - 1);
    if (out->count > 1) {
        // As the script names it: v[2]
        if (key.period != 0) {
            snprintf(out->label, sizeof(out->label), STR_FMT "[%u] (%s)", STR_ARG(ident), out->member + 1, series_period_name(app, key.period));
        } else {
            snprintf(out->label, sizeof(out->label), STR_FMT "[%u]", STR_ARG(ident), out->member + 1);
        }
    }
    const size_t member_values = (size_t)out->dim[0] * out->dim[1] * out->dim[2];

    const float* values = nullptr;
    if (key.period != 0) {
        // Over the period's frames: as its slot holds it, all of the volumes
        uint32_t beg, end;
        if (!key_period_frames(&beg, &end, app, key)) return false;
        const SlotEntry* slot = series_period_slot(app, ident, key.period, beg, end);
        if (!slot || md_array_size(slot->values) != md_attribute_element_count(&attr->format)) return false;
        values = slot->values;
        out->version = slot->version;
    } else {
        values = (const float*)attr->data;
        out->version = attr->version;
    }
    out->data = values + out->member * member_values;
    return true;
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

    char path[SERIES_PATH_CAP];
    const int len = snprintf(path, sizeof(path), STR_FMT STR_FMT, STR_ARG(script_prefix), STR_ARG(name));
    if (len <= 0 || (size_t)len >= sizeof(path)) return false;

    const md_attributes_t* table = series_table(app, source, str_t{path, (size_t)len});
    if (!table) return false;
    if (!md_attributes_find(table, str_t{path, (size_t)len})) return false;

    *out = series_key(source, str_t{path, (size_t)len});
    return true;
}

// ## Subplots

// The key's own period is period 0: the frames of the key
static uint32_t own_period(const PlotSeries& s, uint32_t id) {
    return id == s.key.period ? 0 : id;
}

bool plot_series_has_period(const PlotSeries& s, uint32_t id) {
    id = own_period(s, id);
    for (int i = 0; i < s.num_periods; ++i) {
        if (s.periods[i] == id) return true;
    }
    return false;
}

bool plot_series_add_period(PlotSeries& s, uint32_t id) {
    id = own_period(s, id);
    if (plot_series_has_period(s, id) || s.num_periods >= PLOT_MAX_PERIODS) return false;
    s.periods[s.num_periods++] = id;
    return true;
}

void plot_series_remove_period(PlotSeries& s, uint32_t id) {
    id = own_period(s, id);
    if (s.num_periods <= 1) return;
    for (int i = 0; i < s.num_periods; ++i) {
        if (s.periods[i] != id) continue;
        for (int j = i + 1; j < s.num_periods; ++j) s.periods[j - 1] = s.periods[j];
        s.periods[--s.num_periods] = 0;
        return;
    }
}

int plot_find_series(const PlotSubplot& sp, const SeriesKey& key) {
    for (int i = 0; i < sp.count; ++i) {
        if (series_key_equal(sp.series[i].key, key)) return i;
    }
    return -1;
}

int plot_add_series(const ApplicationState* app, PlotSubplot& sp, const SeriesKey& key) {
    if (key.path[0] == '\0') return -1;
    const int existing = plot_find_series(sp, key);
    if (existing != -1) return existing;
    if (sp.count >= PLOT_MAX_SERIES_PER_SUBPLOT) {
        MD_LOG_INFO("A subplot holds at most %d series", PLOT_MAX_SERIES_PER_SUBPLOT);
        return -1;
    }
    PlotSeries& s = sp.series[sp.count];
    s = PlotSeries{};
    s.key = key;
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

void plot_resolve_legacy_periods(const ApplicationState* app, PlotSubplot* subplots, int count) {
    for (int i = 0; i < count; ++i) {
        PlotSubplot& sp = subplots[i];
        for (int j = 0; j < sp.count;) {
            PlotSeries& s = sp.series[j];
            if (s.key.period != SERIES_PERIOD_LEGACY_FILTER) { ++j; continue; }
            s.key.period = series_legacy_period(app, s.key.period);
            // Its periods again, of which the period it now has of its own is its period 0
            const PlotSeries was = s;
            s.num_periods = 0;
            for (int k = 0; k < was.num_periods; ++k) plot_series_add_period(s, was.periods[k]);
            if (s.num_periods == 0) s.num_periods = 1;
            const int first = plot_find_series(sp, s.key);
            if (first != j) {
                plot_remove_series(sp, j);
            } else {
                ++j;
            }
        }
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
            if (s.key.period != 0) {
                viamd::write_int(state, STR_LIT("Period"), (int)s.key.period);
            }
            viamd::write_vec4(state, STR_LIT("Color"),   vec4_set(s.color.x, s.color.y, s.color.z, s.color.w));
            viamd::write_str (state, STR_LIT("PlotType"), str_from_cstr(plot_type_name[s.plot_type < PlotType_Count ? s.plot_type : 0]));
            viamd::write_int (state, STR_LIT("Marker"),  (int)s.marker);
            viamd::write_flt (state, STR_LIT("MarkerSize"), s.marker_size);
            viamd::write_flt (state, STR_LIT("BarWidth"), s.bar_width);
            viamd::write_int (state, STR_LIT("NumBins"), s.num_bins);
            viamd::write_bool(state, STR_LIT("UseColormap"), s.use_colormap);
            viamd::write_int (state, STR_LIT("Colormap"), (int)s.colormap);
            viamd::write_flt (state, STR_LIT("ColormapAlpha"), s.colormap_alpha);
            if (plot_series_has_periods(s)) {
                int periods[PLOT_MAX_PERIODS];
                for (int k = 0; k < s.num_periods; ++k) periods[k] = (int)s.periods[k];
                viamd::write_int_vec(state, STR_LIT("Periods"), periods, s.num_periods);
            }
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
                series_source_from_name(&s.key.source, str, &s.key.period);
            } else if (str_eq_cstr(ident, "Period") && viamd::extract_int(i, arg)) {
                s.key.period = (uint32_t)MAX(i, 0);
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
            } else if (str_eq_cstr(ident, "Periods")) {
                // Period ids, as many as there are, separated by commas
                int num = 0;
                uint32_t value = 0;
                bool digits = false;
                for (size_t c = 0; c <= arg.len; ++c) {
                    const char ch = c < arg.len ? arg.ptr[c] : ',';
                    if (ch >= '0' && ch <= '9') {
                        value = value * 10 + (uint32_t)(ch - '0');
                        digits = true;
                    } else if (ch == ',') {
                        if (digits && num < PLOT_MAX_PERIODS) s.periods[num++] = value;
                        value = 0;
                        digits = false;
                    }
                }
                if (num > 0) s.num_periods = num;
            }
        }

        if (0 <= subplot && subplot < PLOT_MAX_SUBPLOTS) {
            const int idx = plot_add_series(app, subplots[subplot], s.key);
            if (idx != -1) {
                PlotSeries& dst = subplots[subplot].series[idx];
                const ImVec4 default_color = dst.color;
                dst = s;
                if (!has_color) dst.color = default_color;
            }
        }
        return true;
    }

    return false;
}
