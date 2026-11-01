#include "nci_surface_model.h"

#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <gfx/gl_utils.h>

#include <algorithm>
#include <iterator>
#include <cstdio>
#include <string>

namespace topology {

NciFamilySettings nci_family_default_settings(nci_surface::Family family) {
    const nci_surface::Defaults defaults = nci_surface::family_defaults(nci_surface::family_bit(family));
    NciFamilySettings settings;
    settings.s_iso = (float)defaults.s_iso;
    settings.rho_max = (float)defaults.rho_max;
    settings.radius = (float)defaults.radius;
    return settings;
}

NciUniformSettings nci_uniform_default_settings() {
    const nci_surface::UniformDefaults defaults = nci_surface::uniform_defaults();
    return {(float)defaults.s_iso, (float)defaults.rho_max};
}

ColorScale nci_default_scale() {
    const float range = (float)nci_surface::uniform_defaults().color_range;
    return {.colormap = ImPlotColormap_Jet, .range_beg = -range, .range_end = range,
            .symmetric = true, .auto_range = true, .show_legend = true};
}

namespace {

size_t family_index(const NciSurface& n) { return nci_settings_index(n); }

template <typename F>
void for_each_surface(NciState& state, F&& f) {
    f(state.primary);
    for (auto& n : state.full_system) f(*n);
    for (auto& n : state.previews) f(*n);
    for (auto& n : state.objects) f(*n);
}

void apply_visibility(NciState& state) {
    for_each_surface(state, [&](NciSurface& n) {
        n.layer_visible = !state.layer_hidden && nci_layer_key(n) == state.selected_layer;
    });
}

std::vector<int> saved_layer_keys(const NciState& state) {
    std::vector<int> keys;
    for (const auto& n : state.objects) {
        const int key = nci_layer_key(*n);
        if (std::find(keys.begin(), keys.end(), key) == keys.end()) keys.push_back(key);
    }
    return keys;
}

}  // namespace

int nci_layer_key(const NciSurface& n) {
    if (n.saved) return n.object_group ? n.object_group : n.id;
    if (n.family_preview) return NCI_LAYER_PREVIEW;
    return NCI_LAYER_FULL_SYSTEM;
}

bool nci_scale_member(const NciSurface& n) {
    return n.shown && n.layer_visible;
}

ColorScale nci_effective_scale(const NciState& state, const NciSurface& n) {
    return n.id == -1 ? n.scale : state.scale;
}

void nci_select_layer(NciState& state, int key) {
    state.selected_layer = key;
    state.layer_hidden = false;
    state.view = key == NCI_LAYER_FULL_SYSTEM ? nci_surface::View::Overview
               : key == NCI_LAYER_PREVIEW     ? nci_surface::View::Preview
                                              : nci_surface::View::Saved;
    const auto members = nci_layer_members(state, key);
    state.active_id = members.empty() ? 0 : members.front()->id;
    apply_visibility(state);
}

void nci_select(NciState& state, int id) {
    int key = NCI_LAYER_FULL_SYSTEM;
    for_each_surface(state, [&](NciSurface& n) { if (id != 0 && n.id == id) key = nci_layer_key(n); });
    nci_select_layer(state, key);
    if (id != 0) state.active_id = id;
}

std::vector<NciSurface*> nci_layer_members(NciState& state, int key) {
    std::vector<NciSurface*> members;
    if (key == NCI_LAYER_FULL_SYSTEM) {
        for (auto& n : state.full_system) if (n->shown) members.push_back(n.get());
        if (members.empty() && state.primary.shown) members.push_back(&state.primary);
    } else if (key == NCI_LAYER_PREVIEW) {
        for (auto& n : state.previews) if (n->shown) members.push_back(n.get());
    } else if (key != NCI_LAYER_NONE) {
        for (auto& n : state.objects) if (nci_layer_key(*n) == key) members.push_back(n.get());
    }
    return members;
}

size_t nci_settings_index(const NciSurface& n) {
    return n.families ? (size_t)n.group : (size_t)nci_surface::Family::Unassigned;
}

std::string nci_layer_types(const std::vector<NciSurface*>& members) {
    std::array<uint32_t, nci_surface::family_count> counts = {};
    bool whole_system = false;
    for (const NciSurface* n : members) {
        whole_system |= !n->families;
        for (size_t f = 0; f < counts.size(); ++f) counts[f] += n->type_counts[f];
    }
    std::string out = whole_system ? "Whole system" : "";
    for (size_t f = 0; f < counts.size(); ++f) {
        if (!counts[f]) continue;
        char text[96];
        snprintf(text, sizeof(text), "%s%s %u", out.empty() ? "" : ", ", nci_surface::family_name((nci_surface::Family)f), counts[f]);
        out += text;
    }
    return out.empty() ? std::string("No CPs") : out;
}

std::string nci_family_label(const NciSurface& n) {
    return n.families ? nci_surface::family_name(n.group) : "Whole system";
}

void nci_apply_settings(NciState& state) {
    apply_visibility(state);
    const auto apply = [&](const std::vector<NciSurface*>& members) {
        for (NciSurface* n : members) {
            const NciFamilySettings& settings = state.family_settings[family_index(*n)];
            // Without classified CPs the whole system has no CP type to adapt to
            const bool adaptive = state.mode == NciParameterMode::CPAdaptive && n->families;
            n->s_iso = adaptive ? settings.s_iso : state.uniform.s_iso;
            n->rho_max = adaptive ? settings.rho_max : state.uniform.rho_max;
            n->near_radius = settings.radius;
            n->resolution = state.resolution;
            n->connected = state.connected;
            n->radius_cap = state.radius_cap;
            n->opacity = state.opacity;
            n->editing |= state.settings_editing;
        }
    };
    apply(nci_layer_members(state, NCI_LAYER_FULL_SYSTEM));
    apply(nci_layer_members(state, NCI_LAYER_PREVIEW));
    for (int key : saved_layer_keys(state)) apply(nci_layer_members(state, key));
}

void nci_update_scale(NciState& state) {
    state.span = {};
    for_each_surface(state, [&](NciSurface& n) {
        if (!nci_scale_member(n) || !n.have_key || n.cancelled || n.failed || (!n.saved && n.want_key != n.have_key)) return;
        color_scale_include_span(&state.span, n.span);
    });
    // With no shown surface Auto keeps the last range rather than jumping back to the default
    color_scale_update_range(&state.scale, state.span);
}

bool nci_surface_visible(const NciState& state, const NciSurface& n) {
    return state.enabled && nci_scale_member(n) && n.want_key != 0 && n.want_key == n.have_key &&
           !n.cancelled && !n.failed && n.triangles > 0 && n.mesh.index_count > 0;
}

void nci_surface_release(NciSurface& n) {
    nci_surface_job_destroy(n.job);
    n.job = nullptr;
    n.rho.reset();
    n.rho_key = n.want_key = n.have_key = 0;
    n.triangles = 0;
    n.span = {};
    n.scatter = {};
    n.source_mesh = {};
    n.outline.clear();
    n.cancelled = n.rho_failed = n.failed = false;
    isomesh::mesh_free(&n.mesh);
    gl::free_texture(&n.colormap);
    n.uploaded_colormap = -1;
}

void nci_state_free(NciState& state) {
    nci_surface_release(state.primary);
    nci_surface_release(state.context);
    for (auto& n : state.objects) nci_surface_release(*n);
    state.objects.clear();
    for (auto& n : state.previews) nci_surface_release(*n);
    state.previews.clear();
    for (auto& n : state.full_system) nci_surface_release(*n);
    state.full_system.clear();
    state.layer_outlines.clear();
    state.scatter_layers.clear();
    state.view = nci_surface::View::Overview;
    state.classified_overview = true;
    state.selection_key = 0;
    state.full_selection_key = 0;
    state.force_all_cps = false;
    state.next_group_id = 1;
    state.active_id = 0;
    state.selected_layer = NCI_LAYER_FULL_SYSTEM;
    state.layer_hidden = false;
    state.layer_order.clear();
    state.next_id = 1;
    state.primary.families = 0;
    state.primary.environments.clear();
    state.span = {};
}

void nci_sync_preview(NciState& state, const NciPreviewInput& input) {
    const uint64_t cp_selection_key = input.selection_key;
    if (state.force_all_cps && cp_selection_key != state.full_selection_key) state.force_all_cps = false;
    const bool selected = input.enabled && input.has_selection && !state.force_all_cps;
    uint64_t key = input.selection_key;
    key = md_hash64(&input.enabled, sizeof(input.enabled), key);
    key = md_hash64(&state.force_all_cps, sizeof(state.force_all_cps), key);

    std::vector<nci_surface::CriticalPointAnchor> anchors;
    anchors.reserve(input.enabled ? input.critical_points.size() : 0);
    // Grouped by parameters: the CPs of one interaction region make one surface, cut at one threshold
    std::array<std::vector<std::vector<std::string>>, nci_surface::family_count> all_environments, selected_environments;
    using TypeCounts = std::array<uint32_t, nci_surface::family_count>;
    std::array<TypeCounts, nci_surface::family_count> all_types = {}, selected_types = {};
    for (const NciCriticalPoint& cp : input.critical_points) {
        if (!nci_surface::has_nci(cp.parameter_family)) continue;
        const size_t group = (size_t)cp.parameter_family;
        key = md_hash64(&cp.family, sizeof(cp.family), key);
        key = md_hash64(&cp.parameter_family, sizeof(cp.parameter_family), key);
        anchors.push_back({cp.position, cp.parameter_family, cp.selected});
        all_types[group][(size_t)cp.family]++;
        if (cp.selected) selected_types[group][(size_t)cp.family]++;
        if (!cp.environment.empty()) {
            all_environments[group].push_back(cp.environment);
            if (cp.selected) selected_environments[group].push_back(cp.environment);
        }
    }

    const nci_surface::CriticalPointGroups anchor_groups = nci_surface::group_critical_points(anchors);
    const auto& all_anchors = anchor_groups.all_anchors;
    const auto& selected_anchors = anchor_groups.selected_anchors;
    for (const auto& near : all_anchors)
        if (!near.empty()) key = md_hash64(near.data(), near.size() * sizeof(double), key);
    for (const auto& near : selected_anchors)
        if (!near.empty()) key = md_hash64(near.data(), near.size() * sizeof(double), key);
    key |= 1;
    if (!nci_surface::follow_selection(&state.view, &state.selection_key, key, selected)) return;

    state.classified_overview = false;
    for (size_t f = 0; f < all_anchors.size(); ++f) {
        state.classified_overview |= !all_anchors[f].empty();
        const auto update_family = [&](std::vector<std::unique_ptr<NciSurface>>& surfaces,
                                       const std::array<std::vector<double>, nci_surface::family_count>& family_anchors,
                                       const std::array<std::vector<std::vector<std::string>>, nci_surface::family_count>& environments,
                                       const std::array<TypeCounts, nci_surface::family_count>& types, int id, bool full_system) {
            NciSurface* surface = nullptr;
            for (auto& candidate : surfaces) if (candidate->id == id) { surface = candidate.get(); break; }
            if (!surface && family_anchors[f].empty()) return;
            if (!surface) {
                auto created = std::make_unique<NciSurface>();
                created->family_preview = !full_system;
                created->full_system_family = full_system;
                created->id = id;
                created->shown = true;
                created->scope = NciSurface::Scope_Selected;
                created->group = (nci_surface::Family)f;
                created->resolution = state.resolution;
                created->generation = state.primary.generation;
                created->chip = nci_layer_colour(full_system ? 0 : 1);
                snprintf(created->name, sizeof(created->name), "%s", full_system ? "Full system" : "Preview");
                surface = created.get();
                surfaces.push_back(std::move(created));
            }
            surface->shown = !family_anchors[f].empty();
            surface->anchors = family_anchors[f];
            surface->environments = environments[f];
            surface->type_counts = types[f];
            surface->families = 0;
            for (size_t t = 0; t < types[f].size(); ++t) if (types[f][t]) surface->families |= 1u << t;
        };
        update_family(state.full_system, all_anchors, all_environments, all_types, -100 - (int)f, true);
        update_family(state.previews, selected_anchors, selected_environments, selected_types, -2 - (int)f, false);
    }

    state.primary.shown = !state.classified_overview;
    state.primary.scope = NciSurface::Scope_Everywhere;
    // A new CP selection shows its preview; clearing it returns from the preview to the full system.
    // A saved layer stays shown while the selection is cleared.
    // A layer the user hid stays hidden.
    if (selected) {
        nci_select_layer(state, NCI_LAYER_PREVIEW);
    } else if (state.selected_layer == NCI_LAYER_PREVIEW || (state.selected_layer > 0 && state.layer_hidden)) {
        nci_select_layer(state, NCI_LAYER_FULL_SYSTEM);
    } else {
        const auto members = nci_layer_members(state, state.selected_layer);
        if (std::none_of(members.begin(), members.end(), [&](const NciSurface* n) { return n->id == state.active_id; }))
            state.active_id = members.empty() ? 0 : members.front()->id;
        apply_visibility(state);
    }
}

// A computed result, current or, while the layer is hidden (want_key 0), the last one
static bool outline_member(const NciSurface& n) {
    return n.shown && n.have_key && (!n.want_key || n.have_key == n.want_key) && !n.cancelled && !n.failed;
}

const char* nci_scatter_area_name(NciScatterArea area) {
    switch (area) {
    case NciScatterArea::None: return "No area";
    case NciScatterArea::Filled: return "Filled";
    case NciScatterArea::Heatmap: return "Heatmap";
    default: return "";
    }
}

const char* nci_point_colouring_name(NciPointColouring colouring) {
    switch (colouring) {
    case NciPointColouring::Scale: return "Colour scale";
    case NciPointColouring::Interaction: return "Interaction";
    case NciPointColouring::Layer: return "Layer colour";
    case NciPointColouring::Region: return "Region";
    default: return "";
    }
}

ImVec4 nci_type_colour(nci_surface::Family family) {
    using nci_surface::Family;
    switch (family) {
    case Family::HydrogenBond: return ImVec4(0.20f, 0.45f, 0.90f, 1.0f);  // blue
    case Family::PiStacking: return ImVec4(0.88f, 0.20f, 0.20f, 1.0f);    // red
    case Family::Ring: return ImVec4(0.20f, 0.68f, 0.28f, 1.0f);          // green
    case Family::HalogenBond: return ImVec4(0.58f, 0.32f, 0.80f, 1.0f);   // purple
    case Family::OtherContact: return ImVec4(0.96f, 0.55f, 0.10f, 1.0f);  // orange
    case Family::Cage: return ImVec4(0.58f, 0.38f, 0.20f, 1.0f);          // brown
    default: return ImVec4(0.50f, 0.50f, 0.50f, 1.0f);                    // grey
    }
}

ImVec4 nci_layer_colour(int index) {
    static const ImVec4 colours[] = {
        ImVec4(0.00f, 0.80f, 0.85f, 1.0f),  // Full system: cyan
        ImVec4(0.95f, 0.85f, 0.15f, 1.0f),  // Preview: yellow
        ImVec4(0.92f, 0.30f, 0.78f, 1.0f),  // magenta
        ImVec4(0.60f, 0.90f, 0.20f, 1.0f),  // lime
        ImVec4(1.00f, 0.65f, 0.72f, 1.0f),  // pink
        ImVec4(0.10f, 0.55f, 0.55f, 1.0f),  // teal
        ImVec4(0.62f, 0.62f, 0.10f, 1.0f),  // olive
        ImVec4(0.92f, 0.92f, 0.92f, 1.0f),  // white
    };
    constexpr int count = (int)(sizeof(colours) / sizeof(colours[0]));
    if (index < 2) return colours[MAX(index, 0)];
    return colours[2 + (index - 2) % (count - 2)];
}

void nci_update_layer_outlines(NciState& state) {
    std::vector<int> keys = {NCI_LAYER_FULL_SYSTEM, NCI_LAYER_PREVIEW};
    keys.insert(keys.end(), state.layer_order.begin(), state.layer_order.end());
    for (auto it = state.layer_outlines.begin(); it != state.layer_outlines.end();)
        it = std::find(keys.begin(), keys.end(), it->first) == keys.end() ? state.layer_outlines.erase(it) : std::next(it);
    for (int layer : keys) {
        const std::vector<NciSurface*> members = nci_layer_members(state, layer);
        uint64_t key = 0;
        size_t sample_count = 0;
        float rho_max = 0.0f;
        for (const NciSurface* n : members) {
            key = md_hash64(&n->id, sizeof(n->id), key);
            key = md_hash64(&n->shown, sizeof(n->shown), key);
            key = md_hash64(&n->want_key, sizeof(n->want_key), key);
            key = md_hash64(&n->have_key, sizeof(n->have_key), key);
            key = md_hash64(&n->rho_max, sizeof(n->rho_max), key);
            if (outline_member(*n)) {
                sample_count += n->scatter.s.size();
                rho_max = MAX(rho_max, n->rho_max);
            }
        }
        key |= 1;
        NciState::LayerOutline& outline = state.layer_outlines[layer];
        if (key == outline.key) continue;
        outline.key = key;
        outline.segments.clear();
        if (!sample_count) continue;

        // Around the interactions' own regions with a one-bin margin, so it never runs on top of their thin outlines
        std::vector<const nci_surface::Scatter*> scatters;
        for (const NciSurface* n : members)
            if (outline_member(*n)) scatters.push_back(&n->scatter);
        if (!nci_surface::scatter_union_outline(&outline.segments, scatters, rho_max, 1))
            MD_LOG_ERROR("NCI layer: invalid combined scatter samples for density outline");
    }
}

void nci_update_surfaces(NciState& state, const NciSurfaceRuntimeInput& runtime,
                         const NciPreviewInput& preview, bool enabled, bool show_context) {
    state.enabled = enabled;
    if (!enabled) {
        const auto stop = [](NciSurface& n) {
            nci_surface_job_destroy(n.job);
            n.job = nullptr;
            n.want_key = 0;
        };
        stop(state.primary);
        stop(state.context);
        for (auto& n : state.previews) stop(*n);
        for (auto& n : state.full_system) stop(*n);
        for (auto& n : state.objects) stop(*n);
        return;
    }

    nci_sync_preview(state, preview);
    nci_apply_settings(state);
    update_nci_surface(state.primary, runtime);
    for (auto& n : state.previews) update_nci_surface(*n, runtime);
    for (auto& n : state.full_system) update_nci_surface(*n, runtime);
    for (auto& n : state.objects) update_nci_surface(*n, runtime);
    // While a setting is dragged, the context cloud and the combined outline keep their last result.
    if (state.settings_editing) {
        nci_update_scale(state);
        return;
    }
    nci_update_layer_outlines(state);

    NciSurface& context = state.context;
    context.shown = show_context;
    context.id = -1;
    context.scope = NciSurface::Scope_Everywhere;
    context.connected = false;
    context.resolution = state.resolution;
    context.rho_max = state.primary.rho_max;
    for (const auto& n : state.previews) if (nci_scale_member(*n)) context.rho_max = MAX(context.rho_max, n->rho_max);
    for (const auto& n : state.full_system) if (nci_scale_member(*n)) context.rho_max = MAX(context.rho_max, n->rho_max);
    for (const auto& n : state.objects) if (nci_scale_member(*n)) context.rho_max = MAX(context.rho_max, n->rho_max);
    context.generation = state.primary.generation;
    update_nci_surface(context, runtime, true);
    nci_update_scale(state);
}

void nci_end_panel_frame(NciState& state) {
    state.settings_editing = false;
    state.primary.editing = false;
    for (auto& n : state.previews) n->editing = false;
    for (auto& n : state.full_system) n->editing = false;
    for (auto& n : state.objects) n->editing = false;
}

}  // namespace topology
