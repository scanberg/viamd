#pragma once

#include "nci_surface_core.h"
#include "nci_surface_job.h"

#include <color_scale.h>
#include <gfx/isomesh_utils.h>
#include <imgui.h>
#include <implot.h>

#include <memory>
#include <unordered_map>
#include <string>
#include <vector>

namespace topology {

// Uniform: every CP type uses the same isovalue and density cutoff (NciState::uniform).
// CPAdaptive: each CP type uses its own isovalue and density cutoff (NciState::family_settings).
// In both modes the extent around the CPs is the radius of each CP type, and one colour scale applies.
enum class NciParameterMode : uint8_t { Uniform, CPAdaptive };

constexpr int NCI_LAYER_FULL_SYSTEM = -1000;
constexpr int NCI_LAYER_PREVIEW = -1001;
constexpr int NCI_LAYER_NONE = 0;

struct NciSurface {
    bool shown = false;
    bool layer_visible = true;
    float s_iso = 0.5f;
    float rho_max = 0.05f;
    ColorScale scale = {.colormap = ImPlotColormap_Jet, .range_beg = -0.04f, .range_end = 0.04f,
                        .symmetric = true, .auto_range = true, .show_legend = true};
    ColorScaleSpan span = {};
    // The CP types of the anchors, as shown, and how many of each; one surface holds the CPs that share parameters
    nci_surface::Families families = 0;
    std::array<uint32_t, nci_surface::family_count> type_counts = {};
    // Internal: the type whose settings the surface uses (see nci_surface::parameter_family_of)
    nci_surface::Family group = nci_surface::Family::Unassigned;
    std::vector<std::vector<std::string>> environments;
    ImVec4 chip = ImVec4(0.25f, 0.65f, 0.95f, 1.0f);
    nci_surface::Scatter scatter;
    bool connected = true;
    bool radius_cap = true;
    bool saved = false;
    bool family_preview = false;
    bool full_system_family = false;
    int object_group = 0;
    int id = 0;
    char name[64] = "Full system";
    std::vector<double> anchors;
    nci_surface::Mesh source_mesh;
    std::vector<std::array<float, 4>> outline;
    float opacity = 1.0f;
    enum Resolution : int { Res_Coarse, Res_Medium, Res_Fine, Res_Count };
    int resolution = Res_Medium;
    enum Scope : int { Scope_Everywhere, Scope_Selected };
    int scope = Scope_Selected;
    float near_radius = 2.7f;
    bool editing = false;

    std::shared_ptr<const std::vector<float>> rho;
    nci_surface::Grid grid = {};
    uint64_t rho_key = 0;
    double rho_seconds = 0.0;
    bool rho_failed = false;

    NciSurfaceJob* job = nullptr;
    uint64_t want_key = 0;
    uint64_t have_key = 0;
    double seconds = 0.0;
    size_t triangles = 0;
    bool cancelled = false;
    bool failed = false;
    uint64_t generation = 1;

    isomesh::Mesh mesh = {};
    uint32_t colormap = 0;
    int uploaded_colormap = -1;
};

// The settings of one CP type, applied automatically to every layer and tuned in the NCI settings.
// The radius applies in both modes; the isovalue and density cutoff only in CP-adaptive mode.
struct NciFamilySettings {
    float s_iso = 0.5f;
    float rho_max = 0.05f;
    float radius = 2.7f;
};

struct NciUniformSettings {
    float s_iso = 0.5f;
    float rho_max = 0.05f;
};

NciFamilySettings nci_family_default_settings(nci_surface::Family family);
NciUniformSettings nci_uniform_default_settings();
ColorScale nci_default_scale();


// The area of a layer in the scatter plot: none, its outline's inside filled, or its sample density
enum class NciScatterArea : uint8_t { None, Filled, Heatmap, Count };
const char* nci_scatter_area_name(NciScatterArea area);

// What the scatter plot draws of a layer, as the Ramachandran layers: visible is the layer's eye, apart from the layer
// shown in the viewport
struct NciScatterLayer {
    bool visible = true;
    bool points = true;
    bool lines = true;
    NciScatterArea area = NciScatterArea::None;
};

// An interaction in the scatter plot, the same in every layer: its isolines and area in its colour; visible toggles them
// and its expected region. Not its points: they are the layer's, as in the viewport.
struct NciTypeStyle {
    bool visible = true;
    bool lines = true;
    NciScatterArea area = NciScatterArea::None;
};

struct NciState {
    NciSurface primary;
    std::array<NciFamilySettings, nci_surface::family_count> family_settings = [] {
        std::array<NciFamilySettings, nci_surface::family_count> settings;
        for (size_t f = 0; f < settings.size(); ++f) settings[f] = nci_family_default_settings((nci_surface::Family)f);
        return settings;
    }();
    NciUniformSettings uniform = nci_uniform_default_settings();
    NciParameterMode mode = NciParameterMode::Uniform;
    // The one colour scale of every NCI surface; Auto follows the shown layer
    ColorScale scale = nci_default_scale();
    ColorScaleSpan span = {};
    int resolution = NciSurface::Res_Medium;
    bool connected = true;
    bool radius_cap = true;
    float opacity = 1.0f;
    bool settings_editing = false;

    NciSurface context;
    std::vector<std::unique_ptr<NciSurface>> objects;
    std::vector<std::unique_ptr<NciSurface>> previews;
    std::vector<std::unique_ptr<NciSurface>> full_system;
    // Per layer key: the outline of all the layer's samples together, and the key of the samples it was made from
    struct LayerOutline {
        uint64_t key = 0;
        std::vector<std::array<float, 4>> segments;
    };
    std::unordered_map<int, LayerOutline> layer_outlines;
    nci_surface::View view = nci_surface::View::Overview;
    bool classified_overview = true;
    uint64_t selection_key = 0;
    uint64_t full_selection_key = 0;
    bool force_all_cps = false;
    int next_group_id = 1;
    int next_id = 1;
    int active_id = 0;
    // The one layer shown at a time; NCI_LAYER_NONE hides every layer.
    int selected_layer = NCI_LAYER_FULL_SYSTEM;
    bool layer_hidden = false;
    int pending_delete = 0;
    std::vector<int> layer_order;   // saved layers only
    bool show_guides = true;
    // The expected region of each interaction in the scatter plot (nci_surface::expected_region)
    bool show_reference = true;
    // The context cloud: every grid sample of the whole density, under the layers
    bool show_context_samples = true;
    // Per layer key, what the scatter plot draws of it; a layer without an entry has the defaults
    std::unordered_map<int, NciScatterLayer> scatter_layers;
    // Interactions toggled in the legend: their outlines and expected regions
    std::array<NciTypeStyle, nci_surface::family_count> type_styles = {};
    // Scatter plot look, in pixels: a layer's outline (thick), its interactions' outlines (thin), the samples
    static constexpr float default_layer_outline_width = 4.0f, default_type_outline_width = 1.5f, default_point_size = 1.8f;
    float layer_outline_width = default_layer_outline_width;
    float type_outline_width = default_type_outline_width;
    float point_size = default_point_size;
    bool capacity_reported = false;
    bool enabled = false;
};

struct NciSurfaceRuntimeInput {
    nci_surface::DensityProvider density_provider = {};
    std::array<double, 3> bounds_min = {};
    std::array<double, 3> bounds_max = {};
};

struct NciCriticalPoint {
    std::array<double, 3> position = {};
    nci_surface::Family family = nci_surface::Family::Unassigned;
    nci_surface::Family parameter_family = nci_surface::Family::Unassigned;
    bool selected = false;
    std::vector<std::string> environment;
};

struct NciPreviewInput {
    uint64_t selection_key = 0;
    bool enabled = false;
    bool has_selection = false;
    std::vector<NciCriticalPoint> critical_points;
};

int nci_layer_key(const NciSurface& surface);
ColorScale nci_effective_scale(const NciState& state, const NciSurface& surface);
void nci_select_layer(NciState& state, int key);
void nci_select(NciState& state, int id);
bool nci_scale_member(const NciSurface& surface);
std::vector<NciSurface*> nci_layer_members(NciState& state, int key);
// The CP types of the members with their CP counts, e.g. "Hydrogen bonds 2, Rings 1".
std::string nci_layer_types(const std::vector<NciSurface*>& members);
// The interaction a surface is computed for, e.g. "\xCF\x80-stacking" (its rings closed by \xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80 included)
std::string nci_family_label(const NciSurface& surface);
// The index of the settings the surface uses
size_t nci_settings_index(const NciSurface& surface);
void nci_apply_settings(NciState& state);
void nci_update_scale(NciState& state);
bool nci_surface_visible(const NciState& state, const NciSurface& surface);
void nci_surface_release(NciSurface& surface);
void nci_state_free(NciState& state);
void nci_sync_preview(NciState& state, const NciPreviewInput& input);
void nci_update_layer_outlines(NciState& state);
// Scatter colours. An interaction always has the same colour; layers get colours apart from those.
// Layer index 0 is the Full system, 1 the Preview, 2 and up the saved layers.
ImVec4 nci_type_colour(nci_surface::Family family);
ImVec4 nci_layer_colour(int index);
void nci_update_surfaces(NciState& state, const NciSurfaceRuntimeInput& runtime,
                         const NciPreviewInput& preview, bool enabled, bool show_context);
void nci_end_panel_frame(NciState& state);
bool upload_nci_surface_mesh(NciSurface& surface);
void update_nci_surface(NciSurface& surface, const NciSurfaceRuntimeInput& input, bool context = false);

}  // namespace topology
