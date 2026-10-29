// RAMACHANDRAN
//
// The backbone dihedrals (phi, psi) of the protein residues, plotted in the four classes MolProbity validates them in
// (general, glycine, proline, pre-proline), or all in one plot.
//
// What is plotted is a stack of LAYERS the user defines. A layer is either the reference distribution (the favoured and
// allowed regions of high quality structures), or a subset of the system's residues (a filter, or the selection as it
// changes) over some of its frames: the current one, all of them, or those of a timeline period. A layer of the system
// is drawn as points, as contours of its density (at levels of the fraction of the samples they enclose), or as a heat
// map of it. The layers are drawn bottom up, the first one on top.
//
// The density of a layer is estimated off the main thread (density_map.h): the angles of its residues over its frames are
// binned on a periodic grid, smoothed with a gaussian of the layer's bandwidth, and normalized per class. What its
// contours are drawn at is found from the cumulative distribution of the density (the highest density regions): a level
// of 90% is the line that encloses 90% of the samples, whatever the number of frames or the bandwidth. The jobs are keyed
// by everything they depend on, so a layer is recomputed exactly when its result would differ.
//
// The plots of the layers are composited into one texture per plot, drawn again only when something they show
// changes.
//
//   shift-drag             selects the residues of the points in the box (right button: deselects them)
//   shift-click            on a point, selects its residue (right button: deselects it); right click elsewhere clears
//   hover                  highlights the residue of the point under the pointer, and says what it is
//   double-click           fits the view
//   a period dropped       from the timeline on a plot, a new layer of the density over it; on a layer, its frames

#include "density_gen.inl"
#include "density_gly.inl"
#include "density_pro.inl"
#include "density_pre.inl"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <md_system.h>
#include <md_util.h>

#include <viamd_event.h>
#include <viamd.h>
#include <serialization_utils.h>
#include <timeline_periods_ui.h>
#include <filter_input.h>
#include <property_grid.h>
#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <task_system.h>
#include <density_map.h>
#include <layer_groups.h>

#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <imgui_internal.h>
#include <implot_internal.h>

#include <atomic>
#include <memory>
#include <math.h>
#include <string.h>
#include <stdlib.h>

namespace ramachandran {

// ## Constants

// The classes of residues, in the order of md_ramachandran_type_t (less its UNKNOWN), and the plot of all of them
enum : int {
    Type_General = 0,
    Type_Glycine,
    Type_Proline,
    Type_PreProline,
    NUM_TYPES,
};

constexpr int PANEL_ALL  = NUM_TYPES;       // Every residue in one plot
constexpr int NUM_PANELS = NUM_TYPES + 1;

static const char* panel_names[NUM_PANELS]       = { "General", "Glycine", "Proline", "Pre-proline", "All residues" };

constexpr int MAX_LAYERS        = 16;
constexpr int MAX_LEVELS        = DENSITY_MAP_MAX_LEVELS;
constexpr size_t MAX_POINTS     = 100000;   // Of a layer of points over frames, in one plot: frames are skipped beyond

// The reference: densities normalized to a maximum of 1 on grids of 2 degree bins centred on -179, -177, ..., 179,
// indexed [phi][psi]; and MolProbity's levels of them: favoured encloses 98% of its residues, allowed 99.95% (general)
// or 99.8% (the others)
typedef double density_map_t[180][180];
static const density_map_t* ref_maps[NUM_TYPES] = { &density_gen, &density_gly, &density_pro, &density_pre };
static const float ref_favored[NUM_TYPES] = { 0.02f,   0.02f,  0.02f,  0.02f  };
static const float ref_allowed[NUM_TYPES] = { 0.0005f, 0.002f, 0.002f, 0.002f };

// The fills of the reference, allowed and favoured, by class (the plot of all residues takes the general's)
static const vec4_t ref_fill_allowed[NUM_TYPES] = {
    {0.70f, 0.91f, 1.00f, 1.0f},
    {1.00f, 0.91f, 0.77f, 1.0f},
    {0.82f, 1.00f, 0.77f, 1.0f},
    {0.70f, 0.91f, 1.00f, 1.0f},
};
static const vec4_t ref_fill_favored[NUM_TYPES] = {
    {0.50f, 0.85f, 1.00f, 1.0f},
    {1.00f, 0.80f, 0.50f, 1.0f},
    {0.50f, 1.00f, 0.55f, 1.0f},
    {0.50f, 0.85f, 1.00f, 1.0f},
};

enum Region : int {
    Region_Favored = 0,
    Region_Allowed,
    Region_Outlier,
    Region_Count,
};

static const char* region_names[Region_Count] = { "Favored", "Allowed", "Outlier" };
static const vec4_t region_colors[Region_Count] = {
    {0.95f, 0.95f, 0.95f, 1.0f},   // Favoured: the layer's own colour, this only where a colour is asked for by itself
    {1.00f, 0.76f, 0.20f, 1.0f},
    {1.00f, 0.27f, 0.27f, 1.0f},
};

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

// ## The reference

// The reference density of a class at (phi, psi) in degrees, interpolated
static double ref_density(int type, double phi, double psi) {
    const density_map_t& m = *ref_maps[type];
    const double x = (phi + 180.0) * 0.5 - 0.5;
    const double y = (psi + 180.0) * 0.5 - 0.5;
    const double fx = floor(x);
    const double fy = floor(y);
    const int x0 = (((int)fx % 180) + 180) % 180;
    const int y0 = (((int)fy % 180) + 180) % 180;
    const int x1 = (x0 + 1) % 180;
    const int y1 = (y0 + 1) % 180;
    const double tx = x - fx;
    const double ty = y - fy;
    return lerp(lerp(m[x0][y0], m[x1][y0], tx), lerp(m[x0][y1], m[x1][y1], tx), ty);
}

static inline Region ref_region(int type, float phi_rad, float psi_rad) {
    const double d = ref_density(type, RAD_TO_DEG(phi_rad), RAD_TO_DEG(psi_rad));
    return d >= ref_favored[type] ? Region_Favored : d >= ref_allowed[type] ? Region_Allowed : Region_Outlier;
}

// Angles of (0, 0) are those a segment has none of: the ends of a chain, frames not computed yet
static inline bool has_angles(md_backbone_angles_t a) {
    return !(a.phi == 0.0f && a.psi == 0.0f);
}


// The plane of the angles, in degrees: periodic
static const DensityMapDomain rama_domain = { -180.0f, 180.0f, -180.0f, 180.0f, true };

// ## Layers

enum LayerSource : int {
    LayerSource_Reference = 0,  // The reference distribution
    LayerSource_System,         // Residues of the system
    LayerSource_Count,
};

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
    PointColoring_Region,       // Allowed and outliers stand out
    PointColoring_Time,         // The frame, through the layer's colour map
    PointColoring_Count,
};

static const char* source_names[LayerSource_Count] = { "Reference", "System" };
static const char* style_names[LayerStyle_Count] = { "Points", "Contours", "Heatmap" };
static const char* point_coloring_names[PointColoring_Count] = { "Layer color", "Region", "Time" };

// What a density job found
// A density job's: the density's distribution by class (a channel each, and all of them in PANEL_ALL), and where its
// samples fall in the reference
struct DensityResult {
    DensityMapResult map = {};
    uint64_t regions[NUM_PANELS][Region_Count] = {};
    uint32_t num_frames = 0;
};

struct Layer {
    uint32_t id = 0;            // Unique while it lives, never 0
    RepresentationNode node = {};   // Where it is in the list: the group it is in (layer_groups.h). node.id is id.
    char     name[64] = "";     // Empty: named after what it shows (layer_label)
    bool     visible = true;

    LayerSource source = LayerSource_System;
    LayerFrames frames = LayerFrames_Current;
    uint32_t    period = 0;     // Of LayerFrames_Period, by id
    AtomFilter  filter = {};    // The residues: those whose CA it picks, within its group's. Empty: all of them.
    bool        follow_selection = false;  // The residues of the selection instead of the filter

    LayerStyle style = LayerStyle_Points;
    vec4_t   color = {1, 1, 1, 1};
    float    opacity = 1.0f;

    // Contours
    bool     fill = false;
    bool     lines = true;
    float    line_width = 1.5f;            // Pixels
    int      num_levels = 3;
    float    levels[MAX_LEVELS] = { 0.5f, 0.8f, 0.95f };   // Fractions of the samples enclosed, increasing
    float    bandwidth = 5.0f;             // Degrees: the sigma of the gaussian the density is smoothed with

    // Heatmap (and points by time)
    ImPlotColormap colormap = ImPlotColormap_Viridis;
    bool     log_scale = false;

    // Points
    PointColoring point_coloring = PointColoring_Layer;
    float    point_size = 4.0f;            // Pixels across

    // ### Derived, not saved

    // The residues, as segments of the backbone by class, and what they were found from
    md_bitset_t atoms = {};
    md_array(uint32_t) segs[NUM_TYPES] = {};
    uint64_t segs_key = 0;          // Of what they were found from
    uint64_t segs_hash = 0;         // Of the segments themselves
    bool     segs_valid = false;    // The filter evaluated (or there is none)

    // The density, on the GPU, and what it is of
    uint32_t tex = 0;
    DensityResult result = {};
    bool     has_result = false;
    uint64_t shown_hash = 0;        // Of the inputs of the density in tex
    uint64_t want_hash = 0;         // Of the inputs it should be of, 0 for none (as for points)
    DensityMapJob* job = nullptr;   // Under way
    uint64_t job_serial = 0;        // Of the job launched last: the callbacks of cancelled ones are not its

    // Of the points shown: how many fall in each region, by plot
    uint32_t point_regions[NUM_PANELS][Region_Count] = {};
    uint32_t point_count[NUM_PANELS] = {};
    uint32_t point_stride = 1;      // Frames skipped to keep the points within MAX_POINTS
};

// ## Masks

// Adds the atoms a backbone segment stands for to mask: those of its residue, or for a segment without a residue (a
// backbone found from the bonds alone, as for an XYZ file) those of its backbone
static void mask_add_segment(md_bitmask_t* mask, const md_system_t& sys, size_t seg_idx) {
    const md_residue_idx_t res_idx = sys.protein_backbone.segment.res_idx[seg_idx];
    if (res_idx >= 0) {
        const md_urange_t range = md_system_residue_atom_range(&sys, res_idx);
        if (range.beg != range.end) {
            modify_field(mask, range, SelectionOperator::Or);
        }
        return;
    }
    const md_amino_acid_atoms_t& a = sys.protein_backbone.segment.atoms[seg_idx];
    const md_atom_idx_t atoms[] = { a.n, a.ca, a.c, a.o };
    for (md_atom_idx_t i : atoms) {
        if (i >= 0 && (size_t)i < sys.atom.count) md_bitmask_set_bit(mask, (uint32_t)i);
    }
}

// The atom that tells whether a segment is in a set of atoms: its CA (else its N)
static inline md_atom_idx_t segment_atom(const md_system_t& sys, size_t seg_idx) {
    const md_amino_acid_atoms_t& a = sys.protein_backbone.segment.atoms[seg_idx];
    return a.ca >= 0 ? a.ca : a.n;
}

// ## The component

struct Ramachandran : viamd::EventHandler {
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
    static constexpr uint32_t LIST_ID = 1;  // What its rows are dragged as (LayerNodeRef)

    // The class of each segment of the backbone, or -1 where it has no angles to plot (the ends of chains, traces)
    md_array(int8_t) seg_type = 0;
    uint64_t seg_type_key = 0;      // Of the topology it was found for
    uint32_t seg_count[NUM_PANELS] = {};
    uint64_t system_generation = 1;  // Moves with every system loaded or freed: a density is of one of them

    // How it is shown
    int  view = -1;                 // The plot shown alone (a class or PANEL_ALL), -1 for the four classes side by side
    bool show_sidebar = true;
    bool show_ticks = true;
    ImPlotRect viewrect = ImPlotRect(-180, 180, -180, 180);
    // Fitting the view: the whole domain, centred, extended along the longer side of the plots. Kept fitted as the
    // plots change shape until it is panned or zoomed.
    bool  fit_view = true;

    // The highlight is whoever sets it to clear: what this window highlighted and no longer does is cleared once
    bool  owns_highlight = false;
    bool  highlights_now = false;
    bool  view_fitted = false;
    float plot_aspect = 0.0f;       // Width over height of the plots, the last frame

    // GL: the reference's densities, and the densities of the layers composited by plot
    uint32_t ref_tex = 0;
    DensityMapCanvas canvas[NUM_PANELS] = {};

    // Box selection, in the plot it started in
    struct {
        bool active = false;
        int  panel = -1;
        SelectionOperator op = SelectionOperator::Or;
        ImPlotPoint p0, p1;
        int64_t click_seg = -1;     // The point pressed on, which a click (a box of no size) takes
    } box;

    // Of a workspace from before there were layers: its settings, turned into layers once it is read
    struct {
        bool  present = false;
        float alpha[3] = { 0.85f, 0.85f, 0.85f };
        int   mode[3] = { 0, 1, 1 };    // IsoLevels, IsoLines, Colormap
        int   colormap[3] = { ImPlotColormap_Hot, ImPlotColormap_Plasma, ImPlotColormap_Viridis };
        vec4_t isoline_color[3] = { {1,1,1,1}, {1,1,1,1}, {1,1,1,1} };
        int   show[4] = { 1, 1, 1, 1 };
        int   layout = 0;
        int   period = 0;
        float blur_sigma = 0.0f;
    } legacy;
    bool read_layers = false;       // While a workspace is read: whether it had any
    // and the groups read, by the ids the file has for them
    int      num_read_groups = 0;
    int      read_group_file_id[LAYER_MAX_GROUPS] = {};
    uint32_t read_group_id[LAYER_MAX_GROUPS] = {};
    bool layers_declared = false;   // and whether it said how many (possibly none)

    Ramachandran() { viamd::event_system_register_handler(*this); }

    // ## Events

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app = (ApplicationState*)e.payload;
                initialize();
                workspace_register_window("Ramachandran", &show_window);
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
                legacy = {};
                read_layers = false;
                layers_declared = false;
                num_read_groups = 0;
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                const str_t section = viamd::section_header(state);
                if (str_eq(section, STR_LIT("Ramachandran"))) {
                    deserialize_settings(state);
                } else if (str_eq(section, STR_LIT("RamachandranLayer"))) {
                    deserialize_layer(state);
                } else if (str_eq(section, STR_LIT("RamachandranGroup"))) {
                    deserialize_group(state);
                }
                break;
            }
            case viamd::EventType_ViamdDeserializeEnd:
                if (!read_layers && !layers_declared) {
                    if (legacy.present) {
                        create_legacy_layers();
                    } else {
                        create_default_layers();
                    }
                }
                break;
            case viamd::EventType_ViamdSystemInit:
                // Classified again, whatever its topology's version says
                seg_type_key = 0;
                system_generation += 1;
                break;
            case viamd::EventType_ViamdSystemFree:
            case viamd::EventType_ViamdTrajectoryFree:
                // What the densities were of is gone: nothing of it is shown again, even where the next system's
                // attributes happen to have the same versions
                cancel_all_jobs(true);
                for (int i = 0; i < num_layers; ++i) {
                    layers[i].has_result = false;
                    layers[i].shown_hash = 0;
                    layers[i].want_hash = 0;
                }
                system_generation += 1;
                panel_hash_reset();
                break;
            case viamd::EventType_ViamdFrameTick:
                update();
                draw();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Ramachandran", &show_window);
                break;
            default:
                break;
            }
        }
    }

    void initialize() {
        ref_tex = density_map_texture_create(true);
        init_reference(ref_tex);

        create_default_layers();
    }

    void shutdown() {
        cancel_all_jobs(true);
        clear_layers();
        md_array_free(seg_type, md_get_heap_allocator());
        density_map_texture_free(&ref_tex);
        for (int i = 0; i < NUM_PANELS; ++i) density_map_canvas_free(&canvas[i]);
        density_map_shutdown();
    }

    // The reference densities, resampled onto the grid of the densities (texel centres) and smoothed a little: they are
    // tabulated at 2 degrees
    static void init_reference(uint32_t tex) {
        STATIC_ASSERT(NUM_TYPES == DENSITY_MAP_CHANNELS, "a channel for each class");
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        const int dim = DENSITY_MAP_DIM;
        float* data = (float*)md_temp_alloc_aligned(temp, sizeof(float) * 4 * dim * dim, 16);
        for (int y = 0; y < dim; ++y) {
            const double psi = (y + 0.5) * 360.0 / dim - 180.0;
            for (int x = 0; x < dim; ++x) {
                const double phi = (x + 0.5) * 360.0 / dim - 180.0;
                float* v = data + 4 * (y * dim + x);
                for (int t = 0; t < NUM_TYPES; ++t) {
                    v[t] = (float)ref_density(t, phi, psi);
                }
            }
        }
        density_map_blur(data, dim, 1.5f, true);
        density_map_texture_upload(tex, data);
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
                if (layers[i].source == LayerSource_System && vec4_equal(layers[i].color, layer_palette[p])) ++uses;
            }
            if (uses < best_uses) {
                best = p;
                best_uses = uses;
            }
        }
        return layer_palette[best];
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
                const Layer tmp = layers[j];
                layers[j] = layers[j - 1];
                layers[j - 1] = tmp;
                moved = true;
            }
        }
        if (moved) panel_hash_reset();
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

    // Adds a layer: first in the group 'parent' (or at the top) for 'at' 0, last for -1, or before what is at position
    // 'at' of the list. NULL when there is no room. In place in the list once the tree is applied (apply_tree).
    Layer* add_layer(LayerSource source, LayerFrames frames, LayerStyle style, int at = -1, uint32_t parent = 0) {
        if (num_layers >= MAX_LAYERS) return nullptr;
        Layer& l = layers[num_layers++];
        l = Layer{};
        l.id = next_layer_id++;
        l.node.id = l.id;
        l.node.parent = parent;
        if (at < 0) {
            l.node.order = LAYER_TREE_LAST;
        } else if (at == 0) {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            l.node.order = layer_tree_first_order(nodes, tree_nodes(nodes), parent);
        } else {
            l.node.order = 2 * at - 1;  // Between those at at - 1 and at, the orders being two apart
        }
        l.source = source;
        l.frames = frames;
        l.style = style;
        l.color = source == LayerSource_Reference ? vec4_t{1, 1, 1, 1} : next_layer_color();
        if (source == LayerSource_Reference) {
            l.fill = true;
            l.lines = false;
            l.opacity = 0.85f;
        } else if (style == LayerStyle_Contours) {
            l.fill = false;
            l.lines = true;
        }
        if (frames == LayerFrames_Current && style == LayerStyle_Points) {
            l.point_coloring = PointColoring_Region;
        }
        tree_changed();
        return &l;
    }

    void free_layer(Layer& l) {
        md_allocator_i* heap = md_get_heap_allocator();
        cancel_job(l, true);
        if (l.tex) glDeleteTextures(1, &l.tex);
        l.tex = 0;
        md_bitset_free(&l.atoms, heap);
        for (int t = 0; t < NUM_TYPES; ++t) md_array_free(l.segs[t], heap);
    }

    void remove_layer(uint32_t id) {
        const int idx = layer_index(id);
        if (idx < 0) return;
        free_layer(layers[idx]);
        for (int i = idx; i + 1 < num_layers; ++i) layers[i] = layers[i + 1];
        num_layers -= 1;
        layers[num_layers] = Layer{};
        if (selected_layer == id) {
            selected_layer = num_layers > 0 ? layers[MIN(idx, num_layers - 1)].id : 0;
        }
        tree_changed();
        panel_hash_reset();
    }

    // A copy of what is set, without what is derived from it (found again), at 'order' in the group 'parent'
    Layer* copy_layer(const Layer& src, uint32_t parent, int32_t order) {
        if (num_layers >= MAX_LAYERS) return nullptr;
        const Layer copy = src;
        Layer& l = layers[num_layers++];
        l = copy;
        l.id = next_layer_id++;
        l.node = {l.id, parent, order, false};
        l.atoms = {};
        for (int t = 0; t < NUM_TYPES; ++t) l.segs[t] = nullptr;
        l.segs_key = 0;
        l.segs_hash = 0;
        l.segs_valid = false;
        l.tex = 0;
        l.has_result = false;
        l.shown_hash = 0;
        l.want_hash = 0;
        l.job = nullptr;
        l.job_serial = 0;
        tree_changed();
        return &l;
    }

    // Just above it, in its group
    void duplicate_layer(uint32_t id) {
        const Layer* src = find_layer(id);
        if (!src) return;
        if (Layer* l = copy_layer(*src, src->node.parent, src->node.order - 1)) {
            selected_layer = l->id;
            selected_group = 0;
        }
    }

    void clear_layers() {
        for (int i = 0; i < num_layers; ++i) {
            free_layer(layers[i]);
            layers[i] = Layer{};
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
        panel_hash_reset();
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

    // A group of the layer, where it was, with the layer's residues as its subset: the layers added to it show them
    void group_layer(uint32_t id) {
        Layer* l = find_layer(id);
        if (!l || l->node.parent != 0) return;
        LayerGroup* g = add_group(l->node.order - 1);
        if (!g) return;
        MEMCPY(g->filter.text, l->filter.text, sizeof(g->filter.text));
        g->follow_selection = l->follow_selection;
        l->filter.text[0] = '\0';
        l->follow_selection = false;
        l->segs_key = 0;
        l->node.parent = g->node.id;
        l->node.order = LAYER_TREE_LAST;
        tree_changed();
        selected_group = g->node.id;
        selected_layer = 0;
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

    // In a group now: a filter that is the group's own says nothing more (and within it would pick the same again, or
    // less, as residue(2:8) within residue(2:8) is its first seven residues), so it goes
    void entered_group(Layer& l) {
        l.segs_key = 0;
        const LayerGroup* g = layer_group(l);
        if (g && !g->follow_selection && g->filter.text[0] && strcmp(g->filter.text, l.filter.text) == 0) l.filter.text[0] = '\0';
        if (g && g->follow_selection && l.follow_selection) l.follow_selection = false;
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
            const LayerGroup src = *g;
            LayerGroup* copy = add_group(src.node.order + 1);
            if (!copy) return;
            const uint32_t copy_id = copy->node.id;
            const RepresentationNode node = copy->node;
            *copy = LayerGroup{};
            copy->node = node;
            MEMCPY(copy->name, src.name, sizeof(copy->name));
            copy->visible = src.visible;
            copy->open = src.open;
            copy->filter.text[0] = '\0';
            MEMCPY(copy->filter.text, src.filter.text, sizeof(copy->filter.text));
            copy->follow_selection = src.follow_selection;
            // In the order they are in, each last in the copy
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t n = tree_nodes(nodes);
            int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t nk = representation_tree_children(kids, ARRAY_SIZE(kids), nodes, n, id);
            for (size_t k = 0; k < nk; ++k) {
                if (kids[k] < num_layers) copy_layer(layers[kids[k]], copy_id, LAYER_TREE_LAST);
            }
            selected_group = copy_id;
            selected_layer = 0;
            return;
        }
        if (action == LayerGroupAction_Remove) {
            for (int i = num_layers - 1; i >= 0; --i) {
                if (layers[i].node.parent == id) remove_layer(layers[i].id);
            }
        } else if (action == LayerGroupAction_Ungroup) {
            const int32_t at = g->node.order + 1;   // Where the group was, in the order they were in it
            for (int i = 0; i < num_layers; ++i) {
                Layer& l = layers[i];
                if (l.node.parent != id) continue;
                // What can be kept is (layer_group_keeps_filter); what cannot is left as it is, as the menu warned
                if (l.source == LayerSource_System && layer_group_keeps_filter(l.filter.text, l.follow_selection, *g)) {
                    if (g->follow_selection) l.follow_selection = true;
                    else layer_group_compose_filter(l.filter.text, sizeof(l.filter.text), *g);
                }
                l.node.parent = 0;
                l.node.order = at;
                l.segs_key = 0;
            }
        }
        const int gi = (int)(g - groups);
        layer_group_free(g);
        for (int i = gi; i + 1 < num_groups; ++i) groups[i] = groups[i + 1];
        num_groups -= 1;
        groups[num_groups] = LayerGroup{};
        if (selected_group == id) selected_group = 0;
        tree_changed();
    }

    void reset_settings() {
        view = -1;
        show_sidebar = true;
        show_ticks = true;
        viewrect = ImPlotRect(-180, 180, -180, 180);
        fit_view = true;
    }

    // The reference, the current frame as points, and the whole trajectory as contours (shown when there is one)
    //
    // The two of the system in a group, as an example of one: its subset is where to narrow them all at once, and what
    // is added to it shows the same residues. All of them to begin with, which every system with a backbone has (one
    // read from an XYZ file, whose backbone is found from its bonds, has no protein to name).
    void create_default_layers() {
        clear_layers();
        uint32_t group = 0;
        if (LayerGroup* g = add_group(LAYER_TREE_LAST)) {
            snprintf(g->name, sizeof(g->name), "Residues");
            snprintf(g->filter.text, sizeof(g->filter.text), "all");
            group = g->node.id;
        }
        Layer* cur = add_layer(LayerSource_System, LayerFrames_Current, LayerStyle_Points, -1, group);
        Layer* all = add_layer(LayerSource_System, LayerFrames_All, LayerStyle_Contours, -1, group);
        add_layer(LayerSource_Reference, LayerFrames_All, LayerStyle_Contours);
        if (cur) {
            cur->color = {1, 1, 1, 1};
            selected_layer = cur->id;
        }
        if (all) {
            all->color = {1, 1, 1, 1};     // Lines that stand out on the reference's regions and the points
        }
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

    // The density of the angles of the layer's residues over frames [frame_beg, frame_end): a channel for each class.
    // Where the samples fall in the reference is counted along, by chunk.
    void launch_job(Layer& l, uint64_t hash, uint32_t frame_beg, uint32_t frame_end) {
        ASSERT(!l.job);
        const md_backbone_angles_t* angles = app->trajectory_data.backbone_angles.data;
        ASSERT(angles && frame_end > frame_beg);

        struct Input {
            const md_backbone_angles_t* angles = nullptr;
            size_t   stride = 0;
            uint32_t frame_beg = 0, frame_end = 0;
            md_array(uint32_t) segs[NUM_TYPES] = {};
            uint64_t regions[8][NUM_TYPES][Region_Count] = {};     // By chunk
            ~Input() { for (int t = 0; t < NUM_TYPES; ++t) md_array_free(segs[t], md_get_heap_allocator()); }
        };
        std::shared_ptr<Input> in = std::make_shared<Input>();
        in->angles = angles;
        in->stride = app->trajectory_data.backbone_angles.stride;
        in->frame_beg = frame_beg;
        in->frame_end = frame_end;
        uint64_t num_samples = 0;
        for (int t = 0; t < NUM_TYPES; ++t) {
            md_array_push_array(in->segs[t], l.segs[t], md_array_size(l.segs[t]), md_get_heap_allocator());
            num_samples += md_array_size(l.segs[t]);
        }
        num_samples *= (frame_end - frame_beg);
        // A chunk for every 2M samples, no more than there are frames
        const uint32_t num_chunks = (uint32_t)MIN((uint64_t)MIN(8u, frame_end - frame_beg), (num_samples >> 21) + 1);

        auto source = [in](uint32_t chunk, uint32_t chunks, DensityMapBins& bins, const std::atomic<int32_t>& cancel) {
            const uint32_t n = in->frame_end - in->frame_beg;
            const uint32_t beg = in->frame_beg + (uint32_t)((uint64_t)n * chunk / chunks);
            const uint32_t end = in->frame_beg + (uint32_t)((uint64_t)n * (chunk + 1) / chunks);
            for (uint32_t f = beg; f < end; ++f) {
                if (cancel.load(std::memory_order_relaxed)) return;
                const md_backbone_angles_t* a = in->angles + (size_t)f * in->stride;
                for (int t = 0; t < NUM_TYPES; ++t) {
                    const uint32_t* segs = in->segs[t];
                    const size_t num = md_array_size(segs);
                    for (size_t i = 0; i < num; ++i) {
                        const md_backbone_angles_t ang = a[segs[i]];
                        if (!has_angles(ang)) continue;
                        bins.add(t, RAD_TO_DEG(ang.phi), RAD_TO_DEG(ang.psi));
                        in->regions[chunk][t][ref_region(t, ang.phi, ang.psi)] += 1;
                    }
                }
            }
        };
        const uint32_t id = l.id;
        const uint64_t serial = ++l.job_serial;
        auto done = [this, in, id, hash, serial](const DensityMapResult* result, const float* grid) {
            Layer* layer = find_layer(id);
            if (!layer || !layer->job || layer->job_serial != serial || layer->want_hash != hash) return;
            layer->job = nullptr;
            if (!result) {
                // Interrupted: tried again
                layer->want_hash = 0;
                return;
            }
            if (!layer->tex) layer->tex = density_map_texture_create(true);
            density_map_texture_upload(layer->tex, grid);
            DensityResult& r = layer->result;
            r = {};
            r.map = *result;
            r.num_frames = in->frame_end - in->frame_beg;
            for (int c = 0; c < 8; ++c) {
                for (int t = 0; t < NUM_TYPES; ++t) {
                    for (int k = 0; k < Region_Count; ++k) {
                        r.regions[t][k] += in->regions[c][t][k];
                        r.regions[PANEL_ALL][k] += in->regions[c][t][k];
                    }
                }
            }
            layer->has_result = true;
            layer->shown_hash = hash;
        };
        l.want_hash = hash;
        l.job = density_map_launch(STR_LIT("Ramachandran density"), rama_domain, l.bandwidth, num_chunks, source, done);
    }

    // ## Update

    // The class of each segment: rebuilt when the topology changes
    void classify_segments() {
        const md_system_t& sys = app->mold.sys;
        const size_t num_segments = sys.protein_backbone.segment.count;
        const uint64_t key = md_hash64_combine(sys.topology_version, num_segments);
        if (key == seg_type_key && md_array_size(seg_type) == num_segments) return;

        // The jobs under way index the old segments
        cancel_all_jobs(true);

        md_allocator_i* heap = md_get_heap_allocator();
        md_array_resize(seg_type, num_segments, heap);
        MEMSET(seg_count, 0, sizeof(seg_count));
        if (num_segments > 0) {
            MEMSET(seg_type, -1, num_segments);
            md_temp_scope_t temp = md_temp_begin();
            md_ramachandran_type_t* types = md_temp_alloc_array(temp, md_ramachandran_type_t, num_segments);
            md_util_backbone_ramachandran_classify(types, num_segments, &sys);
            // The first and the last segment of a backbone have no phi and no psi (md_util_backbone_angles_compute)
            const md_protein_backbone_data_t& bb = sys.protein_backbone;
            for (size_t r = 0; r < bb.range.count; ++r) {
                const uint32_t beg = bb.range.offset[r];
                const uint32_t end = bb.range.offset[r + 1];
                if (end - beg < 4) continue;
                for (uint32_t i = beg + 1; i + 1 < end && i < num_segments; ++i) {
                    if (types[i] == MD_RAMACHANDRAN_TYPE_UNKNOWN) continue;
                    const int t = (int)types[i] - (int)MD_RAMACHANDRAN_TYPE_GENERAL;
                    if (t < 0 || t >= NUM_TYPES) continue;
                    seg_type[i] = (int8_t)t;
                    seg_count[t] += 1;
                    seg_count[PANEL_ALL] += 1;
                }
            }
            md_temp_end(temp);
        }
        seg_type_key = key;
        // What the layers had indexes the segments of before
        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            l.segs_key = 0;
            l.segs_valid = false;
            for (int t = 0; t < NUM_TYPES; ++t) md_array_shrink(l.segs[t], 0);
        }
    }

    // The residues of the layer, as segments by class: found again when what they are found from changes
    void update_layer_segments(Layer& l) {
        if (l.source != LayerSource_System) return;
        const md_system_t& sys = app->mold.sys;
        md_allocator_i* heap = md_get_heap_allocator();

        // What decides them: its own, within its group's
        const LayerGroup* g = layer_group(l);
        const bool in_context = g && layer_group_has_context(*g);
        uint64_t key = md_hash64_combine(seg_type_key, l.follow_selection ? 1 : 2);
        if (in_context) key = md_hash64_combine(key, g->hash);
        bool evaluate = false;
        if (l.follow_selection) {
            const md_bitset_t sel = md_bitmask_view(&app->selection.selection_mask);
            key = md_hash64_combine(key, md_bitset_hash64(sel, 0));
        } else if (l.filter.text[0] != '\0') {
            evaluate = true;
            key = md_hash64_combine(key, md_hash64(l.filter.text, strlen(l.filter.text), app->script.ir_fingerprint));
            // What depends on the coordinates follows the frame as points of it, but is taken once for frames that
            // are many: the density over them would be found again at every frame
            if (l.filter.result.is_dynamic && l.frames == LayerFrames_Current) {
                key = md_hash64_combine(key, (uint64_t)(app->animation.frame * 4096.0));
            }
        }
        if (key == l.segs_key) return;
        l.segs_key = key;

        md_bitset_t atoms = {};
        bool all = false;
        l.segs_valid = true;
        if (in_context && !g->valid) {
            // Nothing within a group whose filter does not evaluate
            l.segs_valid = false;
        } else if (l.follow_selection) {
            atoms = md_bitmask_view(&app->selection.selection_mask);
            if (in_context) {
                md_bitset_free(&l.atoms, heap);
                l.atoms = md_bitset_and(atoms, g->atoms, heap);
                atoms = l.atoms;
            }
        } else if (evaluate) {
            l.segs_valid = filter_evaluate(&l.filter, &l.atoms, layer_group_filter_context(app, in_context ? g : nullptr), heap);
            atoms = l.atoms;
        } else if (in_context) {
            atoms = g->atoms;
        } else {
            all = true;
        }

        for (int t = 0; t < NUM_TYPES; ++t) md_array_shrink(l.segs[t], 0);
        if (l.segs_valid) {
            const size_t n = md_array_size(seg_type);
            for (size_t i = 0; i < n; ++i) {
                const int t = seg_type[i];
                if (t < 0) continue;
                if (!all) {
                    const md_atom_idx_t a = segment_atom(sys, i);
                    if (a < 0 || !md_bitset_test(atoms, (uint32_t)a)) continue;
                }
                md_array_push(l.segs[t], (uint32_t)i, heap);
            }
        }
        uint64_t h = 0x51ED27;
        for (int t = 0; t < NUM_TYPES; ++t) {
            h = md_hash64(l.segs[t], md_array_bytes(l.segs[t]), h ^ (uint64_t)(t + 1));
        }
        l.segs_hash = h;
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

    bool angles_available() const {
        const md_system_t& sys = app->mold.sys;
        return app->trajectory_data.backbone_angles.data && app->trajectory_data.backbone_angles.stride == sys.protein_backbone.segment.count &&
               md_array_size(seg_type) == sys.protein_backbone.segment.count;
    }

    // A layer of the system drawn as a density: whether it is to be shown, with a trajectory to show it of
    bool layer_draws_density(const Layer& l) const {
        if (l.source != LayerSource_System || l.style == LayerStyle_Points) return false;
        // Of a single structure, the frames are the frame: points say it better than a smoothing of them
        if (l.frames != LayerFrames_Current && run_num_frames(app) < 2) return false;
        return true;
    }

    void update() {
        if (!app) return;
        classify_segments();
        if (!show_window) return;

        char path_buf[256];
        const uint64_t angle_version = md_attributes_version(&app->mold.sys.attributes,
            md_attributes_id_from_path(run_attribute_path(path_buf, sizeof(path_buf), app, STR_LIT("backbone/angle"))));
        const bool available = angles_available();

        if (tree_dirty) apply_tree();

        // The groups' subsets, before what is within them. One that depends on the coordinates follows the frame only
        // when every layer in it is of the frame shown: a density over frames could not follow it.
        for (int gi = 0; gi < num_groups; ++gi) {
            LayerGroup& g = groups[gi];
            bool follow_frame = true;
            for (int i = 0; i < num_layers; ++i) {
                const Layer& l = layers[i];
                if (l.node.parent == g.node.id && l.source == LayerSource_System && l.frames != LayerFrames_Current) follow_frame = false;
            }
            layer_group_update(&g, app, md_hash64_combine(seg_type_key, system_generation), follow_frame, false);
        }

        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            if (l.source != LayerSource_System) continue;
            update_layer_segments(l);
            if (!layer_shown(l)) continue;

            if (!layer_draws_density(l) || !available || !l.segs_valid) continue;

            uint32_t beg = 0, end = 0;
            if (!layer_frame_range(l, &beg, &end)) continue;

            uint64_t h = md_hash64_combine(angle_version, l.segs_hash);
            h = md_hash64_combine(h, md_hash64_combine(seg_type_key, system_generation));
            const uint32_t frames[2] = { beg, end };
            h = md_hash64(frames, sizeof(frames), h);
            h = md_hash64(&l.bandwidth, sizeof(l.bandwidth), h);
            h |= 1;

            if (h == l.shown_hash || h == l.want_hash) continue;
            // One at a time: a job under way is let finish (its result is shown meanwhile), the latest is asked for after
            if (l.job) continue;
            launch_job(l, h, beg, end);
        }
    }

    // ## Workspace

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("Ramachandran"));
        viamd::write_int(state, STR_LIT("View"), view);
        viamd::write_bool(state, STR_LIT("ShowLayers"), show_sidebar);
        viamd::write_bool(state, STR_LIT("ShowTicks"), show_ticks);
        // Said even when there are none, so that none stays none (rather than the defaults) once read again
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
            layer_group_serialize(state, STR_LIT("RamachandranGroup"), g);
            int kids[MAX_LAYERS + LAYER_MAX_GROUPS];
            const size_t nk = representation_tree_children(kids, ARRAY_SIZE(kids), nodes, n, g.node.id);
            for (size_t c = 0; c < nk; ++c) {
                if (kids[c] < num_layers) serialize_layer(state, layers[kids[c]]);
            }
        }
    }

    void serialize_layer(viamd::serialization_state_t& state, const Layer& l) {
        {
            viamd::write_section_header(state, STR_LIT("RamachandranLayer"));
            if (l.node.parent) viamd::write_int(state, STR_LIT("Group"), (int64_t)l.node.parent);
            if (l.name[0]) viamd::write_str(state, STR_LIT("Name"), str_from_cstr(l.name));
            viamd::write_bool(state, STR_LIT("Visible"), l.visible);
            viamd::write_str(state, STR_LIT("Source"), l.source == LayerSource_Reference ? STR_LIT("reference") : STR_LIT("system"));
            if (l.source == LayerSource_System) {
                static const str_t frames[LayerFrames_Count] = { STR_LIT("current"), STR_LIT("all"), STR_LIT("period") };
                viamd::write_str(state, STR_LIT("Frames"), frames[l.frames]);
                if (l.frames == LayerFrames_Period) viamd::write_int(state, STR_LIT("Period"), (int)l.period);
                if (l.filter.text[0]) viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(l.filter.text));
                viamd::write_bool(state, STR_LIT("FollowSelection"), l.follow_selection);
            }
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
            static const str_t colorings[PointColoring_Count] = { STR_LIT("layer"), STR_LIT("region"), STR_LIT("time") };
            viamd::write_str(state, STR_LIT("PointColoring"), colorings[l.point_coloring]);
            viamd::write_flt(state, STR_LIT("PointSize"), l.point_size);
        }
    }

    void deserialize_settings(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "View")) {
                int v = -1;
                if (viamd::extract_int(v, arg)) view = CLAMP(v, -1, NUM_PANELS - 1);
            } else if (str_eq_cstr(ident, "ShowLayers")) {
                viamd::extract_bool(show_sidebar, arg);
            } else if (str_eq_cstr(ident, "ShowTicks")) {
                viamd::extract_bool(show_ticks, arg);
            } else if (str_eq_cstr(ident, "Layers")) {
                layers_declared = true;
            }
            // From before there were layers
            else if (str_eq_cstr(ident, "BlurSigma")) {
                legacy.present = true;
                viamd::extract_flt(legacy.blur_sigma, arg);
            } else if (str_eq_cstr(ident, "LayerAlpha")) {
                legacy.present = true;
                viamd::extract_flt_vec(legacy.alpha, 3, arg);
            } else if (str_eq_cstr(ident, "LayerDisplayMode")) {
                legacy.present = true;
                viamd::extract_int_vec(legacy.mode, 3, arg);
            } else if (str_eq_cstr(ident, "LayerColormap")) {
                legacy.present = true;
                viamd::extract_int_vec(legacy.colormap, 3, arg);
            } else if (str_eq_cstr(ident, "ShowLayer")) {
                legacy.present = true;
                viamd::extract_int_vec(legacy.show, 4, arg);
            } else if (str_eq_cstr(ident, "LayoutMode")) {
                legacy.present = true;
                viamd::extract_int(legacy.layout, arg);
            } else if (str_eq_cstr(ident, "LayerPeriod")) {
                legacy.present = true;
                viamd::extract_int(legacy.period, arg);
            } else {
                for (int i = 0; i < 3; ++i) {
                    char key[32];
                    snprintf(key, sizeof(key), "LayerIsolineColor%d", i);
                    if (str_eq_cstr(ident, key)) {
                        legacy.present = true;
                        viamd::extract_vec4(legacy.isoline_color[i], arg);
                        break;
                    }
                }
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
            // The first of the workspace: the ones there are go
            clear_layers();
            read_layers = true;
        }
        Layer* l = add_layer(LayerSource_System, LayerFrames_Current, LayerStyle_Points);
        str_t ident, arg;
        if (!l) {
            while (viamd::next_entry(ident, arg, state)) {}
            return;
        }
        l->point_coloring = PointColoring_Layer;
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
            } else if (str_eq_cstr(ident, "Source") && viamd::extract_str(s, arg)) {
                l->source = str_eq_cstr(s, "reference") ? LayerSource_Reference : LayerSource_System;
            } else if (str_eq_cstr(ident, "Frames") && viamd::extract_str(s, arg)) {
                l->frames = str_eq_cstr(s, "all") ? LayerFrames_All : str_eq_cstr(s, "period") ? LayerFrames_Period : LayerFrames_Current;
            } else if (str_eq_cstr(ident, "Period")) {
                int id = 0;
                if (viamd::extract_int(id, arg)) l->period = (uint32_t)MAX(id, 0);
            } else if (str_eq_cstr(ident, "Filter")) {
                viamd::extract_to_char_buf(l->filter.text, sizeof(l->filter.text), arg);
            } else if (str_eq_cstr(ident, "FollowSelection")) {
                viamd::extract_bool(l->follow_selection, arg);
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
                l->point_coloring = str_eq_cstr(s, "region") ? PointColoring_Region : str_eq_cstr(s, "time") ? PointColoring_Time : PointColoring_Layer;
            } else if (str_eq_cstr(ident, "PointSize")) {
                viamd::extract_flt(l->point_size, arg);
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
        l.source = (LayerSource)CLAMP((int)l.source, 0, LayerSource_Count - 1);
        l.frames = (LayerFrames)CLAMP((int)l.frames, 0, LayerFrames_Count - 1);
        l.style = (LayerStyle)CLAMP((int)l.style, 0, LayerStyle_Count - 1);
        if (l.source == LayerSource_Reference && l.style == LayerStyle_Points) l.style = LayerStyle_Contours;
        l.opacity = CLAMP(l.opacity, 0.0f, 1.0f);
        l.line_width = CLAMP(l.line_width, 0.5f, 8.0f);
        l.bandwidth = CLAMP(l.bandwidth, 0.5f, 30.0f);
        l.point_size = CLAMP(l.point_size, 1.0f, 16.0f);
        l.num_levels = CLAMP(l.num_levels, 1, MAX_LEVELS);
        for (int i = 0; i < l.num_levels; ++i) l.levels[i] = CLAMP(l.levels[i], 0.01f, 0.999f);
        sort_levels(l);
        if (l.colormap < 0 || l.colormap >= ImPlot::GetColormapCount()) l.colormap = ImPlotColormap_Viridis;
    }

    // The three layers and the points from before there were layers, as they were shown
    void create_legacy_layers() {
        clear_layers();
        Layer* cur = add_layer(LayerSource_System, LayerFrames_Current, LayerStyle_Points);
        if (cur) {
            cur->visible = legacy.show[3] != 0;
            cur->color = {1, 1, 1, 1};
            selected_layer = cur->id;
        }

        uint32_t period = (uint32_t)MAX(legacy.period, 0);
        if (period == 0 && app) period = app->timeline.periods.legacy_filter_id;
        const float bandwidth = legacy.blur_sigma > 0.0f ? CLAMP(legacy.blur_sigma * 1.4f, 1.0f, 20.0f) : 5.0f;

        // Each under the points, over the ones before it: the reference at the bottom, the period on top
        for (int i = 0; i < 3; ++i) {
            if (i == 2 && period == 0) continue;
            const LayerSource source = i == 0 ? LayerSource_Reference : LayerSource_System;
            const LayerFrames frames = i == 2 ? LayerFrames_Period : LayerFrames_All;
            const int mode = CLAMP(legacy.mode[i], 0, 2);
            Layer* l = add_layer(source, frames, mode == 2 ? LayerStyle_Heatmap : LayerStyle_Contours, i == 0 ? -1 : 1);
            if (!l) continue;
            l->visible = legacy.show[i] != 0;
            l->opacity = CLAMP(legacy.alpha[i], 0.0f, 1.0f);
            l->period = period;
            l->bandwidth = bandwidth;
            l->colormap = legacy.colormap[i] >= 0 && legacy.colormap[i] < ImPlot::GetColormapCount() ? legacy.colormap[i] : ImPlotColormap_Viridis;
            if (mode == 0) {
                l->fill = true;
                l->lines = false;
            } else if (mode == 1) {
                l->fill = false;
                l->lines = true;
                l->color = legacy.isoline_color[i];
            }
        }
        view = legacy.layout >= 1 && legacy.layout <= 4 ? legacy.layout - 1 : -1;
    }

    // ## Drawing

    // What a layer shows, as its name when it has none of its own
    void layer_label(char* buf, size_t cap, const Layer& l) const {
        if (l.name[0]) {
            snprintf(buf, cap, "%s", l.name);
            return;
        }
        if (l.source == LayerSource_Reference) {
            snprintf(buf, cap, "Reference");
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
        if (l.follow_selection) {
            if (l.frames == LayerFrames_Current) snprintf(buf, cap, "Selection");
            else snprintf(buf, cap, "Selection \xc2\xb7 %s", frames);
        } else if (l.filter.text[0]) {
            snprintf(buf, cap, "%s \xc2\xb7 %s", frames, l.filter.text);
        } else {
            snprintf(buf, cap, "%s", frames);
        }
    }

    void panel_hash_reset() {
        for (int i = 0; i < NUM_PANELS; ++i) density_map_canvas_invalidate(&canvas[i]);
    }

    // Fills what layer l is drawn with in the plot, false when there is nothing of it to draw there
    bool layer_draw(DensityMapDraw* d, const Layer& l, int panel) const {
        density_map_draw_init(d);
        const int type = panel == PANEL_ALL ? -1 : panel;
        d->opacity = l.opacity;
        d->line_width = l.line_width;
        MEMCPY(d->line_color, l.color.elem, sizeof(d->line_color));
        d->mode = l.style == LayerStyle_Heatmap ? 1 : 0;
        d->fill = l.fill;
        d->lines = l.lines;
        if (l.style == LayerStyle_Heatmap) {
            density_map_draw_colormap(d, l.colormap, l.log_scale);
        }

        if (l.source == LayerSource_Reference) {
            // The plot of all residues shows the general reference
            const int t = type < 0 ? Type_General : type;
            d->tex = ref_tex;
            d->weights[t] = 1.0f;
            d->num_levels = 2;
            d->levels[0] = ref_favored[t];
            d->levels[1] = ref_allowed[t];
            MEMCPY(d->fill_colors[0], ref_fill_favored[t].elem, sizeof(float) * 4);
            MEMCPY(d->fill_colors[1], ref_fill_allowed[t].elem, sizeof(float) * 4);
            d->cmap_max = 1.0f;
            d->content = 1;
            return true;
        }

        if (!layer_draws_density(l) || !l.has_result || !l.tex || !l.segs_valid || !angles_available()) return false;
        const DensityMapResult& r = l.result.map;
        if (r.samples[panel] == 0 || r.max[panel] <= 0.0f) return false;
        d->tex = l.tex;
        density_map_draw_levels(d, r, panel, l.levels, l.num_levels);
        density_map_draw_fills_from(d, vec_cast(l.color));
        d->content = l.shown_hash;
        return true;
    }

    // The densities of the layers in the plot, composited (density_map_composite): 0 when there are none
    uint32_t render_panel(int panel) {
        DensityMapDraw draws[MAX_LAYERS];
        int num_draws = 0;
        // Bottom up: the first layer is on top
        for (int i = num_layers - 1; i >= 0; --i) {
            const Layer& l = layers[i];
            if (!layer_shown(l) || l.style == LayerStyle_Points || l.opacity <= 0.0f) continue;
            if (layer_draw(&draws[num_draws], l, panel)) ++num_draws;
        }
        int w = 0, h = 0;
        density_map_plot_canvas_size(&w, &h);
        if (num_draws == 0 || w <= 0 || h <= 0) return 0;
        float view_uv[4];
        density_map_plot_view_uv(view_uv, rama_domain);
        return density_map_composite(&canvas[panel], w, h, view_uv, draws, num_draws);
    }

    // ### Points

    // Where the points of a layer come from
    struct PointSource {
        const md_backbone_angles_t* base = nullptr;
        size_t   stride = 0;        // Between frames
        uint32_t frame_beg = 0;
        uint32_t frame_end = 0;
        uint32_t frame_step = 1;
        bool     current = false;   // Of the state shown (interpolated), not of a frame
    };

    bool point_source(PointSource* ps, Layer& l, int panel) const {
        if (l.source != LayerSource_System || l.style != LayerStyle_Points || !l.segs_valid) return false;
        size_t num_segs = 0;
        for (int t = 0; t < NUM_TYPES; ++t) {
            if (panel == PANEL_ALL || panel == t) num_segs += md_array_size(l.segs[t]);
        }
        if (num_segs == 0) return false;

        // The frame shown, or the only one there is
        if (l.frames == LayerFrames_Current || run_num_frames(app) < 2) {
            ps->base = md_util_state_backbone_angles(&app->mold.state, &app->mold.sys);
            if (!ps->base || md_array_size(seg_type) != app->mold.sys.protein_backbone.segment.count) return false;
            ps->stride = 0;
            ps->frame_beg = 0;
            ps->frame_end = 1;
            ps->frame_step = 1;
            ps->current = true;
            l.point_stride = 1;
            return true;
        }
        if (!angles_available()) return false;
        uint32_t beg = 0, end = 0;
        if (!layer_frame_range(l, &beg, &end)) return false;
        ps->base = app->trajectory_data.backbone_angles.data;
        ps->stride = app->trajectory_data.backbone_angles.stride;
        ps->frame_beg = beg;
        ps->frame_end = end;
        const uint64_t total = (uint64_t)num_segs * (end - beg);
        ps->frame_step = (uint32_t)MAX((uint64_t)1, (total + MAX_POINTS - 1) / MAX_POINTS);
        ps->current = false;
        l.point_stride = ps->frame_step;
        return true;
    }

    template <typename F>
    void for_each_point(const Layer& l, const PointSource& ps, int panel, F&& f) const {
        for (uint32_t frame = ps.frame_beg; frame < ps.frame_end; frame += ps.frame_step) {
            const md_backbone_angles_t* a = ps.base + (size_t)frame * ps.stride;
            for (int t = 0; t < NUM_TYPES; ++t) {
                if (panel != PANEL_ALL && panel != t) continue;
                const uint32_t* segs = l.segs[t];
                const size_t n = md_array_size(segs);
                for (size_t i = 0; i < n; ++i) {
                    const md_backbone_angles_t ang = a[segs[i]];
                    if (!has_angles(ang)) continue;
                    f(t, segs[i], frame, ang);
                }
            }
        }
    }

    vec4_t point_color(const Layer& l, const PointSource& ps, int type, uint32_t frame, md_backbone_angles_t ang) const {
        switch (l.point_coloring) {
        case PointColoring_Region: {
            const Region r = ref_region(type, ang.phi, ang.psi);
            return r == Region_Favored ? l.color : region_colors[r];
        }
        case PointColoring_Time: {
            if (ps.frame_end - ps.frame_beg <= 1) return l.color;
            const float t = (float)(frame - ps.frame_beg) / (float)(ps.frame_end - ps.frame_beg - 1);
            return vec_cast(ImPlot::SampleColormap(t, l.colormap));
        }
        default:
            return l.color;
        }
    }

    // What the pointer is on
    struct Hover {
        int      layer = -1;
        uint32_t seg = 0;
        uint32_t frame = 0;
        int      type = 0;
        bool     current = true;
        md_backbone_angles_t ang = {};
    };

    void draw_tooltip(const Hover& h) const {
        const md_system_t& sys = app->mold.sys;
        const Layer& l = layers[h.layer];
        const md_residue_idx_t res_idx = sys.protein_backbone.segment.res_idx[h.seg];
        ImGui::BeginTooltip();
        if (res_idx >= 0) {
            const str_t name = md_system_residue_name(&sys, res_idx);
            const md_sequence_id_t seq_id = md_system_residue_seq_id(&sys, res_idx);
            const md_chain_idx_t chain = md_system_chain_find_by_residue_idx(&sys, (size_t)res_idx);
            if (chain >= 0) {
                const str_t chain_id = md_system_chain_id(&sys, (size_t)chain);
                ImGui::Text("%.*s %d  \xc2\xb7  chain %.*s", (int)name.len, name.ptr, seq_id, (int)chain_id.len, chain_id.ptr);
            } else {
                ImGui::Text("%.*s %d", (int)name.len, name.ptr, seq_id);
            }
        } else {
            ImGui::Text("CA atom[%d]", sys.protein_backbone.segment.atoms[h.seg].ca + 1);
        }
        ImGui::Text("\xcf\x86 %.1f\xc2\xb0   \xcf\x88 %.1f\xc2\xb0", RAD_TO_DEG(h.ang.phi), RAD_TO_DEG(h.ang.psi));
        const Region r = ref_region(h.type, h.ang.phi, h.ang.psi);
        ImGui::TextDisabled("%s \xc2\xb7", panel_names[h.type]);
        ImGui::SameLine();
        ImGui::TextColored(r == Region_Favored ? ImGui::GetStyleColorVec4(ImGuiCol_Text) : vec_cast(region_colors[r]), "%s", region_names[r]);
        if (!h.current) {
            ImGui::TextDisabled("frame %u", h.frame);
        }
        char label[128];
        layer_label(label, sizeof(label), l);
        ImGui::TextDisabled("%s", label);
        ImGui::EndTooltip();
    }

    void apply_to_selection(const md_bitmask_t* atoms, SelectionOperator op) {
        modify_field(&app->selection.selection_mask, md_bitmask_view(atoms), op);
    }

    // The points of the visible layers in the plot: hovered and boxed, then drawn
    void draw_points(int panel, bool hovered) {
        const md_system_t& sys = app->mold.sys;
        md_bitmask_t* selection_mask = &app->selection.selection_mask;
        md_bitmask_t* highlight_mask = &app->selection.highlight_mask;

        const ImPlotRect lim = ImPlot::GetPlotLimits();
        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 size = ImPlot::GetPlotSize();
        const double sx = size.x / lim.Size().x;
        const double sy = size.y / lim.Size().y;
        const double mid_x = (lim.X.Min + lim.X.Max) * 0.5;
        const double mid_y = (lim.Y.Min + lim.Y.Max) * 0.5;
        // Where the point is drawn: in the period around the middle of the view, and where the view is wider than the
        // domain also in the periods beside it
        const double margin_x = lim.Size().x * 0.02;
        const double margin_y = lim.Size().y * 0.02;
        const bool images_x = lim.Size().x > 360.0;
        const bool images_y = lim.Size().y > 360.0;
        auto for_images = [&](md_backbone_angles_t a, auto&& fn) {
            const double x = deperiodize_ortho(RAD_TO_DEG(a.phi), mid_x, 360.0);
            const double y = deperiodize_ortho(RAD_TO_DEG(a.psi), mid_y, 360.0);
            for (int ix = images_x ? -1 : 0; ix <= (images_x ? 1 : 0); ++ix) {
                const double xx = x + 360.0 * ix;
                if (xx < lim.X.Min - margin_x || xx > lim.X.Max + margin_x) continue;
                for (int iy = images_y ? -1 : 0; iy <= (images_y ? 1 : 0); ++iy) {
                    const double yy = y + 360.0 * iy;
                    if (yy < lim.Y.Min - margin_y || yy > lim.Y.Max + margin_y) continue;
                    fn(ImVec2((float)(pos.x + (xx - lim.X.Min) * sx), (float)(pos.y + (lim.Y.Max - yy) * sy)));
                }
            }
        };

        PointSource sources[MAX_LAYERS];
        bool has_source[MAX_LAYERS] = {};
        for (int i = 0; i < num_layers; ++i) {
            Layer& l = layers[i];
            if (!layer_shown(l)) continue;
            has_source[i] = point_source(&sources[i], l, panel);
        }

        // What the pointer is on, the nearest within reach, of the layer on top
        Hover hover;
        if (hovered && !box.active) {
            const ImVec2 mouse = ImGui::GetMousePos();
            const float reach = 8.0f;
            float best = reach * reach;
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                const float r = MAX(reach, layers[i].point_size);
                const float r2 = r * r;
                for_each_point(layers[i], sources[i], panel, [&](int t, uint32_t seg, uint32_t frame, md_backbone_angles_t ang) {
                    for_images(ang, [&](ImVec2 p) {
                        const float dx = p.x - mouse.x;
                        const float dy = p.y - mouse.y;
                        const float d2 = dx * dx + dy * dy;
                        if (d2 < best && d2 < r2) {
                            best = d2;
                            hover = Hover{ i, seg, frame, t, sources[i].current, ang };
                        }
                    });
                });
                // The layer on top has it, where it has a point at all
                if (hover.layer >= 0) break;
            }
        }

        // The highlight: what the box covers, or what is under the pointer
        if (hovered || (box.active && box.panel == panel)) {
            md_bitmask_clear(highlight_mask);
            highlights_now = true;
        }
        // Not moved from where it was pressed: a click, on the point pressed on (if any)
        bool box_is_click = false;
        if (box.active && box.panel == panel) {
            const ImVec2 a = ImPlot::PlotToPixels(box.p0);
            const ImVec2 m = ImGui::GetMousePos();
            box_is_click = fabsf(a.x - m.x) < 3.0f && fabsf(a.y - m.y) < 3.0f;
        }
        if (box.active && box.panel == panel && box_is_click) {
            if (box.click_seg >= 0) {
                mask_add_segment(highlight_mask, sys, (size_t)box.click_seg);
                grow_mask_by_selection_granularity(highlight_mask, app->selection.granularity, sys);
            }
        } else if (box.active && box.panel == panel) {
            const double x0 = MIN(box.p0.x, box.p1.x), x1 = MAX(box.p0.x, box.p1.x);
            const double y0 = MIN(box.p0.y, box.p1.y), y1 = MAX(box.p0.y, box.p1.y);
            const double cx = (x0 + x1) * 0.5, cy = (y0 + y1) * 0.5;
            for (int i = 0; i < num_layers; ++i) {
                if (!has_source[i]) continue;
                for_each_point(layers[i], sources[i], panel, [&](int, uint32_t seg, uint32_t, md_backbone_angles_t ang) {
                    const double x = deperiodize_ortho(RAD_TO_DEG(ang.phi), cx, 360.0);
                    const double y = deperiodize_ortho(RAD_TO_DEG(ang.psi), cy, 360.0);
                    if (x0 <= x && x <= x1 && y0 <= y && y <= y1) {
                        mask_add_segment(highlight_mask, sys, seg);
                    }
                });
            }
            grow_mask_by_selection_granularity(highlight_mask, app->selection.granularity, sys);
        } else if (hover.layer >= 0) {
            mask_add_segment(highlight_mask, sys, hover.seg);
            grow_mask_by_selection_granularity(highlight_mask, app->selection.granularity, sys);
        }

        // Drawn bottom up, then what is selected and highlighted over all of them
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        density_map_marks_begin(dl);
        const vec4_t sel_col = app->selection.color.selection.outline;
        const vec4_t hl_col = app->selection.color.highlight.outline;
        for (int i = num_layers - 1; i >= 0; --i) {
            if (!has_source[i]) continue;
            Layer& l = layers[i];
            const PointSource& ps = sources[i];
            const float r = l.point_size * 0.5f;
            const bool many = !ps.current;   // Over frames: small squares, no outlines
            const float a = l.opacity;
            for_each_point(l, ps, panel, [&](int t, uint32_t, uint32_t frame, md_backbone_angles_t ang) {
                vec4_t c = point_color(l, ps, t, frame, ang);
                c.w *= a;
                const ImU32 col = ImGui::ColorConvertFloat4ToU32(vec_cast(c));
                // Counted by class and in all, whichever plots draw them
                const Region reg = ref_region(t, ang.phi, ang.psi);
                l.point_regions[t][reg] += 1;
                l.point_count[t] += 1;
                l.point_regions[PANEL_ALL][reg] += 1;
                l.point_count[PANEL_ALL] += 1;
                for_images(ang, [&](ImVec2 p) {
                    if (many) {
                        density_map_mark(dl, p, r, DensityMapMark_Square, col);
                    } else {
                        density_map_mark(dl, p, r + 1.0f, DensityMapMark_Disc, col);
                    }
                });
            });
        }
        // Rings around the points of what is selected and highlighted: of every layer, as the residue is the same
        const bool any_selected = !md_bitmask_empty(selection_mask);
        const bool any_highlighted = !md_bitmask_empty(highlight_mask);
        if (any_selected || any_highlighted) {
            for (int i = num_layers - 1; i >= 0; --i) {
                if (!has_source[i]) continue;
                const Layer& l = layers[i];
                const float r = MAX(l.point_size * 0.5f + 2.5f, 4.0f);
                for_each_point(l, sources[i], panel, [&](int, uint32_t seg, uint32_t, md_backbone_angles_t ang) {
                    const md_atom_idx_t atom = segment_atom(sys, seg);
                    if (atom < 0) return;
                    const bool hl = any_highlighted && md_bitmask_test_bit(highlight_mask, (uint32_t)atom);
                    const bool sel = any_selected && md_bitmask_test_bit(selection_mask, (uint32_t)atom);
                    if (!hl && !sel) return;
                    const ImU32 c = ImGui::ColorConvertFloat4ToU32(vec_cast(hl ? hl_col : sel_col));
                    for_images(ang, [&](ImVec2 p) {
                        density_map_mark(dl, p, r + 1.0f, DensityMapMark_Ring, c);
                    });
                });
            }
        }
        density_map_marks_end(dl);
        ImPlot::PopPlotClipRect();

        if (hover.layer >= 0) {
            draw_tooltip(hover);
        }

        // Selecting: shift-drag a box, shift-click a point; the right button takes away; a right click on nothing clears
        const bool shift = ImGui::GetIO().KeyShift;
        const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
        if (!box.active && hovered && shift && (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
            box.active = true;
            box.panel = panel;
            box.op = ImGui::IsMouseClicked(ImGuiMouseButton_Left) ? SelectionOperator::Or : SelectionOperator::AndNot;
            box.p0 = box.p1 = mouse;
            box.click_seg = hover.layer >= 0 ? (int64_t)hover.seg : -1;
        }
        if (box.active && box.panel == panel) {
            box.p1 = mouse;
            const ImVec2 a = ImPlot::PlotToPixels(box.p0);
            const ImVec2 b = ImPlot::PlotToPixels(box.p1);
            const bool is_click = box_is_click;
            if (!is_click) {
                const vec4_t c = box.op == SelectionOperator::Or ? sel_col : vec4_t{1.0f, 0.4f, 0.4f, 0.9f};
                ImPlot::PushPlotClipRect();
                dl->AddRectFilled(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 0.12f)));
                dl->AddRect(ImMin(a, b), ImMax(a, b), ImGui::ColorConvertFloat4ToU32(vec_cast(c)), 0.0f, 0, 1.0f);
                ImPlot::PopPlotClipRect();
            }
            const bool released = box.op == SelectionOperator::Or ? ImGui::IsMouseReleased(ImGuiMouseButton_Left) : ImGui::IsMouseReleased(ImGuiMouseButton_Right);
            const bool held = box.op == SelectionOperator::Or ? ImGui::IsMouseDown(ImGuiMouseButton_Left) : ImGui::IsMouseDown(ImGuiMouseButton_Right);
            if (released || !held) {
                if (!is_click || box.click_seg >= 0) {
                    apply_to_selection(highlight_mask, box.op);
                } else if (box.op == SelectionOperator::AndNot) {
                    md_bitmask_clear(selection_mask);
                }
                box.active = false;
                box.panel = -1;
                box.click_seg = -1;
            }
        }
    }

    // Ticks at whole multiples of a step that fits the view, labelled within (-180, 180]
    // Ticks at whole multiples of a step that fits the view (the grid lines), returns how many
    static int setup_ticks(double* values, int cap, ImAxis axis, double lo, double hi) {
        const double span = hi - lo;
        const double step = span > 200.0 ? 90.0 : span > 90.0 ? 45.0 : span > 40.0 ? 15.0 : 5.0;
        int n = 0;
        for (double v = ceil(lo / step) * step; v <= hi && n < cap; v += step) {
            values[n++] = v;
        }
        if (n > 0) ImPlot::SetupAxisTicks(axis, values, n, nullptr, false);
        return n;
    }

    // An angle as it is labelled: within [-180, 180), but 180 where the domain ends on the right
    static int format_angle(double v, char* buf, int size, void*) {
        double w = v - 360.0 * floor((v + 180.0) / 360.0);
        if (v > 0.0 && w == -180.0) w = 180.0;
        return snprintf(buf, size, "%g\xc2\xb0", w);
    }

    void draw_panel(int panel, ImVec2 size, bool small) {
        constexpr ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
        // No labels outside the plot: they are drawn within it (DrawInsetAxes)
        const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_Foreground | ImPlot::InsetAxisFlags;

        char id[32];
        snprintf(id, sizeof(id), "##rama_%d", panel);
        if (!ImPlot::BeginPlot(id, size, flags)) return;

        ImPlot::SetupAxesLimits(-180, 180, -180, 180, ImPlotCond_Once);
        ImPlot::SetupAxisLinks(ImAxis_X1, &viewrect.X.Min, &viewrect.X.Max);
        ImPlot::SetupAxisLinks(ImAxis_Y1, &viewrect.Y.Min, &viewrect.Y.Max);
        ImPlot::SetupAxes("\xcf\x86", "\xcf\x88", axis_flags, axis_flags);
        ImPlot::SetupAxisFormat(ImAxis_X1, format_angle);
        ImPlot::SetupAxisFormat(ImAxis_Y1, format_angle);
        double ticks[64];
        setup_ticks(ticks, 64, ImAxis_X1, viewrect.X.Min, viewrect.X.Max);
        setup_ticks(ticks, 64, ImAxis_Y1, viewrect.Y.Min, viewrect.Y.Max);
        ImPlot::SetupFinish();

        const ImPlotRect lim = ImPlot::GetPlotLimits();
        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 psize = ImPlot::GetPlotSize();
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        if (psize.x > 1.0f && psize.y > 1.0f) plot_aspect = psize.x / psize.y;

        // The densities, composited, under everything else
        density_map_plot_canvas(render_panel(panel));

        // The edges of the periodic domain, where the view goes past them
        {
            ImPlot::PushPlotClipRect();
            const ImU32 edge = ImGui::GetColorU32(ImGuiCol_TextDisabled, 0.5f);
            for (double e = ceil((lim.X.Min - 180.0) / 360.0) * 360.0 + 180.0; e < lim.X.Max; e += 360.0) {
                if (e <= lim.X.Min) continue;
                const float x = ImPlot::PlotToPixels(e, 0).x;
                dl->AddLine(ImVec2(x, pos.y), ImVec2(x, pos.y + psize.y), edge, 1.0f);
            }
            for (double e = ceil((lim.Y.Min - 180.0) / 360.0) * 360.0 + 180.0; e < lim.Y.Max; e += 360.0) {
                if (e <= lim.Y.Min) continue;
                const float y = ImPlot::PlotToPixels(0, e).y;
                dl->AddLine(ImVec2(pos.x, y), ImVec2(pos.x + psize.x, y), edge, 1.0f);
            }
            ImPlot::PopPlotClipRect();
        }

        const bool hovered = ImPlot::IsPlotHovered();
        draw_points(panel, hovered);

        // The class and how many residues of it there are, in the corner; the labels of the axes along the edges
        ImVec2 title_size = ImVec2(0, 0);
        {
            char text[64];
            snprintf(text, sizeof(text), "%s  %u", panel_names[panel], seg_count[panel]);
            const float fs = ImGui::GetFontSize() * (small ? 0.85f : 1.0f);
            const ImVec2 p = ImVec2(pos.x + 6, pos.y + 4);
            const ImVec2 ts = ImGui::GetFont()->CalcTextSizeA(fs, FLT_MAX, 0.0f, text);
            ImPlot::PushPlotClipRect();
            dl->AddRectFilled(ImVec2(p.x - 3, p.y - 1), ImVec2(p.x + ts.x + 3, p.y + ts.y + 1), IM_COL32(0, 0, 0, 120), 3.0f);
            dl->AddText(ImGui::GetFont(), fs, p, IM_COL32(255, 255, 255, 225), text);
            ImPlot::PopPlotClipRect();
            title_size = ImVec2(ts.x + 9, ts.y + 6);
        }
        if (psize.x > 120.0f && psize.y > 120.0f) {
            ImPlot::DrawInsetAxes(show_ticks ? 0 : ImPlotInsetFlags_NoTickLabels, title_size, 0.78f);
        }

        // A period from the timeline: a new layer of the density over it
        if (ImPlot::BeginDragDropTargetPlot()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                ASSERT(payload->DataSize == sizeof(uint32_t));
                add_period_layer(*(const uint32_t*)payload->Data);
            }
            ImPlot::EndDragDropTarget();
        }

        // A double-click fits the view, as the plots have it
        if (ImPlot::FitThisFrame()) {
            fit_view = true;
        }

        ImPlot::EndPlot();
    }

    // In the group of what is selected, or in 'group' when given
    Layer* add_period_layer(uint32_t period_id, uint32_t group = UINT32_MAX) {
        const TimelinePeriod* p = timeline_periods_find(&app->timeline.periods, period_id);
        if (!p) return nullptr;
        Layer* l = add_layer(LayerSource_System, LayerFrames_Period, LayerStyle_Contours, 0, group == UINT32_MAX ? current_group() : group);
        if (!l) return nullptr;
        l->period = period_id;
        l->color = p->color;
        l->color.w = 1.0f;
        selected_layer = l->id;
        selected_group = 0;
        return l;
    }

    // The view: no wider than the domain and no narrower than 10 degrees, as the plots pan and zoom it
    void constrain_view() {
        if (viewrect.X.Max < viewrect.X.Min) ImSwap(viewrect.X.Min, viewrect.X.Max);
        if (viewrect.Y.Max < viewrect.Y.Min) ImSwap(viewrect.Y.Min, viewrect.Y.Max);
        const double mid_x = (viewrect.X.Min + viewrect.X.Max) * 0.5;
        const double mid_y = (viewrect.Y.Min + viewrect.Y.Max) * 0.5;
        // The shorter side within [10, 360] degrees, the other in proportion: all of the domain fits along the shorter
        // side, and the longer one shows some of it twice
        double w = MAX(viewrect.Size().x, 1e-3);
        double h = MAX(viewrect.Size().y, 1e-3);
        const double short_side = MIN(w, h);
        const double scl = short_side > 360.0 ? 360.0 / short_side : short_side < 10.0 ? 10.0 / short_side : 1.0;
        w *= scl;
        h *= scl;
        viewrect.X.Min = mid_x - w * 0.5;
        viewrect.X.Max = mid_x + w * 0.5;
        viewrect.Y.Min = mid_y - h * 0.5;
        viewrect.Y.Max = mid_y + h * 0.5;
        // Back into the domain once it has gone a whole period away
        if (mid_x > 360.0 || mid_x < -360.0) {
            const double s = 360.0 * floor((mid_x + 180.0) / 360.0);
            viewrect.X.Min -= s;
            viewrect.X.Max -= s;
        }
        if (mid_y > 360.0 || mid_y < -360.0) {
            const double s = 360.0 * floor((mid_y + 180.0) / 360.0);
            viewrect.Y.Min -= s;
            viewrect.Y.Max -= s;
        }
    }

    // The zoom of the wheel this frame, so that the view ends exactly at its limits (the whole domain along the
    // shorter side, out; 10 degrees, in) rather than going past them to be pulled back after: the plots would show the
    // view past the limit for a frame, and pulled back about its middle it would no longer be where the pointer zoomed.
    //
    // ImPlot zooms one step per frame whatever the wheel's amount: with equal axes, both sides by (1 + r/2)^2 for a rate
    // r, out, and by (1 - r/(2(1 + 2r)))^2 in (ImPlot's UpdateInput). Solved here for the rate that reaches the limit.
    void limit_zoom_rate(float* rate) const {
        const float wheel = ImGui::GetIO().MouseWheel;
        if (wheel == 0.0f || *rate <= 0.0f) return;
        const double short_side = MIN(viewrect.Size().x, viewrect.Size().y);
        if (short_side <= 0.0) return;
        double r = *rate;
        if (wheel < 0.0f) {
            const double f = 360.0 / short_side;                // What is left to zoom out
            r = f <= 1.0 + 1e-9 ? 0.0 : MIN(r, 2.0 * (sqrt(f) - 1.0));
        } else {
            const double f = 10.0 / short_side;                 // What is left to zoom in
            const double q = f >= 1.0 - 1e-9 ? 0.0 : 1.0 - sqrt(f);
            r = q <= 0.0 ? 0.0 : q < 0.25 ? MIN(r, 2.0 * q / (1.0 - 4.0 * q)) : r;
        }
        *rate = (float)r;
    }

    void draw_plots() {
        // Counted again as the points are drawn
        for (int i = 0; i < num_layers; ++i) {
            MEMSET(layers[i].point_regions, 0, sizeof(layers[i].point_regions));
            MEMSET(layers[i].point_count, 0, sizeof(layers[i].point_count));
        }

        // Fitted to the shape of the plots of the last frame
        const float aspect_before = plot_aspect;
        if (fit_view && plot_aspect > 0.0f) {
            const double w = plot_aspect >= 1.0f ? 360.0 * plot_aspect : 360.0;
            const double h = plot_aspect >= 1.0f ? 360.0 : 360.0 / plot_aspect;
            viewrect = ImPlotRect(-w * 0.5, w * 0.5, -h * 0.5, h * 0.5);
            fit_view = false;
            view_fitted = true;
        }
        const ImPlotRect view_before = viewrect;

        ImPlotInputMap& input = ImPlot::GetInputMap();
        const ImPlotInputMap saved_input = input;
        // Shift is for selecting: the plots do not pan or zoom meanwhile
        input.OverrideMod = ImGuiMod_Shift;
        limit_zoom_rate(&input.ZoomRate);

        // Edge to edge: the labels are within the plots. Side by side, a few pixels between them to tell them apart.
        ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, view < 0 ? ImVec2(2, 2) : ImVec2(1, 1));
        ImPlot::PushStyleVar(ImPlotStyleVar_LabelPadding, ImVec2(0, 0));
        // Opaque, whatever the windows are: what is behind the window is not to show through the data
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImVec4(0.11f, 0.11f, 0.12f, 1.0f));
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        if (view < 0) {
            const bool small = MIN(avail.x, avail.y) * 0.5f < 260.0f;
            constexpr ImPlotSubplotFlags subplot_flags = ImPlotSubplotFlags_NoResize | ImPlotSubplotFlags_NoMenus | ImPlotSubplotFlags_NoLegend;
            if (ImPlot::BeginSubplots("##rama_grid", 2, 2, ImVec2(-1, -1), subplot_flags)) {
                for (int t = 0; t < NUM_TYPES; ++t) {
                    draw_panel(t, ImVec2(), small);
                }
                ImPlot::EndSubplots();
            }
        } else {
            draw_panel(view, ImVec2(-1, -1), MIN(avail.x, avail.y) < 220.0f);
        }
        ImPlot::PopStyleVar(2);
        ImPlot::PopStyleColor();

        // Panned or zoomed: no longer fitted. Fitted and changed shape: fitted again.
        if (view_fitted) {
            const bool moved = fabs(viewrect.X.Min - view_before.X.Min) > 1e-6 * 360 || fabs(viewrect.Y.Min - view_before.Y.Min) > 1e-6 * 360;
            const bool reshaped = fabsf(plot_aspect - aspect_before) > 1e-3f;
            if (reshaped) {
                fit_view = true;
            } else if (moved) {
                view_fitted = false;
            }
        }
        if (aspect_before == 0.0f) fit_view = true;
        constrain_view();

        input = saved_input;
    }

    // ### The window

    void draw_toolbar() {
        static const char* labels[] = { ICON_FA_TABLE_CELLS_LARGE, "All", "General", "Gly", "Pro", "Pre-Pro" };
        static const int values[] = { -1, PANEL_ALL, Type_General, Type_Glycine, Type_Proline, Type_PreProline };
        static const char* tooltips[] = {
            "The four classes side by side",
            "Every residue in one plot (with the general reference)",
            "Residues other than glycine and proline, and not before a proline",
            "Glycine",
            "Proline",
            "Residues before a proline",
        };

        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6, 2));
        if (ImGui::Button(show_sidebar ? ICON_FA_LAYER_GROUP "##layers" : ICON_FA_LAYER_GROUP "##layers_off")) {
            show_sidebar = !show_sidebar;
        }
        ImGui::SetItemTooltip(show_sidebar ? "Hide the layers" : "Show the layers");
        ImGui::SameLine();
        density_map_segmented("view", labels, values, (int)ARRAY_SIZE(labels), &view, tooltips);

        // On the right: fitting the view and the settings
        ImGui::SameLine();
        const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + ImGui::CalcTextSize(ICON_FA_GEAR).x + ImGui::GetStyle().FramePadding.x * 4 + ImGui::GetStyle().ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
        if (ImGui::Button(ICON_FA_EXPAND)) {
            fit_view = true;
        }
        ImGui::SetItemTooltip("Fit the view (or double-click a plot)");
        ImGui::SameLine();
        if (ImGui::Button(ICON_FA_GEAR)) ImGui::OpenPopup("##settings");
        if (ImGui::BeginPopup("##settings")) {
            ImGui::Checkbox("Tick labels", &show_ticks);
            ImGui::Separator();
            ImGui::TextDisabled("Shift-drag to select residues, with the right button to deselect them.");
            ImGui::TextDisabled("Drop a period from the timeline on a plot to add its density.");
            ImGui::EndPopup();
        }
        ImGui::PopStyleVar();
    }

    // Why a layer of the system has no residues
    const char* layer_error(const Layer& l) const {
        const LayerGroup* g = layer_group(l);
        if (g && layer_group_has_context(*g) && !g->valid) return "the group's subset does not evaluate";
        if (!l.follow_selection && l.filter.text[0] && l.filter.result.error[0]) return l.filter.result.error;
        return "no residues";
    }

    // Statistics of a layer of the system as one line: how many residues, and where they are
    void layer_stats(char* buf, size_t cap, const Layer& l) const {
        buf[0] = '\0';
        if (l.source == LayerSource_Reference) {
            snprintf(buf, cap, "favored and allowed regions");
            return;
        }
        if (!l.segs_valid) {
            snprintf(buf, cap, "%s", layer_error(l));
            return;
        }
        // Of the plot shown, or of all of them
        const int panel = view < 0 ? PANEL_ALL : view;
        uint32_t num_res = 0;
        for (int t = 0; t < NUM_TYPES; ++t) {
            if (panel == PANEL_ALL || panel == t) num_res += (uint32_t)md_array_size(l.segs[t]);
        }
        uint64_t count = 0;
        uint64_t regions[Region_Count] = {};
        if (l.style != LayerStyle_Points && !layer_draws_density(l)) {
            snprintf(buf, cap, "%u res \xc2\xb7 no trajectory", num_res);
            return;
        }
        if (l.style == LayerStyle_Points) {
            count = l.point_count[panel];
            for (int k = 0; k < Region_Count; ++k) regions[k] = l.point_regions[panel][k];
        } else if (l.has_result) {
            count = l.result.map.samples[panel];
            for (int k = 0; k < Region_Count; ++k) regions[k] = l.result.regions[panel][k];
        }
        if (count == 0) {
            snprintf(buf, cap, "%u residues", num_res);
            return;
        }
        const double inv = 100.0 / (double)count;
        if (regions[Region_Outlier] > 0) {
            snprintf(buf, cap, "%u res \xc2\xb7 %.0f%% fav \xc2\xb7 %.1f%% out", num_res, regions[Region_Favored] * inv, regions[Region_Outlier] * inv);
        } else {
            snprintf(buf, cap, "%u res \xc2\xb7 %.0f%% favored", num_res, regions[Region_Favored] * inv);
        }
    }

    // What a layer holds, in full: by class, and where its samples fall
    void draw_layer_tooltip(const Layer& l) const {
        char label[128];
        layer_label(label, sizeof(label), l);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(label);
        if (l.source == LayerSource_Reference) {
            ImGui::TextDisabled("Favored and allowed regions of high quality structures, as MolProbity validates them");
            ImGui::EndTooltip();
            return;
        }
        if (!l.segs_valid) {
            ImGui::TextColored(vec_cast(region_colors[Region_Outlier]), "%s", layer_error(l));
            ImGui::EndTooltip();
            return;
        }
        const bool points = l.style == LayerStyle_Points;
        if (!points && l.has_result) {
            ImGui::TextDisabled("%u frames, smoothed by %.1f\xc2\xb0", l.result.num_frames, l.bandwidth);
        }
        if (ImGui::BeginTable("##stats", 5, ImGuiTableFlags_SizingFixedFit)) {
            ImGui::TableSetupColumn("");
            ImGui::TableSetupColumn("residues");
            ImGui::TableSetupColumn(region_names[Region_Favored]);
            ImGui::TableSetupColumn(region_names[Region_Allowed]);
            ImGui::TableSetupColumn(region_names[Region_Outlier]);
            ImGui::TableHeadersRow();
            for (int p = 0; p < NUM_PANELS; ++p) {
                const uint32_t num_res = p == PANEL_ALL
                    ? (uint32_t)(md_array_size(l.segs[0]) + md_array_size(l.segs[1]) + md_array_size(l.segs[2]) + md_array_size(l.segs[3]))
                    : (uint32_t)md_array_size(l.segs[p]);
                if (num_res == 0 && p != PANEL_ALL) continue;
                uint64_t count = 0;
                uint64_t regions[Region_Count] = {};
                if (points) {
                    count = l.point_count[p];
                    for (int k = 0; k < Region_Count; ++k) regions[k] = l.point_regions[p][k];
                } else if (l.has_result) {
                    count = l.result.map.samples[p];
                    for (int k = 0; k < Region_Count; ++k) regions[k] = l.result.regions[p][k];
                }
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(panel_names[p]);
                ImGui::TableNextColumn();
                ImGui::Text("%u", num_res);
                for (int k = 0; k < Region_Count; ++k) {
                    ImGui::TableNextColumn();
                    if (count) ImGui::Text("%.1f%%", 100.0 * regions[k] / count);
                    else ImGui::TextDisabled("-");
                }
            }
            ImGui::EndTable();
        }
        if (points && l.point_stride > 1) {
            ImGui::TextDisabled("Every %u frames are shown", l.point_stride);
        }
        ImGui::EndTooltip();
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
        const bool selected = selected_layer == l.id && !selected_group;
        const ImRect row_rect(row_min, ImVec2(row_min.x + row_w, row_min.y + row_h));

        char label[128];
        layer_label(label, sizeof(label), l);

        // The row, a selectable under the eye, the swatch and the text
        ImGui::SetNextItemAllowOverlap();
        if (ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h))) {
            selected_layer = l.id;
            selected_group = 0;
        }
        const bool row_hovered = ImGui::IsItemHovered();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal | ImGuiHoveredFlags_NoSharedDelay) && !ImGui::IsMouseDragging(0)) {
            draw_layer_tooltip(l);
        }
        layer_tree_drag_source({LIST_ID, l.id, false}, label);
        if (ImGui::BeginPopupContextItem("##menu")) {
            selected_layer = l.id;
            selected_group = 0;
            if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false, num_layers < MAX_LAYERS)) act->duplicate = l.id;
            if (l.node.parent == 0) {
                if (ImGui::MenuItem(ICON_FA_FOLDER " Group", nullptr, false, num_groups < LAYER_MAX_GROUPS && l.source == LayerSource_System)) act->group = l.id;
                ImGui::SetItemTooltip("A new group of it, with its residues as the group's: the layers added to the group show them");
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
        // A period from the timeline: this layer over its frames
        if (l.source == LayerSource_System && ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                l.frames = LayerFrames_Period;
                l.period = *(const uint32_t*)payload->Data;
                selected_layer = l.id;
                selected_group = 0;
            }
            ImGui::EndDragDropTarget();
        }
        // A layer or group of the list: moved here
        LayerNodeRef src;
        RepresentationDrop drop;
        if (layer_tree_drop_target(row_rect, row_min.x, {LIST_ID, l.id, false}, false, nodes, num_nodes, &src, &drop)) {
            act->move = true;
            act->move_src = src;
            act->move_target = {LIST_ID, l.id, false};
            act->move_drop = drop;
        }
        if (row_hovered && l.source == LayerSource_System && l.frames == LayerFrames_Period) {
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
        const ImVec2 sw_mid = ImVec2((sw_min.x + sw_max.x) * 0.5f, (sw_min.y + sw_max.y) * 0.5f);
        const float vis_alpha = l.visible ? 1.0f : 0.35f;
        if (l.source == LayerSource_Reference) {
            const vec4_t a = ref_fill_allowed[0];
            const vec4_t f = ref_fill_favored[0];
            dl->AddRectFilled(sw_min, sw_max, ImGui::ColorConvertFloat4ToU32(ImVec4(a.x, a.y, a.z, a.w * vis_alpha)), 3.0f);
            dl->AddCircleFilled(sw_mid, line_h * 0.28f, ImGui::ColorConvertFloat4ToU32(ImVec4(f.x * 0.8f, f.y * 0.8f, f.z * 0.9f, f.w * vis_alpha)));
        } else {
            const DensityMapSwatch kind = l.style == LayerStyle_Points ? DensityMapSwatch_Points : l.style == LayerStyle_Contours ? DensityMapSwatch_Contours : DensityMapSwatch_Heatmap;
            density_map_swatch(dl, sw_min, sw_max, kind, vec_cast(l.color), l.fill, l.colormap, vis_alpha);
        }

        // The name, and below it what it holds
        char stats[128];
        layer_stats(stats, sizeof(stats), l);
        const ImU32 text_col = ImGui::GetColorU32(l.visible ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        // On the selection's colour the disabled one is too faint to read
        const ImU32 dim_col = selected ? ImGui::GetColorU32(ImGuiCol_Text, 0.7f) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        const ImVec4 clip = ImVec4(text_col_min.x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), text_col_min, text_col, label, nullptr, 0.0f, &clip);
        const bool invalid = l.source == LayerSource_System && !l.segs_valid;
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(text_col_min.x, text_col_min.y + line_h + 1),
                    invalid ? IM_COL32(255, 110, 110, 255) : dim_col, stats, nullptr, 0.0f, &clip);

        // Working on its density
        if (l.job) {
            const float t = (float)ImGui::GetTime();
            const ImVec2 c = ImVec2(row_min.x + row_w - pad - 6, row_min.y + 9);
            dl->PathArcTo(c, 4.5f, t * 6.0f, t * 6.0f + 4.0f, 10);
            dl->PathStroke(dim_col, 0, 1.5f);
        }

        ImGui::SetCursorScreenPos(ImVec2(row_min.x - indent, row_min.y + row_h + 2));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
        ImGui::PopID();
    }

    // How many residues of the system are within a group's subset
    uint32_t group_residue_count(const LayerGroup& g) const {
        const md_system_t& sys = app->mold.sys;
        uint32_t n = 0;
        const size_t count = md_array_size(seg_type);
        for (size_t i = 0; i < count; ++i) {
            if (seg_type[i] < 0) continue;
            const md_atom_idx_t a = segment_atom(sys, i);
            if (a >= 0 && md_bitset_test(g.atoms, (uint32_t)a)) ++n;
        }
        return n;
    }

    void draw_group(LayerGroup& g, RepresentationNode* const nodes[], size_t num_nodes, SidebarActions* act) {
        int num_in = 0;
        for (int i = 0; i < num_layers; ++i) num_in += layers[i].node.parent == g.node.id ? 1 : 0;
        const bool context = layer_group_has_context(g);
        const bool invalid = context && !g.valid;
        char summary[64];
        if (invalid)      snprintf(summary, sizeof(summary), "does not evaluate");
        else if (context) snprintf(summary, sizeof(summary), "%u res", group_residue_count(g));
        else              snprintf(summary, sizeof(summary), "%d layer%s", num_in, num_in == 1 ? "" : "s");

        // What ungrouping would not keep
        bool exact = true;
        for (int i = 0; i < num_layers; ++i) {
            const Layer& l = layers[i];
            if (l.node.parent == g.node.id && l.source == LayerSource_System) exact = exact && layer_group_keeps_filter(l.filter.text, l.follow_selection, g);
        }
        const char* warning = exact ? nullptr :
            "Some of its layers cannot be given a filter that keeps what they show: one following the selection within a filter, "
            "or one of their own within the selection or within a filter that depends on the coordinates. Those keep their own, "
            "over the whole system.";
        const LayerGroupRow row = layer_group_row(g, LIST_ID, selected_group == g.node.id, summary, invalid, warning);
        if (row.clicked) {
            selected_group = g.node.id;
            selected_layer = 0;
        }
        if (row.action != LayerGroupAction_None) {
            act->group_id = g.node.id;
            act->group_action = row.action;
        }
        // Its subset lit up
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

    // What is added goes first in the group of what is selected: a layer there shows the group's residues
    void draw_add_menu() {
        const uint32_t group = current_group();
        const LayerGroup* g = find_group(group);
        if (g) {
            char gl[96];
            layer_group_label(gl, sizeof(gl), *g);
            ImGui::TextDisabled("In %s", gl);
        }
        const bool room = num_layers < MAX_LAYERS;
        auto select = [this](Layer* l) {
            if (!l) return;
            selected_layer = l->id;
            selected_group = 0;
        };
        if (ImGui::MenuItem("Current frame", nullptr, false, room)) {
            select(add_layer(LayerSource_System, LayerFrames_Current, LayerStyle_Points, 0, group));
        }
        if (ImGui::MenuItem("All frames", nullptr, false, room)) {
            select(add_layer(LayerSource_System, LayerFrames_All, LayerStyle_Contours, 0, group));
        }
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
        if (ps.count == 0) ImGui::SetItemTooltip("Mark periods in the timeline (shift-drag) to show the residues over them");
        if (ImGui::MenuItem("Selection", nullptr, false, room)) {
            if (Layer* l = add_layer(LayerSource_System, LayerFrames_Current, LayerStyle_Points, 0, group)) {
                l->follow_selection = true;
                l->point_coloring = PointColoring_Layer;
                l->point_size = 6.0f;
                select(l);
            }
        }
        ImGui::SetItemTooltip(g ? "The residues selected within the group's, as the selection changes" : "The residues selected, as the selection changes");
        ImGui::Separator();
        if (ImGui::MenuItem("Reference", nullptr, false, room)) {
            select(add_layer(LayerSource_Reference, LayerFrames_All, LayerStyle_Contours, -1, 0));
        }
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_FA_FOLDER " Group", nullptr, false, num_groups < LAYER_MAX_GROUPS)) {
            RepresentationNode* nodes[MAX_LAYERS + LAYER_MAX_GROUPS];
            if (LayerGroup* ng = add_group(layer_tree_first_order(nodes, tree_nodes(nodes), 0))) {
                selected_group = ng->node.id;
                selected_layer = 0;
            }
        }
        ImGui::SetItemTooltip("A group of layers that show the same residues: its subset, which the layers in it are within");
    }

    void draw_levels(Layer& l) {
        prop::row("levels", "The contours enclose these fractions of the samples (highest density regions)");
        density_map_levels_editor(l.levels, &l.num_levels);
    }

    void draw_inspector(Layer& l) {
        prop::begin("ramachandran");
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

        prop::row("source");
        {
            int s = (int)l.source;
            static const int values[] = { LayerSource_Reference, LayerSource_System };
            if (density_map_segmented("##source", source_names, values, LayerSource_Count, &s, nullptr, true)) {
                l.source = (LayerSource)s;
                if (l.source == LayerSource_Reference && l.style == LayerStyle_Points) {
                    l.style = LayerStyle_Contours;
                    l.fill = true;
                    l.lines = false;
                }
                l.segs_key = 0;
            }
        }

        if (l.source == LayerSource_System) {
            prop::row("residues", "The residues whose CA the filter picks, or those selected");
            {
                const float btn_w = ImGui::GetFrameHeight() + 6;
                ImGui::SetNextItemWidth(MAX(1.0f, prop::control_width() - btn_w - 2));
                ImGui::BeginDisabled(l.follow_selection);
                const bool invalid = !l.follow_selection && !l.segs_valid && l.filter.text[0];
                if (invalid) ImGui::PushInvalid();
                const LayerGroup* g = layer_group(l);
                const bool in_context = g && layer_group_has_context(*g) && g->valid;
                char none[96];
                if (in_context) {
                    char gl[64];
                    layer_group_label(gl, sizeof(gl), *g);
                    snprintf(none, sizeof(none), "all of the group's (%s)", gl);
                } else {
                    snprintf(none, sizeof(none), "all residues");
                }
                if (filter_input("##filter", &l.filter, layer_group_filter_context(app, in_context ? g : nullptr), 0, nullptr, l.follow_selection ? "the selection" : none)) {
                    l.segs_key = 0;
                }
                if (invalid) ImGui::PopInvalid();
                ImGui::EndDisabled();
                ImGui::SameLine(0, 2);
                const bool pushed = l.follow_selection;
                if (pushed) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
                if (ImGui::Button(ICON_FA_ARROW_POINTER "##follow", ImVec2(btn_w, 0))) {
                    l.follow_selection = !l.follow_selection;
                    l.segs_key = 0;
                }
                if (pushed) ImGui::PopStyleColor();
                ImGui::SetItemTooltip(l.follow_selection ? "Following the selection: click to use the filter" : "Follow the selection instead of the filter");
                if (!l.follow_selection && l.filter.result.is_dynamic && l.frames != LayerFrames_Current) {
                    prop::hint("Depends on the coordinates: taken in the frame it was evaluated in");
                }
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
                if (l.frames != LayerFrames_Current && run_num_frames(app) < 2) {
                    prop::hint(l.style == LayerStyle_Points ? "No trajectory: the frame there is" : "No trajectory to find a density over");
                }
            }
        }

        prop::row("style");
        {
            int s = (int)l.style;
            if (l.source == LayerSource_Reference) {
                static const int values[] = { LayerStyle_Contours, LayerStyle_Heatmap };
                density_map_segmented("##style", style_names + 1, values, 2, &s, nullptr, true);
            } else {
                static const int values[] = { LayerStyle_Points, LayerStyle_Contours, LayerStyle_Heatmap };
                density_map_segmented("##style", style_names, values, LayerStyle_Count, &s, nullptr, true);
            }
            l.style = (LayerStyle)s;
        }

        // A heat map is coloured by its colour map, points by time too
        if (l.style != LayerStyle_Heatmap && !(l.style == LayerStyle_Points && l.point_coloring == PointColoring_Time)) {
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
            if (l.source == LayerSource_System) {
                draw_levels(l);
            } else {
                prop::hint("Favored: 98%%\nAllowed: 99.95%%, or 99.8%% for Gly, Pro and pre-Pro");
            }
        }
        if (l.style == LayerStyle_Heatmap) {
            prop::row("colormap");
            ImPlot::ColormapSelection("##cmap", &l.colormap);
            prop::row("log scale", "Four decades below the maximum");
            ImGui::Checkbox("##log", &l.log_scale);
        }
        if (l.source == LayerSource_System && l.style != LayerStyle_Points) {
            prop::row("bandwidth", "The standard deviation of the gaussian the density is smoothed with");
            ImGui::SliderFloat("##bandwidth", &l.bandwidth, 1.0f, 20.0f, "%.1f\xc2\xb0", ImGuiSliderFlags_AlwaysClamp);
        }
        if (l.source == LayerSource_System && l.style == LayerStyle_Points) {
            prop::row("point size");
            ImGui::SliderFloat("##point_size", &l.point_size, 1.0f, 12.0f, "%.0f px");
            prop::row("color by");
            int c = (int)l.point_coloring;
            if (ImGui::Combo("##coloring", &c, point_coloring_names, PointColoring_Count)) l.point_coloring = (PointColoring)c;
            if (l.point_coloring == PointColoring_Time) {
                prop::row("colormap");
                ImPlot::ColormapSelection("##cmap_time", &l.colormap);
            }
            if (l.frames != LayerFrames_Current && l.point_stride > 1) {
                prop::hint("Every %u frames, to keep to %u points", l.point_stride, (uint32_t)MAX_POINTS);
            }
        }

        ImGui::PopID();
        prop::end();
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
        if (num_layers == 0 && num_groups == 0) {
            ImGui::TextDisabled("No layers: add one with +");
        }
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
            prop::begin("ramachandran");
            layer_group_inspector(*g, app, "residues");
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
        draw_window();
        if (app && owns_highlight && !highlights_now) md_bitmask_clear(&app->selection.highlight_mask);
        owns_highlight = highlights_now;
    }

    void draw_window() {
        if (!show_window || !app) return;

        ImGui::SetNextWindowSize(ImVec2(760, 480), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Ramachandran", &show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar)) {
            draw_toolbar();
            if (show_sidebar) {
                ImGui::BeginChild("##sidebar", ImVec2(320, 0), ImGuiChildFlags_ResizeX | ImGuiChildFlags_Borders);
                draw_sidebar();
                ImGui::EndChild();
                ImGui::SameLine();
            }
            ImGui::BeginChild("##plots", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
            if (app->mold.sys.protein_backbone.segment.count == 0) {
                ImGui::TextDisabled("The system has no protein backbone: only the reference is shown");
            }
            draw_plots();
            ImGui::EndChild();
        }
        ImGui::End();
    }
};

static Ramachandran instance = {};

}  // namespace ramachandran
