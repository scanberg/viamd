#define IMGUI_DEFINE_MATH_OPERATORS

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_grid.h>
#include <core/md_os.h>
#include <core/md_spatial_acc.h>
#include <core/md_coord_stream.h>

#include <md_system.h>
#include <md_util.h>

#include <viamd.h>
#include <viamd_event.h>
#include <event.h>

#include <imgui_widgets.h>
#include <implot_internal.h>
#include <implot3d.h>

#include "void_analysis_core.h"

#include <float.h>

/*
    Void and nanopore characterization: the component.

    The computation lives in void_analysis_core.h - the distance field, the profile it is summarized into and the
    reductions of it, and connectivity through z - and is tested on its own under tests/. This file owns what needs
    an application: the window and its parameters, reading the system, spreading the field pass over the task
    system, and the 3D overlay.

    Residency: statistics are accumulated per tile while the tile is still in cache, so the full field is only
    materialized when explicitly asked for. At 0.5 nm voxels a 670 x 670 x 175 nm box is 2.4e9 voxels, which is
    9.6 GB as float - streaming is the default for a reason.
*/

namespace {

constexpr int   NUM_BINS  = 512;                      // Distance bins of the profile histogram
constexpr int   MAX_SLABS = 512;                      // Upper bound on the z resolution of the profile
constexpr int   SASA_MAX_POINTS = 1024;
constexpr float ANGSTROM_PER_NM = 10.0f;
constexpr double PI_D = 3.14159265358979323846;

// Angstrom^2 per Dalton to m^2 per gram
constexpr double SPECIFIC_AREA_SCALE = 1.0e-20 / 1.66053906892e-24;

// Opacity of the channel ribbon. The overlay pass blends on source alpha, so this is how much of the
// structure is still read through the tube.
constexpr uint32_t RIBBON_ALPHA = 0x60;

// Colors are ABGR. The chain is cool and the throat warm, so the one tight place on the route reads
// against the rest of it rather than as another bead.
constexpr uint32_t ROUTE_LINE_COLOR   = 0xff50f0a0;
constexpr uint32_t ROUTE_SPHERE_COLOR = 0xffe0b060;
constexpr uint32_t THROAT_COLOR       = 0xff3c96ff;

// Centre to centre spacing of the probe spheres, as a multiple of the local radius. At 1.5 two
// equal spheres overlap by half a radius.
constexpr float SPHERE_SPACING = 1.5f;
constexpr size_t MAX_ROUTE_SPHERES = 256;

enum RadiusSource {
    RadiusSource_Vdw = 0,       // Per atom radius as reported by the system
    RadiusSource_Uniform,       // One radius for every bead, for coarse grained models with no meaningful element
    RadiusSource_Count,
};

const char* radius_source_lbl[RadiusSource_Count] = {
    "Van der Waals",
    "Uniform",
};

// Which part of the box the scalar results are reported over. A film simulated with an open z axis
// sits in a box with vacuum above and below it, and a porosity averaged over that box measures how
// much vacuum the box was given rather than anything about the film. There is no safe default, so
// the region is a stated choice and the profiles always show the whole grid regardless of it.
enum Region {
    Region_Box = 0,     // The whole grid
    Region_Film,        // The film slab, located by the half density convention
    Region_Manual,      // A z range set by hand
    Region_Count,
};

const char* region_lbl[Region_Count] = {
    "Whole box",
    "Film slab (auto)",
    "Manual z range",
};

// How the V(R,z) surface is shown. The two views carry the same numbers: the surface is what the
// quantity is, the heatmap is what it is easier to read a value off.
enum SurfaceView {
    SurfaceView_Surface = 0,
    SurfaceView_Heatmap,
    SurfaceView_Hidden,
    SurfaceView_Count,
};

const char* surface_view_lbl[SurfaceView_Count] = {
    "Surface (3D)",
    "Heatmap (2D)",
    "Hidden",
};

// A slab volume is the arbitrary thickness the profile was cut at, so the normalized form is the one
// that compares across slabs and across runs. The absolute one is still worth having, because it is
// what integrates back to the accessible volume of the region.
enum SurfaceValue {
    SurfaceValue_Fraction = 0,
    SurfaceValue_Volume,
    SurfaceValue_Count,
};

const char* surface_value_lbl[SurfaceValue_Count] = {
    "V(R,z) / V(z)",
    "V(R,z) in nm^3",
};

// The surface is a quad per grid cell, drawn every frame. The slab and radius counts are the user's
// to raise and this is not the place to pay for it: past this many points per axis the grid is
// sampled rather than drawn whole, which changes nothing visible and keeps the plot interactive.
constexpr uint32_t SURF_MAX_DIM = 128;

struct Stats {
    uint64_t num_grid      = 0;      // Voxels in the grid
    uint64_t num_voxels    = 0;      // Of those, the ones the statistics considered - see the in cell mask
    uint64_t num_solid     = 0;      // Voxels inside a bead, i.e. d <= 0
    uint64_t num_clamped   = 0;      // Voxels with no bead within max_dist, reported at that limit
    double   d_min         = 0.0;
    double   d_max         = 0.0;

    // The field summarized as a histogram of the distance, resolved along z. Porosity, accessible
    // volume and the coarea surface area are all reductions of this one array, so they cannot
    // disagree with each other the way three separate passes could. See void_analysis_core.h.
    md_array(uint64_t) hist       = 0;   // [slab * NUM_BINS + bin], void voxels only
    md_array(uint64_t) slab_solid = 0;   // [slab]
    md_array(uint64_t) slab_total = 0;   // [slab]
    uint32_t num_slabs     = 0;

    double   bin_width     = 0.0;
    double   z_min         = 0.0;
    double   slab_height   = 0.0;
    double   voxel_volume  = 0.0;
    double   seconds       = 0.0;
};

struct SasaResult {
    double total_area    = 0.0;      // Angstrom^2, extrapolated to every bead
    double std_err       = 0.0;      // Angstrom^2, from the spread across the sampled beads
    double specific_area = 0.0;      // m^2/g
    double area_density  = 0.0;      // Angstrom^2 per Angstrom^3
    double probe         = 0.0;      // The probe radius the result belongs to
    double t_layer       = 0.0;
    size_t num_sampled   = 0;
    md_array(double) type_area = 0;  // Per atom type, extrapolated
    double seconds       = 0.0;
};

}  // namespace

struct VoidAnalysis : viamd::EventHandler {
    bool show_window = false;

    // Parameters, all in Angstrom to match the rest of mdlib. The design document is written in nm, so the UI
    // presents nm and converts on the way in.
    float voxel_spacing = 5.0f;      // 0.5 nm
    float t_layer       = 0.0f;      // Bound water layer added to every bead radius, swept rather than fitted
    float probe_radius  = 1.4f;      // Water sized by default
    float max_dist      = 320.0f;    // 32 nm, the range the field is resolved over
    float uniform_radius = 5.0f;
    float cell_ext      = 20.0f;     // Cell extent handed to the acceleration structure

    RadiusSource radius_source = RadiusSource_Vdw;

    bool  materialize = false;
    md_array(float) field = 0;
    md_grid_t grid = {};

    // Percolation: r_c, what it costs a probe to be larger than it, and where the route is limited
    float  channel_r_min     = 2.0f;   // Also what bounds the memory: voxels below it never enter
    int    channel_num_radii = 64;
    bool   has_perc          = false;
    double perc_seconds      = 0.0;
    channel_percolation_t perc = {};
    md_array(double) pc_r      = 0;    // The curves, in nm and volume fractions, for plotting
    md_array(double) pc_top    = 0;
    md_array(double) pc_bot    = 0;
    md_array(double) pc_void   = 0;
    md_array(double) pc_open   = 0;
    md_array(double) pc_span   = 0;
    md_array(double) pc_closed = 0;

    // The route the critical radius belongs to
    bool   has_route    = false;
    bool   show_route   = true;
    double route_length = 0.0;
    double route_radius = 0.0;
    md_array(vec4_t) route = 0;

    // The sweep tree. Kept because a branching pore is worth seeing, but it is no longer where any
    // number comes from - on a real network it is one component with thousands of branches.
    bool  has_channels      = false;
    double channel_seconds  = 0.0;
    channel_tree_t channels = {};
    md_array(float)    node_slot    = 0;   // Dendrogram column per branch, negative when not drawn
    float    num_slots     = 0.0f;
    uint32_t hovered_node  = CHANNEL_INVALID_INDEX;

    // Porosity and accessible volume
    int      z_slabs     = 128;      // z resolution of the profile, capped at the grid and MAX_SLABS
    Region   region      = Region_Box;
    float    film_frac   = 0.5f;     // Solid fraction, as a share of the interior value, at the film edge
    float    manual_z_lo = 0.0f;
    float    manual_z_hi = 0.0f;
    uint32_t region_beg  = 0;        // Resolved slab range the scalars are reported over
    uint32_t region_end  = 0;
    uint32_t film_beg    = 0;
    uint32_t film_end    = 0;
    double   film_solid  = 0.0;
    bool     has_film    = false;

    // Cached reductions. The curve and the heatmap are a few million bin sums, so they are rebuilt
    // when the result or the region changes rather than every frame; the per slab accessible
    // fraction is one sweep of the bins and can follow the probe slider.
    int    acc_num_radii = 64;
    SurfaceView  surface_view  = SurfaceView_Surface;
    SurfaceValue surface_value = SurfaceValue_Fraction;
    bool   curves_dirty  = true;
    bool   probe_dirty   = true;
    double surf_max      = 1.0;      // Largest value on the surface, for the colour scale
    uint32_t surf_nx     = 0;        // Sampled grid of the surface: slabs along x, radii along y
    uint32_t surf_ny     = 0;
    double surf_x0 = 0.0, surf_x1 = 1.0;   // World bounds of the sampled grid, nm
    double surf_y0 = 0.0, surf_y1 = 1.0;
    md_array(double) acc_r    = 0;   // Probe radius samples, nm
    md_array(double) acc_frac = 0;   // V(R)/V over the region
    md_array(float)  surf_x   = 0;   // [ix + iy * surf_nx], z of the slab
    md_array(float)  surf_y   = 0;   // ... probe radius
    md_array(float)  surf_z   = 0;   // ... the value
    md_array(float)  acc_heat = 0;   // [(surf_ny - 1 - iy) * surf_nx + ix], row 0 is the largest radius
    md_array(double) prof_z   = 0;   // Slab centres, nm
    md_array(double) prof_phi = 0;   // Porosity per slab
    md_array(double) prof_acc = 0;   // V(R,z)/V(z) per slab at the current probe
    md_array(double) hist_x   = 0;   // Distance histogram of the region, nm
    md_array(double) hist_y   = 0;

    // Surface area
    int   sasa_points   = 256;
    float sasa_fraction = 1.0f;      // Fraction of beads sampled; the estimate is extrapolated from them

    Stats stats = {};
    bool  has_result = false;
    SasaResult sasa = {};
    bool  has_sasa = false;
    char  error[256] = "";

    md_allocator_i* arena = nullptr;
    ApplicationState* app_state = nullptr;

    VoidAnalysis() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(app_state->allocator.persistent, MEGABYTES(1));
                break;
            }
            case viamd::EventType_ViamdShutdown:
                // The percolation result and the traced route are on the heap rather than the arena,
                // so destroying the arena does not take them with it.
                clear_percolation();
                md_arena_allocator_destroy(arena);
                arena = nullptr;
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Void Analysis", &show_window);
                break;
            case viamd::EventType_ViamdSystemFree:
                clear_result();
                clear_sasa();
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                if (!show_window) break;
                if (e.payload_type != viamd::EventPayloadType_ApplicationState) break;
                const bool want_route  = show_route && (has_route || (has_perc && perc.has_r_c));
                const bool want_branch = has_channels && hovered_node != CHANNEL_INVALID_INDEX;
                if (!want_route && !want_branch) break;
                const ApplicationState& state = *(ApplicationState*)e.payload;
                immediate::Scope scope(state.gfx.overlay, "void_channel_path");
                // Camera Z in world space. The overlay is rendered with an identity model matrix, so
                // this is the same space the grid is in.
                const vec3_t cam_axis = vec3_normalize(vec3_from_vec4(state.view.param.matrix.inv.view.col[2]));
                if (want_route)  draw_route_3d(scope);
                if (want_branch) draw_channel_path_3d(scope, hovered_node, cam_axis);
                break;
            }
            default:
                break;
            }
        }
    }

    void clear_result() {
        md_array_free(field, arena);
        field = 0;
        md_array_free(stats.hist, arena);
        md_array_free(stats.slab_solid, arena);
        md_array_free(stats.slab_total, arena);
        stats = {};
        md_array_free(acc_r, arena);
        md_array_free(acc_frac, arena);
        md_array_free(surf_x, arena);
        md_array_free(surf_y, arena);
        md_array_free(surf_z, arena);
        md_array_free(acc_heat, arena);
        md_array_free(prof_z, arena);
        md_array_free(prof_phi, arena);
        md_array_free(prof_acc, arena);
        md_array_free(hist_x, arena);
        md_array_free(hist_y, arena);
        acc_r = 0; acc_frac = 0; acc_heat = 0;
        surf_x = 0; surf_y = 0; surf_z = 0;
        surf_nx = surf_ny = 0;
        prof_z = 0; prof_phi = 0; prof_acc = 0;
        hist_x = 0; hist_y = 0;
        region_beg = region_end = 0;
        has_film = false;
        curves_dirty = true;
        probe_dirty  = true;
        has_result = false;
        clear_percolation();
        clear_channels();
    }

    void clear_channels() {
        channel_tree_free(&channels);
        md_array_free(node_slot, arena);
        node_slot    = 0;
        num_slots    = 0.0f;
        hovered_node = CHANNEL_INVALID_INDEX;
        has_channels = false;
    }

    void clear_percolation() {
        channel_percolation_free(&perc);
        md_array_free(pc_r, arena);
        md_array_free(pc_top, arena);
        md_array_free(pc_bot, arena);
        md_array_free(pc_void, arena);
        md_array_free(pc_open, arena);
        md_array_free(pc_span, arena);
        md_array_free(pc_closed, arena);
        pc_r = 0; pc_top = 0; pc_bot = 0; pc_void = 0; pc_open = 0; pc_span = 0; pc_closed = 0;
        has_perc = false;
        perc_seconds = 0.0;
        clear_route();
    }

    void clear_route() {
        md_array_free(route, md_get_heap_allocator());
        route = 0;
        has_route = false;
        route_length = 0.0;
        route_radius = 0.0;
    }

    void clear_sasa() {
        md_array_free(sasa.type_area, arena);
        sasa = {};
        has_sasa = false;
    }

    // Radii the geometry is actually built from: the bead, its bound water layer, and optionally a probe.
    // A colloid cannot displace the bound layer, so the excluded surface is bead plus layer.
    void fill_radii(float* out_radii, const md_system_t& sys, size_t count, double extra) const {
        if (radius_source == RadiusSource_Vdw) {
            md_atom_extract_radii(out_radii, 0, count, &sys.atom);
        } else {
            for (size_t i = 0; i < count; ++i) out_radii[i] = uniform_radius;
        }
        for (size_t i = 0; i < count; ++i) out_radii[i] += (float)((double)t_layer + extra);
    }

    // Axis aligned grid covering the unit cell, or the atoms when there is no cell.
    // The voxel count is chosen so the grid tiles the box exactly, which keeps a periodic field seamless.
    bool setup_grid(md_grid_t* out_grid, const md_system_state_t& state) {
        return void_field_grid(out_grid, &state.unitcell, state.x, state.y, state.z, state.num_atoms, voxel_spacing);
    }

    void compute() {
        clear_result();
        error[0] = '\0';

        const md_system_t&       sys   = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;

        if (state.num_atoms == 0) {
            snprintf(error, sizeof(error), "No system loaded");
            return;
        }

        if (!setup_grid(&grid, state)) {
            snprintf(error, sizeof(error), "Could not derive a grid from the system");
            return;
        }

        const size_t num_voxels = md_grid_num_points(&grid);
        if (materialize) {
            const size_t bytes = num_voxels * sizeof(float);
            if (bytes > GIGABYTES(4)) {
                snprintf(error, sizeof(error), "Materialized field would need %.1f GB, refusing", (double)bytes / (double)GIGABYTES(1));
                return;
            }
            md_array_resize(field, num_voxels, arena);
        }

        const md_tick_t t0 = md_tick_now();

        md_temp_scope_t temp_scope = md_temp_begin();
        defer { md_temp_end(temp_scope); };

        // Bead radii, offset by the hydration layer. A colloid cannot displace the bound layer, so the excluded
        // surface is the bead plus the layer, and every reported quantity is a function of it.
        float* radii = (float*)md_temp_alloc(temp_scope, state.num_atoms * sizeof(float));
        fill_radii(radii, sys, state.num_atoms, 0.0);

        md_coord_stream_t coords = md_coord_stream_from_soa(state.x, state.y, state.z, NULL, state.num_atoms);

        md_spatial_acc_t acc = {};
        acc.alloc = md_temp_allocator(temp_scope);
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &coords;
        desc.radii    = radii;
        desc.cell_ext = cell_ext;
        desc.unitcell = &state.unitcell;
        md_spatial_acc_init_desc(&acc, &desc);
        defer { md_spatial_acc_free(&acc); };

        const uint32_t num_tiles = void_field_num_tiles(&grid);

        // z resolution of the profile. Never finer than the grid: a slab thinner than a plane of
        // voxels is either empty or a duplicate of its neighbour, and neither is a measurement.
        const uint32_t num_slabs = (uint32_t)CLAMP(z_slabs, 1, MIN(grid.dim[2], MAX_SLABS));

        // The bins are the resolution in R of every accessible volume reported later, so they are
        // deliberately finer than the voxel spacing: R is a continuous parameter and the voxelization,
        // not the binning, is what should be limiting. The coarea derivative is the one reader that
        // needs a coarser window, and it widens its own rather than making everyone else share it.
        void_field_desc_t fdesc = {};
        fdesc.acc       = &acc;
        fdesc.cell      = &state.unitcell;
        fdesc.grid      = &grid;
        fdesc.max_dist  = (double)max_dist;
        fdesc.num_slabs = num_slabs;
        fdesc.num_bins  = NUM_BINS;
        fdesc.field     = field;

        // Per thread accumulators, merged once the range task has completed. Nothing is shared while it runs.
        const size_t num_threads = MAX((size_t)1, task_system::pool_num_threads() + 1);
        const size_t hist_stride = (size_t)num_slabs * NUM_BINS;
        void_field_accum_t* accum = (void_field_accum_t*)md_temp_alloc(temp_scope, num_threads * sizeof(void_field_accum_t));
        for (size_t i = 0; i < num_threads; ++i) {
            accum[i].hist  = (uint64_t*)md_temp_alloc(temp_scope, hist_stride * sizeof(uint64_t));
            accum[i].solid = (uint64_t*)md_temp_alloc(temp_scope, num_slabs * sizeof(uint64_t));
            accum[i].total = (uint64_t*)md_temp_alloc(temp_scope, num_slabs * sizeof(uint64_t));
            void_field_accum_reset(&accum[i], num_slabs, NUM_BINS);
        }

        task_system::ID task = task_system::create_pool_task(STR_LIT("Void distance field"), num_tiles,
            [&](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                const uint32_t ti = MIN((uint32_t)(num_threads - 1), thread_num);
                void_field_eval_tiles(&accum[ti], &fdesc, range_beg, range_end);
            }, 1);

        task_system::enqueue_task(task);
        task_system::task_wait_for(task);

        md_array_resize(stats.hist, hist_stride, arena);
        md_array_resize(stats.slab_solid, num_slabs, arena);
        md_array_resize(stats.slab_total, num_slabs, arena);

        void_field_accum_t merged = {};
        merged.hist  = stats.hist;
        merged.solid = stats.slab_solid;
        merged.total = stats.slab_total;
        void_field_accum_reset(&merged, num_slabs, NUM_BINS);
        for (size_t t = 0; t < num_threads; ++t) {
            void_field_accum_merge(&merged, &accum[t], num_slabs, NUM_BINS);
        }

        const void_profile_t prof = void_field_profile(&merged, &fdesc);

        stats.num_grid    = num_voxels;
        stats.num_voxels  = 0;
        stats.num_solid   = 0;
        for (uint32_t sl = 0; sl < num_slabs; ++sl) {
            stats.num_voxels += merged.total[sl];
            stats.num_solid  += merged.solid[sl];
        }
        stats.num_clamped = merged.num_clamped;
        if (merged.d_min <= merged.d_max) {
            stats.d_min = (double)merged.d_min;
            stats.d_max = (double)merged.d_max;
        } else {
            stats.d_min = 0.0;
            stats.d_max = 0.0;
        }

        stats.num_slabs    = num_slabs;
        stats.bin_width    = prof.bin_width;
        stats.z_min        = prof.z_min;
        stats.slab_height  = prof.slab_height;
        stats.voxel_volume = prof.voxel_volume;
        stats.seconds      = md_tick_to_seconds(md_tick_now() - t0);

        has_result = true;

        // A manual range the user has already set survives a recompute; an unset one starts as the
        // whole grid so dragging it narrows rather than starting from nothing.
        if (!(manual_z_hi > manual_z_lo)) {
            manual_z_lo = grid.origin.z;
            manual_z_hi = grid.origin.z + grid.spacing.z * (float)grid.dim[2];
        }
        update_region();
    }

    // A read only view of the accumulated histogram. Every scalar and every curve below goes through
    // it, so porosity, V(R) and the coarea area are the same reduction evaluated at different R
    // rather than three computations which have to be kept in agreement by hand.
    void_profile_t make_profile() const {
        void_profile_t p = {};
        if (!has_result) return p;
        p.hist         = stats.hist;
        p.solid        = stats.slab_solid;
        p.total        = stats.slab_total;
        p.num_slabs    = stats.num_slabs;
        p.num_bins     = NUM_BINS;
        p.bin_width    = stats.bin_width;
        p.z_min        = stats.z_min;
        p.slab_height  = stats.slab_height;
        p.voxel_volume = stats.voxel_volume;
        return p;
    }

    // Resolve the reporting region to a slab range. The film extent is located whichever region is
    // selected, because it is worth showing next to a whole box number as the reason that number may
    // not mean what it looks like.
    void update_region() {
        has_film   = false;
        region_beg = 0;
        region_end = stats.num_slabs;
        if (!has_result) return;

        const void_profile_t p = make_profile();
        has_film = void_profile_film_extent(&p, (double)film_frac, &film_beg, &film_end, &film_solid);

        if (region == Region_Film && has_film) {
            region_beg = film_beg;
            region_end = film_end;
        } else if (region == Region_Manual) {
            const double lo = MIN(manual_z_lo, manual_z_hi);
            const double hi = MAX(manual_z_lo, manual_z_hi);
            region_beg = void_profile_slab_at(&p, lo);
            region_end = MIN(void_profile_slab_at(&p, hi) + 1, stats.num_slabs);
            if (region_end <= region_beg) region_end = MIN(region_beg + 1, stats.num_slabs);
        }

        curves_dirty = true;
        probe_dirty  = true;
    }

    // V(R)/V over the region, the same per slab as a function of R, and the distance histogram of
    // the region. All of it is sums over bins - a few million of them - which is cheap enough to do
    // on demand and far too much to do every frame.
    void update_curves() {
        curves_dirty = false;
        if (!has_result) return;

        const void_profile_t p = make_profile();
        const size_t nr = (size_t)CLAMP(acc_num_radii, 4, 256);

        // From zero, so the left endpoint of the curve is the porosity, out to the largest clearance
        // in the box, where it reaches zero. The curve between them is the whole of what this
        // section reports; the scalars are two points read off it.
        const double r_hi_nm = MAX(stats.d_max, stats.bin_width) / (double)ANGSTROM_PER_NM;

        md_array_resize(acc_r, nr, arena);
        md_array_resize(acc_frac, nr, arena);
        for (size_t i = 0; i < nr; ++i) {
            acc_r[i]    = r_hi_nm * (double)i / (double)(nr - 1);
            acc_frac[i] = void_profile_accessible_fraction(&p, region_beg, region_end, acc_r[i] * (double)ANGSTROM_PER_NM);
        }

        const uint32_t ns = stats.num_slabs;
        md_array_resize(prof_z,   ns, arena);
        md_array_resize(prof_phi, ns, arena);
        for (uint32_t sl = 0; sl < ns; ++sl) {
            prof_z[sl]   = 0.5 * (void_profile_z_lo(&p, sl) + void_profile_z_hi(&p, sl)) / (double)ANGSTROM_PER_NM;
            prof_phi[sl] = void_profile_porosity(&p, sl, sl + 1);
        }

        // The V(R,z) grid: a histogram of the accessible volume per z slab over every radius. Both
        // views read from this - a surface over (z, R) and, transposed and flipped, a heatmap - so
        // rotating one to check a value against the other is comparing a quantity with itself.
        const uint32_t stride_x = (ns + SURF_MAX_DIM - 1) / SURF_MAX_DIM;
        const uint32_t stride_y = ((uint32_t)nr + SURF_MAX_DIM - 1) / SURF_MAX_DIM;
        surf_nx = (ns + stride_x - 1) / stride_x;
        surf_ny = ((uint32_t)nr + stride_y - 1) / stride_y;

        md_array_resize(surf_x, (size_t)surf_nx * surf_ny, arena);
        md_array_resize(surf_y, (size_t)surf_nx * surf_ny, arena);
        md_array_resize(surf_z, (size_t)surf_nx * surf_ny, arena);
        md_array_resize(acc_heat, (size_t)surf_nx * surf_ny, arena);

        const double nm3 = (double)ANGSTROM_PER_NM * (double)ANGSTROM_PER_NM * (double)ANGSTROM_PER_NM;
        surf_max = 0.0;
        for (uint32_t iy = 0; iy < surf_ny; ++iy) {
            const size_t i = MIN(nr - 1, (size_t)iy * stride_y);
            const double r = acc_r[i] * (double)ANGSTROM_PER_NM;
            for (uint32_t ix = 0; ix < surf_nx; ++ix) {
                const uint32_t sl = MIN(ns - 1, ix * stride_x);
                const double frac = void_profile_accessible_fraction(&p, sl, sl + 1, r);
                const double v    = (surface_value == SurfaceValue_Volume)
                                  ? frac * (double)void_profile_num_total(&p, sl, sl + 1) * stats.voxel_volume / nm3
                                  : frac;
                const size_t k = (size_t)ix + (size_t)iy * surf_nx;
                surf_x[k] = (float)prof_z[sl];
                surf_y[k] = (float)acc_r[i];
                surf_z[k] = (float)v;
                // Row 0 of a heatmap is drawn at the top of the bounds, so its rows run from the
                // largest radius down while the surface runs the other way.
                acc_heat[(size_t)(surf_ny - 1 - iy) * surf_nx + ix] = (float)v;
                if (v > surf_max) surf_max = v;
            }
        }
        if (!(surf_max > 0.0)) surf_max = 1.0;
        surf_x0 = surf_x[0];
        surf_x1 = surf_x[(size_t)surf_nx - 1];
        surf_y0 = surf_y[0];
        surf_y1 = surf_y[(size_t)(surf_ny - 1) * surf_nx];

        md_array_resize(hist_x, NUM_BINS, arena);
        md_array_resize(hist_y, NUM_BINS, arena);
        const double n_tot = (double)MAX((uint64_t)1, void_profile_num_total(&p, region_beg, region_end));
        for (int b = 0; b < NUM_BINS; ++b) {
            uint64_t c = 0;
            for (uint32_t sl = region_beg; sl < region_end; ++sl) {
                c += stats.hist[(size_t)sl * NUM_BINS + b];
            }
            hist_x[b] = ((double)b + 0.5) * stats.bin_width / (double)ANGSTROM_PER_NM;
            hist_y[b] = (double)c / n_tot;
        }
    }

    // V(R,z)/V(z) at the current probe radius. One sweep of the bins per slab, so this can follow
    // the slider rather than waiting for a button.
    void update_probe_profile() {
        probe_dirty = false;
        if (!has_result) return;
        const void_profile_t p = make_profile();
        md_array_resize(prof_acc, stats.num_slabs, arena);
        for (uint32_t sl = 0; sl < stats.num_slabs; ++sl) {
            prof_acc[sl] = void_profile_accessible_fraction(&p, sl, sl + 1, (double)probe_radius);
        }
    }

    // A(r) = -dV_acc/dr, and with |grad d| = 1 almost everywhere the coarea formula makes that
    // derivative the density of the distance histogram at r. One field therefore already carries the
    // accessible surface area at every probe radius, for the cost of a division.
    //
    // The window is widened to the voxel spacing: a window narrower than that only splits the same
    // samples into noisier buckets, and A(r) is not flat across it. It is still a discretized
    // estimate - expect a few percent, against a few tenths of a percent for the sampling estimator
    // - so it earns its place as an independent cross check rather than as the number to quote.
    //
    // Always over the whole grid, never the reporting region: the sampling estimator it is checked
    // against counts every bead in the system, and comparing it against a slab of the box would
    // report the difference between two regions as a disagreement between two methods.
    double coarea_area(double r) const {
        if (!has_result) return 0.0;
        const void_profile_t p = make_profile();
        const double voxel_min = (double)MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        return void_profile_density(&p, 0, stats.num_slabs, r, voxel_min) * stats.voxel_volume;
    }

    // The z profile, one row per slab, at the current probe radius.
    void export_profile_csv() const {
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("csv"))) return;
        str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};

        md_file_t file = {};
        if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            MD_LOG_ERROR("Void analysis: could not open '%.*s' for writing", (int)path.len, path.ptr);
            return;
        }
        defer { md_file_close(&file); };

        const void_profile_t p = make_profile();
        const double nm  = (double)ANGSTROM_PER_NM;
        const double nm3 = nm * nm * nm;

        md_file_printf(file, "# VIAMD void analysis, z profile\n");
        md_file_printf(file, "# hydration_layer_nm,%.6f\n", (double)t_layer / nm);
        md_file_printf(file, "# probe_radius_nm,%.6f\n",    (double)probe_radius / nm);
        md_file_printf(file, "# region_slabs,%u,%u\n", region_beg, region_end);
        md_file_printf(file, "z_lo_nm,z_hi_nm,n_voxels,n_solid,porosity,accessible_fraction,void_volume_nm3,accessible_volume_nm3\n");
        for (uint32_t sl = 0; sl < stats.num_slabs; ++sl) {
            const double v_slab = (double)void_profile_num_total(&p, sl, sl + 1) * stats.voxel_volume;
            md_file_printf(file, "%.6f,%.6f,%llu,%llu,%.8f,%.8f,%.8f,%.8f\n",
                void_profile_z_lo(&p, sl) / nm,
                void_profile_z_hi(&p, sl) / nm,
                (unsigned long long)void_profile_num_total(&p, sl, sl + 1),
                (unsigned long long)void_profile_num_solid(&p, sl, sl + 1),
                void_profile_porosity(&p, sl, sl + 1),
                void_profile_accessible_fraction(&p, sl, sl + 1, (double)probe_radius),
                v_slab * void_profile_porosity(&p, sl, sl + 1) / nm3,
                void_profile_accessible_volume(&p, sl, sl + 1, (double)probe_radius) / nm3);
        }
    }

    // The full V(R,z) surface, long form: one row per (radius, slab) pair, which is what a plotting
    // tool wants and what a matrix written as a grid of numbers is not.
    void export_accessible_csv() const {
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("csv"))) return;
        str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};

        md_file_t file = {};
        if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            MD_LOG_ERROR("Void analysis: could not open '%.*s' for writing", (int)path.len, path.ptr);
            return;
        }
        defer { md_file_close(&file); };

        const void_profile_t p = make_profile();
        const double nm  = (double)ANGSTROM_PER_NM;
        const double nm3 = nm * nm * nm;
        const size_t nr  = md_array_size(acc_r);

        md_file_printf(file, "# VIAMD void analysis, accessible volume V(R,z)\n");
        md_file_printf(file, "# hydration_layer_nm,%.6f\n", (double)t_layer / nm);
        md_file_printf(file, "# z_slab -1 is the reporting region as a whole\n");
        md_file_printf(file, "z_slab,z_center_nm,R_nm,accessible_fraction,accessible_volume_nm3\n");
        for (size_t i = 0; i < nr; ++i) {
            const double r = acc_r[i] * nm;
            md_file_printf(file, "-1,%.6f,%.6f,%.8f,%.8f\n",
                0.5 * (void_profile_z_lo(&p, region_beg) + void_profile_z_hi(&p, region_end - 1)) / nm,
                acc_r[i],
                void_profile_accessible_fraction(&p, region_beg, region_end, r),
                void_profile_accessible_volume(&p, region_beg, region_end, r) / nm3);
        }
        for (uint32_t sl = 0; sl < stats.num_slabs; ++sl) {
            for (size_t i = 0; i < nr; ++i) {
                const double r = acc_r[i] * nm;
                md_file_printf(file, "%u,%.6f,%.6f,%.8f,%.8f\n", sl,
                    0.5 * (void_profile_z_lo(&p, sl) + void_profile_z_hi(&p, sl)) / nm,
                    acc_r[i],
                    void_profile_accessible_fraction(&p, sl, sl + 1, r),
                    void_profile_accessible_volume(&p, sl, sl + 1, r) / nm3);
            }
        }
    }

    // Shrake-Rupley through the same nearest query: scatter points over a bead's expanded sphere and count those
    // nothing else buries. No voxels involved, and it converges as 1/sqrt(N) rather than with the grid.
    //
    // The exposure test is "the query returned this bead", not "the distance came back positive". A sample point
    // lies exactly on its own bead's expanded surface, so its own distance is zero up to rounding - and at the
    // coordinates of a 670 nm box a float carries about 1e-3 A of it there, which swamps any small outward offset
    // one might add to force the sign. The index test does not depend on the sign at all: the bead wins unless
    // something genuinely closer buries the point, and a periodic image of the bead carries the same index.
    void compute_sasa() {
        clear_sasa();
        error[0] = '\0';

        const md_system_t&       sys   = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;

        if (state.num_atoms == 0) {
            snprintf(error, sizeof(error), "No system loaded");
            return;
        }

        const md_tick_t t0 = md_tick_now();

        md_temp_scope_t temp_scope = md_temp_begin();
        defer { md_temp_end(temp_scope); };

        float* radii = (float*)md_temp_alloc(temp_scope, state.num_atoms * sizeof(float));
        fill_radii(radii, sys, state.num_atoms, (double)probe_radius);

        md_coord_stream_t coords = md_coord_stream_from_soa(state.x, state.y, state.z, NULL, state.num_atoms);

        md_spatial_acc_t acc = {};
        acc.alloc = md_temp_allocator(temp_scope);
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &coords;
        desc.radii    = radii;
        desc.cell_ext = cell_ext;
        desc.unitcell = &state.unitcell;
        md_spatial_acc_init_desc(&acc, &desc);
        defer { md_spatial_acc_free(&acc); };

        // Unit sphere directions, generated once and shared read only. A Fibonacci spiral spreads the points evenly
        // enough that the quadrature error is the 1/sqrt(N) term rather than a pattern in the lattice.
        const int NP = CLAMP(sasa_points, 8, SASA_MAX_POINTS);
        float* dir = (float*)md_temp_alloc(temp_scope, (size_t)NP * 3 * sizeof(float));
        {
            const double golden = PI_D * (3.0 - sqrt(5.0));
            for (int k = 0; k < NP; ++k) {
                const double cz  = 1.0 - 2.0 * ((double)k + 0.5) / (double)NP;
                const double rho = sqrt(MAX(0.0, 1.0 - cz * cz));
                const double th  = golden * (double)k;
                dir[k * 3 + 0] = (float)(rho * cos(th));
                dir[k * 3 + 1] = (float)(rho * sin(th));
                dir[k * 3 + 2] = (float)cz;
            }
        }

        // Every bead is equally likely to be sampled regardless of where it sits, so a stride is unbiased for the
        // total; correlation along a fibril costs variance, not accuracy, and the reported error bar sees it.
        const size_t stride      = MAX((size_t)1, (size_t)llround(1.0 / (double)CLAMP(sasa_fraction, 0.001f, 1.0f)));
        const size_t num_sampled = (state.num_atoms + stride - 1) / stride;
        const size_t num_types   = MAX((size_t)1, sys.atom.type.count);

        const size_t num_threads = MAX((size_t)1, task_system::pool_num_threads() + 1);
        double* acc_area  = (double*)md_temp_alloc(temp_scope, num_threads * sizeof(double));
        double* acc_area2 = (double*)md_temp_alloc(temp_scope, num_threads * sizeof(double));
        double* acc_mass  = (double*)md_temp_alloc(temp_scope, num_threads * sizeof(double));
        double* acc_type  = (double*)md_temp_alloc(temp_scope, num_threads * num_types * sizeof(double));
        MEMSET(acc_area,  0, num_threads * sizeof(double));
        MEMSET(acc_area2, 0, num_threads * sizeof(double));
        MEMSET(acc_mass,  0, num_threads * sizeof(double));
        MEMSET(acc_type,  0, num_threads * num_types * sizeof(double));

        // A sample point sits on its own bead, so the search never has to travel: anything which could bury it is
        // within two radii. Bounding the query here keeps the traversal to the immediate neighborhood.
        const double query_max_dist = 2.0 * (double)acc.max_rad + 1.0;

        task_system::ID task = task_system::create_pool_task(STR_LIT("Accessible surface area"), (uint32_t)num_sampled,
            [&](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                const size_t ti = MIN(num_threads - 1, (size_t)thread_num);
                double* my_type = acc_type + ti * num_types;

                float px[SASA_MAX_POINTS], py[SASA_MAX_POINTS], pz[SASA_MAX_POINTS];
                uint32_t idx[SASA_MAX_POINTS];

                for (uint32_t s = range_beg; s < range_end; ++s) {
                    const size_t i = (size_t)s * stride;
                    if (i >= state.num_atoms) continue;

                    const float a = radii[i];
                    if (!(a > 0.0f)) continue;

                    for (int k = 0; k < NP; ++k) {
                        px[k] = state.x[i] + a * dir[k * 3 + 0];
                        py[k] = state.y[i] + a * dir[k * 3 + 1];
                        pz[k] = state.z[i] + a * dir[k * 3 + 2];
                    }

                    md_coord_stream_t pts = md_coord_stream_from_soa(px, py, pz, NULL, (size_t)NP);
                    md_spatial_acc_query_nearest(&acc, &pts, query_max_dist, idx, NULL);

                    int exposed = 0;
                    for (int k = 0; k < NP; ++k) exposed += (idx[k] == (uint32_t)i);

                    const double area = 4.0 * PI_D * (double)a * (double)a * (double)exposed / (double)NP;
                    acc_area [ti] += area;
                    acc_area2[ti] += area * area;
                    acc_mass [ti] += (double)md_atom_mass(&sys.atom, i);

                    const size_t type = (size_t)md_atom_type_idx(&sys.atom, i);
                    if (type < num_types) my_type[type] += area;
                }
            }, 64);

        task_system::enqueue_task(task);
        task_system::task_wait_for(task);

        md_array_resize(sasa.type_area, num_types, arena);
        MEMSET(sasa.type_area, 0, num_types * sizeof(double));

        double sum_area = 0.0, sum_area2 = 0.0, sum_mass = 0.0;
        for (size_t t = 0; t < num_threads; ++t) {
            sum_area  += acc_area[t];
            sum_area2 += acc_area2[t];
            sum_mass  += acc_mass[t];
            for (size_t k = 0; k < num_types; ++k) sasa.type_area[k] += acc_type[t * num_types + k];
        }

        const double scale = (double)state.num_atoms / (double)num_sampled;
        for (size_t k = 0; k < num_types; ++k) sasa.type_area[k] *= scale;

        sasa.total_area  = sum_area * scale;
        sasa.num_sampled = num_sampled;

        // The sampled beads are a sample of the system, so their spread is the uncertainty of the extrapolation.
        // At full sampling this is zero, which is honest: what remains then is the quadrature error of the point
        // count, not a sampling error.
        if (num_sampled > 1 && stride > 1) {
            const double mean = sum_area / (double)num_sampled;
            const double var  = MAX(0.0, sum_area2 / (double)num_sampled - mean * mean);
            sasa.std_err = sqrt(var / (double)num_sampled) * (double)state.num_atoms;
        }

        // A ratio of two sums over the same sample, so the subsampling cancels rather than needing the scale factor
        sasa.specific_area = (sum_mass > 0.0) ? (sum_area / sum_mass) * SPECIFIC_AREA_SCALE : 0.0;

        double A[3][3];
        md_unitcell_A_extract_double(A, &state.unitcell);
        const double vol = fabs(A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
                              - A[1][0] * (A[0][1] * A[2][2] - A[2][1] * A[0][2])
                              + A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]));
        sasa.area_density = (vol > 0.0) ? sasa.total_area / vol : 0.0;

        sasa.probe   = (double)probe_radius;
        sasa.t_layer = (double)t_layer;
        sasa.seconds = md_tick_to_seconds(md_tick_now() - t0);

        has_sasa = true;
    }

    // --- Channels -----------------------------------------------------------------------------

    bool build_channel_field(channel_field_t* out) const {
        if (!field || !has_result) return false;
        MEMSET(out, 0, sizeof(*out));
        out->data = field;
        const uint32_t flags = md_unitcell_flags(&app_state->mold.state.unitcell);
        for (int a = 0; a < 3; ++a) {
            out->dim[a]     = grid.dim[a];
            out->spacing[a] = grid.spacing.elem[a];
            out->origin[a]  = grid.origin.elem[a];
            out->pbc[a]     = (flags & (MD_UNITCELL_PBC_X << a)) != 0;
        }
        // The sweep axis needs two faces to connect, so it is treated as open whatever the cell says
        out->pbc[2] = false;
        return true;
    }

    // One pass in order of decreasing clearance. It answers what the bisection and the per radius
    // sweeps answered - r_c, and how many components get through at each radius - and several things
    // they could not: where the route is limited, how deep a probe too large to cross still reaches,
    // and which of the pore space is open to the outside rather than merely large enough.
    void compute_percolation() {
        clear_percolation();
        error[0] = '\0';

        channel_field_t f;
        if (!build_channel_field(&f)) {
            snprintf(error, sizeof(error), "Compute the field first, with 'Materialize field' enabled");
            return;
        }

        const md_tick_t t0 = md_tick_now();

        if (!channel_percolate(&perc, &f, (double)channel_r_min, (uint32_t)CLAMP(channel_num_radii, 8, 512), md_get_heap_allocator())) {
            snprintf(error, sizeof(error), "Nothing at or above %.2f nm to analyse", channel_r_min / ANGSTROM_PER_NM);
            channel_percolation_free(&perc);
            return;
        }
        has_perc = true;

        // At exactly r_c the set is barely connected, so any route through it passes the throat the
        // pass just located. Tracing there is therefore tracing the route r_c belongs to, not a
        // route that merely happens to fit.
        if (perc.has_r_c) trace_route(perc.r_c);

        perc_seconds = md_tick_to_seconds(md_tick_now() - t0);
        build_perc_curves();
    }

    void trace_route(double r) {
        clear_route();
        channel_field_t f;
        if (!build_channel_field(&f)) return;
        has_route = channel_trace_path(&route, &route_length, &f, r, md_get_heap_allocator());
        route_radius = r;
        if (!has_route) {
            snprintf(error, sizeof(error), "No route at %.2f nm", r / ANGSTROM_PER_NM);
        }
    }

    void build_perc_curves() {
        const size_t n = md_array_size(perc.radius);
        md_array_resize(pc_r,      n, arena);
        md_array_resize(pc_top,    n, arena);
        md_array_resize(pc_bot,    n, arena);
        md_array_resize(pc_void,   n, arena);
        md_array_resize(pc_open,   n, arena);
        md_array_resize(pc_span,   n, arena);
        md_array_resize(pc_closed, n, arena);
        for (size_t i = 0; i < n; ++i) {
            pc_r[i]      = perc.radius[i]        / (double)ANGSTROM_PER_NM;
            pc_top[i]    = perc.z_from_top[i]    / (double)ANGSTROM_PER_NM;
            pc_bot[i]    = perc.z_from_bottom[i] / (double)ANGSTROM_PER_NM;
            pc_void[i]   = perc.frac_void[i];
            pc_open[i]   = perc.frac_open[i];
            pc_span[i]   = perc.frac_spanning[i];
            pc_closed[i] = perc.frac_void[i] - perc.frac_open[i];
        }
    }

    // Read a curve at an arbitrary radius, which is what the scalars below the plots quote.
    double perc_at(const md_array(double) y, double r_nm) const {
        const size_t n = md_array_size(pc_r);
        if (n == 0 || !y) return 0.0;
        if (r_nm <= pc_r[0])     return y[0];
        if (r_nm >= pc_r[n - 1]) return y[n - 1];
        size_t i = 0;
        while (i + 2 < n && pc_r[i + 1] < r_nm) ++i;
        const double t = (r_nm - pc_r[i]) / MAX(1.0e-12, pc_r[i + 1] - pc_r[i]);
        return (1.0 - t) * y[i] + t * y[i + 1];
    }

    // Nearest sample to a radius, for the integer curves where interpolating between two counts
    // would invent a component that is not there.
    size_t perc_index_at(double r_nm) const {
        const size_t n = md_array_size(pc_r);
        if (n == 0) return 0;
        size_t best = 0;
        double best_d = fabs(pc_r[0] - r_nm);
        for (size_t i = 1; i < n; ++i) {
            const double d = fabs(pc_r[i] - r_nm);
            if (d < best_d) { best_d = d; best = i; }
        }
        return best;
    }

    void compute_channels() {
        clear_channels();
        error[0] = '\0';

        channel_field_t f;
        if (!build_channel_field(&f)) {
            snprintf(error, sizeof(error), "Compute the field first, with 'Materialize field' enabled");
            return;
        }

        const md_tick_t t0 = md_tick_now();
        channel_sweep(&channels, &f, (double)probe_radius, true, arena);
        layout_channel_tree();
        channel_seconds = md_tick_to_seconds(md_tick_now() - t0);
        has_channels = true;
    }

    void layout_channel_tree() {
        const size_t num = md_array_size(channels.nodes);
        md_array_resize(node_slot, num, arena);
        num_slots = channel_tree_layout(node_slot, &channels, arena);
    }

    // Everything below is drawn into the overlay queue, which is rendered last with the depth test
    // off, so a route stays in front of the structure rather than being buried by it.
    // One stretch of a route, as a camera facing ribbon whose half width is the clearance there with
    // an opaque line on top. A line alone is one pixel wide at any zoom, which in a dense network is
    // easy to lose and thin enough for temporal AA to eat; the ribbon carries the width of the
    // channel and is translucent so the structure still reads through it.
    void draw_clearance_segment(immediate::Queue* q, vec4_t a, vec4_t b, uint32_t col, vec3_t cam_axis) const {
        // A route may wrap through the periodic faces; it must not be drawn straight back across the
        // box when it does.
        const float wrap_x = 0.5f * grid.spacing.x * (float)grid.dim[0];
        const float wrap_y = 0.5f * grid.spacing.y * (float)grid.dim[1];
        if (fabsf(a.x - b.x) > wrap_x || fabsf(a.y - b.y) > wrap_y) return;

        const vec3_t p0 = {a.x, a.y, a.z};
        const vec3_t p1 = {b.x, b.y, b.z};

        const vec3_t d   = vec3_sub(p1, p0);
        const float  len = vec3_length(d);
        if (len > 1.0e-6f) {
            // Perpendicular to both the segment and the line of sight. It only degenerates when the
            // segment points straight at the camera, where a ribbon has nothing to show.
            vec3_t side = vec3_cross(vec3_mul1(d, 1.0f / len), cam_axis);
            const float side_len = vec3_length(side);
            if (side_len > 1.0e-4f) {
                side = vec3_mul1(side, 1.0f / side_len);
                const vec3_t e0 = vec3_mul1(side, a.w);
                const vec3_t e1 = vec3_mul1(side, b.w);
                const uint32_t fill = (col & 0x00FFFFFFu) | (RIBBON_ALPHA << 24);
                immediate::triangle(q, vec3_sub(p0, e0), vec3_add(p0, e0), vec3_add(p1, e1), fill);
                immediate::triangle(q, vec3_sub(p0, e0), vec3_add(p1, e1), vec3_sub(p1, e1), fill);
            }
        }
        immediate::line(q, p0, p1, col);
    }

    // The route a probe of the traced radius can follow, and the throat that stops a larger one.
    //
    // A centre line plus the probe itself, placed along it. A clearance scaled ribbon was tried
    // first and is wrong here for a reason worth recording: the route is straightened, so its
    // anchors can be sixty voxels apart, and a ribbon interpolates the clearance between them.
    // Through a wide cavity that draws a few enormous flat billboards whose width was never measured
    // anywhere along them. Spheres read as what they are - where a probe of that size fits - and
    // every one is sized by the field at its own centre, which is why the route is resampled.
    //
    // The throat is the same thing at the one place it is tight: "a sphere this big fits exactly
    // here and no bigger one gets past" is the whole content of the critical radius.
    void draw_route_3d(immediate::Queue* q) const {
        const size_t n = md_array_size(route);
        const float wrap_x = 0.5f * grid.spacing.x * (float)grid.dim[0];
        const float wrap_y = 0.5f * grid.spacing.y * (float)grid.dim[1];

        // Minimum image, so a step across a periodic face counts as the half voxel it is rather than
        // the width of the box. The same test says whether the segment may be drawn at all.
        auto step = [&](const vec4_t& a, const vec4_t& b, bool* wrapped) {
            float dx = b.x - a.x;
            float dy = b.y - a.y;
            *wrapped = false;
            if (fabsf(dx) > wrap_x) { dx -= copysignf(2.0f * wrap_x, dx); *wrapped = true; }
            if (fabsf(dy) > wrap_y) { dy -= copysignf(2.0f * wrap_y, dy); *wrapped = true; }
            const float dz = b.z - a.z;
            return sqrtf(dx * dx + dy * dy + dz * dz);
        };

        for (size_t i = 1; i < n; ++i) {
            bool wrapped = false;
            step(route[i - 1], route[i], &wrapped);
            if (wrapped) continue;
            const vec3_t p0 = {route[i - 1].x, route[i - 1].y, route[i - 1].z};
            const vec3_t p1 = {route[i].x, route[i].y, route[i].z};
            immediate::line(q, p0, p1, ROUTE_LINE_COLOR);
        }

        if (n == 0) return;

        // Spheres are spaced at SPHERE_SPACING times the local radius, so consecutive ones overlap by
        // about half a radius: enough to read as one channel, not so much that they fuse into a
        // sausage. The requirement is taken as the largest radius seen since the last sphere, which
        // is what keeps the overlap bounded where the channel opens out.
        const float min_step = MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        auto walk = [&](float scale, auto&& emit) {
            size_t count = 1;
            emit(route[0]);
            float travelled = 0.0f;
            float need = scale * SPHERE_SPACING * MAX(route[0].w, min_step);
            for (size_t i = 1; i < n; ++i) {
                bool wrapped = false;
                travelled += step(route[i - 1], route[i], &wrapped);
                need = MAX(need, scale * SPHERE_SPACING * route[i].w);
                if (travelled >= need) {
                    emit(route[i]);
                    count += 1;
                    travelled = 0.0f;
                    need = scale * SPHERE_SPACING * MAX(route[i].w, min_step);
                }
            }
            return count;
        };

        // Count first: a long route through a tight channel would otherwise fill the overlay with
        // thousands of spheres, so the spacing widens rather than the frame rate collapsing.
        const size_t want = walk(1.0f, [](const vec4_t&) {});
        const float  scale = (want > MAX_ROUTE_SPHERES) ? (float)want / (float)MAX_ROUTE_SPHERES : 1.0f;
        const int    detail = (want > 64) ? 6 : 8;

        walk(scale, [&](const vec4_t& p) {
            immediate::sphere_wireframe(q, vec3_t{p.x, p.y, p.z}, p.w, ROUTE_SPHERE_COLOR, detail, detail + 2);
        });

        if (has_perc && perc.has_r_c) {
            const vec3_t c = {perc.throat[0], perc.throat[1], perc.throat[2]};
            immediate::sphere_wireframe(q, c, (float)perc.r_c, THROAT_COLOR, 12, 16);
        }
    }

    void draw_channel_path_3d(immediate::Queue* q, uint32_t node, vec3_t cam_axis) const {
        auto segment = [&](vec4_t a, vec4_t b, uint32_t col) {
            draw_clearance_segment(q, a, b, col, cam_axis);
        };

        uint32_t n = node;
        bool first = true;
        while (n != CHANNEL_INVALID_INDEX) {
            const channel_node_t& cn = channels.nodes[n];
            const uint32_t col = first ? 0xff3060ff : 0xffb0a060;   // ABGR: the branch, then its route down

            for (size_t i = 1; i < md_array_size(cn.path); ++i) {
                segment(cn.path[i - 1], cn.path[i], col);
            }
            if (cn.parent != CHANNEL_INVALID_INDEX && md_array_size(cn.path) > 0) {
                const channel_node_t& pn = channels.nodes[cn.parent];
                if (md_array_size(pn.path) > 0) {
                    segment(*md_array_last(cn.path), pn.path[0], col);
                }
            }
            first = false;
            n = cn.parent;
        }
    }

    void draw_channel_tree(ImVec2 size) {
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        ImGui::InvisibleButton("##channel_tree", size);
        const bool active = ImGui::IsItemHovered();
        const ImVec2 mouse = ImGui::GetIO().MousePos;

        dl->AddRectFilled(p0, ImVec2(p0.x + size.x, p0.y + size.y), IM_COL32(24, 24, 28, 255));

        if (num_slots <= 0.0f) {
            const char* msg = "No channel gets through at this probe radius";
            const ImVec2 ts = ImGui::CalcTextSize(msg);
            dl->AddText(ImVec2(p0.x + 0.5f * (size.x - ts.x), p0.y + 0.5f * (size.y - ts.y)), IM_COL32(150, 150, 150, 255), msg);
            hovered_node = CHANNEL_INVALID_INDEX;
            return;
        }

        const float z_lo = grid.origin.z;
        const float z_hi = grid.origin.z + grid.spacing.z * (float)grid.dim[2];
        const float span = MAX(1.0e-6f, z_hi - z_lo);

        auto to_y = [&](float z) { return p0.y + (z_hi - z) / span * size.y; };
        auto to_x = [&](float s) { return p0.x + (s + 0.5f) / MAX(1.0f, num_slots) * size.x; };

        float max_clear = 1.0e-6f;
        for (size_t i = 0; i < md_array_size(channels.nodes); ++i) {
            if (node_slot[i] >= 0.0f) max_clear = MAX(max_clear, channels.nodes[i].max_clearance);
        }

        uint32_t best = CHANNEL_INVALID_INDEX;
        float    best_dist = 6.0f;

        for (size_t i = 0; i < md_array_size(channels.nodes); ++i) {
            if (node_slot[i] < 0.0f) continue;
            const channel_node_t& n = channels.nodes[i];

            const float x  = to_x(node_slot[i]);
            const float y0 = to_y(n.z_top);
            const float y1 = to_y(n.z_bot);

            // Width carries the tightest point of the branch: a thin line is a branch a probe barely fits through
            const float t = CLAMP(n.min_clearance / max_clear, 0.0f, 1.0f);
            const float w = 1.0f + 4.0f * t;

            const bool is_hovered = (hovered_node == (uint32_t)i);
            const ImU32 col = is_hovered ? IM_COL32(255, 150, 60, 255)
                                         : IM_COL32(90 + (int)(140 * t), 150 + (int)(80 * t), 230, 255);

            dl->AddLine(ImVec2(x, y0), ImVec2(x, y1), col, is_hovered ? w + 2.0f : w);

            if (n.parent != CHANNEL_INVALID_INDEX && node_slot[n.parent] >= 0.0f) {
                const float px = to_x(node_slot[n.parent]);
                const float py = to_y(channels.nodes[n.parent].z_top);
                dl->AddLine(ImVec2(x, y1), ImVec2(px, py), col, 1.0f);
            }

            if (active) {
                const float dy = (mouse.y < MIN(y0, y1)) ? MIN(y0, y1) - mouse.y
                               : (mouse.y > MAX(y0, y1)) ? mouse.y - MAX(y0, y1) : 0.0f;
                const float d = sqrtf((mouse.x - x) * (mouse.x - x) + dy * dy);
                if (d < best_dist) {
                    best_dist = d;
                    best = (uint32_t)i;
                }
            }
        }

        if (active) hovered_node = best;

        if (hovered_node != CHANNEL_INVALID_INDEX && hovered_node < md_array_size(channels.nodes)) {
            const channel_node_t& n = channels.nodes[hovered_node];
            ImGui::SetTooltip("z %.1f to %.1f nm\nTightest point %.2f nm\nWidest point %.2f nm\n%u voxels",
                              n.z_bot / ANGSTROM_PER_NM, n.z_top / ANGSTROM_PER_NM,
                              n.min_clearance / ANGSTROM_PER_NM, n.max_clearance / ANGSTROM_PER_NM,
                              n.num_voxels);
        }
    }

    // Porosity and accessible volume. Both are V(R) = Vol[d > R] over a slab range: porosity is the
    // R = 0 end of the same curve, which is why they share a section and a plot rather than sitting
    // in two places able to disagree.
    void draw_porosity_section() {
        const void_profile_t p = make_profile();
        const double nm  = (double)ANGSTROM_PER_NM;
        const double nm3 = nm * nm * nm;

        ImGui::SeparatorText("Porosity and accessible volume");

        if (ImGui::BeginCombo("Report over", region_lbl[region])) {
            for (int i = 0; i < Region_Count; ++i) {
                if (ImGui::Selectable(region_lbl[i], region == i)) {
                    region = (Region)i;
                    update_region();
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("A film simulated with an open z axis sits in a box with vacuum above and below it, and a\n"
                              "porosity averaged over that box measures how much vacuum the box was given. The scalars\n"
                              "below are reported over one region; the profiles always show the whole grid.");

        if (region == Region_Film) {
            float frac_pct = film_frac * 100.0f;
            if (ImGui::SliderFloat("Film edge", &frac_pct, 5.0f, 95.0f, "%.0f%% of interior")) {
                film_frac = frac_pct * 0.01f;
                update_region();
            }
            ImGui::SetItemTooltip("The film ends where the solid fraction has fallen to this share of its interior value.\n"
                                  "50%% is the half density convention; for a symmetric surface that is the Gibbs\n"
                                  "dividing surface, which is what a thickness means when the surface is rough.");
            if (!has_film) {
                ImGui::TextDisabled("No solid anywhere in the grid - falling back to the whole box");
            }
        } else if (region == Region_Manual) {
            float lo = manual_z_lo / (float)nm;
            float hi = manual_z_hi / (float)nm;
            const float z0 = grid.origin.z / (float)nm;
            const float z1 = (grid.origin.z + grid.spacing.z * (float)grid.dim[2]) / (float)nm;
            if (ImGui::DragFloatRange2("z range (nm)", &lo, &hi, 0.05f * (z1 - z0), z0, z1, "%.2f", "%.2f", ImGuiSliderFlags_AlwaysClamp)) {
                manual_z_lo = lo * (float)nm;
                manual_z_hi = hi * (float)nm;
                update_region();
            }
        }

        if (has_film) {
            const double f_lo = void_profile_z_lo(&p, film_beg) / nm;
            const double f_hi = void_profile_z_hi(&p, film_end - 1) / nm;
            ImGui::TextDisabled("Film slab: %.2f to %.2f nm, %.2f nm thick, interior solid fraction %.3f",
                                f_lo, f_hi, f_hi - f_lo, film_solid);
        }

        const uint64_t n_region = void_profile_num_total(&p, region_beg, region_end);
        const double   v_region = (double)n_region * stats.voxel_volume;
        const double   phi      = void_profile_porosity(&p, region_beg, region_end);
        const double   v_acc    = void_profile_accessible_volume(&p, region_beg, region_end, (double)probe_radius);
        const double   phi_acc  = (v_region > 0.0) ? v_acc / v_region : 0.0;

        ImGui::Text("Region: z %.2f to %.2f nm, %u of %u slabs, %.4g nm^3",
                    void_profile_z_lo(&p, region_beg) / nm,
                    void_profile_z_hi(&p, region_end - 1) / nm,
                    region_end - region_beg, stats.num_slabs, v_region / nm3);

        ImGui::Text("Porosity: %.4f", phi);
        ImGui::SetItemTooltip("Void volume fraction of the region: the share of it a point sized probe can occupy.\n"
                              "The solid it is measured against is the bead plus its hydration layer, so this moves\n"
                              "with that layer like everything else here - sweep it rather than fitting it once.");
        ImGui::SameLine();
        ImGui::TextDisabled("(solid %.4f, void volume %.4g nm^3)", 1.0 - phi, phi * v_region / nm3);

        ImGui::Text("V(R) at R = %.2f nm: %.4g nm^3, V(R)/V = %.4f", probe_radius / nm, v_acc / nm3, phi_acc);
        ImGui::SetItemTooltip("Volume whose distance to the nearest bead surface exceeds R, i.e. where the CENTRE of a\n"
                              "probe of radius R fits. Not the volume the probe body fills, which is that set dilated\n"
                              "by R and is larger, and not a statement that the volume can be reached from outside -\n"
                              "a closed cavity counts here and is invisible to infiltration. The channel sweep below\n"
                              "is what answers reachability.");
        ImGui::SameLine();
        ImGui::TextDisabled("(%.1f%% of the pore space)", (phi > 0.0) ? 100.0 * phi_acc / phi : 0.0);

        ImGui::SliderInt("Radius samples", &acc_num_radii, 8, 256, "%d", ImGuiSliderFlags_AlwaysClamp);
        if (ImGui::IsItemDeactivatedAfterEdit()) curves_dirty = true;

        // V(R)/V against R. Porosity is its left endpoint and the probe radius reads off it, so the
        // two scalars above are two points on this curve and the curve is the actual result.
        if (md_array_size(acc_r) > 1 && ImPlot::BeginPlot("##acc_curve", ImVec2(-1, 180))) {
            ImPlot::SetupAxes("Probe radius R (nm)", "V(R) / V");
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, MAX(1.0e-4, phi * 1.05), ImPlotCond_Always);
            ImPlot::PlotLine("V(R)/V", acc_r, acc_frac, (int)md_array_size(acc_r));
            double r_nm = probe_radius / nm;
            ImPlot::DragLineX(0, &r_nm, ImVec4(1, 0.6f, 0.2f, 1), 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }

        // The same against z. The R = 0 curve is the porosity profile, so the gap between the two
        // lines is the pore volume the probe is too large for, slab by slab.
        if (md_array_size(prof_z) > 1 && md_array_size(prof_acc) == md_array_size(prof_z) &&
            ImPlot::BeginPlot("##z_profile", ImVec2(-1, 180))) {
            ImPlot::SetupAxes("z (nm)", "Volume fraction");
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Once);
            ImPlot::PlotLine("Porosity", prof_z, prof_phi, (int)md_array_size(prof_z));
            ImPlot::PlotLine("V(R,z)/V(z)", prof_z, prof_acc, (int)md_array_size(prof_z));
            if (region_end > region_beg && (region_beg > 0 || region_end < stats.num_slabs)) {
                double rz[2] = { void_profile_z_lo(&p, region_beg) / nm, void_profile_z_hi(&p, region_end - 1) / nm };
                ImPlot::DragLineX(1, &rz[0], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                ImPlot::DragLineX(2, &rz[1], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
            }
            ImPlot::EndPlot();
        }

        // V(R,z): the accessible volume of every z slab over every probe radius, in one picture. A
        // cut at constant z is the curve above for that slab; a cut at constant R is the profile.
        // What the surface shows that neither cut does is where the two run out together - the
        // radius at which the interior of the film goes flat while the surface layers still have
        // room is the shape of an infiltration limit.
        if (ImGui::BeginCombo("V(R,z)", surface_view_lbl[surface_view])) {
            for (int i = 0; i < SurfaceView_Count; ++i) {
                if (ImGui::Selectable(surface_view_lbl[i], surface_view == i)) surface_view = (SurfaceView)i;
            }
            ImGui::EndCombo();
        }
        if (surface_view != SurfaceView_Hidden) {
            if (ImGui::BeginCombo("Height", surface_value_lbl[surface_value])) {
                for (int i = 0; i < SurfaceValue_Count; ++i) {
                    if (ImGui::Selectable(surface_value_lbl[i], surface_value == i)) {
                        surface_value = (SurfaceValue)i;
                        curves_dirty  = true;
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("The normalized form compares across slabs and across runs; the absolute one is what\n"
                                  "sums back to the accessible volume of the region. A slab thickness is a choice, so\n"
                                  "the absolute surface moves with the slab count and the normalized one does not.");
        }

        if (surface_view == SurfaceView_Surface && surf_nx > 1 && surf_ny > 1) {
            ImPlot3D::PushColormap(ImPlot3DColormap_Viridis);
            if (ImPlot3D::BeginPlot("##acc_surface", ImVec2(-1, 320))) {
                ImPlot3D::SetupAxes("z (nm)", "R (nm)", surface_value_lbl[surface_value]);
                ImPlot3D::SetupAxesLimits(surf_x0, surf_x1, surf_y0, surf_y1, 0.0, surf_max, ImPlot3DCond_Always);
                ImPlot3DSpec spec;
                spec.Marker    = ImPlot3DMarker_None;
                spec.FillAlpha = 1.0f;
                spec.LineColor = ImVec4(0.0f, 0.0f, 0.0f, 0.25f);
                spec.Flags     = ImPlot3DSurfaceFlags_NoMarkers;
                ImPlot3D::PlotSurface("V(R,z)", surf_x, surf_y, surf_z, (int)surf_nx, (int)surf_ny, 0.0, surf_max, spec);
                ImPlot3D::EndPlot();
            }
            ImPlot3D::PopColormap();
            ImGui::TextDisabled("Drag to rotate. %u slabs x %u radii%s",
                                surf_nx, surf_ny,
                                (surf_nx < stats.num_slabs || (int)surf_ny < acc_num_radii) ? ", sampled from a finer grid" : "");
        } else if (surface_view == SurfaceView_Heatmap && surf_nx > 1 && surf_ny > 1) {
            ImPlot::PushColormap(ImPlotColormap_Viridis);
            if (ImPlot::BeginPlot("##acc_heat", ImVec2(-76, 220), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
                ImPlot::SetupAxes("z (nm)", "Probe radius R (nm)");
                ImPlot::SetupAxesLimits(surf_x0, surf_x1, surf_y0, surf_y1, ImPlotCond_Always);
                ImPlot::PlotHeatmap("V(R,z)", acc_heat, (int)surf_ny, (int)surf_nx, 0.0, surf_max, nullptr,
                                    ImPlotPoint(surf_x0, surf_y0), ImPlotPoint(surf_x1, surf_y1));
                ImPlot::EndPlot();
            }
            ImGui::SameLine();
            ImPlot::ColormapScale("##acc_heat_scale", 0.0, surf_max, ImVec2(60, 220), "%.3g");
            ImPlot::PopColormap();
        }

        // Distance histogram of the region. V(R)/V above is its reverse cumulative, so this is the
        // same data read as a density: where the pore space sits in clearance, rather than how much
        // of it survives a given probe.
        if (md_array_size(hist_x) == NUM_BINS && ImPlot::BeginPlot("##Void distance", ImVec2(-1, 160))) {
            ImPlot::SetupAxes("Distance to nearest surface (nm)", "Volume fraction");
            ImPlot::PlotBars("Void", hist_x, hist_y, NUM_BINS, stats.bin_width / nm);
            ImPlot::EndPlot();
        }

        if (ImGui::Button("Export z profile")) {
            export_profile_csv();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export V(R,z)")) {
            export_accessible_csv();
        }
        ImGui::SetItemTooltip("Long form: one row per radius and slab, plus the region as a whole at slab -1.");
    }

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({420, 520}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Void Analysis", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            const bool has_system = app_state && app_state->mold.state.num_atoms > 0;

            ImGui::SeparatorText("Field");

            float spacing_nm = voxel_spacing / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Voxel spacing (nm)", &spacing_nm, 0.1f, 2.0f, "%.2f")) {
                voxel_spacing = spacing_nm * ANGSTROM_PER_NM;
            }

            float max_dist_nm = max_dist / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Max distance (nm)", &max_dist_nm, 1.0f, 100.0f, "%.1f")) {
                max_dist = max_dist_nm * ANGSTROM_PER_NM;
            }
            ImGui::SetItemTooltip("Distances are resolved out to this range; beyond it a voxel is reported at the limit.\n"
                                  "It is also the span of the distance binning, so a range far beyond the largest pore\n"
                                  "spends bins on nothing.");

            ImGui::SliderInt("z profile slabs", &z_slabs, 8, MAX_SLABS);
            ImGui::SetItemTooltip("How finely the profiles resolve z. Capped at one slab per plane of voxels -\n"
                                  "a slab thinner than that is empty or a copy of its neighbour.");

            float cell_ext_nm = cell_ext / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Acc cell extent (nm)", &cell_ext_nm, 0.3f, 10.0f, "%.2f")) {
                cell_ext = cell_ext_nm * ANGSTROM_PER_NM;
            }
            ImGui::SetItemTooltip("Cell size of the spatial acceleration structure. Smaller cells search less but cost memory:\n"
                                  "the structure allocates one offset per cell, and it is capped at 1024 cells per axis.");

            ImGui::SeparatorText("Radii");

            if (ImGui::BeginCombo("Source", radius_source_lbl[radius_source])) {
                for (int i = 0; i < RadiusSource_Count; ++i) {
                    if (ImGui::Selectable(radius_source_lbl[i], radius_source == i)) {
                        radius_source = (RadiusSource)i;
                    }
                }
                ImGui::EndCombo();
            }
            if (radius_source == RadiusSource_Uniform) {
                float r_nm = uniform_radius / ANGSTROM_PER_NM;
                if (ImGui::SliderFloat("Bead radius (nm)", &r_nm, 0.05f, 5.0f, "%.2f")) {
                    uniform_radius = r_nm * ANGSTROM_PER_NM;
                }
            }

            float t_layer_nm = t_layer / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Hydration layer (nm)", &t_layer_nm, 0.0f, 1.0f, "%.2f")) {
                t_layer = t_layer_nm * ANGSTROM_PER_NM;
            }
            ImGui::SetItemTooltip("Added to every bead radius. Sweep it rather than fitting it once:\n"
                                  "in a network with 2 nm gaps, 0.3 nm per side removes a third of the gap.");

            float probe_nm = probe_radius / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Probe radius (nm)", &probe_nm, 0.0f, 10.0f, "%.2f")) {
                probe_radius = probe_nm * ANGSTROM_PER_NM;
                probe_dirty  = true;
            }
            ImGui::SetItemTooltip("Water sized (0.14 nm) for sorption, colloid sized for exclusion.");

            ImGui::SeparatorText("Output");

            ImGui::Checkbox("Materialize field", &materialize);
            ImGui::SetItemTooltip("Keep the full voxel grid in memory. Statistics alone do not need it.");

            if (has_system) {
                md_grid_t preview = {};
                if (setup_grid(&preview, app_state->mold.state)) {
                    const size_t num_voxels = md_grid_num_points(&preview);
                    ImGui::Text("Grid: %i x %i x %i  (%.3g voxels)", preview.dim[0], preview.dim[1], preview.dim[2], (double)num_voxels);
                    if (materialize) {
                        ImGui::Text("Field: %.2f GB", (double)(num_voxels * sizeof(float)) / (double)GIGABYTES(1));
                    }
                }
            }

            ImGui::BeginDisabled(!has_system);
            if (ImGui::Button("Compute")) {
                compute();
            }
            ImGui::EndDisabled();

            if (error[0]) {
                ImGui::TextColored({1.0f, 0.4f, 0.4f, 1.0f}, "%s", error);
            }

            if (has_result) {
                if (curves_dirty) update_curves();
                if (probe_dirty)  update_probe_profile();

                ImGui::SeparatorText("Result");
                ImGui::Text("Computed in %.2f s over %.4g voxels", stats.seconds, (double)stats.num_voxels);
                if (stats.num_voxels < stats.num_grid) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("(%.4g outside the cell, skipped)", (double)(stats.num_grid - stats.num_voxels));
                    ImGui::SetItemTooltip("The grid covers the bounding box of a triclinic cell. Voxels outside the cell\n"
                                          "itself are periodic images of voxels already counted and are left out.");
                }
                ImGui::Text("Distance range: %.2f to %.2f nm", stats.d_min / ANGSTROM_PER_NM, stats.d_max / ANGSTROM_PER_NM);
                if (stats.num_clamped > 0) {
                    ImGui::TextColored({1.0f, 0.8f, 0.35f, 1.0f}, "%.2f%% of voxels have no bead within %.1f nm",
                                       100.0 * (double)stats.num_clamped / (double)MAX((uint64_t)1, stats.num_voxels),
                                       max_dist / ANGSTROM_PER_NM);
                    ImGui::SetItemTooltip("Their distance is that limit, not a measurement. Porosity is unaffected - they\n"
                                          "are void either way - but an accessible volume asked for at a radius near or\n"
                                          "above the limit reads them as fitting the probe. Raise the max distance.");
                }

                draw_porosity_section();
            }

            ImGui::SeparatorText("Percolation through z");

            ImGui::TextDisabled("Along z, which is treated as open at both ends.");

            float r_min_nm = channel_r_min / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Smallest radius (nm)", &r_min_nm, 0.05f, 5.0f, "%.2f")) {
                channel_r_min = r_min_nm * ANGSTROM_PER_NM;
            }
            ImGui::SetItemTooltip("The lowest radius the curves reach, and what bounds the memory: a voxel with\n"
                                  "less clearance than this never enters the structure. Below it a network is\n"
                                  "usually one connected pore and there is nothing left to resolve.");

            ImGui::SliderInt("Radius samples", &channel_num_radii, 8, 512);

            ImGui::BeginDisabled(!has_result || !field);
            if (ImGui::Button("Analyse percolation")) {
                compute_percolation();
            }
            ImGui::EndDisabled();
            if (!field) {
                ImGui::SameLine();
                ImGui::TextDisabled("(needs a materialized field)");
            }

            if (has_perc) {
                const double nm  = (double)ANGSTROM_PER_NM;
                const size_t np  = md_array_size(pc_r);
                const double r_p = probe_radius / nm;

                ImGui::Text("Computed in %.2f s over %.3g M active voxels, %.2f GB",
                            perc_seconds, (double)perc.num_active / 1.0e6,
                            (double)perc.bytes / (double)GIGABYTES(1));

                if (perc.has_r_c) {
                    ImGui::Text("Critical radius r_c: %.3f nm", perc.r_c / nm);
                    ImGui::SetItemTooltip("The largest probe which still gets from one z face to the other. It is the\n"
                                          "tightest throat along the widest route, not a cavity anywhere on it, and it\n"
                                          "is read off the insertion order rather than bracketed by bisection.");

                    ImGui::Text("Limiting throat at %.1f, %.1f, %.1f nm",
                                (double)perc.throat[0] / nm, (double)perc.throat[1] / nm, (double)perc.throat[2] / nm);
                    ImGui::SetItemTooltip("The one place a probe of r_c only just fits: the voxel whose insertion joined\n"
                                          "the two faces. Drawn in 3D as a sphere of that radius, which is the probe.");

                    if (has_route && md_array_size(route) > 1) {
                        const vec4_t a = route[0];
                        const vec4_t b = *md_array_last(route);
                        const double straight = MAX(1.0e-6, fabs((double)b.z - (double)a.z));
                        ImGui::Text("Route at %.2f nm: %.1f nm long, tortuosity %.2f",
                                    route_radius / nm, route_length / nm, route_length / straight);
                        ImGui::SetItemTooltip("A route a probe of that radius can actually follow, not a representative\n"
                                              "centreline: every point on it clears the probe and each is reachable from\n"
                                              "the last. Tortuosity is its length over the straight-line depth.");
                    }
                } else {
                    ImGui::Text("Nothing gets through at %.2f nm or above", channel_r_min / nm);
                    ImGui::SetItemTooltip("Either the network is closed to a probe this size, or the smallest radius is\n"
                                          "set above r_c. Lower it to find out which.");
                }

                ImGui::Checkbox("Show route and throat in 3D", &show_route);
                ImGui::SameLine();
                ImGui::BeginDisabled(!perc.has_r_c || (double)probe_radius > perc.r_c);
                if (ImGui::Button("Trace at the probe radius")) {
                    trace_route((double)probe_radius);
                }
                ImGui::EndDisabled();
                ImGui::SetItemTooltip("Retrace the route for the current probe rather than for r_c. Only possible below\n"
                                      "r_c, since above it there is no route to trace.");

                // How far in a probe gets, from each face. The two fronts approach each other as the
                // probe shrinks and meet exactly at r_c, so this is the percolation threshold arrived
                // at from the other direction - and unlike the single number, the band between them
                // says how much of the film a probe too large to cross can still get into.
                if (np > 1 && ImPlot::BeginPlot("##penetration", ImVec2(-1, 200))) {
                    ImPlot::SetupAxes("Probe radius R (nm)", "z (nm)");
                    ImPlot::PlotShaded("Reached by neither", pc_r, pc_bot, pc_top, (int)np);
                    ImPlot::PlotLine("From the top",    pc_r, pc_top, (int)np);
                    ImPlot::PlotLine("From the bottom", pc_r, pc_bot, (int)np);
                    if (has_film) {
                        const void_profile_t vp = make_profile();
                        double fz[2] = { void_profile_z_lo(&vp, film_beg) / nm, void_profile_z_hi(&vp, film_end - 1) / nm };
                        ImPlot::DragLineY(10, &fz[0], ImVec4(0.6f, 0.6f, 0.6f, 0.6f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImPlot::DragLineY(11, &fz[1], ImVec4(0.6f, 0.6f, 0.6f, 0.6f), 1.0f, ImPlotDragToolFlags_NoInputs);
                    }
                    if (perc.has_r_c) {
                        double rc_nm = perc.r_c / nm;
                        ImPlot::DragLineX(12, &rc_nm, ImVec4(1, 0.6f, 0.2f, 1), 1.0f, ImPlotDragToolFlags_NoInputs);
                    }
                    ImPlot::EndPlot();
                }
                ImGui::TextDisabled("The band is what a probe of that radius cannot enter from either side.");

                // Fits, reachable, spanning. The gap between the first two is the closed porosity -
                // pore space large enough for the probe and sealed off from it - which is the
                // exterior flood fill, arriving out of the same ordering rather than a second pass.
                if (np > 1 && ImPlot::BeginPlot("##perc_fractions", ImVec2(-1, 190))) {
                    ImPlot::SetupAxes("Probe radius R (nm)", "Volume fraction");
                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, MAX(1.0e-4, pc_void[0] * 1.05), ImPlotCond_Always);
                    ImPlot::PlotShaded("Closed", pc_r, pc_open, pc_void, (int)np);
                    ImPlot::PlotLine("Fits",      pc_r, pc_void, (int)np);
                    ImPlot::PlotLine("Reachable", pc_r, pc_open, (int)np);
                    ImPlot::PlotLine("Spanning",  pc_r, pc_span, (int)np);
                    if (perc.has_r_c) {
                        double rc_nm = perc.r_c / nm;
                        ImPlot::DragLineX(13, &rc_nm, ImVec4(1, 0.6f, 0.2f, 1), 1.0f, ImPlotDragToolFlags_NoInputs);
                    }
                    ImPlot::EndPlot();
                }

                {
                    const double v_fit  = perc_at(pc_void, r_p);
                    const double v_open = perc_at(pc_open, r_p);
                    const double v_span = perc_at(pc_span, r_p);
                    ImGui::Text("At probe %.2f nm - fits %.4f, reachable %.4f, spanning %.4f", r_p, v_fit, v_open, v_span);
                    ImGui::Text("Closed to it: %.4f of the box, %.1f%% of what it would otherwise fit",
                                v_fit - v_open, (v_fit > 0.0) ? 100.0 * (1.0 - v_open / v_fit) : 0.0);
                    ImGui::SetItemTooltip("Pore space wide enough for the probe with no route to the outside. The\n"
                                          "accessible volume above counts it; nothing that has to get there can.");
                    const size_t pi = perc_index_at(r_p);
                    ImGui::Text("Components at that radius: %u, of which %u get through",
                                perc.num_components[pi], perc.num_spanning[pi]);
                    ImGui::SetItemTooltip("Counting components is what this section used to report on its own. Below r_c\n"
                                          "a percolating network is one component and the count is 1, which is why the\n"
                                          "curves above are here instead.");
                }
            }

            if (ImGui::CollapsingHeader("Branch tree at the probe radius")) {
                ImGui::TextDisabled("One component with many branches on a real network - a shape, not a number.");

                ImGui::BeginDisabled(!has_result || !field);
                if (ImGui::Button("Build tree")) {
                    compute_channels();
                }
                ImGui::EndDisabled();

                if (has_channels) {
                    ImGui::Text("Built in %.2f s. %u components, %u of them through.",
                                channel_seconds, channels.num_components, channels.num_spanning);
                    ImGui::TextDisabled("Hover a branch to trace it in 3D");
                    draw_channel_tree(ImVec2(ImGui::GetContentRegionAvail().x, 220.0f));
                }
            }

            ImGui::SeparatorText("Accessible surface area");

            ImGui::SliderInt("Points per bead", &sasa_points, 32, SASA_MAX_POINTS);
            ImGui::SetItemTooltip("Quadrature over each bead's sphere. The error falls as 1/sqrt(N):\n"
                                  "256 points lands within a few tenths of a percent on an analytic case.");

            ImGui::SliderFloat("Bead sample fraction", &sasa_fraction, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Area is an average over beads, so sampling a fraction of them is unbiased.\n"
                                  "The reported uncertainty is the spread across the beads actually sampled.");

            ImGui::BeginDisabled(!has_system);
            if (ImGui::Button("Compute surface area")) {
                compute_sasa();
            }
            ImGui::EndDisabled();

            if (has_sasa) {
                const double area_nm2 = sasa.total_area / (double)(ANGSTROM_PER_NM * ANGSTROM_PER_NM);
                const double err_nm2  = sasa.std_err    / (double)(ANGSTROM_PER_NM * ANGSTROM_PER_NM);

                ImGui::Text("Computed in %.2f s from %zu beads", sasa.seconds, sasa.num_sampled);
                if (sasa.std_err > 0.0) {
                    ImGui::Text("Area: %.4g +/- %.2g nm^2", area_nm2, err_nm2);
                } else {
                    ImGui::Text("Area: %.4g nm^2", area_nm2);
                }
                if (sasa.area_density > 0.0) {
                    ImGui::Text("Area per volume: %.4g nm^-1", sasa.area_density * (double)ANGSTROM_PER_NM);
                }
                if (sasa.specific_area > 0.0) {
                    ImGui::Text("Specific area: %.4g m^2/g", sasa.specific_area);
                }
                ImGui::TextDisabled("At probe %.2f nm, hydration layer %.2f nm",
                                    sasa.probe / (double)ANGSTROM_PER_NM, sasa.t_layer / (double)ANGSTROM_PER_NM);

                // Independent cross check. The two estimators share almost no code path, so agreement within a few
                // percent is evidence; a larger gap means one of them is wrong, and the coarea one is the discretized
                // of the two.
                if (has_result) {
                    const double coarea = coarea_area(probe_radius) / (double)(ANGSTROM_PER_NM * ANGSTROM_PER_NM);
                    if (coarea > 0.0) {
                        ImGui::Text("Coarea cross check: %.4g nm^2  (%+.1f%%)", coarea,
                                    100.0 * (coarea / MAX(area_nm2, 1.0e-12) - 1.0));
                        ImGui::SetItemTooltip("Read off the distance histogram of the field, if one has been computed.\n"
                                              "Discretized, so a few percent is expected.");
                    }
                }

                // Per bead type. This is the differential hydration question: if the trajectory does not preserve
                // face identity, every bead lands in one bucket and that shows up here first.
                const md_atom_type_data_t& type = app_state->mold.sys.atom.type;
                if (md_array_size(sasa.type_area) > 1 && type.name) {
                    if (ImGui::BeginTable("##type_area", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchProp)) {
                        ImGui::TableSetupColumn("Type");
                        ImGui::TableSetupColumn("Area (nm^2)");
                        ImGui::TableSetupColumn("Share");
                        ImGui::TableHeadersRow();
                        for (size_t k = 0; k < md_array_size(sasa.type_area); ++k) {
                            if (sasa.type_area[k] <= 0.0) continue;
                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);
                            ImGui::Text("%.*s", (int)type.name[k].len, type.name[k].buf);
                            ImGui::TableSetColumnIndex(1);
                            ImGui::Text("%.4g", sasa.type_area[k] / (double)(ANGSTROM_PER_NM * ANGSTROM_PER_NM));
                            ImGui::TableSetColumnIndex(2);
                            ImGui::Text("%.1f%%", 100.0 * sasa.type_area[k] / MAX(sasa.total_area, 1.0e-12));
                        }
                        ImGui::EndTable();
                    }
                }
            }
        }
        ImGui::End();
    }
};

static VoidAnalysis instance = {};
