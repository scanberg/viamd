#define IMGUI_DEFINE_MATH_OPERATORS

#include <core/md_log.h>
#include <core/md_allocator.h>
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

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>

#include "void_analysis_core.h"

#include <core/md_str_builder.h>

#include <float.h>
#include <atomic>
#include <memory>

/*
    Void and nanopore characterization: the component.

    The computation lives in void_analysis_core.h - the distance field, the profile it is summarized into and the
    reductions of it, the surface topography, and the pore network - and is tested on its own under tests/.
    This file owns what needs an application: the window and its parameters, reading the system, spreading the
    field and topography passes over the task system, the volume rendering of the field, and the pore network in
    the 3D view, with its picking and tooltips.

    Residency: statistics are accumulated per tile while the tile is still in cache, so the full field is only
    materialized when explicitly asked for. At 0.5 nm voxels a 670 x 670 x 175 nm box is 2.4e9 voxels, which is
    9.6 GB as float - streaming is the default for a reason.
*/

namespace {

constexpr int   NUM_BINS  = 512;                      // Distance bins of the profile histogram
constexpr int   MAX_SLABS = 512;                      // Above this many voxel planes, slabs take several
constexpr int   SASA_MAX_POINTS = 1024;
constexpr float ANGSTROM_PER_NM = 10.0f;
constexpr double PI_D = 3.14159265358979323846;

// Angstrom^2 per Dalton to m^2 per gram
constexpr double SPECIFIC_AREA_SCALE = 1.0e-20 / 1.66053906892e-24;

// Dalton per Angstrom^3 to kg/m^3: 1.66053906892e-27 kg over 1e-30 m^3
constexpr double DA_PER_A3_TO_KG_PER_M3 = 1660.53906892;

// Largest 3D texture edge the field is uploaded at for volume rendering. A larger field is averaged
// down by a whole number of voxels per axis, which keeps the texture within what any GL 4 driver
// accepts and within a few hundred megabytes of video memory.
constexpr int VOL_MAX_DIM = 512;
constexpr int VOL_TF_RES  = 1024;    // Fine enough that the cut at R is sharp to a small fraction of a voxel

// Topography display: past this many columns per axis the heat map is averaged down. The statistics
// are always taken over every column.
constexpr uint32_t HEIGHT_MAX_DIM = 256;

// Pore network. Colours are ABGR, one per PORE_CLASS_, from the Okabe-Ito set: distinct under
// the common colour vision deficiencies, and dark enough to read against a white background as
// well as a black one. Vermillion for what gets through the film, blue and green for what is
// reachable from one side only, purple for what fits the probe and leads nowhere, grey for what
// the probe does not fit.
constexpr uint32_t PORE_CLASS_COLOR[5] = {
    0xff8c8c8c,     // PORE_CLASS_SMALL     #8C8C8C
    0xff7733aa,     // PORE_CLASS_CLOSED    #AA3377
    0xffb27200,     // PORE_CLASS_TOP       #0072B2
    0xff739e00,     // PORE_CLASS_BOTTOM    #009E73
    0xff005ed5,     // PORE_CLASS_SPANNING  #D55E00
};
const char* pore_class_lbl[5] = {
    "Too narrow",
    "Closed",
    "Reachable from the top",
    "Reachable from the bottom",
    "Through the film",
};
constexpr uint32_t PORE_ROUTE_COLOR         = 0xff003d9e;     // #9E3D00, a darker vermillion: the route is spanning by definition
constexpr uint32_t PORE_FOCUS_ON_LIGHT      = 0xff1a1a1a;
constexpr uint32_t PORE_FOCUS_ON_DARK       = 0xfff2f2f2;
constexpr float    PORE_POINT_SIZE          = 7.0f;     // Pixels, the smallest a pore is drawn
constexpr float    PORE_POINT_SIZE_FOCUS    = 12.0f;

// Past this many selected pores only their points are drawn in focus. A box over a dense network can
// take thousands, and each one in full focus is its sphere and a sphere per throat.
constexpr uint32_t PORE_MAX_FOCUS_DRAWN     = 2000;

// The skeleton is rebuilt every frame, so a network past this many throats draws the widest of them
// - which, since the list runs widest first, is a prefix - and says so in the window.
constexpr size_t   PORE_MAX_DRAWN_THROATS   = 250000;

constexpr PickingDomainID PickingDomain_PoreNetwork = HASH_STR_LIT64("picking domain void pore network");

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
// much vacuum the box was given rather than anything about the film. The film slab is therefore the
// default. On a box filled to its faces it covers the whole box unless an edge slab falls below the
// threshold, and the region line says which slabs were used. The profiles always show the whole grid.
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

// Which face of the film the topography map shows. The thickness is the one of the three which does
// not depend on where the film happens to sit in the box.
enum HeightView {
    HeightView_Top = 0,
    HeightView_Bottom,
    HeightView_Thickness,
    HeightView_Count,
};

const char* height_view_lbl[HeightView_Count] = {
    "Top surface",
    "Bottom surface",
    "Thickness",
};

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
    md_array(double)   slab_mass  = 0;   // [slab], Dalton, binned on the same slab boundaries
    md_array(double)   slab_count = 0;   // [slab], atoms
    double   mass_total    = 0.0;        // Of the atoms that landed in a slab
    double   count_total   = 0.0;
    uint32_t num_slabs     = 0;

    double   bin_width     = 0.0;
    double   z_min         = 0.0;
    double   z_max         = 0.0;
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
    size_t num_sampled   = 0;
    md_array(double) type_area = 0;  // Per atom type, extrapolated
    double seconds       = 0.0;
};

}  // namespace

// --- Background passes -----------------------------------------------------------------------------
//
// The field and the topography run on the task pool, never on the render thread. A pass owns
// everything it reads - a copy of the beads as they were when it was asked for - so playback, a new
// frame or a freed system cannot change what a running pass sees, and it owns everything it writes
// until the main thread installs it at the end. A pass that is cancelled, or superseded by a newer
// one, is told to stop and frees itself when its last task has run.

// The cutoff the acceleration structure is built for. Every query here is a box which grows as far
// as it has to, so this only sets the scale of the cells; the field pass costs the same within 10%
// from 10 to 40 A.
constexpr double ACC_CUTOFF = 20.0;

enum PassState { Pass_Running = 0, Pass_Done, Pass_Failed };
enum PassPhase { Phase_Setup = 0, Phase_Work, Phase_Finish };

void stats_free(Stats& s) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_array_free(s.hist, heap);
    md_array_free(s.slab_solid, heap);
    md_array_free(s.slab_total, heap);
    md_array_free(s.slab_mass, heap);
    md_array_free(s.slab_count, heap);
    s = {};
}

// The beads as a pass sees them: a copy, and the structure built over it on a pool thread.
struct BeadInput {
    md_array(vec3_t)  xyz   = 0;
    md_array(float)   radii = 0;
    md_unitcell_t     cell  = {};
    md_coord_stream_t coords = {};
    md_spatial_acc_t  acc   = {};
    void_beads_t      beads = {};
    bool              built = false;

    void build() {
        coords = md_coord_stream_from_aos((const float*)xyz, sizeof(vec3_t), NULL, md_array_size(xyz));
        acc = {};
        acc.alloc = md_get_heap_allocator();
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &coords;
        desc.cutoff   = ACC_CUTOFF;
        desc.unitcell = &cell;
        md_spatial_acc_init(&acc, &desc);
        beads = void_beads(&acc, radii, md_array_size(radii));
        built = true;
    }

    void release() {
        if (built) md_spatial_acc_free(&acc);
        built = false;
        beads = {};
        md_array_free(xyz,   md_get_heap_allocator());
        md_array_free(radii, md_get_heap_allocator());
        xyz = 0;
        radii = 0;
    }
};

struct Pass {
    std::atomic<bool> cancel{false};
    std::atomic<int>  state{Pass_Running};
    std::atomic<int>  phase{Phase_Setup};
    task_system::ID   task_work = task_system::INVALID_ID;
    std::atomic<uint32_t> work_done{0};         // Tiles or patches, counted one at a time: the pool
    uint32_t          work_total = 0;           // hands out ranges of hundreds, too coarse to report
    md_tick_t         t0 = 0;
    double            seconds = 0.0;
    char              error[256] = "";

    bool cancelled() const { return cancel.load(std::memory_order_relaxed); }

    // Called from a pool thread; the main thread only reads error once the pass is finished
    void fail(const char* msg) {
        snprintf(error, sizeof(error), "%s", msg);
        state  = Pass_Failed;
        cancel = true;
    }

    // Setup builds the structure and is not divisible; the work is the bulk and reports per range
    float progress() const {
        switch (phase.load()) {
        case Phase_Setup: return 0.0f;
        case Phase_Work:  return work_total ? (float)work_done.load(std::memory_order_relaxed) / (float)work_total : 0.0f;
        default:          return 1.0f;
        }
    }
};

struct FieldPass : Pass {
    BeadInput in;
    md_array(float) masses = 0;
    md_grid_t grid = {};
    double    max_dist = 0.0;
    uint32_t  planes_per_slab = 1;
    uint32_t  num_slabs = 1;
    bool      materialize = false;
    bool      pbc_z = false;
    size_t    num_threads = 1;

    void_field_desc_t desc = {};
    md_array(void_field_accum_t) accum = 0;     // One per thread, merged in the finish
    md_array(uint64_t) accum_mem = 0;

    md_array(float) field = 0;                  // Handed over on install
    Stats stats = {};

    ~FieldPass() {
        md_allocator_i* heap = md_get_heap_allocator();
        in.release();
        md_array_free(masses, heap);
        md_array_free(accum, heap);
        md_array_free(accum_mem, heap);
        md_array_free(field, heap);
        stats_free(stats);
    }
};

// The materialized field, shared between the window and any pass reading it. Replacing or clearing
// the field drops the window's reference only; the memory goes when the last reader is done.
struct FieldBuffer {
    md_array(float) data = 0;
    ~FieldBuffer() { md_array_free(data, md_get_heap_allocator()); }
};

struct NetworkPass : Pass {
    std::shared_ptr<FieldBuffer> keep;          // Alive for as long as the build reads it
    channel_field_t f = {};
    double   r_min = 0.0;
    double   merge = 0.0;
    uint64_t result_gen = 0;                    // The field it was built from

    pore_network_t net = {};                    // Handed over on install
    md_array(uint32_t) route = 0;
    double   bottleneck = 0.0;
    double   route_length = 0.0;

    ~NetworkPass() {
        pore_network_free(&net);
        md_array_free(route, md_get_heap_allocator());
    }
};

// The build reports from inside its sweep, which is also where it can be told to stop
bool network_progress(float fraction, void* user) {
    NetworkPass* p = (NetworkPass*)user;
    p->work_done.store((uint32_t)(fraction * (float)p->work_total), std::memory_order_relaxed);
    return !p->cancelled();
}

struct HeightPass : Pass {
    BeadInput in;
    md_grid_t grid = {};
    double    probe = 0.0;
    double    max_dist = 0.0;
    uint64_t  result_gen = 0;                   // The field it was made for
    void_heightmap_desc_t desc = {};

    md_array(float) top = 0;                    // Handed over on install
    md_array(float) bot = 0;
    md_array(float) thk = 0;
    void_heightmap_stats_t st[HeightView_Count] = {};

    ~HeightPass() {
        md_allocator_i* heap = md_get_heap_allocator();
        in.release();
        md_array_free(top, heap);
        md_array_free(bot, heap);
        md_array_free(thk, heap);
    }
};

struct VoidAnalysis : viamd::EventHandler {
    bool show_window = false;

    // Parameters, all in Angstrom to match the rest of mdlib. The design document is written in nm, so the UI
    // presents nm and converts on the way in.
    float voxel_spacing = 5.0f;      // 0.5 nm
    float probe_radius  = 1.4f;      // Water sized by default
    float max_dist      = 320.0f;    // 32 nm, the range the field is resolved over
    float uniform_radius = 5.0f;

    RadiusSource radius_source = RadiusSource_Vdw;

    bool  materialize = false;
    md_array(float) field = 0;                  // The data of field_owner, or null
    std::shared_ptr<FieldBuffer> field_owner;
    md_grid_t grid = {};

    // Pore network: the clearance field as pores and the throats between them, drawn in 3D.
    float    net_r_min       = 2.0f;    // Smallest throat and pore; also what bounds the memory
    float    net_merge       = 0.5f;    // Persistence below which two pores are one, in voxel spacings
    bool     has_net         = false;
    bool     show_net        = true;
    bool     show_net_route  = true;
    bool     net_show_small  = false;   // Pores and throats too narrow for the probe
    double   net_seconds     = 0.0;
    pore_network_t net       = {};
    md_array(uint8_t)  net_class = 0;   // PORE_CLASS_ per pore at net_class_r
    float    net_class_r     = -1.0f;
    uint32_t net_class_count[5] = {};
    md_array(uint32_t) net_route = 0;   // Widest route, top face first
    double   net_route_bottleneck = 0.0;
    double   net_route_length     = 0.0;
    uint32_t net_hovered     = PORE_INVALID;
    PickingRange net_picking = {};
    size_t   net_drawn_throats = 0;
    md_array(uint8_t) net_selected = 0; // Per pore, set by clicking it in the 3D view
    uint32_t net_num_selected = 0;
    md_array(uint8_t) net_region = 0;   // Per pore, inside the box being dragged out right now
    uint32_t net_num_region = 0;
    bool     net_region_removing = false;

    // Porosity and accessible volume
    Region   region      = Region_Film;
    float    film_frac   = 0.5f;     // Solid fraction, as a share of the interior value, at the film edge
    float    manual_z_lo = 0.0f;
    float    manual_z_hi = 0.0f;
    uint32_t region_beg  = 0;        // Resolved slab range the scalars are reported over
    uint32_t region_end  = 0;
    uint32_t film_beg    = 0;
    uint32_t film_end    = 0;
    double   film_solid  = 0.0;
    bool     has_film    = false;

    // Cached reductions. The curves are rebuilt when the result or the region changes rather than
    // every frame; the per slab accessible fraction is one sweep of the bins and can follow the
    // probe slider.
    bool   curves_dirty  = true;
    bool   probe_dirty   = true;
    md_array(double) acc_r    = 0;   // Probe radius at the histogram bin edges, nm: the knots of V(R)
    md_array(double) acc_frac = 0;   // V(R)/V over the region
    md_array(double) prof_z   = 0;   // Slab centres, nm
    md_array(double) prof_phi = 0;   // Porosity per slab
    md_array(double) prof_acc = 0;   // V(R,z)/V(z) per slab at the current probe
    md_array(double) prof_rho = 0;   // Density per slab, in density_unit()
    double z_link[2] = {0.0, 1.0};   // Shared z axis of the two profiles against z, nm
    md_array(double) hist_x   = 0;   // Distance histogram of the region, nm
    md_array(double) hist_y   = 0;

    // Surface topography, as a probe of the current radius finds it from above and below. Computed
    // along with the field, and again when the probe radius is changed, since unlike the profile it
    // cannot be read off a histogram: it is a property of the geometry at one R.
    bool   has_height       = false;
    bool   height_pending   = false;    // The probe radius changed since the map was made
    double height_probe     = 0.0;      // R the map belongs to
    double height_seconds   = 0.0;
    HeightView height_view  = HeightView_Top;
    md_array(float) height_top  = 0;    // [x + y * dim_x], world z of the probe apex, NaN where open
    md_array(float) height_bot  = 0;
    md_array(float) height_thk  = 0;    // top - bottom
    void_heightmap_stats_t height_stats[HeightView_Count] = {};
    bool   height_disp_dirty = true;
    uint32_t height_nx = 0;             // Display grid, averaged down from the columns
    uint32_t height_ny = 0;
    md_array(float) height_disp = 0;    // [(ny - 1 - j) * nx + i], nm, row 0 is the largest y
    double height_lo = 0.0, height_hi = 1.0;   // Colour scale, nm

    // The distance field as a volume. The transfer function is transparent below the probe radius,
    // so what is drawn is exactly the set a probe centre of that size can occupy, coloured by how
    // much room it has there.
    bool     show_volume   = false;
    float    vol_opacity   = 0.1f;
    float    vol_upper     = 0.0f;      // Transparent above this clearance, Angstrom; 0 until a result sets it
    bool     vol_clip_region = true;    // Clip z to the reporting region
    int      vol_colormap  = ImPlotColormap_Viridis;
    bool     vol_dirty     = false;     // Field changed, texture needs uploading
    bool     vol_tf_dirty  = true;
    uint32_t vol_tex       = 0;
    uint32_t vol_tf_tex    = 0;
    int      vol_dim[3]    = {0, 0, 0};
    int      vol_stride    = 1;
    vec3_t   vol_spacing   = {};
    mat4_t   vol_model     = {};
    float    vol_z0        = 0.0f;      // World z extent the texture covers, for the region clip
    float    vol_z1        = 1.0f;

    // Surface area
    int   sasa_points   = 256;
    float sasa_fraction = 1.0f;      // Fraction of beads sampled; the estimate is extrapolated from them

    Stats stats = {};
    bool  has_result = false;
    SasaResult sasa = {};
    bool  has_sasa = false;
    char  error[256] = "";

    md_allocator_i* arena = nullptr;    // The heap; see EventType_ViamdInitialize

    FieldPass*  field_pass  = nullptr;  // Running, or null. Owned by its tasks, see install_field_pass
    HeightPass* height_pass = nullptr;
    NetworkPass* net_pass   = nullptr;
    uint64_t    result_gen  = 0;        // Bumped whenever a new field is installed
    ApplicationState* app_state = nullptr;

    VoidAnalysis() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                // Everything the window keeps is freed when it is replaced. An arena never frees, so
                // it would keep every field ever computed until shutdown.
                arena = md_get_heap_allocator();
                break;
            }
            case viamd::EventType_ViamdShutdown:
                // Running passes stop and free themselves; what the window holds is freed here.
                cancel_field_pass();
                cancel_height_pass();
                cancel_network_pass();
                clear_result();
                clear_sasa();
                free_volume();
                arena = nullptr;
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Void Analysis", &show_window);
                break;
            case viamd::EventType_ViamdSystemFree:
                cancel_field_pass();
                cancel_height_pass();
                clear_result();
                clear_sasa();
                break;
            case viamd::EventType_ViamdRenderOpaque: {
                if (e.payload_type != viamd::EventPayloadType_ApplicationState) break;
                // Queued into the world pass, which follows this event and depth tests. Like the
                // volume, an explicit toggle that stays up with the window closed.
                if (show_net && has_net) {
                    draw_network_3d(*(const ApplicationState*)e.payload);
                }
                break;
            }
            case viamd::EventType_ViamdRenderTransparent: {
                if (e.payload_type != viamd::EventPayloadType_ApplicationState) break;
                // The volume is an explicit toggle and stays up with the window closed, like any
                // other representation.
                if (show_volume && has_result) {
                    draw_volume(*(const ApplicationState*)e.payload);
                }
                break;
            }
            case viamd::EventType_ViamdPickingRangeReserve: {
                if (e.payload_type != viamd::EventPayloadType_PickingSpace) break;
                net_picking = {};
                if (show_net && has_net && md_array_size(net.vertices) > 0) {
                    if (!picking_range_reserve(&net_picking, (PickingSpace*)e.payload, PickingDomain_PoreNetwork, md_array_size(net.vertices))) {
                        net_picking = {};
                    }
                }
                break;
            }
            case viamd::EventType_ViamdInteractionSurface: {
                if (e.payload_type != viamd::EventPayloadType_InteractionSurfaceEvent) break;
                network_interaction(*(const InteractionSurfaceEvent*)e.payload);
                break;
            }
            case viamd::EventType_ViamdPickingTooltipTextRequest: {
                if (e.payload_type != viamd::EventPayloadType_PickingTooltipTextRequest) break;
                PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)e.payload;
                if (req->hit.domain == PickingDomain_PoreNetwork) network_tooltip(&req->sb, req->hit.local_idx);
                break;
            }
            default:
                break;
            }
        }
    }

    void clear_result() {
        field_owner.reset();
        field = 0;
        md_array_free(stats.hist, arena);
        md_array_free(stats.slab_solid, arena);
        md_array_free(stats.slab_total, arena);
        md_array_free(stats.slab_mass, arena);
        md_array_free(stats.slab_count, arena);
        stats = {};
        md_array_free(acc_r, arena);
        md_array_free(acc_frac, arena);
        md_array_free(prof_z, arena);
        md_array_free(prof_phi, arena);
        md_array_free(prof_acc, arena);
        md_array_free(prof_rho, arena);
        md_array_free(hist_x, arena);
        md_array_free(hist_y, arena);
        acc_r = 0; acc_frac = 0;
        prof_z = 0; prof_phi = 0; prof_acc = 0; prof_rho = 0;
        hist_x = 0; hist_y = 0;
        region_beg = region_end = 0;
        has_film = false;
        curves_dirty = true;
        probe_dirty  = true;
        has_result = false;
        clear_height();
        vol_dirty = false;
        clear_network();
    }

    void clear_height() {
        md_array_free(height_top, arena);
        md_array_free(height_bot, arena);
        md_array_free(height_thk, arena);
        md_array_free(height_disp, arena);
        height_top = 0; height_bot = 0; height_thk = 0; height_disp = 0;
        height_nx = height_ny = 0;
        MEMSET(height_stats, 0, sizeof(height_stats));
        has_height = false;
        height_pending = false;
        height_disp_dirty = true;
    }

    void free_volume() {
        if (vol_tex)    gl::free_texture(&vol_tex);
        if (vol_tf_tex) gl::free_texture(&vol_tf_tex);
        vol_tex = 0;
        vol_tf_tex = 0;
        vol_dim[0] = vol_dim[1] = vol_dim[2] = 0;
        vol_tf_dirty = true;
    }

    void clear_sasa() {
        md_array_free(sasa.type_area, arena);
        sasa = {};
        has_sasa = false;
    }

    // Radii the geometry is actually built from: the bead, and optionally a probe.
    // returns max radius, for the acceleration structure
    double fill_radii(float* out_radii, const md_system_t& sys, size_t count, double extra) const {
        double max = 0.0;
        if (radius_source == RadiusSource_Vdw) {
            md_atom_extract_radii(out_radii, 0, count, &sys.atom);
            for (size_t i = 0; i < count; ++i) {
                float rad = out_radii[i] + extra;
                max = MAX(max, rad);
                out_radii[i] = rad;
            }
        } else {
            float rad = uniform_radius + extra;
            max = rad;
            for (size_t i = 0; i < count; ++i) {
                out_radii[i] = rad;
            }
        }
        return max;
    }

    // Axis aligned grid covering the unit cell, or the atoms when there is no cell.
    // The voxel count is chosen so the grid tiles the box exactly, which keeps a periodic field seamless.
    bool setup_grid(md_grid_t* out_grid, const md_system_state_t& state) {
        return void_field_grid(out_grid, &state.unitcell, state.xyz, state.num_atoms, voxel_spacing);
    }

    // A copy of the beads, as every pass wants them. Main thread, since it reads the system.
    void snapshot_beads(BeadInput* in) const {
        const md_system_t&       sys   = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;
        md_allocator_i* heap = md_get_heap_allocator();
        md_array_resize(in->xyz,   state.num_atoms, heap);
        md_array_resize(in->radii, state.num_atoms, heap);
        MEMCPY(in->xyz, state.xyz, state.num_atoms * sizeof(vec3_t));
        // Bead radii. Every reported quantity is a function of them, so the radius source is a
        // choice worth stating next to any number quoted from here.
        fill_radii(in->radii, sys, state.num_atoms, 0.0);
        in->cell = state.unitcell;
    }

    bool height_running() const { return height_pass != nullptr; }

    // Abandon the running pass, if any. It stops at the next tile and frees itself.
    void cancel_field_pass() {
        if (!field_pass) return;
        field_pass->cancel = true;
        task_system::task_interrupt(field_pass->task_work);
        field_pass = nullptr;
    }

    void cancel_height_pass() {
        if (!height_pass) return;
        height_pass->cancel = true;
        task_system::task_interrupt(height_pass->task_work);
        height_pass = nullptr;
    }

    // Snapshot the system, then hand everything else to the pool: building the structure, the
    // field over the tiles, and the reduction of the per thread accumulators. The previous result
    // stays up until the new one is installed.
    void compute() {
        error[0] = '\0';

        const md_system_t&       sys   = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;

        if (state.num_atoms == 0) {
            snprintf(error, sizeof(error), "No system loaded");
            return;
        }

        md_grid_t g = {};
        if (!setup_grid(&g, state)) {
            snprintf(error, sizeof(error), "Could not derive a grid from the system");
            return;
        }

        const size_t num_voxels = md_grid_num_points(&g);
        if (materialize && num_voxels * sizeof(float) > GIGABYTES(4)) {
            snprintf(error, sizeof(error), "Materialized field would need %.1f GB, refusing", (double)(num_voxels * sizeof(float)) / (double)GIGABYTES(1));
            return;
        }

        cancel_field_pass();

        FieldPass* p = new FieldPass();
        p->t0 = md_tick_now();
        snapshot_beads(&p->in);
        md_array_resize(p->masses, state.num_atoms, md_get_heap_allocator());
        md_atom_extract_masses(p->masses, 0, state.num_atoms, &sys.atom);

        p->grid        = g;
        p->max_dist    = (double)max_dist;
        p->materialize = materialize;
        p->pbc_z       = (md_unitcell_flags(&state.unitcell) & MD_UNITCELL_PBC_Z) != 0;

        // z resolution of the profile: one plane of voxels per slab, which is the resolution the
        // field has. Whole planes only, see void_profile_planes_per_slab for why.
        p->planes_per_slab = void_profile_planes_per_slab(g.dim[2], MAX_SLABS);
        p->num_slabs       = void_profile_num_slabs(g.dim[2], p->planes_per_slab);

        // Thread 0 is the main thread, which does not take part here but is counted so that any
        // thread index the pool reports has an accumulator of its own.
        p->num_threads = MAX((size_t)1, task_system::pool_num_threads() + 1);

        const uint32_t num_tiles = void_field_num_tiles(&g);

        const task_system::ID setup = task_system::create_pool_task(STR_LIT("Void field setup"), [p, num_voxels]() {
            if (p->cancelled()) return;
            md_allocator_i* heap = md_get_heap_allocator();
            p->in.build();
            if (p->materialize) {
                md_array_resize(p->field, num_voxels, heap);
                if (!p->field) { p->fail("Could not allocate the materialized field"); return; }
            }

            const size_t hist_stride = (size_t)p->num_slabs * NUM_BINS;
            const size_t per_thread  = hist_stride + 2 * (size_t)p->num_slabs;
            md_array_resize(p->accum_mem, p->num_threads * per_thread, heap);
            md_array_resize(p->accum, p->num_threads, heap);
            for (size_t i = 0; i < p->num_threads; ++i) {
                uint64_t* base = p->accum_mem + i * per_thread;
                p->accum[i] = {};
                p->accum[i].hist  = base;
                p->accum[i].solid = base + hist_stride;
                p->accum[i].total = base + hist_stride + p->num_slabs;
                void_field_accum_reset(&p->accum[i], p->num_slabs, NUM_BINS);
            }

            // The bins are the resolution in R of every accessible volume reported later, so they
            // are deliberately finer than the voxel spacing: R is a continuous parameter and the
            // voxelization, not the binning, is what should be limiting.
            p->desc = {};
            p->desc.beads           = p->in.beads;
            p->desc.cell            = &p->in.cell;
            p->desc.grid            = &p->grid;
            p->desc.max_dist        = p->max_dist;
            p->desc.planes_per_slab = p->planes_per_slab;
            p->desc.num_bins        = NUM_BINS;
            p->desc.field           = p->field;
            p->phase = Phase_Work;
        });

        // Per thread accumulators: nothing is shared while the tiles run
        // One tile at a time within a range, so a cancel lands within a tile rather than a range
        p->work_total = num_tiles;
        p->task_work = task_system::create_pool_task(STR_LIT("Void distance field"), num_tiles, [p](uint32_t beg, uint32_t end, uint32_t thread_num) {
            const size_t ti = MIN(p->num_threads - 1, (size_t)thread_num);
            for (uint32_t t = beg; t < end && !p->cancelled(); ++t) {
                void_field_eval_tiles(&p->accum[ti], &p->desc, t, t + 1);
                p->work_done.fetch_add(1, std::memory_order_relaxed);
            }
        }, 1);

        const task_system::ID finish = task_system::create_pool_task(STR_LIT("Void field reduce"), [p, num_voxels]() {
            p->phase = Phase_Finish;
            if (!p->cancelled()) finish_field_pass(p, num_voxels);
            p->in.release();
            md_array_free(p->accum, md_get_heap_allocator());
            md_array_free(p->accum_mem, md_get_heap_allocator());
            p->accum = 0;
            p->accum_mem = 0;
            p->seconds = md_tick_to_seconds(md_tick_now() - p->t0);
            if (p->state == Pass_Running) p->state = p->cancelled() ? Pass_Failed : Pass_Done;
        });

        const task_system::ID install = task_system::create_main_task(STR_LIT("Void field install"), [this, p]() {
            install_field_pass(p);
        });

        task_system::set_task_dependency(p->task_work, setup);
        task_system::set_task_dependency(finish, p->task_work);
        task_system::set_task_dependency(install, finish);

        field_pass = p;
        task_system::enqueue_task(setup);
    }

    // Pool thread. Everything here reads and writes the pass only.
    static void finish_field_pass(FieldPass* p, size_t num_voxels) {
        md_allocator_i* heap = md_get_heap_allocator();
        Stats& s = p->stats;
        const uint32_t ns = p->num_slabs;

        md_array_resize(s.hist, (size_t)ns * NUM_BINS, heap);
        md_array_resize(s.slab_solid, ns, heap);
        md_array_resize(s.slab_total, ns, heap);

        void_field_accum_t merged = {};
        merged.hist  = s.hist;
        merged.solid = s.slab_solid;
        merged.total = s.slab_total;
        void_field_accum_reset(&merged, ns, NUM_BINS);
        for (size_t t = 0; t < p->num_threads; ++t) {
            void_field_accum_merge(&merged, &p->accum[t], ns, NUM_BINS);
        }

        const void_profile_t prof = void_field_profile(&merged, &p->desc);

        s.num_grid   = num_voxels;
        s.num_voxels = 0;
        s.num_solid  = 0;
        for (uint32_t sl = 0; sl < ns; ++sl) {
            s.num_voxels += merged.total[sl];
            s.num_solid  += merged.solid[sl];
        }
        s.num_clamped = merged.num_clamped;
        if (merged.d_min <= merged.d_max) {
            s.d_min = (double)merged.d_min;
            s.d_max = (double)merged.d_max;
        }

        s.num_slabs    = ns;
        s.bin_width    = prof.bin_width;
        s.z_min        = prof.z_min;
        s.z_max        = prof.z_max;
        s.slab_height  = prof.slab_height;
        s.voxel_volume = prof.voxel_volume;

        // Mass along z, on the same slabs as the histogram so that a density and a porosity quoted
        // for a slab range are about the same volume. z is wrapped only when the cell says it is
        // periodic; on an open axis an atom outside the grid is outside what is being measured.
        const size_t n = md_array_size(p->in.xyz);
        md_array_resize(s.slab_mass,  ns, heap);
        md_array_resize(s.slab_count, ns, heap);
        s.mass_total  = void_profile_bin_mass(s.slab_mass,  &prof, p->in.xyz, p->masses, n, p->pbc_z);
        s.count_total = void_profile_bin_mass(s.slab_count, &prof, p->in.xyz, NULL,      n, p->pbc_z);
    }

    // Main thread, when the pass's last task has run. Only the pass the window is still waiting for
    // is installed; anything else was cancelled or superseded and is just freed.
    void install_field_pass(FieldPass* p) {
        if (p != field_pass) {
            delete p;
            return;
        }
        field_pass = nullptr;

        if (p->state != Pass_Done) {
            if (p->error[0]) snprintf(error, sizeof(error), "%s", p->error);
            delete p;
            return;
        }

        // Once asked for, the topography follows the field it belongs to
        const bool want_height = has_height || height_pass != nullptr;
        cancel_height_pass();
        clear_result();

        grid  = p->grid;
        if (p->field) {
            field_owner = std::make_shared<FieldBuffer>();
            field_owner->data = p->field;
            p->field = 0;
            field = field_owner->data;
        }
        stats = p->stats;
        p->stats = {};
        stats.seconds = p->seconds;
        delete p;

        has_result  = true;
        result_gen += 1;

        // A new field needs uploading, and the colour range follows the new distances unless the
        // user had narrowed it to something the new result still covers.
        vol_dirty    = (field != nullptr);
        vol_tf_dirty = true;
        if (!(vol_upper > 0.0f) || vol_upper > (float)stats.d_max) vol_upper = (float)stats.d_max;

        z_link[0] = grid.origin.z / (double)ANGSTROM_PER_NM;
        z_link[1] = (grid.origin.z + grid.spacing.z * (float)grid.dim[2]) / (double)ANGSTROM_PER_NM;

        // A manual range the user has already set survives a recompute; an unset one starts as the
        // whole grid so dragging it narrows rather than starting from nothing.
        if (!(manual_z_hi > manual_z_lo)) {
            manual_z_lo = grid.origin.z;
            manual_z_hi = grid.origin.z + grid.spacing.z * (float)grid.dim[2];
        }
        update_region();
        if (want_height) start_height();
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
        p.z_max        = stats.z_max;
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

    // V(R)/V over the region, the porosity and density profiles against z, and the distance
    // histogram of the region. All of it is sums over bins, cheap enough to do on demand and too
    // much to do every frame.
    void update_curves() {
        curves_dirty = false;
        if (!has_result) return;

        const void_profile_t p = make_profile();
        const double nm = (double)ANGSTROM_PER_NM;

        const uint32_t ns = stats.num_slabs;
        md_array_resize(prof_z,   ns, arena);
        md_array_resize(prof_phi, ns, arena);
        for (uint32_t sl = 0; sl < ns; ++sl) {
            prof_z[sl]   = 0.5 * (void_profile_z_lo(&p, sl) + void_profile_z_hi(&p, sl)) / nm;
            prof_phi[sl] = void_profile_porosity(&p, sl, sl + 1);
        }

        md_array_resize(prof_rho, ns, arena);
        for (uint32_t sl = 0; sl < ns; ++sl) {
            prof_rho[sl] = density_over(sl, sl + 1);
        }

        // The distance histogram of the region, summed over its slabs once. The bars are this as it
        // stands and the V(R)/V curve is its reverse cumulative.
        uint64_t region_hist[NUM_BINS] = {};
        for (uint32_t sl = region_beg; sl < region_end; ++sl) {
            const uint64_t* h = stats.hist + (size_t)sl * NUM_BINS;
            for (int b = 0; b < NUM_BINS; ++b) region_hist[b] += h[b];
        }
        const double n_tot = (double)MAX((uint64_t)1, void_profile_num_total(&p, region_beg, region_end));

        md_array_resize(hist_x, NUM_BINS, arena);
        md_array_resize(hist_y, NUM_BINS, arena);
        for (int b = 0; b < NUM_BINS; ++b) {
            hist_x[b] = ((double)b + 0.5) * stats.bin_width / nm;
            hist_y[b] = (double)region_hist[b] / n_tot;
        }

        // V(R)/V. void_profile_count_above splits the bin holding R linearly, so the curve is
        // piecewise linear with its knots at the bin edges and is exact when evaluated there - there
        // is no sampling resolution to choose. From R = 0, where it is the porosity, to the first
        // edge past the largest clearance, where it reaches zero.
        const double   d_hi      = MAX(stats.d_max, 0.0);
        const uint32_t num_edges = (stats.bin_width > 0.0) ? (uint32_t)MIN((double)NUM_BINS, floor(d_hi / stats.bin_width) + 1.0) : 1;
        const uint32_t num_knots = num_edges + 1;
        md_array_resize(acc_r,    num_knots, arena);
        md_array_resize(acc_frac, num_knots, arena);
        uint64_t above = 0;
        for (int b = NUM_BINS; b >= 0; --b) {
            if (b < NUM_BINS) above += region_hist[b];
            if ((uint32_t)b < num_knots) {
                acc_r[b]    = (double)b * stats.bin_width / nm;
                acc_frac[b] = (double)above / n_tot;
            }
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

    // --- Density ------------------------------------------------------------------------------
    //
    // Mass over the volume the profile considered, on the same slabs as the porosity. A coarse
    // grained model without masses still has a meaningful number density, so that is what is shown
    // when there is no mass to speak of rather than a column of zeros.

    bool has_mass() const { return stats.mass_total > 0.0; }

    const char* density_unit() const { return has_mass() ? "kg/m^3" : "atoms/nm^3"; }

    double density_over(uint32_t slab_beg, uint32_t slab_end) const {
        if (!has_result) return 0.0;
        const void_profile_t p = make_profile();
        if (has_mass()) {
            return void_profile_mass_density(&p, stats.slab_mass, slab_beg, slab_end) * DA_PER_A3_TO_KG_PER_M3;
        }
        return void_profile_mass_density(&p, stats.slab_count, slab_beg, slab_end) * 1000.0;
    }

    // --- Surface topography -------------------------------------------------------------------

    // Heights of the probe apex over every column of the grid, from above and from below, at the
    // current probe radius. beads are the ones the field was built from, or built the same way.
    // The topography is asked for, not computed with every field: a pass of its own over every
    // column, which on a sparse network costs as much as the field again. Same shape as the field
    // pass - snapshot here, everything else on the pool.
    void start_height() {
        height_pending = false;
        if (!has_result || !app_state || app_state->mold.state.num_atoms == 0) return;

        const size_t n = (size_t)grid.dim[0] * (size_t)grid.dim[1];
        if (n == 0) return;
        if (!((double)probe_radius < (double)max_dist)) {
            snprintf(error, sizeof(error), "Probe radius must be below the max distance for the topography");
            return;
        }

        cancel_height_pass();

        HeightPass* p = new HeightPass();
        p->t0         = md_tick_now();
        p->grid       = grid;
        p->probe      = (double)probe_radius;
        p->max_dist   = (double)max_dist;
        p->result_gen = result_gen;
        snapshot_beads(&p->in);

        const task_system::ID setup = task_system::create_pool_task(STR_LIT("Void topography setup"), [p, n]() {
            if (p->cancelled()) return;
            md_allocator_i* heap = md_get_heap_allocator();
            p->in.build();
            md_array_resize(p->top, n, heap);
            md_array_resize(p->bot, n, heap);
            md_array_resize(p->thk, n, heap);
            p->desc = {};
            p->desc.beads        = p->in.beads;
            p->desc.grid         = &p->grid;
            p->desc.probe_radius = p->probe;
            p->desc.max_dist     = p->max_dist;
            p->phase = Phase_Work;
        });

        p->work_total = void_heightmap_num_patches(&p->grid);
        p->task_work = task_system::create_pool_task(STR_LIT("Void surface topography"), p->work_total,
            [p](uint32_t beg, uint32_t end, uint32_t) {
                for (uint32_t i = beg; i < end && !p->cancelled(); ++i) {
                    void_heightmap_eval_patches(p->top, p->bot, &p->desc, i, i + 1);
                    p->work_done.fetch_add(1, std::memory_order_relaxed);
                }
            }, 1);

        const task_system::ID finish = task_system::create_pool_task(STR_LIT("Void topography reduce"), [p, n]() {
            p->phase = Phase_Finish;
            if (!p->cancelled()) {
                for (size_t i = 0; i < n; ++i) {
                    const float t = p->top[i];
                    const float b = p->bot[i];
                    p->thk[i] = (isfinite(t) && isfinite(b)) ? t - b : NAN;
                }
                void_heightmap_stats(&p->st[HeightView_Top],       p->top, n);
                void_heightmap_stats(&p->st[HeightView_Bottom],    p->bot, n);
                void_heightmap_stats(&p->st[HeightView_Thickness], p->thk, n);
            }
            p->in.release();
            p->seconds = md_tick_to_seconds(md_tick_now() - p->t0);
            if (p->state == Pass_Running) p->state = p->cancelled() ? Pass_Failed : Pass_Done;
        });

        const task_system::ID install = task_system::create_main_task(STR_LIT("Void topography install"), [this, p]() {
            install_height_pass(p);
        });

        task_system::set_task_dependency(p->task_work, setup);
        task_system::set_task_dependency(finish, p->task_work);
        task_system::set_task_dependency(install, finish);

        height_pass = p;
        task_system::enqueue_task(setup);
    }

    void install_height_pass(HeightPass* p) {
        // Superseded, cancelled, or made for a field that has since been replaced
        if (p != height_pass || p->state != Pass_Done || p->result_gen != result_gen) {
            if (p == height_pass) height_pass = nullptr;
            delete p;
            return;
        }
        height_pass = nullptr;

        const bool pending = height_pending;    // The probe moved while this ran; it reruns on release
        clear_height();
        height_top = p->top; p->top = 0;
        height_bot = p->bot; p->bot = 0;
        height_thk = p->thk; p->thk = 0;
        MEMCPY(height_stats, p->st, sizeof(height_stats));
        height_nx         = 0;
        height_ny         = 0;
        height_probe      = p->probe;
        height_seconds    = p->seconds;
        has_height        = true;
        height_disp_dirty = true;
        height_pending    = pending;
        delete p;
    }

    // A progress bar for a running pass, with a way to stop it. Returns true when it was stopped.
    bool draw_pass_progress(const Pass* p, const char* what) {
        static const char* phase_lbl[] = { "building the neighbour structure", "", "reducing" };
        const int ph = p->phase.load();
        const float frac = p->progress();
        char overlay[128];
        if (ph == Phase_Work) snprintf(overlay, sizeof(overlay), "%s: %.0f%%", what, 100.0f * frac);
        else                  snprintf(overlay, sizeof(overlay), "%s: %s", what, phase_lbl[ph]);
        ImGui::PushID(what);
        ImGui::ProgressBar(frac, ImVec2(-80.0f, 0.0f), overlay);
        ImGui::SameLine();
        const bool stop = ImGui::Button("Cancel");
        ImGui::PopID();
        return stop;
    }

    const float* height_source(HeightView v) const {
        switch (v) {
        case HeightView_Top:       return height_top;
        case HeightView_Bottom:    return height_bot;
        case HeightView_Thickness: return height_thk;
        default:                   return nullptr;
        }
    }

    // Average the selected map down to at most HEIGHT_MAX_DIM per axis for the heat map. An open
    // column has no height, and ImPlot has no way to leave a cell empty, so it is drawn at the bottom
    // of the colour scale: the deepest thing on the map, which is what a hole is. For the thickness
    // that is zero, which is exactly right.
    void update_height_display() {
        height_disp_dirty = false;
        if (!has_height) return;
        const float* src = height_source(height_view);
        if (!src) return;

        const double nm = (double)ANGSTROM_PER_NM;
        const void_heightmap_stats_t& st = height_stats[height_view];
        if (st.num_valid > 0) {
            height_lo = st.min / nm;
            height_hi = st.max / nm;
        } else {
            height_lo = 0.0;
            height_hi = 1.0;
        }
        if (height_view == HeightView_Thickness) height_lo = 0.0;
        if (!(height_hi > height_lo)) height_hi = height_lo + 1.0e-3;

        const uint32_t dx = (uint32_t)grid.dim[0];
        const uint32_t dy = (uint32_t)grid.dim[1];
        const uint32_t sx = (dx + HEIGHT_MAX_DIM - 1) / HEIGHT_MAX_DIM;
        const uint32_t sy = (dy + HEIGHT_MAX_DIM - 1) / HEIGHT_MAX_DIM;
        height_nx = (dx + sx - 1) / sx;
        height_ny = (dy + sy - 1) / sy;
        md_array_resize(height_disp, (size_t)height_nx * height_ny, arena);

        for (uint32_t j = 0; j < height_ny; ++j) {
            for (uint32_t i = 0; i < height_nx; ++i) {
                double sum = 0.0;
                uint32_t cnt = 0;
                for (uint32_t y = j * sy; y < MIN(dy, (j + 1) * sy); ++y) {
                    for (uint32_t x = i * sx; x < MIN(dx, (i + 1) * sx); ++x) {
                        const float v = src[(size_t)y * dx + x];
                        if (isfinite(v)) { sum += (double)v; cnt += 1; }
                    }
                }
                const double v = cnt ? (sum / (double)cnt) / nm : height_lo;
                height_disp[(size_t)(height_ny - 1 - j) * height_nx + i] = (float)v;
            }
        }
    }

    // --- Volume rendering of the field --------------------------------------------------------

    // Upload the materialized field as a 3D texture, averaged down to VOL_MAX_DIM per axis. The
    // texture covers a whole number of blocks, so when the stride does not divide the grid it runs
    // a fraction of a block past the far faces; those texels average only the voxels that exist.
    void upload_volume() {
        vol_dirty = false;
        if (!field || !has_result) return;

        const int* d = grid.dim;
        const int maxd = MAX(d[0], MAX(d[1], d[2]));
        const int s = MAX(1, (maxd + VOL_MAX_DIM - 1) / VOL_MAX_DIM);
        const int dd[3] = { (d[0] + s - 1) / s, (d[1] + s - 1) / s, (d[2] + s - 1) / s };

        if (!gl::init_texture_3D(&vol_tex, dd[0], dd[1], dd[2], GL_R16F)) {
            MD_LOG_ERROR("Void analysis: could not create a %i x %i x %i volume texture", dd[0], dd[1], dd[2]);
            return;
        }

        if (s == 1) {
            gl::set_texture_3D_data(vol_tex, 0, field, GL_R32F);
        } else {
            const size_t count = (size_t)dd[0] * dd[1] * dd[2];
            const size_t bytes = count * sizeof(float);
            float* buf = (float*)md_alloc(md_get_heap_allocator(), bytes);
            defer { md_free(md_get_heap_allocator(), buf, bytes); };

            const float* src = field;
            task_system::ID task = task_system::create_pool_task(STR_LIT("Void field downsample"), (uint32_t)dd[2],
                [buf, src, d, dd, s](uint32_t range_beg, uint32_t range_end, uint32_t) {
                    for (uint32_t k = range_beg; k < range_end; ++k) {
                        for (int j = 0; j < dd[1]; ++j) {
                            for (int i = 0; i < dd[0]; ++i) {
                                double sum = 0.0;
                                int cnt = 0;
                                for (int z = (int)k * s; z < MIN(d[2], ((int)k + 1) * s); ++z) {
                                    for (int y = j * s; y < MIN(d[1], (j + 1) * s); ++y) {
                                        const float* row = src + ((size_t)z * d[1] + y) * d[0];
                                        for (int x = i * s; x < MIN(d[0], (i + 1) * s); ++x) {
                                            sum += (double)row[x];
                                            cnt += 1;
                                        }
                                    }
                                }
                                buf[((size_t)k * dd[1] + j) * dd[0] + i] = cnt ? (float)(sum / (double)cnt) : 0.0f;
                            }
                        }
                    }
                }, 1);
            task_system::enqueue_task(task);
            task_system::task_wait_for(task);

            gl::set_texture_3D_data(vol_tex, 0, buf, GL_R32F);
        }

        vol_stride = s;
        MEMCPY(vol_dim, dd, sizeof(dd));
        vol_spacing = vec3_t{ grid.spacing.x * (float)s, grid.spacing.y * (float)s, grid.spacing.z * (float)s };
        const vec3_t min_aabb = grid.origin;
        const vec3_t max_aabb = {
            grid.origin.x + vol_spacing.x * (float)dd[0],
            grid.origin.y + vol_spacing.y * (float)dd[1],
            grid.origin.z + vol_spacing.z * (float)dd[2],
        };
        vol_model = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
        vol_z0 = min_aabb.z;
        vol_z1 = max_aabb.z;
    }

    // Clearance range the colours span. The lower end is always zero, so a colour means the same
    // clearance whatever the probe is; the probe only decides where the transparency ends.
    double vol_tf_max() const {
        return MAX((double)vol_upper, 1.0e-3);
    }

    // Transparent below the probe radius - those are places a probe centre of that size cannot be -
    // and above the upper limit, which by default is the query range, so vacuum no bead is within
    // reach of is not drawn as the largest void in the box. In between the opacity rises with the
    // clearance, so wide cavities read as the body of the pore space and the margins as a haze.
    void update_tf() {
        vol_tf_dirty = false;
        uint32_t px[VOL_TF_RES];
        const double tf_max = vol_tf_max();
        const double R = (double)probe_radius;
        for (int i = 0; i < VOL_TF_RES; ++i) {
            const float  t = (float)i / (float)(VOL_TF_RES - 1);
            const double d = (double)t * tf_max;
            ImVec4 c = ImPlot::SampleColormap(t, vol_colormap);
            float a = 0.0f;
            if (d >= R && i < VOL_TF_RES - 1) {
                const double u = (tf_max > R) ? (d - R) / (tf_max - R) : 1.0;
                a = vol_opacity * (float)(0.25 + 0.75 * u);
            }
            c.w = CLAMP(a, 0.0f, 1.0f);
            px[i] = ImGui::ColorConvertFloat4ToU32(c);
        }
        gl::init_texture_2D(&vol_tf_tex, VOL_TF_RES, 1, GL_RGBA8);
        gl::set_texture_2D_data(vol_tf_tex, 0, px, GL_RGBA8);
    }

    void draw_volume(const ApplicationState& state) {
        // The texture belongs to the last materialized field; without one it would be a stale picture
        if (!field) return;
        if (vol_dirty)    upload_volume();
        if (vol_tf_dirty) update_tf();
        if (!vol_tex || !vol_tf_tex) return;

        // The reporting region, as a clip in the texture's z. The texture may run past the grid by
        // part of a block, which is why this is taken against its own extent.
        vec3_t clip_min = {0, 0, 0};
        vec3_t clip_max = {1, 1, 1};
        if (vol_clip_region && region_end > region_beg && (region_beg > 0 || region_end < stats.num_slabs) && vol_z1 > vol_z0) {
            const void_profile_t p = make_profile();
            clip_min.z = CLAMP((float)((void_profile_z_lo(&p, region_beg)     - vol_z0) / (vol_z1 - vol_z0)), 0.0f, 1.0f);
            clip_max.z = CLAMP((float)((void_profile_z_hi(&p, region_end - 1) - vol_z0) / (vol_z1 - vol_z0)), 0.0f, 1.0f);
        }

        volume::RenderDesc desc = {};
        desc.render_target.depth  = state.gbuffer.tex.depth;
        desc.render_target.color  = state.gbuffer.tex.transparency;
        desc.render_target.width  = state.gbuffer.width;
        desc.render_target.height = state.gbuffer.height;

        desc.texture.density_volume    = vol_tex;
        desc.texture.transfer_function = vol_tf_tex;

        desc.matrix.model    = vol_model;
        desc.matrix.view     = state.view.param.matrix.curr.view;
        desc.matrix.proj     = state.view.param.matrix.curr.proj;
        desc.matrix.inv_proj = state.view.param.matrix.inv.proj;

        desc.clip_volume.min = clip_min;
        desc.clip_volume.max = clip_max;

        desc.temporal.enabled = state.visuals.temporal_aa.enabled;

        desc.dvr.enabled      = true;
        desc.dvr.min_tf_value = 0.0f;
        desc.dvr.max_tf_value = (float)vol_tf_max();

        desc.shading.env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f;
        desc.shading.roughness    = 0.3f;
        desc.shading.dir_radiance = {10, 10, 10};
        desc.shading.ior          = 1.5f;
        desc.shading.exposure     = state.visuals.tonemapping.exposure;
        desc.shading.gamma        = state.visuals.tonemapping.gamma;

        desc.voxel_spacing = vol_spacing;

        volume::render_volume(desc);
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
        md_file_printf(file, "# probe_radius_nm,%.6f\n",    (double)probe_radius / nm);
        md_file_printf(file, "# region_slabs,%u,%u\n", region_beg, region_end);
        md_file_printf(file, "# density_unit,%s\n", density_unit());
        md_file_printf(file, "z_lo_nm,z_hi_nm,n_voxels,n_solid,porosity,accessible_fraction,void_volume_nm3,accessible_volume_nm3,density\n");
        for (uint32_t sl = 0; sl < stats.num_slabs; ++sl) {
            const double v_slab = (double)void_profile_num_total(&p, sl, sl + 1) * stats.voxel_volume;
            md_file_printf(file, "%.6f,%.6f,%llu,%llu,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                void_profile_z_lo(&p, sl) / nm,
                void_profile_z_hi(&p, sl) / nm,
                (unsigned long long)void_profile_num_total(&p, sl, sl + 1),
                (unsigned long long)void_profile_num_solid(&p, sl, sl + 1),
                void_profile_porosity(&p, sl, sl + 1),
                void_profile_accessible_fraction(&p, sl, sl + 1, (double)probe_radius),
                v_slab * void_profile_porosity(&p, sl, sl + 1) / nm3,
                void_profile_accessible_volume(&p, sl, sl + 1, (double)probe_radius) / nm3,
                density_over(sl, sl + 1));
        }
    }

    // V(R,z) in long form: one row per (radius, slab) pair, which is what a plotting tool wants and
    // what a matrix written as a grid of numbers is not. The radii are the knots of the curve, so a
    // linear interpolation between rows reproduces it exactly.
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
        double max_rad = fill_radii(radii, sys, state.num_atoms, (double)probe_radius);

        md_coord_stream_t coords = md_coord_stream_from_aos((const float*)state.xyz, sizeof(vec3_t), NULL, state.num_atoms);

        md_spatial_acc_t acc = {};
        acc.alloc = md_temp_allocator(temp_scope);
        md_spatial_acc_desc_t desc = {};
        desc.coords   = &coords;
        desc.cutoff   = ACC_CUTOFF;
        desc.unitcell = &state.unitcell;
        md_spatial_acc_init(&acc, &desc);
        defer { md_spatial_acc_free(&acc); };
        const void_beads_t beads = void_beads(&acc, radii, state.num_atoms);

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
        const double query_max_dist = 2.0 * max_rad + 1.0;

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
                        px[k] = state.xyz[i].x + a * dir[k * 3 + 0];
                        py[k] = state.xyz[i].y + a * dir[k * 3 + 1];
                        pz[k] = state.xyz[i].z + a * dir[k * 3 + 2];
                    }

                    void_beads_nearest(&beads, px, py, pz, (size_t)NP, query_max_dist, idx, NULL);

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
        sasa.seconds = md_tick_to_seconds(md_tick_now() - t0);

        has_sasa = true;
    }

    // --- Pore network -------------------------------------------------------------------------

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

    void cancel_network_pass() {
        if (!net_pass) return;
        net_pass->cancel = true;
        net_pass = nullptr;
    }

    // Also abandons a build in progress: whatever it would produce belongs to what is being cleared
    void clear_network() {
        cancel_network_pass();
        pore_network_free(&net);
        md_array_free(net_class, md_get_heap_allocator());
        md_array_free(net_route, md_get_heap_allocator());
        md_array_free(net_selected, md_get_heap_allocator());
        md_array_free(net_region, md_get_heap_allocator());
        net_class = 0;
        net_route = 0;
        net_selected = 0;
        net_num_selected = 0;
        net_region = 0;
        net_num_region = 0;
        net_drawn_throats = 0;
        net_class_r = -1.0f;
        MEMSET(net_class_count, 0, sizeof(net_class_count));
        net_route_bottleneck = 0.0;
        net_route_length = 0.0;
        net_hovered = PORE_INVALID;
        net_picking = {};
        has_net = false;
        net_seconds = 0.0;
    }

    // The build is one ordered sweep and does not divide over threads, so it runs as a single pool
    // task. It holds a reference to the field, which keeps it alive if a new one is installed while
    // the build runs; the network that would have come of it is then discarded.
    void compute_network() {
        error[0] = '\0';

        channel_field_t f;
        if (!build_channel_field(&f)) {
            snprintf(error, sizeof(error), "Compute the field first, with 'Materialize field' enabled");
            return;
        }

        cancel_network_pass();

        NetworkPass* p = new NetworkPass();
        p->t0         = md_tick_now();
        p->keep       = field_owner;
        p->f          = f;
        p->r_min      = (double)net_r_min;
        p->merge      = (double)net_merge * (double)MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        p->result_gen = result_gen;
        p->work_total = 1000;

        p->task_work = task_system::create_pool_task(STR_LIT("Void pore network"), [p]() {
            md_allocator_i* heap = md_get_heap_allocator();
            p->phase = Phase_Work;
            if (!p->cancelled() && !pore_network_build(&p->net, &p->f, p->r_min, p->merge, heap, network_progress, p)) {
                if (!p->cancelled()) {
                    snprintf(p->error, sizeof(p->error), "Nothing at or above %.2f nm to build a network from", p->r_min / ANGSTROM_PER_NM);
                }
                p->cancel = true;
            }
            p->keep.reset();

            p->phase = Phase_Finish;
            if (!p->cancelled() && pore_network_widest_route(&p->route, &p->bottleneck, &p->net, heap)) {
                // Pore to throat to pore along the route, plus the straight legs in from each face.
                // Over the depth of the box that is the tortuosity of the route at the resolution of
                // the graph.
                const pore_network_t& net = p->net;
                const size_t n = md_array_size(p->route);
                double len = 0.0;
                for (size_t i = 0; i + 1 < n; ++i) {
                    const uint32_t e = pore_network_find_edge(&net, p->route[i], p->route[i + 1]);
                    if (e == PORE_INVALID) continue;
                    vec3_t a, s, b;
                    pore_network_edge_points(&a, &s, &b, &net, e);
                    len += vec3_length(vec3_sub(s, a)) + vec3_length(vec3_sub(b, s));
                }
                const float z_top = net.box_min[2] + net.box_ext[2];
                len += fabs((double)z_top - (double)net.vertices[p->route[0]].pos[2]);
                len += fabs((double)net.vertices[p->route[n - 1]].pos[2] - (double)net.box_min[2]);
                p->route_length = len;
            }

            p->seconds = md_tick_to_seconds(md_tick_now() - p->t0);
            p->state   = p->cancelled() ? Pass_Failed : Pass_Done;
        });

        const task_system::ID install = task_system::create_main_task(STR_LIT("Void pore network install"), [this, p]() {
            install_network_pass(p);
        });
        task_system::set_task_dependency(install, p->task_work);

        net_pass = p;
        task_system::enqueue_task(p->task_work);
    }

    // Main thread. The previous network stays up until this replaces it.
    void install_network_pass(NetworkPass* p) {
        if (p != net_pass || p->state != Pass_Done || p->result_gen != result_gen) {
            if (p == net_pass) {
                net_pass = nullptr;
                if (p->error[0]) snprintf(error, sizeof(error), "%s", p->error);
            }
            delete p;
            return;
        }
        net_pass = nullptr;

        clear_network();
        net = p->net;
        p->net = {};
        net_route = p->route;
        p->route = 0;
        net_route_bottleneck = p->bottleneck;
        net_route_length     = p->route_length;
        net_seconds          = p->seconds;
        delete p;

        has_net = true;
        md_array_resize(net_selected, md_array_size(net.vertices), md_get_heap_allocator());
        md_array_resize(net_region,   md_array_size(net.vertices), md_get_heap_allocator());
        clear_network_selection();
        clear_network_region();
        update_network_classes();
    }

    // What every pore is to a probe of the current radius. A sweep of the throats, so it follows the
    // probe slider rather than waiting for a button.
    void update_network_classes() {
        net_class_r = probe_radius;
        MEMSET(net_class_count, 0, sizeof(net_class_count));
        const size_t V = md_array_size(net.vertices);
        md_array_resize(net_class, V, md_get_heap_allocator());
        if (V == 0) return;
        pore_network_classify(net_class, &net, (double)probe_radius, md_get_heap_allocator());
        for (size_t i = 0; i < V; ++i) net_class_count[net_class[i]] += 1;
    }

    bool net_pore_visible(uint32_t i) const {
        return net_show_small || net_class[i] != PORE_CLASS_SMALL;
    }

    // Contrast against the background, for what is in focus: near black on a light background and
    // near white on a dark one. No hue survives both, and every hue is already spent on the classes.
    static uint32_t focus_color(const ApplicationState& state) {
        const vec3_t c = state.visuals.background.color;
        const float  lum = (0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z) * state.visuals.background.intensity;
        return (lum > 0.5f) ? PORE_FOCUS_ON_LIGHT : PORE_FOCUS_ON_DARK;
    }

    bool net_pore_selected(uint32_t i) const {
        return i < md_array_size(net_selected) && net_selected[i];
    }

    void clear_network_selection() {
        for (size_t i = 0; i < md_array_size(net_selected); ++i) net_selected[i] = 0;
        net_num_selected = 0;
    }

    void clear_network_region() {
        if (net_num_region == 0) return;
        for (size_t i = 0; i < md_array_size(net_region); ++i) net_region[i] = 0;
        net_num_region = 0;
    }

    // The pores whose centres project inside the box, out of those drawn - or, when removing, out of
    // those selected - exactly as the atoms are chosen.
    void update_network_region(const InteractionSurfaceEvent& ev) {
        clear_network_region();
        const size_t V = md_array_size(net.vertices);
        if (md_array_size(net_region) != V || md_array_size(net_class) != V) return;
        const bool removing = ev.selection_mode == InteractionSelectionMode::Remove;
        net_region_removing = removing;
        for (uint32_t i = 0; i < (uint32_t)V; ++i) {
            if (removing ? !net_selected[i] : !net_pore_visible(i)) continue;
            const pore_vertex_t& v = net.vertices[i];
            const vec4_t c = mat4_mul_vec4(ev.world_to_clip, vec4_set(v.pos[0], v.pos[1], v.pos[2], 1.0f));
            if (!(c.w > 0.0f)) continue;    // Behind the camera, where the divide would mirror it into view
            const float sx = ( c.x / c.w * 0.5f + 0.5f) * ev.surface_size.x;
            const float sy = (-c.y / c.w * 0.5f + 0.5f) * ev.surface_size.y;
            if (ev.region_min.x <= sx && sx <= ev.region_max.x && ev.region_min.y <= sy && sy <= ev.region_max.y) {
                net_region[i] = 1;
                net_num_region += 1;
            }
        }
    }

    void set_network_selected(uint32_t i, bool on) {
        if (i >= md_array_size(net_selected) || (net_selected[i] != 0) == on) return;
        net_selected[i] = on ? 1 : 0;
        if (on) net_num_selected += 1; else net_num_selected -= 1;
    }

    // Selection follows the atoms: a click picks one pore, shift adds and shift with the right
    // button removes, a shift drag adds or removes everything in the box, and a plain click on
    // nothing clears. Hover only feeds the tooltip and the focus drawing.
    void network_interaction(const InteractionSurfaceEvent& ev) {
        if (ev.surface_id != interaction_surface_main) return;
        const size_t V = md_array_size(net.vertices);
        const bool live = has_net && show_net;
        const bool ours = live && ev.hit.domain == PickingDomain_PoreNetwork && ev.hit.local_idx < V;

        if (ev.kind != InteractionSurfaceEventKind::RegionSelect) {
            clear_network_region();
        }

        if (ev.kind == InteractionSurfaceEventKind::RegionSelect) {
            net_hovered = PORE_INVALID;
            if (!live) return;
            update_network_region(ev);
            if (ev.region_phase == InteractionSurfaceEventPhase::Commit) {
                const bool on = ev.selection_mode != InteractionSelectionMode::Remove;
                for (uint32_t i = 0; i < (uint32_t)V; ++i) {
                    if (net_region[i]) set_network_selected(i, on);
                }
                clear_network_region();
            }
        } else if (ev.kind == InteractionSurfaceEventKind::Hover) {
            net_hovered = ours ? ev.hit.local_idx : PORE_INVALID;
        } else if (ev.kind == InteractionSurfaceEventKind::Click && live) {
            if (ours) {
                switch (ev.selection_mode) {
                case InteractionSelectionMode::None:
                    clear_network_selection();
                    set_network_selected(ev.hit.local_idx, true);
                    break;
                case InteractionSelectionMode::Append: set_network_selected(ev.hit.local_idx, true);  break;
                case InteractionSelectionMode::Remove: set_network_selected(ev.hit.local_idx, false); break;
                }
            } else if (ev.hit.domain == 0 && ev.selection_mode != InteractionSelectionMode::Append) {
                clear_network_selection();
            }
        }
    }

    // Drawn into the world queue, which is rendered with the depth test into the G-buffer: lit by the
    // same pass as the structure, hidden by what is in front of it, and picked by depth rather than
    // by draw order.
    void draw_network_3d(const ApplicationState& state) {
        if (net_class_r != probe_radius || md_array_size(net_class) != md_array_size(net.vertices)) {
            update_network_classes();
        }

        const size_t V = md_array_size(net.vertices);
        const size_t E = md_array_size(net.edges);
        if (V == 0) return;

        md_allocator_i* frame = state.allocator.frame;

        // Points and lines have no surface to light, so their normal is made to face the camera:
        // the deferred pass then shades them at their full colour from any direction.
        const vec3_t facing = vec3_normalize(vec3_from_vec4(state.view.param.matrix.inv.view.col[2]));
        const uint32_t focus = focus_color(state);

        auto pore_pos = [&](uint32_t i) {
            return vec3_t{ net.vertices[i].pos[0], net.vertices[i].pos[1], net.vertices[i].pos[2] };
        };

        md_array(immediate::Vertex) lines = 0;
        md_array_ensure(lines, 4 * MIN(E, PORE_MAX_DRAWN_THROATS) + 1024, frame);
        auto seg = [&](vec3_t a, vec3_t b, uint32_t col) {
            const immediate::Vertex va = { a, col, facing, 0xFFFFFFFFu };
            const immediate::Vertex vb = { b, col, facing, 0xFFFFFFFFu };
            md_array_push(lines, va, frame);
            md_array_push(lines, vb, frame);
        };
        auto throat = [&](uint32_t e, uint32_t col) {
            vec3_t a, s, b;
            pore_network_edge_points(&a, &s, &b, &net, e);
            seg(a, s, col);
            seg(s, b, col);
        };
        auto wire_sphere = [&](vec3_t c, float r, uint32_t col, int stacks, int slices) {
            const float pi = 3.14159265358979f;
            for (int i = 1; i < stacks; ++i) {          // Parallels
                const float t = pi * (float)i / (float)stacks;
                const float z = r * cosf(t), rr = r * sinf(t);
                for (int j = 0; j < slices; ++j) {
                    const float p0 = 2.0f * pi * (float)j / (float)slices, p1 = 2.0f * pi * (float)(j + 1) / (float)slices;
                    seg(vec3_t{ c.x + rr * cosf(p0), c.y + rr * sinf(p0), c.z + z }, vec3_t{ c.x + rr * cosf(p1), c.y + rr * sinf(p1), c.z + z }, col);
                }
            }
            for (int j = 0; j < slices; j += 2) {       // Meridians, every other slice
                const float p = 2.0f * pi * (float)j / (float)slices;
                for (int i = 0; i < stacks; ++i) {
                    const float t0 = pi * (float)i / (float)stacks, t1 = pi * (float)(i + 1) / (float)stacks;
                    seg(vec3_t{ c.x + r * sinf(t0) * cosf(p), c.y + r * sinf(t0) * sinf(p), c.z + r * cosf(t0) },
                        vec3_t{ c.x + r * sinf(t1) * cosf(p), c.y + r * sinf(t1) * sinf(p), c.z + r * cosf(t1) }, col);
                }
            }
        };

        // A pore in focus - hovered or selected - shows its largest sphere and each of its throats at
        // its own width. First, so that where it coincides with the ordinary drawing of the same
        // throat, the depth test keeps the focus colour.
        auto draw_focus = [&](uint32_t i) {
            wire_sphere(pore_pos(i), net.vertices[i].radius, focus, 10, 20);
            for (uint32_t j = net.adj_offset[i]; j < net.adj_offset[i + 1]; ++j) {
                const uint32_t e = net.adj[j];
                const pore_edge_t& edge = net.edges[e];
                throat(e, focus);
                wire_sphere(vec3_t{ edge.pos[0], edge.pos[1], edge.pos[2] }, edge.radius, focus, 6, 12);
            }
        };
        const bool hovered = net_hovered < V;
        if (hovered) draw_focus(net_hovered);
        if (net_num_selected > 0 && net_num_selected <= PORE_MAX_FOCUS_DRAWN) {
            for (uint32_t i = 0; i < (uint32_t)V; ++i) {
                if (net_selected[i] && i != net_hovered) draw_focus(i);
            }
        }

        // The voxel that limits the widest route, at r_c
        const size_t nr = md_array_size(net_route);
        const bool route = show_net_route && nr > 0;
        if (route && net.has_r_c) {
            wire_sphere(vec3_t{ net.throat[0], net.throat[1], net.throat[2] }, (float)net.r_c, PORE_ROUTE_COLOR, 8, 16);
        }

        // Throats. Widest first, so the ones a probe of this radius passes are a prefix. An open
        // throat joins two pores of the same component, so either end gives its class.
        net_drawn_throats = 0;
        for (size_t e = 0; e < E && net_drawn_throats < PORE_MAX_DRAWN_THROATS; ++e, ++net_drawn_throats) {
            const pore_edge_t& edge = net.edges[e];
            const bool open = edge.radius >= probe_radius;
            if (!open && !net_show_small) break;
            throat((uint32_t)e, open ? PORE_CLASS_COLOR[net_class[edge.a]] : PORE_CLASS_COLOR[PORE_CLASS_SMALL]);
        }

        immediate::Scope scope(state.gfx.world, "void_pore_network");
        immediate::lines(scope, lines, md_array_size(lines));

        // The widest route through the film as a tube, so it reads as a path rather than as one more
        // line among thousands. Sized to the grid and to r_c, not to the screen.
        if (route) {
            const float spacing = MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
            const float tube = MAX(0.5f * spacing, net.has_r_c ? 0.15f * (float)net.r_c : 0.0f);
            for (size_t i = 0; i + 1 < nr; ++i) {
                const uint32_t e = pore_network_find_edge(&net, net_route[i], net_route[i + 1]);
                if (e == PORE_INVALID) continue;
                vec3_t a, s, b;
                pore_network_edge_points(&a, &s, &b, &net, e);
                immediate::cylinder(scope, a, s, tube, PORE_ROUTE_COLOR, 0xFFFFFFFFu, 8);
                immediate::cylinder(scope, s, b, tube, PORE_ROUTE_COLOR, 0xFFFFFFFFu, 8);
            }
            const vec3_t first = pore_pos(net_route[0]);
            const vec3_t last  = pore_pos(net_route[nr - 1]);
            immediate::cylinder(scope, first, vec3_t{ first.x, first.y, net.box_min[2] + net.box_ext[2] }, tube, PORE_ROUTE_COLOR, 0xFFFFFFFFu, 8);
            immediate::cylinder(scope, last,  vec3_t{ last.x,  last.y,  net.box_min[2] }, tube, PORE_ROUTE_COLOR, 0xFFFFFFFFu, 8);
        }

        // Pores, as points that keep their size on screen. Those in focus go first and larger, so the
        // depth test keeps their colour where the ordinary point of the same pore lands on it.
        const bool pickable = net_picking.domain == PickingDomain_PoreNetwork && net_picking.end - net_picking.beg == (uint32_t)V;
        if (pickable) immediate::set_picking_base_idx(scope, net_picking.beg);
        const uint32_t no_pick = 0xFFFFFFFFu;

        md_array(immediate::Vertex) pts = 0;
        md_array_ensure(pts, V, frame);
        auto pore_vertex = [&](uint32_t i, uint32_t col) {
            const immediate::Vertex v = { pore_pos(i), col, facing, pickable ? i : no_pick };
            md_array_push(pts, v, frame);
        };
        // A box being dragged previews what it will take: what it adds, or, when removing, what it
        // leaves selected.
        const bool removing_region = net_num_region > 0 && net_region_removing;
        if (hovered) pore_vertex(net_hovered, focus);
        for (uint32_t i = 0; i < (uint32_t)V && (net_num_selected > 0 || net_num_region > 0); ++i) {
            if (i == net_hovered) continue;
            const bool sel = net_selected[i] != 0;
            const bool reg = net_region[i] != 0;
            if (removing_region ? (sel && !reg) : (sel || reg)) pore_vertex(i, focus);
        }
        if (md_array_size(pts) > 0) {
            immediate::set_point_size(scope, PORE_POINT_SIZE_FOCUS);
            immediate::points(scope, pts, md_array_size(pts));
            md_array_shrink(pts, 0);
        }
        for (uint32_t i = 0; i < (uint32_t)V; ++i) {
            if (net_pore_visible(i)) pore_vertex(i, PORE_CLASS_COLOR[net_class[i]]);
        }
        immediate::set_point_size(scope, PORE_POINT_SIZE);
        immediate::points(scope, pts, md_array_size(pts));
    }

    void network_tooltip(md_strb_t* sb, uint32_t i) const {
        if (!has_net || i >= md_array_size(net.vertices)) return;
        const pore_vertex_t& v = net.vertices[i];
        const double nm  = (double)ANGSTROM_PER_NM;
        const double nm3 = nm * nm * nm;

        md_strb_fmt(sb, "Pore %u\n", i);
        md_strb_fmt(sb, "Radius: %.2f nm (largest sphere that fits)\n", v.radius / nm);
        md_strb_fmt(sb, "Volume: %.4g nm^3\n", (double)v.num_voxels * net.voxel_volume / nm3);
        md_strb_fmt(sb, "Centre: %.1f, %.1f, %.1f nm\n", v.pos[0] / nm, v.pos[1] / nm, v.pos[2] / nm);

        const uint32_t deg = net.adj_offset[i + 1] - net.adj_offset[i];
        if (deg > 0) {
            float t_min = FLT_MAX, t_max = 0.0f;
            uint32_t passable = 0;
            for (uint32_t j = net.adj_offset[i]; j < net.adj_offset[i + 1]; ++j) {
                const float r = net.edges[net.adj[j]].radius;
                t_min = MIN(t_min, r);
                t_max = MAX(t_max, r);
                passable += (r >= probe_radius);
            }
            md_strb_fmt(sb, "Throats: %u, %.2f to %.2f nm, %u wide enough for the probe\n", deg, t_min / nm, t_max / nm, passable);
        } else {
            md_strb_fmt(sb, "Throats: none\n");
        }
        if (v.face_top    >= 0.0f) md_strb_fmt(sb, "Opens on the top face, %.2f nm wide\n", v.face_top / nm);
        if (v.face_bottom >= 0.0f) md_strb_fmt(sb, "Opens on the bottom face, %.2f nm wide\n", v.face_bottom / nm);

        if (i < md_array_size(net_class)) {
            md_strb_fmt(sb, "Probe %.2f nm: %s", probe_radius / nm, pore_class_lbl[net_class[i]]);
        }
        for (size_t k = 0; k < md_array_size(net_route); ++k) {
            if (net_route[k] == i) {
                md_strb_fmt(sb, "\nOn the widest route, pore %zu of %zu from the top", k + 1, md_array_size(net_route));
                break;
            }
        }
    }

    void draw_network_section() {
        const double nm = (double)ANGSTROM_PER_NM;

        if (!ImGui::CollapsingHeader("Pore network", ImGuiTreeNodeFlags_DefaultOpen)) return;
        ImGui::TextDisabled("The void as pores and the throats between them, in the 3D view. Through z, open at both ends.");

        float r_min_nm = net_r_min / ANGSTROM_PER_NM;
        if (ImGui::SliderFloat("Smallest clearance (nm)", &r_min_nm, 0.05f, 5.0f, "%.2f")) {
            net_r_min = r_min_nm * ANGSTROM_PER_NM;
        }
        ImGui::SetItemTooltip("The narrowest throat and the smallest pore the network can have, and what bounds\n"
                              "the memory: a voxel with less clearance never enters the pass.");

        const float spacing = MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        ImGui::SliderFloat("Merge depth (voxels)", &net_merge, 0.0f, 4.0f, "%.2f");
        ImGui::SetItemTooltip("Two pores are one when the throat between them is less than this far below the\n"
                              "smaller one's radius - a shallow dip rather than a constriction. Every bump the\n"
                              "voxelization leaves on a wall is a dip of a fraction of a voxel, so below about\n"
                              "half a voxel the network fills with pores that are not there. Currently %.2f nm.",
                              (spacing > 0.0f) ? net_merge * spacing / ANGSTROM_PER_NM : 0.0f);

        ImGui::BeginDisabled(!has_result || !field);
        if (ImGui::Button(net_pass ? "Restart##net" : "Build network")) {
            compute_network();
        }
        ImGui::EndDisabled();
        if (!field) {
            ImGui::SameLine();
            ImGui::TextDisabled("(needs a materialized field)");
        } else if (net_pass) {
            ImGui::SameLine();
            if (draw_pass_progress(net_pass, "Pore network")) cancel_network_pass();
        }

        if (!has_net) return;

        const size_t V = md_array_size(net.vertices);
        const size_t E = md_array_size(net.edges);
        ImGui::Text("Built in %.2f s: %zu pores, %zu throats, mean coordination %.2f",
                    net_seconds, V, E, V > 0 ? 2.0 * (double)E / (double)V : 0.0);

        if (net.has_r_c) {
            ImGui::Text("Critical radius r_c: %.3f nm", net.r_c / nm);
            ImGui::SetItemTooltip("The largest probe which still gets from one z face to the other: the tightest\n"
                                  "throat on the widest route. Exact - read off the voxels in the same pass, not off\n"
                                  "the graph, so merging pores does not move it.");
            const size_t nr = md_array_size(net_route);
            if (nr > 0) {
                ImGui::Text("Widest route: %zu pores, tortuosity %.2f", nr, net_route_length / MAX(1.0e-6, (double)net.box_ext[2]));
                ImGui::SetItemTooltip("Through pore centres and throats, with straight legs in from each face, over the\n"
                                      "depth of the box. A route through the graph, so it is as coarse as the pores are.%s",
                                      (fabs(net_route_bottleneck - net.r_c) > 1.0e-3 * nm)
                                      ? "\nIts narrowest throat reads wider than r_c because it passes a merged pore." : "");
            }
        } else {
            ImGui::Text("Nothing gets through at %.2f nm or above", net_r_min / nm);
        }

        ImGui::Checkbox("Show in 3D", &show_net);
        ImGui::SameLine();
        ImGui::Checkbox("Widest route", &show_net_route);
        ImGui::SameLine();
        ImGui::Checkbox("Too narrow for the probe", &net_show_small);
        if (net_drawn_throats >= PORE_MAX_DRAWN_THROATS) {
            ImGui::TextColored({1.0f, 0.8f, 0.35f, 1.0f}, "Drawing the %zu widest throats only", PORE_MAX_DRAWN_THROATS);
        }

        if (net_class_r != probe_radius) update_network_classes();
        ImGui::Text("At probe %.2f nm:", probe_radius / nm);
        static const int order[] = { PORE_CLASS_SPANNING, PORE_CLASS_TOP, PORE_CLASS_BOTTOM, PORE_CLASS_CLOSED, PORE_CLASS_SMALL };
        for (int c : order) {
            ImGui::ColorButton(pore_class_lbl[c], ImGui::ColorConvertU32ToFloat4(PORE_CLASS_COLOR[c] | 0xFF000000u),
                               ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(10, 10));
            ImGui::SameLine();
            ImGui::Text("%s: %u", pore_class_lbl[c], net_class_count[c]);
        }
        ImGui::TextDisabled("Hover a pore in the 3D view for its details. Click to select, shift click or drag to add,\n"
                            "shift right click or drag to remove.");
        if (net_num_selected > 0) {
            ImGui::Text("%u selected%s", net_num_selected,
                        net_num_selected > PORE_MAX_FOCUS_DRAWN ? " (spheres drawn up to 2000)" : "");
            ImGui::SameLine();
            if (ImGui::SmallButton("Clear selection")) clear_network_selection();
        }
    }

    // Porosity and accessible volume. Both are V(R) = Vol[d > R] over a slab range: porosity is the
    // R = 0 end of the same curve, which is why they share a section and a plot rather than sitting
    // in two places able to disagree.
    void draw_porosity_section() {
        const void_profile_t p = make_profile();
        const double nm  = (double)ANGSTROM_PER_NM;
        const double nm3 = nm * nm * nm;

        if (!ImGui::CollapsingHeader("Porosity and accessible volume", ImGuiTreeNodeFlags_DefaultOpen)) return;

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
                              "The solid it is measured against is the union of the bead spheres, so this moves with\n"
                              "the radius source like everything else here.");
        ImGui::SameLine();
        ImGui::TextDisabled("(solid %.4f, void volume %.4g nm^3)", 1.0 - phi, phi * v_region / nm3);

        ImGui::Text("V(R) at R = %.2f nm: %.4g nm^3, V(R)/V = %.4f", probe_radius / nm, v_acc / nm3, phi_acc);
        ImGui::SetItemTooltip("Volume whose distance to the nearest bead surface exceeds R, i.e. where the CENTRE of a\n"
                              "probe of radius R fits. Not the volume the probe body fills, which is that set dilated\n"
                              "by R and is larger, and not a statement that the volume can be reached from outside -\n"
                              "a closed cavity counts here and is invisible to infiltration.");
        ImGui::SameLine();
        ImGui::TextDisabled("(%.1f%% of the pore space)", (phi > 0.0) ? 100.0 * phi_acc / phi : 0.0);

        const double rho = density_over(region_beg, region_end);
        ImGui::Text("Density: %.4g %s", rho, density_unit());
        if (has_mass()) {
            ImGui::SetItemTooltip("Mass of the atoms in the region over the volume of the region - the same volume the\n"
                                  "porosity is a fraction of. With a periodic z the atoms are wrapped into the box;\n"
                                  "with an open z an atom outside the grid is outside the region and not counted.");
            // Bulk density is the skeletal density diluted by the pore space, rho = (1 - phi) rho_s,
            // so the skeletal one is what the solid itself weighs per volume it occupies.
            ImGui::SameLine();
            ImGui::TextDisabled("(skeletal %.4g kg/m^3)", (phi < 1.0) ? rho / (1.0 - phi) : 0.0);
            ImGui::SetItemTooltip("Mass over the solid volume alone: rho / (1 - porosity). For a bead model this says\n"
                                  "how heavy the chosen radii make the solid, which is a check on the radii.");
        } else {
            ImGui::SetItemTooltip("No atom carries a mass, so this is the number density instead.");
        }

        // Three views side by side. V(R)/V against R, whose left endpoint is the porosity and which the
        // probe radius reads off; the same against z, where the R = 0 curve is the porosity profile
        // and the gap to the probe curve is the pore volume the probe is too large for; and the
        // density against z on the same slabs, which is what says where the film actually is. The
        // two z plots share their axis, so zooming one zooms the other.
        const bool have_rows = md_array_size(prof_z) > 1 && md_array_size(prof_acc) == md_array_size(prof_z) &&
                               md_array_size(prof_rho) == md_array_size(prof_z);
        if (ImGui::BeginTable("##profiles", 3, ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_Resizable | ImGuiTableFlags_NoSavedSettings)) {
            const float plot_h = 220.0f;
            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            if (md_array_size(acc_r) > 1 && ImPlot::BeginPlot("##acc_curve", ImVec2(-1, plot_h))) {
                ImPlot::SetupAxes("Probe radius R (nm)", "V(R) / V");
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, MAX(1.0e-4, phi * 1.05), ImPlotCond_Always);
                ImPlot::PlotLine("V(R)/V", acc_r, acc_frac, (int)md_array_size(acc_r));
                double r_nm = probe_radius / nm;
                ImPlot::DragLineX(0, &r_nm, ImVec4(1, 0.6f, 0.2f, 1), 1.0f, ImPlotDragToolFlags_NoInputs);
                ImPlot::EndPlot();
            }

            double rz[2] = { 0.0, 0.0 };
            const bool show_region = region_end > region_beg && (region_beg > 0 || region_end < stats.num_slabs);
            if (show_region) {
                rz[0] = void_profile_z_lo(&p, region_beg) / nm;
                rz[1] = void_profile_z_hi(&p, region_end - 1) / nm;
            }

            ImGui::TableSetColumnIndex(1);
            if (have_rows && ImPlot::BeginPlot("##z_profile", ImVec2(-1, plot_h))) {
                ImPlot::SetupAxes("z (nm)", "Volume fraction");
                ImPlot::SetupAxisLinks(ImAxis_X1, &z_link[0], &z_link[1]);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Once);
                ImPlot::PlotLine("Porosity", prof_z, prof_phi, (int)md_array_size(prof_z));
                ImPlot::PlotLine("V(R,z)/V(z)", prof_z, prof_acc, (int)md_array_size(prof_z));
                if (show_region) {
                    ImPlot::DragLineX(1, &rz[0], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                    ImPlot::DragLineX(2, &rz[1], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                }
                ImPlot::EndPlot();
            }

            ImGui::TableSetColumnIndex(2);
            if (have_rows && ImPlot::BeginPlot("##density_profile", ImVec2(-1, plot_h))) {
                char y_lbl[64];
                snprintf(y_lbl, sizeof(y_lbl), "Density (%s)", density_unit());
                ImPlot::SetupAxes("z (nm)", y_lbl);
                ImPlot::SetupAxisLinks(ImAxis_X1, &z_link[0], &z_link[1]);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Once);
                ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.25f);
                ImPlot::PlotShaded("##density_fill", prof_z, prof_rho, (int)md_array_size(prof_z), 0.0);
                ImPlot::PlotLine("Density", prof_z, prof_rho, (int)md_array_size(prof_z));
                if (show_region) {
                    ImPlot::DragLineX(3, &rz[0], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                    ImPlot::DragLineX(4, &rz[1], ImVec4(1, 0.6f, 0.2f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                }
                // The region's value as a level: the global density is the mean of this curve over
                // the region's slabs, weighted by their volume.
                double rho_mean = rho;
                ImPlot::DragLineY(5, &rho_mean, ImVec4(0.7f, 0.7f, 0.7f, 0.7f), 1.0f, ImPlotDragToolFlags_NoInputs);
                ImPlot::EndPlot();
            }

            ImGui::EndTable();
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

    // The topography of each face as a heat map over (x, y), with the roughness numbers that
    // summarize it. The map is of the probe apex, so it is the surface a probe of the current radius
    // would feel, not the bead surface - they coincide only at R = 0.
    void draw_height_section() {
        const double nm = (double)ANGSTROM_PER_NM;

        if (!ImGui::CollapsingHeader("Surface topography", ImGuiTreeNodeFlags_DefaultOpen)) return;

        if (height_pass) {
            if (draw_pass_progress(height_pass, "Topography")) cancel_height_pass();
        }
        if (!has_height) {
            if (!height_pass) {
                if (ImGui::Button("Compute topography")) start_height();
                ImGui::SetItemTooltip("The height a probe of the current radius finds over every column, from above and\n"
                                      "from below. A pass of its own, about as costly as the field on a sparse network.");
            }
            return;
        }

        if (ImGui::BeginCombo("Map", height_view_lbl[height_view])) {
            for (int i = 0; i < HeightView_Count; ++i) {
                if (ImGui::Selectable(height_view_lbl[i], height_view == i)) {
                    height_view = (HeightView)i;
                    height_disp_dirty = true;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("Top: where a probe lowered from above first touches, as the height of its apex.\n"
                              "Bottom: the same from below. Thickness: the two subtracted, column by column.");

        if (height_disp_dirty) update_height_display();

        ImGui::TextDisabled("At probe %.2f nm, %.2f s over %i x %i columns",
                            height_probe / nm, height_seconds, grid.dim[0], grid.dim[1]);
        if (height_pending && !height_pass) {
            ImGui::SameLine();
            ImGui::TextColored({1.0f, 0.8f, 0.35f, 1.0f}, "(updates on release)");
        }

        const void_heightmap_stats_t& st = height_stats[height_view];
        const uint64_t n_cols = st.num_valid + st.num_open;
        if (st.num_valid > 0) {
            ImGui::Text("Mean %.3f nm, range %.3f to %.3f nm", st.mean / nm, st.min / nm, st.max / nm);
            ImGui::Text("Roughness Rq %.3f nm, Ra %.3f nm, peak to valley %.3f nm",
                        st.rq / nm, st.ra / nm, (st.max - st.min) / nm);
            ImGui::SetItemTooltip("Rq is the RMS deviation from the mean height, Ra the mean absolute deviation.\n"
                                  "Both are of the surface the probe feels, so they fall as the probe grows and\n"
                                  "stops reaching into the narrower dips.");
        }
        if (st.num_open > 0) {
            ImGui::Text("Open columns: %.2f%%", 100.0 * (double)st.num_open / (double)MAX((uint64_t)1, n_cols));
            ImGui::SetItemTooltip("Columns the probe falls straight through without touching anything: a pore that\n"
                                  "runs the full height of the box at this radius. Drawn at the bottom of the\n"
                                  "colour scale and left out of every statistic above.");
        }

        if (height_nx > 0 && height_ny > 0) {
            const double x0 = grid.origin.x / nm;
            const double y0 = grid.origin.y / nm;
            const double x1 = (grid.origin.x + grid.spacing.x * (float)grid.dim[0]) / nm;
            const double y1 = (grid.origin.y + grid.spacing.y * (float)grid.dim[1]) / nm;

            const float h = 320.0f;
            ImPlot::PushColormap(ImPlotColormap_Viridis);
            if (ImPlot::BeginPlot("##height_map", ImVec2(-80, h), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText | ImPlotFlags_Equal)) {
                ImPlot::SetupAxes("x (nm)", "y (nm)");
                ImPlot::SetupAxesLimits(x0, x1, y0, y1, ImPlotCond_Once);
                ImPlot::PlotHeatmap("##height", height_disp, (int)height_ny, (int)height_nx, height_lo, height_hi, nullptr,
                                    ImPlotPoint(x0, y0), ImPlotPoint(x1, y1));

                // Read the full resolution map under the cursor, not the averaged display
                if (ImPlot::IsPlotHovered()) {
                    const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                    const int ix = (int)floor((mp.x * nm - grid.origin.x) / grid.spacing.x);
                    const int iy = (int)floor((mp.y * nm - grid.origin.y) / grid.spacing.y);
                    const float* src = height_source(height_view);
                    if (src && ix >= 0 && iy >= 0 && ix < grid.dim[0] && iy < grid.dim[1]) {
                        const float v = src[(size_t)iy * grid.dim[0] + ix];
                        if (isfinite(v)) {
                            ImGui::SetTooltip("x %.2f, y %.2f nm\n%s %.3f nm", mp.x, mp.y,
                                              height_view == HeightView_Thickness ? "thickness" : "height", v / nm);
                        } else {
                            ImGui::SetTooltip("x %.2f, y %.2f nm\nopen: the probe falls through", mp.x, mp.y);
                        }
                    }
                }
                ImPlot::EndPlot();
            }
            ImGui::SameLine();
            ImPlot::ColormapScale("##height_scale", height_lo, height_hi, ImVec2(70, h), "%.2f");
            ImPlot::PopColormap();
            if (height_nx < (uint32_t)grid.dim[0] || height_ny < (uint32_t)grid.dim[1]) {
                ImGui::TextDisabled("Shown averaged to %u x %u; the statistics and the tooltip use every column.", height_nx, height_ny);
            }
        }
    }

    void draw_volume_section() {
        if (!ImGui::CollapsingHeader("Volume rendering", ImGuiTreeNodeFlags_DefaultOpen)) return;

        ImGui::Checkbox("Show distance field", &show_volume);
        ImGui::SetItemTooltip("Ray casts the distance field in the 3D view. It is transparent wherever the clearance\n"
                              "is below the probe radius, so what remains is exactly where a probe centre of that\n"
                              "size fits, coloured by how much room it has there.");
        if (!show_volume) return;

        if (!has_result || !field) {
            ImGui::TextDisabled("Needs a result computed with 'Materialize field' enabled.");
            return;
        }

        const double nm = (double)ANGSTROM_PER_NM;

        if (ImGui::BeginCombo("Colormap", ImPlot::GetColormapName(vol_colormap))) {
            for (int i = 0; i < ImPlot::GetColormapCount(); ++i) {
                if (ImGui::Selectable(ImPlot::GetColormapName(i), vol_colormap == i)) {
                    vol_colormap = i;
                    vol_tf_dirty = true;
                }
            }
            ImGui::EndCombo();
        }

        if (ImGui::SliderFloat("Opacity", &vol_opacity, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
            vol_tf_dirty = true;
        }
        ImGui::SetItemTooltip("Opacity at the widest clearance shown, per reference step. It falls to a quarter\n"
                              "of that at the probe radius, so the margins of the pore space read as a haze.");

        float upper_nm = vol_upper / (float)nm;
        const float lo_nm = probe_radius / (float)nm;
        const float hi_nm = MAX(lo_nm + 0.01f, (float)(stats.d_max / nm));
        if (ImGui::SliderFloat("Upper limit (nm)", &upper_nm, lo_nm, hi_nm, "%.2f")) {
            vol_upper = upper_nm * (float)nm;
            vol_tf_dirty = true;
        }
        ImGui::SetItemTooltip("Transparent above this clearance too, and the top of the colour scale. The default\n"
                              "is the largest distance in the field, which hides only what no bead was found within\n"
                              "the max distance of - the vacuum around a film, typically.");

        ImGui::Checkbox("Clip to the reporting region", &vol_clip_region);

        ImGui::TextDisabled("Colour spans 0 to %.2f nm, transparent below %.2f nm", vol_tf_max() / nm, probe_radius / nm);
        if (vol_stride > 1) {
            ImGui::TextDisabled("Uploaded at %i x %i x %i, averaged over %i^3 voxels", vol_dim[0], vol_dim[1], vol_dim[2], vol_stride);
        }
    }

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({900, 720}, ImGuiCond_FirstUseEver);
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

            float probe_nm = probe_radius / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Probe radius (nm)", &probe_nm, 0.0f, 10.0f, "%.2f")) {
                probe_radius = probe_nm * ANGSTROM_PER_NM;
                probe_dirty  = true;
                vol_tf_dirty = true;
                if (has_height || height_running()) height_pending = true;
            }
            // The profile and the volume follow the slider; the topography is a geometric pass of its
            // own and waits until the value is let go.
            if (ImGui::IsItemDeactivatedAfterEdit() && (has_height || height_running())) {
                start_height();
            }
            ImGui::SetItemTooltip("Water sized (0.14 nm) for sorption, colloid sized for exclusion.\n"
                                  "The accessible volume, the volume rendering and the topography all follow it.");

            ImGui::SeparatorText("Output");

            ImGui::Checkbox("Materialize field", &materialize);
            ImGui::SetItemTooltip("Keep the full voxel grid in memory. Statistics and the topography do not need it;\n"
                                  "the volume rendering does.");

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
            if (ImGui::Button(field_pass ? "Restart" : "Compute")) {
                compute();
            }
            ImGui::EndDisabled();
            if (field_pass) {
                ImGui::SetItemTooltip("Abandon the running pass and start again with the current settings.");
                ImGui::SameLine();
                if (draw_pass_progress(field_pass, "Distance field")) cancel_field_pass();
            }

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
                draw_height_section();
            }

            draw_volume_section();

            draw_network_section();

#if 0
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
                ImGui::TextDisabled("At probe %.2f nm", sasa.probe / (double)ANGSTROM_PER_NM);

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
#endif
        }
        ImGui::End();
    }
};

static VoidAnalysis instance = {};
