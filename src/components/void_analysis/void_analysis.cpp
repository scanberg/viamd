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

#include "channels.h"

#include <float.h>

/*
    Void and nanopore characterization: the distance field stage.

    For every voxel of a regular grid this computes the additively weighted distance to the nearest bead surface,

        d(p) = min_i ( |p - c_i| - R_i ),

    which is the radius of the largest probe sphere that fits at p. The work is done by the weighted nearest query
    on md_spatial_acc, so the field is exact everywhere rather than banded, and periodicity is whatever the unit cell
    of the system says it is - periodic in x and y and open in z for a supported film, for instance.

    Everything downstream of the field - thresholding, connected components, the percolation threshold radius, pore
    size distributions - is deliberately absent here. This stage is what those need, and it is what the analytic
    test cases in mdlib validate.

    Residency: statistics are accumulated per tile while the tile is still in cache, so the full field is only
    materialized when explicitly asked for. At 0.5 nm voxels a 670 x 670 x 175 nm box is 2.4e9 voxels, which is
    9.6 GB as float - streaming is the default for a reason.
*/

namespace {

constexpr int   TILE_DIM  = 8;                        // Voxels per tile along each axis
constexpr int   TILE_SIZE = TILE_DIM * TILE_DIM * TILE_DIM;
constexpr int   NUM_BINS  = 256;
constexpr int   SASA_MAX_POINTS = 1024;
constexpr float ANGSTROM_PER_NM = 10.0f;
constexpr double PI_D = 3.14159265358979323846;

// Angstrom^2 per Dalton to m^2 per gram
constexpr double SPECIFIC_AREA_SCALE = 1.0e-20 / 1.66053906892e-24;

// Opacity of the channel ribbon. The overlay pass blends on source alpha, so this is how much of the
// structure is still read through the tube.
constexpr uint32_t RIBBON_ALPHA = 0x60;

enum RadiusSource {
    RadiusSource_Vdw = 0,       // Per atom radius as reported by the system
    RadiusSource_Uniform,       // One radius for every bead, for coarse grained models with no meaningful element
    RadiusSource_Count,
};

const char* radius_source_lbl[RadiusSource_Count] = {
    "Van der Waals",
    "Uniform",
};

struct Stats {
    uint64_t num_voxels    = 0;
    uint64_t num_solid     = 0;      // Voxels inside a bead, i.e. d <= 0
    double   d_min         = 0.0;
    double   d_max         = 0.0;
    md_array(uint64_t) histogram = 0; // Void voxels binned over [0, hist_max]
    double   hist_max      = 0.0;
    double   bin_width     = 0.0;
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

    // Channels
    float channel_r_min     = 2.0f;   // Radii below this are not worth sweeping for a channel
    int   channel_num_radii = 32;
    bool  has_channels      = false;
    double r_c              = 0.0;    // Largest probe radius which still gets through
    double channel_seconds  = 0.0;
    channel_tree_t channels = {};
    md_array(double)   curve_radius = 0;
    md_array(uint32_t) curve_count  = 0;
    md_array(float)    node_slot    = 0;   // Dendrogram column per branch, negative when not drawn
    md_array(double)   plot_x       = 0;   // Scratch for the channel count plot, in nm
    md_array(double)   plot_y       = 0;
    float    num_slots     = 0.0f;
    uint32_t hovered_node  = CHANNEL_INVALID_INDEX;

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
                if (!show_window || !has_channels || hovered_node == CHANNEL_INVALID_INDEX) break;
                if (e.payload_type != viamd::EventPayloadType_ApplicationState) break;
                const ApplicationState& state = *(ApplicationState*)e.payload;
                immediate::Scope scope(state.gfx.overlay, "void_channel_path");
                // Camera Z in world space. The overlay is rendered with an identity model matrix, so
                // this is the same space the grid is in.
                const vec3_t cam_axis = vec3_normalize(vec3_from_vec4(state.view.param.matrix.inv.view.col[2]));
                draw_channel_path_3d(scope, hovered_node, cam_axis);
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
        md_array_free(stats.histogram, arena);
        stats = {};
        has_result = false;
        clear_channels();
    }

    void clear_channels() {
        channel_tree_free(&channels);
        md_array_free(curve_radius, arena);
        md_array_free(curve_count, arena);
        md_array_free(node_slot, arena);
        md_array_free(plot_x, arena);
        md_array_free(plot_y, arena);
        plot_x = 0;
        plot_y = 0;
        curve_radius = 0;
        curve_count  = 0;
        node_slot    = 0;
        num_slots    = 0.0f;
        hovered_node = CHANNEL_INVALID_INDEX;
        has_channels = false;
        r_c = 0.0;
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
        double A[3][3];
        md_unitcell_A_extract_double(A, &state.unitcell);
        const uint32_t flags = md_unitcell_flags(&state.unitcell);

        vec3_t origin = {0, 0, 0};
        vec3_t extent = {0, 0, 0};

        if (flags != MD_UNITCELL_NONE) {
            // Cartesian bounds of the cell parallelepiped: the sum of the positive components of the basis vectors
            for (int r = 0; r < 3; ++r) {
                double lo = 0.0, hi = 0.0;
                for (int c = 0; c < 3; ++c) {
                    const double v = A[c][r];
                    if (v < 0.0) lo += v; else hi += v;
                }
                origin.elem[r] = (float)lo;
                extent.elem[r] = (float)(hi - lo);
            }
        } else {
            vec3_t aabb_min = { FLT_MAX,  FLT_MAX,  FLT_MAX};
            vec3_t aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
            for (size_t i = 0; i < state.num_atoms; ++i) {
                const vec3_t p = {state.x[i], state.y[i], state.z[i]};
                aabb_min = vec3_min(aabb_min, p);
                aabb_max = vec3_max(aabb_max, p);
            }
            const float margin = 2.0f * voxel_spacing;
            for (int a = 0; a < 3; ++a) {
                origin.elem[a] = aabb_min.elem[a] - margin;
                extent.elem[a] = (aabb_max.elem[a] - aabb_min.elem[a]) + 2.0f * margin;
            }
        }

        for (int a = 0; a < 3; ++a) {
            const int dim = (int)round((double)extent.elem[a] / (double)voxel_spacing);
            out_grid->dim[a] = MAX(1, dim);
            out_grid->spacing.elem[a] = extent.elem[a] / (float)out_grid->dim[a];
        }
        out_grid->origin = origin;
        out_grid->orientation = mat3_ident();

        return md_grid_num_points(out_grid) > 0;
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

        const md_timestamp_t t0 = md_time_now();

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

        const int tiles[3] = {
            (grid.dim[0] + TILE_DIM - 1) / TILE_DIM,
            (grid.dim[1] + TILE_DIM - 1) / TILE_DIM,
            (grid.dim[2] + TILE_DIM - 1) / TILE_DIM,
        };
        const uint32_t num_tiles = (uint32_t)tiles[0] * (uint32_t)tiles[1] * (uint32_t)tiles[2];

        // Per thread accumulators, merged once the range task has completed. Nothing is shared while it runs.
        const size_t num_threads = MAX((size_t)1, task_system::pool_num_threads() + 1);
        uint64_t* hist   = (uint64_t*)md_temp_alloc(temp_scope, num_threads * NUM_BINS * sizeof(uint64_t));
        uint64_t* solid  = (uint64_t*)md_temp_alloc(temp_scope, num_threads * sizeof(uint64_t));
        float*    d_lo   = (float*)   md_temp_alloc(temp_scope, num_threads * sizeof(float));
        float*    d_hi   = (float*)   md_temp_alloc(temp_scope, num_threads * sizeof(float));
        MEMSET(hist,  0, num_threads * NUM_BINS * sizeof(uint64_t));
        MEMSET(solid, 0, num_threads * sizeof(uint64_t));
        for (size_t i = 0; i < num_threads; ++i) {
            d_lo[i] =  FLT_MAX;
            d_hi[i] = -FLT_MAX;
        }

        // Never bin finer than the voxel sampling. A bin narrower than a voxel only splits the same samples into
        // noisier buckets, and the histogram density is exactly what the surface area estimate below reads.
        const double voxel_min  = (double)MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        const double bin_width  = MAX(voxel_min, (double)max_dist / (double)NUM_BINS);
        const double hist_max   = bin_width * (double)NUM_BINS;
        const double bin_scale  = 1.0 / bin_width;

        task_system::ID task = task_system::create_pool_task(STR_LIT("Void distance field"), num_tiles,
            [&](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                const uint32_t ti = MIN((uint32_t)(num_threads - 1), thread_num);
                uint64_t* my_hist = hist + (size_t)ti * NUM_BINS;

                float qx[TILE_SIZE], qy[TILE_SIZE], qz[TILE_SIZE];
                float dist[TILE_SIZE];
                int   vi[TILE_SIZE], vj[TILE_SIZE], vk[TILE_SIZE];

                for (uint32_t t = range_beg; t < range_end; ++t) {
                    const int tx = (int)(t % (uint32_t)tiles[0]);
                    const int ty = (int)((t / (uint32_t)tiles[0]) % (uint32_t)tiles[1]);
                    const int tz = (int)(t / ((uint32_t)tiles[0] * (uint32_t)tiles[1]));

                    // A tile is a compact block of voxels, which is what the batched query is efficient at
                    int n = 0;
                    for (int k = 0; k < TILE_DIM; ++k) {
                        const int z = tz * TILE_DIM + k;
                        if (z >= grid.dim[2]) break;
                        for (int j = 0; j < TILE_DIM; ++j) {
                            const int y = ty * TILE_DIM + j;
                            if (y >= grid.dim[1]) break;
                            for (int i = 0; i < TILE_DIM; ++i) {
                                const int x = tx * TILE_DIM + i;
                                if (x >= grid.dim[0]) break;
                                qx[n] = grid.origin.x + ((float)x + 0.5f) * grid.spacing.x;
                                qy[n] = grid.origin.y + ((float)y + 0.5f) * grid.spacing.y;
                                qz[n] = grid.origin.z + ((float)z + 0.5f) * grid.spacing.z;
                                vi[n] = x;
                                vj[n] = y;
                                vk[n] = z;
                                n += 1;
                            }
                        }
                    }
                    if (n == 0) continue;

                    md_coord_stream_t pts = md_coord_stream_from_soa(qx, qy, qz, NULL, (size_t)n);
                    md_spatial_acc_query_nearest(&acc, &pts, (double)max_dist, NULL, dist);

                    for (int p = 0; p < n; ++p) {
                        const float d = dist[p];
                        d_lo[ti] = MIN(d_lo[ti], d);
                        d_hi[ti] = MAX(d_hi[ti], d);

                        if (d <= 0.0f) {
                            solid[ti] += 1;
                        } else {
                            int bin = (int)((double)d * bin_scale);
                            bin = CLAMP(bin, 0, NUM_BINS - 1);
                            my_hist[bin] += 1;
                        }

                        if (field) {
                            // Output contract: unsigned distance in the void, clamped inside the solid
                            const size_t idx = ((size_t)vk[p] * grid.dim[1] + (size_t)vj[p]) * grid.dim[0] + (size_t)vi[p];
                            field[idx] = MAX(0.0f, d);
                        }
                    }
                }
            }, 1);

        task_system::enqueue_task(task);
        task_system::task_wait_for(task);

        md_array_resize(stats.histogram, NUM_BINS, arena);
        MEMSET(stats.histogram, 0, NUM_BINS * sizeof(uint64_t));

        stats.num_voxels = num_voxels;
        stats.num_solid  = 0;
        stats.d_min      =  DBL_MAX;
        stats.d_max      = -DBL_MAX;
        for (size_t t = 0; t < num_threads; ++t) {
            stats.num_solid += solid[t];
            if (d_lo[t] !=  FLT_MAX) stats.d_min = MIN(stats.d_min, (double)d_lo[t]);
            if (d_hi[t] != -FLT_MAX) stats.d_max = MAX(stats.d_max, (double)d_hi[t]);
            for (int b = 0; b < NUM_BINS; ++b) {
                stats.histogram[b] += hist[t * NUM_BINS + b];
            }
        }
        if (stats.d_min == DBL_MAX) {
            stats.d_min = 0.0;
            stats.d_max = 0.0;
        }
        stats.hist_max     = hist_max;
        stats.bin_width    = bin_width;
        stats.voxel_volume = (double)grid.spacing.x * (double)grid.spacing.y * (double)grid.spacing.z;
        stats.seconds      = md_time_as_seconds(md_time_now() - t0);

        has_result = true;
    }

    // Fraction of the box a probe of the supplied radius can occupy the centre of. This is a volume measure only,
    // it says nothing about whether that volume is reachable from outside - that needs the connectivity stage.
    double accessible_fraction(double r_probe) const {
        if (!has_result || stats.num_voxels == 0) return 0.0;
        const double bin_width = stats.bin_width;
        uint64_t count = 0;
        for (int b = 0; b < NUM_BINS; ++b) {
            if ((double)(b + 1) * bin_width > r_probe) count += stats.histogram[b];
        }
        return (double)count / (double)stats.num_voxels;
    }

    // A(r) = -dV_acc/dr, and with |grad d| = 1 almost everywhere the coarea formula makes that derivative the
    // density of the distance histogram at r. One field therefore already carries the accessible surface area at
    // every probe radius, for the cost of a division.
    //
    // Read at bin centres and interpolate: taking the bin which contains r instead biases the answer by half a bin,
    // and A(r) is not flat. It is still a discretized estimate - expect a few percent, against a few tenths of a
    // percent for the sampling estimator below - so it earns its place as an independent cross check rather than as
    // the number to quote. The two share almost no code path.
    double coarea_area(double r) const {
        if (!has_result || stats.bin_width <= 0.0) return 0.0;
        const double u  = r / stats.bin_width - 0.5;
        const int    b0 = (int)floor(u);
        const double t  = u - (double)b0;
        const double d0 = (b0     >= 0 && b0     < NUM_BINS) ? (double)stats.histogram[b0]     : 0.0;
        const double d1 = (b0 + 1 >= 0 && b0 + 1 < NUM_BINS) ? (double)stats.histogram[b0 + 1] : 0.0;
        return ((1.0 - t) * d0 + t * d1) * stats.voxel_volume / stats.bin_width;
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

        const md_timestamp_t t0 = md_time_now();

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
        sasa.seconds = md_time_as_seconds(md_time_now() - t0);

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

    void compute_channels() {
        clear_channels();
        error[0] = '\0';

        channel_field_t f;
        if (!build_channel_field(&f)) {
            snprintf(error, sizeof(error), "Compute the field first, with 'Materialize field' enabled");
            return;
        }

        const md_timestamp_t t0 = md_time_now();

        const double r_lo  = (double)channel_r_min;
        const double d_max = MAX(stats.d_max, r_lo + 1.0);
        const double tol   = 0.05 * (double)MIN(grid.spacing.x, MIN(grid.spacing.y, grid.spacing.z));
        const size_t n     = (size_t)CLAMP(channel_num_radii, 2, 128);

        // r_c first, because it is also where the count curve has to stop. Nothing percolates above
        // it by definition, so sweeping up to the largest clearance in the box - which measures a
        // cavity, not a throat, and is just the max_dist clamp wherever a voxel has no bead in range
        // - spends nearly every sample on a radius which cannot possibly get through, and leaves the
        // one region that carries information sampled more coarsely than the voxel grid.
        r_c = channel_critical_radius(&f, r_lo, d_max, tol, md_get_heap_allocator());

        // A little past r_c so the curve is seen to fall to zero rather than ending on a cliff.
        const double r_hi = (r_c >= r_lo) ? MIN(d_max, MAX(1.1 * r_c, r_lo + 10.0 * tol)) : d_max;

        md_array_resize(curve_radius, n, arena);
        md_array_resize(curve_count,  n, arena);
        for (size_t i = 0; i < n; ++i) {
            curve_radius[i] = r_lo + (r_hi - r_lo) * (double)i / (double)(n - 1);
        }

        // Each radius is an independent sweep, so this is a plain fan out
        task_system::ID task = task_system::create_pool_task(STR_LIT("Channel percolation"), (uint32_t)n,
            [&](uint32_t beg, uint32_t end, uint32_t) {
                channel_spanning_counts(curve_count + beg, curve_radius + beg, end - beg, &f, md_get_heap_allocator());
            }, 1);
        task_system::enqueue_task(task);
        task_system::task_wait_for(task);

        channel_sweep(&channels, &f, (double)probe_radius, true, arena);
        layout_channel_tree();

        channel_seconds = md_time_as_seconds(md_time_now() - t0);
        has_channels = true;
    }

    void layout_channel_tree() {
        const size_t num = md_array_size(channels.nodes);
        md_array_resize(node_slot, num, arena);
        num_slots = channel_tree_layout(node_slot, &channels, arena);
    }

    // The hovered branch, then the route it takes on down to the bottom. Drawn into the overlay
    // queue, which is rendered last with the depth test off, so it stays in front of the structure
    // rather than being buried by it.
    //
    // A line is one pixel wide at any zoom, which in a dense network is easy to lose and thin enough
    // for temporal AA to eat, so each segment is also drawn as a camera facing ribbon whose half
    // width is the clearance there. The ribbon carries the width of the channel and is translucent so
    // the structure still reads through it; the line stays opaque on top of it, which is what keeps a
    // tight channel visible when the ribbon is below a pixel across.
    void draw_channel_path_3d(immediate::Queue* q, uint32_t node, vec3_t cam_axis) const {
        const float wrap_x = 0.5f * grid.spacing.x * (float)grid.dim[0];
        const float wrap_y = 0.5f * grid.spacing.y * (float)grid.dim[1];

        auto segment = [&](vec4_t a, vec4_t b, uint32_t col) {
            // A channel may wrap through the periodic faces; the centreline must not be drawn
            // straight back across the box when it does.
            if (fabsf(a.x - b.x) > wrap_x || fabsf(a.y - b.y) > wrap_y) return;

            const vec3_t p0 = {a.x, a.y, a.z};
            const vec3_t p1 = {b.x, b.y, b.z};

            const vec3_t d   = vec3_sub(p1, p0);
            const float  len = vec3_length(d);
            if (len > 1.0e-6f) {
                // Perpendicular to both the segment and the line of sight. It only degenerates when
                // the segment points straight at the camera, where a ribbon has nothing to show.
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
            ImGui::SetItemTooltip("Distances are resolved out to this range; beyond it a voxel is reported at the limit.");

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
                ImGui::SeparatorText("Result");
                ImGui::Text("Computed in %.2f s", stats.seconds);
                ImGui::Text("Solid fraction: %.4f", (double)stats.num_solid / (double)stats.num_voxels);
                ImGui::Text("Distance range: %.2f to %.2f nm", stats.d_min / ANGSTROM_PER_NM, stats.d_max / ANGSTROM_PER_NM);

                ImGui::Text("Accessible volume fraction: %.4f", accessible_fraction(probe_radius));
                ImGui::TextDisabled("Volume only - reachability from outside needs the connectivity stage.");

                if (ImPlot::BeginPlot("##Void distance", ImVec2(-1, 200))) {
                    const double bin_width = stats.hist_max / (double)NUM_BINS;
                    double centers[NUM_BINS];
                    double values[NUM_BINS];
                    for (int b = 0; b < NUM_BINS; ++b) {
                        centers[b] = ((double)b + 0.5) * bin_width / (double)ANGSTROM_PER_NM;
                        values[b]  = (double)stats.histogram[b] / (double)stats.num_voxels;
                    }
                    ImPlot::SetupAxes("Distance to nearest surface (nm)", "Volume fraction");
                    ImPlot::PlotBars("Void", centers, values, NUM_BINS, bin_width / (double)ANGSTROM_PER_NM);
                    ImPlot::EndPlot();
                }
            }

            ImGui::SeparatorText("Channels through z");

            ImGui::TextDisabled("Sweeps along z, which is treated as open at both ends.");

            float r_min_nm = channel_r_min / ANGSTROM_PER_NM;
            if (ImGui::SliderFloat("Smallest radius (nm)", &r_min_nm, 0.05f, 5.0f, "%.2f")) {
                channel_r_min = r_min_nm * ANGSTROM_PER_NM;
            }
            ImGui::SetItemTooltip("The radius sweep starts here. Below it a network is usually one\n"
                                  "connected pore and the channel count stops meaning much.");

            ImGui::SliderInt("Radius samples", &channel_num_radii, 4, 128);

            ImGui::BeginDisabled(!has_result || !field);
            if (ImGui::Button("Find channels")) {
                compute_channels();
            }
            ImGui::EndDisabled();
            if (!field) {
                ImGui::SameLine();
                ImGui::TextDisabled("(needs a materialized field)");
            }

            if (has_channels) {
                ImGui::Text("Computed in %.2f s", channel_seconds);
                if (r_c >= channel_r_min) {
                    ImGui::Text("Critical radius r_c: %.3f nm", r_c / (double)ANGSTROM_PER_NM);
                    ImGui::SetItemTooltip("The largest probe which still gets from one z face to the other.\n"
                                          "Set by the tightest throat along the best route, not by any cavity on it.");
                } else {
                    ImGui::Text("Nothing gets through at %.2f nm or above", channel_r_min / ANGSTROM_PER_NM);
                }
                ImGui::Text("Channels at probe %.2f nm: %u", probe_radius / ANGSTROM_PER_NM, channels.num_spanning);
                ImGui::SetItemTooltip("Independent channels: two routes which join anywhere along the way\n"
                                      "are one channel, not two.");

                if (md_array_size(curve_radius) > 1 && ImPlot::BeginPlot("##Channel count", ImVec2(-1, 140))) {
                    const size_t n = md_array_size(curve_radius);
                    ImPlot::SetupAxes("Probe radius (nm)", "Channels");
                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                    md_array_resize(plot_x, n, arena);
                    md_array_resize(plot_y, n, arena);
                    for (size_t i = 0; i < n; ++i) {
                        plot_x[i] = curve_radius[i] / (double)ANGSTROM_PER_NM;
                        plot_y[i] = (double)curve_count[i];
                    }
                    ImPlot::PlotStairs("Channels", plot_x, plot_y, (int)n);
                    if (r_c >= channel_r_min) {
                        double rc_nm = r_c / (double)ANGSTROM_PER_NM;
                        ImPlot::DragLineX(0, &rc_nm, ImVec4(1, 0.6f, 0.2f, 1), 1.0f, ImPlotDragToolFlags_NoInputs);
                    }
                    ImPlot::EndPlot();
                }

                ImGui::TextDisabled("Sweep tree at probe %.2f nm - hover a branch to trace it in 3D",
                                    probe_radius / ANGSTROM_PER_NM);
                draw_channel_tree(ImVec2(ImGui::GetContentRegionAvail().x, 220.0f));
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
