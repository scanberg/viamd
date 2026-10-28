#pragma once

// DENSITY MAPS
// Densities of samples in a plane (the backbone angles of residues, the shapes of structures), estimated off the main
// thread and drawn into ImPlot plots as contours or heat maps, with marks for the samples themselves. What the
// Ramachandran and the Shapespace windows share.
//
// ## Estimation
// A job bins samples on a grid of DENSITY_MAP_DIM^2 texels over a domain, in up to DENSITY_MAP_CHANNELS channels (classes
// of samples, which plots show apart or together), in chunks on the task pool. Each channel is then normalized into a
// probability density of its own, smoothed with a gaussian, and its cumulative distribution tabulated: what density
// encloses a fraction of its samples (density_map_threshold), so that a contour at 90% is the line around 90% of them
// whatever their number or the smoothing (highest density regions). One more channel, DENSITY_MAP_ALL, is all of them:
// the channels weighted by their shares of the samples.
//
// A domain is periodic (angles: what leaves one edge comes in at the other) or bounded (there is nothing outside it). A
// bounded one may have a triangle the samples lie within (the shapes of structures): the smoothing would spread the
// density of samples near its edges out of it, where there are none, and underestimate it there (at a corner of 60
// degrees by up to 5/6). What leaves is reflected back in instead (the reflection method), and the densities and the
// fractions enclosed are of what is within the triangle.
//
// A job is the caller's to keep while it runs; it calls back on the main thread once done, cancelled or not, and is
// freed after. A job cancelled is not to be touched again by the caller (but its callback still comes).
//
// ## Drawing
// The density layers of one plot are composited into one texture of the plot's size (a canvas), drawn again only when
// what they show changes, and put into the plot premultiplied. A layer is drawn as contours (fills of the regions
// above levels, lines along them, of a constant width in pixels) or a heat map, of a weighted sum of the channels of its
// grid, sampled with a cubic B-spline so that contours stay smooth however far the view is zoomed in. A canvas may be
// clipped to a triangle (the triangle of shapes), with antialiased edges.
//
// Marks: points as textured quads (a disc with a dark rim, a ring, a square), thousands of them for the cost of as
// many quads.

#include <core/md_str.h>
#include <task_system.h>

#include <imgui.h>
#include <implot.h>

#include <atomic>
#include <functional>
#include <math.h>
#include <stdint.h>

constexpr int DENSITY_MAP_DIM        = 512;    // Texels along each side of the grid, a power of two
constexpr int DENSITY_MAP_CHANNELS   = 4;
constexpr int DENSITY_MAP_ALL        = DENSITY_MAP_CHANNELS;    // The channel of all of them
constexpr int DENSITY_MAP_TABLE      = 256;    // Entries of the cumulative distribution of a channel
constexpr int DENSITY_MAP_MAX_LEVELS = 8;

struct DensityMapDomain {
    float x_min = 0.0f, x_max = 1.0f;
    float y_min = 0.0f, y_max = 1.0f;
    bool  periodic = false;
    // An equilateral triangle the samples lie within (x, y of each corner, counter-clockwise), bounded domains only:
    // with room around it for the reflections, some 4 sigma, or what they would put there is lost
    bool  triangle = false;
    float tri[6] = {};
};

// The grid of a chunk of a job, that its samples are added to (on the task pool)
struct DensityMapBins {
    float*   data = nullptr;    // DENSITY_MAP_DIM^2 texels of DENSITY_MAP_CHANNELS
    uint64_t samples[DENSITY_MAP_CHANNELS] = {};
    DensityMapDomain domain;
    float    scale_x = 1.0f, scale_y = 1.0f;    // Texels per unit of the domain

    // Of the triangle, if any: the inward normal and offset of each edge (from corner k to k + 1), and how near an edge
    // or a corner a sample is to have images across it
    float    edge_n[3][2] = {};
    float    edge_c[3] = {};
    float    reach = 0.0f;

    // A sample of channel c at (x, y), in the domain's coordinates: wrapped into a periodic domain, dropped outside a
    // bounded one; and its images across the edges and about the corners of the triangle near it, if there is one
    inline void add(int c, float x, float y) {
        if (!put(c, x, y)) return;
        samples[c] += 1;
        if (domain.triangle) add_images(c, x, y);
    }

    inline bool put(int c, float x, float y) {
        int ix = (int)floorf((x - domain.x_min) * scale_x);
        int iy = (int)floorf((y - domain.y_min) * scale_y);
        if (domain.periodic) {
            ix &= DENSITY_MAP_DIM - 1;
            iy &= DENSITY_MAP_DIM - 1;
        } else if (ix < 0 || iy < 0 || ix >= DENSITY_MAP_DIM || iy >= DENSITY_MAP_DIM) {
            return false;
        }
        data[DENSITY_MAP_CHANNELS * (iy * DENSITY_MAP_DIM + ix) + c] += 1.0f;
        return true;
    }

    void add_images(int c, float x, float y);
};

struct DensityMapResult {
    float    table[DENSITY_MAP_CHANNELS + 1][DENSITY_MAP_TABLE] = {};   // See density_map_threshold
    float    max[DENSITY_MAP_CHANNELS + 1] = {};                        // The highest density of each channel
    float    weights_all[DENSITY_MAP_CHANNELS] = {};                    // Of each channel in DENSITY_MAP_ALL
    uint64_t samples[DENSITY_MAP_CHANNELS + 1] = {};
};

// Adds the samples of chunk 'chunk' (of 'num_chunks') to bins, on the task pool. Returns early once cancel is set.
using DensityMapSource = std::function<void(uint32_t chunk, uint32_t num_chunks, DensityMapBins& bins, const std::atomic<int32_t>& cancel)>;

// On the main thread, once the job is done: its result and its grid (DENSITY_MAP_DIM^2 texels of DENSITY_MAP_CHANNELS),
// both NULL when it was cancelled. Valid only during the call.
using DensityMapDone = std::function<void(const DensityMapResult* result, const float* grid)>;

struct DensityMapJob;

// Starts a job. sigma: the standard deviation of the gaussian it is smoothed with, in the domain's units of x (the grid's
// texels are taken to be square). num_chunks: how many parts the samples are binned in, at most what the pool has
// threads for (each takes a grid of its own); the source is called once for each.
DensityMapJob* density_map_launch(str_t label, const DensityMapDomain& domain, float sigma, uint32_t num_chunks, DensityMapSource source, DensityMapDone done);

// Its result is not wanted: done is called with NULL. With wait, returns once the source is no longer called (what it
// reads may then go).
void density_map_cancel(DensityMapJob* job, bool wait);

// The density at which the region of everything denser encloses the fraction 'level' of the samples of a channel
float density_map_threshold(const float table[DENSITY_MAP_TABLE], float level);

// ## CPU helpers, exposed to be tested

// Fills the cumulative distribution of a density (every stride'th value of v, of those mask has nonzero if given) and
// its maximum
void density_map_mass_table(float table[DENSITY_MAP_TABLE], float* out_max, const float* v, size_t count, size_t stride, const uint8_t* mask = nullptr);

// The bins of a domain (data and samples left to the caller), with what a triangle's images need for a sigma (domain
// units)
void density_map_bins_init(DensityMapBins* bins, const DensityMapDomain& domain, float sigma);

// Of the texels whose centres are within the domain's triangle, 1 (DENSITY_MAP_DIM^2 of them)
void density_map_triangle_mask(uint8_t* mask, const DensityMapDomain& domain);

// A gaussian blur of a dim x dim grid of DENSITY_MAP_CHANNELS, in place: periodic, or with nothing outside it
void density_map_blur(float* data, int dim, float sigma_texels, bool periodic);

// ## GL

// A texture for the grid of a job: wrapped where the domain is periodic, nothing (zero) outside it where it is bounded
uint32_t density_map_texture_create(bool periodic);
void     density_map_texture_upload(uint32_t tex, const float* grid);
void     density_map_texture_free(uint32_t* tex);

// What one layer is drawn with in one plot. Zeroed whole (density_map_draw_init) before it is filled: a canvas is drawn
// again only when the bytes of its draws change.
struct DensityMapDraw {
    uint32_t tex;
    int32_t  mode;          // 0: contours, 1: heat map
    int32_t  fill;
    int32_t  lines;
    int32_t  num_levels;
    int32_t  cmap_count;
    int32_t  cmap_log;
    float    line_width;    // Pixels
    float    weights[DENSITY_MAP_CHANNELS];
    float    levels[DENSITY_MAP_MAX_LEVELS];       // Densities, the highest first
    float    fill_colors[DENSITY_MAP_MAX_LEVELS][4];
    float    line_color[4];
    float    cmap[16][4];
    float    cmap_max;
    float    opacity;
    uint64_t content;       // Of what is in tex: another result, another drawing
};

void density_map_draw_init(DensityMapDraw* d);
// Contours at the levels (fractions of the samples enclosed, increasing) of a channel of a result, or of the channels
// weighted as in it for DENSITY_MAP_ALL
void density_map_draw_levels(DensityMapDraw* d, const DensityMapResult& r, int channel, const float* levels, int num_levels);
// Fills of the regions between levels, each a share of color so that they add up to most of it at the middle
void density_map_draw_fills_from(DensityMapDraw* d, ImVec4 color);
// A heat map through an ImPlot colour map: up to 16 of its colours, or 16 samples of it where it has more
void density_map_draw_colormap(DensityMapDraw* d, ImPlotColormap cmap, bool log_scale);

// The composite of the density layers of one plot
struct DensityMapCanvas {
    uint32_t tex = 0;
    int      width = 0, height = 0;
    uint64_t hash = 0;
};

// Composites the draws, the first at the bottom, into the canvas, of width x height texels, when what they show
// changed. view_uv: the uv of the grid at the plot's lower left corner, and the extent of the plot in uv. clip_uv, if
// given: a triangle (3 uv) outside which nothing is drawn. Returns the texture, 0 when there is nothing to show.
uint32_t density_map_composite(DensityMapCanvas* canvas, int width, int height, const float view_uv[4], const DensityMapDraw* draws, int num_draws, const float* clip_uv = nullptr);
void     density_map_canvas_free(DensityMapCanvas* canvas);
void     density_map_canvas_invalidate(DensityMapCanvas* canvas);

// The canvas over the whole of the current plot (between BeginPlot and EndPlot, after SetupFinish)
void density_map_plot_canvas(uint32_t tex);

// The size in texels a canvas for the current plot is to have (its size times the framebuffer's scale)
void density_map_plot_canvas_size(int* width, int* height);

// The uv rect (see density_map_composite) the current plot shows of a domain
void density_map_plot_view_uv(float out_view_uv[4], const DensityMapDomain& domain);

// Frees what is shared by everything drawn (the program, the marks). Safe to call more than once.
void density_map_shutdown();

// ## Marks

enum DensityMapMark : int {
    DensityMapMark_Disc = 0,    // Tinted, with a dark rim
    DensityMapMark_Ring,        // Tinted: around what is selected or highlighted
    DensityMapMark_Square,      // Plain: many points
};

// Between begin and end, marks of radius r (pixels) at p into the draw list (the plot's, within its clip rect)
void density_map_marks_begin(ImDrawList* dl);
void density_map_mark(ImDrawList* dl, ImVec2 p, float r, DensityMapMark mark, ImU32 col);
void density_map_marks_end(ImDrawList* dl);

// ## UI shared by the windows of layers

// A row of buttons of which one is on. With fill, they share the width left on the row.
bool density_map_segmented(const char* id, const char* const* labels, const int* values, int count, int* value, const char* const* tooltips = nullptr, bool fill = false);

// The levels of contours as fractions enclosed, a field for each (in %, dragged, right-click to remove) and + to add
// one. Kept sorted. Returns true when they change.
bool density_map_levels_editor(float* levels, int* num_levels);

// A swatch of how a layer looks, within [min, max]: points, contours (filled or not), or a heat map
enum DensityMapSwatch : int {
    DensityMapSwatch_Points = 0,
    DensityMapSwatch_Contours,
    DensityMapSwatch_Heatmap,
};
void density_map_swatch(ImDrawList* dl, ImVec2 min, ImVec2 max, DensityMapSwatch kind, ImVec4 color, bool fill, ImPlotColormap cmap, float alpha);
