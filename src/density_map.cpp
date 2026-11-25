#include <density_map.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_vec_math.h>
#include <core/md_hash.h>
#include <core/md_log.h>

#include <gfx/profiler.h>
#include <app/imgui_impl_md_gpu.h>

#include "density_map_shaders.inl"

#include <implot_internal.h>

#include <math.h>
#include <string.h>

// ## Estimation

static constexpr int MAX_CHUNKS = 8;   // Each holds a grid of its own (16 MB of DENSITY_MAP_DIM = 512 and 4 channels)

static size_t grid_bytes() { return sizeof(float) * DENSITY_MAP_CHANNELS * DENSITY_MAP_DIM * DENSITY_MAP_DIM; }

// Entry k of the cumulative distribution of a density is the density where the region of everything denser encloses
// mass_level(k) of it. The levels crowd towards 1 (up to 99.9%), where contours of rare samples are drawn.
static inline double mass_level(int k) {
    return 1.0 - pow(10.0, -3.0 * k / (DENSITY_MAP_TABLE - 1));
}

float density_map_threshold(const float table[DENSITY_MAP_TABLE], float level) {
    level = MD_CLAMP(level, 0.0f, 0.999f);
    const double x = -log10(1.0 - level) / 3.0 * (DENSITY_MAP_TABLE - 1);
    const int k = (int)x;
    if (k >= DENSITY_MAP_TABLE - 1) return table[DENSITY_MAP_TABLE - 1];
    return (float)lerp((double)table[k], (double)table[k + 1], x - k);
}

// The mass sorted into bins over 8 decades below the maximum, walked from the densest
void density_map_mass_table(float table[DENSITY_MAP_TABLE], float* out_max, const float* v, size_t count, size_t stride, const uint8_t* mask) {
    float vmax = 0.0f;
    double total = 0.0;
    for (size_t i = 0; i < count; ++i) {
        if (mask && !mask[i]) continue;
        const float x = v[i * stride];
        vmax = MD_MAX(vmax, x);
        total += x;
    }
    *out_max = vmax;
    if (vmax <= 0.0f || total <= 0.0) {
        md_memset(table, 0, sizeof(float) * DENSITY_MAP_TABLE);
        return;
    }

    constexpr int NUM_BINS = 4096;
    constexpr double DECADES = 8.0;
    const double bins_per_decade = NUM_BINS / DECADES;
    md_allocator_i* heap = md_get_heap_allocator();
    double* hist = (double*)md_alloc(heap, sizeof(double) * NUM_BINS);
    md_memset(hist, 0, sizeof(double) * NUM_BINS);
    const double inv_max = 1.0 / vmax;
    for (size_t i = 0; i < count; ++i) {
        if (mask && !mask[i]) continue;
        const float x = v[i * stride];
        if (x <= 0.0f) continue;
        const int b = (int)(-log10(x * inv_max) * bins_per_decade);
        hist[MD_CLAMP(b, 0, NUM_BINS - 1)] += x;
    }

    table[0] = vmax;
    int b = 0;
    double cum = 0.0;
    for (int k = 1; k < DENSITY_MAP_TABLE; ++k) {
        const double m = mass_level(k) * total;
        while (b < NUM_BINS && cum + hist[b] < m) {
            cum += hist[b];
            ++b;
        }
        if (b >= NUM_BINS) {
            table[k] = (float)(vmax * pow(10.0, -DECADES));
            continue;
        }
        const double frac = hist[b] > 0.0 ? (m - cum) / hist[b] : 0.0;
        table[k] = (float)(vmax * pow(10.0, -(b + frac) / bins_per_decade));
    }
    md_free(heap, hist, sizeof(double) * NUM_BINS);
}

// A running box of half-width r along each row: wrapped around a periodic grid, or with nothing outside a bounded one
static void blur_rows(md_float4_t* out, const md_float4_t* in, int dim, int r, bool periodic) {
    const int mask = dim - 1;
    const float scl = 1.0f / (2 * r + 1);
    for (int row = 0; row < dim; ++row) {
        const md_float4_t* src = in + dim * row;
        md_float4_t* dst = out + dim * row;
        md_float4_t acc = md_float4_zero();
        if (periodic) {
            for (int x = -r; x <= r; ++x) acc = acc + src[x & mask];
            for (int x = 0; x < dim; ++x) {
                dst[x] = acc * scl;
                acc = md_float4_max(md_float4_zero(), acc - src[(x - r) & mask] + src[(x + r + 1) & mask]);
            }
        } else {
            for (int x = 0; x <= r && x < dim; ++x) acc = acc + src[x];
            for (int x = 0; x < dim; ++x) {
                dst[x] = acc * scl;
                if (x + r + 1 < dim) acc = acc + src[x + r + 1];
                if (x - r >= 0)      acc = acc - src[x - r];
                acc = md_float4_max(md_float4_zero(), acc);
            }
        }
    }
}

static void transpose(md_float4_t* dst, const md_float4_t* src, int n) {
    const int block = 8;
    MD_ASSERT(n % block == 0);
    for (int i = 0; i < n; i += block) {
        for (int j = 0; j < n; j += block) {
            for (int k = i; k < i + block; ++k) {
                for (int l = j; l < j + block; ++l) {
                    dst[k + l * n] = src[l + k * n];
                }
            }
        }
    }
}

// The half-widths of n boxes whose succession approximates a gaussian of sigma
static void boxes_for_gauss(int* half_w, int n, float sigma) {
    const float w_ideal = sqrtf((12 * sigma * sigma / n) + 1);
    int wl = (int)w_ideal;
    if (wl % 2 == 0) wl--;
    const int wu = wl + 2;
    const float m_ideal = (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
    const int m = (int)(m_ideal + 0.5f);
    for (int i = 0; i < n; i++) {
        const int w = i < m ? wl : wu;
        half_w[i] = MD_MAX(0, (w - 1) / 2);
    }
}

void density_map_blur(float* data, int dim, float sigma, bool periodic) {
    MD_STATIC_ASSERT(DENSITY_MAP_CHANNELS == 4, "the blur works on md_float4_t");
    MD_ASSERT(dim > 0 && (dim & (dim - 1)) == 0);
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_float4_t* grid = (md_float4_t*)data;
    md_float4_t* tmp = (md_float4_t*)md_temp_alloc_aligned(temp, dim * dim * sizeof(md_float4_t), 16);

    int r[3];
    boxes_for_gauss(r, 3, sigma);
    for (int pass = 0; pass < 2; ++pass) {
        blur_rows(tmp, grid, dim, r[0], periodic);
        blur_rows(grid, tmp, dim, r[1], periodic);
        blur_rows(tmp, grid, dim, r[2], periodic);
        transpose(grid, tmp, dim);
    }
}

// ## The triangle's images
//
// The reflection method: the density within a region of samples near its boundary is estimated as that of the samples
// and of their mirror images across the boundary, which puts back what the kernels of the samples spread out of it.
// At an edge, the image across it. At a corner of 60 degrees, the images of the wedge's group (of the six symmetries of
// the corner's two edges and the line through it at 60 degrees from both): the edges' two, the rotations by 120 degrees
// either way about the corner, and the reflection across the third line, which is parallel to the opposite edge.

void density_map_bins_init(DensityMapBins* bins, const DensityMapDomain& domain, float sigma) {
    bins->domain = domain;
    bins->scale_x = DENSITY_MAP_DIM / (domain.x_max - domain.x_min);
    bins->scale_y = DENSITY_MAP_DIM / (domain.y_max - domain.y_min);
    if (!domain.triangle) return;
    const float* t = domain.tri;
    // Counter-clockwise or not, the normals inward: towards the corner opposite each edge
    for (int k = 0; k < 3; ++k) {
        const float ax = t[2 * k], ay = t[2 * k + 1];
        const float bx = t[2 * ((k + 1) % 3)], by = t[2 * ((k + 1) % 3) + 1];
        const float ox = t[2 * ((k + 2) % 3)], oy = t[2 * ((k + 2) % 3) + 1];
        float nx = -(by - ay), ny = bx - ax;
        const float len = sqrtf(nx * nx + ny * ny);
        nx /= MD_MAX(len, 1e-12f);
        ny /= MD_MAX(len, 1e-12f);
        if (nx * (ox - ax) + ny * (oy - ay) < 0.0f) {
            nx = -nx;
            ny = -ny;
        }
        bins->edge_n[k][0] = nx;
        bins->edge_n[k][1] = ny;
        bins->edge_c[k] = nx * ax + ny * ay;
    }
    // Beyond 4 sigma, an image puts less than 0.04% of the kernel's peak within the triangle
    bins->reach = 4.0f * sigma;
}

void DensityMapBins::add_images(int c, float x, float y) {
    const float* t = domain.tri;
    float d[3];
    for (int k = 0; k < 3; ++k) {
        d[k] = edge_n[k][0] * x + edge_n[k][1] * y - edge_c[k];
        if (d[k] < reach) put(c, x - 2.0f * d[k] * edge_n[k][0], y - 2.0f * d[k] * edge_n[k][1]);
    }
    constexpr float cs = -0.5f, sn = 0.86602540378f;    // 120 degrees
    for (int v = 0; v < 3; ++v) {
        const float vx = t[2 * v], vy = t[2 * v + 1];
        const float px = x - vx, py = y - vy;
        if (px * px + py * py >= reach * reach) continue;
        put(c, vx + cs * px - sn * py, vy + sn * px + cs * py);
        put(c, vx + cs * px + sn * py, vy - sn * px + cs * py);
        // The edge opposite corner v is edge v + 1 (from corner v + 1 to v + 2)
        const float* n = edge_n[(v + 1) % 3];
        const float e = n[0] * px + n[1] * py;
        put(c, x - 2.0f * e * n[0], y - 2.0f * e * n[1]);
    }
}

void density_map_triangle_mask(uint8_t* mask, const DensityMapDomain& domain) {
    DensityMapBins b;
    density_map_bins_init(&b, domain, 0.0f);
    const float dx = (domain.x_max - domain.x_min) / DENSITY_MAP_DIM;
    const float dy = (domain.y_max - domain.y_min) / DENSITY_MAP_DIM;
    for (int iy = 0; iy < DENSITY_MAP_DIM; ++iy) {
        const float y = domain.y_min + (iy + 0.5f) * dy;
        for (int ix = 0; ix < DENSITY_MAP_DIM; ++ix) {
            const float x = domain.x_min + (ix + 0.5f) * dx;
            bool in = true;
            for (int k = 0; k < 3; ++k) in = in && (b.edge_n[k][0] * x + b.edge_n[k][1] * y - b.edge_c[k] >= 0.0f);
            mask[iy * DENSITY_MAP_DIM + ix] = in ? 1 : 0;
        }
    }
}

struct DensityMapJob {
    DensityMapDomain domain;
    float    sigma = 0.0f;                  // Domain units
    float    sigma_texels = 1.0f;
    uint32_t num_chunks = 0;
    uint32_t num_bufs = 0;              // The chunks', and at least one more for the channel of all (job_finalize)
    float*   bins[MAX_CHUNKS] = {};
    uint64_t samples[MAX_CHUNKS][DENSITY_MAP_CHANNELS] = {};

    DensityMapSource source;
    DensityMapDone   done;

    std::atomic<int32_t>  cancel{0};
    std::atomic<uint32_t> chunks_done{0};
    bool ok = false;
    DensityMapResult result;

    // Bins the chunks and merges them on the pool, then hands the result over on the main thread, where the job
    // ends itself (job_free)
    task_system::Job task;
};

static void job_free(DensityMapJob* job) {
    md_allocator_i* heap = md_get_heap_allocator();
    for (uint32_t i = 0; i < job->num_bufs; ++i) {
        if (job->bins[i]) md_free(heap, job->bins[i], grid_bytes());
    }
    delete job;
}

static void job_bin_chunk(DensityMapJob* job, uint32_t c) {
    DensityMapBins bins;
    bins.data = job->bins[c];
    density_map_bins_init(&bins, job->domain, job->sigma);
    job->source(c, job->num_chunks, bins, job->cancel);
    md_memcpy(job->samples[c], bins.samples, sizeof(bins.samples));
    if (!job->cancel.load()) job->chunks_done.fetch_add(1);
}

// Merges the chunks, normalizes each channel, smooths, and finds the distribution of each
static void job_finalize(DensityMapJob* job) {
    if (job->cancel.load() || job->chunks_done.load() != job->num_chunks) return;

    const size_t num_texels = (size_t)DENSITY_MAP_DIM * DENSITY_MAP_DIM;
    const int C = DENSITY_MAP_CHANNELS;
    float* d = job->bins[0];
    for (uint32_t c = 1; c < job->num_chunks; ++c) {
        const float* s = job->bins[c];
        for (size_t i = 0; i < num_texels * C; ++i) d[i] += s[i];
    }

    DensityMapResult& r = job->result;
    uint64_t total = 0;
    for (int t = 0; t < C; ++t) {
        for (uint32_t c = 0; c < job->num_chunks; ++c) r.samples[t] += job->samples[c][t];
        total += r.samples[t];
    }
    r.samples[DENSITY_MAP_ALL] = total;

    // Each channel a probability density of its own: plots of the channels compare their shapes whatever their counts
    float inv[C];
    for (int t = 0; t < C; ++t) {
        inv[t] = r.samples[t] ? 1.0f / (float)r.samples[t] : 0.0f;
        r.weights_all[t] = total ? (float)((double)r.samples[t] / (double)total) : 0.0f;
    }
    for (size_t i = 0; i < num_texels; ++i) {
        for (int t = 0; t < C; ++t) d[C * i + t] *= inv[t];
    }

    if (job->cancel.load()) return;
    density_map_blur(d, DENSITY_MAP_DIM, job->sigma_texels, job->domain.periodic);
    if (job->cancel.load()) return;

    // Of what is within the triangle only, if there is one: outside it are the images
    uint8_t* mask = nullptr;
    md_allocator_i* heap = md_get_heap_allocator();
    if (job->domain.triangle) {
        mask = (uint8_t*)md_alloc(heap, num_texels);
        density_map_triangle_mask(mask, job->domain);
    }
    for (int t = 0; t < C; ++t) {
        density_map_mass_table(r.table[t], &r.max[t], d + t, num_texels, C, mask);
    }
    // All of them: the channels weighted by their shares, in the second buffer (merged into the first already)
    MD_ASSERT(job->num_bufs >= 2);
    float* all = job->bins[1];
    for (size_t i = 0; i < num_texels; ++i) {
        float v = 0.0f;
        for (int t = 0; t < C; ++t) v += d[C * i + t] * r.weights_all[t];
        all[i] = v;
    }
    density_map_mass_table(r.table[DENSITY_MAP_ALL], &r.max[DENSITY_MAP_ALL], all, num_texels, 1, mask);
    if (mask) md_free(heap, mask, num_texels);

    job->ok = true;
}

DensityMapJob* density_map_launch(md_str_t label, const DensityMapDomain& domain, float sigma, uint32_t num_chunks, DensityMapSource source, DensityMapDone done) {
    MD_ASSERT(source && done);
    md_allocator_i* heap = md_get_heap_allocator();
    DensityMapJob* job = new DensityMapJob();
    job->domain = domain;
    job->domain.triangle = domain.triangle && !domain.periodic;
    job->sigma = sigma;
    job->sigma_texels = MD_MAX(0.25f, sigma * DENSITY_MAP_DIM / MD_MAX(1e-12f, domain.x_max - domain.x_min));
    job->source = std::move(source);
    job->done = std::move(done);

    const uint32_t max_chunks = (uint32_t)MD_CLAMP(task_system::num_threads(), 1u, (uint32_t)MAX_CHUNKS);
    job->num_chunks = MD_CLAMP(num_chunks, 1u, max_chunks);
    job->num_bufs = MD_MAX(job->num_chunks, 2u);
    for (uint32_t i = 0; i < job->num_bufs; ++i) {
        job->bins[i] = (float*)md_alloc(heap, grid_bytes());
        md_memset(job->bins[i], 0, grid_bytes());
    }

    job->task.start(label, [job]() {
        task_system::parallel_for(job->num_chunks, 1, [job](uint32_t beg, uint32_t end) {
            for (uint32_t c = beg; c < end; ++c) job_bin_chunk(job, c);
        });
        job_finalize(job);
    }, [job](task_system::Status status) {
        const bool ok = status == task_system::Status::Completed && job->ok && !job->cancel.load();
        job->done(ok ? &job->result : nullptr, ok ? job->bins[0] : nullptr);
        job_free(job);
    });
    return job;
}

void density_map_cancel(DensityMapJob* job, bool wait) {
    if (!job) return;
    job->cancel.store(1);
    job->task.cancel();
    if (wait) job->task.wait();
}

// ## GPU

// density_map.slang's parameters of a draw (PARAMS float4, at P_*) and Args
struct GpuDraw {
    md_gpu_float4 weights;
    md_gpu_float4 levels[2];
    md_gpu_float4 fill[DENSITY_MAP_MAX_LEVELS];
    md_gpu_float4 line_color;
    md_gpu_float4 cmap[16];
    uint32_t mode, fill_on, lines_on, num_levels;
    uint32_t cmap_count, cmap_log, _p0, _p1;
    float    line_width, cmap_max, opacity, _p2;
    md_gpu_float4 _p3;
};
static_assert(sizeof(GpuDraw) == 32 * 16, "density_map.slang's PARAMS");

struct GpuArgs {
    md_gpu_float4 uv_rect, clip01, clip2;
    uint32_t width, height, count, periodic;
    md_gpu_addr_t params;
    md_gpu_addr_t grids;        // md_gpu_sampled_tex_t per draw
    md_gpu_storage_tex_t dst;
    uint32_t _p0, _p1;
};
static_assert(sizeof(GpuArgs) == 96, "density_map.slang's Args");

static struct {
    md_gpu_stream_t  stream = nullptr;
    md_gpu_kernel_t  kernel = nullptr;
    md_gpu_texture_t marks = nullptr;
} gpu;

void density_map_initialize(md_gpu_stream_t stream) {
    gpu.stream = stream;
}

md_gpu_texture_t density_map_texture_create() {
    if (!gpu.stream) return nullptr;
    md_gpu_texture_desc_t d = {};
    d.type   = MD_GPU_TEX_2D;
    d.format = MD_GPU_FORMAT_RGBA32_FLOAT;
    d.usage  = MD_GPU_TEX_SAMPLED;
    d.width  = DENSITY_MAP_DIM;
    d.height = DENSITY_MAP_DIM;
    d.label  = "density map";
    md_gpu_texture_t tex = md_gpu_texture_create(gpu.stream, &d);
    if (!tex) MD_LOG_ERROR("A density map's texture could not be made: %s", md_gpu_last_error());
    return tex;
}

void density_map_texture_upload(md_gpu_texture_t tex, const float* grid) {
    if (!tex || !gpu.stream) return;
    if (!md_gpu_upload_texture(gpu.stream, tex, nullptr, grid, grid_bytes())) {
        MD_LOG_ERROR("A density map could not be uploaded: %s", md_gpu_last_error());
    }
}

void density_map_texture_free(md_gpu_texture_t* tex) {
    if (tex && *tex) {
        md_gpu_texture_destroy(*tex);
        *tex = nullptr;
    }
}

void density_map_shutdown() {
    if (gpu.kernel) md_gpu_kernel_destroy(gpu.kernel);
    if (gpu.marks)  md_gpu_texture_destroy(gpu.marks);
    gpu.kernel = nullptr;
    gpu.marks  = nullptr;
}

void density_map_draw_init(DensityMapDraw* d) {
    md_memset(d, 0, sizeof(DensityMapDraw));
    d->opacity = 1.0f;
    d->line_width = 1.5f;
}

void density_map_draw_levels(DensityMapDraw* d, const DensityMapResult& r, int channel, const float* levels, int num_levels) {
    MD_ASSERT(0 <= channel && channel <= DENSITY_MAP_ALL);
    md_memset(d->weights, 0, sizeof(d->weights));
    if (channel == DENSITY_MAP_ALL) {
        for (int t = 0; t < DENSITY_MAP_CHANNELS; ++t) d->weights[t] = r.weights_all[t];
    } else {
        d->weights[channel] = 1.0f;
    }
    d->num_levels = MD_CLAMP(num_levels, 0, DENSITY_MAP_MAX_LEVELS);
    // Levels increase in what they enclose, so their densities decrease: the densest first
    for (int i = 0; i < d->num_levels; ++i) d->levels[i] = density_map_threshold(r.table[channel], levels[i]);
    d->cmap_max = r.max[channel];
}

void density_map_draw_fills_from(DensityMapDraw* d, ImVec4 color) {
    const int n = MD_MAX(1, d->num_levels);
    const float a = 1.0f - powf(0.2f, 1.0f / (float)n);
    for (int i = 0; i < DENSITY_MAP_MAX_LEVELS; ++i) {
        d->fill_colors[i][0] = color.x;
        d->fill_colors[i][1] = color.y;
        d->fill_colors[i][2] = color.z;
        d->fill_colors[i][3] = i < d->num_levels ? color.w * a : 0.0f;
    }
}

void density_map_draw_colormap(DensityMapDraw* d, ImPlotColormap cmap, bool log_scale) {
    const int size = ImPlot::GetColormapSize(cmap);
    const int n = MD_MIN(size, 16);
    d->cmap_count = n;
    d->cmap_log = log_scale;
    for (int i = 0; i < n; ++i) {
        const ImVec4 c = size <= 16 ? ImPlot::GetColormapColor(i, cmap) : ImPlot::SampleColormap((float)i / (n - 1), cmap);
        d->cmap[i][0] = c.x;
        d->cmap[i][1] = c.y;
        d->cmap[i][2] = c.z;
        d->cmap[i][3] = c.w;
    }
}

static md_gpu_float4 f4(float x, float y, float z, float w) {
    md_gpu_float4 v = {x, y, z, w};
    return v;
}

md_gpu_texture_t density_map_composite(DensityMapCanvas* canvas, int width, int height, const float view_uv[4], const DensityMapDraw* draws, int num_draws,
                                       bool periodic, const float* clip_uv) {
    MD_ASSERT(canvas);
    if (num_draws <= 0 || !gpu.stream) return nullptr;
    if (!gpu.kernel) {
        md_gpu_kernel_desc_t kd = viamd_shader_density_map_composite_kernel();
        gpu.kernel = md_gpu_kernel_create(md_gpu_stream_device(gpu.stream), &kd);
        if (!gpu.kernel) {
            MD_LOG_ERROR("The density maps' kernel could not be made: %s", md_gpu_last_error());
            return nullptr;
        }
    }

    width  = MD_CLAMP(width,  1, 4096);
    height = MD_CLAMP(height, 1, 4096);
    if (!canvas->tex || canvas->width != width || canvas->height != height) {
        if (canvas->tex) md_gpu_texture_destroy(canvas->tex);
        md_gpu_texture_desc_t d = {};
        d.type   = MD_GPU_TEX_2D;
        d.format = MD_GPU_FORMAT_RGBA8_UNORM;
        d.usage  = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED;
        d.width  = (uint32_t)width;
        d.height = (uint32_t)height;
        d.label  = "density map canvas";
        canvas->tex = md_gpu_texture_create(gpu.stream, &d);
        canvas->width = width;
        canvas->height = height;
        canvas->hash = 0;
        if (!canvas->tex) {
            MD_LOG_ERROR("A density map's canvas (%d x %d) could not be made: %s", width, height, md_gpu_last_error());
            *canvas = {};
            return nullptr;
        }
    }

    uint64_t h = md_hash64(draws, sizeof(DensityMapDraw) * num_draws, (uint64_t)(uintptr_t)canvas->tex);
    h = md_hash64(view_uv, sizeof(float) * 4, h ^ (periodic ? 0x9E71 : 0));
    if (clip_uv) h = md_hash64(clip_uv, sizeof(float) * 6, h ^ 0xC11B);
    if (h == canvas->hash) return canvas->tex;
    canvas->hash = h;

    PUSH_GPU_SECTION("Density map");
    md_gpu_temp_t scope = md_gpu_temp_begin(gpu.stream);
    GpuDraw* g = (GpuDraw*)md_alloc(md_get_heap_allocator(), sizeof(GpuDraw) * num_draws);
    md_gpu_sampled_tex_t* grids = (md_gpu_sampled_tex_t*)md_alloc(md_get_heap_allocator(), sizeof(md_gpu_sampled_tex_t) * num_draws);
    md_memset(g, 0, sizeof(GpuDraw) * num_draws);
    for (int i = 0; i < num_draws; ++i) {
        const DensityMapDraw& d = draws[i];
        GpuDraw& o = g[i];
        o.weights = f4(d.weights[0], d.weights[1], d.weights[2], d.weights[3]);
        o.levels[0] = f4(d.levels[0], d.levels[1], d.levels[2], d.levels[3]);
        o.levels[1] = f4(d.levels[4], d.levels[5], d.levels[6], d.levels[7]);
        for (int k = 0; k < DENSITY_MAP_MAX_LEVELS; ++k) o.fill[k] = f4(d.fill_colors[k][0], d.fill_colors[k][1], d.fill_colors[k][2], d.fill_colors[k][3]);
        o.line_color = f4(d.line_color[0], d.line_color[1], d.line_color[2], d.line_color[3]);
        for (int k = 0; k < 16; ++k) o.cmap[k] = f4(d.cmap[k][0], d.cmap[k][1], d.cmap[k][2], d.cmap[k][3]);
        o.mode       = (uint32_t)d.mode;
        o.fill_on    = d.fill ? 1u : 0u;
        o.lines_on   = d.lines ? 1u : 0u;
        o.num_levels = (uint32_t)MD_CLAMP(d.num_levels, 0, DENSITY_MAP_MAX_LEVELS);
        o.cmap_count = (uint32_t)MD_CLAMP(d.cmap_count, 1, 16);
        o.cmap_log   = d.cmap_log ? 1u : 0u;
        o.line_width = d.line_width;
        o.cmap_max   = d.cmap_max;
        o.opacity    = d.opacity;
        grids[i]     = d.tex ? md_gpu_texture_sampled(d.tex) : md_gpu_sampled_tex_t{};
        if (!d.tex) o.opacity = 0.0f;
    }
    // A layer without a texture is drawn as nothing: a texture of the others stands in for it
    md_gpu_sampled_tex_t any = {};
    for (int i = 0; i < num_draws; ++i) if (grids[i].handle) any = grids[i];
    for (int i = 0; i < num_draws; ++i) if (!grids[i].handle) grids[i] = any;

    GpuArgs a = {};
    a.uv_rect  = f4(view_uv[0], view_uv[1], view_uv[2], view_uv[3]);
    a.clip01   = clip_uv ? f4(clip_uv[0], clip_uv[1], clip_uv[2], clip_uv[3]) : f4(0, 0, 1, 0);
    a.clip2    = clip_uv ? f4(clip_uv[4], clip_uv[5], 1, 0) : f4(0, 1, 0, 0);
    a.width    = (uint32_t)width;
    a.height   = (uint32_t)height;
    a.count    = any.handle ? (uint32_t)num_draws : 0u;
    a.periodic = periodic ? 1u : 0u;
    a.params   = md_gpu_temp_upload(scope, g, sizeof(GpuDraw) * num_draws);
    a.grids    = md_gpu_temp_upload(scope, grids, sizeof(md_gpu_sampled_tex_t) * num_draws);
    a.dst      = md_gpu_texture_storage(canvas->tex, 0);
    md_free(md_get_heap_allocator(), g, sizeof(GpuDraw) * num_draws);
    md_free(md_get_heap_allocator(), grids, sizeof(md_gpu_sampled_tex_t) * num_draws);
    const md_gpu_ordering_t prev = md_gpu_stream_set_ordering(gpu.stream, MD_GPU_ORDER_IMPLICIT);
    const bool ok = a.params && a.grids && md_gpu_launch(gpu.stream, gpu.kernel, md_gpu_grid_for(gpu.kernel, a.width, a.height, 1), &a, sizeof(a));
    md_gpu_stream_set_ordering(gpu.stream, prev);
    md_gpu_temp_end(scope);
    if (!ok) {
        MD_LOG_ERROR("A density map could not be drawn: %s", md_gpu_last_error());
        canvas->hash = 0;
    }
    POP_GPU_SECTION();

    return canvas->tex;
}

void density_map_canvas_free(DensityMapCanvas* canvas) {
    if (canvas && canvas->tex) md_gpu_texture_destroy(canvas->tex);
    if (canvas) *canvas = {};
}

void density_map_canvas_invalidate(DensityMapCanvas* canvas) {
    if (canvas) canvas->hash = 0;
}

void density_map_plot_canvas(md_gpu_texture_t tex) {
    if (!tex) return;
    const ImVec2 pos = ImPlot::GetPlotPos();
    const ImVec2 size = ImPlot::GetPlotSize();
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlot::PushPlotClipRect();
    ImGui_ImplMdGpu_PushState(dl, ImGui_ImplMdGpu_Premultiplied);
    dl->AddImage(ImGui_ImplMdGpu_TextureID(tex), pos, ImVec2(pos.x + size.x, pos.y + size.y), ImVec2(0, 0), ImVec2(1, 1));
    ImGui_ImplMdGpu_PopState(dl);
    ImPlot::PopPlotClipRect();
}

void density_map_plot_canvas_size(int* width, int* height) {
    // Of the viewport the window is in, where the backend says (multiple viewports), else of the display
    ImVec2 scale = ImGui::GetWindowViewport()->FramebufferScale;
    if (scale.x <= 0.0f || scale.y <= 0.0f) scale = ImGui::GetIO().DisplayFramebufferScale;
    const ImVec2 size = ImPlot::GetPlotSize();
    *width  = (int)(size.x * scale.x + 0.5f);
    *height = (int)(size.y * scale.y + 0.5f);
}

void density_map_plot_view_uv(float out[4], const DensityMapDomain& domain) {
    const ImPlotRect lim = ImPlot::GetPlotLimits();
    const double w = domain.x_max - domain.x_min;
    const double h = domain.y_max - domain.y_min;
    out[0] = (float)((lim.X.Min - domain.x_min) / w);
    out[1] = (float)((lim.Y.Min - domain.y_min) / h);
    out[2] = (float)(lim.Size().x / w);
    out[3] = (float)(lim.Size().y / h);
}

// ## Marks

// Side by side: a disc with a dark rim (tinted, the rim staying dark) and a ring (tinted)
// The marks' texture, with its mips (each a box filter of the one above it): what a mark a few pixels across samples
static md_gpu_texture_t create_marks() {
    constexpr int S = 64;
    static uint32_t pixels[S * 2 * S];
    for (int y = 0; y < S; ++y) {
        for (int x = 0; x < 2 * S; ++x) {
            const int cx = x < S ? S / 2 : S + S / 2;
            const float dx = (x + 0.5f - cx) / (S * 0.5f - 2.0f);
            const float dy = (y + 0.5f - S / 2) / (S * 0.5f - 2.0f);
            const float r = sqrtf(dx * dx + dy * dy);
            const float aa = 1.5f / (S * 0.5f);
            float lum = 1.0f, alpha = 0.0f;
            if (x < S) {
                const float rim = MD_CLAMP((r - 0.72f) / aa + 0.5f, 0.0f, 1.0f);
                lum = 1.0f - rim;
                alpha = MD_CLAMP((1.0f - r) / aa + 0.5f, 0.0f, 1.0f) * (1.0f - 0.25f * rim);
            } else {
                alpha = MD_CLAMP((1.0f - r) / aa + 0.5f, 0.0f, 1.0f) * MD_CLAMP((r - 0.7f) / aa + 0.5f, 0.0f, 1.0f);
            }
            const uint32_t l8 = (uint32_t)(lum * 255.0f + 0.5f);
            const uint32_t a8 = (uint32_t)(alpha * 255.0f + 0.5f);
            pixels[y * 2 * S + x] = l8 | (l8 << 8) | (l8 << 16) | (a8 << 24);
        }
    }
    if (!gpu.stream) return nullptr;
    constexpr uint32_t MIPS = 8;    // 128 x 64 down to 1 x 1
    md_gpu_texture_desc_t d = {};
    d.type       = MD_GPU_TEX_2D;
    d.format     = MD_GPU_FORMAT_RGBA8_UNORM;
    d.usage      = MD_GPU_TEX_SAMPLED;
    d.width      = 2 * S;
    d.height     = S;
    d.mip_levels = MIPS;
    d.label      = "density map marks";
    md_gpu_texture_t tex = md_gpu_texture_create(gpu.stream, &d);
    if (!tex) {
        MD_LOG_ERROR("The marks' texture could not be made: %s", md_gpu_last_error());
        return nullptr;
    }
    static uint32_t mip[S * S];     // the next level, at most half of the first
    uint32_t* src = pixels;
    uint32_t w = 2 * S, h = S;
    for (uint32_t m = 0; m < MIPS; ++m) {
        md_gpu_tex_region_t r = {};
        r.extent[0] = w;
        r.extent[1] = h;
        r.extent[2] = 1;
        r.mip = m;
        md_gpu_upload_texture(gpu.stream, tex, &r, src, (size_t)w * h * 4);
        const uint32_t nw = MD_MAX(w / 2, 1u), nh = MD_MAX(h / 2, 1u);
        for (uint32_t y = 0; y < nh; ++y) {
            for (uint32_t x = 0; x < nw; ++x) {
                uint32_t sum[4] = {};
                for (uint32_t k = 0; k < 4; ++k) {
                    const uint32_t sx = MD_MIN(2 * x + (k & 1), w - 1), sy = MD_MIN(2 * y + (k >> 1), h - 1);
                    const uint32_t p = src[sy * w + sx];
                    for (int c = 0; c < 4; ++c) sum[c] += (p >> (8 * c)) & 0xFF;
                }
                uint32_t o = 0;
                for (int c = 0; c < 4; ++c) o |= ((sum[c] + 2) / 4) << (8 * c);
                mip[y * nw + x] = o;
            }
        }
        if (src == pixels) src = mip;   // in place from here on: each texel reads only texels at or after its own
        w = nw;
        h = nh;
    }
    return tex;
}

void density_map_marks_begin(ImDrawList* dl) {
    if (!gpu.marks) gpu.marks = create_marks();
    dl->PushTexture(ImTextureRef(ImGui_ImplMdGpu_TextureID(gpu.marks)));
}

void density_map_mark(ImDrawList* dl, ImVec2 p, float r, DensityMapMark mark, ImU32 col) {
    static const ImVec2 uv[3][2] = {
        { ImVec2(0.0f,  0.0f), ImVec2(0.5f,  1.0f) },   // Disc
        { ImVec2(0.5f,  0.0f), ImVec2(1.0f,  1.0f) },   // Ring
        { ImVec2(0.25f, 0.5f), ImVec2(0.25f, 0.5f) },   // Square: the middle of the disc
    };
    dl->PrimReserve(6, 4);
    dl->PrimRectUV(ImVec2(p.x - r, p.y - r), ImVec2(p.x + r, p.y + r), uv[mark][0], uv[mark][1], col);
}

void density_map_marks_end(ImDrawList* dl) {
    dl->PopTexture();
}

// ## UI

bool density_map_segmented(const char* id, const char* const* labels, const int* values, int count, int* value, const char* const* tooltips, bool fill) {
    bool changed = false;
    ImGui::PushID(id);
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(1, ImGui::GetStyle().ItemSpacing.y));
    const float w = fill ? MD_MAX(1.0f, (ImGui::GetContentRegionAvail().x - (count - 1)) / count) : 0.0f;
    for (int i = 0; i < count; ++i) {
        const bool on = *value == values[i];
        if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
        if (ImGui::Button(labels[i], ImVec2(w, 0))) {
            *value = values[i];
            changed = true;
        }
        if (on) ImGui::PopStyleColor();
        if (tooltips && tooltips[i]) ImGui::SetItemTooltip("%s", tooltips[i]);
        if (i + 1 < count) ImGui::SameLine();
    }
    ImGui::PopStyleVar();
    ImGui::PopID();
    return changed;
}

static void sort_levels(float* levels, int n) {
    for (int i = 1; i < n; ++i) {
        for (int j = i; j > 0 && levels[j - 1] > levels[j]; --j) {
            const float tmp = levels[j];
            levels[j] = levels[j - 1];
            levels[j - 1] = tmp;
        }
    }
}

bool density_map_levels_editor(float* levels, int* num_levels) {
    bool changed = false;
    int n = *num_levels;
    const float w = ImGui::CalcTextSize("99.9%").x + ImGui::GetStyle().FramePadding.x * 2;
    int remove = -1;
    ImGui::PushID("##levels");
    for (int i = 0; i < n; ++i) {
        ImGui::PushID(i);
        if (i > 0) ImGui::SameLine(0, 2);
        if (ImGui::GetContentRegionAvail().x < w) ImGui::NewLine();
        float pct = levels[i] * 100.0f;
        ImGui::SetNextItemWidth(w);
        if (ImGui::DragFloat("##lvl", &pct, 0.2f, 1.0f, 99.9f, pct < 99.0f ? "%.0f%%" : "%.1f%%", ImGuiSliderFlags_AlwaysClamp)) {
            levels[i] = pct / 100.0f;
            changed = true;
        }
        if (ImGui::IsItemDeactivatedAfterEdit()) sort_levels(levels, n);
        if (ImGui::BeginPopupContextItem("##lvl_menu")) {
            if (ImGui::MenuItem("Remove", nullptr, false, n > 1)) remove = i;
            ImGui::EndPopup();
        }
        ImGui::SetItemTooltip("Drag to change, right-click to remove");
        ImGui::PopID();
    }
    ImGui::SameLine(0, 2);
    ImGui::BeginDisabled(n >= DENSITY_MAP_MAX_LEVELS);
    if (ImGui::Button("+##add_level")) {
        // Halfway between the outermost and all of them
        const float last = n > 0 ? levels[n - 1] : 0.5f;
        levels[n++] = MD_MIN(0.999f, last + (1.0f - last) * 0.5f);
        changed = true;
    }
    ImGui::EndDisabled();
    if (remove >= 0 && n > 1) {
        for (int i = remove; i + 1 < n; ++i) levels[i] = levels[i + 1];
        n -= 1;
        changed = true;
    }
    ImGui::PopID();
    *num_levels = n;
    return changed;
}

void density_map_swatch(ImDrawList* dl, ImVec2 mn, ImVec2 mx, DensityMapSwatch kind, ImVec4 color, bool fill, ImPlotColormap cmap, float alpha) {
    const float s = mx.y - mn.y;
    const ImVec2 mid = ImVec2((mn.x + mx.x) * 0.5f, (mn.y + mx.y) * 0.5f);
    ImVec4 c = color;
    c.w *= alpha;
    const ImU32 col = ImGui::ColorConvertFloat4ToU32(c);
    switch (kind) {
    case DensityMapSwatch_Points:
        dl->AddCircleFilled(ImVec2(mn.x + s * 0.3f, mn.y + s * 0.35f), 2.5f, col);
        dl->AddCircleFilled(ImVec2(mn.x + s * 0.7f, mn.y + s * 0.45f), 2.5f, col);
        dl->AddCircleFilled(ImVec2(mn.x + s * 0.45f, mn.y + s * 0.75f), 2.5f, col);
        break;
    case DensityMapSwatch_Contours:
        if (fill) {
            ImVec4 f = c;
            f.w *= 0.45f;
            dl->AddCircleFilled(mid, s * 0.46f, ImGui::ColorConvertFloat4ToU32(f));
        }
        dl->AddCircle(mid, s * 0.46f, col, 0, 1.5f);
        dl->AddCircle(mid, s * 0.22f, col, 0, 1.5f);
        break;
    case DensityMapSwatch_Heatmap: {
        const int n = 5;
        for (int k = 0; k < n; ++k) {
            ImVec4 cc = ImPlot::SampleColormap((float)k / (n - 1), cmap);
            cc.w *= alpha;
            const float x0 = mn.x + (mx.x - mn.x) * k / n;
            const float x1 = mn.x + (mx.x - mn.x) * (k + 1) / n;
            dl->AddRectFilled(ImVec2(x0, mn.y), ImVec2(x1, mx.y), ImGui::ColorConvertFloat4ToU32(cc));
        }
        break;
    }
    default:
        break;
    }
}
