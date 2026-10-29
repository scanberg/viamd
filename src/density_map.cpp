#include <density_map.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_vec_math.h>
#include <core/md_hash.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/profiler.h>

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
    level = CLAMP(level, 0.0f, 0.999f);
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
        vmax = MAX(vmax, x);
        total += x;
    }
    *out_max = vmax;
    if (vmax <= 0.0f || total <= 0.0) {
        MEMSET(table, 0, sizeof(float) * DENSITY_MAP_TABLE);
        return;
    }

    constexpr int NUM_BINS = 4096;
    constexpr double DECADES = 8.0;
    const double bins_per_decade = NUM_BINS / DECADES;
    md_allocator_i* heap = md_get_heap_allocator();
    double* hist = (double*)md_alloc(heap, sizeof(double) * NUM_BINS);
    MEMSET(hist, 0, sizeof(double) * NUM_BINS);
    const double inv_max = 1.0 / vmax;
    for (size_t i = 0; i < count; ++i) {
        if (mask && !mask[i]) continue;
        const float x = v[i * stride];
        if (x <= 0.0f) continue;
        const int b = (int)(-log10(x * inv_max) * bins_per_decade);
        hist[CLAMP(b, 0, NUM_BINS - 1)] += x;
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
static void blur_rows(vec4_t* out, const vec4_t* in, int dim, int r, bool periodic) {
    const int mask = dim - 1;
    const float scl = 1.0f / (2 * r + 1);
    for (int row = 0; row < dim; ++row) {
        const vec4_t* src = in + dim * row;
        vec4_t* dst = out + dim * row;
        vec4_t acc = vec4_zero();
        if (periodic) {
            for (int x = -r; x <= r; ++x) acc = acc + src[x & mask];
            for (int x = 0; x < dim; ++x) {
                dst[x] = acc * scl;
                acc = vec4_max(vec4_zero(), acc - src[(x - r) & mask] + src[(x + r + 1) & mask]);
            }
        } else {
            for (int x = 0; x <= r && x < dim; ++x) acc = acc + src[x];
            for (int x = 0; x < dim; ++x) {
                dst[x] = acc * scl;
                if (x + r + 1 < dim) acc = acc + src[x + r + 1];
                if (x - r >= 0)      acc = acc - src[x - r];
                acc = vec4_max(vec4_zero(), acc);
            }
        }
    }
}

static void transpose(vec4_t* dst, const vec4_t* src, int n) {
    const int block = 8;
    ASSERT(n % block == 0);
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
        half_w[i] = MAX(0, (w - 1) / 2);
    }
}

void density_map_blur(float* data, int dim, float sigma, bool periodic) {
    STATIC_ASSERT(DENSITY_MAP_CHANNELS == 4, "the blur works on vec4_t");
    ASSERT(dim > 0 && (dim & (dim - 1)) == 0);
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    vec4_t* grid = (vec4_t*)data;
    vec4_t* tmp = (vec4_t*)md_temp_alloc_aligned(temp, dim * dim * sizeof(vec4_t), 16);

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
        nx /= MAX(len, 1e-12f);
        ny /= MAX(len, 1e-12f);
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

    task_system::ID bin_task = task_system::INVALID_ID;
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
    MEMCPY(job->samples[c], bins.samples, sizeof(bins.samples));
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
    ASSERT(job->num_bufs >= 2);
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

DensityMapJob* density_map_launch(str_t label, const DensityMapDomain& domain, float sigma, uint32_t num_chunks, DensityMapSource source, DensityMapDone done) {
    ASSERT(source && done);
    md_allocator_i* heap = md_get_heap_allocator();
    DensityMapJob* job = new DensityMapJob();
    job->domain = domain;
    job->domain.triangle = domain.triangle && !domain.periodic;
    job->sigma = sigma;
    job->sigma_texels = MAX(0.25f, sigma * DENSITY_MAP_DIM / MAX(1e-12f, domain.x_max - domain.x_min));
    job->source = std::move(source);
    job->done = std::move(done);

    const uint32_t max_chunks = (uint32_t)CLAMP(task_system::pool_num_threads(), (size_t)1, (size_t)MAX_CHUNKS);
    job->num_chunks = CLAMP(num_chunks, 1u, max_chunks);
    job->num_bufs = MAX(job->num_chunks, 2u);
    for (uint32_t i = 0; i < job->num_bufs; ++i) {
        job->bins[i] = (float*)md_alloc(heap, grid_bytes());
        MEMSET(job->bins[i], 0, grid_bytes());
    }

    task_system::ID bin_task = task_system::create_pool_task(label, job->num_chunks, [job](uint32_t beg, uint32_t end, uint32_t) {
        for (uint32_t c = beg; c < end; ++c) job_bin_chunk(job, c);
    });
    task_system::ID finalize_task = task_system::create_pool_task(STR_LIT("##Density finalize"), [job]() {
        job_finalize(job);
    });
    task_system::ID done_task = task_system::create_main_task(STR_LIT("##Density done"), [job]() {
        const bool ok = job->ok && !job->cancel.load();
        job->done(ok ? &job->result : nullptr, ok ? job->bins[0] : nullptr);
        job_free(job);
    });

    job->bin_task = bin_task;
    task_system::set_task_dependency(finalize_task, bin_task);
    task_system::set_task_dependency(done_task, finalize_task);
    task_system::enqueue_task(bin_task);
    return job;
}

void density_map_cancel(DensityMapJob* job, bool wait) {
    if (!job) return;
    job->cancel.store(1);
    if (wait) task_system::task_interrupt_and_wait_for(job->bin_task);
}

// ## GL

uint32_t density_map_texture_create(bool periodic) {
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    const GLint wrap = periodic ? GL_REPEAT : GL_CLAMP_TO_BORDER;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap);
    const float border[4] = { 0, 0, 0, 0 };
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, DENSITY_MAP_DIM, DENSITY_MAP_DIM);
    glBindTexture(GL_TEXTURE_2D, 0);
    return tex;
}

void density_map_texture_upload(uint32_t tex, const float* grid) {
    gl::set_texture_2D_data(tex, 0, grid, GL_RGBA32F);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void density_map_texture_free(uint32_t* tex) {
    if (tex && *tex) {
        glDeleteTextures(1, tex);
        *tex = 0;
    }
}

static constexpr str_t v_fs_triangle_src = STR_LIT(R"(
#version 150 core

void main() {
    uint idx = uint(gl_VertexID) % 3U;
    gl_Position = vec4(
        (float( idx     &1U)) * 4.0 - 1.0,
        (float((idx>>1U)&1U)) * 4.0 - 1.0,
        0, 1.0);
}
)");

// A layer over the canvas, premultiplied: contours (fills of the regions above the levels, lines along them) or a heat
// map, of the channels weighted. Sampled with a cubic B-spline, from four bilinear fetches.
static constexpr str_t f_layer_src = STR_LIT(R"(
#version 330 core

out vec4 out_color;

uniform sampler2D u_tex;
uniform vec2  u_tex_size;
uniform vec4  u_uv_rect;        // The uv at the lower left corner, and the extent
uniform vec2  u_inv_res;
uniform vec4  u_weights;

uniform int   u_mode;           // 0: contours, 1: heat map
uniform int   u_fill;
uniform int   u_lines;
uniform int   u_num_levels;
uniform float u_levels[8];      // Densities, the highest first
uniform vec4  u_fill_colors[8];
uniform vec4  u_line_color;
uniform float u_line_width;     // Pixels

uniform vec4  u_cmap[16];
uniform int   u_cmap_count;
uniform float u_cmap_max;
uniform int   u_cmap_log;

uniform float u_opacity;

uniform int   u_clip;           // Within the triangle u_clip_tri (uv) only
uniform vec2  u_clip_tri[3];

float sample_density(vec2 uv) {
    vec2 st = uv * u_tex_size - 0.5;
    vec2 i  = floor(st);
    vec2 f  = st - i;
    vec2 f2 = f * f;
    vec2 f3 = f2 * f;
    vec2 w0 = (1.0 / 6.0) * (-f3 + 3.0 * f2 - 3.0 * f + 1.0);
    vec2 w1 = (1.0 / 6.0) * (3.0 * f3 - 6.0 * f2 + 4.0);
    vec2 w2 = (1.0 / 6.0) * (-3.0 * f3 + 3.0 * f2 + 3.0 * f + 1.0);
    vec2 w3 = (1.0 / 6.0) * f3;
    vec2 g0 = w0 + w1;
    vec2 g1 = w2 + w3;
    vec2 h0 = (i - 1.0 + w1 / g0 + 0.5) / u_tex_size;
    vec2 h1 = (i + 1.0 + w3 / g1 + 0.5) / u_tex_size;
    vec4 t00 = texture(u_tex, vec2(h0.x, h0.y));
    vec4 t10 = texture(u_tex, vec2(h1.x, h0.y));
    vec4 t01 = texture(u_tex, vec2(h0.x, h1.y));
    vec4 t11 = texture(u_tex, vec2(h1.x, h1.y));
    vec4 t = g0.y * (g0.x * t00 + g1.x * t10) + g1.y * (g0.x * t01 + g1.x * t11);
    return dot(t, u_weights);
}

vec4 over(vec4 src, vec4 dst) {
    return src + dst * (1.0 - src.a);
}

vec4 premultiply(vec4 c) {
    return vec4(c.rgb * c.a, c.a);
}

// How much of the pixel is within the triangle: its distance in pixels to the nearest edge, inside
float clip_coverage(vec2 uv) {
    vec2 a = u_clip_tri[0], b = u_clip_tri[1], c = u_clip_tri[2];
    float s = sign((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
    float cov = 1.0;
    for (int k = 0; k < 3; ++k) {
        vec2 p = u_clip_tri[k];
        vec2 q = u_clip_tri[(k + 1) % 3];
        float e = s * ((q.x - p.x) * (uv.y - p.y) - (q.y - p.y) * (uv.x - p.x));
        float g = max(length(vec2(dFdx(e), dFdy(e))), 1.0e-30);
        cov = min(cov, clamp(e / g + 0.5, 0.0, 1.0));
    }
    return cov;
}

void main() {
    vec2 uv = u_uv_rect.xy + gl_FragCoord.xy * u_inv_res * u_uv_rect.zw;
    float d = sample_density(uv);
    float g = max(length(vec2(dFdx(d), dFdy(d))), 1.0e-30);
    float clip = clip_coverage(uv);

    vec4 c = vec4(0.0);
    if (u_mode == 0) {
        if (u_fill != 0) {
            for (int i = u_num_levels - 1; i >= 0; --i) {
                float cov = clamp((d - u_levels[i]) / g + 0.5, 0.0, 1.0);
                c = over(premultiply(u_fill_colors[i]) * cov, c);
            }
        }
        if (u_lines != 0) {
            float a = 0.0;
            for (int i = 0; i < u_num_levels; ++i) {
                float dist = abs(d - u_levels[i]) / g;
                a = max(a, clamp(0.5 * u_line_width + 0.5 - dist, 0.0, 1.0));
            }
            c = over(premultiply(u_line_color) * a, c);
        }
    } else {
        float v = d / max(u_cmap_max, 1.0e-30);
        if (u_cmap_log != 0) v = 1.0 + log(max(v, 1.0e-4)) / log(10.0) * 0.25;
        v = clamp(v, 0.0, 1.0);
        float s = v * float(u_cmap_count - 1);
        int i0 = int(floor(s));
        int i1 = min(i0 + 1, u_cmap_count - 1);
        vec4 col = mix(u_cmap[i0], u_cmap[i1], s - float(i0));
        col.a *= smoothstep(0.0, 0.06, v);
        c = premultiply(col);
    }
    out_color = c * u_opacity * (u_clip != 0 ? clip : 1.0);
}
)");

static struct {
    GLuint program = 0;
    struct {
        GLint tex, tex_size, uv_rect, inv_res, weights, mode, fill, lines, num_levels, levels, fill_colors, line_color,
              line_width, cmap, cmap_count, cmap_max, cmap_log, opacity, clip, clip_tri;
    } loc = {};
    GLuint fbo = 0;
    GLuint vao = 0;
    GLuint marks = 0;
} gl_state;

static void gl_init() {
    if (gl_state.program) return;
    GLuint v_shader = gl::compile_shader_from_source(v_fs_triangle_src, GL_VERTEX_SHADER);
    GLuint f_shader = gl::compile_shader_from_source(f_layer_src, GL_FRAGMENT_SHADER);
    gl_state.program = glCreateProgram();
    const GLuint shaders[] = { v_shader, f_shader };
    gl::attach_link_detach(gl_state.program, shaders, (int)ARRAY_SIZE(shaders));
    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    const GLuint p = gl_state.program;
    auto& l = gl_state.loc;
    l.tex         = glGetUniformLocation(p, "u_tex");
    l.tex_size    = glGetUniformLocation(p, "u_tex_size");
    l.uv_rect     = glGetUniformLocation(p, "u_uv_rect");
    l.inv_res     = glGetUniformLocation(p, "u_inv_res");
    l.weights     = glGetUniformLocation(p, "u_weights");
    l.mode        = glGetUniformLocation(p, "u_mode");
    l.fill        = glGetUniformLocation(p, "u_fill");
    l.lines       = glGetUniformLocation(p, "u_lines");
    l.num_levels  = glGetUniformLocation(p, "u_num_levels");
    l.levels      = glGetUniformLocation(p, "u_levels");
    l.fill_colors = glGetUniformLocation(p, "u_fill_colors");
    l.line_color  = glGetUniformLocation(p, "u_line_color");
    l.line_width  = glGetUniformLocation(p, "u_line_width");
    l.cmap        = glGetUniformLocation(p, "u_cmap");
    l.cmap_count  = glGetUniformLocation(p, "u_cmap_count");
    l.cmap_max    = glGetUniformLocation(p, "u_cmap_max");
    l.cmap_log    = glGetUniformLocation(p, "u_cmap_log");
    l.opacity     = glGetUniformLocation(p, "u_opacity");
    l.clip        = glGetUniformLocation(p, "u_clip");
    l.clip_tri    = glGetUniformLocation(p, "u_clip_tri");

    glGenFramebuffers(1, &gl_state.fbo);
    glGenVertexArrays(1, &gl_state.vao);
}

void density_map_shutdown() {
    if (gl_state.program) glDeleteProgram(gl_state.program);
    if (gl_state.fbo)     glDeleteFramebuffers(1, &gl_state.fbo);
    if (gl_state.vao)     glDeleteVertexArrays(1, &gl_state.vao);
    if (gl_state.marks)   glDeleteTextures(1, &gl_state.marks);
    gl_state = {};
}

void density_map_draw_init(DensityMapDraw* d) {
    MEMSET(d, 0, sizeof(DensityMapDraw));
    d->opacity = 1.0f;
    d->line_width = 1.5f;
}

void density_map_draw_levels(DensityMapDraw* d, const DensityMapResult& r, int channel, const float* levels, int num_levels) {
    ASSERT(0 <= channel && channel <= DENSITY_MAP_ALL);
    MEMSET(d->weights, 0, sizeof(d->weights));
    if (channel == DENSITY_MAP_ALL) {
        for (int t = 0; t < DENSITY_MAP_CHANNELS; ++t) d->weights[t] = r.weights_all[t];
    } else {
        d->weights[channel] = 1.0f;
    }
    d->num_levels = CLAMP(num_levels, 0, DENSITY_MAP_MAX_LEVELS);
    // Levels increase in what they enclose, so their densities decrease: the densest first
    for (int i = 0; i < d->num_levels; ++i) d->levels[i] = density_map_threshold(r.table[channel], levels[i]);
    d->cmap_max = r.max[channel];
}

void density_map_draw_fills_from(DensityMapDraw* d, ImVec4 color) {
    const int n = MAX(1, d->num_levels);
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
    const int n = MIN(size, 16);
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

uint32_t density_map_composite(DensityMapCanvas* canvas, int width, int height, const float view_uv[4], const DensityMapDraw* draws, int num_draws, const float* clip_uv) {
    ASSERT(canvas);
    if (num_draws <= 0) return 0;
    gl_init();

    width  = CLAMP(width,  1, 4096);
    height = CLAMP(height, 1, 4096);
    if (!canvas->tex || canvas->width != width || canvas->height != height) {
        if (canvas->tex) glDeleteTextures(1, &canvas->tex);
        glGenTextures(1, &canvas->tex);
        glBindTexture(GL_TEXTURE_2D, canvas->tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, width, height);
        glBindTexture(GL_TEXTURE_2D, 0);
        canvas->width = width;
        canvas->height = height;
        canvas->hash = 0;
    }

    uint64_t h = md_hash64(draws, sizeof(DensityMapDraw) * num_draws, (uint64_t)canvas->tex);
    h = md_hash64(view_uv, sizeof(float) * 4, h);
    if (clip_uv) h = md_hash64(clip_uv, sizeof(float) * 6, h ^ 0xC11B);
    if (h == canvas->hash) return canvas->tex;
    canvas->hash = h;

    PUSH_GPU_SECTION("Density map");
    GLint last_fbo = 0;
    GLint last_viewport[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl_state.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, canvas->tex, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, width, height);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_SCISSOR_TEST);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFuncSeparate(GL_ONE, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    const auto& l = gl_state.loc;
    glUseProgram(gl_state.program);
    glBindVertexArray(gl_state.vao);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(l.tex, 0);
    glUniform2f(l.tex_size, (float)DENSITY_MAP_DIM, (float)DENSITY_MAP_DIM);
    glUniform4fv(l.uv_rect, 1, view_uv);
    glUniform2f(l.inv_res, 1.0f / width, 1.0f / height);
    glUniform1i(l.clip, clip_uv ? 1 : 0);
    const float no_clip[6] = { 0, 0, 1, 0, 0, 1 };
    glUniform2fv(l.clip_tri, 3, clip_uv ? clip_uv : no_clip);

    for (int i = 0; i < num_draws; ++i) {
        const DensityMapDraw& d = draws[i];
        glBindTexture(GL_TEXTURE_2D, d.tex);
        glUniform4fv(l.weights, 1, d.weights);
        glUniform1i(l.mode, d.mode);
        glUniform1i(l.fill, d.fill);
        glUniform1i(l.lines, d.lines);
        glUniform1i(l.num_levels, d.num_levels);
        glUniform1fv(l.levels, DENSITY_MAP_MAX_LEVELS, d.levels);
        glUniform4fv(l.fill_colors, DENSITY_MAP_MAX_LEVELS, d.fill_colors[0]);
        glUniform4fv(l.line_color, 1, d.line_color);
        glUniform1f(l.line_width, d.line_width);
        glUniform4fv(l.cmap, 16, d.cmap[0]);
        glUniform1i(l.cmap_count, MAX(d.cmap_count, 1));
        glUniform1f(l.cmap_max, d.cmap_max);
        glUniform1i(l.cmap_log, d.cmap_log);
        glUniform1f(l.opacity, d.opacity);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_BLEND);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, (GLuint)last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    POP_GPU_SECTION();

    return canvas->tex;
}

void density_map_canvas_free(DensityMapCanvas* canvas) {
    if (canvas && canvas->tex) glDeleteTextures(1, &canvas->tex);
    if (canvas) *canvas = {};
}

void density_map_canvas_invalidate(DensityMapCanvas* canvas) {
    if (canvas) canvas->hash = 0;
}

static void imgui_blend_premultiplied(const ImDrawList*, const ImDrawCmd*) {
    glBlendFuncSeparate(GL_ONE, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
}

void density_map_plot_canvas(uint32_t tex) {
    if (!tex) return;
    const ImVec2 pos = ImPlot::GetPlotPos();
    const ImVec2 size = ImPlot::GetPlotSize();
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlot::PushPlotClipRect();
    dl->AddCallback(imgui_blend_premultiplied, nullptr);
    dl->AddImage((ImTextureID)(intptr_t)tex, pos, ImVec2(pos.x + size.x, pos.y + size.y), ImVec2(0, 1), ImVec2(1, 0));
    dl->AddCallback(ImDrawCallback_ResetRenderState, nullptr);
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
static GLuint create_marks() {
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
                const float rim = CLAMP((r - 0.72f) / aa + 0.5f, 0.0f, 1.0f);
                lum = 1.0f - rim;
                alpha = CLAMP((1.0f - r) / aa + 0.5f, 0.0f, 1.0f) * (1.0f - 0.25f * rim);
            } else {
                alpha = CLAMP((1.0f - r) / aa + 0.5f, 0.0f, 1.0f) * CLAMP((r - 0.7f) / aa + 0.5f, 0.0f, 1.0f);
            }
            const uint32_t l8 = (uint32_t)(lum * 255.0f + 0.5f);
            const uint32_t a8 = (uint32_t)(alpha * 255.0f + 0.5f);
            pixels[y * 2 * S + x] = l8 | (l8 << 8) | (l8 << 16) | (a8 << 24);
        }
    }
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 2 * S, S, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
    return tex;
}

void density_map_marks_begin(ImDrawList* dl) {
    if (!gl_state.marks) gl_state.marks = create_marks();
    dl->PushTexture(ImTextureRef((ImTextureID)(intptr_t)gl_state.marks));
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
    const float w = fill ? MAX(1.0f, (ImGui::GetContentRegionAvail().x - (count - 1)) / count) : 0.0f;
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
        levels[n++] = MIN(0.999f, last + (1.0f - last) * 0.5f);
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
