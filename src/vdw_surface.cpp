#include <vdw_surface.h>
#include <glyphs.h>
#include <attribute_access.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/view_param.h>

#include <md_system.h>
#include <md_system_state.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_vec_math.h>

#include <imgui.h>
#include <implot.h>

#include <shaders.inl>

#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>

static const int GRID_MAX_DIM = 128;        // cells a side of the points' grid at most
static const float KERNEL_CUT = 9.21f;      // zeta d^2 past which a point's Gaussian is taken as 0 (e^-9.21 = 1e-4)
static const double AREA_EPS = 1e-6;        // of the largest area: a point below it has none (buried)

enum { BUF_SPHERE, BUF_SPHERE_COLOR, BUF_SPHERE_ATOM, BUF_POINT, BUF_POINT_W, BUF_CELL, BUF_COUNT };

struct VdwSurfaceData {
    // The spheres: of the atoms shown, compacted
    std::vector<uint32_t> atom_colors;      // per atom of the system, as the representation computed them
    uint64_t colors_version = 0;
    uint64_t sphere_hash = 0;
    uint32_t num_spheres = 0;

    // The points colouring it, sorted by the cell of the grid they are in
    uint64_t point_hash = 0;
    uint32_t num_points = 0;
    ColorScaleSpan span = {};
    float grid_origin[3] = {};
    float grid_cell = 1.0f;
    int   grid_dim[3] = {};

    GLuint buf[BUF_COUNT] = {};
    GLuint tex[BUF_COUNT] = {};
    GLuint colormap_tex = 0;
    int    colormap = -1;
    uint32_t version = 0;
};

static const md_attribute_t* get(const md_system_t& sys, md_attribute_id_t id) {
    return id != MD_ATTRIBUTE_INVALID ? md_attributes_get(&sys.attributes, id) : nullptr;
}

// G/leaf of the group a is in
static const md_attribute_t* sibling(const md_system_t& sys, const md_attribute_t* a, const char* leaf) {
    if (!a) return nullptr;
    char path[256];
    const md_str_t g = md_attribute_group(a);
    snprintf(path, sizeof(path), "%.*s/%s", (int)g.len, g.ptr, leaf);
    return md_attributes_find(&sys.attributes, md_str_from_cstr(path));
}

// The point set of a quantity: G/coordinate of its group G
static md_attribute_id_t points_of(const md_system_t& sys, md_attribute_id_t quantity) {
    const md_attribute_t* c = sibling(sys, get(sys, quantity), "coordinate");
    return c && glyphs_is_point_set(sys, c->id) ? c->id : MD_ATTRIBUTE_INVALID;
}

static bool is_length_series(const md_attribute_t* a, size_t count) {
    return attribute_has_shape(a, 1, 1) && (count == 0 || a->format.shape[0] == count) && md_unit_base_equal(a->unit, md_unit_angstrom());
}

size_t vdw_surface_radii_sets(md_attribute_id_t out[], size_t cap, const md_system_t& sys) {
    const size_t num_qm = attribute_count(sys, MD_STR_LIT("qm/atom/coordinate"));
    size_t n = 0;
    if (num_qm == 0) return 0;
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, MD_STR_LIT("")); md_attributes_next(&it);) {
        const md_attribute_t* a = it.attr;
        if (md_attribute_is_alias(a) || !md_str_eq(md_attribute_leaf(a), MD_STR_LIT("sphere_radius")) || !is_length_series(a, num_qm)) continue;
        if (n < cap) out[n] = a->id;
        n++;
    }
    return n;
}

size_t vdw_surface_quantities(md_attribute_id_t out[], size_t cap, const md_system_t& sys) {
    md_attribute_id_t sets[64];
    const size_t ns = MD_MIN(glyphs_point_sets(sets, MD_COUNTOF(sets), sys), MD_COUNTOF(sets));
    size_t n = 0;
    for (size_t s = 0; s < ns; ++s) {
        md_attribute_id_t q[256];
        const size_t nq = MD_MIN(glyphs_quantities(q, MD_COUNTOF(q), sys, sets[s]), MD_COUNTOF(q));
        for (size_t i = 0; i < nq; ++i) {
            const md_attribute_t* a = get(sys, q[i]);
            const md_str_t leaf = md_attribute_leaf(a);
            if (md_str_eq(leaf, MD_STR_LIT("atom_index")) || md_str_eq(leaf, MD_STR_LIT("system_index"))) continue;
            if (a->flags & MD_ATTRIBUTE_FLAG_EXTENSIVE) continue;
            if (n < cap) out[n] = q[i];
            n++;
        }
    }
    return n;
}

void vdw_surface_select_quantity(VdwSurfaceRepresentation* v, const md_system_t& sys, md_attribute_id_t quantity) {
    v->quantity = points_of(sys, quantity) != MD_ATTRIBUTE_INVALID ? quantity : MD_ATTRIBUTE_INVALID;
    if (v->quantity != MD_ATTRIBUTE_INVALID) glyphs_scale_for(&v->color_scale, sys, v->quantity);
}

bool vdw_surface_select_cavity(VdwSurfaceRepresentation* v, const md_system_t& sys, md_attribute_id_t radii) {
    const md_attribute_t* r = get(sys, radii);
    if (!r || !is_length_series(r, attribute_count(sys, MD_STR_LIT("qm/atom/coordinate")))) return false;
    v->radii = radii;
    // coloured by the principal quantity of a point set beside it (G/sphere_radius beside G/tessera/...)
    md_attribute_id_t q[256];
    const size_t nq = MD_MIN(vdw_surface_quantities(q, MD_COUNTOF(q), sys), MD_COUNTOF(q));
    const md_str_t group = md_attribute_group(r);
    md_attribute_id_t pick = MD_ATTRIBUTE_INVALID;
    for (size_t i = 0; i < nq && pick == MD_ATTRIBUTE_INVALID; ++i) {
        const md_attribute_t* a = get(sys, q[i]);
        if ((a->flags & MD_ATTRIBUTE_FLAG_PRINCIPAL) && md_str_begins_with(a->path, group)) pick = q[i];
    }
    vdw_surface_select_quantity(v, sys, pick);
    v->color_scale.show_legend = pick != MD_ATTRIBUTE_INVALID;   // what the colours mean, beside them
    return true;
}

void vdw_surface_set_colors(VdwSurfaceRepresentation* v, const uint32_t* colors, size_t count) {
    if (!v->data) v->data = new VdwSurfaceData();
    v->data->atom_colors.assign(colors, colors + count);
    v->data->colors_version++;
}

// ---------------------------------------------------------------------------
// GPU buffers
// ---------------------------------------------------------------------------

static void upload_buffer(VdwSurfaceData* d, int i, GLenum format, const void* data, size_t bytes) {
    if (!d->buf[i]) glGenBuffers(1, &d->buf[i]);
    if (!d->tex[i]) glGenTextures(1, &d->tex[i]);
    glBindBuffer(GL_TEXTURE_BUFFER, d->buf[i]);
    // A buffer texture of nothing is not allowed to be bound: one element at least
    static const uint32_t zero[4] = {};
    glBufferData(GL_TEXTURE_BUFFER, bytes ? bytes : sizeof(zero), bytes ? data : zero, GL_STATIC_DRAW);
    glBindTexture(GL_TEXTURE_BUFFER, d->tex[i]);
    glTexBuffer(GL_TEXTURE_BUFFER, format, d->buf[i]);
    glBindTexture(GL_TEXTURE_BUFFER, 0);
    glBindBuffer(GL_TEXTURE_BUFFER, 0);
}

static void free_gpu(VdwSurfaceData* d) {
    for (int i = 0; i < BUF_COUNT; ++i) {
        if (d->tex[i]) glDeleteTextures(1, &d->tex[i]);
        if (d->buf[i]) glDeleteBuffers(1, &d->buf[i]);
        d->tex[i] = d->buf[i] = 0;
    }
    if (d->colormap_tex) glDeleteTextures(1, &d->colormap_tex);
    d->colormap_tex = 0;
    d->colormap = -1;
}

static GLuint colormap_texture(VdwSurfaceData* d, int colormap) {
    if (d->colormap_tex && d->colormap == colormap) return d->colormap_tex;
    uint32_t pixels[256];
    for (int i = 0; i < 256; ++i) {
        ImVec4 col = ImPlot::SampleColormap((float)i / 255.0f, colormap);
        col.w = 1.0f;
        pixels[i] = ImGui::ColorConvertFloat4ToU32(col);
    }
    if (!d->colormap_tex) {
        gl::init_texture_2D(&d->colormap_tex, 256, 1, GL_RGBA8);
        glBindTexture(GL_TEXTURE_2D, d->colormap_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    gl::set_texture_2D_data(d->colormap_tex, 0, pixels, GL_RGBA8);
    d->colormap = colormap;
    return d->colormap_tex;
}

// ---------------------------------------------------------------------------
// The spheres
// ---------------------------------------------------------------------------

// Each atom's radius: of the QM atoms' series (by qm/atom/system_index), or the van der Waals radius times the scale;
// 0 for none
static std::vector<float> radii(const VdwSurfaceRepresentation& v, const md_system_t& sys) {
    const size_t n = sys.atom.count;
    std::vector<float> r(n, 0.0f);
    if (const md_attribute_t* a = get(sys, v.radii)) {
        const std::vector<double> qr = attribute_values(sys, a, md_unit_angstrom());
        const std::vector<double> map = attribute_values(sys, MD_STR_LIT("qm/atom/system_index"), md_unit_none(), 1, qr.size());
        for (size_t i = 0; i < qr.size(); ++i) {
            const double s = map.empty() ? (double)i : map[i];
            if (s >= 0 && s < (double)n && isfinite(qr[i])) r[(size_t)s] = (float)qr[i];
        }
        return r;
    }
    for (size_t i = 0; i < n; ++i) r[i] = md_atom_radius(&sys.atom, i) * v.scale;
    return r;
}

static void update_spheres(VdwSurfaceRepresentation* v, const md_system_t& sys, const md_system_state_t& state) {
    VdwSurfaceData* d = v->data;
    const size_t n = MD_MIN((size_t)sys.atom.count, d->atom_colors.size());
    uint64_t h = md_hash64(&d->colors_version, sizeof(d->colors_version), 0);
    h = md_hash64(&v->scale, sizeof(v->scale), h);
    h = md_hash64(&v->radii, sizeof(v->radii), h);
    const uint64_t rv = v->radii != MD_ATTRIBUTE_INVALID ? md_attributes_version(&sys.attributes, v->radii) : 0;
    h = md_hash64(&rv, sizeof(rv), h);
    if (state.xyz) {
        for (size_t i = 0; i < n && i < state.num_atoms; ++i) {
            if (d->atom_colors[i] >> 24) h = md_hash64(&state.xyz[i], sizeof(md_float3_t), h);
        }
    }
    if (h == d->sphere_hash && d->tex[BUF_SPHERE]) return;
    d->sphere_hash = h;

    const std::vector<float> r = radii(*v, sys);
    std::vector<float> s;
    std::vector<uint32_t> col, atom;
    s.reserve(4 * n);
    for (size_t i = 0; i < n && state.xyz && i < state.num_atoms; ++i) {
        if (!(d->atom_colors[i] >> 24) || !(r[i] > 0.0f)) continue;
        s.push_back(state.xyz[i].x);
        s.push_back(state.xyz[i].y);
        s.push_back(state.xyz[i].z);
        s.push_back(r[i]);
        col.push_back(d->atom_colors[i]);
        atom.push_back((uint32_t)i);
    }
    d->num_spheres = (uint32_t)atom.size();
    upload_buffer(d, BUF_SPHERE,       GL_RGBA32F, s.data(),    s.size() * sizeof(float));
    upload_buffer(d, BUF_SPHERE_COLOR, GL_RGBA8,   col.data(),  col.size() * sizeof(uint32_t));
    upload_buffer(d, BUF_SPHERE_ATOM,  GL_R32UI,   atom.data(), atom.size() * sizeof(uint32_t));
    d->version++;
}

// ---------------------------------------------------------------------------
// The points
// ---------------------------------------------------------------------------

struct PointGrid {
    float origin[3] = {};
    float cell = 1.0f;
    int   dim[3] = {1, 1, 1};
    std::vector<uint32_t> start;    // per cell, into order; one past the last cell
    std::vector<uint32_t> order;    // the points by cell
    int cell_of(const float* p) const {
        int c[3];
        for (int k = 0; k < 3; ++k) c[k] = MD_CLAMP((int)floorf((p[k] - origin[k]) / cell), 0, dim[k] - 1);
        return (c[2] * dim[1] + c[1]) * dim[0] + c[0];
    }
};

// Of the points (xyz, 3 each) in cells of at least `cell` a side, at most GRID_MAX_DIM of them a side
static PointGrid grid_build(const std::vector<float>& xyz, float cell) {
    PointGrid g;
    const size_t n = xyz.size() / 3;
    float lo[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (size_t j = 0; j < n; ++j) for (int k = 0; k < 3; ++k) { lo[k] = MD_MIN(lo[k], xyz[3 * j + k]); hi[k] = MD_MAX(hi[k], xyz[3 * j + k]); }
    if (n == 0) for (int k = 0; k < 3; ++k) lo[k] = hi[k] = 0.0f;
    float ext = 0.0f;
    for (int k = 0; k < 3; ++k) ext = MD_MAX(ext, hi[k] - lo[k]);
    g.cell = MD_MAX(cell, ext / (float)(GRID_MAX_DIM - 1));
    g.cell = MD_MAX(g.cell, 1e-3f);
    for (int k = 0; k < 3; ++k) {
        g.origin[k] = lo[k];
        g.dim[k] = MD_MAX(1, (int)floorf((hi[k] - lo[k]) / g.cell) + 1);
    }
    const size_t nc = (size_t)g.dim[0] * g.dim[1] * g.dim[2];
    std::vector<uint32_t> count(nc + 1, 0), of(n);
    for (size_t j = 0; j < n; ++j) { of[j] = (uint32_t)g.cell_of(&xyz[3 * j]); count[of[j]]++; }
    g.start.assign(nc + 1, 0);
    for (size_t c = 0; c < nc; ++c) g.start[c + 1] = g.start[c] + count[c];
    std::vector<uint32_t> fill(g.start.begin(), g.start.end() - 1);
    g.order.resize(n);
    for (size_t j = 0; j < n; ++j) g.order[fill[of[j]]++] = (uint32_t)j;
    return g;
}

// For each point the others within r (each pair once, j < k), by the grid
template <typename F>
static void grid_pairs(const PointGrid& g, const std::vector<float>& xyz, float r, F&& f) {
    const int reach = (int)ceilf(r / g.cell);
    const float r2 = r * r;
    const size_t n = xyz.size() / 3;
    for (size_t j = 0; j < n; ++j) {
        const float* p = &xyz[3 * j];
        int c[3];
        for (int k = 0; k < 3; ++k) c[k] = MD_CLAMP((int)floorf((p[k] - g.origin[k]) / g.cell), 0, g.dim[k] - 1);
        for (int z = MD_MAX(c[2] - reach, 0); z <= MD_MIN(c[2] + reach, g.dim[2] - 1); ++z)
        for (int y = MD_MAX(c[1] - reach, 0); y <= MD_MIN(c[1] + reach, g.dim[1] - 1); ++y)
        for (int x = MD_MAX(c[0] - reach, 0); x <= MD_MIN(c[0] + reach, g.dim[0] - 1); ++x) {
            const size_t cell = ((size_t)z * g.dim[1] + y) * g.dim[0] + x;
            for (uint32_t o = g.start[cell]; o < g.start[cell + 1]; ++o) {
                const uint32_t k = g.order[o];
                if (k <= j) continue;
                const float dx = xyz[3 * k] - p[0], dy = xyz[3 * k + 1] - p[1], dz = xyz[3 * k + 2] - p[2];
                const float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < r2) f((uint32_t)j, k, d2);
            }
        }
    }
}

// The values diffused over the surface for the diffusion length l (the heat kernel's width sqrt(2 t), D = 1), on the
// points of an area (a > 0): a point cloud Laplacian (Belkin, Sun and Wang's, of the kernel width eps = h^2 / 2 of the
// mean spacing h, out to 2.5 h) in its flux form, a_j dv_j/dt = sum_k W_jk (v_k - v_j), W symmetric: sum a_j v_j is
// kept. Explicit steps, as many as are stable.
static void diffuse(std::vector<double>* values, const std::vector<float>& xyz, const std::vector<double>& area, double h, double l) {
    const size_t n = values->size();
    if (n < 2 || !(l > 0.0) || !(h > 0.0)) return;
    const double eps = 0.5 * h * h;
    const float r = (float)(2.5 * h);
    const PointGrid g = grid_build(xyz, r);
    struct Edge { uint32_t j, k; double w; };
    std::vector<Edge> edges;
    const double norm = 1.0 / (4.0 * MD_PI * eps * eps);
    grid_pairs(g, xyz, r, [&](uint32_t j, uint32_t k, float d2) {
        if (!(area[j] > 0.0) || !(area[k] > 0.0)) return;
        edges.push_back({j, k, norm * area[j] * area[k] * exp(-(double)d2 / (4.0 * eps))});
    });
    if (edges.empty()) return;
    // stable while dt sum_k W_jk / a_j <= 1 for every j
    std::vector<double> deg(n, 0.0);
    for (const Edge& e : edges) { deg[e.j] += e.w; deg[e.k] += e.w; }
    double rate = 0.0;
    for (size_t j = 0; j < n; ++j) if (area[j] > 0.0) rate = MD_MAX(rate, deg[j] / area[j]);
    const double T = 0.5 * l * l;
    const int steps = MD_MIN((int)ceil(T * rate / 0.9), 10000);
    const double dt = T / steps;
    std::vector<double> flux(n);
    std::vector<double>& v = *values;
    double before = 0.0;
    for (size_t j = 0; j < n; ++j) if (area[j] > 0.0) before += area[j] * v[j];
    for (int s = 0; s < steps; ++s) {
        std::fill(flux.begin(), flux.end(), 0.0);
        for (const Edge& e : edges) {
            const double f = e.w * (v[e.k] - v[e.j]);
            flux[e.j] += f;
            flux[e.k] -= f;
        }
        for (size_t j = 0; j < n; ++j) if (area[j] > 0.0) v[j] += dt * flux[j] / area[j];
    }
    double after = 0.0;
    for (size_t j = 0; j < n; ++j) if (area[j] > 0.0) after += area[j] * v[j];
    MD_LOG_DEBUG("VdW Surface: the values diffused for %.2f Angstrom over %zu points (%zu pairs, %d steps); sum a v kept to %.1e",
                 l, n, edges.size(), steps, before != 0.0 ? fabs(after / before - 1.0) : fabs(after));
}

static void update_points(VdwSurfaceRepresentation* v, const md_system_t& sys, const md_system_state_t& state) {
    VdwSurfaceData* d = v->data;
    const md_attribute_t* q = get(sys, v->quantity);
    const md_attribute_id_t points = points_of(sys, v->quantity);
    uint64_t h = md_hash64(&v->quantity, sizeof(v->quantity), 0);
    if (q) {
        h = md_hash64(&v->smoothing, sizeof(v->smoothing), h);
        const uint64_t qv = md_attributes_version(&sys.attributes, v->quantity);
        h = md_hash64(&qv, sizeof(qv), h);
        const uint64_t ph = glyphs_placement_hash(sys, state, points);
        h = md_hash64(&ph, sizeof(ph), h);
    }
    if (h == d->point_hash && d->tex[BUF_POINT]) return;
    d->point_hash = h;
    d->num_points = 0;
    d->span = {};

    std::vector<float> xyz = q ? glyphs_place_points(sys, state, points) : std::vector<float>{};
    const size_t n = xyz.size() / 3;
    std::vector<double> value = n ? attribute_values(sys, q, q->unit) : std::vector<double>{};
    if (value.size() != n) { xyz.clear(); value.clear(); }
    const md_attribute_t* c = get(sys, points);
    std::vector<double> area = n ? attribute_values(sys, sibling(sys, c, "area"), md_unit_pow(md_unit_angstrom(), 2)) : std::vector<double>{};
    std::vector<double> zeta = n ? attribute_values(sys, sibling(sys, c, "exponent"), md_unit_inv(md_unit_pow(md_unit_angstrom(), 2))) : std::vector<double>{};
    if (area.size() != n) area.assign(n, 1.0);
    // a point of no area (buried), or of no value, takes no part
    double amax = 0.0;
    for (double a : area) if (isfinite(a)) amax = MD_MAX(amax, a);
    for (size_t j = 0; j < n; ++j) if (!(area[j] > AREA_EPS * amax) || !isfinite(value[j])) area[j] = 0.0;

    // the spacing: the square root of the mean area, else the mean distance to the nearest point
    double spacing = 0.0;
    if (n) {
        if (sibling(sys, c, "area")) {
            double sum = 0.0; size_t m = 0;
            for (double a : area) if (a > 0.0) { sum += a; m++; }
            spacing = m ? sqrt(sum / m) : 0.0;
        } else {
            std::vector<float> nn(n, FLT_MAX);
            const PointGrid g = grid_build(xyz, 1.0f);
            grid_pairs(g, xyz, 4.0f, [&](uint32_t j, uint32_t k, float d2) { nn[j] = MD_MIN(nn[j], d2); nn[k] = MD_MIN(nn[k], d2); });
            double sum = 0.0; size_t m = 0;
            for (float x : nn) if (x < FLT_MAX) { sum += sqrt(x); m++; }
            spacing = m ? sum / m : 1.0;
        }
    }
    if (zeta.size() != n) zeta.assign(n, spacing > 0.0 ? 1.0 / (2.0 * spacing * spacing) : 1.0);

    diffuse(&value, xyz, area, spacing, v->smoothing);
    std::vector<double> shown;
    for (size_t j = 0; j < n; ++j) if (area[j] > 0.0) shown.push_back(value[j]);

    // the grid the pixels look up: cells as wide as the widest Gaussian reaches
    float reach = 0.0f;
    for (size_t j = 0; j < n; ++j) if (area[j] > 0.0 && zeta[j] > 0.0) reach = MD_MAX(reach, sqrtf(KERNEL_CUT / (float)zeta[j]));
    reach = MD_CLAMP(reach, 0.1f, 8.0f);
    const PointGrid g = grid_build(xyz, reach);
    std::vector<float> p(4 * n), w(2 * n);
    for (size_t o = 0; o < n; ++o) {
        const uint32_t j = g.order[o];
        for (int k = 0; k < 3; ++k) p[4 * o + k] = xyz[3 * j + k];
        p[4 * o + 3] = (float)value[j];
        w[2 * o + 0] = (float)area[j];
        w[2 * o + 1] = (float)MD_MAX(zeta[j], (double)(KERNEL_CUT / (reach * reach)));   // never past the grid's reach
    }
    std::vector<uint32_t> cells(2 * (g.start.size() - 1));
    for (size_t cc = 0; cc + 1 < g.start.size(); ++cc) { cells[2 * cc] = g.start[cc]; cells[2 * cc + 1] = g.start[cc + 1] - g.start[cc]; }
    upload_buffer(d, BUF_POINT,   GL_RGBA32F, p.data(), p.size() * sizeof(float));
    upload_buffer(d, BUF_POINT_W, GL_RG32F,   w.data(), w.size() * sizeof(float));
    upload_buffer(d, BUF_CELL,    GL_RG32UI,  cells.data(), cells.size() * sizeof(uint32_t));
    for (int k = 0; k < 3; ++k) { d->grid_origin[k] = g.origin[k]; d->grid_dim[k] = g.dim[k]; }
    d->grid_cell = g.cell;
    d->num_points = (uint32_t)n;

    // the span of what is shown: of the values at the points (diffused), the 1st to 99th percentile
    if (!shown.empty()) {
        std::sort(shown.begin(), shown.end());
        d->span.valid = true;
        d->span.min = (float)shown.front();
        d->span.max = (float)shown.back();
        d->span.lo  = (float)shown[(size_t)(0.01 * (shown.size() - 1))];
        d->span.hi  = (float)shown[(size_t)(0.99 * (shown.size() - 1))];
    }
    d->version++;
}

bool vdw_surface_update(VdwSurfaceRepresentation* v, const md_system_t& sys, const md_system_state_t& state) {
    if (!v->data) v->data = new VdwSurfaceData();
    if (v->data->atom_colors.empty()) return false;
    update_spheres(v, sys, state);
    update_points(v, sys, state);
    if (v->quantity != MD_ATTRIBUTE_INVALID) color_scale_update_range(&v->color_scale, v->data->span);
    return v->data->num_spheres > 0;
}

void vdw_surface_invalidate(VdwSurfaceRepresentation* v) {
    if (!v->data) return;
    v->data->sphere_hash = 0;
    v->data->point_hash = 0;
}

void vdw_surface_free(VdwSurfaceRepresentation* v) {
    if (!v->data) return;
    free_gpu(v->data);
    delete v->data;
    v->data = nullptr;
}

ColorScaleSpan vdw_surface_span(const VdwSurfaceRepresentation& v) {
    return v.data && v.quantity != MD_ATTRIBUTE_INVALID ? v.data->span : ColorScaleSpan{};
}

uint64_t vdw_surface_hash(const VdwSurfaceRepresentation& v) {
    uint64_t h = md_hash64(&v.scale, sizeof(v.scale), 0);
    h = md_hash64(&v.radii, sizeof(v.radii), h);
    h = md_hash64(&v.quantity, sizeof(v.quantity), h);
    h = md_hash64(&v.color_scale, sizeof(v.color_scale), h);
    h = md_hash64(&v.smoothing, sizeof(v.smoothing), h);
    if (v.data) h = md_hash64(&v.data->version, sizeof(v.data->version), h);
    return h;
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

static struct {
    GLuint program = 0;
    GLuint vao = 0;
    bool   failed = false;
} gl_state;

static bool ensure_program() {
    if (gl_state.failed) return false;
    if (!gl_state.program) {
        GLuint vs = gl::compile_shader_from_source({(const char*)vdw_surface_vert, vdw_surface_vert_size}, GL_VERTEX_SHADER);
        GLuint fs = gl::compile_shader_from_source({(const char*)vdw_surface_frag, vdw_surface_frag_size}, GL_FRAGMENT_SHADER);
        if (vs && fs) {
            gl_state.program = glCreateProgram();
            const GLuint shaders[] = {vs, fs};
            if (!gl::attach_link_detach(gl_state.program, shaders, 2)) {
                MD_LOG_ERROR("VdW Surface: the program did not link");
                glDeleteProgram(gl_state.program);
                gl_state.program = 0;
            }
        }
        if (vs) glDeleteShader(vs);
        if (fs) glDeleteShader(fs);
        if (!gl_state.program) {
            gl_state.failed = true;
            return false;
        }
    }
    if (!gl_state.vao) glGenVertexArrays(1, &gl_state.vao);
    return true;
}

static void jitter_uv(float out[2], const md_float4x4_t& P) {
    // As md_gl reads it off the projection
    if (P.elem[2][3] == 0.0f) {
        out[0] = -P.elem[3][0] * 0.5f;
        out[1] = -P.elem[3][1] * 0.5f;
    } else {
        out[0] = P.elem[2][0] * 0.5f;
        out[1] = P.elem[2][1] * 0.5f;
    }
}

void vdw_surface_draw(const VdwSurfaceRepresentation& v, const ViewParam& view, uint32_t atom_base) {
    VdwSurfaceData* d = v.data;
    if (!d || d->num_spheres == 0 || !d->tex[BUF_SPHERE] || !ensure_program()) return;
    const GLuint prog = gl_state.program;
    glUseProgram(prog);
    const char* names[BUF_COUNT] = {"u_sphere", "u_sphere_color", "u_sphere_atom", "u_point", "u_point_w", "u_cell"};
    for (int i = 0; i < BUF_COUNT; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, d->tex[i]);
        glUniform1i(glGetUniformLocation(prog, names[i]), i);
    }
    const bool by_points = v.quantity != MD_ATTRIBUTE_INVALID && d->num_points > 0 && d->tex[BUF_POINT];
    glActiveTexture(GL_TEXTURE0 + BUF_COUNT);
    glBindTexture(GL_TEXTURE_2D, colormap_texture(d, v.color_scale.colormap));
    glUniform1i(glGetUniformLocation(prog, "u_colormap"), BUF_COUNT);

    const md_float4x4_t& P = view.matrix.curr.proj;
    const md_float4x4_t curr_view_to_prev_clip = view.matrix.prev.proj * view.matrix.prev.view * view.matrix.inv.view;
    float jit[4];
    jitter_uv(jit + 0, view.matrix.curr.proj);
    jitter_uv(jit + 2, view.matrix.prev.proj);
    const bool ortho = P.elem[2][3] == 0.0f;
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_world_to_view"), 1, GL_FALSE, &view.matrix.curr.view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_to_world"), 1, GL_FALSE, &view.matrix.inv.view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_to_clip"), 1, GL_FALSE, &P.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip"), 1, GL_FALSE, &curr_view_to_prev_clip.elem[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, jit);
    glUniform1i(glGetUniformLocation(prog, "u_ortho"), ortho ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_by_points"), by_points ? 1 : 0);
    glUniform2f(glGetUniformLocation(prog, "u_range"), v.color_scale.range_beg, v.color_scale.range_end);
    glUniform3fv(glGetUniformLocation(prog, "u_grid_origin"), 1, d->grid_origin);
    glUniform1f(glGetUniformLocation(prog, "u_grid_cell"), d->grid_cell);
    glUniform3i(glGetUniformLocation(prog, "u_grid_dim"), d->grid_dim[0], d->grid_dim[1], d->grid_dim[2]);
    glUniform1ui(glGetUniformLocation(prog, "u_pick_base"), atom_base);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)d->num_spheres);
    glBindVertexArray(0);
    glEnable(GL_CULL_FACE);

    for (int i = 0; i < BUF_COUNT; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
    }
    glActiveTexture(GL_TEXTURE0 + BUF_COUNT);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
}
