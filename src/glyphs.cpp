#include <glyphs.h>
#include <attribute_access.h>
#include <qm_placement.h>
#include <display_units.h>

#include <md_system.h>
#include <md_system_state.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>

static const double REF_PERCENTILE = 0.99;  // the magnitude a quantity that is not a length or an area maps to its scale at

struct GlyphData {
    md_attribute_id_t points = MD_ATTRIBUTE_INVALID;
    uint64_t source_hash = 0;
    uint32_t num = 0;
    std::vector<double> pos, nrm;           // 3 per point as published: Angstrom, unit (empty: not oriented)

    // Where they are: per point a system atom (system_index), or the QM atoms' frame (atom_index), or as published
    enum Place { AsIs, PerAtom, QmFrame } place = AsIs;
    std::vector<uint32_t> atom_of;          // PerAtom: per point; QmFrame: per QM atom, its system atom
    std::vector<double> qm_xyz;             // QmFrame: 3 per QM atom, Angstrom
    uint64_t placement_hash = 0;

    uint64_t size_hash = 0;
    ColorScaleSpan span = {};
    std::vector<float> color_value;         // per point, the colour's quantity (NAN-free: 0 where it has none)

    std::vector<float> size;                // 2 per point: the radius, the height (Angstrom, signed)
    bool sized = false;

    // What md_render draws (glyphs_draw_data)
    bool cylinders = false;
    std::vector<float>    draw_points;      // 4 per point of a glyph (x, y, z, radius): a sphere's one, a cylinder's two
    std::vector<float>    draw_radius;      // per glyph
    std::vector<uint32_t> draw_colors;      // per glyph, RGBA8
    bool     placed = false;
    uint64_t color_hash = 0;
    uint32_t version = 0, color_version = 0;
};

// ---------------------------------------------------------------------------
// The point sets and their quantities
// ---------------------------------------------------------------------------

static const md_attribute_t* get(const md_system_t& sys, md_attribute_id_t id) {
    return id != MD_ATTRIBUTE_INVALID ? md_attributes_get(&sys.attributes, id) : nullptr;
}

// A series of `components` values each, of exactly count of them (0: any)
static bool is_series(const md_attribute_t* a, uint32_t components, uint32_t count) {
    return attribute_has_shape(a, 1, components) && (count == 0 || a->format.shape[0] == count);
}

static bool is_length(md_unit_t u) { return md_unit_base_equal(u, md_unit_angstrom()); }
static bool is_area(md_unit_t u)   { return md_unit_base_equal(u, md_unit_pow(md_unit_angstrom(), 2)); }

static const md_attribute_t* coordinate(const md_system_t& sys, md_attribute_id_t points) {
    const md_attribute_t* c = get(sys, points);
    if (!is_series(c, 3, 0) || !is_length(c->unit) || c->format.shape[0] == 0) return nullptr;
    if (!md_str_eq(md_attribute_leaf(c), MD_STR_LIT("coordinate"))) return nullptr;
    return c;
}

static const md_attribute_t* sibling(const md_system_t& sys, const md_attribute_t* c, const char* leaf) {
    return md_attributes_sibling(&sys.attributes, c, md_str_from_cstr(leaf));
}

bool glyphs_is_point_set(const md_system_t& sys, md_attribute_id_t points) {
    return coordinate(sys, points) != nullptr;
}

bool glyphs_is_oriented(const md_system_t& sys, md_attribute_id_t points) {
    const md_attribute_t* c = coordinate(sys, points);
    return c && is_series(sibling(sys, c, "normal"), 3, c->format.shape[0]);
}

md_str_t glyphs_group(const md_system_t& sys, md_attribute_id_t points) {
    const md_attribute_t* c = coordinate(sys, points);
    return c ? md_attribute_group(c) : md_str_t{};
}

size_t glyphs_point_sets(md_attribute_id_t out[], size_t cap, const md_system_t& sys) {
    size_t n = 0;
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, MD_STR_LIT("")); md_attributes_next(&it);) {
        if (md_attribute_is_alias(it.attr) || !coordinate(sys, it.attr->id)) continue;
        if (n < cap) out[n] = it.attr->id;
        n++;
    }
    return n;
}

size_t glyphs_quantities(md_attribute_id_t out[], size_t cap, const md_system_t& sys, md_attribute_id_t points) {
    const md_attribute_t* c = coordinate(sys, points);
    if (!c) return 0;
    const uint32_t count = c->format.shape[0];
    const md_str_t group = md_attribute_group(c);
    size_t n = 0;
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, group); md_attributes_next(&it);) {
        const md_attribute_t* a = it.attr;
        if (md_attribute_is_alias(a) || !md_str_eq(md_attribute_group(a), group) || !is_series(a, 1, count)) continue;
        if (n < cap) out[n] = a->id;
        n++;
    }
    return n;
}

bool glyphs_unit_print(char* buf, size_t cap, double* factor, md_unit_t unit) {
    md_unit_base_t b = unit.base;
    b.dim.length = b.dim.angle = b.dim.count = 0;
    const bool geometric = b.raw_bits == 0;
    double f = 1.0;
    if (geometric) {
        f = display_units::factor_print(buf, cap, unit);
    } else {
        display_units::print(buf, cap, unit);
    }
    if (factor) *factor = f;
    return geometric;
}

bool glyphs_scale_is_factor(const md_system_t& sys, md_attribute_id_t quantity) {
    const md_attribute_t* a = get(sys, quantity);
    return a && (is_length(a->unit) || is_area(a->unit));
}

float glyphs_default_scale(const md_system_t& sys, md_attribute_id_t quantity, bool radius) {
    if (!get(sys, quantity)) return radius ? 0.15f : 0.1f;
    return glyphs_scale_is_factor(sys, quantity) ? 1.0f : (radius ? 0.5f : 0.3f);
}

// The values of a quantity at every point, in its own unit (or, for a size, in Angstrom or Angstrom^2); empty when
// it is none
static std::vector<double> values(const md_system_t& sys, md_attribute_id_t quantity, uint32_t count, md_unit_t unit) {
    const md_attribute_t* a = get(sys, quantity);
    return is_series(a, 1, count) ? attribute_values(sys, a, unit) : std::vector<double>{};
}

static ColorScaleSpan span_of(const std::vector<double>& v) {
    ColorScaleSpan s = {};
    std::vector<double> w;
    w.reserve(v.size());
    for (double x : v) if (isfinite(x)) w.push_back(x);
    if (w.empty()) return s;
    std::sort(w.begin(), w.end());
    s.valid = true;
    s.min = (float)w.front();
    s.max = (float)w.back();
    s.lo  = (float)w[(size_t)(0.01 * (w.size() - 1))];
    s.hi  = (float)w[(size_t)(0.99 * (w.size() - 1))];
    return s;
}

void glyphs_select_color(GlyphRepresentation* g, const md_system_t& sys, md_attribute_id_t quantity) {
    g->color = get(sys, quantity) ? quantity : MD_ATTRIBUTE_INVALID;
    const md_attribute_t* c = coordinate(sys, g->points);
    if (!c || g->color == MD_ATTRIBUTE_INVALID) return;
    const md_attribute_t* a = get(sys, g->color);
    const ColorScaleSpan span = span_of(values(sys, g->color, c->format.shape[0], a->unit));
    const bool signed_values = color_scale_span_signed(span);
    g->scale.symmetric  = signed_values;
    g->scale.auto_range = true;
    if (signed_values && g->scale.colormap == COLOR_SCALE_COLORMAP_VIRIDIS) g->scale.colormap = COLOR_SCALE_COLORMAP_RDBU;
    else if (!signed_values && g->scale.colormap == COLOR_SCALE_COLORMAP_RDBU) g->scale.colormap = COLOR_SCALE_COLORMAP_VIRIDIS;
    color_scale_update_range(&g->scale, span);
}

void glyphs_select(GlyphRepresentation* g, const md_system_t& sys, md_attribute_id_t points) {
    if (!coordinate(sys, points)) {
        // the first oriented one, else the first
        md_attribute_id_t sets[64];
        const size_t n = MD_MIN(glyphs_point_sets(sets, MD_COUNTOF(sets), sys), MD_COUNTOF(sets));
        points = MD_ATTRIBUTE_INVALID;
        for (size_t i = 0; i < n && points == MD_ATTRIBUTE_INVALID; ++i) if (glyphs_is_oriented(sys, sets[i])) points = sets[i];
        if (points == MD_ATTRIBUTE_INVALID && n > 0) points = sets[0];
    }
    g->points = points;
    g->shape = glyphs_is_oriented(sys, points) ? GlyphShape::Cylinder : GlyphShape::Sphere;
    g->radius = {MD_ATTRIBUTE_INVALID, 0.15f};
    g->height = {MD_ATTRIBUTE_INVALID, 0.1f};
    md_attribute_id_t q[256];
    const size_t nq = MD_MIN(glyphs_quantities(q, MD_COUNTOF(q), sys, points), MD_COUNTOF(q));
    for (size_t i = 0; i < nq; ++i) {
        if (is_area(get(sys, q[i])->unit)) { g->radius = {q[i], 1.0f}; break; }
    }
    // coloured by what the set is about, as its publisher says (MD_ATTRIBUTE_FLAG_PRINCIPAL), else uniform
    md_attribute_id_t color = MD_ATTRIBUTE_INVALID;
    if (const md_attribute_t* p = md_attributes_principal(&sys.attributes, glyphs_group(sys, points))) {
        for (size_t i = 0; i < nq; ++i) if (q[i] == p->id) color = p->id;
    }
    glyphs_select_color(g, sys, color);
}

uint32_t glyphs_point_atom(const md_system_t& sys, md_attribute_id_t points, uint32_t point) {
    const md_attribute_t* c = coordinate(sys, points);
    if (!c || point >= c->format.shape[0]) return UINT32_MAX;
    double v = NAN;
    if (const md_attribute_t* s = sibling(sys, c, "system_index")) {
        if (is_series(s, 1, c->format.shape[0])) md_attribute_extract_f64(&v, 1, &sys.attributes, s, md_attribute_slice_1(point), s->unit);
        return isfinite(v) && v >= 0 && v < (double)sys.atom.count ? (uint32_t)v : UINT32_MAX;
    }
    const md_attribute_t* a = sibling(sys, c, "atom_index");
    if (!is_series(a, 1, c->format.shape[0])) return UINT32_MAX;
    md_attribute_extract_f64(&v, 1, &sys.attributes, a, md_attribute_slice_1(point), a->unit);
    if (!isfinite(v) || v < 0) return UINT32_MAX;
    const uint32_t qm = (uint32_t)v;
    const md_attribute_t* map = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/system_index"));
    double s = qm;
    if (map && is_series(map, 1, 0) && qm < map->format.shape[0]) md_attribute_extract_f64(&s, 1, &sys.attributes, map, md_attribute_slice_1(qm), map->unit);
    return isfinite(s) && s >= 0 && s < (double)sys.atom.count ? (uint32_t)s : UINT32_MAX;
}

double glyphs_value(const md_system_t& sys, md_attribute_id_t quantity, uint32_t point) {
    const md_attribute_t* a = get(sys, quantity);
    if (!is_series(a, 1, 0) || point >= a->format.shape[0]) return NAN;
    double v = NAN;
    if (md_attribute_extract_f64(&v, 1, &sys.attributes, a, md_attribute_slice_1(point), a->unit) != 1) return NAN;
    return v;
}

// ---------------------------------------------------------------------------
// Reading and placing
// ---------------------------------------------------------------------------

// What the points are read from: the group's attributes and the QM atoms', by id and version
static uint64_t source_hash(const md_system_t& sys, md_attribute_id_t points) {
    uint64_t h = md_hash64(&points, sizeof(points), 0x9e3779b97f4a7c15ull);
    const md_attribute_t* c = coordinate(sys, points);
    if (!c) return h;
    const char* leaves[] = {"coordinate", "normal", "system_index", "atom_index"};
    for (const char* leaf : leaves) {
        const md_attribute_t* a = sibling(sys, c, leaf);
        const uint64_t v[2] = {a ? (uint64_t)a->id : 0, a ? md_attributes_version(&sys.attributes, a->id) : 0};
        h = md_hash64(v, sizeof(v), h);
    }
    const md_str_t qm[] = {MD_STR_LIT("qm/atom/coordinate"), MD_STR_LIT("qm/atom/system_index")};
    for (const md_str_t& p : qm) {
        const md_attribute_t* a = md_attributes_find(&sys.attributes, p);
        const uint64_t v[2] = {a ? (uint64_t)a->id : 0, a ? md_attributes_version(&sys.attributes, a->id) : 0};
        h = md_hash64(v, sizeof(v), h);
    }
    return h;
}

static bool read_points(GlyphData* d, const md_system_t& sys, md_attribute_id_t points) {
    const md_attribute_t* c = coordinate(sys, points);
    if (!c) return false;
    const uint32_t n = c->format.shape[0];
    d->points = points;
    d->num = n;
    d->pos = attribute_values(sys, c, md_unit_angstrom());
    if (d->pos.size() != 3 * (size_t)n) return false;
    if (const md_attribute_t* nr = sibling(sys, c, "normal"); is_series(nr, 3, n)) d->nrm = attribute_values(sys, nr, md_unit_none());
    d->place = GlyphData::AsIs;
    if (const md_attribute_t* s = sibling(sys, c, "system_index"); is_series(s, 1, n)) {
        std::vector<double> v = values(sys, s->id, n, s->unit);
        d->atom_of.resize(n);
        for (uint32_t j = 0; j < n; ++j) d->atom_of[j] = j < v.size() && v[j] >= 0 ? (uint32_t)v[j] : UINT32_MAX;
        d->place = GlyphData::PerAtom;
    } else if (const md_attribute_t* a = sibling(sys, c, "atom_index"); is_series(a, 1, n)) {
        const md_attribute_t* q = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));
        if (is_series(q, 3, 0)) {
            const uint32_t m = q->format.shape[0];
            d->qm_xyz = attribute_values(sys, q, md_unit_angstrom());
            if (d->qm_xyz.size() == 3 * (size_t)m) {
                const md_attribute_t* map = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/system_index"));
                std::vector<double> v = is_series(map, 1, m) ? values(sys, map->id, m, map->unit) : std::vector<double>{};
                d->atom_of.resize(m);
                for (uint32_t i = 0; i < m; ++i) d->atom_of[i] = v.empty() ? i : (v[i] >= 0 ? (uint32_t)v[i] : UINT32_MAX);
                d->place = GlyphData::QmFrame;
            }
        }
    }
    return true;
}

// Where the glyphs are now (world, Angstrom): a sphere on each point, or a cylinder from it along its normal by its
// height (inward for a negative one, never shorter than a thousandth of an Angstrom)
static void place(GlyphData* d, const md_system_state_t& state) {
    const uint32_t n = d->num;
    QmPlacement pl;
    if (d->place == GlyphData::QmFrame) {
        qm_place(&pl, d->qm_xyz.data(), (uint32_t)d->atom_of.size(), state.xyz, state.num_atoms, d->atom_of.data(), 1.0);
    }
    const uint32_t per = d->cylinders ? 2 : 1;
    d->draw_points.assign(4 * (size_t)per * n, 0.0f);
    d->draw_radius.assign(n, 0.0f);
    for (uint32_t j = 0; j < n; ++j) {
        double p[3] = {d->pos[3 * j], d->pos[3 * j + 1], d->pos[3 * j + 2]};
        double q[3] = {0, 0, 0};
        if (!d->nrm.empty()) for (int k = 0; k < 3; ++k) q[k] = d->nrm[3 * j + k];
        if (d->place == GlyphData::QmFrame) {
            double pp[3], qq[3];
            qm_place_point(pl, p, pp);
            qm_place_dir(pl, q, qq);
            for (int k = 0; k < 3; ++k) { p[k] = pp[k]; q[k] = qq[k]; }
        } else if (d->place == GlyphData::PerAtom) {
            const uint32_t i = d->atom_of[j];
            if (state.xyz && i < state.num_atoms) { p[0] = state.xyz[i].x; p[1] = state.xyz[i].y; p[2] = state.xyz[i].z; }
        }
        const float r = d->size[2 * j];
        float* o = &d->draw_points[4 * (size_t)per * j];
        for (int k = 0; k < 3; ++k) o[k] = (float)p[k];
        o[3] = r;
        d->draw_radius[j] = r;
        if (per == 2) {
            const double len = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
            const double u[3] = {len > 0 ? q[0] / len : 0.0, len > 0 ? q[1] / len : 0.0, len > 0 ? q[2] / len : 1.0};
            const double h = d->size[2 * j + 1];
            const double L = h < 0.0 ? MD_MIN(h, -1e-3) : MD_MAX(h, 1e-3);
            for (int k = 0; k < 3; ++k) o[4 + k] = (float)(p[k] + u[k] * L);
            o[7] = r;
        }
    }
}

// The magnitude of the q-quantile of the values
static double quantile_magnitude(const std::vector<double>& v, double q) {
    std::vector<double> w;
    w.reserve(v.size());
    for (double x : v) if (isfinite(x)) w.push_back(fabs(x));
    if (w.empty()) return 0.0;
    const size_t k = MD_MIN(w.size() - 1, (size_t)(q * (double)w.size()));
    std::nth_element(w.begin(), w.begin() + k, w.end());
    return w[k];
}

// A size (Angstrom) of each point from a mapping, by the quantity's unit (glyphs.h); a height keeps its sign
static void sizes(std::vector<float>* out, const md_system_t& sys, uint32_t n, const GlyphMapping& m, bool radius) {
    out->assign(n, m.scale);
    const md_attribute_t* a = get(sys, m.quantity);
    if (!is_series(a, 1, n)) return;
    if (is_length(a->unit) || is_area(a->unit)) {
        const bool area = is_area(a->unit);
        const std::vector<double> v = values(sys, m.quantity, n, area ? md_unit_pow(md_unit_angstrom(), 2) : md_unit_angstrom());
        for (uint32_t j = 0; j < n && j < v.size(); ++j) {
            const double x = isfinite(v[j]) ? v[j] : 0.0;
            const double s = area ? copysign(sqrt(fabs(x) / MD_PI), x) : x;
            (*out)[j] = (float)(m.scale * (radius ? fabs(s) : s));
        }
        return;
    }
    const std::vector<double> v = values(sys, m.quantity, n, a->unit);
    const double ref = quantile_magnitude(v, REF_PERCENTILE);
    for (uint32_t j = 0; j < n && j < v.size(); ++j) {
        const double x = ref > 0 && isfinite(v[j]) ? v[j] / ref : 0.0;
        (*out)[j] = (float)(m.scale * (radius ? sqrt(fabs(x)) : x));
    }
}

static void set_sizes(GlyphData* d, const md_system_t& sys, const GlyphRepresentation& g) {
    const uint32_t n = d->num;
    std::vector<float> r, h;
    sizes(&r, sys, n, g.radius, true);
    sizes(&h, sys, n, g.height, false);
    d->size.assign(2 * (size_t)n, 0.0f);
    for (uint32_t j = 0; j < n; ++j) { d->size[2 * j + 0] = r[j]; d->size[2 * j + 1] = h[j]; }

    d->span = {};
    d->color_value.assign(n, 0.0f);
    if (const md_attribute_t* c = get(sys, g.color); is_series(c, 1, n)) {
        const std::vector<double> v = values(sys, g.color, n, c->unit);
        for (uint32_t j = 0; j < n && j < v.size(); ++j) d->color_value[j] = isfinite(v[j]) ? (float)v[j] : 0.0f;
        d->span = span_of(v);
    }
}

// ---------------------------------------------------------------------------
// Update
// ---------------------------------------------------------------------------

static uint64_t size_hash(const md_system_t& sys, const GlyphRepresentation& g) {
    uint64_t h = md_hash64(&g.radius, sizeof(g.radius), 0);
    h = md_hash64(&g.height, sizeof(g.height), h);
    h = md_hash64(&g.color, sizeof(g.color), h);
    const md_attribute_id_t ids[3] = {g.radius.quantity, g.height.quantity, g.color};
    for (md_attribute_id_t id : ids) {
        const uint64_t v = id != MD_ATTRIBUTE_INVALID ? md_attributes_version(&sys.attributes, id) : 0;
        h = md_hash64(&v, sizeof(v), h);
    }
    return h;
}

bool glyphs_update(GlyphRepresentation* g, const md_system_t& sys, const md_system_state_t& state) {
    if (!coordinate(sys, g->points)) return false;
    if (!g->data) g->data = new GlyphData();
    GlyphData* d = g->data;

    const uint64_t h = source_hash(sys, g->points);
    if (h != d->source_hash) {
        *d = GlyphData{};
        d->source_hash = h;
        if (!read_points(d, sys, g->points)) return false;
    }
    if (d->num == 0) return false;

    // Sizes and colours when their mapping (or what it maps) changed; then the points placed again, beside which the
    // colours are
    bool replace = false;
    const uint64_t sh = size_hash(sys, *g);
    if (sh != d->size_hash || !d->sized) {
        d->size_hash = sh;
        set_sizes(d, sys, *g);
        d->sized = true;
        replace = true;
    }
    // Cylinders where the points are oriented, else spheres
    const bool cylinders = g->shape == GlyphShape::Cylinder && !d->nrm.empty();
    if (cylinders != d->cylinders) {
        d->cylinders = cylinders;
        replace = true;
    }
    // Placed on the atoms as they are
    uint64_t ph = (uint64_t)d->place + 1;
    if (d->place != GlyphData::AsIs && state.xyz) {
        for (uint32_t i : d->atom_of) {
            if (i < state.num_atoms) ph = md_hash64(&state.xyz[i], sizeof(md_float3_t), ph);
        }
    }
    if (replace || ph != d->placement_hash || !d->placed) {
        d->placement_hash = ph;
        place(d, state);
        d->placed = true;
        d->version++;
    }
    color_scale_update_range(&g->scale, d->span);
    // The colours, by the scale over the values (or a uniform grey), opaque
    uint64_t ch = md_hash64(&g->scale, sizeof(g->scale), d->size_hash);
    ch = md_hash64(&g->color, sizeof(g->color), ch);
    if (ch != d->color_hash || d->draw_colors.size() != d->num) {
        d->color_hash = ch;
        d->draw_colors.resize(d->num);
        const bool uniform = g->color == MD_ATTRIBUTE_INVALID;
        for (uint32_t j = 0; j < d->num; ++j) {
            d->draw_colors[j] = uniform ? 0xFFBFBFBFu : (color_scale_color_u32(g->scale, d->color_value[j]) | 0xFF000000u);
        }
        d->color_version++;
    }
    return true;
}

void glyphs_invalidate(GlyphRepresentation* g) {
    if (g->data) g->data->source_hash = 0;
}

void glyphs_free(GlyphRepresentation* g) {
    if (!g->data) return;
    delete g->data;
    g->data = nullptr;
}

ColorScaleSpan glyphs_span(const GlyphRepresentation& g) {
    return g.data ? g.data->span : ColorScaleSpan{};
}

uint64_t glyphs_hash(const GlyphRepresentation& g) {
    uint64_t h = md_hash64(&g.points, sizeof(g.points), 0);
    h = md_hash64(&g.radius, sizeof(g.radius), h);
    h = md_hash64(&g.height, sizeof(g.height), h);
    h = md_hash64(&g.color, sizeof(g.color), h);
    h = md_hash64(&g.shape, sizeof(g.shape), h);
    h = md_hash64(&g.scale, sizeof(g.scale), h);
    h = md_hash64(&g.unlit, sizeof(g.unlit), h);
    if (g.data) h = md_hash64(&g.data->version, sizeof(g.data->version), h);
    return h;
}

uint32_t glyphs_count(const GlyphRepresentation& g) { return g.data ? g.data->num : 0; }

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

bool glyphs_draw_data(GlyphDraw* out, const GlyphRepresentation& g) {
    *out = {};
    const GlyphData* d = g.data;
    if (!d || d->num == 0 || !d->placed || d->draw_colors.size() != d->num) return false;
    out->points = d->draw_points.data();
    out->radius = d->draw_radius.data();
    out->colors = d->draw_colors.data();
    out->count = d->num;
    out->cylinders = d->cylinders;
    out->unlit = g.unlit;
    out->version = d->version;
    out->color_version = d->color_version;
    return true;
}
