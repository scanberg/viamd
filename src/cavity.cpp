#include <cavity.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/gbuffer.h>
#include <gfx/view_param.h>

#include <md_system.h>
#include <md_system_state.h>
#include <md_attributes.h>
#include <md_util_geom.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_vec_math.h>
#include <core/md_allocator.h>

#include <imgui.h>
#include <implot.h>

#include <shaders.inl>

#include <math.h>
#include <float.h>
#include <string.h>
#include <vector>
#include <algorithm>

// The attributes it is read from (md_vlx.h)
namespace cavity_path {
static const md_str_t coordinate   = MD_STR_INIT("vlx/scf/gostshyp/tessera/coordinate");
static const md_str_t normal       = MD_STR_INIT("vlx/scf/gostshyp/tessera/normal");
static const md_str_t area         = MD_STR_INIT("vlx/scf/gostshyp/tessera/area");
static const md_str_t amplitude    = MD_STR_INIT("vlx/scf/gostshyp/tessera/amplitude");
static const md_str_t exponent     = MD_STR_INIT("vlx/scf/gostshyp/tessera/exponent");
static const md_str_t decay_length = MD_STR_INIT("vlx/scf/gostshyp/tessera/decay_length");
static const md_str_t atom_index   = MD_STR_INIT("vlx/scf/gostshyp/tessera/atom_index");
static const md_str_t radius       = MD_STR_INIT("vlx/scf/gostshyp/sphere_radius");
static const md_str_t pressure     = MD_STR_INIT("vlx/scf/gostshyp/pressure");
static const md_str_t discretization = MD_STR_INIT("vlx/scf/gostshyp/discretization");
static const md_str_t num_lebedev  = MD_STR_INIT("vlx/scf/gostshyp/num_lebedev_points");
static const md_str_t r_ext        = MD_STR_INIT("vlx/scf/gostshyp/r_ext");
static const md_str_t qm_coord     = MD_STR_INIT("qm/atom/coordinate");
static const md_str_t qm_map       = MD_STR_INIT("qm/atom/system_index");
}

// Ångström and bohr, and the energy of P a lambda in kJ/mol for P in GPa and a lambda in Å^3: 1e9 Pa 1e-30 m^3 N_A
static const double BOHR_A        = 0.52917721092;
static const double GPA_A3_KJMOL  = 0.602214076;
static const double AU_PRESSURE_GPA = 29421.02648438959;

// The surface is integrated on points this far apart (Å), and no kernel is narrower than SIGMA_FLOOR, so that the
// quadrature resolves every one of them: the Gaussians of the tesserae at the seams of the spheres, whose areas the
// switching took nearly to nothing, are as narrow as their areas are small.
static const float QUAD_SPACING = 0.06f;
static const float SIGMA_FLOOR  = 0.12f;
static const float KERNEL_CUT   = 3.0f;     // in sigmas: the same cut in the normalization and in the shader
static const float SCREEN_PERIOD = 6.0f;    // pixels: of the halftone of the weight missing (cavity_resolve.frag)
static const float POWER_REACH  = 2.0f;     // Å: a tessera further than this from a point is never its patch by power, in
                                            // the fit as in the shader

struct CellGrid {
    float    min[3] = {};
    float    cell = 1.0f;
    int      dim[3] = {};
    std::vector<uint32_t> offset;   // [cells + 1]
    std::vector<uint32_t> items;

    int cell_of(const float p[3], int c[3]) const {
        for (int k = 0; k < 3; ++k) {
            c[k] = (int)floorf((p[k] - min[k]) / cell);
            if (c[k] < 0 || c[k] >= dim[k]) return -1;
        }
        return c[0] + dim[0] * (c[1] + dim[1] * c[2]);
    }
};

struct QuadPoint {
    float p[3];
    float n[3];
    float da;
    float w;        // GOSTSHYP's weight of the area here: its switching, 1 where the cavity is not switched
    float d_nn;     // to the nearest tessera facing the same way
    uint32_t nn;
    float p_nn;     // the least power, |p - c_j|^2 - w_j, of a tessera facing the same way
};

struct CavityData {
    uint64_t source_hash = 0;

    // The tesserae, in the frame of the file (Å)
    uint32_t num_tess = 0;
    std::vector<float>    pos;      // 3 per
    std::vector<float>    nrm;      // 3 per
    std::vector<double>   area, amp, zeta, lambda;
    std::vector<uint32_t> atom;
    bool   has_lambda = false;
    double pressure = 0.0;          // GPa, 0 when the file does not say

    // The spheres, one per QM atom, in the frame of the file: radius 0 where there is none
    uint32_t num_spheres = 0;
    std::vector<float>    sph;      // 4 per: centre, radius
    std::vector<uint32_t> sys_idx;  // the system's atom of each

    // The switching of VeloxChem's 'swig' (cavity_geometry.frag): available when the file says the cavity was made so,
    // and the spheres around each sphere within reach of its switching
    bool  has_switching = false;
    float sw_gamma = 0, sw_alpha = 0, sw_r_ext = 0;
    std::vector<uint32_t> nbr_off;  // [spheres + 1]
    std::vector<uint32_t> nbr;
    bool  switched = false;         // whether the surface below is weighted by the switching (QuadPoint::w)

    // The tesserae on the surface: those whose area is at least min_area of a full point (QuadPoint aside, a point of the
    // Lebedev grid that the switching left whole), measured on its own sphere. The rest are left out of everything
    // below: of the patches, the kernels and the grid the shader reads. They stay in the data.
    float min_area = -1.0f;
    double full_share = 0;          // the largest a_j / 4 pi R^2: a whole point of the largest orbit of the grid
    std::vector<uint8_t> shown;

    // The weights of the power diagram the surface is split by, |x - c_j|^2 - w_j (Å^2): fitted so that each patch
    // has the area of its tessera, relative to the surface drawn (fit_weights). Made once per cavity.
    std::vector<float> weight;

    // Made for a kernel width
    float kernel_scale = -1.0f;
    std::vector<float>  sigma;      // per tessera
    std::vector<double> norm;       // N_j
    std::vector<double> share;      // s_j: the scale of each kernel in the partition of unity whose parts have the areas
                                    // of the tesserae (fit_partition)
    double share_mismatch = 0;      // how far the parts are from those areas, relative to the surface they cover
    std::vector<QuadPoint> quad;
    double exposed_area = 0.0;      // of the union of the spheres
    double weighted_area = 0.0;     // the same weighted by the switching: what the tesserae's areas are measured in
    CellGrid lists;                 // per cell: every tessera that matters to a point of the surface in it

    // Made for a field, and whether the kernel's parts are fitted to the areas
    int field = -1;
    bool conserve = true;
    std::vector<float> nearest_value;   // per tessera
    std::vector<float> coef;            // 2 per: numerator, denominator
    bool conservative = false;          // the field is the numerator alone (Area coverage), else their ratio
    ColorScaleSpan span_nearest = {};
    ColorScaleSpan span_kernel = {};
    float kernel_error = 0;             // of a fitted field: how far its integral is from the sum of the values over the
                                        // areas of the tesserae

    // Where it is: world = R file + t
    float rot[9] = {1,0,0, 0,1,0, 0,0,1};   // column major
    float trans[3] = {0,0,0};
    uint64_t placement_hash = 0;

    // GPU
    GLuint buf[9] = {};             // tess_a, tess_b, tess_c, spheres, cells, items, nbr_range, nbr_items, glyph
    GLuint tex[9] = {};
    GLuint colormap_tex = 0;
    int    colormap = -1;
    uint64_t version = 0;
};

// ---------------------------------------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------------------------------------

static const md_attribute_t* find(const md_system_t& sys, md_str_t path) {
    return md_attributes_find(&sys.attributes, path);
}

bool cavity_available(const md_system_t& sys) {
    const md_attribute_t* c = find(sys, cavity_path::coordinate);
    return c && md_attribute_view(c, MD_ATTRIBUTE_TYPE_F64, 3, 1) && find(sys, cavity_path::atom_index) && find(sys, cavity_path::radius);
}

static uint64_t source_hash(const md_system_t& sys) {
    const md_str_t paths[] = {cavity_path::coordinate, cavity_path::normal, cavity_path::area, cavity_path::amplitude, cavity_path::exponent,
                              cavity_path::decay_length, cavity_path::atom_index, cavity_path::radius, cavity_path::pressure, cavity_path::qm_coord, cavity_path::qm_map,
                              cavity_path::discretization, cavity_path::num_lebedev, cavity_path::r_ext};
    uint64_t h = 0x9e3779b97f4a7c15ull;
    for (const md_str_t& p : paths) {
        const md_attribute_t* a = find(sys, p);
        const uint64_t v[2] = {a ? (uint64_t)a->id : 0, a ? a->version : 0};
        h = md_hash64(v, sizeof(v), h);
    }
    return h;
}

template <typename T>
static const T* view(const md_system_t& sys, md_str_t path, md_attribute_type_t type, uint32_t comp, uint32_t count) {
    const md_attribute_t* a = find(sys, path);
    if (!a || a->format.rank != 1 || a->format.shape[0] != count) return nullptr;
    return (const T*)md_attribute_view(a, type, comp, 1);
}

static bool read_source(CavityData* d, const md_system_t& sys) {
    const md_attribute_t* c = find(sys, cavity_path::coordinate);
    if (!c || c->format.rank != 1) return false;
    const uint32_t n = c->format.shape[0];
    const double* xyz  = view<double>(sys, cavity_path::coordinate, MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* nrm  = view<double>(sys, cavity_path::normal,     MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* area = view<double>(sys, cavity_path::area,       MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* amp  = view<double>(sys, cavity_path::amplitude,  MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* zeta = view<double>(sys, cavity_path::exponent,   MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* lam  = view<double>(sys, cavity_path::decay_length, MD_ATTRIBUTE_TYPE_F64, 1, n);
    const uint32_t* atom = view<uint32_t>(sys, cavity_path::atom_index, MD_ATTRIBUTE_TYPE_U32, 1, n);
    const md_attribute_t* r = find(sys, cavity_path::radius);
    const double* rad = r && r->format.rank == 1 ? (const double*)md_attribute_view(r, MD_ATTRIBUTE_TYPE_F64, 1, 1) : nullptr;
    if (!xyz || !nrm || !area || !amp || !zeta || !atom || !rad) {
        MD_LOG_ERROR("Cavity: the tesserae are incomplete");
        return false;
    }
    const uint32_t num_atoms = r->format.shape[0];

    d->num_tess = n;
    d->pos.resize(3 * n); d->nrm.resize(3 * n);
    d->area.assign(area, area + n); d->amp.assign(amp, amp + n); d->zeta.assign(zeta, zeta + n);
    d->has_lambda = lam != nullptr;
    d->lambda.assign(n, 0.0);
    if (lam) d->lambda.assign(lam, lam + n);
    d->atom.assign(atom, atom + n);
    for (uint32_t j = 0; j < n; ++j) {
        if (d->atom[j] >= num_atoms) {
            MD_LOG_ERROR("Cavity: tessera %u names atom %u of %u", j, d->atom[j], num_atoms);
            return false;
        }
        for (int k = 0; k < 3; ++k) {
            d->pos[3 * j + k] = (float)xyz[3 * j + k];
            d->nrm[3 * j + k] = (float)nrm[3 * j + k];
        }
    }

    d->pressure = 0.0;
    if (const md_attribute_t* p = find(sys, cavity_path::pressure)) {
        if (const double* v = (const double*)md_attribute_view(p, MD_ATTRIBUTE_TYPE_F64, 1, 0)) d->pressure = *v;
    }

    // The spheres: at the QM atoms where the file put them, and the system's atom each one is
    const double* qm = view<double>(sys, cavity_path::qm_coord, MD_ATTRIBUTE_TYPE_F64, 3, num_atoms);
    const uint32_t* map = view<uint32_t>(sys, cavity_path::qm_map, MD_ATTRIBUTE_TYPE_U32, 1, num_atoms);
    if (!qm) {
        MD_LOG_ERROR("Cavity: the QM atoms are not where the cavity is (qm/atom/coordinate)");
        return false;
    }
    d->num_spheres = num_atoms;
    d->sph.resize(4 * num_atoms);
    d->sys_idx.resize(num_atoms);
    for (uint32_t i = 0; i < num_atoms; ++i) {
        for (int k = 0; k < 3; ++k) d->sph[4 * i + k] = (float)qm[3 * i + k];
        d->sph[4 * i + 3] = isfinite(rad[i]) ? (float)rad[i] : 0.0f;
        d->sys_idx[i] = map ? map[i] : i;
        if (d->sys_idx[i] >= sys.atom.count) {
            MD_LOG_ERROR("Cavity: QM atom %u is atom %u of a system of %zu", i, d->sys_idx[i], sys.atom.count);
            return false;
        }
    }

    // The switching, as VeloxChem's tessellation.py makes it for 'swig'. 'iswig' switches with error functions of a
    // width per point, which is not done here: such a cavity is drawn as the union of its spheres.
    d->has_switching = false;
    const md_attribute_t* disc = find(sys, cavity_path::discretization);
    const md_str_t disc_str = disc && disc->format.type == MD_ATTRIBUTE_TYPE_STR ? md_attribute_str(&sys.attributes, disc, 0) : md_str_t{};
    double n_leb = 0, rext = 0;
    if (const md_attribute_t* a = find(sys, cavity_path::num_lebedev)) {
        if (const double* v = (const double*)md_attribute_view(a, MD_ATTRIBUTE_TYPE_F64, 1, 0)) n_leb = *v;
    }
    if (const md_attribute_t* a = find(sys, cavity_path::r_ext)) {
        if (const double* v = (const double*)md_attribute_view(a, MD_ATTRIBUTE_TYPE_F64, 1, 0)) rext = *v;
    }
    if (md_str_eq_cstr_ignore_case(disc_str, "swig") && n_leb > 0) {
        const double gamma = sqrt(14.0 / n_leb);
        d->has_switching = true;
        d->sw_gamma = (float)gamma;
        d->sw_alpha = (float)(0.5 + 1.0 / gamma - sqrt(1.0 / (gamma * gamma) - 1.0 / 28.0));
        d->sw_r_ext = (float)rext;
    }
    // The spheres within reach of each one's switching: of sphere j it ends (1 + gamma (1 - alpha)) r_j from its centre
    // (r_j grown by r_ext), and a point of sphere i is r_i + r_ext from its own
    d->nbr_off.assign(num_atoms + 1, 0);
    d->nbr.clear();
    for (uint32_t i = 0; i < num_atoms; ++i) {
        d->nbr_off[i] = (uint32_t)d->nbr.size();
        const float Ri = d->sph[4 * i + 3];
        if (Ri <= 0) continue;
        for (uint32_t j = 0; j < num_atoms; ++j) {
            const float Rj = d->sph[4 * j + 3];
            if (j == i || Rj <= 0) continue;
            const float rj = Rj + d->sw_r_ext;
            const float reach = (Ri + d->sw_r_ext) + rj * (1.0f + d->sw_gamma * (1.0f - d->sw_alpha));
            const float dx = d->sph[4 * i] - d->sph[4 * j], dy = d->sph[4 * i + 1] - d->sph[4 * j + 1], dz = d->sph[4 * i + 2] - d->sph[4 * j + 2];
            if (dx * dx + dy * dy + dz * dz < reach * reach) d->nbr.push_back(j);
        }
    }
    d->nbr_off[num_atoms] = (uint32_t)d->nbr.size();
    return true;
}

static inline float smoothstep5(float t) {
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    return t * t * t * (10.0f - 15.0f * t + 6.0f * t * t);
}

// The switching of the point of sphere i in the direction n (cavity_geometry.frag)
static float switching(const CavityData& d, uint32_t i, const float n[3]) {
    const float* s = &d.sph[4 * i];
    const float pe[3] = {s[0] + n[0] * (s[3] + d.sw_r_ext), s[1] + n[1] * (s[3] + d.sw_r_ext), s[2] + n[2] * (s[3] + d.sw_r_ext)};
    float F = 1.0f;
    for (uint32_t k = d.nbr_off[i]; k < d.nbr_off[i + 1] && F > 0.0f; ++k) {
        const float* sj = &d.sph[4 * d.nbr[k]];
        const float rj = sj[3] + d.sw_r_ext;
        const float dx = pe[0] - sj[0], dy = pe[1] - sj[1], dz = pe[2] - sj[2];
        F *= smoothstep5((sqrtf(dx * dx + dy * dy + dz * dz) - (rj - d.sw_alpha * d.sw_gamma * rj)) / (d.sw_gamma * rj));
    }
    return F;
}


// ---------------------------------------------------------------------------------------------------------
// The kernel's normalization and the lists of the cells
// ---------------------------------------------------------------------------------------------------------

static inline float dist2(const float* a, const float* b) {
    const float x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
    return x * x + y * y + z * z;
}

// A coarse grid of the tesserae on the CPU, for the queries made while building: cell >= the largest reach
static void build_coarse(CellGrid* g, const CavityData& d, float cell) {
    float lo[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (uint32_t i = 0; i < d.num_spheres; ++i) {
        const float R = d.sph[4 * i + 3];
        if (R <= 0) continue;
        for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], d.sph[4 * i + k] - R); hi[k] = fmaxf(hi[k], d.sph[4 * i + k] + R); }
    }
    for (uint32_t j = 0; j < d.num_tess; ++j) {
        for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], d.pos[3 * j + k]); hi[k] = fmaxf(hi[k], d.pos[3 * j + k]); }
    }
    g->cell = cell;
    for (int k = 0; k < 3; ++k) {
        g->min[k] = lo[k] - cell;
        g->dim[k] = (int)ceilf((hi[k] - lo[k]) / cell) + 3;
    }
    const size_t cells = (size_t)g->dim[0] * g->dim[1] * g->dim[2];
    std::vector<uint32_t> count(cells + 1, 0);
    std::vector<int> cell_of(d.num_tess);
    for (uint32_t j = 0; j < d.num_tess; ++j) {
        int c[3];
        cell_of[j] = d.shown[j] ? g->cell_of(&d.pos[3 * j], c) : -1;
        if (cell_of[j] >= 0) count[cell_of[j] + 1]++;
    }
    for (size_t i = 0; i < cells; ++i) count[i + 1] += count[i];
    g->offset = count;
    g->items.assign(count[cells], 0);
    std::vector<uint32_t> fill(count.begin(), count.end() - 1);
    for (uint32_t j = 0; j < d.num_tess; ++j) {
        if (cell_of[j] >= 0) g->items[fill[cell_of[j]]++] = j;
    }
}

// Calls f(j) for every tessera in the cells within reach of the box [lo, hi] grown by r
template <typename F>
static void query(const CellGrid& g, const float lo[3], const float hi[3], float r, F f) {
    int a[3], b[3];
    for (int k = 0; k < 3; ++k) {
        a[k] = MD_MAX(0, (int)floorf((lo[k] - r - g.min[k]) / g.cell));
        b[k] = MD_MIN(g.dim[k] - 1, (int)floorf((hi[k] + r - g.min[k]) / g.cell));
    }
    for (int z = a[2]; z <= b[2]; ++z)
        for (int y = a[1]; y <= b[1]; ++y)
            for (int x = a[0]; x <= b[0]; ++x) {
                const size_t c = (size_t)x + g.dim[0] * ((size_t)y + g.dim[1] * (size_t)z);
                for (uint32_t k = g.offset[c]; k < g.offset[c + 1]; ++k) f(g.items[k]);
            }
}

// The weights of the power diagram whose cells, on the exposed surface, have the areas of the tesserae: the dual of a
// semi-discrete optimal transport problem, concave in the weights, climbed by its gradient - the target area of a
// tessera less the area of its cell, both by the quadrature. Areas are measured as GOSTSHYP measures them, weighted by
// the switching (QuadPoint::w): a tessera's area is its Lebedev share of its sphere times the switching, so its cell is
// sized by the same measure, and a cell in a band along a seam, where the weight is low, is wide on the surface. The
// targets are the areas of the tesserae scaled to the weighted area of the union: the switching keeps part of the area
// inside other spheres, out of sight, so only the relative areas can be matched. A tessera no point of the surface can
// reach keeps a weight of 0 and no cell.
static void fit_weights(CavityData* d, const CellGrid& coarse) {
    const uint32_t n = d->num_tess;
    const float REACH = POWER_REACH;    // Å: the candidates of a point; the cells stay within it
    const float W_MAX = 3.0f;           // Å^2
    const int   ITERATIONS = 200;
    const float STEP = 0.25f;

    // The candidates of each point and their squared distances
    const size_t Q = d->quad.size();
    std::vector<uint32_t> off(Q + 1, 0), cand;
    std::vector<float> dist;
    cand.reserve(Q * 32);
    dist.reserve(Q * 32);
    std::vector<uint8_t> reached(n, 0);
    for (size_t i = 0; i < Q; ++i) {
        const QuadPoint& q = d->quad[i];
        query(coarse, q.p, q.p, REACH, [&](uint32_t j) {
            const float* nj = &d->nrm[3 * j];
            if (nj[0] * q.n[0] + nj[1] * q.n[1] + nj[2] * q.n[2] <= 0) return;
            const float dd = dist2(q.p, &d->pos[3 * j]);
            if (dd >= REACH * REACH) return;
            cand.push_back(j);
            dist.push_back(dd);
            reached[j] = 1;
        });
        off[i + 1] = (uint32_t)cand.size();
    }

    double sum_a = 0;
    for (uint32_t j = 0; j < n; ++j) if (reached[j]) sum_a += d->area[j];
    const double scale = sum_a > 0 ? d->weighted_area / sum_a : 0.0;
    std::vector<double> target(n), cell(n);
    for (uint32_t j = 0; j < n; ++j) target[j] = reached[j] ? d->area[j] * scale : 0.0;

    d->weight.assign(n, 0.0f);
    double mismatch = 0;
    for (int it = 0; it <= ITERATIONS; ++it) {
        std::fill(cell.begin(), cell.end(), 0.0);
        for (size_t i = 0; i < Q; ++i) {
            float best = FLT_MAX;
            uint32_t bj = UINT32_MAX;
            for (uint32_t k = off[i]; k < off[i + 1]; ++k) {
                const float v = dist[k] - d->weight[cand[k]];
                if (v < best) { best = v; bj = cand[k]; }
            }
            if (bj != UINT32_MAX) cell[bj] += d->quad[i].da * d->quad[i].w;
        }
        double err = 0;
        for (uint32_t j = 0; j < n; ++j) err += fabs(cell[j] - target[j]);
        mismatch = d->weighted_area > 0 ? err / d->weighted_area : 0.0;
        if (it == 0) MD_LOG_DEBUG("Cavity: the nearest tesserae split %.1f%% of the surface otherwise than their areas", 100.0 * mismatch);
        if (it == ITERATIONS) break;
        for (uint32_t j = 0; j < n; ++j) {
            if (!reached[j]) continue;
            d->weight[j] = MD_CLAMP(d->weight[j] + STEP * (float)(target[j] - cell[j]), -W_MAX, W_MAX);
        }
    }
    MD_LOG_DEBUG("Cavity: the power diagram splits %.1f%% of the surface otherwise than the areas of the tesserae", 100.0 * mismatch);
}

// The kernels as a partition of unity, phi_j(x) = s_j K_j(x) / sum_k s_k K_k(x), whose parts have the areas of the
// tesserae: the integral of phi_j F dA over the surface is a_j, scaled to the surface the kernels cover (the same
// relative areas as the power diagram, of which this is the smooth version, the kernel width its smoothing). The
// scales are fitted by Sinkhorn's iteration, s_j <- s_j a_j / (the integral of phi_j F dA). A field made of them,
// f = sum_j v_j phi_j, is everywhere a weighted mean of the values of the tesserae around it - it never leaves their
// range, and a field of one value is that value everywhere - and integrates over the surface to sum_j v_j a_j.
static void fit_partition(CavityData* d, const CellGrid& coarse) {
    const uint32_t n = d->num_tess;
    const int ITERATIONS = 300;
    const double TOLERANCE = 1e-3;

    // The kernels at each point
    const size_t Q = d->quad.size();
    std::vector<uint32_t> off(Q + 1, 0), cand;
    std::vector<float> kv;
    cand.reserve(Q * 32);
    kv.reserve(Q * 32);
    double covered = 0;
    for (size_t i = 0; i < Q; ++i) {
        const QuadPoint& q = d->quad[i];
        query(coarse, q.p, q.p, coarse.cell, [&](uint32_t j) {
            if (d->norm[j] <= 0) return;
            const float* nj = &d->nrm[3 * j];
            if (nj[0] * q.n[0] + nj[1] * q.n[1] + nj[2] * q.n[2] <= 0) return;
            const float dd = dist2(q.p, &d->pos[3 * j]);
            const float s = d->sigma[j];
            if (dd >= KERNEL_CUT * KERNEL_CUT * s * s) return;
            cand.push_back(j);
            kv.push_back(expf(-0.5f * dd / (s * s)));
        });
        off[i + 1] = (uint32_t)cand.size();
        if (off[i + 1] > off[i]) covered += (double)q.da * q.w;
    }

    double sum_a = 0;
    for (uint32_t j = 0; j < n; ++j) if (d->norm[j] > 0) sum_a += d->area[j];
    const double scale = sum_a > 0 ? covered / sum_a : 0.0;
    std::vector<double> target(n), part(n);
    d->share.assign(n, 0.0);
    for (uint32_t j = 0; j < n; ++j) {
        target[j] = d->norm[j] > 0 ? d->area[j] * scale : 0.0;
        d->share[j] = d->norm[j] > 0 ? target[j] / d->norm[j] : 0.0;   // the area-weighted kernels, to start from
    }

    int it = 0;
    for (;; ++it) {
        std::fill(part.begin(), part.end(), 0.0);
        for (size_t i = 0; i < Q; ++i) {
            double D = 0;
            for (uint32_t k = off[i]; k < off[i + 1]; ++k) D += d->share[cand[k]] * kv[k];
            if (D <= 0) continue;
            const double m = (double)d->quad[i].da * d->quad[i].w / D;
            for (uint32_t k = off[i]; k < off[i + 1]; ++k) part[cand[k]] += d->share[cand[k]] * kv[k] * m;
        }
        double err = 0;
        for (uint32_t j = 0; j < n; ++j) err += fabs(part[j] - target[j]);
        d->share_mismatch = covered > 0 ? err / covered : 0.0;
        if (it == 0) MD_LOG_DEBUG("Cavity: the area-weighted kernels give %.1f%% of the surface otherwise than the areas", 100.0 * d->share_mismatch);
        if (d->share_mismatch < TOLERANCE || it == ITERATIONS) break;
        for (uint32_t j = 0; j < n; ++j) {
            if (part[j] > 0 && target[j] > 0) d->share[j] *= target[j] / part[j];
        }
    }
    MD_LOG_DEBUG("Cavity: the kernels' partition of unity, %d iterations: %.2f%% of the surface otherwise than the areas", it, 100.0 * d->share_mismatch);
}

static void build_kernel(CavityData* d, float scale) {
    const uint32_t n = d->num_tess;
    d->kernel_scale = scale;
    d->sigma.resize(n);
    float sigma_max = 0;
    for (uint32_t j = 0; j < n; ++j) {
        const double s = d->zeta[j] > 0 ? 1.0 / sqrt(2.0 * d->zeta[j]) : SIGMA_FLOOR;
        d->sigma[j] = scale * fmaxf((float)s, SIGMA_FLOOR);
        sigma_max = fmaxf(sigma_max, d->sigma[j]);
    }

    // The nearest tessera of a point of the surface is never further than the spacing of the tesserae, which the
    // Lebedev grid of the largest sphere bounds: about 1.2 Å. The coarse grid reaches both that and the kernels.
    CellGrid coarse;
    build_coarse(&coarse, *d, fmaxf(KERNEL_CUT * sigma_max, POWER_REACH));

    // The exposed surface on points: a Fibonacci sphere per atom, the points inside another sphere left out, each
    // weighted by the switching where the cavity is switched
    d->quad.clear();
    d->exposed_area = 0;
    d->weighted_area = 0;
    for (uint32_t i = 0; i < d->num_spheres; ++i) {
        const float* c = &d->sph[4 * i];
        const float R = c[3];
        if (R <= 0) continue;
        std::vector<uint32_t> nb;
        for (uint32_t k = 0; k < d->num_spheres; ++k) {
            const float Rk = d->sph[4 * k + 3];
            if (k == i || Rk <= 0) continue;
            if (dist2(c, &d->sph[4 * k]) < (R + Rk) * (R + Rk)) nb.push_back(k);
        }
        const double area = 4.0 * M_PI * R * R;
        const uint32_t m = (uint32_t)ceil(area / (QUAD_SPACING * QUAD_SPACING));
        const double golden = M_PI * (3.0 - sqrt(5.0));
        for (uint32_t s = 0; s < m; ++s) {
            const double z = 1.0 - (2.0 * s + 1.0) / m;
            const double r = sqrt(fmax(0.0, 1.0 - z * z));
            const double phi = golden * s;
            QuadPoint q;
            q.n[0] = (float)(r * cos(phi)); q.n[1] = (float)(r * sin(phi)); q.n[2] = (float)z;
            for (int k = 0; k < 3; ++k) q.p[k] = c[k] + R * q.n[k];
            bool buried = false;
            for (uint32_t k : nb) {
                const float Rk = d->sph[4 * k + 3];
                if (dist2(q.p, &d->sph[4 * k]) < Rk * Rk) { buried = true; break; }
            }
            if (buried) continue;
            q.da = (float)(area / m);
            q.w = d->switched ? switching(*d, i, q.n) : 1.0f;
            q.d_nn = FLT_MAX;
            q.nn = UINT32_MAX;
            q.p_nn = FLT_MAX;
            d->quad.push_back(q);
            d->exposed_area += q.da;
            d->weighted_area += (double)q.da * q.w;
        }
    }

    if (d->weight.size() != n) fit_weights(d, coarse);

    // N_j, the integral of each kernel over the exposed surface by the weighted measure, and the nearest tessera of each
    // point, plain and by power
    d->norm.assign(n, 0.0);
    for (QuadPoint& q : d->quad) {
        query(coarse, q.p, q.p, coarse.cell, [&](uint32_t j) {
            const float* nj = &d->nrm[3 * j];
            if (nj[0] * q.n[0] + nj[1] * q.n[1] + nj[2] * q.n[2] <= 0) return;
            const float dd = dist2(q.p, &d->pos[3 * j]);
            if (dd < q.d_nn) { q.d_nn = dd; q.nn = j; }
            if (dd < POWER_REACH * POWER_REACH) q.p_nn = fminf(q.p_nn, dd - d->weight[j]);
            const float s = d->sigma[j];
            if (dd < KERNEL_CUT * KERNEL_CUT * s * s) d->norm[j] += exp(-0.5 * dd / (s * s)) * q.da * q.w;
        });
        q.d_nn = q.nn != UINT32_MAX ? sqrtf(q.d_nn) : FLT_MAX;
    }
    // A kernel that only reaches points of no weight (the very seam) spreads nothing either
    uint32_t num_unreached = 0;
    for (uint32_t j = 0; j < n; ++j) {
        if (d->norm[j] < 1e-9) d->norm[j] = 0.0;
        num_unreached += d->norm[j] <= 0.0 && d->shown[j];
    }
    if (num_unreached) {
        // Buried: the switching kept it, but no point of the exposed surface is within its reach. It colours nothing.
        MD_LOG_DEBUG("Cavity: %u of %u tesserae reach no exposed point of the surface", num_unreached, n);
    }
    fit_partition(d, coarse);

    // The lists of the cells. A cell of the surface gets every tessera within the furthest a point of it is from its
    // nearest tessera (the quadrature's, grown by its spacing), and every tessera whose kernel reaches into it. A cell
    // next to one of the surface gets the reach of its neighbours, for the points the quadrature missed in it.
    CellGrid& g = d->lists;
    float ext = 0;
    for (int k = 0; k < 3; ++k) ext = fmaxf(ext, (coarse.dim[k] - 3) * coarse.cell);
    const float cell = fmaxf(0.4f, ext / 96.0f);
    g.cell = cell;
    for (int k = 0; k < 3; ++k) {
        g.min[k] = coarse.min[k];
        g.dim[k] = (int)ceilf(coarse.dim[k] * coarse.cell / cell) + 1;
    }
    const size_t cells = (size_t)g.dim[0] * g.dim[1] * g.dim[2];
    // And by power: every tessera whose power can be below the least power of a point in it, grown by how much the
    // power of a tessera changes between the points of the quadrature (2 d h, d within the reach of the fit)
    std::vector<float> reach(cells, -1.0f), preach(cells, -FLT_MAX);
    for (const QuadPoint& q : d->quad) {
        int c[3];
        const int ci = g.cell_of(q.p, c);
        if (ci >= 0 && q.d_nn < FLT_MAX) reach[ci] = fmaxf(reach[ci], q.d_nn + QUAD_SPACING);
        if (ci >= 0 && q.p_nn < FLT_MAX) preach[ci] = fmaxf(preach[ci], q.p_nn + 4.0f * 2.0f * QUAD_SPACING);
    }
    float w_max = 0;
    for (float w : d->weight) w_max = fmaxf(w_max, w);
    std::vector<float> grown = reach, pgrown = preach;
    for (int z = 0; z < g.dim[2]; ++z)
        for (int y = 0; y < g.dim[1]; ++y)
            for (int x = 0; x < g.dim[0]; ++x) {
                const size_t ci = (size_t)x + g.dim[0] * ((size_t)y + g.dim[1] * (size_t)z);
                if (reach[ci] >= 0) continue;
                float r = -1.0f, pr = -FLT_MAX;
                for (int dz = -1; dz <= 1; ++dz)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            const int X = x + dx, Y = y + dy, Z = z + dz;
                            if (X < 0 || Y < 0 || Z < 0 || X >= g.dim[0] || Y >= g.dim[1] || Z >= g.dim[2]) continue;
                            const size_t ni = (size_t)X + g.dim[0] * ((size_t)Y + g.dim[1] * (size_t)Z);
                            r  = fmaxf(r, reach[ni]);
                            pr = fmaxf(pr, preach[ni]);
                        }
                grown[ci] = r;
                pgrown[ci] = pr;
            }

    g.offset.assign(cells + 1, 0);
    g.items.clear();
    for (int z = 0; z < g.dim[2]; ++z)
        for (int y = 0; y < g.dim[1]; ++y)
            for (int x = 0; x < g.dim[0]; ++x) {
                const size_t ci = (size_t)x + g.dim[0] * ((size_t)y + g.dim[1] * (size_t)z);
                g.offset[ci] = (uint32_t)g.items.size();
                const float r = grown[ci];
                if (r < 0) continue;
                const float pr = pgrown[ci];
                const float pr_reach = pr > -FLT_MAX ? sqrtf(fmaxf(pr + w_max, 0.0f)) : 0.0f;
                const float lo[3] = {g.min[0] + x * cell, g.min[1] + y * cell, g.min[2] + z * cell};
                const float hi[3] = {lo[0] + cell, lo[1] + cell, lo[2] + cell};
                query(coarse, lo, hi, fmaxf(fmaxf(r, pr_reach), KERNEL_CUT * sigma_max), [&](uint32_t j) {
                    // Its distance from the box
                    const float* p = &d->pos[3 * j];
                    float dd = 0;
                    for (int k = 0; k < 3; ++k) {
                        const float e = fmaxf(fmaxf(lo[k] - p[k], 0.0f), p[k] - hi[k]);
                        dd += e * e;
                    }
                    const float kr = KERNEL_CUT * d->sigma[j];
                    if (dd <= r * r || dd - d->weight[j] <= pr || (dd < kr * kr && d->norm[j] > 0)) g.items.push_back(j);
                });
            }
    g.offset[cells] = (uint32_t)g.items.size();
    MD_LOG_DEBUG("Cavity: %zu points on %.1f Å^2 of exposed surface (%.1f weighted); %d x %d x %d cells of %.2f Å, %zu entries",
                 d->quad.size(), d->exposed_area, d->weighted_area, g.dim[0], g.dim[1], g.dim[2], cell, g.items.size());
}

// ---------------------------------------------------------------------------------------------------------
// The field
// ---------------------------------------------------------------------------------------------------------

md_unit_t cavity_field_unit(CavityField field) {
    return field == CavityField::DecayLength ? md_unit_angstrom() : md_unit_none();
}

// Per tessera: the value of its patch, which the kernel averages too; and for Area coverage, the total it spreads
static void tessera_values(const CavityData& d, CavityField field, uint32_t j, double* nearest, double* v, bool* conservative) {
    const double a = d.area[j], lam = d.lambda[j], p = d.amp[j], P = d.pressure;
    *conservative = false;
    switch (field) {
    case CavityField::WallEnergy:
        // Per area, P lambda_j: the energy of the wall at it, P a_j lambda_j, over its area
        *nearest = *v = GPA_A3_KJMOL * P * lam;
        break;
    case CavityField::DecayLength:
        *nearest = *v = lam;
        break;
    case CavityField::WallDensity:
        // <rho|g_j> = P a_j lambda_j / p_j, in atomic units
        *nearest = *v = p > 0 ? (P / AU_PRESSURE_GPA) * (a / (BOHR_A * BOHR_A)) * (lam / BOHR_A) / p : 0;
        break;
    case CavityField::Amplitude:
        *nearest = *v = a > 0 ? p / a : 0;
        break;
    case CavityField::AreaCoverage:
    default:
        // The areas spread by the kernels each normalized to 1 over the surface, sum_j a_j K_j / N_j: 1 where the
        // tesserae account for the surface. What the kernels of the old conservative spreading made of a field of one
        // value, and why it was given up for the partition of unity
        *nearest = 1.0;
        *v = a; *conservative = true;
        break;
    }
}

static void build_field(CavityData* d, CavityField field, bool conserve) {
    const uint32_t n = d->num_tess;
    d->field = (int)field;
    d->conserve = conserve;
    d->nearest_value.resize(n);
    d->coef.resize(2 * n);
    for (uint32_t j = 0; j < n; ++j) {
        double nv, v;
        bool cons;
        tessera_values(*d, field, j, &nv, &v, &cons);
        d->conservative = cons;
        const double N = d->norm[j];
        d->nearest_value[j] = (float)nv;
        // The scale of the kernel in the partition of unity: fitted to the areas, or the area-weighted kernel a_j / N_j
        const double sj = (conserve && d->share.size() == n) ? d->share[j] : (N > 0 ? d->area[j] / N : 0.0);
        if (N > 0 && cons) {
            d->coef[2 * j + 0] = (float)(v / N);
            d->coef[2 * j + 1] = 0.0f;
        } else if (N > 0) {
            d->coef[2 * j + 0] = (float)(v * sj);
            d->coef[2 * j + 1] = (float)sj;
        } else {
            d->coef[2 * j + 0] = d->coef[2 * j + 1] = 0.0f;
        }
    }
}

// The field at a point of the surface (the frame of the file) with its normal, as the shader makes it
static bool eval_field(float* out, const CavityData& d, const float p[3], const float n[3], bool kernel, bool power, uint32_t* nearest = nullptr) {
    const CellGrid& g = d.lists;
    int c[3];
    const int ci = g.cell_of(p, c);
    if (ci < 0) return false;
    float best = FLT_MAX, num = 0, den = 0;
    uint32_t bj = UINT32_MAX;
    for (uint32_t k = g.offset[ci]; k < g.offset[ci + 1]; ++k) {
        const uint32_t j = g.items[k];
        const float* nj = &d.nrm[3 * j];
        if (nj[0] * n[0] + nj[1] * n[1] + nj[2] * n[2] <= 0) continue;
        const float dd = dist2(p, &d.pos[3 * j]);
        const float v = power ? (dd < POWER_REACH * POWER_REACH ? dd - d.weight[j] : FLT_MAX) : dd;
        if (v < best) { best = v; bj = j; }
        const float s = d.sigma[j];
        if (kernel && dd < KERNEL_CUT * KERNEL_CUT * s * s) {
            const float K = expf(-0.5f * dd / (s * s));
            num += K * d.coef[2 * j];
            den += K * d.coef[2 * j + 1];
        }
    }
    if (bj == UINT32_MAX) return false;
    if (nearest) *nearest = bj;
    if (!kernel) *out = d.nearest_value[bj];
    else if (d.conservative) *out = num;
    else *out = den > 0 ? num / den : d.nearest_value[bj];
    return true;
}

// What the field spans over the exposed surface: the 1st to 99th percentile, by area
// Of a field of the fitted partition: its integral over the surface, by the quadrature at its own points through the
// lists of the cells as the shader reads them, against the sum of the values over the parts' areas. A list that left
// out a tessera shows here, and so does how far the fit is from the areas. Of Area coverage: its integral against the
// sum of the areas it spreads.
static float kernel_error(const CavityData& d) {
    if (d.field < 0 || (!d.conserve && !d.conservative)) return 0.0f;
    double integral = 0, total = 0, covered = 0;
    for (const QuadPoint& q : d.quad) {
        float f;
        if (eval_field(&f, d, q.p, q.n, true, true)) {
            integral += (double)f * q.da * q.w;
            covered += (double)q.da * q.w;
        }
    }
    double sum_a = 0;
    for (uint32_t j = 0; j < d.num_tess; ++j) if (d.norm[j] > 0) sum_a += d.area[j];
    for (uint32_t j = 0; j < d.num_tess; ++j) {
        double nv, v;
        bool cons;
        tessera_values(d, (CavityField)d.field, j, &nv, &v, &cons);
        if (d.norm[j] > 0) total += cons ? v : v * d.area[j] * (sum_a > 0 ? covered / sum_a : 0.0);
    }
    return total != 0 ? (float)fabs(integral / total - 1.0) : 0.0f;
}

static ColorScaleSpan field_span(const CavityData& d, bool kernel) {
    std::vector<float> v;
    const size_t stride = MD_MAX((size_t)1, d.quad.size() / 40000);
    for (size_t i = 0; i < d.quad.size(); i += stride) {
        float f;
        if (eval_field(&f, d, d.quad[i].p, d.quad[i].n, kernel, true)) v.push_back(f);
    }
    ColorScaleSpan span = {};
    if (v.empty()) return span;
    std::sort(v.begin(), v.end());
    span.valid = true;
    span.min = v.front();
    span.max = v.back();
    span.lo = v[(size_t)(0.01 * (v.size() - 1))];
    span.hi = v[(size_t)(0.99 * (v.size() - 1))];
    return span;
}

// ---------------------------------------------------------------------------------------------------------
// GPU
// ---------------------------------------------------------------------------------------------------------

static void upload_buffer(CavityData* d, int i, GLenum format, const void* data, size_t bytes) {
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

enum { BUF_TESS_A, BUF_TESS_B, BUF_TESS_C, BUF_SPHERES, BUF_CELLS, BUF_ITEMS, BUF_NBR_RANGE, BUF_NBR_ITEMS, BUF_GLYPH, BUF_COUNT };

static void upload_kernel(CavityData* d) {
    const uint32_t n = d->num_tess;
    std::vector<float> a(4 * n);
    for (uint32_t j = 0; j < n; ++j) {
        a[4 * j + 0] = d->pos[3 * j + 0]; a[4 * j + 1] = d->pos[3 * j + 1]; a[4 * j + 2] = d->pos[3 * j + 2];
        a[4 * j + 3] = d->sigma[j];
    }
    upload_buffer(d, BUF_TESS_A, GL_RGBA32F, a.data(), a.size() * sizeof(float));
    upload_buffer(d, BUF_SPHERES, GL_RGBA32F, d->sph.data(), d->sph.size() * sizeof(float));
    const CellGrid& g = d->lists;
    const size_t cells = g.offset.size() - 1;
    std::vector<uint32_t> range(2 * cells);
    for (size_t c = 0; c < cells; ++c) {
        range[2 * c + 0] = g.offset[c];
        range[2 * c + 1] = g.offset[c + 1] - g.offset[c];
    }
    upload_buffer(d, BUF_CELLS, GL_RG32UI, range.data(), range.size() * sizeof(uint32_t));
    upload_buffer(d, BUF_ITEMS, GL_R32UI, g.items.data(), g.items.size() * sizeof(uint32_t));
    std::vector<uint32_t> nr(2 * d->num_spheres);
    for (uint32_t i = 0; i < d->num_spheres; ++i) {
        nr[2 * i + 0] = d->nbr_off[i];
        nr[2 * i + 1] = d->nbr_off[i + 1] - d->nbr_off[i];
    }
    upload_buffer(d, BUF_NBR_RANGE, GL_RG32UI, nr.data(), nr.size() * sizeof(uint32_t));
    upload_buffer(d, BUF_NBR_ITEMS, GL_R32UI, d->nbr.data(), d->nbr.size() * sizeof(uint32_t));

    // Of the glyphs: the radius where the Gaussian is at half its peak, exp(-zeta r^2) = 1/2, which with zeta a = pi
    // ln 2 is sqrt(a / pi); the decay length; its weight; whether it is on the surface. The weight is the switching at
    // its point where the cavity is switched: the file's areas are the Lebedev weights times it, a_j = w_j 4 pi R^2 F_j
    // (on 5helicene.h5 to 4e-13), so this is F alone, without the Lebedev weight. Else the share of a whole point.
    std::vector<float> gl(4 * n);
    for (uint32_t j = 0; j < n; ++j) {
        const double R = d->sph[4 * d->atom[j] + 3];
        const double full = d->full_share * 4.0 * M_PI * R * R;
        gl[4 * j + 0] = d->zeta[j] > 0 ? (float)sqrt(log(2.0) / d->zeta[j]) : (float)sqrt(d->area[j] / M_PI);
        gl[4 * j + 1] = (float)d->lambda[j];
        if (d->has_switching) {
            const float* s = &d->sph[4 * d->atom[j]];
            float n[3] = {d->pos[3 * j] - s[0], d->pos[3 * j + 1] - s[1], d->pos[3 * j + 2] - s[2]};
            const float len = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            for (float& v : n) v /= len > 0 ? len : 1.0f;
            gl[4 * j + 2] = switching(*d, d->atom[j], n);
        } else {
            gl[4 * j + 2] = full > 0 ? (float)(d->area[j] / full) : 0.0f;
        }
        gl[4 * j + 3] = d->shown[j] ? 1.0f : 0.0f;
    }
    upload_buffer(d, BUF_GLYPH, GL_RGBA32F, gl.data(), gl.size() * sizeof(float));
}

static void upload_field(CavityData* d) {
    const uint32_t n = d->num_tess;
    std::vector<float> b(4 * n);
    for (uint32_t j = 0; j < n; ++j) {
        b[4 * j + 0] = d->nrm[3 * j + 0]; b[4 * j + 1] = d->nrm[3 * j + 1]; b[4 * j + 2] = d->nrm[3 * j + 2];
        b[4 * j + 3] = d->nearest_value[j];
    }
    upload_buffer(d, BUF_TESS_B, GL_RGBA32F, b.data(), b.size() * sizeof(float));
    std::vector<float> c(4 * n);
    for (uint32_t j = 0; j < n; ++j) {
        c[4 * j + 0] = d->coef[2 * j + 0];
        c[4 * j + 1] = d->coef[2 * j + 1];
        c[4 * j + 2] = d->weight.size() == n ? d->weight[j] : 0.0f;
        c[4 * j + 3] = 0.0f;
    }
    upload_buffer(d, BUF_TESS_C, GL_RGBA32F, c.data(), c.size() * sizeof(float));
}

static void free_gpu(CavityData* d) {
    for (int i = 0; i < BUF_COUNT; ++i) {
        if (d->tex[i]) glDeleteTextures(1, &d->tex[i]);
        if (d->buf[i]) glDeleteBuffers(1, &d->buf[i]);
        d->tex[i] = d->buf[i] = 0;
    }
    if (d->colormap_tex) glDeleteTextures(1, &d->colormap_tex);
    d->colormap_tex = 0;
    d->colormap = -1;
}

static GLuint colormap_texture(CavityData* d, int colormap) {
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

// ---------------------------------------------------------------------------------------------------------
// Placement: the cavity of the file's geometry on the atoms as they are
// ---------------------------------------------------------------------------------------------------------

static void place(CavityData* d, const md_system_state_t& state) {
    const uint32_t m = d->num_spheres;
    if (!state.xyz || state.num_atoms == 0) return;
    uint64_t h = 0;
    for (uint32_t i = 0; i < m; ++i) h = md_hash64(&state.xyz[d->sys_idx[i]], sizeof(md_f32x3_t), h);
    if (h == d->placement_hash) return;
    d->placement_hash = h;

    // The rigid motion taking the QM atoms of the file onto where they are now (Kabsch)
    double cf[3] = {0, 0, 0}, cw[3] = {0, 0, 0};
    for (uint32_t i = 0; i < m; ++i) {
        const md_f32x3_t w = state.xyz[d->sys_idx[i]];
        for (int k = 0; k < 3; ++k) cf[k] += d->sph[4 * i + k] / m;
        cw[0] += w.x / m; cw[1] += w.y / m; cw[2] += w.z / m;
    }
    std::vector<md_f32x4_t> ref(m), trg(m);
    float dev = 0;
    for (uint32_t i = 0; i < m; ++i) {
        const md_f32x3_t w = state.xyz[d->sys_idx[i]];
        ref[i] = {(float)(d->sph[4 * i + 0] - cf[0]), (float)(d->sph[4 * i + 1] - cf[1]), (float)(d->sph[4 * i + 2] - cf[2]), 1.0f};
        trg[i] = {(float)(w.x - cw[0]), (float)(w.y - cw[1]), (float)(w.z - cw[2]), 1.0f};
        dev = fmaxf(dev, fabsf(ref[i].x - trg[i].x) + fabsf(ref[i].y - trg[i].y) + fabsf(ref[i].z - trg[i].z));
    }
    // R maps the target onto the reference, R (q - cw) = p - cf; the cavity goes the other way: q = R^T (p - cf) + cw
    float R[9] = {1,0,0, 0,1,0, 0,0,1};
    if (m >= 3 && dev > 1e-4f) {
        const md_f32x3x3_t M = md_util_optimal_rotation_rel_vec4(ref.data(), trg.data(), m);
        for (int c = 0; c < 3; ++c)
            for (int r = 0; r < 3; ++r) R[3 * c + r] = M.col[r].elem[c];   // the transpose, column major
    }
    memcpy(d->rot, R, sizeof(R));
    for (int k = 0; k < 3; ++k) {
        d->trans[k] = (float)(cw[k] - (R[k] * cf[0] + R[3 + k] * cf[1] + R[6 + k] * cf[2]));
    }

    d->version++;
}

// world -> the frame of the file, as a column major 4x4
static md_f32x4x4_t world_to_file(const CavityData& d) {
    md_f32x4x4_t M = {};
    // R^T, and -R^T t
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r) M.elem[c][r] = d.rot[3 * r + c];
    for (int r = 0; r < 3; ++r) {
        M.elem[3][r] = -(d.rot[3 * r + 0] * d.trans[0] + d.rot[3 * r + 1] * d.trans[1] + d.rot[3 * r + 2] * d.trans[2]);
    }
    M.elem[3][3] = 1.0f;
    return M;
}

// ---------------------------------------------------------------------------------------------------------
// Public
// ---------------------------------------------------------------------------------------------------------

bool cavity_update(CavityRepresentation* cav, const md_system_t& sys, const md_system_state_t& state) {
    if (!cavity_available(sys)) return false;
    if (!cav->data) cav->data = new CavityData();
    CavityData* d = cav->data;

    const uint64_t h = source_hash(sys);
    if (h != d->source_hash) {
        free_gpu(d);
        CavityData fresh;
        *d = std::move(fresh);
        if (!read_source(d, sys)) {
            d->num_tess = 0;
            d->source_hash = h;
            return false;
        }
        d->source_hash = h;
    }
    if (d->num_tess == 0) return false;

    // The surface: everything below is made over it, the weights of the patches too
    const bool switched = cav->switched && d->has_switching;
    if (switched != d->switched) {
        d->switched = switched;
        d->weight.clear();
        d->kernel_scale = -1.0f;
    }
    const float min_area = MD_CLAMP(cav->min_area, 0.0f, 0.5f);
    if (min_area != d->min_area) {
        d->min_area = min_area;
        if (d->full_share <= 0) {
            for (uint32_t j = 0; j < d->num_tess; ++j) {
                const double R = d->sph[4 * d->atom[j] + 3];
                if (R > 0) d->full_share = fmax(d->full_share, d->area[j] / (4.0 * M_PI * R * R));
            }
        }
        d->shown.assign(d->num_tess, 1);
        uint32_t hidden = 0;
        for (uint32_t j = 0; j < d->num_tess; ++j) {
            const double R = d->sph[4 * d->atom[j] + 3];
            const double full = d->full_share * 4.0 * M_PI * R * R;
            if (d->area[j] < min_area * full) { d->shown[j] = 0; hidden++; }
        }
        MD_LOG_DEBUG("Cavity: %u of %u tesserae below %.3g of a whole point, left off the surface", hidden, d->num_tess, min_area);
        d->weight.clear();
        d->kernel_scale = -1.0f;
    }
    const float scale = MD_CLAMP(cav->kernel_scale, 0.25f, 4.0f);
    if (scale != d->kernel_scale) {
        build_kernel(d, scale);
        upload_kernel(d);
        d->field = -1;
        d->version++;
    }
    if ((int)cav->field != d->field || cav->conserve != d->conserve) {
        build_field(d, cav->field, cav->conserve);
        upload_field(d);
        d->span_nearest = field_span(*d, false);
        d->span_kernel  = field_span(*d, true);
        d->kernel_error = kernel_error(*d);
        MD_LOG_DEBUG("Cavity: %s, the integral of the kernel off its sum by %.2e", cavity_field_str[d->field], d->kernel_error);
        d->version++;
    }
    place(d, state);
    color_scale_update_range(&cav->scale, cavity_span(*cav));
    return true;
}

void cavity_invalidate(CavityRepresentation* cav) {
    if (cav->data) cav->data->source_hash = 0;
}

void cavity_free(CavityRepresentation* cav) {
    if (!cav->data) return;
    free_gpu(cav->data);
    delete cav->data;
    cav->data = nullptr;
}

ColorScaleSpan cavity_span(const CavityRepresentation& cav) {
    if (!cav.data || cav.data->field != (int)cav.field) return {};
    return cav.coloring == CavityColoring::Kernel || cav.field == CavityField::AreaCoverage ? cav.data->span_kernel : cav.data->span_nearest;
}

uint64_t cavity_hash(const CavityRepresentation& cav) {
    uint64_t h = md_hash64(&cav.field, sizeof(cav.field), 0);
    h = md_hash64(&cav.coloring, sizeof(cav.coloring), h);
    h = md_hash64(&cav.kernel_scale, sizeof(cav.kernel_scale), h);
    h = md_hash64(&cav.unlit, sizeof(cav.unlit), h);
    h = md_hash64(&cav.conserve, sizeof(cav.conserve), h);
    h = md_hash64(&cav.fit_areas, sizeof(cav.fit_areas), h);
    h = md_hash64(&cav.switched, sizeof(cav.switched), h);
    h = md_hash64(&cav.min_area, sizeof(cav.min_area), h);
    h = md_hash64(&cav.edges, sizeof(cav.edges), h);
    h = md_hash64(&cav.edge_width, sizeof(cav.edge_width), h);
    h = md_hash64(&cav.weight_outline, sizeof(cav.weight_outline), h);
    h = md_hash64(&cav.weight_level, sizeof(cav.weight_level), h);
    h = md_hash64(&cav.scale, sizeof(cav.scale), h);
    h = md_hash64(&cav.surface, sizeof(cav.surface), h);
    h = md_hash64(&cav.glyphs, sizeof(cav.glyphs), h);
    h = md_hash64(&cav.glyph_shape, sizeof(cav.glyph_shape), h);
    h = md_hash64(&cav.glyph_set, sizeof(cav.glyph_set), h);
    h = md_hash64(&cav.glyph_length, sizeof(cav.glyph_length), h);
    h = md_hash64(&cav.glyph_color, sizeof(cav.glyph_color), h);
    h = md_hash64(&cav.glyph_length_scale, sizeof(cav.glyph_length_scale), h);
    h = md_hash64(&cav.glyph_radius_halfmax, sizeof(cav.glyph_radius_halfmax), h);
    h = md_hash64(&cav.glyph_radius, sizeof(cav.glyph_radius), h);
    if (cav.data) h = md_hash64(&cav.data->version, sizeof(cav.data->version), h);
    return h;
}

uint32_t cavity_sphere_count(const CavityRepresentation& cav)  { return cav.data ? cav.data->num_spheres : 0; }
uint32_t cavity_tessera_count(const CavityRepresentation& cav) { return cav.data ? cav.data->num_tess : 0; }

bool cavity_switched(const CavityRepresentation& cav) { return cav.data && cav.data->switched; }
bool cavity_has_switching(const CavityRepresentation& cav) { return cav.data && cav.data->has_switching; }

// ---------------------------------------------------------------------------------------------------------
// Passes
// ---------------------------------------------------------------------------------------------------------

static const char* v_shader_fullscreen = R"(
#version 410 core
void main() {
    vec2 uv = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    gl_Position = vec4(uv * 2.0 - 1.0, 0.0, 1.0);
}
)";

static struct {
    GLuint geometry = 0;
    GLuint glyph = 0;
    GLuint resolve = 0;
    GLuint marks = 0;
    GLuint vao = 0;
    GLuint pick_copy = 0;
    GLuint pick_fbo = 0;
    int    pick_w = 0, pick_h = 0;
    bool   failed = false;
} gl_state;

static GLuint make_program(const char* label, const unsigned char* frag, size_t frag_size) {
    GLuint vs = gl::compile_shader_from_source(md_str_from_cstr(v_shader_fullscreen), GL_VERTEX_SHADER);
    GLuint fs = gl::compile_shader_from_source({(const char*)frag, frag_size}, GL_FRAGMENT_SHADER);
    GLuint prog = 0;
    if (vs && fs) {
        prog = glCreateProgram();
        const GLuint shaders[] = {vs, fs};
        if (!gl::attach_link_detach(prog, shaders, 2)) {
            MD_LOG_ERROR("Cavity: the %s program did not link", label);
            glDeleteProgram(prog);
            prog = 0;
        }
    }
    if (vs) glDeleteShader(vs);
    if (fs) glDeleteShader(fs);
    return prog;
}

static GLuint make_program_vf(const char* label, const unsigned char* vert, size_t vert_size, const unsigned char* frag, size_t frag_size) {
    GLuint vs = gl::compile_shader_from_source({(const char*)vert, vert_size}, GL_VERTEX_SHADER);
    GLuint fs = gl::compile_shader_from_source({(const char*)frag, frag_size}, GL_FRAGMENT_SHADER);
    GLuint prog = 0;
    if (vs && fs) {
        prog = glCreateProgram();
        const GLuint shaders[] = {vs, fs};
        if (!gl::attach_link_detach(prog, shaders, 2)) {
            MD_LOG_ERROR("Cavity: the %s program did not link", label);
            glDeleteProgram(prog);
            prog = 0;
        }
    }
    if (vs) glDeleteShader(vs);
    if (fs) glDeleteShader(fs);
    return prog;
}

static bool ensure_programs() {
    if (gl_state.failed) return false;
    if (!gl_state.geometry) gl_state.geometry = make_program_vf("geometry", cavity_geometry_vert, cavity_geometry_vert_size, cavity_geometry_frag, cavity_geometry_frag_size);
    if (!gl_state.glyph)    gl_state.glyph    = make_program_vf("glyph", cavity_glyph_vert, cavity_glyph_vert_size, cavity_glyph_frag, cavity_glyph_frag_size);
    if (!gl_state.resolve) gl_state.resolve = make_program("resolve", cavity_resolve_frag, cavity_resolve_frag_size);
    if (!gl_state.marks)   gl_state.marks   = make_program("marks", cavity_marks_frag, cavity_marks_frag_size);
    if (!gl_state.vao)     glGenVertexArrays(1, &gl_state.vao);
    if (!gl_state.geometry || !gl_state.glyph || !gl_state.resolve || !gl_state.marks) {
        gl_state.failed = true;
        return false;
    }
    return true;
}

// The textures every pass reads, from unit 0: tess_a, tess_b, tess_c, spheres, cells, items
static void bind_buffers(GLuint prog, const CavityData& d) {
    static const char* names[8] = {"u_tess_a", "u_tess_b", "u_tess_c", "u_spheres", "u_cells", "u_items", "u_nbr_range", "u_nbr_items"};
    for (int i = 0; i < 8; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, d.tex[i]);
        glUniform1i(glGetUniformLocation(prog, names[i]), i);
    }
    const CellGrid& g = d.lists;
    glUniform3f(glGetUniformLocation(prog, "u_grid_min"), g.min[0], g.min[1], g.min[2]);
    glUniform1f(glGetUniformLocation(prog, "u_inv_cell"), 1.0f / g.cell);
    glUniform3i(glGetUniformLocation(prog, "u_grid_dim"), g.dim[0], g.dim[1], g.dim[2]);
}

static void unbind_buffers() {
    for (int i = 0; i < 8; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
    }
    glActiveTexture(GL_TEXTURE0);
}

void cavity_resolve(const CavityRepresentation& cav, const GBuffer& gbuf, const ViewParam& view, uint32_t sphere_base, uint32_t tessera_base) {
    const CavityData* d = cav.data;
    if (!cav.surface || !d || d->num_tess == 0 || d->field != (int)cav.field || !ensure_programs()) return;

    // The picking target as it is, read while the G-buffer's is written
    if (!gl_state.pick_copy || gl_state.pick_w != (int)gbuf.width || gl_state.pick_h != (int)gbuf.height) {
        if (!gl_state.pick_copy) glGenTextures(1, &gl_state.pick_copy);
        glBindTexture(GL_TEXTURE_2D, gl_state.pick_copy);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gbuf.width, gbuf.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
        if (!gl_state.pick_fbo) glGenFramebuffers(1, &gl_state.pick_fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl_state.pick_fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl_state.pick_copy, 0);
        gl_state.pick_w = gbuf.width;
        gl_state.pick_h = gbuf.height;
    }
    glBindFramebuffer(GL_READ_FRAMEBUFFER, gbuf.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl_state.pick_fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glBlitFramebuffer(0, 0, gbuf.width, gbuf.height, 0, 0, gbuf.width, gbuf.height, GL_COLOR_BUFFER_BIT, GL_NEAREST);

    // Into the colour and the picking of the G-buffer; the depth is attached but neither tested nor written nor read
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
    const GLenum bufs[2] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_PICKING};
    glDrawBuffers(2, bufs);
    glViewport(0, 0, gbuf.width, gbuf.height);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);

    const GLuint prog = gl_state.resolve;
    glUseProgram(prog);
    bind_buffers(prog, *d);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, gl_state.pick_copy);
    glUniform1i(glGetUniformLocation(prog, "u_tex_picking"), 8);
    glActiveTexture(GL_TEXTURE9);
    glBindTexture(GL_TEXTURE_2D, colormap_texture((CavityData*)d, cav.scale.colormap));
    glUniform1i(glGetUniformLocation(prog, "u_colormap"), 9);

    // The rays of the pixels in the frame of the file, by the same (jittered) projection the spheres were drawn with
    const md_f32x4x4_t clip_to_world = md_f32x4x4_inverse(view.matrix.curr.proj * view.matrix.curr.view);
    const md_f32x4x4_t clip_to_file = world_to_file(*d) * clip_to_world;
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_clip_to_cavity"), 1, GL_FALSE, &clip_to_file.elem[0][0]);
    glUniform2f(glGetUniformLocation(prog, "u_inv_size"), 1.0f / gbuf.width, 1.0f / gbuf.height);
    glUniform1ui(glGetUniformLocation(prog, "u_sphere_base"), sphere_base);
    glUniform1ui(glGetUniformLocation(prog, "u_sphere_count"), d->num_spheres);
    glUniform1ui(glGetUniformLocation(prog, "u_tessera_base"), tessera_base);
    const int mode = cav.coloring == CavityColoring::Nearest && cav.field != CavityField::AreaCoverage ? 0 : (d->conservative ? 2 : 1);
    glUniform1i(glGetUniformLocation(prog, "u_mode"), mode);
    glUniform2f(glGetUniformLocation(prog, "u_range"), cav.scale.range_beg, cav.scale.range_end);
    glUniform1i(glGetUniformLocation(prog, "u_unlit"), cav.unlit ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_edges"), cav.edges && mode == 0 ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_edge_width"), MD_CLAMP(cav.edge_width, 0.001f, 0.2f));
    glUniform1f(glGetUniformLocation(prog, "u_level"), MD_CLAMP(cav.weight_level, 0.01f, 1.0f));
    glUniform1i(glGetUniformLocation(prog, "u_outline"), cav.weight_outline ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_power"), cav.fit_areas ? 1.0f : 0.0f);
    glUniform1f(glGetUniformLocation(prog, "u_power_reach2"), POWER_REACH * POWER_REACH);
    glUniform1i(glGetUniformLocation(prog, "u_switched"), d->switched ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_gamma"), d->sw_gamma);
    glUniform1f(glGetUniformLocation(prog, "u_alpha"), d->sw_alpha);
    glUniform1f(glGetUniformLocation(prog, "u_r_ext"), d->sw_r_ext);
    glUniform1f(glGetUniformLocation(prog, "u_screen_period"), SCREEN_PERIOD);

    glBindVertexArray(gl_state.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    unbind_buffers();
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE9);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glEnable(GL_CULL_FACE);
}

static void draw_glyph_marks(const CavityRepresentation& cav, const CavityData& d, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view);

void cavity_draw_marks(const CavityRepresentation& cav, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view) {
    const CavityData* d = cav.data;
    if (!d || tessera >= d->num_tess || d->field != (int)cav.field || !ensure_programs()) return;

    // The sphere of the tessera first, then those that overlap it: the patch is on them, or within them
    const uint32_t own = d->atom[tessera];
    const float* c = &d->sph[4 * own];
    int local[32];
    int num_local = 0;
    local[num_local++] = (int)own;
    for (uint32_t k = 0; k < d->num_spheres && num_local < 32; ++k) {
        const float Rk = d->sph[4 * k + 3];
        if (k == own || Rk <= 0) continue;
        if (dist2(c, &d->sph[4 * k]) < (c[3] + Rk) * (c[3] + Rk)) local[num_local++] = (int)k;
    }

    // The marks are of the image shown: unjittered, at its size (gbuffer_marks_begin bound the target and blending)
    glViewport(0, 0, gbuf.out_width, gbuf.out_height);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    if (cav.glyphs) draw_glyph_marks(cav, *d, tessera, gbuf, view);
    if (!cav.surface) {
        glEnable(GL_CULL_FACE);
        return;
    }

    const GLuint prog = gl_state.marks;
    glUseProgram(prog);
    bind_buffers(prog, *d);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, gbuf.tex.depth);
    glUniform1i(glGetUniformLocation(prog, "u_tex_depth"), 8);

    const md_f32x4x4_t to_file = world_to_file(*d);
    const md_f32x4x4_t file_to_world = md_f32x4x4_inverse(to_file);
    const md_f32x4x4_t clip_to_world = md_f32x4x4_inverse(view.matrix.curr.proj_no_jitter * view.matrix.curr.view);
    const md_f32x4x4_t clip_to_file = to_file * clip_to_world;
    const md_f32x4x4_t file_to_view = view.matrix.curr.view * file_to_world;
    const md_f32x4x4_t clip_to_view = md_f32x4x4_inverse(view.matrix.curr.proj);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_clip_to_cavity"), 1, GL_FALSE, &clip_to_file.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_cavity_to_view"), 1, GL_FALSE, &file_to_view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_clip_to_view"), 1, GL_FALSE, &clip_to_view.elem[0][0]);
    glUniform1f(glGetUniformLocation(prog, "u_proj_11"), view.matrix.curr.proj.elem[1][1]);
    glUniform1i(glGetUniformLocation(prog, "u_ortho"), view.matrix.curr.proj.elem[3][3] == 1.0f ? 1 : 0);
    glUniform2f(glGetUniformLocation(prog, "u_inv_size"), 1.0f / gbuf.out_width, 1.0f / gbuf.out_height);
    glUniform1iv(glGetUniformLocation(prog, "u_local"), num_local, local);
    glUniform1i(glGetUniformLocation(prog, "u_num_local"), num_local);
    glUniform1ui(glGetUniformLocation(prog, "u_tessera"), tessera);
    glUniform1f(glGetUniformLocation(prog, "u_power"), cav.fit_areas ? 1.0f : 0.0f);
    glUniform1f(glGetUniformLocation(prog, "u_power_reach2"), POWER_REACH * POWER_REACH);

    glBindVertexArray(gl_state.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    unbind_buffers();
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
    glEnable(GL_CULL_FACE);
}

static void jitter_uv(float out[2], const md_f32x4x4_t& P) {
    // As md_gl reads it off the projection
    if (P.elem[2][3] == 0.0f) {
        out[0] = -P.elem[3][0] * 0.5f;
        out[1] = -P.elem[3][1] * 0.5f;
    } else {
        out[0] = P.elem[2][0] * 0.5f;
        out[1] = P.elem[2][1] * 0.5f;
    }
}

void cavity_draw_geometry(const CavityRepresentation& cav, const ViewParam& view, uint32_t sphere_base) {
    const CavityData* d = cav.data;
    if (!cav.surface || !d || d->num_spheres == 0 || !ensure_programs()) return;

    const GLuint prog = gl_state.geometry;
    glUseProgram(prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, d->tex[BUF_SPHERES]);
    glUniform1i(glGetUniformLocation(prog, "u_spheres"), 0);
    const md_f32x4x4_t to_file = world_to_file(*d);
    const md_f32x4x4_t file_to_view = view.matrix.curr.view * md_f32x4x4_inverse(to_file);
    const md_f32x4x4_t curr_view_to_prev_clip = view.matrix.prev.proj * view.matrix.prev.view * view.matrix.inv.view;
    float jit[4];
    jitter_uv(jit + 0, view.matrix.curr.proj);
    jitter_uv(jit + 2, view.matrix.prev.proj);
    const bool ortho = view.matrix.curr.proj.elem[2][3] == 0.0f;

    glUniformMatrix4fv(glGetUniformLocation(prog, "u_cavity_to_view"), 1, GL_FALSE, &file_to_view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_to_clip"), 1, GL_FALSE, &view.matrix.curr.proj.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip"), 1, GL_FALSE, &curr_view_to_prev_clip.elem[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, jit);
    glUniform1i(glGetUniformLocation(prog, "u_ortho"), ortho ? 1 : 0);
    glUniform1ui(glGetUniformLocation(prog, "u_sphere_base"), sphere_base);

    // Into the G-buffer bound, as the opaque representations draw: depth tested and written
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)d->num_spheres);
    glBindVertexArray(0);
    glEnable(GL_CULL_FACE);

    glBindTexture(GL_TEXTURE_BUFFER, 0);
    glUseProgram(0);
}

static const float GLYPH_MIN_PX = 1.5f;     // pixels: the radius a glyph is never thinner than, to be seen and hovered

// The program of the glyphs with its buffers and uniforms, for the G-buffer (proj, jittered) or the marks (unjittered)
static void glyph_program(const CavityRepresentation& cav, const CavityData& d, const ViewParam& view, const md_f32x4x4_t& proj, bool marks) {
    const GLuint prog = gl_state.glyph;
    glUseProgram(prog);
    const GLuint texs[3] = {d.tex[BUF_TESS_A], d.tex[BUF_TESS_B], d.tex[BUF_GLYPH]};
    const char* names[3] = {"u_tess_a", "u_tess_b", "u_glyph"};
    for (int i = 0; i < 3; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, texs[i]);
        glUniform1i(glGetUniformLocation(prog, names[i]), i);
    }
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, colormap_texture((CavityData*)&d, cav.scale.colormap));
    glUniform1i(glGetUniformLocation(prog, "u_colormap"), 3);
    // The depth of the marks on a unit of its own even when unused: two types of sampler on one unit fail the draw
    glUniform1i(glGetUniformLocation(prog, "u_tex_depth"), 8);

    const md_f32x4x4_t to_file = world_to_file(d);
    const md_f32x4x4_t file_to_view = view.matrix.curr.view * md_f32x4x4_inverse(to_file);
    const md_f32x4x4_t curr_view_to_prev_clip = view.matrix.prev.proj * view.matrix.prev.view * view.matrix.inv.view;
    float jit[4];
    jitter_uv(jit + 0, view.matrix.curr.proj);
    jitter_uv(jit + 2, view.matrix.prev.proj);
    const bool ortho = proj.elem[2][3] == 0.0f;
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    const float p11 = proj.elem[1][1];

    glUniformMatrix4fv(glGetUniformLocation(prog, "u_cavity_to_view"), 1, GL_FALSE, &file_to_view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_to_clip"), 1, GL_FALSE, &proj.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip"), 1, GL_FALSE, &curr_view_to_prev_clip.elem[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, jit);
    glUniform1i(glGetUniformLocation(prog, "u_ortho"), ortho ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_only"), -1);
    glUniform1i(glGetUniformLocation(prog, "u_set"), (int)cav.glyph_set);
    glUniform1i(glGetUniformLocation(prog, "u_length_mode"), (int)cav.glyph_length);
    glUniform1f(glGetUniformLocation(prog, "u_length_scale"), cav.glyph_length_scale);
    glUniform1i(glGetUniformLocation(prog, "u_radius_halfmax"), cav.glyph_radius_halfmax ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_radius"), cav.glyph_radius);
    glUniform2f(glGetUniformLocation(prog, "u_range"), cav.scale.range_beg, cav.scale.range_end);
    glUniform1f(glGetUniformLocation(prog, "u_min_px"), GLYPH_MIN_PX);
    glUniform1f(glGetUniformLocation(prog, "u_pixel"), 2.0f / (p11 * (float)MD_MAX(vp[3], 1)));
    glUniform1i(glGetUniformLocation(prog, "u_capsule"), cav.glyph_shape == CavityGlyphShape::Capsule ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_color_mode"), (int)cav.glyph_color);
    glUniform1i(glGetUniformLocation(prog, "u_unlit"), cav.unlit ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_marks"), marks ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_proj_11"), p11);
}

static void glyph_unbind() {
    for (int i = 0; i < 3; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
    }
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
}

void cavity_draw_glyphs(const CavityRepresentation& cav, const ViewParam& view, uint32_t tessera_base) {
    const CavityData* d = cav.data;
    if (!cav.glyphs || !d || d->num_tess == 0 || d->field != (int)cav.field || !ensure_programs()) return;
    glyph_program(cav, *d, view, view.matrix.curr.proj, false);
    glUniform1ui(glGetUniformLocation(gl_state.glyph, "u_tessera_base"), tessera_base);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)d->num_tess);
    glBindVertexArray(0);
    glEnable(GL_CULL_FACE);
    glyph_unbind();
}

// The glyph of a tessera into the marks bound, as cavity_draw_marks marks its patch
static void draw_glyph_marks(const CavityRepresentation& cav, const CavityData& d, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view) {
    glyph_program(cav, d, view, view.matrix.curr.proj_no_jitter, true);
    const GLuint prog = gl_state.glyph;
    glUniform1i(glGetUniformLocation(prog, "u_only"), (int)tessera);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, gbuf.tex.depth);
    glUniform1i(glGetUniformLocation(prog, "u_tex_depth"), 8);
    const md_f32x4x4_t clip_to_view = md_f32x4x4_inverse(view.matrix.curr.proj);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_clip_to_view"), 1, GL_FALSE, &clip_to_view.elem[0][0]);
    glUniform2f(glGetUniformLocation(prog, "u_inv_size"), 1.0f / gbuf.out_width, 1.0f / gbuf.out_height);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
    glBindVertexArray(0);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, 0);
    glyph_unbind();
}

bool cavity_tessera_info(CavityTesseraInfo* info, const CavityRepresentation& cav, uint32_t tessera, const float* world_pos) {
    const CavityData* d = cav.data;
    if (!d || tessera >= d->num_tess) return false;
    *info = {};
    const double a = d->area[tessera], lam = d->lambda[tessera], p = d->amp[tessera], P = d->pressure;
    info->atom = d->atom[tessera];
    info->sys_atom = d->sys_idx[d->atom[tessera]];
    info->area = (float)a;
    info->decay_length = (float)lam;
    info->wall_energy = (float)(GPA_A3_KJMOL * P * lam);
    info->wall_density = p > 0 ? (float)((P / AU_PRESSURE_GPA) * (a / (BOHR_A * BOHR_A)) * (lam / BOHR_A) / p) : 0.0f;
    info->amplitude = (float)p;

    if (world_pos && d->field == (int)cav.field) {
        // On the sphere of the tessera, under the cursor: in the frame of the file
        const md_f32x4x4_t M = world_to_file(*d);
        float x[3];
        for (int r = 0; r < 3; ++r) x[r] = M.elem[0][r] * world_pos[0] + M.elem[1][r] * world_pos[1] + M.elem[2][r] * world_pos[2] + M.elem[3][r];
        // The sphere it is on: the patch of a tessera may reach onto a neighbouring sphere
        uint32_t on = d->atom[tessera];
        float best = FLT_MAX;
        for (uint32_t i = 0; i < d->num_spheres; ++i) {
            const float R = d->sph[4 * i + 3];
            if (R <= 0) continue;
            const float e = fabsf(sqrtf(dist2(x, &d->sph[4 * i])) - R);
            if (e < best) { best = e; on = i; }
        }
        const float* s = &d->sph[4 * on];
        float n[3] = {x[0] - s[0], x[1] - s[1], x[2] - s[2]};
        const float len = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        // On the surface, not on a glyph standing off it (within what the depth buffer resolves)
        if (len > 0 && best < 0.08f) {
            for (int k = 0; k < 3; ++k) { n[k] /= len; x[k] = s[k] + s[3] * n[k]; }
            const bool kernel = cav.coloring == CavityColoring::Kernel || cav.field == CavityField::AreaCoverage;
            info->has_field_here = eval_field(&info->field_here, *d, x, n, kernel, cav.fit_areas);
            if (d->switched) {
                info->switching_here = switching(*d, on, n);
                info->has_switching_here = true;
            }
        }
    }
    return true;
}

bool cavity_summary(CavitySummary* out, const CavityRepresentation& cav) {
    const CavityData* d = cav.data;
    if (!d || d->num_tess == 0) return false;
    *out = {};
    out->pressure_gpa = (float)d->pressure;
    out->num_tesserae = d->num_tess;
    for (uint32_t i = 0; i < d->num_spheres; ++i) out->num_spheres += d->sph[4 * i + 3] > 0;
    double A = 0, E = 0;
    for (uint32_t j = 0; j < d->num_tess; ++j) {
        A += d->area[j];
        E += GPA_A3_KJMOL * d->pressure * d->area[j] * d->lambda[j];
    }
    out->tessera_area = (float)A;
    out->exposed_area = (float)d->exposed_area;
    out->weighted_area = (float)d->weighted_area;
    out->switched = d->switched;
    out->wall_energy = (float)E;
    out->max_kernel_error = d->kernel_error;

    double buried_area = 0, unreached = 0, hidden_area = 0, hidden_energy = 0;
    for (uint32_t j = 0; j < d->num_tess; ++j) {
        const float* p = &d->pos[3 * j];
        for (uint32_t k = 0; k < d->num_spheres; ++k) {
            const float Rk = d->sph[4 * k + 3];
            if (k == d->atom[j] || Rk <= 0) continue;
            if (dist2(p, &d->sph[4 * k]) < Rk * Rk * 0.9999f) {
                out->num_buried += 1;
                buried_area += d->area[j];
                break;
            }
        }
        const double e = GPA_A3_KJMOL * d->pressure * d->area[j] * d->lambda[j];
        if (d->shown.size() == d->num_tess && !d->shown[j]) {
            out->num_hidden += 1;
            hidden_area += d->area[j];
            hidden_energy += e;
        } else if (d->norm.size() == d->num_tess && d->norm[j] <= 0) {
            unreached += e;
        }
    }
    out->hidden_area = A > 0 ? (float)(hidden_area / A) : 0.0f;
    out->hidden_energy = E > 0 ? (float)(hidden_energy / E) : 0.0f;
    out->buried_area = A > 0 ? (float)(buried_area / A) : 0.0f;
    out->unreached_energy = E > 0 ? (float)(unreached / E) : 0.0f;
    return true;
}
