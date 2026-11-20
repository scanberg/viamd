#include <gostshyp.h>

#include <md_system.h>
#include <md_attributes.h>
#include <core/md_log.h>
#include <core/md_common.h>

#include <math.h>
#include <float.h>
#include <stdlib.h>

// The attributes it is read from (md_vlx.h)
namespace gostshyp_path {
static const md_str_t coordinate   = MD_STR_INIT("vlx/scf/gostshyp/tessera/coordinate");
static const md_str_t normal       = MD_STR_INIT("vlx/scf/gostshyp/tessera/normal");
static const md_str_t amplitude    = MD_STR_INIT("vlx/scf/gostshyp/tessera/amplitude");
static const md_str_t exponent     = MD_STR_INIT("vlx/scf/gostshyp/tessera/exponent");
static const md_str_t decay_length = MD_STR_INIT("vlx/scf/gostshyp/tessera/decay_length");
static const md_str_t qm_coord     = MD_STR_INIT("qm/atom/coordinate");
}

static const double BOHR_A = 0.52917721092;

static const md_attribute_t* find(const md_system_t& sys, md_str_t path) {
    return md_attributes_find(&sys.attributes, path);
}

template <typename T>
static const T* view(const md_system_t& sys, md_str_t path, md_attribute_type_t type, uint32_t comp, uint32_t count) {
    const md_attribute_t* a = find(sys, path);
    if (!a || a->format.rank != 1 || a->format.shape[0] != count) return nullptr;
    return (const T*)md_attribute_view(a, type, comp, 1);
}

bool gostshyp_available(const md_system_t& sys) {
    const md_attribute_t* c = find(sys, gostshyp_path::coordinate);
    return c && md_attribute_view(c, MD_ATTRIBUTE_TYPE_F64, 3, 1) && find(sys, gostshyp_path::amplitude);
}

bool gostshyp_read(GostshypTesserae* d, const md_system_t& sys) {
    *d = GostshypTesserae{};
    const md_attribute_t* c = find(sys, gostshyp_path::coordinate);
    if (!c || c->format.rank != 1) return false;
    const uint32_t n = c->format.shape[0];
    const double* xyz  = view<double>(sys, gostshyp_path::coordinate,   MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* nrm  = view<double>(sys, gostshyp_path::normal,       MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* amp  = view<double>(sys, gostshyp_path::amplitude,    MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* zeta = view<double>(sys, gostshyp_path::exponent,     MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* lam  = view<double>(sys, gostshyp_path::decay_length, MD_ATTRIBUTE_TYPE_F64, 1, n);
    if (!xyz || !nrm || !amp || !zeta || !lam) {
        MD_LOG_ERROR("GOSTSHYP: the tesserae are incomplete");
        return false;
    }
    const md_attribute_t* q = find(sys, gostshyp_path::qm_coord);
    const uint32_t num_atoms = q && q->format.rank == 1 ? q->format.shape[0] : 0;
    const double* qm = num_atoms ? view<double>(sys, gostshyp_path::qm_coord, MD_ATTRIBUTE_TYPE_F64, 3, num_atoms) : nullptr;
    if (!qm) {
        MD_LOG_ERROR("GOSTSHYP: the QM atoms are not where the cavity is (qm/atom/coordinate)");
        return false;
    }
    d->num = n;
    d->num_atoms = num_atoms;
    d->pos.assign(xyz, xyz + 3 * n);
    d->nrm.assign(nrm, nrm + 3 * n);
    d->amp.assign(amp, amp + n);
    d->zeta.assign(zeta, zeta + n);
    d->lambda.assign(lam, lam + n);
    d->qm_xyz.assign(qm, qm + 3 * num_atoms);
    return true;
}

// ---------------------------------------------------------------------------
// The wall
// ---------------------------------------------------------------------------

// A Gaussian is dropped where exp(-zeta r^2) is below this: its contribution to the potential is then below 1e-9
// of its peak
static const double CUTOFF_EXP = 20.7;      // -ln(1e-9)
static const double CELL_BOHR  = 1.5;

struct GostshypWall {
    uint32_t num = 0;
    std::vector<double> c, n;           // 3 per, bohr, unit
    std::vector<double> zeta, coef;     // bohr^-2; p (zeta / pi)^(3/2), Hartree
    std::vector<double> lambda2z;       // 2 zeta lambda, bohr^-1
    std::vector<double> cut2;           // bohr^2
    double lo[3] = {0, 0, 0};
    int dim[3] = {0, 0, 0};
    std::vector<uint32_t> cell_off;     // [cells + 1]
    std::vector<uint32_t> items;
};

GostshypWall* gostshyp_wall_create(const GostshypTesserae& t, const GostshypPlacement& pl) {
    GostshypWall* w = new GostshypWall();
    const uint32_t n = t.num;
    w->num = n;
    w->c.resize(3 * n); w->n.resize(3 * n);
    w->zeta.resize(n); w->coef.resize(n); w->lambda2z.resize(n); w->cut2.resize(n);
    const double* R = pl.R;
    const double to_bohr = 1.0 / BOHR_A;
    double hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (int k = 0; k < 3; ++k) w->lo[k] = DBL_MAX;
    for (uint32_t j = 0; j < n; ++j) {
        // The file's position, scaled as the placement wants it, moved, and in bohr
        double p[3], q[3];
        for (int k = 0; k < 3; ++k) p[k] = t.pos[3 * j + k] * pl.scale;
        for (int k = 0; k < 3; ++k) q[k] = (R[k] * p[0] + R[3 + k] * p[1] + R[6 + k] * p[2] + pl.t[k]) / pl.scale * to_bohr;
        for (int k = 0; k < 3; ++k) {
            w->c[3 * j + k] = q[k];
            w->n[3 * j + k] = R[k] * t.nrm[3 * j] + R[3 + k] * t.nrm[3 * j + 1] + R[6 + k] * t.nrm[3 * j + 2];
        }
        const double z = t.zeta[j] * BOHR_A * BOHR_A;
        w->zeta[j] = z;
        w->coef[j] = t.amp[j] * pow(z / M_PI, 1.5);
        w->lambda2z[j] = 2.0 * z * t.lambda[j] * to_bohr;
        w->cut2[j] = CUTOFF_EXP / z;
        const double cut = sqrt(w->cut2[j]);
        for (int k = 0; k < 3; ++k) {
            w->lo[k] = fmin(w->lo[k], q[k] - cut);
            hi[k] = fmax(hi[k], q[k] + cut);
        }
    }
    if (n == 0) return w;
    size_t cells = 1;
    for (int k = 0; k < 3; ++k) {
        w->dim[k] = MD_MAX(1, (int)ceil((hi[k] - w->lo[k]) / CELL_BOHR));
        cells *= (size_t)w->dim[k];
    }
    // Each Gaussian into every cell its cutoff's box overlaps
    std::vector<uint32_t> count(cells, 0);
    auto for_cells = [&](uint32_t j, auto f) {
        const double cut = sqrt(w->cut2[j]);
        int a[3], b[3];
        for (int k = 0; k < 3; ++k) {
            a[k] = MD_CLAMP((int)floor((w->c[3 * j + k] - cut - w->lo[k]) / CELL_BOHR), 0, w->dim[k] - 1);
            b[k] = MD_CLAMP((int)floor((w->c[3 * j + k] + cut - w->lo[k]) / CELL_BOHR), 0, w->dim[k] - 1);
        }
        for (int z = a[2]; z <= b[2]; ++z)
            for (int y = a[1]; y <= b[1]; ++y)
                for (int x = a[0]; x <= b[0]; ++x) f(((size_t)z * w->dim[1] + y) * w->dim[0] + x);
    };
    for (uint32_t j = 0; j < n; ++j) for_cells(j, [&](size_t ci) { count[ci]++; });
    w->cell_off.assign(cells + 1, 0);
    for (size_t i = 0; i < cells; ++i) w->cell_off[i + 1] = w->cell_off[i] + count[i];
    w->items.resize(w->cell_off[cells]);
    std::vector<uint32_t> fill(w->cell_off.begin(), w->cell_off.end() - 1);
    for (uint32_t j = 0; j < n; ++j) for_cells(j, [&](size_t ci) { w->items[fill[ci]++] = j; });
    return w;
}

void gostshyp_wall_destroy(GostshypWall* w) {
    delete w;
}

static inline long cell_of(const GostshypWall* w, const float p[3]) {
    int c[3];
    for (int k = 0; k < 3; ++k) {
        c[k] = (int)floor((p[k] - w->lo[k]) / CELL_BOHR);
        if (c[k] < 0 || c[k] >= w->dim[k]) return -1;
    }
    return ((long)c[2] * w->dim[1] + c[1]) * w->dim[0] + c[0];
}

void gostshyp_wall_evaluate(const GostshypWall* w, GostshypField field, const float* xyz, size_t n, double* out) {
    const bool eff = field == GostshypField::EffectivePotential;
    for (size_t i = 0; i < n; ++i) {
        const float* p = xyz + 3 * i;
        double v = 0.0;
        const long ci = w ? cell_of(w, p) : -1;
        if (ci >= 0) {
            for (uint32_t k = w->cell_off[ci]; k < w->cell_off[ci + 1]; ++k) {
                const uint32_t j = w->items[k];
                const double dx = p[0] - w->c[3 * j], dy = p[1] - w->c[3 * j + 1], dz = p[2] - w->c[3 * j + 2];
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > w->cut2[j]) continue;
                double g = w->coef[j] * exp(-w->zeta[j] * r2);
                if (eff) g *= 1.0 + w->lambda2z[j] * (w->n[3 * j] * dx + w->n[3 * j + 1] * dy + w->n[3 * j + 2] * dz);
                v += g;
            }
        }
        out[i] = v;
    }
}

bool gostshyp_wall_reaches(const GostshypWall* w, const float xyz[3]) {
    const long ci = w ? cell_of(w, xyz) : -1;
    if (ci < 0) return false;
    for (uint32_t k = w->cell_off[ci]; k < w->cell_off[ci + 1]; ++k) {
        const uint32_t j = w->items[k];
        const double dx = xyz[0] - w->c[3 * j], dy = xyz[1] - w->c[3 * j + 1], dz = xyz[2] - w->c[3 * j + 2];
        if (dx * dx + dy * dy + dz * dz <= w->cut2[j]) return true;
    }
    return false;
}
