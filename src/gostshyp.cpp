#include <gostshyp.h>

#include <md_system.h>
#include <md_attributes.h>
#include <core/md_geom.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_common.h>

#include <math.h>
#include <float.h>
#include <stdlib.h>

// The attributes it is read from (md_vlx.h)
namespace gostshyp_path {
static const md_str_t coordinate   = MD_STR_INIT("vlx/scf/gostshyp/tessera/coordinate");
static const md_str_t normal       = MD_STR_INIT("vlx/scf/gostshyp/tessera/normal");
static const md_str_t area         = MD_STR_INIT("vlx/scf/gostshyp/tessera/area");
static const md_str_t amplitude    = MD_STR_INIT("vlx/scf/gostshyp/tessera/amplitude");
static const md_str_t exponent     = MD_STR_INIT("vlx/scf/gostshyp/tessera/exponent");
static const md_str_t decay_length = MD_STR_INIT("vlx/scf/gostshyp/tessera/decay_length");
static const md_str_t atom_index   = MD_STR_INIT("vlx/scf/gostshyp/tessera/atom_index");
static const md_str_t switching    = MD_STR_INIT("vlx/scf/gostshyp/tessera/switching");
static const md_str_t lebedev      = MD_STR_INIT("vlx/scf/gostshyp/tessera/lebedev_weight");
static const md_str_t radius       = MD_STR_INIT("vlx/scf/gostshyp/sphere_radius");
static const md_str_t pressure     = MD_STR_INIT("vlx/scf/gostshyp/pressure");
static const md_str_t qm_coord     = MD_STR_INIT("qm/atom/coordinate");
static const md_str_t qm_map       = MD_STR_INIT("qm/atom/system_index");
}

static const double BOHR_A          = 0.52917721092;
static const double GPA_A3_KJMOL    = 0.602214076;          // 1 GPa Angstrom^3 in kJ/mol
static const double AU_PRESSURE_GPA = 29421.02648438959;    // the atomic unit of pressure, Eh / bohr^3

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
    return c && md_attribute_view(c, MD_ATTRIBUTE_TYPE_F64, 3, 1) && find(sys, gostshyp_path::atom_index) && find(sys, gostshyp_path::radius);
}

uint64_t gostshyp_source_hash(const md_system_t& sys) {
    const md_str_t paths[] = {gostshyp_path::coordinate, gostshyp_path::normal, gostshyp_path::area, gostshyp_path::amplitude,
                              gostshyp_path::exponent, gostshyp_path::decay_length, gostshyp_path::atom_index, gostshyp_path::switching,
                              gostshyp_path::lebedev, gostshyp_path::radius, gostshyp_path::pressure, gostshyp_path::qm_coord, gostshyp_path::qm_map};
    uint64_t h = 0x9e3779b97f4a7c15ull;
    for (const md_str_t& p : paths) {
        const md_attribute_t* a = find(sys, p);
        const uint64_t v[2] = {a ? (uint64_t)a->id : 0, a ? a->version : 0};
        h = md_hash64(v, sizeof(v), h);
    }
    return h;
}

bool gostshyp_read(GostshypTesserae* d, const md_system_t& sys) {
    *d = GostshypTesserae{};
    const md_attribute_t* c = find(sys, gostshyp_path::coordinate);
    if (!c || c->format.rank != 1) return false;
    const uint32_t n = c->format.shape[0];
    const double* xyz  = view<double>(sys, gostshyp_path::coordinate,   MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* nrm  = view<double>(sys, gostshyp_path::normal,       MD_ATTRIBUTE_TYPE_F64, 3, n);
    const double* area = view<double>(sys, gostshyp_path::area,         MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* amp  = view<double>(sys, gostshyp_path::amplitude,    MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* zeta = view<double>(sys, gostshyp_path::exponent,     MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* lam  = view<double>(sys, gostshyp_path::decay_length, MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* F    = view<double>(sys, gostshyp_path::switching,    MD_ATTRIBUTE_TYPE_F64, 1, n);
    const double* w    = view<double>(sys, gostshyp_path::lebedev,      MD_ATTRIBUTE_TYPE_F64, 1, n);
    const uint32_t* atom = view<uint32_t>(sys, gostshyp_path::atom_index, MD_ATTRIBUTE_TYPE_U32, 1, n);
    const md_attribute_t* r = find(sys, gostshyp_path::radius);
    const double* rad = r && r->format.rank == 1 ? (const double*)md_attribute_view(r, MD_ATTRIBUTE_TYPE_F64, 1, 1) : nullptr;
    if (!xyz || !nrm || !area || !amp || !zeta || !lam || !F || !w || !atom || !rad) {
        MD_LOG_ERROR("GOSTSHYP: the tesserae are incomplete");
        return false;
    }
    const uint32_t num_atoms = r->format.shape[0];
    const double* qm = view<double>(sys, gostshyp_path::qm_coord, MD_ATTRIBUTE_TYPE_F64, 3, num_atoms);
    const uint32_t* map = view<uint32_t>(sys, gostshyp_path::qm_map, MD_ATTRIBUTE_TYPE_U32, 1, num_atoms);
    if (!qm) {
        MD_LOG_ERROR("GOSTSHYP: the QM atoms are not where the cavity is (qm/atom/coordinate)");
        return false;
    }

    d->num = n;
    d->num_atoms = num_atoms;
    d->pos.assign(xyz, xyz + 3 * n);
    d->nrm.assign(nrm, nrm + 3 * n);
    d->area.assign(area, area + n);
    d->amp.assign(amp, amp + n);
    d->zeta.assign(zeta, zeta + n);
    d->lambda.assign(lam, lam + n);
    d->switching.assign(F, F + n);
    d->lebedev.assign(w, w + n);
    d->atom.assign(atom, atom + n);
    for (uint32_t j = 0; j < n; ++j) {
        if (d->atom[j] >= num_atoms) {
            MD_LOG_ERROR("GOSTSHYP: tessera %u names atom %u of %u", j, d->atom[j], num_atoms);
            *d = GostshypTesserae{};
            return false;
        }
    }
    d->qm_xyz.assign(qm, qm + 3 * num_atoms);
    d->radius.assign(rad, rad + num_atoms);
    d->sys_idx.resize(num_atoms);
    for (uint32_t i = 0; i < num_atoms; ++i) {
        d->sys_idx[i] = map ? map[i] : i;
        if (d->sys_idx[i] >= sys.atom.count) {
            MD_LOG_ERROR("GOSTSHYP: QM atom %u is atom %u of a system of %zu", i, d->sys_idx[i], sys.atom.count);
            *d = GostshypTesserae{};
            return false;
        }
    }
    if (const md_attribute_t* p = find(sys, gostshyp_path::pressure)) {
        if (const double* v = (const double*)md_attribute_view(p, MD_ATTRIBUTE_TYPE_F64, 1, 0)) d->pressure_gpa = *v;
    }
    return true;
}

double gostshyp_energy_per_area(const GostshypTesserae& t, uint32_t j) {
    return GPA_A3_KJMOL * t.pressure_gpa * t.lambda[j];
}

double gostshyp_energy(const GostshypTesserae& t, uint32_t j) {
    return GPA_A3_KJMOL * t.pressure_gpa * t.area[j] * t.lambda[j];
}

double gostshyp_wall_density(const GostshypTesserae& t, uint32_t j) {
    // <rho|g> = P a lambda / p, in atomic units
    const double p = t.amp[j];
    return p > 0 ? (t.pressure_gpa / AU_PRESSURE_GPA) * (t.area[j] / (BOHR_A * BOHR_A)) * (t.lambda[j] / BOHR_A) / p : 0.0;
}

bool gostshyp_place(GostshypPlacement* out, const GostshypTesserae& t, const md_float3_t* xyz, size_t count, const uint32_t* idx, double scale) {
    *out = GostshypPlacement{};
    out->scale = scale;
    const uint32_t m = t.num_atoms;
    if (!xyz || m == 0) return false;
    for (uint32_t i = 0; i < m; ++i) {
        if ((idx ? idx[i] : i) >= count) return false;
    }
    // The rigid motion taking the QM atoms of the file (scaled) onto where they are (Kabsch)
    double cf[3] = {0, 0, 0}, cw[3] = {0, 0, 0};
    for (uint32_t i = 0; i < m; ++i) {
        const md_float3_t w = xyz[idx ? idx[i] : i];
        for (int k = 0; k < 3; ++k) cf[k] += t.qm_xyz[3 * i + k] * scale / m;
        cw[0] += w.x / m; cw[1] += w.y / m; cw[2] += w.z / m;
    }
    md_float4_t* ref = (md_float4_t*)malloc(sizeof(md_float4_t) * m);
    md_float4_t* trg = (md_float4_t*)malloc(sizeof(md_float4_t) * m);
    if (!ref || !trg) { free(ref); free(trg); return false; }
    double dev = 0;
    for (uint32_t i = 0; i < m; ++i) {
        const md_float3_t w = xyz[idx ? idx[i] : i];
        ref[i] = {(float)(t.qm_xyz[3 * i + 0] * scale - cf[0]), (float)(t.qm_xyz[3 * i + 1] * scale - cf[1]), (float)(t.qm_xyz[3 * i + 2] * scale - cf[2]), 1.0f};
        trg[i] = {(float)(w.x - cw[0]), (float)(w.y - cw[1]), (float)(w.z - cw[2]), 1.0f};
        dev = fmax(dev, fabs(ref[i].x - trg[i].x) + fabs(ref[i].y - trg[i].y) + fabs(ref[i].z - trg[i].z));
    }
    // R maps the target onto the reference, R (q - cw) = p - cf; the cavity goes the other way: q = R^T (p - cf) + cw
    if (m >= 3 && dev > 1e-4 * scale) {
        const md_float3x3_t M = md_geom_optimal_rotation_rel(ref, trg, m);
        for (int c = 0; c < 3; ++c)
            for (int r = 0; r < 3; ++r) out->R[3 * c + r] = M.col[r].elem[c];   // the transpose, column major
    }
    free(ref);
    free(trg);
    const double* R = out->R;
    for (int k = 0; k < 3; ++k) out->t[k] = cw[k] - (R[k] * cf[0] + R[3 + k] * cf[1] + R[6 + k] * cf[2]);
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
