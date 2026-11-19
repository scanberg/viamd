#include <grid_field.h>

#include <md_system.h>
#include <core/md_common.h>
#include <core/md_geom.h>
#include <core/md_log.h>
#include <core/md_unit.h>
#include <core/md_str.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

static const double ANGSTROM_TO_BOHR = 1.0 / 0.52917721092;

bool grid_field_is_grid(const md_attribute_t* attr) {
    if (!attr || attr->format.type != MD_ATTRIBUTE_TYPE_F32 || attr->format.components != 1) return false;
    if (attr->format.rank != 3 && attr->format.rank != 4) return false;
    // cube/<name>: one segment below the group
    const md_str_t prefix = MD_STR_LIT(GRID_FIELD_PREFIX "/");
    if (!md_str_begins_with(attr->path, prefix)) return false;
    return md_str_find_char(md_str_substr(attr->path, prefix.len), '/') == MD_STR_NPOS;
}

md_attribute_id_t grid_field_first(const md_system_t& sys) {
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, MD_STR_LIT(GRID_FIELD_PREFIX)); md_attributes_next(&it);) {
        if (grid_field_is_grid(it.attr)) return it.attr->id;
    }
    return MD_ATTRIBUTE_INVALID;
}

static const md_attribute_t* grid_attr(const md_system_t& sys, md_attribute_id_t key) {
    const md_attribute_t* attr = md_attributes_get(&sys.attributes, key);
    return grid_field_is_grid(attr) ? attr : nullptr;
}

uint32_t grid_field_count(const md_system_t& sys, md_attribute_id_t key) {
    const md_attribute_t* attr = grid_attr(sys, key);
    if (!attr) return 0;
    return attr->format.rank == 4 ? attr->format.shape[0] : 1;
}

// A sibling of the grid's values: <path>/<leaf>
static const md_attribute_t* sibling(const md_system_t& sys, const md_attribute_t* attr, const char* leaf) {
    char path[256];
    const int n = snprintf(path, sizeof(path), "%.*s/%s", (int)attr->path.len, attr->path.ptr, leaf);
    if (n <= 0 || n >= (int)sizeof(path)) return nullptr;
    return md_attributes_find(&sys.attributes, md_str_t{path, (size_t)n});
}

static size_t extract(double* dst, size_t count, const md_system_t& sys, const md_attribute_t* attr, md_unit_t unit) {
    if (!attr || md_attribute_element_count(&attr->format) != count) return 0;
    return md_attribute_extract_f64(dst, count, &sys.attributes, attr, md_attribute_slice_all(), unit);
}

int grid_field_member_id(const md_system_t& sys, md_attribute_id_t key, uint32_t member) {
    const md_attribute_t* attr = grid_attr(sys, key);
    if (!attr) return -1;
    const uint32_t m = grid_field_count(sys, key);
    const md_attribute_t* ids = sibling(sys, attr, "id");
    if (!ids || member >= m) return -1;
    std::vector<double> v(m);
    if (extract(v.data(), m, sys, ids, md_unit_none()) != m) return -1;
    return (int)v[member];
}

void grid_field_free(GridField* f) {
    if (!f) return;
    free(f->values);
    *f = GridField{};
}

// The rigid motion taking the cube's atoms (p) onto where the system's are (q): q = R p + t, by mass alike (unweighted,
// as the atoms are the same). None where they are where the cube has them.
static void place(GridField* f, const md_double3_t* p, const md_float3_t* q, size_t n, md_float3x3_t* R, md_float3_t* t) {
    *R = md_float3x3_ident();
    *t = {};
    if (n == 0) return;
    double cp[3] = {0, 0, 0}, cq[3] = {0, 0, 0};
    for (size_t i = 0; i < n; ++i) {
        cp[0] += p[i].x / n; cp[1] += p[i].y / n; cp[2] += p[i].z / n;
        cq[0] += q[i].x / n; cq[1] += q[i].y / n; cq[2] += q[i].z / n;
    }
    std::vector<md_float4_t> ref(n), trg(n);
    double dev = 0;
    for (size_t i = 0; i < n; ++i) {
        ref[i] = {(float)(q[i].x - cq[0]), (float)(q[i].y - cq[1]), (float)(q[i].z - cq[2]), 1.0f};
        trg[i] = {(float)(p[i].x - cp[0]), (float)(p[i].y - cp[1]), (float)(p[i].z - cp[2]), 1.0f};
        dev = fmax(dev, fabs(ref[i].x - trg[i].x) + fabs(ref[i].y - trg[i].y) + fabs(ref[i].z - trg[i].z));
    }
    // R (p - cp) = q - cq: where the atoms are the cube's, to its precision, translated at most, not rotated by noise
    if (n >= 3 && dev > 1e-4) *R = md_geom_optimal_rotation_rel(ref.data(), trg.data(), n);
    const md_float3_t rc = md_float3x3_mul_float3(*R, md_float3_t{(float)cp[0], (float)cp[1], (float)cp[2]});
    *t = md_float3_t{(float)cq[0] - rc.x, (float)cq[1] - rc.y, (float)cq[2] - rc.z};
    double ss = 0;
    for (size_t i = 0; i < n; ++i) {
        const md_float3_t x = md_float3_add(md_float3x3_mul_float3(*R, md_float3_t{(float)p[i].x, (float)p[i].y, (float)p[i].z}), *t);
        const md_float3_t d = md_float3_sub(x, q[i]);
        ss += md_float3_dot(d, d);
    }
    f->rmsd = sqrt(ss / n);
    f->num_placed_by = n;
}

bool grid_field_get(GridField* out, const md_system_t& sys, const md_system_state_t& state, md_attribute_id_t key, uint32_t member, bool values) {
    *out = GridField{};
    const md_attribute_t* attr = grid_attr(sys, key);
    if (!attr || !attr->data) return false;
    const uint32_t m = grid_field_count(sys, key);
    if (member >= m) return false;
    const uint32_t* shape = attr->format.shape + (attr->format.rank - 3);
    const size_t num = (size_t)shape[0] * shape[1] * shape[2];

    double o[3], b[9];
    if (extract(o, 3, sys, sibling(sys, attr, "origin"), md_unit_angstrom()) != 3 ||
        extract(b, 9, sys, sibling(sys, attr, "basis"), md_unit_angstrom()) != 9) {
        MD_LOG_ERROR("Grid '" MD_STR_FMT "': no origin and basis to place it by", MD_STR_ARG(attr->path));
        return false;
    }

    // Its atoms and the system's they are, as they are now
    md_float3x3_t R;
    md_float3_t t;
    const md_attribute_t* ac = sibling(sys, attr, "atom/coordinate");
    const md_attribute_t* ai = sibling(sys, attr, "atom/system_index");
    const size_t n = ac ? md_attribute_element_count(&ac->format) / 3 : 0;
    std::vector<md_double3_t> p(n);
    std::vector<double> idx(n);
    std::vector<md_float3_t> q(n);
    size_t placed = 0;
    if (n && ai && state.xyz && extract(&p[0].x, 3 * n, sys, ac, md_unit_angstrom()) == 3 * n &&
        extract(idx.data(), n, sys, ai, md_unit_none()) == n) {
        placed = n;
        for (size_t i = 0; i < n; ++i) {
            const size_t j = (size_t)idx[i];
            if (j >= state.num_atoms) { placed = 0; break; }
            q[i] = state.xyz[j];
        }
    }
    place(out, p.data(), q.data(), placed, &R, &t);

    if (values) {
        out->values = (float*)malloc(sizeof(float) * num);
        if (!out->values) return false;
        md_memcpy(out->values, (const float*)attr->data + member * num, sizeof(float) * num);
    }
    for (int a = 0; a < 3; ++a) out->dim[a] = (int)shape[a];
    for (int a = 0; a < 3; ++a) {
        out->axes.col[a] = md_float3x3_mul_float3(R, md_float3_t{(float)b[3 * a], (float)b[3 * a + 1], (float)b[3 * a + 2]});
    }
    out->corner = md_float3_add(md_float3x3_mul_float3(R, md_float3_t{(float)o[0], (float)o[1], (float)o[2]}), t);
    out->world_to_index = md_float3x3_inverse(out->axes);
    return true;
}

float grid_field_sample_bohr(const GridField& f, const float xyz[3]) {
    const float s = (float)(1.0 / ANGSTROM_TO_BOHR);
    const md_float3_t w = {xyz[0] * s - f.corner.x, xyz[1] * s - f.corner.y, xyz[2] * s - f.corner.z};
    const md_float3_t u = md_float3x3_mul_float3(f.world_to_index, w);
    // whole at voxel centres
    const float c[3] = {u.x - 0.5f, u.y - 0.5f, u.z - 0.5f};
    int i0[3];
    float fr[3];
    for (int a = 0; a < 3; ++a) {
        // half a voxel beyond the outermost centres is the grid's edge: clamped to them there, 0 beyond
        if (c[a] < -0.5f || c[a] > (float)f.dim[a] - 0.5f) return 0.0f;
        const float cc = MD_CLAMP(c[a], 0.0f, (float)(f.dim[a] - 1));
        i0[a] = MD_MIN((int)cc, f.dim[a] - 2 < 0 ? 0 : f.dim[a] - 2);
        fr[a] = f.dim[a] > 1 ? cc - (float)i0[a] : 0.0f;
    }
    const size_t sy = (size_t)f.dim[2], sx = (size_t)f.dim[1] * f.dim[2];
    const int d[3] = {f.dim[0] > 1 ? 1 : 0, f.dim[1] > 1 ? 1 : 0, f.dim[2] > 1 ? 1 : 0};
    const float* v = f.values + (size_t)i0[0] * sx + (size_t)i0[1] * sy + (size_t)i0[2];
    float r = 0.0f;
    for (int a = 0; a <= d[0]; ++a) {
        for (int b = 0; b <= d[1]; ++b) {
            for (int k = 0; k <= d[2]; ++k) {
                const float w8 = (a ? fr[0] : 1.0f - fr[0]) * (b ? fr[1] : 1.0f - fr[1]) * (k ? fr[2] : 1.0f - fr[2]);
                r += w8 * v[(size_t)a * sx + (size_t)b * sy + (size_t)k];
            }
        }
    }
    return r;
}

bool grid_field_md_grid(md_grid_t* grid, const GridField& f) {
    md_float3x3_t o = {};
    float len[3];
    for (int a = 0; a < 3; ++a) {
        len[a] = md_float3_length(f.axes.col[a]);
        if (len[a] <= 0.0f) return false;
        o.col[a] = md_float3_mul1(f.axes.col[a], 1.0f / len[a]);
    }
    bool orthogonal = true;
    for (int a = 0; a < 3; ++a) {
        orthogonal &= fabsf(md_float3_dot(o.col[a], o.col[(a + 1) % 3])) < 1.0e-3f;
    }
    if (orthogonal) {
        grid->orientation = o;
        grid->origin  = md_float3_mul1(f.corner, (float)ANGSTROM_TO_BOHR);
        grid->spacing = md_float3_t{len[0] * (float)ANGSTROM_TO_BOHR, len[1] * (float)ANGSTROM_TO_BOHR, len[2] * (float)ANGSTROM_TO_BOHR};
        for (int a = 0; a < 3; ++a) grid->dim[a] = f.dim[a];
        return true;
    }
    // Skewed axes: the box of its corners, along the world's axes, at its finest step
    md_float3_t lo = { FLT_MAX, FLT_MAX, FLT_MAX }, hi = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
    for (int c = 0; c < 8; ++c) {
        md_float3_t p = f.corner;
        for (int a = 0; a < 3; ++a) {
            if (c & (1 << a)) p = md_float3_add(p, md_float3_mul1(f.axes.col[a], (float)f.dim[a]));
        }
        lo = md_float3_min(lo, p);
        hi = md_float3_max(hi, p);
    }
    const float h = MD_MIN(len[0], MD_MIN(len[1], len[2]));
    grid->orientation = md_float3x3_ident();
    grid->origin = md_float3_mul1(lo, (float)ANGSTROM_TO_BOHR);
    for (int a = 0; a < 3; ++a) {
        grid->dim[a] = MD_MAX(1, (int)ceilf((hi.elem[a] - lo.elem[a]) / h));
        grid->spacing.elem[a] = (hi.elem[a] - lo.elem[a]) / grid->dim[a] * (float)ANGSTROM_TO_BOHR;
    }
    return true;
}
