#include "density_volume_core.h"

#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <md_mcs.h>
#include <md_system.h>
#include <md_element.h>
#include <md_infer.h>
#include <core/md_geom.h>
#include <md_xyz.h>

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace density_volume {

static const double BOHR_IN_ANGSTROM = 0.529177210903;

static bool fail(char* err, size_t cap, const char* fmt, ...) {
    if (err && cap) {
        va_list args;
        va_start(args, fmt);
        vsnprintf(err, cap, fmt, args);
        va_end(args);
    }
    return false;
}

// ### CUBE ###

struct Reader {
    const char* p;
    const char* end;
};

static md_str_t next_line(Reader* r) {
    const char* b = r->p;
    while (r->p < r->end && *r->p != '\n') ++r->p;
    md_str_t s = {b, (size_t)(r->p - b)};
    if (r->p < r->end) ++r->p;
    if (s.len && s.ptr[s.len - 1] == '\r') s.len -= 1;
    return s;
}

// The next number of what is left, across lines
static bool next_number(Reader* r, double* out) {
    while (r->p < r->end && (*r->p == ' ' || *r->p == '\t' || *r->p == '\n' || *r->p == '\r')) ++r->p;
    if (r->p >= r->end) return false;
    char buf[64];
    size_t n = 0;
    while (r->p < r->end && *r->p != ' ' && *r->p != '\t' && *r->p != '\n' && *r->p != '\r') {
        if (n + 1 < sizeof(buf)) buf[n++] = *r->p;
        ++r->p;
    }
    buf[n] = 0;
    char* e = nullptr;
    *out = strtod(buf, &e);
    return e && *e == 0 && n > 0;
}

// The numbers of one line
static size_t line_numbers(md_str_t line, double* out, size_t cap) {
    Reader r = {line.ptr, line.ptr + line.len};
    size_t n = 0;
    double x;
    while (n < cap && next_number(&r, &x)) out[n++] = x;
    return n;
}

bool cube_parse(Cube* out, md_str_t text, md_allocator_i* alloc, char* err, size_t err_cap) {
    *out = Cube{};
    out->bulk = NAN;
    out->alloc = alloc;
    Reader r = {text.ptr, text.ptr + text.len};
    const md_str_t c1 = next_line(&r);
    next_line(&r);
    // viamd's own: its bulk density in the first line
    for (size_t i = 0; i + 5 <= c1.len; ++i) {
        if (strncmp(c1.ptr + i, "BULK ", 5) != 0) continue;
        double v[1];
        if (line_numbers({c1.ptr + i + 5, c1.len - i - 5}, v, 1) == 1) out->bulk = v[0];
        break;
    }
    double h[5];
    if (line_numbers(next_line(&r), h, 5) < 4) return fail(err, err_cap, "not a cube file: no atom count and origin on line 3");
    const long long natoms_signed = (long long)h[0];
    const size_t natoms = (size_t)llabs(natoms_signed);
    const bool has_ids = natoms_signed < 0;
    md_float3_t origin_bohr = {(float)h[1], (float)h[2], (float)h[3]};
    double scale[3];
    for (int a = 0; a < 3; ++a) {
        double v[4];
        if (line_numbers(next_line(&r), v, 4) != 4) return fail(err, err_cap, "not a cube file: line %d is not a grid axis", 4 + a);
        const long long n = (long long)v[0];
        if (n == 0 || llabs(n) > (1 << 16)) return fail(err, err_cap, "a grid axis of %lld points", n);
        out->grid.dim[a] = (uint32_t)llabs(n);
        // A negative count: the axis is in Å, else in bohr
        scale[a] = n < 0 ? 1.0 : BOHR_IN_ANGSTROM;
        out->grid.basis[a] = md_float3_t{(float)(v[1] * scale[a]), (float)(v[2] * scale[a]), (float)(v[3] * scale[a])};
    }
    // Positions in the units of the first axis, as readers do
    const double len_scale = scale[0];
    const md_float3_t origin = md_float3_mul1(origin_bohr, (float)len_scale);
    out->num_atoms = natoms;
    out->atomic_numbers = natoms ? (int*)md_alloc(alloc, natoms * sizeof(int)) : nullptr;
    out->xyz = natoms ? (md_float3_t*)md_alloc(alloc, natoms * sizeof(md_float3_t)) : nullptr;
    for (size_t i = 0; i < natoms; ++i) {
        double v[5];
        if (line_numbers(next_line(&r), v, 5) != 5) return fail(err, err_cap, "atom %zu is not 'Z charge x y z'", i + 1);
        out->atomic_numbers[i] = (int)v[0];
        out->xyz[i] = md_float3_t{(float)(v[2] * len_scale), (float)(v[3] * len_scale), (float)(v[4] * len_scale)};
    }
    if (has_ids) {
        // The count of values per point and their ids: one is all that is read
        double v[64];
        const size_t n = line_numbers(next_line(&r), v, 64);
        if (n < 1 || v[0] < 1) return fail(err, err_cap, "no count of values after the atoms");
        if (v[0] > 1) return fail(err, err_cap, "a cube of %d values per point; one is read", (int)v[0]);
    }
    const size_t count = (size_t)out->grid.dim[0] * out->grid.dim[1] * out->grid.dim[2];
    out->values = (float*)md_alloc(alloc, count * sizeof(float));
    for (size_t i = 0; i < count; ++i) {
        double x;
        if (!next_number(&r, &x)) return fail(err, err_cap, "%zu values of %zu", i, count);
        out->values[i] = (float)x;
    }
    // A cube's origin is the first point; a voxel is centred on its point
    out->grid.origin = origin;
    for (int a = 0; a < 3; ++a) out->grid.origin = md_float3_sub(out->grid.origin, md_float3_mul1(out->grid.basis[a], 0.5f));
    out->grid.values = out->values;
    return true;
}

bool cube_read(Cube* out, md_str_t path, md_allocator_i* alloc, char* err, size_t err_cap) {
    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    const md_str_t text = md_load_textfile(path, temp);
    bool ok = false;
    if (!text.len) {
        *out = Cube{};
        fail(err, err_cap, "could not read '%.*s'", (int)path.len, path.ptr);
    } else {
        ok = cube_parse(out, text, alloc, err, err_cap);
    }
    md_arena_allocator_destroy(temp);
    return ok;
}

void cube_free(Cube* cube) {
    if (!cube || !cube->alloc) return;
    if (cube->values) md_free(cube->alloc, cube->values, (size_t)cube->grid.dim[0] * cube->grid.dim[1] * cube->grid.dim[2] * sizeof(float));
    if (cube->atomic_numbers) md_free(cube->alloc, cube->atomic_numbers, cube->num_atoms * sizeof(int));
    if (cube->xyz) md_free(cube->alloc, cube->xyz, cube->num_atoms * sizeof(md_float3_t));
    *cube = Cube{};
}

// ### ALIGNMENT ###

bool align(Alignment* out, const md_system_t* sys, const int32_t* atoms, size_t count, const md_float3_t* xyz, const Cube& other,
           md_allocator_i* alloc, char* err, size_t err_cap) {
    *out = Alignment{};
    if (!sys || !atoms || !count || !xyz) return fail(err, err_cap, "this volume has no reference structure");
    if (!other.num_atoms) return fail(err, err_cap, "the cube has no atoms to place it by");
    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));

    // The cube's atoms as a system of their own: elements and the bonds their distances give
    md_array(char) text = nullptr;
    {
        char line[128];
        int n = snprintf(line, sizeof(line), "%zu\ncube\n", other.num_atoms);
        md_array_push_array(text, line, (size_t)n, temp);
        for (size_t i = 0; i < other.num_atoms; ++i) {
            const md_str_t sym = md_atomic_number_symbol((md_atomic_number_t)other.atomic_numbers[i]);
            n = snprintf(line, sizeof(line), "%.*s %.6f %.6f %.6f\n", (int)sym.len, sym.ptr, other.xyz[i].x, other.xyz[i].y, other.xyz[i].z);
            md_array_push_array(text, line, (size_t)n, temp);
        }
    }
    md_system_t osys = {.alloc = temp};
    md_system_state_t ostate = {.alloc = temp};
    bool ok = md_xyz_system_init_from_str(&osys, &ostate, {text, md_array_size(text)}, NULL) &&
              md_system_infer(&osys, &ostate, MD_SYSTEM_INFER_BOND_BIT, NULL);
    if (!ok) {
        md_arena_allocator_destroy(temp);
        return fail(err, err_cap, "the cube's atoms could not be read as a structure");
    }
    int32_t* oatoms = (int32_t*)md_alloc(temp, other.num_atoms * sizeof(int32_t));
    for (size_t i = 0; i < other.num_atoms; ++i) oatoms[i] = (int32_t)i;

    md_mcs_desc_t desc = {};
    desc.a = {sys, atoms, count, xyz};
    desc.b = {&osys, oatoms, other.num_atoms, (const md_float3_t*)ostate.xyz};
    md_mcs_result_t res = {};
    ok = md_mcs(&res, &desc, temp);
    if (ok) {
        out->matched = res.count;
        out->this_heavy = res.a_heavy;
        out->other_heavy = res.b_heavy;
        out->incomplete = res.incomplete;
    }
    if (ok && res.count < 3) {
        ok = fail(err, err_cap, "%zu heavy atoms in common: too few to place the other volume by", res.count);
    } else if (!ok) {
        fail(err, err_cap, "the common substructure could not be searched");
    }
    if (ok) {
        // Kabsch by mass over the pairs: the other's centre onto this one's, then the rotation
        const size_t n = res.count;
        md_float4_t* p = (md_float4_t*)md_alloc(temp, n * sizeof(md_float4_t));
        md_float4_t* q = (md_float4_t*)md_alloc(temp, n * sizeof(md_float4_t));
        double wsum = 0, cp[3] = {0}, cq[3] = {0};
        for (size_t i = 0; i < n; ++i) {
            const float w = md_atom_mass(&sys->atom, (size_t)res.a[i]);
            const md_float3_t a = xyz[res.a[i]], b = ostate.xyz[res.b[i]];
            p[i] = md_float4_t{a.x, a.y, a.z, w};
            q[i] = md_float4_t{b.x, b.y, b.z, w};
            wsum += w;
            cp[0] += w * a.x; cp[1] += w * a.y; cp[2] += w * a.z;
            cq[0] += w * b.x; cq[1] += w * b.y; cq[2] += w * b.z;
        }
        const md_float3_t comp = {(float)(cp[0] / wsum), (float)(cp[1] / wsum), (float)(cp[2] / wsum)};
        const md_float3_t comq = {(float)(cq[0] / wsum), (float)(cq[1] / wsum), (float)(cq[2] / wsum)};
        for (size_t i = 0; i < n; ++i) {
            p[i].x -= comp.x; p[i].y -= comp.y; p[i].z -= comp.z;
            q[i].x -= comq.x; q[i].y -= comq.y; q[i].z -= comq.z;
        }
        out->R = md_geom_optimal_rotation_rel(p, q, n);
        out->t = md_float3_sub(comp, md_float3x3_mul_float3(out->R, comq));
        double ss = 0;
        for (size_t i = 0; i < n; ++i) {
            const md_float3_t d = md_float3_sub(md_float3x3_mul_float3(out->R, md_float3_t{q[i].x, q[i].y, q[i].z}), md_float3_t{p[i].x, p[i].y, p[i].z});
            ss += p[i].w * md_float3_dot(d, d);
        }
        out->rmsd = sqrt(ss / wsum);
    }
    md_mcs_result_free(&res);
    md_arena_allocator_destroy(temp);
    (void)alloc;
    return ok;
}

// ### RESAMPLING ###

double resample(float* out, const Grid& to, const Grid& from, const Alignment& al) {
    const size_t n = (size_t)to.dim[0] * to.dim[1] * to.dim[2];
    if (!out || !n) return 0.0;
    if (!from.values || !from.dim[0] || !from.dim[1] || !from.dim[2]) {
        memset(out, 0, n * sizeof(float));
        return 0.0;
    }
    // A point of this frame in from's continuous index, whole numbers at voxel centres:
    // f = B^-1 (R^T (p - t) - origin) - 1/2, with B's columns from's steps
    md_float3x3_t B = {};
    for (int a = 0; a < 3; ++a) { B.elem[a][0] = from.basis[a].x; B.elem[a][1] = from.basis[a].y; B.elem[a][2] = from.basis[a].z; }
    const md_float3x3_t Binv = md_float3x3_inverse(B);
    const md_float3x3_t M = md_float3x3_mul(Binv, md_float3x3_transpose(al.R));
    const md_float3_t f0 = md_float3_sub(md_float3x3_mul_float3(Binv, md_float3_sub(md_float3x3_mul_float3(md_float3x3_transpose(al.R), md_float3_mul1(al.t, -1.0f)), from.origin)), md_float3_set1(0.5f));
    // A step of to's index in from's
    md_float3_t step[3];
    for (int a = 0; a < 3; ++a) step[a] = md_float3x3_mul_float3(M, to.basis[a]);
    md_float3_t start = md_float3_add(md_float3x3_mul_float3(M, to.origin), f0);
    for (int a = 0; a < 3; ++a) start = md_float3_add(start, md_float3_mul1(step[a], 0.5f));

    const int64_t d0 = from.dim[0], d1 = from.dim[1], d2 = from.dim[2];
    const float* V = from.values;
    size_t reached = 0, o = 0;
    for (uint32_t i = 0; i < to.dim[0]; ++i) {
        for (uint32_t j = 0; j < to.dim[1]; ++j) {
            md_float3_t f = md_float3_add(start, md_float3_add(md_float3_mul1(step[0], (float)i), md_float3_mul1(step[1], (float)j)));
            for (uint32_t k = 0; k < to.dim[2]; ++k, ++o, f = md_float3_add(f, step[2])) {
                const float c[3] = {f.x, f.y, f.z};
                const int64_t dims[3] = {d0, d1, d2};
                bool inside = true;
                int64_t lo[3];
                float w[3];
                for (int a = 0; a < 3 && inside; ++a) {
                    // Within half a voxel of the outer centres: the outer voxels' values
                    if (!(c[a] >= -0.5f && c[a] <= (float)dims[a] - 0.5f)) { inside = false; break; }
                    const float x = fminf(fmaxf(c[a], 0.0f), (float)(dims[a] - 1));
                    lo[a] = (int64_t)floorf(x);
                    if (lo[a] >= dims[a] - 1) lo[a] = dims[a] - 1;
                    w[a] = x - (float)lo[a];
                }
                if (!inside) { out[o] = 0.0f; continue; }
                reached += 1;
                float v = 0.0f;
                for (int corner = 0; corner < 8; ++corner) {
                    int64_t idx[3];
                    float wt = 1.0f;
                    for (int a = 0; a < 3; ++a) {
                        const int up = (corner >> a) & 1;
                        idx[a] = MD_MIN(lo[a] + up, dims[a] - 1);
                        wt *= up ? w[a] : 1.0f - w[a];
                    }
                    if (wt != 0.0f) v += wt * V[(idx[0] * d1 + idx[1]) * d2 + idx[2]];
                }
                out[o] = v;
            }
        }
    }
    return (double)reached / (double)n;
}

}  // namespace density_volume
