#include "utest.h"

#include "../density_volume_core.h"

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>

#include <math.h>
#include <string.h>
#include <string>
#include <vector>

/* Comparing a volume with one from a cube file: reading it, placing it by the common substructure, sampling it. */

using namespace density_volume;

UTEST(density_volume_core, cube_parse) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    // Two atoms in bohr, a 2 x 3 x 4 grid of 0.5 bohr steps, the values counting up, x outermost and z innermost
    std::string text = "EXPORTED DENSITY VOLUME FROM VIAMD, UNITS IN BOHR, BULK 0.0334\n"
                       "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n"
                       "    2     1.000000     2.000000     3.000000\n"
                       "    2     0.500000     0.000000     0.000000\n"
                       "    3     0.000000     0.500000     0.000000\n"
                       "    4     0.000000     0.000000     0.500000\n"
                       "    8     8.000000     1.000000     0.000000     0.000000\n"
                       "    1     1.000000     0.000000     1.000000     0.000000\n";
    for (int i = 0; i < 24; ++i) text += " " + std::to_string(i) + ((i % 6 == 5) ? "\n" : "");
    Cube c;
    char err[128] = "";
    ASSERT_TRUE(cube_parse(&c, {text.data(), text.size()}, alloc, err, sizeof(err)));
    const float b = 0.529177210903f;
    EXPECT_EQ(2u, c.grid.dim[0]);
    EXPECT_EQ(3u, c.grid.dim[1]);
    EXPECT_EQ(4u, c.grid.dim[2]);
    EXPECT_NEAR(0.5f * b, c.grid.basis[0].x, 1e-6f);
    EXPECT_NEAR(0.5f * b, c.grid.basis[2].z, 1e-6f);
    // The cube's origin is the first point; the grid's is that voxel's corner, half a step back
    EXPECT_NEAR((1.0f - 0.25f) * b, c.grid.origin.x, 1e-5f);
    EXPECT_NEAR((2.0f - 0.25f) * b, c.grid.origin.y, 1e-5f);
    EXPECT_NEAR((3.0f - 0.25f) * b, c.grid.origin.z, 1e-5f);
    ASSERT_EQ(2u, c.num_atoms);
    EXPECT_EQ(8, c.atomic_numbers[0]);
    EXPECT_NEAR(1.0f * b, c.xyz[0].x, 1e-5f);
    EXPECT_NEAR(1.0f * b, c.xyz[1].y, 1e-5f);
    // [i][j][k], k fastest
    EXPECT_EQ(23.0f, c.values[(1 * 3 + 2) * 4 + 3]);
    EXPECT_EQ(5.0f, c.values[(0 * 3 + 1) * 4 + 1]);
    EXPECT_NEAR(0.0334, c.bulk, 1e-9);

    // Too few values; not a cube
    std::string cut = text.substr(0, text.size() - 8);
    EXPECT_FALSE(cube_parse(&c, {cut.data(), cut.size()}, alloc, err, sizeof(err)));
    EXPECT_TRUE(strstr(err, "values") != nullptr);
    EXPECT_FALSE(cube_parse(&c, STR_LIT("hello\nworld\n"), alloc, err, sizeof(err)));
    md_arena_allocator_destroy(alloc);
}

namespace {
struct Peptide {
    md_allocator_i* alloc = nullptr;
    md_system_t sys = {};
    md_system_state_t state = {};
};
const Peptide* peptide() {
    static Peptide p;
    static bool ok = false;
    if (!p.alloc) {
        p.alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
        p.sys = {.alloc = p.alloc};
        p.state = {.alloc = p.alloc};
        ok = md_pdb_system_init_from_file(&p.sys, &p.state, STR_LIT(VIAMD_TEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) &&
             md_util_system_infer(&p.sys, &p.state, MD_UTIL_INFER_ALL);
    }
    return ok ? &p : nullptr;
}

mat3_t rotation(float ax, float ay, float az) {
    auto m = [](vec3_t c0, vec3_t c1, vec3_t c2) { mat3_t M; M.col[0] = c0; M.col[1] = c1; M.col[2] = c2; return M; };
    const mat3_t X = m({1, 0, 0}, {0, cosf(ax), sinf(ax)}, {0, -sinf(ax), cosf(ax)});
    const mat3_t Y = m({cosf(ay), 0, -sinf(ay)}, {0, 1, 0}, {sinf(ay), 0, cosf(ay)});
    const mat3_t Z = m({cosf(az), sinf(az), 0}, {-sinf(az), cosf(az), 0}, {0, 0, 1});
    return mat3_mul(Z, mat3_mul(Y, X));
}
}  // namespace

// The other structure is this one turned and moved, with an atom less (a mutation of sorts): placed back exactly
UTEST(density_volume_core, align) {
    const Peptide* p = peptide();
    ASSERT_TRUE(p != nullptr);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    // The first three residues
    std::vector<int32_t> atoms;
    for (size_t r = 0; r < 3; ++r) {
        const md_urange_t ar = md_system_residue_atom_range(&p->sys, r);
        for (uint32_t i = ar.beg; i < ar.end; ++i) atoms.push_back((int32_t)i);
    }
    ASSERT_GT(atoms.size(), 10u);
    const vec3_t* xyz = (const vec3_t*)p->state.xyz;
    const mat3_t R0 = rotation(0.3f, -1.1f, 2.0f);
    const vec3_t t0 = {5.0f, -3.0f, 12.0f};
    // The other: R0 x + t0, without the last heavy atom
    size_t drop = atoms.size();
    for (size_t i = atoms.size(); i-- > 0;) if (md_atom_atomic_number(&p->sys.atom, atoms[i]) > 1) { drop = i; break; }
    std::vector<int> z;
    std::vector<vec3_t> ox;
    std::vector<int32_t> kept;
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (i == drop) continue;
        z.push_back(md_atom_atomic_number(&p->sys.atom, atoms[i]));
        ox.push_back(vec3_add(mat3_mul_vec3(R0, xyz[atoms[i]]), t0));
        kept.push_back(atoms[i]);
    }
    Cube other;
    other.num_atoms = z.size();
    other.atomic_numbers = z.data();
    other.xyz = ox.data();
    Alignment al;
    char err[128] = "";
    ASSERT_TRUE(align(&al, &p->sys, atoms.data(), atoms.size(), xyz, other, alloc, err, sizeof(err)));
    EXPECT_FALSE(al.incomplete);
    EXPECT_EQ(al.this_heavy - 1, al.matched);
    EXPECT_EQ(al.other_heavy, al.matched);
    EXPECT_LT(al.rmsd, 1e-3);
    // Every atom of the other lands on its own
    for (size_t i = 0; i < kept.size(); ++i) {
        const vec3_t back = vec3_add(mat3_mul_vec3(al.R, ox[i]), al.t);
        EXPECT_LT(vec3_length(vec3_sub(back, xyz[kept[i]])), 2e-3f);
    }
    // No atoms, or a structure with nothing in common: nothing placed
    Cube none;
    EXPECT_FALSE(align(&al, &p->sys, atoms.data(), atoms.size(), xyz, none, alloc, err, sizeof(err)));
    int lone_z[2] = {26, 26};
    vec3_t lone_x[2] = {{0, 0, 0}, {10, 0, 0}};
    Cube lone;
    lone.num_atoms = 2;
    lone.atomic_numbers = lone_z;
    lone.xyz = lone_x;
    EXPECT_FALSE(align(&al, &p->sys, atoms.data(), atoms.size(), xyz, lone, alloc, err, sizeof(err)));
    EXPECT_TRUE(strstr(err, "common") != nullptr);
    md_arena_allocator_destroy(alloc);
}

// A field linear in space is sampled exactly by trilinear interpolation: through a turn and a move, from a grid of
// other steps, the values are the field's where the other grid reaches and 0 where it does not
UTEST(density_volume_core, resample) {
    const vec3_t g = {0.7f, -0.2f, 0.4f};
    const float c0 = 3.0f;
    auto field = [&](vec3_t x) { return vec3_dot(g, x) + c0; };

    Grid from;
    from.dim[0] = 20; from.dim[1] = 24; from.dim[2] = 28;
    from.origin = {-5.0f, -6.0f, -7.0f};
    from.basis[0] = {0.5f, 0.0f, 0.0f};
    from.basis[1] = {0.05f, 0.5f, 0.0f};
    from.basis[2] = {0.0f, 0.0f, 0.5f};
    std::vector<float> fv((size_t)from.dim[0] * from.dim[1] * from.dim[2]);
    for (uint32_t i = 0; i < from.dim[0]; ++i) for (uint32_t j = 0; j < from.dim[1]; ++j) for (uint32_t k = 0; k < from.dim[2]; ++k) {
        const vec3_t x = vec3_add(from.origin, vec3_add(vec3_mul1(from.basis[0], i + 0.5f), vec3_add(vec3_mul1(from.basis[1], j + 0.5f), vec3_mul1(from.basis[2], k + 0.5f))));
        fv[((size_t)i * from.dim[1] + j) * from.dim[2] + k] = field(x);
    }
    from.values = fv.data();

    Alignment al;
    al.R = rotation(0.2f, 0.5f, -0.4f);
    al.t = {1.0f, -0.5f, 2.0f};

    Grid to;
    to.dim[0] = 16; to.dim[1] = 16; to.dim[2] = 16;
    to.origin = {-4.0f, -3.0f, -2.0f};
    to.basis[0] = {0.6f, 0.0f, 0.0f};
    to.basis[1] = {0.0f, 0.6f, 0.0f};
    to.basis[2] = {0.0f, 0.0f, 0.6f};
    std::vector<float> out((size_t)16 * 16 * 16);
    const double covered = resample(out.data(), to, from, al);
    EXPECT_GT(covered, 0.2);
    EXPECT_LT(covered, 1.0);
    // Where from reaches, the field at R^T (p - t); elsewhere 0
    mat3_t B = {};
    for (int a = 0; a < 3; ++a) B.col[a] = from.basis[a];
    const mat3_t Binv = mat3_inverse(B);
    size_t inside = 0, checked = 0;
    for (uint32_t i = 0; i < 16; ++i) for (uint32_t j = 0; j < 16; ++j) for (uint32_t k = 0; k < 16; ++k) {
        const vec3_t p = vec3_add(to.origin, vec3_add(vec3_mul1(to.basis[0], i + 0.5f), vec3_add(vec3_mul1(to.basis[1], j + 0.5f), vec3_mul1(to.basis[2], k + 0.5f))));
        const vec3_t q = mat3_mul_vec3(mat3_transpose(al.R), vec3_sub(p, al.t));
        const vec3_t f = vec3_sub(mat3_mul_vec3(Binv, vec3_sub(q, from.origin)), vec3_set1(0.5f));
        const float v = out[((size_t)i * 16 + j) * 16 + k];
        // Within the outer centres, exact; past them by up to half a voxel, the outer voxels' values; beyond, 0
        const bool in = f.x >= 0 && f.y >= 0 && f.z >= 0 && f.x <= from.dim[0] - 1 && f.y <= from.dim[1] - 1 && f.z <= from.dim[2] - 1;
        const bool out_ = f.x < -0.5f || f.y < -0.5f || f.z < -0.5f || f.x > from.dim[0] - 0.5f || f.y > from.dim[1] - 0.5f || f.z > from.dim[2] - 0.5f;
        if (in) { EXPECT_NEAR(field(q), v, 1e-3f); inside += 1; }
        if (out_) EXPECT_EQ(0.0f, v);
        checked += 1;
    }
    EXPECT_GT(inside, 500u);
}
