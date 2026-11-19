#include "utest.h"

#include <grid_field.h>
#include <loader.h>

#include <md_cube.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_file.h>
#include <core/md_str.h>

#include <math.h>
#include <stdio.h>
#include <string>

// A grid of a cube file, placed on the atoms it is of as they are: values linear in the index (exact under trilinear
// interpolation), three atoms to place it by. Bohr in the file.
static std::string cube_text() {
    std::string t = "grid field\nlinear values\n";
    t += "    3    0.000000    0.000000    0.000000\n";
    t += "    4    0.500000    0.000000    0.000000\n";
    t += "    5    0.000000    0.500000    0.000000\n";
    t += "    6    0.000000    0.000000    0.500000\n";
    t += "    8    8.000000    1.000000    1.000000    1.000000\n";
    t += "    1    1.000000    1.800000    1.000000    1.000000\n";
    t += "    1    1.000000    1.000000    1.900000    1.200000\n";
    char buf[64];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 5; ++j)
            for (int k = 0; k < 6; ++k) {
                snprintf(buf, sizeof(buf), " %g", (double)(i + 10 * j + 100 * k));
                t += buf;
                if (k == 5) t += "\n";
            }
    return t;
}

struct GridFieldTest {
    md_allocator_i* alloc = nullptr;
    md_system_t sys = {};
    md_system_state_t state = {};
    md_attribute_id_t key = MD_ATTRIBUTE_INVALID;
    char path[64] = "grid_field_test.cube";
};

static bool grid_field_test_init(GridFieldTest* t) {
    const std::string text = cube_text();
    md_file_t f = {};
    if (!md_file_open(&f, md_str_from_cstr(t->path), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    md_file_write(f, text.data(), text.size());
    md_file_close(&f);
    t->alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    t->sys.alloc = t->alloc;
    t->state.alloc = t->alloc;
    if (!md_cube_system_init_from_file(&t->sys, &t->state, md_str_from_cstr(t->path), nullptr)) return false;
    t->key = md_attributes_id_from_path(MD_STR_LIT("cube/grid_field_test"));
    return true;
}

static void grid_field_test_free(GridFieldTest* t) {
    md_arena_allocator_destroy(t->alloc);
    remove(t->path);
}

// The centre of voxel (i, j, k) in bohr, as the file has it
static void centre_bohr(float out[3], float i, float j, float k) {
    out[0] = 0.5f * i;
    out[1] = 0.5f * j;
    out[2] = 0.5f * k;
}

UTEST(grid_field, as_the_file_has_it) {
    GridFieldTest t;
    ASSERT_TRUE(grid_field_test_init(&t));
    EXPECT_EQ(1u, grid_field_count(t.sys, t.key));
    EXPECT_EQ(t.key, grid_field_first(t.sys));
    GridField g;
    ASSERT_TRUE(grid_field_get(&g, t.sys, t.state, t.key, 0));
    EXPECT_EQ(3u, g.num_placed_by);
    EXPECT_NEAR(0.0, g.rmsd, 1e-5);
    float p[3];
    centre_bohr(p, 2, 3, 4);
    EXPECT_NEAR(2 + 30 + 400, grid_field_sample_bohr(g, p), 1e-3);
    // between the points: linear, so exact
    centre_bohr(p, 1.25f, 2.5f, 3.75f);
    EXPECT_NEAR(1.25 + 25 + 375, grid_field_sample_bohr(g, p), 1e-2);
    // outside: nothing
    centre_bohr(p, -2, 0, 0);
    EXPECT_EQ(0.0f, grid_field_sample_bohr(g, p));

    // its own grid as an md_grid_t: voxel centres at origin + (i + 1/2) spacing
    md_grid_t grid;
    ASSERT_TRUE(grid_field_md_grid(&grid, g));
    EXPECT_EQ(4, grid.dim[0]);
    EXPECT_EQ(6, grid.dim[2]);
    EXPECT_NEAR(-0.25, grid.origin.x, 1e-4);
    EXPECT_NEAR(0.5, grid.spacing.y, 1e-4);
    grid_field_free(&g);
    grid_field_test_free(&t);
}

// The molecule turned 90 degrees about z and moved: the grid goes with it
UTEST(grid_field, follows_its_atoms) {
    GridFieldTest t;
    ASSERT_TRUE(grid_field_test_init(&t));
    const float bohr = 0.52917721092f;
    const md_float3_t shift = {3.0f, -1.0f, 2.0f};    // A
    auto move = [&](md_float3_t a) { return md_float3_t{-a.y + shift.x, a.x + shift.y, a.z + shift.z}; };
    for (size_t i = 0; i < t.state.num_atoms; ++i) t.state.xyz[i] = move(t.state.xyz[i]);

    GridField g;
    ASSERT_TRUE(grid_field_get(&g, t.sys, t.state, t.key, 0));
    EXPECT_NEAR(0.0, g.rmsd, 1e-4);
    float c[3];
    centre_bohr(c, 2, 1, 3);
    const md_float3_t w = move(md_float3_t{c[0] * bohr, c[1] * bohr, c[2] * bohr});
    const float p[3] = {w.x / bohr, w.y / bohr, w.z / bohr};
    EXPECT_NEAR(2 + 10 + 300, grid_field_sample_bohr(g, p), 1e-2);
    // where it was is outside now
    EXPECT_NE(2.0f + 10.0f + 300.0f, grid_field_sample_bohr(g, c));
    grid_field_free(&g);
    grid_field_test_free(&t);
}

// A cube of the loaded atoms supplements them rather than replacing the system
UTEST(grid_field, a_cube_of_the_atoms_supplements) {
    GridFieldTest t;
    ASSERT_TRUE(grid_field_test_init(&t));
    loader::LoaderState s;
    loader::init(&s, md_str_from_cstr(t.path), &t.sys);
    EXPECT_EQ((LoaderType)LoaderType_CUBE, s.type);
    EXPECT_NE(0u, s.flags & LoaderFlag_Supplemental);
    loader::init(&s, md_str_from_cstr(t.path), nullptr);
    EXPECT_EQ(0u, s.flags & LoaderFlag_Supplemental);
    EXPECT_NE(0u, s.flags & LoaderFlag_System);
    grid_field_test_free(&t);
}
