#include <mol/aminoacid.h>
#include <mol/molecule_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>

#include <stdio.h>
#include <stdlib.h>

#define TEST(x) \
    { printf("TEST " x "\n"); }

#define EXPECT_EQ(x, y)          \
    if (x == y)                  \
        printf("    SUCCESS\n"); \
    else                         \
        printf("    FAIL\n");

int main() {
    TEST("com") {
        vec3 pos[8] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
        float mass[8] = {1, 1, 1, 1, 1, 1, 1, 1};
        Element elem[8] = {Element::C, Element::C, Element::C, Element::C, Element::C, Element::C, Element::C, Element::C};

        EXPECT_EQ(compute_com(pos, mass), vec3(0.5));
        EXPECT_EQ(compute_com(pos, elem), vec3(0.5));

        translate_positions(pos, vec3(-1));
        vec3 com = compute_periodic_com(pos, elem, vec3(4, 4, 4));
        printf("com periodic: %f %f %f\n", com.x, com.y, com.z);
    }

    // allocate_and_load_gro_from_file(&md.molecule, PROJECT_SOURCE_DIR "");

    return 0;
}
