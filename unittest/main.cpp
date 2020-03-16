#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/string_utils.h>
#include <core/spatial_hash.h>
#include <core/bitfield.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>

//#include <glm/gtx/io.hpp>

#ifdef __clang__
#define CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS
#endif
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

extern void transform_ref(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp = 1);

constexpr CStringView CAFFINE_PDB = R"(
ATOM      1  N1  BENZ    1       5.040   1.944  -8.324                          
ATOM      2  C2  BENZ    1       6.469   2.092  -7.915                          
ATOM      3  C3  BENZ    1       7.431   0.865  -8.072                          
ATOM      4  C4  BENZ    1       6.916  -0.391  -8.544                          
ATOM      5  N5  BENZ    1       5.532  -0.541  -8.901                          
ATOM      6  C6  BENZ    1       4.590   0.523  -8.394                          
ATOM      7  C11 BENZ    1       4.045   3.041  -8.005                          
ATOM      8  H111BENZ    1       4.453   4.038  -8.264                          
ATOM      9  H112BENZ    1       3.101   2.907  -8.570                          
ATOM     10  H113BENZ    1       3.795   3.050  -6.926                          
ATOM     11  O21 BENZ    1       6.879   3.181  -7.503                          
ATOM     12  C51 BENZ    1       4.907  -1.659  -9.696                          
ATOM     13  H511BENZ    1       4.397  -1.273 -10.599                          
ATOM     14  H512BENZ    1       5.669  -2.391 -10.028                          
ATOM     15  H513BENZ    1       4.161  -2.209  -9.089                          
ATOM     16  O61 BENZ    1       3.470   0.208  -7.986                          
ATOM     17  N1  NSP3    1B      8.807   0.809  -7.799                          
ATOM     18  N1  NSP3    1C      7.982  -1.285  -8.604                          
ATOM     19  C1  CSP3    1D      9.015  -0.500  -8.152                          
ATOM     20  H1  CSP3    1D     10.007  -0.926  -8.079                          
ATOM     21  C1  CSP3    1E      9.756   1.835  -7.299                          
ATOM     22  H11 CSP3    1E     10.776   1.419  -7.199                          
ATOM     23  H12 CSP3    1E      9.437   2.207  -6.309                          
ATOM     24  H13 CSP3    1E      9.801   2.693  -7.994
)";

TEST_CASE("String utils", "[String Utils]") {
    REQUIRE(to_upper('a') == 'A');
    REQUIRE(to_upper('b') == 'B');
    REQUIRE(to_upper('z') == 'Z');
    REQUIRE(to_upper('6') == '6');

    REQUIRE(to_lower('A') == 'a');
    REQUIRE(to_lower('B') == 'b');
    REQUIRE(to_lower('Z') == 'z');
    REQUIRE(to_lower('6') == '6');
}

TEST_CASE("Bitfield", "[Bitfield]") {
    constexpr int size = 257;
    Bitfield field;
    bitfield::init(&field, size);
    defer { bitfield::free(&field); };
    bitfield::set_bit(field, 2);
    bitfield::set_bit(field, 3);
    bitfield::set_bit(field, 5);
    bitfield::set_bit(field, 255);

    REQUIRE(field.size() == size);
    REQUIRE(field.size_in_bytes() == size / 8 + (size % 8 ? 1 : 0));

    REQUIRE(bitfield::get_bit(field, 2) == true);
    REQUIRE(bitfield::get_bit(field, 3) == true);
    REQUIRE(bitfield::get_bit(field, 5) == true);
    REQUIRE(bitfield::get_bit(field, 255) == true);
    REQUIRE(bitfield::number_of_bits_set(field) == 4);

    REQUIRE(field[2] == true);
    REQUIRE(field[3] == true);
    REQUIRE(field[5] == true);
    REQUIRE(field[255] == true);

    REQUIRE(bitfield::find_first_bit_set(field) == 2);
    REQUIRE(bitfield::find_last_bit_set(field) == 255);

    // bitfield::print(field);
    // printf("\n");

    bitfield::clear_all(field);
    for (i64 i = 0; i < field.size(); i++) {
        REQUIRE(bitfield::get_bit(field, i) == false);
    }

    REQUIRE(bitfield::any_bit_set(field) == false);

    // --- RANGE ----
    const auto beg = 33;
    const auto end = 129;
    bitfield::set_range(field, Range<int>(beg, end));
    REQUIRE(bitfield::number_of_bits_set(field) == (end - beg));

    // bitfield::print(field);
    // printf("\n");

    for (i64 i = 0; i < field.size(); i++) {
        if (beg <= i && i < end) {
            REQUIRE(bitfield::get_bit(field, i) == true);
        } else {
            REQUIRE(bitfield::get_bit(field, i) == false);
        }
    }

    REQUIRE(bitfield::any_bit_set_in_range(field, Range<int>(0, beg)) == false);
    REQUIRE(bitfield::any_bit_set_in_range(field, Range<int>(beg, end)) == true);
    REQUIRE(bitfield::all_bits_set_in_range(field, Range<int>(beg, end)) == true);
    REQUIRE(bitfield::any_bit_set_in_range(field, Range<int>(end, size)) == false);

    bitfield::clear_all(field);
    REQUIRE(bitfield::any_bit_set(field) == false);
    REQUIRE(bitfield::find_next_bit_set(field) == -1);

    bitfield::set_bit(field, 211);
    REQUIRE(bitfield::find_next_bit_set(field) == 211);

    // bitfield::print(field);
    // printf("\n");
}

TEST_CASE("Array", "[Array]") {
    static constexpr int data[4] = {1, 2, 4, 5};
    constexpr Array<const int> arr(data);
    STATIC_ASSERT(arr.size() == 4, "Expected length to be 4");
    STATIC_ASSERT(arr[2] == 4, "Expected data[2] to be 4");

    constexpr CStringView cstr = "Cool";
    STATIC_ASSERT(cstr.size() == 4, "Expected length to be 4");
}

TEST_CASE("Testing DynamicArray", "[DynamicArray]") {
    DynamicArray<int> da1;
    DynamicArray<int> da2;

    auto func = []() -> DynamicArray<int> {
        DynamicArray<int> cool;
        cool.push_back(1);
        cool.push_back(2);
        cool.push_back(3);
        return cool;
    };

    da1 = func();

    REQUIRE(da1.size() == 3);
}

TEST_CASE("Molecule Utils", "[molecule_utils]") {
    MoleculeStructure mol;
    pdb::load_molecule_from_string(&mol, CAFFINE_PDB);
    defer { free_molecule_structure(&mol); };

    SECTION("COM: equal mass") {
        vec3 ref = {0, 0, 0};
        for (i64 i = 0; i < mol.atom.count; i++) {
            ref.x += mol.atom.position.x[i];
            ref.y += mol.atom.position.y[i];
            ref.z += mol.atom.position.z[i];
        }
        ref /= (float)mol.atom.count;

        const vec3 com = compute_com(mol.atom.position, mol.atom.count);

        REQUIRE(ref.x == Approx(com.x));
        REQUIRE(ref.y == Approx(com.y));
        REQUIRE(ref.z == Approx(com.z));
    }

    SECTION("COM: individual mass") {
        vec3 ref = {0, 0, 0};
        float sum = 0.0f;
        for (i64 i = 0; i < mol.atom.count; i++) {
            const auto m = mol.atom.mass[i];
            ref.x += mol.atom.position.x[i] * m;
            ref.y += mol.atom.position.y[i] * m;
            ref.z += mol.atom.position.z[i] * m;
            sum += m;
        }
        ref /= sum;

        const vec3 com = compute_com(mol.atom.position, mol.atom.mass, mol.atom.count);

        REQUIRE(ref.x == Approx(com.x));
        REQUIRE(ref.y == Approx(com.y));
        REQUIRE(ref.z == Approx(com.z));
    }

    SECTION("COM Periodic") {
        MoleculeStructureDescriptor desc;
        desc.num_atoms = 10;
        desc.num_residues = 4;
        desc.num_chains = 1;
        
        MoleculeStructure molecule;
        init_molecule_structure(&molecule, desc);
        defer { free_molecule_structure(&molecule); };

        for (int i = 0; i < 10; i++) {
            if (i < 5) {
                molecule.atom.position.x[i] = 5.0f + i * 1.0f;
            } else {
                molecule.atom.position.x[i] = -5.0f + i * 1.0f;
            }

            molecule.atom.position.y[i] = 5.0f;
            molecule.atom.position.z[i] = 5.0f;
            molecule.atom.mass[i] = 1.0f;
        }
        molecule.residue.atom_range[0] = {0, 3};
        molecule.residue.atom_range[1] = {3, 5};
        molecule.residue.atom_range[2] = {5, 8};
        molecule.residue.atom_range[3] = {8, 10};

        molecule.chain.atom_range[0] = {0, 10};
        molecule.chain.residue_range[0] = {0, 4};

        const mat3 box = {10, 0, 0, 0, 10, 0, 0, 0, 10};
        const vec3 com = compute_com_periodic(molecule.atom.position, molecule.atom.mass, molecule.atom.count, box);
        printf("com: %.2f %.2f %.2f\n", com.x, com.y, com.z);
    }

    SECTION("TRANSFORM") {
        const mat4 matrix = math::mat4_cast(math::angle_axis(math::PI / 4.0f, math::normalize(vec3(1, 1, 1))));

        void* mem = TMP_MALLOC(mol.atom.count * sizeof(float) * 6);
        defer { TMP_FREE(mem); };
        soa_vec3 ref = {
            (float*)mem + 0 * mol.atom.count,
            (float*)mem + 1 * mol.atom.count,
            (float*)mem + 2 * mol.atom.count };
        soa_vec3 pos = {
            (float*)mem + 3 * mol.atom.count,
            (float*)mem + 4 * mol.atom.count,
            (float*)mem + 5 * mol.atom.count };

        memcpy(ref.x, mol.atom.position.x, mol.atom.count * sizeof(float));
        memcpy(ref.y, mol.atom.position.y, mol.atom.count * sizeof(float));
        memcpy(ref.z, mol.atom.position.z, mol.atom.count * sizeof(float));
        memcpy(pos.x, mol.atom.position.x, mol.atom.count * sizeof(float));
        memcpy(pos.y, mol.atom.position.y, mol.atom.count * sizeof(float));
        memcpy(pos.z, mol.atom.position.z, mol.atom.count * sizeof(float));

        transform_ref(ref, mol.atom.count, matrix);
        transform(pos, mol.atom.count, matrix);

        for (i64 i = 0; i < mol.atom.count; i++) {
            REQUIRE(ref.x[i] == Approx(pos.x[i]));
            REQUIRE(ref.y[i] == Approx(pos.y[i]));
            REQUIRE(ref.z[i] == Approx(pos.z[i]));
        }
    }
}

TEST_CASE("Testing pdb loader caffine", "[parse_pdb]") {
    MoleculeStructure mol;
    pdb::load_molecule_from_string(&mol, CAFFINE_PDB);
    defer { free_molecule_structure(&mol); };

    REQUIRE(mol.atom.count == 24);
}

TEST_CASE("Testing filter", "[filter]") {
    MoleculeStructure mol;
    pdb::load_molecule_from_string(&mol, CAFFINE_PDB);
    defer { free_molecule_structure(&mol); };

    filter::initialize();
    Bitfield mask;
    bitfield::init(&mask, mol.atom.count);

    SECTION("filter element N") {
        filter::compute_filter_mask(mask, "element N", mol);
        for (i32 i = 0; i < mask.size(); i++) {
            if (i == 0 || i == 4 || i == 16 || i == 17) {
                REQUIRE(bitfield::get_bit(mask, i) == true);
            } else {
                REQUIRE(bitfield::get_bit(mask, i) == false);
            }
        }
    }

    SECTION("filter atom 1:10") {
        filter::compute_filter_mask(mask, "atom 1:10", mol);
        for (i32 i = 0; i < mask.size(); i++) {
            if (0 <= i && i <= 9) {
                REQUIRE(bitfield::get_bit(mask, i) == true);
            } else {
                REQUIRE(bitfield::get_bit(mask, i) == false);
            }
        }
    }

    SECTION("filter atom 10:*") {
        filter::compute_filter_mask(mask, "atom 10:*", mol);
        for (i32 i = 0; i < mask.size(); i++) {
            if (0 <= i && i < 9) {
                REQUIRE(bitfield::get_bit(mask, i) == false);
            } else {
                REQUIRE(bitfield::get_bit(mask, i) == true);
            }
        }
    }

    SECTION("filter atom *:*") {
        filter::compute_filter_mask(mask, "atom *:*", mol);
        REQUIRE(bitfield::all_bits_set(mask) == true);
    }

    SECTION("filter atom *") {
        filter::compute_filter_mask(mask, "atom *", mol);
        REQUIRE(bitfield::all_bits_set(mask) == true);
    }

    SECTION("filter all") {
        filter::compute_filter_mask(mask, "all", mol);
        REQUIRE(bitfield::all_bits_set(mask) == true);
    }

    SECTION("filter not all") {
        filter::compute_filter_mask(mask, "not all", mol);
        REQUIRE(bitfield::any_bit_set(mask) == false);
    }

    SECTION("filter residue *") {
        filter::compute_filter_mask(mask, "residue *", mol);
        REQUIRE(bitfield::all_bits_set(mask) == true);
    }
}

TEST_CASE("Testing string_utils", "[string_utils]") {
    SECTION("find pattern") {
        auto matches = find_patterns_in_file(UNITTEST_DATA_DIR"/test.txt", "KALLE");
        REQUIRE(matches.size() == 4);
    }
}
