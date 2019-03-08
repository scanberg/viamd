#include <core/gl.h>
#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/string_utils.h>
#include <core/volume.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/spatial_hash.h>

#include <glm/gtx/io.hpp>

#ifdef __clang__
#define CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS
#endif
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

constexpr const char* CAFFINE_PDB = R"(
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

#include <chrono>
#define TIME() std::chrono::high_resolution_clock::now()
#define MILLISEC(x, y) std::chrono::duration_cast<std::chrono::milliseconds>(y - x).count()

TEST_CASE("Testing PdbInfo", "[PdbInfo]") {
    const auto t0 = TIME();
    // String pdb_str = allocate_and_read_textfile(VIAMD_DATA_DIR "/alanine/14ns-300K.pdb");
    String pdb_str = allocate_and_read_textfile(VIAMD_DATA_DIR "/haofan/for_VIAMD.pdb");
    defer { free_string(&pdb_str); };
    const auto t1 = TIME();

    const auto t2 = TIME();
    PdbInfo info;
    extract_pdb_info(&info, pdb_str);
    const auto t3 = TIME();

    const auto t4 = TIME();
    /*
    MoleculeDynamic md;
    allocate_and_parse_pdb_from_string(&md, pdb_str);
    defer {
        free_molecule_structure(&md.molecule);
        free_trajectory(&md.trajectory);
    };
    */
    const auto t5 = TIME();

    printf("Time to load dataset: %.2f\n", (double)MILLISEC(t0, t1));
    printf("Time to extract pdb info: %.2f\n", (double)MILLISEC(t2, t3));
    printf("Time to parse full pdb: %.2f\n", (double)MILLISEC(t4, t5));
    printf("PdbInfo:\n num_atoms: %i \n num_residues: %i \n num_chains: %i \n num_frames: %i \n", info.num_atoms, info.num_residues, info.num_chains, info.num_frames);
}

TEST_CASE("Testing DynamicArray", "[DynamicArray]") {
    // MoleculeDynamic md;
    // allocate_and_parse_pdb_from_string(&md, CAFFINE_PDB);
    // defer { free_molecule_structure(&md.molecule); };
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
}

TEST_CASE("Testing pdb loader caffine", "[parse_pdb]") {
    MoleculeDynamic md;
    allocate_and_parse_pdb_from_string(&md, CAFFINE_PDB);
    defer { free_molecule_structure(&md.molecule); };

    REQUIRE(md.molecule.atom.count == 24);
}

TEST_CASE("Testing filter", "[filter]") {
    MoleculeDynamic md;
    allocate_and_parse_pdb_from_string(&md, CAFFINE_PDB);
    defer { free_molecule_structure(&md.molecule); };

    filter::initialize();
    DynamicArray<bool> mask(md.molecule.atom.count);

    SECTION("filter element N") {
        filter::compute_filter_mask(mask, md, "element N");
        for (int32 i = 0; i < mask.size(); i++) {
            if (i == 0 || i == 4 || i == 16 || i == 17) {
                REQUIRE(mask[i] == true);
            } else {
                REQUIRE(mask[i] == false);
            }
        }
    }

    SECTION("filter atom 1:10") {
        filter::compute_filter_mask(mask, md, "atom 1:10");
        for (int32 i = 0; i < mask.size(); i++) {
            if (0 <= i && i <= 9) {
                REQUIRE(mask[i] == true);
            } else {
                REQUIRE(mask[i] == false);
            }
        }
    }

    SECTION("filter atom 10:*") {
        filter::compute_filter_mask(mask, md, "atom 10:*");
        for (int32 i = 0; i < mask.size(); i++) {
            if (0 <= i && i < 9) {
                REQUIRE(mask[i] == true);
            } else {
                REQUIRE(mask[i] == false);
            }
        }
    }

    SECTION("filter atom *:*") {
        filter::compute_filter_mask(mask, md, "atom *:*");
        for (int32 i = 0; i < mask.size(); i++) {
            REQUIRE(mask[i] == true);
        }
    }

    SECTION("filter all") {
        filter::compute_filter_mask(mask, md, "all");
        for (int32 i = 0; i < mask.size(); i++) {
            REQUIRE(mask[i] == true);
        }
    }

    SECTION("filter not all") {
        filter::compute_filter_mask(mask, md, "not all");
        for (int32 i = 0; i < mask.size(); i++) {
            REQUIRE(mask[i] == false);
        }
    }
}
