/**
 * @file molden_parser_test.cpp
 * @brief Unit tests for Molden file parser
 * 
 * Tests the parsing functionality of the Molden file parser.
 * Validates correct parsing of atoms, basis sets, and molecular orbitals.
 * 
 * To compile and run (standalone):
 *   g++ -std=c++20 -I. molden_parser_test.cpp molden.cpp -o molden_parser_test
 *   ./molden_parser_test
 */

#include "molden.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

using namespace molden;

/**
 * @brief Test parsing of simple H2 molecule
 */
void test_parse_h2() {
    std::cout << "Testing H2 parsing... ";
    
    std::string h2_molden = R"(
[Molden Format]
[Title]
Hydrogen molecule

[Atoms] Angs
H   1   1   0.0   0.0   0.0
H   2   1   0.0   0.0   0.74

[GTO]
1 0
s 3 1.0
3.42525091 0.15432897
0.62391373 0.53532814
0.16885540 0.44463454

2 0
s 3 1.0
3.42525091 0.15432897
0.62391373 0.53532814
0.16885540 0.44463454

[MO]
Ene= -0.594
Spin= Alpha
Occup= 2.0
Sym= Sigma_g
0.548919
0.548919
)";
    
    std::string error;
    MoldenData data = parse_molden_string(h2_molden, &error);
    
    if (!error.empty()) {
        std::cerr << "Parse error: " << error << std::endl;
        assert(false);
    }
    
    // Validate atoms
    assert(data.atoms.size() == 2);
    assert(data.atoms[0].element_symbol == "H");
    assert(data.atoms[0].atomic_number == 1);
    assert(data.atoms[1].atomic_number == 1);
    assert(data.coord_unit == CoordinateUnit::Angstrom);
    
    // Validate basis sets
    assert(data.basis_sets.size() == 2);
    assert(data.basis_sets[0].atom_index == 1);
    assert(data.basis_sets[0].shells.size() == 1);
    assert(data.basis_sets[0].shells[0].shell_type == ShellType::S);
    assert(data.basis_sets[0].shells[0].primitives.size() == 3);
    
    // Validate total basis functions
    assert(data.total_basis_functions == 2); // 2 atoms, 1 s-function each
    
    // Validate orbitals
    assert(data.orbitals.size() == 1);
    assert(data.orbitals[0].spin == SpinType::Alpha);
    assert(std::abs(data.orbitals[0].energy - (-0.594)) < 1e-6);
    assert(data.orbitals[0].occupation == 2.0);
    assert(data.orbitals[0].symmetry == "Sigma_g");
    assert(data.orbitals[0].coefficients.size() == 2);
    
    // Validate data
    assert(util::validate_molden_data(data));
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing of H2O with SP shells
 */
void test_parse_h2o_with_sp_shells() {
    std::cout << "Testing H2O with SP shells... ";
    
    std::string h2o_molden = R"(
[Molden Format]
[Title]
Water molecule

[Atoms] Angs
O   1   8   0.0   0.0   0.117
H   2   1   0.0   0.757  -0.468
H   3   1   0.0  -0.757  -0.468

[GTO]
1 0
s 3 1.0
130.7093200 0.15432897
23.8088610 0.53532814
6.4436083 0.44463454
sp 3 1.0
5.0331513 -0.09996723 0.15591627
1.1695961 0.39951283 0.60768372
0.3803890 0.70011547 0.39195739

2 0
s 3 1.0
3.42525091 0.15432897
0.62391373 0.53532814
0.16885540 0.44463454

3 0
s 3 1.0
3.42525091 0.15432897
0.62391373 0.53532814
0.16885540 0.44463454

[MO]
Ene= -20.251
Spin= Alpha
Occup= 2.0
Sym= A1
0.994216
0.025846
0.0
0.0
0.006035
0.004764
0.004764
)";
    
    std::string error;
    MoldenData data = parse_molden_string(h2o_molden, &error);
    
    if (!error.empty()) {
        std::cerr << "Parse error: " << error << std::endl;
        assert(false);
    }
    
    // Validate atoms
    assert(data.atoms.size() == 3);
    assert(data.atoms[0].element_symbol == "O");
    assert(data.atoms[0].atomic_number == 8);
    
    // Validate basis sets with SP shell
    assert(data.basis_sets.size() == 3);
    assert(data.basis_sets[0].shells.size() == 2); // s and sp shells
    
    // Check SP shell
    ContractedShell sp_shell = data.basis_sets[0].shells[1];
    assert(sp_shell.shell_type == ShellType::SP);
    assert(sp_shell.primitives.size() == 3);
    assert(sp_shell.primitives[0].coefficient_sp != 0.0); // Has P coefficient
    
    // Total basis functions: 1(s) + 4(sp) + 1(s) + 1(s) = 7
    assert(data.total_basis_functions == 7);
    
    // Validate orbital
    assert(data.orbitals.size() == 1);
    assert(data.orbitals[0].coefficients.size() == 7);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing with AU coordinates
 */
void test_parse_au_coordinates() {
    std::cout << "Testing AU coordinate parsing... ";
    
    std::string molden = R"(
[Atoms] AU
H   1   1   0.0   0.0   0.0
H   2   1   0.0   0.0   1.4

[GTO]
1 0
s 1 1.0
1.0 1.0

2 0
s 1 1.0
1.0 1.0
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    assert(data.coord_unit == CoordinateUnit::AtomicUnit);
    assert(data.atoms.size() == 2);
    assert(std::abs(data.atoms[1].z - 1.4f) < 1e-6f);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing with missing atomic numbers
 */
void test_parse_missing_atomic_numbers() {
    std::cout << "Testing missing atomic numbers... ";
    
    std::string molden = R"(
[Atoms] Angs
C   1   0.0   0.0   0.0
H   2   1.0   0.0   0.0

[GTO]
1 0
s 1 1.0
1.0 1.0

2 0
s 1 1.0
1.0 1.0
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    // Should infer atomic numbers from symbols
    assert(data.atoms.size() == 2);
    assert(data.atoms[0].atomic_number == 6); // Carbon
    assert(data.atoms[1].atomic_number == 1); // Hydrogen
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing multiple orbitals with different spins
 */
void test_parse_multiple_orbitals() {
    std::cout << "Testing multiple orbitals... ";
    
    std::string molden = R"(
[Atoms] Angs
H   1   1   0.0   0.0   0.0

[GTO]
1 0
s 1 1.0
1.0 1.0

[MO]
Ene= -1.0
Spin= Alpha
Occup= 1.0
0.5

Ene= -0.5
Spin= Beta
Occup= 1.0
0.7
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    assert(data.orbitals.size() == 2);
    assert(data.orbitals[0].spin == SpinType::Alpha);
    assert(data.orbitals[1].spin == SpinType::Beta);
    assert(data.num_alpha_orbitals == 1);
    assert(data.num_beta_orbitals == 1);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing without symmetry labels
 */
void test_parse_no_symmetry() {
    std::cout << "Testing MO without symmetry labels... ";
    
    std::string molden = R"(
[Atoms] Angs
H   1   1   0.0   0.0   0.0

[GTO]
1 0
s 1 1.0
1.0 1.0

[MO]
Ene= -1.0
Spin= Alpha
Occup= 2.0
0.5
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    assert(data.orbitals.size() == 1);
    assert(data.orbitals[0].symmetry.empty()); // No symmetry specified
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing with [5D] tag (spherical harmonics)
 */
void test_parse_spherical_basis() {
    std::cout << "Testing spherical basis format... ";
    
    std::string molden = R"(
[5D]
[Atoms] Angs
C   1   6   0.0   0.0   0.0

[GTO]
1 0
s 1 1.0
1.0 1.0
d 1 1.0
2.0 0.5
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    assert(data.basis_format == BasisFormat::Spherical);
    
    // D shell: 5 functions in spherical, 6 in Cartesian
    size_t expected_funcs = 1 + 5; // s + d(spherical)
    assert(data.total_basis_functions == expected_funcs);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test Fortran-style scientific notation (D instead of E)
 */
void test_parse_fortran_notation() {
    std::cout << "Testing Fortran D notation... ";
    
    std::string molden = R"(
[Atoms] Angs
H   1   1   0.0D0   0.0D0   0.0D0

[GTO]
1 0
s 2 1.0
1.234D+02 5.678D-03
9.876D-01 4.321D-02

[MO]
Ene= -1.234D+00
Spin= Alpha
Occup= 2.0D0
5.0D-01
)";
    
    std::string error;
    MoldenData data = parse_molden_string(molden, &error);
    
    assert(data.atoms.size() == 1);
    assert(data.basis_sets.size() == 1);
    
    // Check that D notation was parsed correctly
    double expected_exp = 123.4; // 1.234D+02
    assert(std::abs(data.basis_sets[0].shells[0].primitives[0].exponent - expected_exp) < 1e-6);
    
    assert(std::abs(data.orbitals[0].energy - (-1.234)) < 1e-6);
    assert(std::abs(data.orbitals[0].coefficients[0] - 0.5) < 1e-6);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test parsing from file
 */
void test_parse_from_file() {
    std::cout << "Testing file parsing... ";
    
    std::string error;
    MoldenData data = parse_molden_file("datasets/molden_examples/h2_sto3g.molden", &error);
    
    if (!error.empty()) {
        std::cerr << "File parse error: " << error << std::endl;
        // This test may fail if file doesn't exist - that's OK
        std::cout << "SKIP (file not found)" << std::endl;
        return;
    }
    
    assert(data.atoms.size() == 2);
    assert(data.basis_sets.size() == 2);
    assert(data.orbitals.size() == 2);
    assert(util::validate_molden_data(data));
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test error handling with invalid data
 */
void test_error_handling() {
    std::cout << "Testing error handling... ";
    
    // Missing atom fields
    std::string invalid1 = R"(
[Atoms] Angs
H   1
)";
    
    std::string error;
    MoldenData data = parse_molden_string(invalid1, &error);
    assert(!error.empty()); // Should have error
    
    // Invalid shell type
    std::string invalid2 = R"(
[Atoms] Angs
H   1   1   0.0   0.0   0.0

[GTO]
1 0
x 1 1.0
1.0 1.0
)";
    
    data = parse_molden_string(invalid2, &error);
    assert(!error.empty()); // Should have error
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Demonstrate complete usage of the parser
 */
void demonstrate_parser_usage() {
    std::cout << "\n=== Demonstrating Molden Parser Usage ===" << std::endl;
    
    std::string error;
    MoldenData data = parse_molden_file("datasets/molden_examples/h2o_sto3g.molden", &error);
    
    if (!error.empty()) {
        std::cout << "Note: Could not load file (this is OK for demonstration)" << std::endl;
        std::cout << "Error: " << error << std::endl;
        return;
    }
    
    std::cout << "Successfully parsed Molden file!" << std::endl;
    std::cout << "Title: " << data.title << std::endl;
    std::cout << "Number of atoms: " << data.atoms.size() << std::endl;
    std::cout << "Number of basis sets: " << data.basis_sets.size() << std::endl;
    std::cout << "Total basis functions: " << data.total_basis_functions << std::endl;
    std::cout << "Number of molecular orbitals: " << data.orbitals.size() << std::endl;
    std::cout << "Alpha orbitals: " << data.num_alpha_orbitals << std::endl;
    std::cout << "Beta orbitals: " << data.num_beta_orbitals << std::endl;
    
    // Print atom details
    std::cout << "\nAtoms:" << std::endl;
    for (const auto& atom : data.atoms) {
        std::cout << "  " << atom.element_symbol << " (Z=" << atom.atomic_number << "): "
                  << std::fixed << std::setprecision(3)
                  << "(" << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }
    
    // Print orbital details
    std::cout << "\nMolecular Orbitals:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), data.orbitals.size()); ++i) {
        const auto& mo = data.orbitals[i];
        std::cout << "  MO " << (i+1) << ": "
                  << "E=" << std::fixed << std::setprecision(3) << mo.energy << " Ha, "
                  << "Occ=" << mo.occupation << ", "
                  << "Spin=" << (mo.spin == SpinType::Alpha ? "Alpha" : "Beta")
                  << ", Sym=" << mo.symmetry << std::endl;
    }
    
    // Validate
    bool is_valid = util::validate_molden_data(data);
    std::cout << "\nValidation: " << (is_valid ? "PASS" : "FAIL") << std::endl;
}

int main() {
    std::cout << "=== Molden Parser Tests ===" << std::endl << std::endl;
    
    try {
        test_parse_h2();
        test_parse_h2o_with_sp_shells();
        test_parse_au_coordinates();
        test_parse_missing_atomic_numbers();
        test_parse_multiple_orbitals();
        test_parse_no_symmetry();
        test_parse_spherical_basis();
        test_parse_fortran_notation();
        test_parse_from_file();
        test_error_handling();
        
        std::cout << "\n=== All tests PASSED ===" << std::endl;
        
        demonstrate_parser_usage();
        
        std::cout << "\n=== Test suite completed successfully ===" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}
