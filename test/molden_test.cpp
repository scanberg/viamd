/**
 * @file molden_test.cpp
 * @brief Simple test demonstrating Molden data structure usage
 * 
 * This file shows how to use the Molden data structures and validates
 * that the design is sound. It's not a comprehensive unit test suite,
 * but rather a proof-of-concept demonstrating the API.
 * 
 * To compile and run (standalone):
 *   g++ -std=c++20 -I. molden_test.cpp molden.cpp -o molden_test
 *   ./molden_test
 */

#include "molden.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace molden;

/**
 * @brief Test basic Atom structure
 */
void test_atom_structure() {
    std::cout << "Testing Atom structure... ";
    
    Atom atom;
    atom.element_symbol = "C";
    atom.atom_index = 1;
    atom.atomic_number = 6;
    atom.x = 0.0f;
    atom.y = 0.0f;
    atom.z = 0.0f;
    
    assert(atom.element_symbol == "C");
    assert(atom.atomic_number == 6);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test basis set structures
 */
void test_basis_structures() {
    std::cout << "Testing basis set structures... ";
    
    // Create a simple S shell with 3 primitives
    ContractedShell s_shell;
    s_shell.shell_type = ShellType::S;
    s_shell.scale_factor = 1.0;
    
    // Add primitives (example from STO-3G)
    PrimitiveGaussian prim1;
    prim1.exponent = 3.42525091;
    prim1.coefficient = 0.15432897;
    prim1.coefficient_sp = 0.0;
    s_shell.primitives.push_back(prim1);
    
    PrimitiveGaussian prim2;
    prim2.exponent = 0.62391373;
    prim2.coefficient = 0.53532814;
    prim2.coefficient_sp = 0.0;
    s_shell.primitives.push_back(prim2);
    
    PrimitiveGaussian prim3;
    prim3.exponent = 0.16885540;
    prim3.coefficient = 0.44463454;
    prim3.coefficient_sp = 0.0;
    s_shell.primitives.push_back(prim3);
    
    assert(s_shell.primitives.size() == 3);
    assert(s_shell.shell_type == ShellType::S);
    
    // Test basis function count
    size_t num_funcs = util::get_num_basis_functions(ShellType::S);
    assert(num_funcs == 1);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test SP shell handling
 */
void test_sp_shell() {
    std::cout << "Testing SP shell handling... ";
    
    ContractedShell sp_shell;
    sp_shell.shell_type = ShellType::SP;
    sp_shell.scale_factor = 1.0;
    
    // SP shells have two coefficients
    PrimitiveGaussian prim;
    prim.exponent = 1.962;
    prim.coefficient = 0.38075;    // S coefficient
    prim.coefficient_sp = 0.15680; // P coefficient
    sp_shell.primitives.push_back(prim);
    
    assert(sp_shell.shell_type == ShellType::SP);
    assert(sp_shell.primitives[0].coefficient_sp != 0.0);
    
    // SP shell counts as 1 (S) + 3 (P) = 4 functions
    size_t num_funcs = util::get_num_basis_functions(ShellType::SP);
    assert(num_funcs == 4);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test molecular orbital structure
 */
void test_molecular_orbital() {
    std::cout << "Testing molecular orbital structure... ";
    
    MolecularOrbital mo;
    mo.energy = -11.2345;
    mo.spin = SpinType::Alpha;
    mo.occupation = 2.0;
    mo.symmetry = "A1g";
    
    // Add some coefficients
    mo.coefficients.push_back(0.1234);
    mo.coefficients.push_back(-0.5678);
    mo.coefficients.push_back(0.9012);
    
    assert(mo.spin == SpinType::Alpha);
    assert(mo.occupation == 2.0);
    assert(mo.coefficients.size() == 3);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test complete MoldenData structure
 */
void test_molden_data() {
    std::cout << "Testing complete MoldenData structure... ";
    
    MoldenData data;
    data.title = "Water molecule test";
    data.coord_unit = CoordinateUnit::Angstrom;
    data.basis_format = BasisFormat::Cartesian;
    
    // Add atoms (H2O)
    Atom o, h1, h2;
    
    o.element_symbol = "O";
    o.atom_index = 1;
    o.atomic_number = 8;
    o.x = 0.0f; o.y = 0.0f; o.z = 0.0f;
    
    h1.element_symbol = "H";
    h1.atom_index = 2;
    h1.atomic_number = 1;
    h1.x = 0.7570f; h1.y = 0.5860f; h1.z = 0.0f;
    
    h2.element_symbol = "H";
    h2.atom_index = 3;
    h2.atomic_number = 1;
    h2.x = -0.7570f; h2.y = 0.5860f; h2.z = 0.0f;
    
    data.atoms.push_back(o);
    data.atoms.push_back(h1);
    data.atoms.push_back(h2);
    
    assert(data.atoms.size() == 3);
    assert(data.coord_unit == CoordinateUnit::Angstrom);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test utility functions
 */
void test_utility_functions() {
    std::cout << "Testing utility functions... ";
    
    // Test coordinate conversions
    double angstrom_val = 1.0;
    double au_val = util::angstrom_to_au(angstrom_val);
    double back_to_angstrom = util::au_to_angstrom(au_val);
    assert(std::abs(back_to_angstrom - angstrom_val) < 1e-10);
    
    // Test shell type conversions
    assert(util::char_to_shell_type('s') == ShellType::S);
    assert(util::char_to_shell_type('p') == ShellType::P);
    assert(util::char_to_shell_type('d') == ShellType::D);
    
    assert(std::string(util::shell_type_to_string(ShellType::S)) == "s");
    assert(std::string(util::shell_type_to_string(ShellType::SP)) == "sp");
    
    // Test element symbol to atomic number
    assert(util::element_symbol_to_atomic_number("H") == 1);
    assert(util::element_symbol_to_atomic_number("C") == 6);
    assert(util::element_symbol_to_atomic_number("N") == 7);
    assert(util::element_symbol_to_atomic_number("O") == 8);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test basis function counting
 */
void test_basis_function_counting() {
    std::cout << "Testing basis function counting... ";
    
    // Cartesian basis
    assert(util::get_num_basis_functions(ShellType::S, BasisFormat::Cartesian) == 1);
    assert(util::get_num_basis_functions(ShellType::P, BasisFormat::Cartesian) == 3);
    assert(util::get_num_basis_functions(ShellType::D, BasisFormat::Cartesian) == 6);
    assert(util::get_num_basis_functions(ShellType::F, BasisFormat::Cartesian) == 10);
    
    // Spherical basis
    assert(util::get_num_basis_functions(ShellType::D, BasisFormat::Spherical) == 5);
    assert(util::get_num_basis_functions(ShellType::F, BasisFormat::Spherical) == 7);
    assert(util::get_num_basis_functions(ShellType::G, BasisFormat::Spherical) == 9);
    
    // SP shell
    assert(util::get_num_basis_functions(ShellType::SP, BasisFormat::Cartesian) == 4);
    assert(util::get_num_basis_functions(ShellType::SP, BasisFormat::Spherical) == 4);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Test validation functions
 */
void test_validation() {
    std::cout << "Testing validation functions... ";
    
    // Create valid data
    MoldenData data;
    data.coord_unit = CoordinateUnit::Angstrom;
    
    Atom atom;
    atom.element_symbol = "H";
    atom.atom_index = 1;
    atom.atomic_number = 1;
    atom.x = atom.y = atom.z = 0.0f;
    data.atoms.push_back(atom);
    
    // Should be valid
    assert(util::validate_molden_data(data) == true);
    
    // Create invalid data (no atoms)
    MoldenData invalid_data;
    assert(util::validate_molden_data(invalid_data) == false);
    
    std::cout << "PASS" << std::endl;
}

/**
 * @brief Demonstrate data structure usage with a realistic example
 */
void demonstrate_usage() {
    std::cout << "\n=== Demonstrating Molden Data Structure Usage ===" << std::endl;
    
    MoldenData data;
    data.title = "Hydrogen molecule (H2) - STO-3G basis";
    data.coord_unit = CoordinateUnit::Angstrom;
    data.basis_format = BasisFormat::Cartesian;
    
    // Add two hydrogen atoms
    Atom h1, h2;
    h1.element_symbol = "H";
    h1.atom_index = 1;
    h1.atomic_number = 1;
    h1.x = 0.0f; h1.y = 0.0f; h1.z = 0.0f;
    
    h2.element_symbol = "H";
    h2.atom_index = 2;
    h2.atomic_number = 1;
    h2.x = 0.0f; h2.y = 0.0f; h2.z = 0.74f; // Bond length ~0.74 Angstrom
    
    data.atoms.push_back(h1);
    data.atoms.push_back(h2);
    
    std::cout << "Created molecule: " << data.title << std::endl;
    std::cout << "Number of atoms: " << data.atoms.size() << std::endl;
    std::cout << "Coordinates in: " << (data.coord_unit == CoordinateUnit::Angstrom ? "Angstroms" : "AU") << std::endl;
    
    // Add basis set for first hydrogen
    AtomBasisSet basis1;
    basis1.atom_index = 1;
    
    ContractedShell s_shell;
    s_shell.shell_type = ShellType::S;
    s_shell.scale_factor = 1.0;
    
    // STO-3G for hydrogen (3 primitives)
    // Using designated initializers for clarity (C++20)
    s_shell.primitives.push_back(PrimitiveGaussian{
        .exponent = 3.42525091,
        .coefficient = 0.15432897,
        .coefficient_sp = 0.0
    });
    s_shell.primitives.push_back(PrimitiveGaussian{
        .exponent = 0.62391373,
        .coefficient = 0.53532814,
        .coefficient_sp = 0.0
    });
    s_shell.primitives.push_back(PrimitiveGaussian{
        .exponent = 0.16885540,
        .coefficient = 0.44463454,
        .coefficient_sp = 0.0
    });
    
    basis1.shells.push_back(s_shell);
    data.basis_sets.push_back(basis1);
    
    // Add basis set for second hydrogen (same as first)
    AtomBasisSet basis2 = basis1;
    basis2.atom_index = 2;
    data.basis_sets.push_back(basis2);
    
    // Calculate total basis functions
    data.total_basis_functions = util::calculate_total_basis_functions(data);
    std::cout << "Total basis functions: " << data.total_basis_functions << std::endl;
    
    // Add a molecular orbital (bonding sigma orbital)
    MolecularOrbital mo;
    mo.energy = -0.594;  // Hartrees
    mo.spin = SpinType::Alpha;
    mo.occupation = 2.0;
    mo.symmetry = "Sigma_g";
    mo.coefficients = {0.707, 0.707}; // Symmetric combination
    
    data.orbitals.push_back(mo);
    
    std::cout << "Number of molecular orbitals: " << data.orbitals.size() << std::endl;
    std::cout << "First MO energy: " << mo.energy << " Hartree" << std::endl;
    std::cout << "First MO occupation: " << mo.occupation << std::endl;
    
    // Validate the data
    bool is_valid = util::validate_molden_data(data);
    std::cout << "Data validation: " << (is_valid ? "PASS" : "FAIL") << std::endl;
}

int main() {
    std::cout << "=== Molden Data Structure Tests ===" << std::endl << std::endl;
    
    try {
        test_atom_structure();
        test_basis_structures();
        test_sp_shell();
        test_molecular_orbital();
        test_molden_data();
        test_utility_functions();
        test_basis_function_counting();
        test_validation();
        
        std::cout << "\n=== All tests PASSED ===" << std::endl;
        
        demonstrate_usage();
        
        std::cout << "\n=== Test suite completed successfully ===" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}
