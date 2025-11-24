/**
 * @file test_molden_loader.cpp
 * @brief Test program to verify Molden system loader integration
 * 
 * This program tests the Molden loader by loading example Molden files
 * and verifying the conversion to md_system_t.
 */

#include "molden.h"
#include "md_molden_loader.h"

#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_str.h>

#include <iostream>
#include <iomanip>

void print_system_info(const md_system_t* sys) {
    std::cout << "\n=== System Information ===" << std::endl;
    std::cout << "Atoms: " << sys->atom.count << std::endl;
    std::cout << "Atom Types: " << sys->atom.type.count << std::endl;
    std::cout << "Bonds: " << sys->bond.count << std::endl;
    
    std::cout << "\n--- Atom Types ---" << std::endl;
    for (size_t i = 0; i < sys->atom.type.count; ++i) {
        str_t name = LBL_TO_STR(sys->atom.type.name[i]);
        std::cout << "  " << i << ": " 
                  << std::string(name.ptr, name.len)
                  << " (Z=" << (int)sys->atom.type.z[i] << ")" << std::endl;
    }
    
    std::cout << "\n--- Atoms (first 10) ---" << std::endl;
    size_t max_atoms = std::min((size_t)10, sys->atom.count);
    for (size_t i = 0; i < max_atoms; ++i) {
        md_atom_type_idx_t type_idx = sys->atom.type_idx[i];
        str_t type_name = LBL_TO_STR(sys->atom.type.name[type_idx]);
        std::cout << "  " << i << ": " 
                  << std::string(type_name.ptr, type_name.len)
                  << " @ (" << std::fixed << std::setprecision(3)
                  << sys->atom.x[i] << ", "
                  << sys->atom.y[i] << ", "
                  << sys->atom.z[i] << ")" << std::endl;
    }
    
    std::cout << "\n--- Bonds (first 10) ---" << std::endl;
    size_t max_bonds = std::min((size_t)10, sys->bond.count);
    for (size_t i = 0; i < max_bonds; ++i) {
        md_atom_idx_t atom1 = sys->bond.atom_idx[2*i];
        md_atom_idx_t atom2 = sys->bond.atom_idx[2*i + 1];
        
        md_atom_type_idx_t type1 = sys->atom.type_idx[atom1];
        md_atom_type_idx_t type2 = sys->atom.type_idx[atom2];
        str_t name1 = LBL_TO_STR(sys->atom.type.name[type1]);
        str_t name2 = LBL_TO_STR(sys->atom.type.name[type2]);
        
        std::cout << "  " << i << ": " 
                  << std::string(name1.ptr, name1.len) << atom1
                  << " - " 
                  << std::string(name2.ptr, name2.len) << atom2
                  << std::endl;
    }
    std::cout << std::endl;
}

bool test_molden_file(const char* filepath) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Testing: " << filepath << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Get the loader
    md_system_loader_i* loader = md_molden_system_loader();
    if (!loader) {
        std::cerr << "ERROR: Failed to get Molden loader" << std::endl;
        return false;
    }
    
    // Create system structure
    md_system_t sys = {0};
    
    // Get allocator
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Load the file
    str_t filename = str_from_cstr(filepath);
    bool success = loader->init_from_file(&sys, filename, nullptr, alloc);
    
    if (!success) {
        std::cerr << "ERROR: Failed to load Molden file: " << filepath << std::endl;
        return false;
    }
    
    // Print system information
    print_system_info(&sys);
    
    // Cleanup
    md_system_free(&sys, alloc);
    
    std::cout << "SUCCESS: Loaded and verified " << filepath << std::endl;
    return true;
}

int main(int argc, char** argv) {
    std::cout << "\n=== Molden Loader Test ===" << std::endl;
    std::cout << "This test verifies the Molden file loader integration" << std::endl;
    
    // Test with example files
    bool all_passed = true;
    
    // Test H2O
    all_passed &= test_molden_file("datasets/molden_examples/h2o_sto3g.molden");
    
    // Test H2
    all_passed &= test_molden_file("datasets/molden_examples/h2_sto3g.molden");
    
    if (all_passed) {
        std::cout << "\n=== ALL TESTS PASSED ===" << std::endl;
        return 0;
    } else {
        std::cout << "\n=== SOME TESTS FAILED ===" << std::endl;
        return 1;
    }
}
