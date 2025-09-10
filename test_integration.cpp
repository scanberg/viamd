#include <viamd.h>
#include "../src/components/builder/builder.cpp"
#include "../src/components/openmm/openmm.cpp"
#include <iostream>

// Simple test to verify that the UFF minimization integration works
int main() {
    std::cout << "Testing UFF minimization integration..." << std::endl;
    
    // Create a minimal application state
    ApplicationState app_state = {};
    
    // Test that the builder clear function doesn't crash
    builder::MoleculeBuilder builder_test;
    builder_test.app_state = &app_state;
    builder_test.clear_molecule_from_viamd();
    
    std::cout << "✓ Clear molecule function works" << std::endl;
    
    // Test that the OpenMM interface function doesn't crash without OpenMM
    #ifdef VIAMD_ENABLE_OPENMM
    openmm_interface::minimize_energy_if_available(app_state);
    std::cout << "✓ OpenMM minimization interface works" << std::endl;
    #else
    std::cout << "✓ OpenMM not enabled - interface test skipped" << std::endl;
    #endif
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}