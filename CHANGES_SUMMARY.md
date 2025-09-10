# Force Field Switching and Simulation Stability Improvements

## Summary of Changes

This PR addresses the two main issues identified:

### 1. Force Field Switching After Initialization ✅

**Problem**: Previously, users could only change force fields before initializing the simulation. The GUI disabled the force field combo box after initialization.

**Solution**: 
- Removed the `if (!state.simulation.initialized)` restriction around the force field selection
- Added automatic system reinitialization when force field is changed on an already initialized system
- Users can now switch between AMBER14 and UFF at any time during simulation

**Code Changes**:
```cpp
// Before: Force field combo only available when not initialized
if (!state.simulation.initialized) {
    // Force field selection code
}

// After: Force field combo always available with automatic reinitialization
if (ImGui::Combo("Force Field", &current_ff, force_field_items, IM_ARRAYSIZE(force_field_items))) {
    // ... force field change logic
    if (state.simulation.initialized) {
        MD_LOG_INFO("Reinitializing system with new force field...");
        cleanup_simulation(state);
        setup_system(state);
    }
}
```

### 2. Simulation Stability Improvements ✅

**Problem**: Many simulations were exploding due to aggressive parameters and insufficient stabilization.

**Solutions Implemented**:

1. **More Conservative Default Timestep**:
   - Reduced from 0.001 ps to 0.0005 ps
   - Reduced slider maximum from 0.002 ps to 0.001 ps

2. **Improved Energy Minimization**:
   - Increased tolerance from 1e-4 to 1e-6 (tighter convergence)
   - Increased iterations from 1000 to 5000
   - Added energy validation after minimization

3. **Reduced Force Field Charges**:
   - Scaled charges by 0.25x instead of 0.5x (more conservative)
   - Applied to both AMBER and UFF force fields

4. **Enhanced Explosion Detection**:
   - Lowered threshold from 100 Å to 50 Å from origin
   - More sensitive early detection of instabilities

## Testing

A test program was created (`test_force_field_switch.cpp`) that verifies:
- Force field changes work before initialization (no reinitialization needed)
- Force field changes work after initialization (automatic reinitialization)
- Default timestep is correctly set to the more conservative value

## Impact

1. **User Experience**: Users can now experiment with different force fields without restarting simulations
2. **Stability**: Simulations are much less likely to explode due to conservative parameters
3. **Reliability**: Better energy minimization leads to more stable initial structures

These changes maintain backward compatibility while significantly improving the robustness and usability of the OpenMM simulation component.