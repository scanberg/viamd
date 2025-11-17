# Phase 5: TREXIO Component Development

## Overview

Phase 5 implements a TREXIO visualization component for VIAMD that reuses existing infrastructure from the VeloxChem component.

## Implementation Status

### Completed ✅

1. **Basic Component Structure**
   - Created `src/components/trexio/trexio.cpp`
   - Follows VeloxChem component pattern
   - Conditional compilation with `#ifdef MD_TREXIO`
   - Proper initialization and shutdown

2. **UI Windows**
   - **TREXIO Summary Window**: Displays system information
     * Number of atoms
     * Number of electrons
     * Number of molecular orbitals
     * HOMO/LUMO indices and energies
   - **TREXIO Orbital Viewer**: Orbital visualization interface
     * MO selection dropdown
     * Quick HOMO/LUMO buttons
     * Rendering settings (iso value, resolution)

3. **Data Integration**
   - Loads TREXIO files via md_trexio API
   - Extracts quantum chemistry data (atoms, electrons, MOs)
   - Calculates HOMO/LUMO indices from occupation numbers
   - Converts energy units (Hartree → eV)

### In Progress 🚧

1. **Orbital Visualization**
   - Grid-based orbital rendering (reusing VeloxChem infrastructure)
   - GTO (Gaussian Type Orbital) extraction from TREXIO data
   - Volume rendering with adjustable iso-values
   - Multiple resolution levels (Low/Mid/High)

2. **Main View Integration**
   - Display atomic structure in main visualization window
   - Auto-show HOMO orbital on file load
   - Integration with existing camera and rendering systems

### Planned 📋

1. **Advanced Orbital Features**
   - Orbital color schemes (phase coloring)
   - Multiple orbital display
   - Orbital density plots
   - Export orbital data

2. **Extended Quantum Chemistry Data**
   - Basis set information display
   - Electron density visualization
   - Excited state data (when available in TREXIO files)
   - Response properties

3. **Performance Optimizations**
   - Asynchronous orbital computation
   - Cached grid data
   - GPU-accelerated GTO evaluation

## Architecture

### Component Design

The TREXIO component follows the same pattern as VeloxChem:

```cpp
struct Trexio {
    md_allocator_i* allocator;
    md_trexio_t* trexio_data;        // TREXIO file data
    
    struct qc_data {                  // Quantum chemistry data
        int num_atoms, num_electrons, num_mo;
        int homo_idx, lumo_idx;
        double* mo_energies;
        double* mo_occupations;
    };
    
    struct orbital {                  // Orbital rendering state
        GLuint tex_id;
        int selected_mo;
        float iso_value;
        VolumeRes resolution;
        md_grid_t grid;
        float* grid_data;
    };
    
    struct ui {                       // UI state
        bool show_summary;
        bool show_orbital;
    };
};
```

### Reusing VeloxChem Infrastructure

The component reuses several VeloxChem features:

1. **Grid Orbital Rendering**
   - Same volume rendering pipeline
   - Identical grid computation functions
   - Compatible texture formats

2. **GTO Evaluation**
   - md_gto.h functions for Gaussian evaluation
   - Same cutoff and sampling logic
   - Shared volume resolution settings

3. **UI Patterns**
   - Summary window format
   - Orbital selection interface
   - Rendering controls

## Integration with VIAMD

### File Loading

TREXIO files are loaded through the existing loader system:

1. User opens .trexio or .h5 file
2. Loader detects TREXIO format
3. md_trexio_parse_file() extracts data
4. TREXIO component initializes
5. Summary window auto-opens

### Main View Display

On file load:
1. Atomic structure displayed in main view
2. HOMO orbital automatically computed and shown
3. User can switch between MOs via summary window

### Menu Integration

Component windows accessible via:
- Window → TREXIO Summary
- Window → TREXIO Orbitals
- Auto-open on TREXIO file load

## Building

The component is controlled by the `VIAMD_ENABLE_TREXIO` CMake option:

```bash
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

When disabled, the component compiles to stubs that display a message about recompiling with TREXIO enabled.

## Testing

### Test Files

Use the generated PySCF test files:
- `test_data/h2_pyscf.h5` - Simple H2 molecule
- `test_data/h2o_pyscf.h5` - Water molecule

### Manual Testing

1. Build with `-DVIAMD_ENABLE_TREXIO=ON`
2. Run: `./viamd ../test_data/h2_pyscf.h5`
3. Verify:
   - Summary window shows correct atom/electron counts
   - HOMO/LUMO indices are calculated
   - MO energies displayed in eV
   - Orbital viewer interface is functional

## Limitations

### Current

- Orbital visualization not yet fully implemented (placeholder UI)
- GTO extraction needs TREXIO basis set data parsing
- No trajectory support (TREXIO stores single structures)
- Read-only (no file writing capability)

### By Design

- Reuses VeloxChem infrastructure (not a standalone implementation)
- Shares same volume rendering limitations
- Same grid resolution constraints

## Future Enhancements

1. **Phase 5.1**: Complete orbital visualization
   - Implement GTO extraction from TREXIO
   - Wire up volume rendering pipeline
   - Add orbital computation tasks

2. **Phase 5.2**: Enhanced visualization
   - Orbital phase coloring
   - Isosurface rendering options
   - Animation between orbitals

3. **Phase 5.3**: Advanced features
   - Excited state visualization
   - Transition density rendering
   - Integration with analysis tools

## References

- VeloxChem component: `src/components/veloxchem/veloxchem.cpp`
- TREXIO API: `ext/mdlib/src/md_trexio.h`
- Grid utilities: `gfx/volumerender_utils.h`
- GTO functions: `ext/mdlib/src/md_gto.h`

## Notes

This implementation prioritizes:
1. Reusing existing, proven code (VeloxChem component)
2. Minimal duplication
3. Consistent UI/UX with existing quantum chemistry features
4. Incremental development (working placeholders → full features)
