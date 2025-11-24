# Molden File Examples

This directory contains example Molden files for testing and demonstration purposes.

## Files

### h2_sto3g.molden
**Description**: Hydrogen molecule (H₂) with STO-3G basis set

**Contents**:
- 2 hydrogen atoms
- STO-3G basis (3 Gaussian primitives per atom)
- 2 molecular orbitals (bonding σ and antibonding σ*)
- Coordinates in Angstroms
- Simple system for basic testing

**Use Case**: Basic parser functionality testing

### h2o_sto3g.molden
**Description**: Water molecule (H₂O) with STO-3G basis set

**Contents**:
- 3 atoms (1 oxygen, 2 hydrogen)
- STO-3G basis with SP shells on oxygen
- 7 molecular orbitals (5 occupied, 2 virtual)
- Coordinates in Angstroms
- Demonstrates SP shell handling

**Use Case**: Testing SP shells (combined S and P shells)

## Using These Files

### With VIAMD (Phase 3 - Visualization)
**Drag-and-Drop Support**:
Simply drag and drop any .molden or .mold file into the VIAMD window to visualize it.

**File Menu**:
Use File → Open and select a Molden file to load it.

**What You'll See**:
- Atoms rendered with CPK coloring (color by element)
- Covalent bonds inferred from atomic distances
- Multiple representation modes:
  - Ball-and-stick (default)
  - Space-fill (van der Waals spheres)
  - Licorice (stick model)
  - And more...

**Camera Controls**:
- Left-click drag: Rotate view
- Right-click drag: Pan view
- Scroll: Zoom in/out
- Double-click atom: Center view on atom

**Supported Features**:
- Automatic coordinate conversion (Bohr ↔ Angstrom)
- Bond inference based on covalent radii
- Atom selection and highlighting
- Multiple color schemes

### With C++
```cpp
#include "molden.h"

// Load a file
molden::MoldenData data = molden::load_molden_file(
    "datasets/molden_examples/h2o_sto3g.molden"
);

// Check if successful
if (!data.atoms.empty()) {
    std::cout << "Loaded " << data.atoms.size() << " atoms" << std::endl;
}
```

### Running Tests
```bash
cd src
g++ -std=c++20 -I. molden_parser_test.cpp molden.cpp -o molden_test
./molden_test
```

## Molden Format Reference

For detailed information about the Molden format:
- See `../src/MOLDEN_FORMAT_NOTES.md` for format specification
- See `../src/MOLDEN_PARSER_USAGE.md` for parser usage
- See `../src/PHASE2_SUMMARY.md` for implementation details

## Creating Your Own Test Files

To create Molden files from quantum chemistry calculations:

**Gaussian**:
```
# Add to Gaussian input
%chk=molecule.chk
# After calculation:
formchk molecule.chk molecule.fchk
molden molecule.fchk
```

**ORCA**:
```
# Add to ORCA input
%output Print[P_MOs] 1 end
# Molden file generated automatically as molecule.molden
```

**GAMESS**:
```
# GAMESS automatically generates .dat file
# Convert to Molden format using:
molden2aim molecule.dat -molden
```

## File Format Overview

### Basic Structure
```
[Molden Format]
[Title]
 Title line

[Atoms] (AU|Angs)
Element  Index  Z  x  y  z

[GTO]
AtomIndex  0
shell_type  num_primitives  scale_factor
exponent  coefficient

[MO]
Ene= energy
Spin= Alpha|Beta
Occup= occupation
Sym= symmetry
coefficient_1
coefficient_2
...
```

## Additional Resources

- **Phase 1 Documentation**: Data structures and format notes
- **Phase 2 Documentation**: Parser implementation and usage
- **Test Suite**: `../src/molden_parser_test.cpp`

## Notes

These files are minimal examples for testing. Real quantum chemistry calculations may produce much larger files with:
- More atoms (hundreds to thousands)
- Larger basis sets (dozens of shells per atom)
- Many more molecular orbitals (hundreds)
- Additional sections ([FREQ], [GEOMETRIES], etc.)

The parser handles all these cases correctly.
