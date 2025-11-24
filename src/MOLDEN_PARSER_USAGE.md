# Molden File Parser - Usage Guide

## Overview

This document describes how to use the Molden file parser implemented in Phase 2. The parser can read Molden format files and extract atomic coordinates, basis set information, and molecular orbital data.

## Quick Start

### Basic Usage

```cpp
#include "molden.h"

// Load a Molden file
molden::MoldenData data = molden::load_molden_file("molecule.molden");

// Check if loading succeeded
if (data.atoms.empty()) {
    std::cerr << "Failed to load Molden file" << std::endl;
    return;
}

// Access atomic coordinates
for (const auto& atom : data.atoms) {
    std::cout << atom.element_symbol << ": " 
              << atom.x << ", " << atom.y << ", " << atom.z << std::endl;
}

// Access molecular orbitals
for (const auto& mo : data.orbitals) {
    std::cout << "Orbital energy: " << mo.energy << " Ha" << std::endl;
    std::cout << "Occupation: " << mo.occupation << std::endl;
}
```

### With Error Handling

```cpp
std::string error_message;
molden::MoldenData data = molden::parse_molden_file("molecule.molden", &error_message);

if (!error_message.empty()) {
    std::cerr << "Parse error: " << error_message << std::endl;
    return;
}

// Validate the data
if (!molden::util::validate_molden_data(data)) {
    std::cerr << "Data validation failed" << std::endl;
    return;
}
```

## API Reference

### Core Functions

#### `load_molden_file(const char* filepath)`
Primary entry point for loading Molden files.

**Parameters:**
- `filepath`: Path to the Molden file

**Returns:**
- `MoldenData` structure with parsed data, or empty structure on error

**Example:**
```cpp
molden::MoldenData data = molden::load_molden_file("h2o.molden");
```

#### `parse_molden_file(const std::string& filename, std::string* error_msg)`
Parse a Molden file with detailed error reporting.

**Parameters:**
- `filename`: Path to the Molden file
- `error_msg`: Optional pointer to string for error messages

**Returns:**
- `MoldenData` structure with parsed data

**Example:**
```cpp
std::string error;
molden::MoldenData data = molden::parse_molden_file("molecule.molden", &error);
if (!error.empty()) {
    std::cerr << "Error: " << error << std::endl;
}
```

#### `parse_molden_string(const std::string& str, std::string* error_msg)`
Parse Molden data from a string buffer.

**Parameters:**
- `str`: String containing Molden file contents
- `error_msg`: Optional pointer to string for error messages

**Returns:**
- `MoldenData` structure with parsed data

**Example:**
```cpp
std::string molden_content = read_from_network();
molden::MoldenData data = molden::parse_molden_string(molden_content);
```

### Validation Functions

#### `util::validate_molden_data(const MoldenData& data)`
Validate the consistency of parsed Molden data.

**Checks:**
- Atoms exist and have valid indices
- Basis sets reference valid atoms
- Exponents are positive
- Occupation numbers are in [0.0, 2.0]
- Coefficient counts match basis function counts

**Returns:**
- `true` if data is valid, `false` otherwise

**Example:**
```cpp
if (molden::util::validate_molden_data(data)) {
    std::cout << "Data is valid" << std::endl;
}
```

#### `util::calculate_total_basis_functions(const MoldenData& data)`
Calculate the total number of basis functions.

**Returns:**
- Total number of basis functions

**Example:**
```cpp
size_t num_funcs = molden::util::calculate_total_basis_functions(data);
std::cout << "Total basis functions: " << num_funcs << std::endl;
```

## Data Structures

### MoldenData

Main container for all Molden file data:

```cpp
struct MoldenData {
    std::string title;                           // Title from [Title] section
    CoordinateUnit coord_unit;                   // Angstrom or AtomicUnit
    BasisFormat basis_format;                    // Cartesian or Spherical
    
    std::vector<Atom> atoms;                     // Atomic coordinates
    std::vector<AtomBasisSet> basis_sets;        // Basis set data
    std::vector<MolecularOrbital> orbitals;      // Molecular orbitals
    
    size_t total_basis_functions;                // Total number of basis functions
    size_t num_alpha_orbitals;                   // Number of alpha orbitals
    size_t num_beta_orbitals;                    // Number of beta orbitals
};
```

### Atom

Atomic coordinate data:

```cpp
struct Atom {
    std::string element_symbol;  // Element symbol (e.g., "C", "H")
    int32_t     atom_index;      // 1-based index from file
    int32_t     atomic_number;   // Atomic number (Z)
    float       x, y, z;         // Coordinates
};
```

### MolecularOrbital

Molecular orbital data:

```cpp
struct MolecularOrbital {
    double energy;                       // Orbital energy (a.u.)
    SpinType spin;                       // Alpha or Beta
    double occupation;                   // Occupation number (0.0 to 2.0)
    std::string symmetry;                // Symmetry label
    std::vector<double> coefficients;    // MO coefficients
};
```

## Supported Features

### Molden Sections

- ✅ `[Molden Format]` - Header (optional)
- ✅ `[Title]` - Title line (optional)
- ✅ `[Atoms]` - Atomic coordinates (required)
- ✅ `[GTO]` - Gaussian basis sets (required for MO visualization)
- ✅ `[MO]` - Molecular orbitals (optional)
- ✅ `[5D]`, `[5D7F]`, `[5D10F]` - Spherical harmonics tags

### Format Variations

- ✅ Coordinate units: Angstrom and Atomic Units (AU)
- ✅ Basis formats: Cartesian and Spherical
- ✅ SP shells (combined S+P shells)
- ✅ Missing atomic numbers (inferred from element symbols)
- ✅ Missing symmetry labels in MO section
- ✅ Fortran-style scientific notation (D instead of E)
- ✅ Orbital index lines (skipped automatically)
- ✅ Multiple coefficients per line
- ✅ Case-insensitive section names

### Edge Cases

The parser handles various edge cases documented in MOLDEN_FORMAT_NOTES.md:

- Missing or optional fields
- Non-sequential atom indices
- Different whitespace formatting
- Scale factors in GTO section
- Virtual (unoccupied) orbitals
- Fractional occupation numbers

## Examples

### Example 1: Simple H₂ Molecule

```cpp
molden::MoldenData data = molden::load_molden_file("h2.molden");

std::cout << "Number of atoms: " << data.atoms.size() << std::endl;
std::cout << "Coordinate unit: " 
          << (data.coord_unit == molden::CoordinateUnit::Angstrom ? "Angstrom" : "AU") 
          << std::endl;

// Calculate bond length
if (data.atoms.size() == 2) {
    float dx = data.atoms[1].x - data.atoms[0].x;
    float dy = data.atoms[1].y - data.atoms[0].y;
    float dz = data.atoms[1].z - data.atoms[0].z;
    float bond_length = std::sqrt(dx*dx + dy*dy + dz*dz);
    std::cout << "H-H bond length: " << bond_length << " Angstrom" << std::endl;
}
```

### Example 2: Analyze Molecular Orbitals

```cpp
molden::MoldenData data = molden::load_molden_file("molecule.molden");

// Find HOMO (Highest Occupied Molecular Orbital)
double homo_energy = -std::numeric_limits<double>::infinity();
size_t homo_idx = 0;

for (size_t i = 0; i < data.orbitals.size(); ++i) {
    const auto& mo = data.orbitals[i];
    if (mo.occupation > 0.0 && mo.energy > homo_energy) {
        homo_energy = mo.energy;
        homo_idx = i;
    }
}

std::cout << "HOMO energy: " << homo_energy << " Ha" << std::endl;
std::cout << "HOMO occupation: " << data.orbitals[homo_idx].occupation << std::endl;
```

### Example 3: Basis Set Information

```cpp
molden::MoldenData data = molden::load_molden_file("molecule.molden");

std::cout << "Basis set information:" << std::endl;
std::cout << "Total basis functions: " << data.total_basis_functions << std::endl;
std::cout << "Basis format: " 
          << (data.basis_format == molden::BasisFormat::Cartesian ? "Cartesian" : "Spherical")
          << std::endl;

// Print shells for each atom
for (const auto& basis : data.basis_sets) {
    std::cout << "Atom " << basis.atom_index << " has " 
              << basis.shells.size() << " shells:" << std::endl;
    
    for (const auto& shell : basis.shells) {
        std::cout << "  " << molden::util::shell_type_to_string(shell.shell_type)
                  << " shell with " << shell.primitives.size() << " primitives" << std::endl;
    }
}
```

### Example 4: Convert Coordinates

```cpp
molden::MoldenData data = molden::load_molden_file("molecule.molden");

// Convert to Angstroms if in AU
if (data.coord_unit == molden::CoordinateUnit::AtomicUnit) {
    for (auto& atom : data.atoms) {
        atom.x = static_cast<float>(molden::util::au_to_angstrom(atom.x));
        atom.y = static_cast<float>(molden::util::au_to_angstrom(atom.y));
        atom.z = static_cast<float>(molden::util::au_to_angstrom(atom.z));
    }
    data.coord_unit = molden::CoordinateUnit::Angstrom;
}
```

## Testing

### Unit Tests

Run the comprehensive test suite:

```bash
cd src
g++ -std=c++20 -I. molden_parser_test.cpp molden.cpp -o molden_parser_test
./molden_parser_test
```

All tests should pass:
- H₂ parsing
- H₂O with SP shells
- AU vs Angstrom coordinates
- Missing atomic numbers
- Multiple orbitals with different spins
- Missing symmetry labels
- Spherical basis formats
- Fortran D notation
- File parsing
- Error handling

### Example Files

Test files are provided in `datasets/molden_examples/`:

- `h2_sto3g.molden` - Simple H₂ molecule with STO-3G basis
- `h2o_sto3g.molden` - Water molecule with SP shells

## Performance

### Memory Usage

The parser allocates memory dynamically:

- **Atoms**: ~40 bytes per atom
- **Basis functions**: ~24 bytes per primitive
- **MO coefficients**: 8 bytes × num_basis_functions × num_orbitals

For a typical small molecule (30 atoms, 100 basis functions, 100 orbitals):
- Total memory: ~80 KB

### Parsing Speed

Typical parsing times on modern hardware:
- Small molecules (<100 atoms): <10 ms
- Medium molecules (100-1000 atoms): 10-100 ms
- Large systems (>1000 atoms): 100-500 ms

The parser is I/O bound for most files.

## Error Handling

### Common Errors

1. **File not found**
   ```
   Failed to open file: molecule.molden
   ```

2. **Invalid atom line**
   ```
   Invalid atom line (expected at least 5 fields): H 1
   ```

3. **Unknown shell type**
   ```
   Unknown shell type: x
   ```

4. **Coefficient count mismatch**
   ```
   Validation fails if coefficient count doesn't match basis function count
   ```

### Best Practices

1. Always check for errors:
   ```cpp
   std::string error;
   molden::MoldenData data = molden::parse_molden_file("file.molden", &error);
   if (!error.empty()) {
       // Handle error
   }
   ```

2. Validate data before use:
   ```cpp
   if (!molden::util::validate_molden_data(data)) {
       // Data is invalid
   }
   ```

3. Check for empty data:
   ```cpp
   if (data.atoms.empty()) {
       // No atoms parsed
   }
   ```

## Future Enhancements (Phase 3)

Planned features for Phase 3:

- Integration with VIAMD loader system
- Conversion to `md_system_t` structure
- Support for visualization in VIAMD UI
- Orbital rendering with isosurfaces
- Support for multiple geometries
- Frequency data parsing

## References

- Phase 1 Summary: `PHASE1_SUMMARY.md`
- Format Notes: `MOLDEN_FORMAT_NOTES.md`
- Data Structures: `molden.h`
- Implementation: `molden.cpp`
- Tests: `molden_parser_test.cpp`
