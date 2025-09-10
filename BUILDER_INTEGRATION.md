# Molecule Builder Integration in VIAMD

This document describes the molecule builder component integration in VIAMD with RDKit support.

## Overview

VIAMD now includes a powerful molecule builder component that allows users to create and visualize molecules from SMILES strings using the RDKit C++ API. This enables on-demand molecular structure generation and seamless integration with VIAMD's visualization capabilities.

## Features

### Molecule Building Capabilities
- **SMILES parsing** with RDKit's robust chemical notation parser
- **3D structure generation** using distance geometry algorithms
- **Geometry optimization** with UFF (Universal Force Field) energy minimization
- **Hydrogen addition** for complete molecular representation
- **Real-time molecule generation** with immediate visualization

### User Interface Features
- **Menu integration** accessible via "Builder → Molecule Builder"
- **Interactive GUI** with dark theme matching VIAMD's aesthetic
- **Example molecule library** with quick-access buttons for common structures
- **Real-time feedback** with error reporting and molecule statistics
- **Seamless VIAMD integration** for direct loading into visualization pipeline

### Built-in Example Molecules
- **Water (H₂O)**: `O`
- **Ethanol**: `CCO`
- **Benzene**: `c1ccccc1`
- **Caffeine**: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- **Aspirin**: `CC(=O)OC1=CC=CC=C1C(=O)O`
- **Glucose**: `C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O`

## Installation

### Prerequisites
- RDKit library and development headers
- C++20 compatible compiler
- CMake 3.20 or higher

### Ubuntu/Debian Installation
```bash
# Install RDKit development libraries
sudo apt install librdkit-dev librdkit1

# Optional: Install additional RDKit components
sudo apt install rdkit-data rdkit-doc
```

### Fedora/CentOS Installation
```bash
# Enable EPEL repository if needed
sudo dnf install epel-release

# Install RDKit
sudo dnf install rdkit-devel rdkit
```

### Building VIAMD with Molecule Builder

#### Standard Build (with Builder enabled by default)
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

#### Custom Build Options
```bash
# Enable/disable builder module
cmake .. -DVIAMD_ENABLE_BUILDER=ON

# Enable/disable RDKit support (requires VIAMD_ENABLE_BUILDER=ON)
cmake .. -DVIAMD_ENABLE_BUILDER=ON -DVIAMD_ENABLE_RDKIT=ON

# Complete build with all options
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DVIAMD_ENABLE_BUILDER=ON \
  -DVIAMD_ENABLE_RDKIT=ON
```

### CMake Configuration Options

| Option | Default | Description |
|--------|---------|-------------|
| `VIAMD_ENABLE_BUILDER` | `ON` | Enable Molecule Builder Module |
| `VIAMD_ENABLE_RDKIT` | `ON` | Enable RDKit for molecule building |

### Graceful Degradation
If RDKit is not available:
- Component compiles with informative messages
- Builder menu item shows installation requirements
- No impact on other VIAMD functionality

## Usage

### 1. Access the Molecule Builder
- Navigate to **Builder → Molecule Builder** in the main menu
- The builder window will open with SMILES input interface

### 2. Enter SMILES String
- **Manual input**: Type any valid SMILES string in the text field
- **Example buttons**: Click pre-configured molecules for quick access
- **Real-time validation**: Invalid SMILES strings show immediate error feedback

### 3. Generate 3D Structure
- Click **"Build Molecule"** to generate 3D coordinates
- VIAMD automatically:
  - Parses the SMILES string using RDKit
  - Generates 3D coordinates with distance geometry
  - Adds implicit hydrogens
  - Optimizes geometry using UFF force field
  - Reports molecule statistics (atoms, bonds, molecular weight)

### 4. Load into VIAMD
- Click **"Load into VIAMD"** to import the generated molecule
- Molecule appears in the main visualization window
- Full VIAMD functionality available (representations, analysis, etc.)

## Technical Details

### RDKit Integration
The molecule builder uses several RDKit libraries:
- **RDKitGraphMol**: Core molecular graph functionality
- **RDKitSmilesParse**: SMILES string parsing and validation
- **RDKitDistGeomHelpers**: 3D coordinate generation using distance geometry
- **RDKitForceFieldHelpers**: Force field setup for geometry optimization
- **RDKitForceField**: UFF force field implementation
- **RDKitDescriptors**: Molecular property calculations

### Coordinate System and Format Conversion
- **RDKit**: Uses Angstrom units (consistent with VIAMD)
- **Automatic conversion**: Seamless translation between RDKit and VIAMD molecular representations
- **Memory management**: Efficient use of VIAMD's custom allocators

### Event-Driven Architecture
- **Component integration**: Follows VIAMD's event-driven component system
- **Menu registration**: Automatic integration with VIAMD's menu system
- **State management**: Proper cleanup and resource management

## Example Workflows

### Building Simple Molecules
1. **Water molecule**:
   ```
   SMILES: O
   Result: H₂O with proper 3D geometry
   ```

2. **Ethanol**:
   ```
   SMILES: CCO
   Result: C₂H₆O with optimized conformation
   ```

### Building Complex Molecules
1. **Caffeine**:
   ```
   SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
   Result: Complete caffeine structure with ~24 atoms
   ```

2. **Aspirin**:
   ```
   SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
   Result: Acetylsalicylic acid with proper stereochemistry
   ```

### Integration with VIAMD Analysis
1. **Generate molecule** using the builder
2. **Apply representations**: Ball & Stick, Space Filling, etc.
3. **Analyze structure**: Use VIAMD's analysis tools
4. **Export results**: Save coordinates or images

## Troubleshooting

### Common Issues

**"RDKit not found" during build**:
```bash
# Ensure RDKit development packages are installed
sudo apt install librdkit-dev librdkit1

# Verify CMake can find RDKit
cmake .. -DVIAMD_ENABLE_RDKIT=ON
```

**"Invalid SMILES string" error**:
- Check SMILES syntax (valence rules, ring closures, etc.)
- Use RDKit documentation for SMILES specification
- Try simpler molecules first to verify functionality

**"Failed to generate 3D coordinates"**:
- Some complex SMILES may fail distance geometry
- Try different conformer generation parameters
- Verify the molecule is chemically reasonable

**Builder menu not visible**:
```bash
# Ensure builder is enabled during compilation
cmake .. -DVIAMD_ENABLE_BUILDER=ON -DVIAMD_ENABLE_RDKIT=ON
```

### Performance Considerations
- **Large molecules**: 3D generation time increases with molecular size
- **Complex rings**: Strained ring systems may require longer optimization
- **Memory usage**: Proportional to molecular size and conformer generation

## Validation and Testing

### Recommended Test Molecules
```bash
# Simple molecules (fast generation)
O              # Water
CCO            # Ethanol
c1ccccc1       # Benzene

# Medium complexity
CC(C)O         # Isopropanol
CC(=O)O        # Acetic acid

# Complex molecules (test full functionality)
CN1C=NC2=C1C(=O)N(C(=O)N2C)C  # Caffeine
```

### Verification Steps
1. **SMILES parsing**: Verify no syntax errors
2. **3D generation**: Check reasonable coordinates
3. **VIAMD loading**: Confirm proper molecule display
4. **Statistics**: Verify atom count, molecular weight

## Limitations

### Current Limitations
- **SMILES format only**: No support for SDF, MOL, or other formats yet
- **Single conformer**: Generates one 3D structure per SMILES
- **UFF optimization only**: Limited to Universal Force Field
- **No fragment building**: Complete SMILES required

### Future Enhancements
- **Multiple input formats**: SDF, MOL, XYZ file import
- **Conformer ensembles**: Multiple 3D structures per molecule
- **Advanced optimization**: MMFF, GAFF force field support
- **Fragment-based building**: Interactive molecular construction
- **Batch processing**: Multiple molecule generation

## API Integration

### Component Architecture
The molecule builder follows VIAMD's component pattern:
- **Location**: `src/components/builder/builder.cpp`
- **Event integration**: Uses VIAMD's event system
- **Memory management**: VIAMD allocator compatibility
- **Menu integration**: Automatic registration with main menu

### Event Handling
The component responds to:
- `ViamdInitialize`: Component initialization and RDKit setup
- `ViamdWindowDrawMenu`: Menu integration and UI rendering
- `ViamdShutdown`: Cleanup and resource deallocation

### Extension Points
Developers can extend the builder by:
- **Adding input formats**: SDF, MOL file parsers
- **Custom force fields**: Alternative optimization methods
- **Specialized builders**: Protein, nucleic acid constructors
- **Analysis integration**: Property calculation, similarity searches

## Installation Verification

### Test RDKit Installation
```bash
# Check RDKit libraries
ldconfig -p | grep rdkit

# Verify headers
ls /usr/include/rdkit/GraphMol/

# Test VIAMD build
cd build
cmake .. -DVIAMD_ENABLE_BUILDER=ON -DVIAMD_ENABLE_RDKIT=ON
make -j$(nproc)
```

### Runtime Verification
1. **Start VIAMD**: `./viamd`
2. **Check menu**: Builder → Molecule Builder should be available
3. **Test building**: Try SMILES `CCO` and verify successful generation
4. **Load molecule**: Confirm structure appears in visualization

## Conclusion

The molecule builder integration provides VIAMD with powerful on-demand molecular structure generation capabilities. By leveraging RDKit's robust chemistry toolkit, users can quickly create molecules from SMILES notation and immediately visualize them using VIAMD's comprehensive analysis tools. This combination makes VIAMD an excellent platform for chemical education, drug discovery, and molecular research.