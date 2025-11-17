# TREXIO Support in VIAMD

## Overview

VIAMD now supports loading quantum chemistry data from TREXIO files. TREXIO is an open-source file format and library for storing and manipulating quantum chemistry data.

## What is TREXIO?

TREXIO (TREX Input/Output) is a file format designed for storing wave function parameters and matrix elements from quantum chemistry calculations. It supports:

- Atomic structures (nucleus coordinates, charges)
- Basis sets (Gaussian, Slater, numerical)
- Molecular orbitals
- Electron density
- Various quantum chemistry properties

The format is used by many quantum chemistry codes including:
- Quantum Package
- PySCF
- FHI-aims
- CP2K
- CHAMP
- And many others

## Supported Features

VIAMD's TREXIO integration currently supports:

### ✅ Implemented
- Loading atomic coordinates from TREXIO files
- Reading atomic charges and element information
- Displaying molecular structure
- Basic file format detection (.trexio and .h5 formats)

### 🚧 In Progress
- Molecular orbital visualization
- Basis set data extraction
- Electron density visualization

### 📋 Planned
- Excited state data
- Response properties
- Advanced quantum chemistry visualizations

## Building with TREXIO Support

### Prerequisites

TREXIO is automatically downloaded and built from source when you enable TREXIO support in VIAMD. No manual installation is required.

**Optional: HDF5 Library** - For HDF5 backend support (recommended)

HDF5 support enables TREXIO to read .h5 files. If HDF5 is not available, TREXIO will build with text backend support only (.trexio files).

```bash
# Ubuntu/Debian
sudo apt-get install libhdf5-dev

# macOS
brew install hdf5

# Conda
conda install -c conda-forge hdf5
```

### Compilation

TREXIO is automatically downloaded and built from source as part of the VIAMD build process:

```bash
# Initialize all submodules
git submodule update --init --recursive

# Apply the mdlib patch to enable TREXIO support
cd ext/mdlib
git apply ../../docs/mdlib_trexio.patch
cd ../..

# Configure and build VIAMD with TREXIO enabled
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

The build system will:
1. Automatically download the TREXIO 2.6.0 release tarball from GitHub
2. Detect if HDF5 is available on your system
3. Build TREXIO with HDF5 support (if available) or text-only backend
4. Link TREXIO statically into VIAMD

For advanced builds:

```bash
# Build with explicit HDF5 paths (if HDF5 is in a non-standard location)
cmake -DVIAMD_ENABLE_TREXIO=ON -DHDF5_ROOT=/path/to/hdf5 ..

# Combined with VeloxChem support
cmake -DVIAMD_ENABLE_TREXIO=ON -DVIAMD_ENABLE_VELOXCHEM=ON ..
```

## Usage

### Loading TREXIO Files

1. Launch VIAMD
2. Use File → Open or drag-and-drop
3. Select a TREXIO file (.trexio or .h5 extension)
4. The molecular structure will be loaded and displayed

### Supported File Extensions

- `.trexio` - TREXIO text format
- `.h5` - TREXIO HDF5 format (requires HDF5 backend)

**Note**: Both VeloxChem and TREXIO use `.h5` extension. VIAMD will try TREXIO format first if TREXIO support is enabled.

### Example Workflow

```bash
# Generate a TREXIO file using PySCF (Python example)
python -c "
from pyscf import gto, scf
import trexio

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='6-31g')
mf = scf.RHF(mol).run()

# Write to TREXIO
with trexio.File('h2_molecule.h5', 'w', trexio.TREXIO_HDF5) as f:
    trexio.write_nucleus_num(f, mol.natm)
    trexio.write_nucleus_charge(f, mol.atom_charges())
    trexio.write_nucleus_coord(f, mol.atom_coords())
    # ... write more data
"

# Load in VIAMD
./viamd h2_molecule.h5
```

## File Format Details

TREXIO files contain hierarchical data organized into groups:

- **nucleus**: Atomic positions, charges, labels
- **basis**: Basis set information (shells, primitives, coefficients)
- **ao**: Atomic orbital data
- **mo**: Molecular orbital coefficients and properties
- **electron**: Electron counts (up/down spin)
- And more...

For full TREXIO specification, see: https://trex-coe.github.io/trexio/

## Troubleshooting

### TREXIO download failed

**Error**: `Failed to download TREXIO` or network errors during CMake configuration

**Solution**: 
1. Check your internet connection
2. If behind a proxy, configure CMake proxy settings
3. Alternatively, manually download the tarball:
   ```bash
   cd /tmp
   wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
   ```
   Then reconfigure CMake (it will use the cached download)

### HDF5 backend not available

**Error**: Cannot open .h5 TREXIO files, only .trexio text files work

**Solution**:
1. Install HDF5 library (see Prerequisites section)
2. Reconfigure and rebuild VIAMD:
   ```bash
   cd build
   cmake -DVIAMD_ENABLE_TREXIO=ON ..
   make
   ```
The build system will automatically detect HDF5 and enable HDF5 support in TREXIO.

### Cannot open .h5 files

**Error**: File opens as VeloxChem instead of TREXIO (or vice versa)

**Solution**: 
- Ensure the file is a valid TREXIO file
- Try using `.trexio` extension for text format files
- Check that both `VIAMD_ENABLE_TREXIO` and `VIAMD_ENABLE_VELOXCHEM` are set correctly

## Limitations

Current limitations (will be addressed in future releases):

1. **Read-only**: Cannot write/export TREXIO files
2. **Limited groups**: Only nucleus, basis, mo, and electron groups are read
3. **No trajectory support**: TREXIO is for single structures only
4. **Basic visualization**: Advanced quantum chemistry features (orbitals, densities) are not yet visualized

## Contributing

To extend TREXIO support:

1. **Core Implementation**: `ext/mdlib/src/md_trexio.c`
2. **Header**: `ext/mdlib/src/md_trexio.h`
3. **Loader Integration**: `src/loader.cpp`
4. **UI Component** (future): `src/components/trexio/`

See `/tmp/TREXIO_INTEGRATION_DESIGN.md` for detailed technical documentation.

## References

- [TREXIO Documentation](https://trex-coe.github.io/trexio/)
- [TREXIO GitHub](https://github.com/TREX-CoE/trexio)
- [TREXIO Paper](https://doi.org/10.1063/5.0148161)
- [VIAMD Wiki](https://github.com/scanberg/viamd/wiki)

## Support

For questions or issues:
- VIAMD Issues: https://github.com/scanberg/viamd/issues
- TREXIO Issues: https://github.com/TREX-CoE/trexio/issues

## License

TREXIO support in VIAMD follows the same license as VIAMD. The TREXIO library itself is licensed under the 3-clause BSD license.
