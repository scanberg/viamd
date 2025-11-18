# TREXIO Support in VIAMD

> **⚠️ IMPORTANT STATUS NOTICE**  
> TREXIO support is currently **under development** and is **not yet fully functional**.  
> Building with `-DVIAMD_ENABLE_TREXIO=ON` will fail due to API compatibility issues.  
> See [TREXIO_KNOWN_ISSUES.md](TREXIO_KNOWN_ISSUES.md) for details and required fixes.

## Overview

VIAMD has partial support for loading quantum chemistry data from TREXIO files. TREXIO is an open-source file format and library for storing and manipulating quantum chemistry data.

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

### Building VIAMD with TREXIO Support

Once TREXIO is installed, build VIAMD with TREXIO support enabled:

```bash
# Initialize all submodules
git submodule update --init --recursive

# Apply the mdlib patch to enable TREXIO support
./scripts/apply_mdlib_trexio_patch.sh

# Configure and build VIAMD with TREXIO enabled
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

**Note:** The patch application script (`./scripts/apply_mdlib_trexio_patch.sh`) can be run multiple times safely. If you encounter issues with the patch, simply run the script again and it will automatically clean up and reapply the patch.

**Alternative: Manual Patch Application**

If you prefer to apply the patch manually or encounter issues with the script:

```bash
# From the repository root
cd ext/mdlib

# If patch was partially applied, clean up first:
git checkout -- CMakeLists.txt
rm -f src/md_trexio.c src/md_trexio.h

# Apply the patch
git apply ../../docs/mdlib_trexio_original.patch
cd ../..
```

The build system will:
1. Detect the installed TREXIO library via pkg-config
2. Link against the TREXIO library
3. Enable TREXIO file format support in VIAMD

**Troubleshooting CMake Configuration:**

If CMake cannot find TREXIO:

```bash
# Ensure TREXIO is in pkg-config path
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH

# Verify pkg-config can find TREXIO
pkg-config --modversion trexio

# On Linux, ensure library cache is updated
sudo ldconfig

# Reconfigure CMake
cd build && rm -rf * && cmake -DVIAMD_ENABLE_TREXIO=ON ..
```

**Advanced Build Options:**

```bash
# Specify custom TREXIO installation path
cmake -DVIAMD_ENABLE_TREXIO=ON -DTREXIO_ROOT=/custom/path ..

# Combined with VeloxChem support
cmake -DVIAMD_ENABLE_TREXIO=ON -DVIAMD_ENABLE_VELOXCHEM=ON ..

# Specify custom HDF5 path (if needed for TREXIO)
cmake -DVIAMD_ENABLE_TREXIO=ON -DHDF5_ROOT=/custom/hdf5/path ..
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

### TREXIO library not found during CMake configuration

**Error**: `TREXIO support requested but TREXIO library not found`

**Solution**:
1. Verify TREXIO is installed:
   ```bash
   pkg-config --modversion trexio
   # Should output: 2.6.0 (or your installed version)
   ```

2. If not installed, install TREXIO following the instructions in the Prerequisites section above.

3. If installed but not found, ensure pkg-config can locate it:
   ```bash
   # Add TREXIO installation to PKG_CONFIG_PATH
   export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
   
   # On Linux, update library cache
   sudo ldconfig
   
   # Verify
   pkg-config --libs --cflags trexio
   ```

4. Reconfigure CMake:
   ```bash
   cd build && rm -rf * && cmake -DVIAMD_ENABLE_TREXIO=ON ..
   ```

### Patch application failed

**Error**: `error: corrupt patch at line XXX` or `Failed to apply patch`

**Solution**:
1. Clean up any partial application:
   ```bash
   cd ext/mdlib
   git checkout -- CMakeLists.txt
   git clean -fd  # Remove any new TREXIO files
   cd ../..
   ```

2. Use the original (working) patch:
   ```bash
   cd ext/mdlib
   git apply ../../docs/mdlib_trexio_original.patch
   cd ../..
   ```

3. If the patch still fails, check that your mdlib submodule is at the correct commit:
   ```bash
   cd ext/mdlib
   git log --oneline -1  # Should show: 06da5a3 vlx fixes
   cd ../..
   ```

### HDF5 backend not available in TREXIO

**Error**: Cannot open .h5 TREXIO files, only .trexio text files work

**Cause**: TREXIO was built without HDF5 support

**Solution**:
1. Install HDF5 library (see Prerequisites section above)
2. Rebuild and reinstall TREXIO:
   ```bash
   cd /tmp/trexio-2.6.0  # Or wherever you built TREXIO
   make clean
   ./configure --prefix=/usr/local
   make -j4
   sudo make install
   sudo ldconfig  # Linux only
   ```
3. Verify HDF5 support is enabled:
   ```bash
   # Check TREXIO configuration output during ./configure
   # Should show: "Compilation with HDF5: yes"
   ```
4. Rebuild VIAMD:
   ```bash
   cd /path/to/viamd/build
   rm -rf *
   cmake -DVIAMD_ENABLE_TREXIO=ON ..
   make
   ```

### Cannot open .h5 files / File format detection issues

**Error**: File opens as VeloxChem instead of TREXIO (or vice versa)

**Cause**: Both TREXIO and VeloxChem use `.h5` extension

**Solution**: 
- VIAMD tries TREXIO format first if TREXIO support is enabled
- Ensure the file is a valid TREXIO file by checking with the TREXIO library
- Use `.trexio` extension for text format files to avoid ambiguity
- Check file validity:
  ```bash
  # Using Python with trexio library
  python3 -c "import trexio; f = trexio.File('your_file.h5', 'r'); print(f'Valid TREXIO file: {trexio.inquire(\"your_file.h5\")}')"
  ```

### Runtime errors when loading TREXIO files

**Error**: Crash or errors when opening a TREXIO file in VIAMD

**Solution**:
1. Verify the file is valid:
   ```bash
   # Check if file can be opened by TREXIO library
   python3 -c "import trexio; f = trexio.File('your_file.h5', 'r'); print('OK')"
   ```

2. Check VIAMD was built correctly:
   ```bash
   # Verify TREXIO support is enabled
   ./viamd --help | grep -i trexio
   
   # Check linked libraries
   ldd ./viamd | grep trexio  # Linux
   otool -L ./viamd | grep trexio  # macOS
   ```

3. Enable verbose logging (if available) to see detailed error messages

### Known Limitations

See the Limitations section below for current feature limitations.

## Limitations

### Current Implementation Limitations

1. **Read-only**: VIAMD cannot write or export TREXIO files, only read them
2. **Limited data groups**: Currently reads:
   - Nucleus data (coordinates, charges, labels)
   - Electron configuration (up/down spin)
   - Basic system metadata
3. **No trajectory support**: TREXIO format stores single molecular structures, not trajectories
4. **No molecular orbital visualization**: MO data is parsed but not yet visualized
5. **No basis set display**: Basis set information is not currently displayed in the UI
6. **Text format compatibility**: Some manually-created text format files may not be compatible with TREXIO 2.6.0 library

### Known Issues

1. **HDF5 vs Text format**: If TREXIO was built without HDF5 support, only `.trexio` text format files can be loaded
2. **File extension ambiguity**: Both TREXIO and VeloxChem use `.h5` extension; VIAMD tries TREXIO first
3. **Large files**: Very large TREXIO files may be slow to load (performance optimizations planned)
4. **Memory usage**: Loading large molecular systems may require significant memory

### Planned Features

Future releases will add:
- Molecular orbital visualization (grid-based rendering)
- Basis set information display
- Electron density visualization
- Excited state data support (when available in files)
- Response properties visualization
- Export/conversion capabilities

## Contributing and Development

### Code Structure

TREXIO support in VIAMD consists of several components:

1. **Core Loader** (`ext/mdlib/src/md_trexio.c`, `md_trexio.h` - added via patch)
   - Parses TREXIO files using the TREXIO library
   - Implements `md_system_loader_i` interface
   - Converts TREXIO data to VIAMD internal structures
   
2. **Loader Integration** (`src/loader.cpp`)
   - Registers TREXIO loader with the file loading system
   - Handles file extension detection (`.trexio`, `.h5`)
   - Conditional compilation with `#if MD_TREXIO`

3. **UI Component** (`src/components/trexio/trexio.cpp`)
   - Displays TREXIO file information
   - Provides orbital and quantum chemistry data visualization
   - Future: advanced QC visualizations

4. **mdlib Patch** (`docs/mdlib_trexio_original.patch`)
   - Adds TREXIO support to the mdlib library
   - Must be applied before building
   - Can be applied multiple times safely via script

### For Developers

To extend TREXIO support:

1. **Adding new TREXIO data groups**:
   - Edit `ext/mdlib/src/md_trexio.c` to read additional groups
   - Update `md_trexio_t` structure to store new data
   - Add accessor functions in `md_trexio.h`
   - Update the patch file: `cd ext/mdlib && git diff > ../../docs/mdlib_trexio_new.patch`

2. **Modifying the UI**:
   - Edit `src/components/trexio/trexio.cpp`
   - Follow existing patterns from VeloxChem component
   - Use ImGui for UI elements

3. **Testing changes**:
   - Use test files in `test_data/` directory
   - Create new test TREXIO files with `test_data/create_test_trexio.py` or `create_pyscf_trexio.py`
   - Unit tests in `ext/mdlib/unittest/test_trexio.c` (if patch includes them)

### Maintaining the Patch

The mdlib TREXIO patch adds TREXIO support to the mdlib submodule. To update it:

```bash
# 1. Make changes to mdlib TREXIO files
cd ext/mdlib
# ... edit src/md_trexio.c, src/md_trexio.h, or CMakeLists.txt ...

# 2. Generate new patch
git diff > ../../docs/mdlib_trexio_updated.patch

# 3. Test the new patch
cd ..
git clone https://github.com/scanberg/mdlib mdlib_test
cd mdlib_test
git checkout 06da5a3  # Or current commit
git apply ../../docs/mdlib_trexio_updated.patch

# 4. If successful, replace the old patch
cd ../../docs
mv mdlib_trexio_updated.patch mdlib_trexio_original.patch
```

### Build System Integration

TREXIO support is controlled by CMake options:

- **VIAMD level**: `VIAMD_ENABLE_TREXIO=ON` (in root `CMakeLists.txt`)
- **mdlib level**: `MD_ENABLE_TREXIO=ON` (propagated automatically)

When enabled:
1. CMake tries to find TREXIO via `find_package(TREXIO)`
2. If found, adds compile definition `MD_TREXIO`
3. Links against TREXIO library
4. Includes TREXIO source files in build

### Testing

Test the TREXIO implementation:

```bash
# Build with TREXIO support
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make

# Test with sample files
./viamd ../test_data/h2_molecule.trexio
./viamd ../test_data/h2o_molecule.trexio
./viamd ../test_data/ch4_molecule.trexio

# Run validation script (if available)
cd ../test_data
./validate_trexio.sh
```

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
