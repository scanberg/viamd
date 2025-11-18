# TREXIO GUI User Guide

## Overview

The TREXIO GUI panel provides a graphical interface for loading and visualizing quantum chemistry data from TREXIO files. It enables users to view molecular metadata and compute molecular orbital grids for visualization.

## Features

- **File Loading**: Load TREXIO files (`.h5` or `.trexio` format)
- **Metadata Display**: View system information, basis set details, and molecular orbital data
- **Orbital Grid Computation**: Evaluate molecular orbitals on 3D grids for visualization
- **Export**: Save computed grids to binary files for external use

## Requirements

### Build Dependencies

- TREXIO library (v2.6.0 or later)
- HDF5 library
- CMake 3.20+
- C++17 compiler

### Installation

1. **Install TREXIO library**:
```bash
wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
tar -xzf trexio-2.6.0.tar.gz
cd trexio-2.6.0
./configure --prefix=/usr/local
make -j4
sudo make install
sudo ldconfig
```

2. **Build viamd with TREXIO support**:
```bash
cd /path/to/viamd
git submodule update --init --recursive
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make -j4
```

## Usage

### Accessing the TREXIO Panel

1. Launch viamd: `./bin/viamd`
2. Navigate to: **Windows → TREXIO Tools**

### Loading a TREXIO File

1. Click **"Open TREXIO File..."** button
2. Select a `.h5` or `.trexio` file from the file picker
3. The summary window will automatically display metadata

### Viewing Summary Information

The **TREXIO Summary** window displays:

- **System Information**:
  - Number of atoms
  - Number of electrons
  
- **Basis Set**:
  - Number of shells
  - Number of primitives
  - Number of atomic orbitals (AOs)
  
- **Molecular Orbitals**:
  - Number of molecular orbitals (MOs)
  
- **Atom List**:
  - First 50 atoms with element labels and coordinates (x, y, z in Angstrom)

### Computing Orbital Grids

1. **Open the Orbital Grid window**:
   - In the Summary window, click **"Open Orbital Grid Tool"**
   - Or click **"Show Orbital Grid"** from the main TREXIO Tools window

2. **Select Molecular Orbital**:
   - Use the **"MO Index"** slider to select the orbital (0 to N-1)
   - For HOMO, use the index corresponding to the highest occupied orbital
   - For LUMO, use the next index after HOMO

3. **Configure Grid Parameters**:
   
   **Bounding Box** (in Angstrom):
   - X Range: `[min_x, max_x]`
   - Y Range: `[min_y, max_y]`
   - Z Range: `[min_z, max_z]`
   - Default: Auto-suggested from molecule geometry with 5Å padding
   
   **Resolution**:
   - Res X, Y, Z: Number of grid points per dimension (8-128)
   - Default: 32 × 32 × 32
   - Higher resolution = better detail but longer computation time

4. **Compute the Grid**:
   - Click **"Compute Grid"**
   - Watch the progress bar (0-100%)
   - Monitor log messages in the log window
   - **To cancel**: Click "Cancel" button during computation

5. **View Results**:
   
   After successful computation, the **Grid Statistics** section displays:
   - **Min value**: Minimum orbital value on the grid
   - **Max value**: Maximum orbital value on the grid
   - **Mean value**: Average orbital value across all grid points
   - **Output file**: Path to exported grid file

### Verbose Logging

Enable **"Verbose Logging"** in the main TREXIO Tools window to:
- Write detailed logs to `tools/trexio/debug-YYYYMMDD-HHMMSS.log`
- Track allocator operations
- Debug file loading and grid computation

**Note**: The debug log directory (`tools/trexio/`) must exist before enabling verbose logging.

## Grid File Format

Computed grids are saved in a custom binary format (`.grid` files):

### File Structure

```
Header (48 bytes):
  - Magic number: 0x4754524F ("GTRO") - 4 bytes
  - Version: 1 - 4 bytes
  - Resolution: nx, ny, nz - 12 bytes (3 × int32)
  - Bounding box: min_x, max_x, min_y, max_y, min_z, max_z - 24 bytes (6 × float)

Grid Data:
  - Values: nx × ny × nz floats (4 bytes each)
  - Stored in row-major order (z varies fastest, then y, then x)
```

### Reading Grid Files (Python Example)

```python
import numpy as np
import struct

def read_grid_file(filename):
    with open(filename, 'rb') as f:
        # Read header
        magic = struct.unpack('I', f.read(4))[0]
        assert magic == 0x4754524F, "Invalid magic number"
        
        version = struct.unpack('I', f.read(4))[0]
        assert version == 1, "Unsupported version"
        
        nx, ny, nz = struct.unpack('iii', f.read(12))
        bbox = struct.unpack('ffffff', f.read(24))
        min_x, max_x, min_y, max_y, min_z, max_z = bbox
        
        # Read grid data
        grid_size = nx * ny * nz
        grid_data = np.fromfile(f, dtype=np.float32, count=grid_size)
        grid = grid_data.reshape((nx, ny, nz))
        
    return {
        'grid': grid,
        'resolution': (nx, ny, nz),
        'bounds': {'x': (min_x, max_x), 'y': (min_y, max_y), 'z': (min_z, max_z)}
    }
```

## Example Workflows

### Computing HOMO for Water (H2O)

1. Load file: `test_data/h2o_molecule.trexio`
2. Summary shows: 3 atoms, 10 electrons, 24 MOs
3. HOMO index = 4 (5th MO, 0-indexed)
4. Open Orbital Grid tool
5. Set MO Index = 4
6. Grid parameters: Default 32³ with auto-suggested bounds
7. Click "Compute Grid"
8. Wait ~5-30 seconds (depending on CPU)
9. View statistics and grid file path

### High-Resolution Grid for Publication

1. Load your TREXIO file
2. Select desired MO
3. Increase resolution: 64 × 64 × 64 or 128 × 128 × 128
4. Adjust bounding box to focus on region of interest
5. Compute grid (may take several minutes)
6. Export grid file for rendering in external tools

## Performance Tips

- **Start with low resolution** (32³) for quick preview
- **Increase resolution** (64³, 128³) only for final production
- **Adjust bounding box** to reduce computation volume
- **Use verbose logging** to diagnose slow operations
- **Cancel and restart** if parameters need adjustment

## Troubleshooting

### "Failed to load TREXIO file"

- **Check file format**: Ensure file is valid TREXIO (.h5 or .trexio)
- **Verify TREXIO library**: Run `pkg-config --modversion trexio`
- **Check file permissions**: Ensure file is readable

### "No GTOs available for MO evaluation"

- **Missing basis set data**: TREXIO file must contain basis set information
- **Missing MO coefficients**: File must have molecular orbital data
- **Regenerate file**: Use PySCF or other quantum chemistry software with TREXIO export

### Grid computation is very slow

- **Reduce resolution**: Try 16³ or 32³ first
- **Reduce bounding box**: Focus on smaller region
- **Check system load**: Close other applications
- **Use smaller molecule**: Large basis sets take longer

### Debug log not created

- **Create directory**: `mkdir -p tools/trexio`
- **Check permissions**: Ensure write access to tools/trexio/
- **Enable logging first**: Toggle verbose logging before loading file

## Developer Notes

### Extending the UI

The TREXIO panel is event-driven and follows viamd's component architecture:

- **Component**: `src/components/trexio/trexio.cpp` - Event handler registration
- **UI Panel**: `src/ui/trexio_panel.cpp` - ImGui interface
- **Grid Evaluator**: `src/trexio/trexio_orbital_grid.cpp` - Computation engine

### Adding New Features

1. **UI Changes**: Edit `src/ui/trexio_panel.cpp`
2. **Computation**: Modify `src/trexio/trexio_orbital_grid.cpp`
3. **Integration**: Update `src/components/trexio/trexio.cpp`
4. **Rebuild**: `make -j4` in build directory

### Testing

See `test_data/TESTING_GUIDE.md` for:
- Generating test TREXIO files
- Validating grid computations
- Running headless tests

## References

- **TREXIO Library**: https://github.com/TREX-CoE/trexio
- **TREXIO Documentation**: https://trex-coe.github.io/trexio/
- **viamd Repository**: https://github.com/scanberg/viamd
- **GTO Evaluation**: See `ext/mdlib/src/md_gto.h`

## License

This TREXIO integration follows viamd's license terms. TREXIO library is BSD-licensed.

## Support

For issues or questions:
- Open an issue on the viamd GitHub repository
- Check TREXIO_GUI_IMPLEMENTATION.md for technical details
- Review test_data/ examples for reference
