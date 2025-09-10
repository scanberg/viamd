# OpenMM Integration in VIAMD

This document describes the OpenMM molecular dynamics simulation integration in VIAMD.

## Overview

VIAMD now includes integrated molecular dynamics simulation capabilities through OpenMM's C++ API. This allows users to run real-time MD simulations directly within VIAMD while leveraging its powerful visualization and analysis tools.

## Features

### Simulation Capabilities
- **Real-time MD simulation** with live coordinate updates
- **Langevin dynamics integration** with temperature control
- **Multiple force field support** including AMBER and UFF (Universal Force Field)
- **Harmonic bonds and angles** with force field-specific parameters
- **Non-bonded interactions** with Lennard-Jones potentials
- **Customizable parameters** (temperature, timestep, friction)
- **Performance monitoring** and simulation statistics

### Integration Features
- **Seamless integration** with VIAMD's existing visualization pipeline
- **Event-driven architecture** using VIAMD's component system
- **Minimal impact** on existing VIAMD functionality
- **Optional compilation** - can be disabled if OpenMM is not available

## Installation

### Prerequisites
- OpenMM library and development headers
- C++20 compatible compiler
- CMake 3.20 or higher

### Ubuntu/Debian Installation
```bash
sudo apt install libopenmm-dev libopenmm8.0t64 libopenmm-plugins
```

### Building VIAMD with OpenMM
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DVIAMD_ENABLE_OPENMM=ON
make -j$(nproc)
```

## Usage

### 1. Load a Molecular Structure
- Use VIAMD's normal file loading mechanisms to load a PDB, XYZ, or other supported molecular structure file
- The molecule should have proper atom types for force field assignment

### 2. Open the Simulation Window
- Go to **Windows → OpenMM Simulation** in the main menu
- This opens the simulation control panel

### 3. Select Force Field
- Choose between **AMBER14** and **UFF** (Universal Force Field) using the dropdown
- **AMBER14**: Traditional biomolecular force field, optimized for proteins and nucleic acids
- **UFF**: Universal force field, works well for organic molecules and diverse chemical systems
- Force field can only be changed before system initialization

### 4. Initialize the Simulation System
- Click **"Initialize System"** to set up the OpenMM simulation
- VIAMD will automatically:
  - Create particles for each atom with appropriate masses
  - Set up harmonic bond forces based on molecular connectivity and selected force field
  - Configure angle forces using force field-specific parameters
  - Set up non-bonded interactions with force field-appropriate parameters
  - Initialize the Langevin integrator

### 5. Configure Simulation Parameters
- **Temperature**: Set the simulation temperature (default: 300 K)
- **Timestep**: Adjust the integration timestep (default: 0.002 ps)
- **Friction**: Control the friction coefficient (default: 1.0 ps⁻¹)
- **Steps per update**: Number of MD steps between visualization updates (default: 10)

### 6. Run the Simulation
- Click **"Start Simulation"** to begin the MD simulation
- Watch molecules move in real-time with full VIAMD visualization
- Use **"Pause"** and **"Resume"** to control simulation flow
- Click **"Stop"** to halt the simulation
- Use **"Reset"** to return to initial coordinates

## Technical Details

### Force Field Implementation
VIAMD supports two force fields for molecular dynamics simulations:

#### AMBER14 Force Field
- **Traditional biomolecular force field** optimized for proteins, nucleic acids, and biological systems
- **Atom types**: Based on chemical environment and hybridization (CA, C, N, O, H, etc.)
- **Parameters**: Derived from extensive experimental data and quantum chemical calculations
- **Best for**: Proteins, DNA/RNA, carbohydrates, and other biomolecules

#### UFF (Universal Force Field)
- **General-purpose force field** suitable for most of the periodic table
- **Atom types**: Based on element, hybridization, and coordination (C_3, C_2, N_3, O_2, etc.)
- **Parameters**: Derived from simple rules and periodic trends
- **Best for**: Organic molecules, inorganic compounds, and diverse chemical systems

Both force fields include:
- **Harmonic bonds**: Applied to all bonds detected in the molecular structure
- **Harmonic angles**: Applied to all three-atom angle terms
- **Non-bonded interactions**: Lennard-Jones potential with force field-specific parameters
- **Electrostatic interactions**: Coulomb potential with partial charges (reduced for stability)

### Coordinate System
- OpenMM uses nanometers; VIAMD uses Angstroms
- Automatic conversion is handled internally
- Real-time coordinate updates are applied to VIAMD's molecular representation

### Performance Considerations
- Simulation performance depends on system size and update frequency
- Reduce "Steps per update" for smoother visualization of fast dynamics
- Increase "Steps per update" for better simulation performance with large systems

### Memory Management
- Uses VIAMD's allocator system for consistent memory management
- OpenMM contexts are properly cleaned up when simulations end
- Automatic cleanup when new molecular structures are loaded

## Example Workflow

1. **Load the test molecule**:
   ```
   Load datasets/1ALA-500.pdb
   ```

2. **Set up visualization**:
   - Choose appropriate representations (Ball & Stick, etc.)
   - Adjust camera view

3. **Initialize simulation**:
   - Open Windows → OpenMM Simulation
   - Click "Initialize System"
   - Adjust temperature to 350 K for more visible motion

4. **Run simulation**:
   - Click "Start Simulation"
   - Watch the alanine dipeptide move in real-time
   - Observe conformational changes and thermal motion

## Troubleshooting

### Common Issues

**"OpenMM not found" during build**:
- Ensure OpenMM development packages are installed
- Check that CMake can find OpenMM headers and libraries

**"No atoms available for simulation setup"**:
- Load a molecular structure file first
- Ensure the structure contains valid atom information

**Simulation fails to start**:
- Check that the molecular structure has reasonable coordinates
- Verify that bonds were properly detected
- Review console output for detailed error messages

**Poor simulation performance**:
- Reduce system size or increase "Steps per update"
- Consider using simpler molecular systems for real-time visualization

## Limitations

### Current Limitations
- **Simplified force field implementation**: Not suitable for quantitative research, but good for visualization and education
- **No periodic boundary conditions**: Systems are simulated in vacuum
- **Limited atom type recognition**: Basic assignment based on connectivity and element type
- **No advanced sampling methods**: Only basic Langevin dynamics

### Future Enhancements
- **Enhanced force field parameter assignment** with better chemical environment recognition
- **Periodic boundary conditions** support for bulk and interfacial systems
- **Additional force fields** (CHARMM, GROMOS, etc.)
- **Advanced simulation methods** (replica exchange, enhanced sampling)
- **Trajectory saving and analysis tools**

## API Integration

### Component Architecture
The OpenMM integration follows VIAMD's component pattern:
- Located in `src/components/openmm/openmm.cpp`
- Uses VIAMD's event system for communication
- Integrates with ApplicationState for state management

### Event Handling
The component responds to:
- `ViamdInitialize`: Component initialization
- `ViamdFrameTick`: Simulation updates
- `ViamdTopologyInit`: New molecule loaded
- `ViamdWindowDrawMenu`: UI menu integration

### Extension Points
Developers can extend the OpenMM integration by:
- **Adding new force field types** (CHARMM, GROMOS, OPLS, etc.)
- **Implementing custom integrators** for specialized sampling methods
- **Adding analysis tools** for energy, structure, and dynamics analysis
- **Enhancing the user interface** with advanced simulation controls

## Conclusion

The OpenMM integration provides VIAMD with professional-grade molecular dynamics simulation capabilities while maintaining its strength in visualization and analysis. With support for both AMBER and UFF force fields, users can simulate a wide range of molecular systems from biological macromolecules to small organic compounds. This combination makes VIAMD a powerful tool for both education and research in molecular simulation.