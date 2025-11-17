# VIAMD
Visual Interactive Analysis of Molecular Dynamics

VIAMD is an interactive analysis tool for molecular dynamics (MD) written in C/C++. VIAMD is developed at the PDC Center for High Performance Computing (KTH, Stockholm). It exposes a rudementary script language that is used to declare operations which are performed over the frames of the trajectory.
The results can then be viewed in the different windows exposed in the application. 
<p align="center">
<img src="https://github.com/scanberg/viamd/assets/38646069/5651ef62-28bc-4f41-8234-75cf9ba85612" alt="This is an overview of the viamd software" width="800"/>
</p>

## Status
[![Windows (MSVC 19)](https://github.com/scanberg/viamd/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/scanberg/viamd/actions/workflows/windows.yml)
[![Ubuntu 22.04 (GCC 11)](https://github.com/scanberg/viamd/actions/workflows/ubuntu22.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu22.yml)
[![Ubuntu 24.04 (GCC 13)](https://github.com/scanberg/viamd/actions/workflows/ubuntu24.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu24.yml)
[![MacOS (Clang)](https://github.com/scanberg/viamd/actions/workflows/macos.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/macos.yml)

## Running VIAMD 

### Windows
For windows, we recommend to use the latest binary available on the [release page](https://github.com/scanberg/viamd/releases/).

### Ubuntu and MacOs
To [build](https://github.com/scanberg/viamd/wiki/0.-Building) VIAMD on your machine, you can follow the procedure described in details in the wiki for [Linux](https://github.com/scanberg/viamd/wiki/0.-Building#linux) and [MacOS](https://github.com/scanberg/viamd/wiki/0.-Building#mac).

## Building with Optional Features

### TREXIO Support

VIAMD can be built with support for reading TREXIO quantum chemistry files. TREXIO is an open-source file format used by many quantum chemistry codes (Quantum Package, PySCF, FHI-aims, CP2K, etc.).

TREXIO is automatically downloaded and built from source during the VIAMD build process. No manual installation is required.

#### Prerequisites
- HDF5 library (optional, for HDF5 backend support)

#### Installing HDF5 (Optional but Recommended)

HDF5 support enables TREXIO to read .h5 files. If HDF5 is not available, TREXIO will build with text backend support only.

```bash
# Ubuntu/Debian
sudo apt-get install libhdf5-dev

# macOS
brew install hdf5

# Conda
conda install -c conda-forge hdf5
```

#### Building VIAMD with TREXIO

```bash
# Initialize submodules
git submodule update --init --recursive

# Apply mdlib patch to enable TREXIO support
cd ext/mdlib
git apply ../../docs/mdlib_trexio.patch
cd ../..

# Configure and build
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

The build system will automatically:
1. Download the TREXIO 2.6.0 release tarball
2. Detect if HDF5 is available on your system
3. Build TREXIO with HDF5 support (if available) or text-only backend
4. Link TREXIO statically into VIAMD

**CMake Options:**
- `-DVIAMD_ENABLE_TREXIO=ON` - Enable TREXIO file format support (default: OFF)
- `-DVIAMD_ENABLE_VELOXCHEM=ON` - Enable VeloxChem module (default: OFF)

See `docs/TREXIO_SUPPORT.md` for detailed documentation on TREXIO support.

### VeloxChem Support

VIAMD can also be built with VeloxChem support for quantum chemistry calculations:

```bash
cmake -DVIAMD_ENABLE_VELOXCHEM=ON ..
make
```

## Documentation
Documentation about VIAMD is available on the github [wiki](https://github.com/scanberg/viamd/wiki). The two first chapters relate to the [visual](https://github.com/scanberg/viamd/wiki/1.-Visual) and [analysis](https://github.com/scanberg/viamd/wiki/2.-Analysis) features respectively, where we highlight the interactive part of software. The third chapter focus on the VIAMD [language](https://github.com/scanberg/viamd/wiki/3.-Language) used for scripting and the fourth chapter propose a serie of [tutorial](https://github.com/scanberg/viamd/wiki/4.-Tutorials) (under construction). 

A series of videos is available on [youtube](https://youtube.com/playlist?list=PLNx9MpJY8ffr9CeK7WefdOnuGRw_E5rSj&si=VatBHEwiL7jWyhPK).

## Update
If you want to stay informed about the latest update of VIAMD, please register your email address to the [form](https://forms.gle/fAxuWob8nMLcrS5h9). 

## Citations:
* General Framework:
  * R Skånberg, I Hotz, A Ynnerman, M Linares, VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics, J. Chem. Inf. Model. 2023, 63, 23, 7382–7391 https://doi.org/10.1021/acs.jcim.3c01033
  * R Skånberg, C König, P Norman, M Linares, D Jönsson, I Hotz, A Ynnerman, VIA-MD: Visual Interactive Analysis of Molecular Dynamics, 2018, Eurographics Proceedings, p. 19–27

* Specific tool:
  * Selection tool: Robin Skånberg, Mathieu Linares, Martin Falk, Ingrid Hotz, Anders Ynnerman, MolFind-Integrated Multi-Selection Schemes for Complex Molecular Structures, 2019, The Eurographics Association, p. 17-21​
  * Shape Space and Spatial Distribution Function: Robin Skånberg, Martin Falk, Mathieu Linares, Anders Ynnerman, Ingrid Hotz, Tracking Internal Frames of Reference for Consistent Molecular Distribution Functions, 2021, IEEE Transactions on Visualization and Computer Graphics, 28 (9), 3126-3137​

## Financial Support
VIAMD has received constant financial support since 2018 from the Swedish e-Research center ([SeRC](https://e-science.se/)) and the [Wallenberg Foundation](https://www.wallenberg.org/en)

VIAMD is supported by [InfraVis](https://infravis.se/) for specific projets:
- Parser for LAMMPS file (2301-5217 / 140 hours)
- Interactive analysis of [VeloxChem](https://veloxchem.org/docs/intro.html) file (interactive analysis of orbitals and spectra plotting) (600 hours)
- Support for [TREXIO](https://github.com/TREX-CoE/trexio) quantum chemistry file format 

<p align="center">
<img src="https://github.com/scanberg/viamd/assets/38646069/e7245119-3ec4-4b84-9056-7197b3d1448b"  height="75" >
<img src="https://github.com/scanberg/viamd/assets/38646069/f1c8493f-9519-4458-87c6-2d57a4071ad7"  height="75" >
<img src="https://github.com/scanberg/viamd/assets/38646069/cfc3feed-728f-45c2-a7db-c3c0707acbb1"  height="75" >
</p>

## Acknowledgements

https://github.com/glfw/glfw

https://github.com/dougbinks/enkiTS

https://github.com/ocornut/imgui

https://github.com/epezent/implot

https://github.com/BalazsJako/ImGuiColorTextEdit

https://github.com/skaslev/gl3w

https://github.com/max0x7ba/atomic_queue

https://github.com/mlabbe/nativefiledialog

https://github.com/nothings/stb

#
<p align="center">
<img src="https://github.com/user-attachments/assets/39b69b10-88a1-43a7-9d69-68513ac4e632"  width="150" alt="This is the VIAMD logo" >
</p>



