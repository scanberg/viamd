#!/usr/bin/env python3
"""
Create sample TREXIO test files for VIAMD testing.
This script creates simple test molecules in TREXIO text format.
"""

import os
import sys

def create_trexio_text_h2():
    """Create a simple H2 molecule in TREXIO text format."""
    
    # Create directory structure for TREXIO text format
    trexio_dir = "h2_molecule.trexio"
    os.makedirs(trexio_dir, exist_ok=True)
    
    # Write nucleus group
    with open(os.path.join(trexio_dir, "nucleus.txt"), "w") as f:
        f.write("num\n")
        f.write("2\n")
        f.write("\n")
        f.write("charge\n")
        f.write("1.0\n")
        f.write("1.0\n")
        f.write("\n")
        f.write("coord\n")
        f.write("0.0 0.0 0.0\n")
        f.write("0.0 0.0 1.4\n")  # ~0.74 Angstrom in Bohr
        f.write("\n")
        f.write("label\n")
        f.write("H\n")
        f.write("H\n")
        f.write("\n")
        f.write("repulsion\n")
        f.write("0.7142857143\n")  # 1/1.4 in atomic units
        f.write("\n")
    
    # Write electron group
    with open(os.path.join(trexio_dir, "electron.txt"), "w") as f:
        f.write("up_num\n")
        f.write("1\n")
        f.write("\n")
        f.write("dn_num\n")
        f.write("1\n")
        f.write("\n")
    
    # Write metadata
    with open(os.path.join(trexio_dir, "metadata.txt"), "w") as f:
        f.write("package_version\n")
        f.write("test_v1.0\n")
        f.write("\n")
        f.write("description\n")
        f.write("H2 test molecule for VIAMD\n")
        f.write("\n")
    
    print(f"Created TREXIO text format file: {trexio_dir}")
    return trexio_dir

def create_trexio_text_h2o():
    """Create a water molecule in TREXIO text format."""
    
    trexio_dir = "h2o_molecule.trexio"
    os.makedirs(trexio_dir, exist_ok=True)
    
    # Water molecule geometry (in Bohr)
    # O at origin, H atoms forming ~104.5 degree angle
    with open(os.path.join(trexio_dir, "nucleus.txt"), "w") as f:
        f.write("num\n")
        f.write("3\n")
        f.write("\n")
        f.write("charge\n")
        f.write("8.0\n")   # Oxygen
        f.write("1.0\n")   # Hydrogen
        f.write("1.0\n")   # Hydrogen
        f.write("\n")
        f.write("coord\n")
        f.write("0.0 0.0 0.0\n")                    # O
        f.write("1.43 1.11 0.0\n")                  # H1
        f.write("-1.43 1.11 0.0\n")                 # H2
        f.write("\n")
        f.write("label\n")
        f.write("O\n")
        f.write("H\n")
        f.write("H\n")
        f.write("\n")
        f.write("repulsion\n")
        f.write("9.195908\n")  # Approximate nuclear repulsion energy
        f.write("\n")
    
    # Write electron group
    with open(os.path.join(trexio_dir, "electron.txt"), "w") as f:
        f.write("up_num\n")
        f.write("5\n")  # 10 electrons total, 5 up
        f.write("\n")
        f.write("dn_num\n")
        f.write("5\n")  # 5 down
        f.write("\n")
    
    # Write metadata
    with open(os.path.join(trexio_dir, "metadata.txt"), "w") as f:
        f.write("package_version\n")
        f.write("test_v1.0\n")
        f.write("\n")
        f.write("description\n")
        f.write("H2O water molecule for VIAMD testing\n")
        f.write("\n")
    
    print(f"Created TREXIO text format file: {trexio_dir}")
    return trexio_dir

def create_trexio_text_ch4():
    """Create a methane molecule in TREXIO text format."""
    
    trexio_dir = "ch4_molecule.trexio"
    os.makedirs(trexio_dir, exist_ok=True)
    
    # Methane tetrahedral geometry (in Bohr)
    with open(os.path.join(trexio_dir, "nucleus.txt"), "w") as f:
        f.write("num\n")
        f.write("5\n")
        f.write("\n")
        f.write("charge\n")
        f.write("6.0\n")   # Carbon
        f.write("1.0\n")   # Hydrogen
        f.write("1.0\n")
        f.write("1.0\n")
        f.write("1.0\n")
        f.write("\n")
        f.write("coord\n")
        f.write("0.0 0.0 0.0\n")                    # C at origin
        f.write("1.18 1.18 1.18\n")                 # H1
        f.write("-1.18 -1.18 1.18\n")               # H2
        f.write("-1.18 1.18 -1.18\n")               # H3
        f.write("1.18 -1.18 -1.18\n")               # H4
        f.write("\n")
        f.write("label\n")
        f.write("C\n")
        f.write("H\n")
        f.write("H\n")
        f.write("H\n")
        f.write("H\n")
        f.write("\n")
        f.write("repulsion\n")
        f.write("13.486\n")  # Approximate nuclear repulsion energy
        f.write("\n")
    
    # Write electron group
    with open(os.path.join(trexio_dir, "electron.txt"), "w") as f:
        f.write("up_num\n")
        f.write("5\n")  # 10 electrons total
        f.write("\n")
        f.write("dn_num\n")
        f.write("5\n")
        f.write("\n")
    
    # Write metadata
    with open(os.path.join(trexio_dir, "metadata.txt"), "w") as f:
        f.write("package_version\n")
        f.write("test_v1.0\n")
        f.write("\n")
        f.write("description\n")
        f.write("CH4 methane molecule for VIAMD testing\n")
        f.write("\n")
    
    print(f"Created TREXIO text format file: {trexio_dir}")
    return trexio_dir

def create_readme():
    """Create README for test data."""
    
    with open("README.md", "w") as f:
        f.write("# TREXIO Test Data\n\n")
        f.write("This directory contains test TREXIO files for validating VIAMD's TREXIO support.\n\n")
        f.write("## Test Files\n\n")
        f.write("### Text Format (.trexio directories)\n\n")
        f.write("- **h2_molecule.trexio** - Simple H2 molecule (2 atoms)\n")
        f.write("- **h2o_molecule.trexio** - Water molecule (3 atoms)\n")
        f.write("- **ch4_molecule.trexio** - Methane molecule (5 atoms)\n\n")
        f.write("### File Format\n\n")
        f.write("TREXIO text format uses a directory structure with .txt files for each data group.\n")
        f.write("This format is human-readable and doesn't require HDF5 library.\n\n")
        f.write("## Usage\n\n")
        f.write("To test VIAMD with these files:\n\n")
        f.write("```bash\n")
        f.write("# Build VIAMD with TREXIO support\n")
        f.write("cd /home/runner/work/viamd/viamd\n")
        f.write("cd ext/mdlib && git apply ../../docs/mdlib_trexio.patch\n")
        f.write("cd ../..\n")
        f.write("mkdir build && cd build\n")
        f.write("cmake -DVIAMD_ENABLE_TREXIO=ON ..\n")
        f.write("make\n\n")
        f.write("# Load a test file\n")
        f.write("./viamd ../test_data/h2_molecule.trexio\n")
        f.write("```\n\n")
        f.write("## Data Validation\n\n")
        f.write("Each test file contains:\n")
        f.write("- Nucleus data (coordinates in Bohr, charges, labels)\n")
        f.write("- Electron configuration (up/down spin numbers)\n")
        f.write("- Metadata (package version, description)\n\n")
        f.write("Coordinates should be converted from Bohr to Angstrom when loaded in VIAMD.\n")
    
    print("Created README.md")

if __name__ == "__main__":
    print("Creating TREXIO test files...")
    print("=" * 60)
    
    # Change to test_data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Create test files
    h2_dir = create_trexio_text_h2()
    h2o_dir = create_trexio_text_h2o()
    ch4_dir = create_trexio_text_ch4()
    
    # Create README
    create_readme()
    
    print("=" * 60)
    print("\nTest files created successfully!")
    print(f"\nLocation: {script_dir}")
    print("\nFiles created:")
    print(f"  - {h2_dir}/")
    print(f"  - {h2o_dir}/")
    print(f"  - {ch4_dir}/")
    print(f"  - README.md")
