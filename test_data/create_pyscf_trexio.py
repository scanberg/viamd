#!/usr/bin/env python3
"""
Create TREXIO test files using PySCF (if available).
Falls back to manual text format if PySCF/TREXIO are not installed.

This script demonstrates how to create TREXIO files from quantum chemistry
calculations, useful for Phase 4 testing.
"""

import os
import sys

def check_dependencies():
    """Check if PySCF and TREXIO are available."""
    has_pyscf = False
    has_trexio = False
    
    try:
        import pyscf
        has_pyscf = True
        print(f"✓ PySCF {pyscf.__version__} found")
    except ImportError:
        print("✗ PySCF not found (pip install pyscf)")
    
    try:
        import trexio
        has_trexio = True
        print(f"✓ TREXIO Python bindings found")
    except ImportError:
        print("✗ TREXIO Python bindings not found (pip install trexio)")
    
    return has_pyscf, has_trexio

def create_h2_with_pyscf(output_file='h2_pyscf.h5'):
    """Create H2 molecule TREXIO file using PySCF."""
    try:
        from pyscf import gto, scf
        import trexio
        import numpy as np
        
        print(f"\nCreating {output_file} with PySCF...")
        
        # Define H2 molecule
        mol = gto.M(
            atom='H 0 0 0; H 0 0 0.74',
            basis='6-31g',
            unit='Angstrom'
        )
        
        # Run SCF calculation
        mf = scf.RHF(mol)
        mf.kernel()
        
        # Write to TREXIO HDF5 format
        with trexio.File(output_file, 'w', trexio.TREXIO_HDF5) as f:
            # Nucleus data
            trexio.write_nucleus_num(f, mol.natm)
            trexio.write_nucleus_charge(f, mol.atom_charges())
            trexio.write_nucleus_coord(f, mol.atom_coords())
            
            # Convert atom labels
            labels = [mol.atom_symbol(i) for i in range(mol.natm)]
            # TREXIO expects fixed-length strings, pad to 4 chars
            labels_padded = [label.ljust(4) for label in labels]
            trexio.write_nucleus_label(f, labels_padded)
            
            # Electron data
            nelec = mol.nelectron
            trexio.write_electron_up_num(f, nelec // 2)
            trexio.write_electron_dn_num(f, nelec - nelec // 2)
            
            # Basis set data (simplified - full basis would be more complex)
            # For testing, we'll write minimal basis info
            
            # AO/MO data
            nao = mf.mo_coeff.shape[0]
            nmo = mf.mo_coeff.shape[1]
            
            trexio.write_ao_num(f, nao)
            trexio.write_mo_num(f, nmo)
            trexio.write_mo_coefficient(f, mf.mo_coeff)
            trexio.write_mo_energy(f, mf.mo_energy)
            
            # MO occupation
            mo_occ = np.zeros(nmo)
            mo_occ[:nelec//2] = 2.0  # Closed shell
            trexio.write_mo_occupation(f, mo_occ)
        
        print(f"✓ Created {output_file}")
        print(f"  Atoms: {mol.natm}")
        print(f"  Electrons: {mol.nelectron}")
        print(f"  Basis functions: {nao}")
        print(f"  MOs: {nmo}")
        print(f"  SCF energy: {mf.e_tot:.6f} Hartree")
        
        return True
        
    except Exception as e:
        print(f"✗ Failed to create {output_file}: {e}")
        return False

def create_h2o_with_pyscf(output_file='h2o_pyscf.h5'):
    """Create H2O molecule TREXIO file using PySCF."""
    try:
        from pyscf import gto, scf
        import trexio
        import numpy as np
        
        print(f"\nCreating {output_file} with PySCF...")
        
        # Define water molecule
        mol = gto.M(
            atom='''
            O  0.0000  0.0000  0.1173
            H  0.0000  0.7572 -0.4692
            H  0.0000 -0.7572 -0.4692
            ''',
            basis='sto-3g',
            unit='Angstrom'
        )
        
        # Run SCF calculation
        mf = scf.RHF(mol)
        mf.kernel()
        
        # Write to TREXIO
        with trexio.File(output_file, 'w', trexio.TREXIO_HDF5) as f:
            # Nucleus data
            trexio.write_nucleus_num(f, mol.natm)
            trexio.write_nucleus_charge(f, mol.atom_charges())
            trexio.write_nucleus_coord(f, mol.atom_coords())
            
            labels = [mol.atom_symbol(i).ljust(4) for i in range(mol.natm)]
            trexio.write_nucleus_label(f, labels)
            
            # Electron data
            nelec = mol.nelectron
            trexio.write_electron_up_num(f, nelec // 2)
            trexio.write_electron_dn_num(f, nelec - nelec // 2)
            
            # AO/MO data
            nao = mf.mo_coeff.shape[0]
            nmo = mf.mo_coeff.shape[1]
            
            trexio.write_ao_num(f, nao)
            trexio.write_mo_num(f, nmo)
            trexio.write_mo_coefficient(f, mf.mo_coeff)
            trexio.write_mo_energy(f, mf.mo_energy)
            
            mo_occ = np.zeros(nmo)
            mo_occ[:nelec//2] = 2.0
            trexio.write_mo_occupation(f, mo_occ)
        
        print(f"✓ Created {output_file}")
        print(f"  Atoms: {mol.natm}")
        print(f"  Electrons: {mol.nelectron}")
        print(f"  Basis functions: {nao}")
        print(f"  SCF energy: {mf.e_tot:.6f} Hartree")
        
        return True
        
    except Exception as e:
        print(f"✗ Failed to create {output_file}: {e}")
        return False

def main():
    print("=" * 70)
    print("TREXIO Test File Generator (PySCF)")
    print("=" * 70)
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Check dependencies
    print("\nChecking dependencies...")
    has_pyscf, has_trexio = check_dependencies()
    
    if not (has_pyscf and has_trexio):
        print("\n" + "=" * 70)
        print("INSTALL INSTRUCTIONS")
        print("=" * 70)
        print("\nTo generate TREXIO files with PySCF, install:")
        print("  pip install pyscf")
        print("  pip install trexio  # or: conda install -c conda-forge trexio")
        print("\nAlternatively, use the basic test files created by:")
        print("  python3 create_test_trexio.py")
        print("=" * 70)
        return 1
    
    # Create test files
    print("\n" + "=" * 70)
    print("Creating TREXIO test files with quantum chemistry data...")
    print("=" * 70)
    
    success = []
    
    # H2 molecule
    if create_h2_with_pyscf('h2_pyscf.h5'):
        success.append('h2_pyscf.h5')
    
    # H2O molecule
    if create_h2o_with_pyscf('h2o_pyscf.h5'):
        success.append('h2o_pyscf.h5')
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nCreated {len(success)} TREXIO files:")
    for f in success:
        print(f"  ✓ {f}")
    
    if success:
        print("\nThese files can be used to test VIAMD's TREXIO support:")
        print("  ./viamd h2_pyscf.h5")
        print("  ./viamd h2o_pyscf.h5")
        print("\nOr run unit tests with:")
        print("  cd build && ctest -V -R trexio")
    
    print("=" * 70)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
