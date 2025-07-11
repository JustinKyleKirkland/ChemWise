#!/usr/bin/env python3
"""
Test script for new basis set functionality.
"""

from chemwise.core.molecule import Molecule, Atom
from chemwise.core.basis import BasisSet, list_available_basis_sets
from chemwise.methods.hf import HartreeFock

def test_basis_sets():
    """Test the new basis set functionality."""
    
    print("=== ChemWise Basis Set Test ===\n")
    
    # List available basis sets
    print("Available basis sets:")
    for basis in list_available_basis_sets():
        print(f"  - {basis}")
    print()
    
    # Create a simple H2 molecule
    import numpy as np
    h1 = Atom('H', np.array([0.0, 0.0, 0.0]))
    h2_atom = Atom('H', np.array([0.0, 0.0, 1.4]))
    h2 = Molecule([h1, h2_atom])
    
    print(f"Testing molecule: {h2}")
    print("Geometry:")
    for i, atom in enumerate(h2.atoms):
        print(f"  {atom.symbol}: {atom.coordinates}")
    print()
    
    # Test different basis sets
    basis_sets_to_test = ['sto-3g', 'cc-pvdz', 'def2-svp', '6-31g']
    
    results = {}
    
    for basis_name in basis_sets_to_test:
        try:
            print(f"--- Testing {basis_name.upper()} basis set ---")
            
            # Create basis set
            basis = BasisSet.from_molecule(h2, basis_name)
            print(f"Basis set: {basis}")
            print(f"Number of basis functions: {len(basis)}")
            
            # Print basis function details
            for i, bf in enumerate(basis.basis_functions[:5]):  # Show first 5
                print(f"  {i:2d}: {bf}")
            if len(basis) > 5:
                print(f"  ... and {len(basis) - 5} more")
            
            # Run HF calculation
            hf = HartreeFock(h2, basis)
            hf_energy = hf.run()
            
            results[basis_name] = {
                'n_basis': len(basis),
                'hf_energy': hf_energy
            }
            
            print(f"HF Energy: {hf_energy:.8f} Hartree")
            print()
            
        except Exception as e:
            print(f"Error with {basis_name}: {e}")
            print()
    
    # Summary
    print("=== Results Summary ===")
    print(f"{'Basis Set':<12} {'N_basis':<8} {'HF Energy (Hartree)':<18}")
    print("-" * 40)
    
    for basis_name, data in results.items():
        print(f"{basis_name.upper():<12} {data['n_basis']:<8} {data['hf_energy']:<18.8f}")
    
    print()
    print("Test completed successfully!")

if __name__ == "__main__":
    test_basis_sets()
