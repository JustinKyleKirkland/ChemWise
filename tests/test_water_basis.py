#!/usr/bin/env python3
"""
Test water molecule with different basis sets.
"""

import numpy as np
from chemwise.core.molecule import Molecule, Atom
from chemwise.core.basis import BasisSet, list_available_basis_sets
from chemwise.methods.hf import HartreeFock

def test_water():
    """Test basis sets with water molecule."""
    
    print("=== Water Molecule Basis Set Test ===\n")
    
    # Create water molecule (H2O)
    o_atom = Atom('O', np.array([0.0, 0.0, 0.0]))
    h1_atom = Atom('H', np.array([0.0, 0.757, 0.587])) 
    h2_atom = Atom('H', np.array([0.0, -0.757, 0.587]))
    water = Molecule([o_atom, h1_atom, h2_atom])
    
    print(f"Testing molecule: {water}")
    print("Geometry:")
    for i, atom in enumerate(water.atoms):
        print(f"  {atom.symbol}: {atom.coordinates}")
    print()
    
    # Test with STO-3G first (should be stable)
    basis_sets_to_test = ['sto-3g', '6-31g']  # Start with simpler ones
    
    for basis_name in basis_sets_to_test:
        try:
            print(f"--- Testing {basis_name.upper()} basis set ---")
            
            # Create basis set
            basis = BasisSet.from_molecule(water, basis_name)
            print(f"Basis set: {basis}")
            print(f"Number of basis functions: {len(basis)}")
            
            # Run HF calculation
            hf = HartreeFock(water, basis)
            hf_energy = hf.run()
            
            print(f"HF Energy: {hf_energy:.8f} Hartree")
            print()
            
        except Exception as e:
            print(f"Error with {basis_name}: {e}")
            print()

if __name__ == "__main__":
    test_water()
