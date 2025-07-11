#!/usr/bin/env python3
"""
Test cc-pVDZ basis set specifically.
"""

import numpy as np
from chemwise.core.molecule import Molecule, Atom
from chemwise.core.basis import BasisSet
from chemwise.methods.hf import HartreeFock

def test_cc_pvdz():
    """Test cc-pVDZ basis set."""
    
    print("=== cc-pVDZ Basis Set Test ===\n")
    
    # Test with H2 molecule
    h1 = Atom('H', np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', np.array([0.0, 0.0, 1.4]))
    h2_mol = Molecule([h1, h2])
    
    print("Testing H2 with cc-pVDZ:")
    try:
        basis = BasisSet.from_molecule(h2_mol, 'cc-pvdz')
        print(f"  Basis set: {basis}")
        
        # Print basis function details
        for i, bf in enumerate(basis.basis_functions):
            print(f"    {i}: {bf}")
        
        print(f"  Number of basis functions: {len(basis)}")
        
        hf = HartreeFock(h2_mol, basis)
        energy = hf.run()
        print(f"  Energy: {energy:.8f} Hartree")
    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
    print()

if __name__ == "__main__":
    test_cc_pvdz()
