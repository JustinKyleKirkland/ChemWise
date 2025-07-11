#!/usr/bin/env python3
"""
Simple test for the integral fix.
"""

import numpy as np
from chemwise.core.molecule import Molecule, Atom
from chemwise.core.basis import BasisSet
from chemwise.methods.hf import HartreeFock

def test_simple_molecules():
    """Test simple molecules with basis sets."""
    
    print("=== Simple Molecule Tests ===\n")
    
    # Test 1: H2 with STO-3G
    print("1. H2 with STO-3G:")
    h1 = Atom('H', np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', np.array([0.0, 0.0, 1.4]))
    h2_mol = Molecule([h1, h2])
    
    try:
        basis = BasisSet.from_molecule(h2_mol, 'sto-3g')
        hf = HartreeFock(h2_mol, basis)
        energy = hf.run()
        print(f"  Energy: {energy:.8f} Hartree")
    except Exception as e:
        print(f"  Error: {e}")
    print()
    
    # Test 2: H2 with 6-31G
    print("2. H2 with 6-31G:")
    try:
        basis = BasisSet.from_molecule(h2_mol, '6-31g')
        hf = HartreeFock(h2_mol, basis)
        energy = hf.run()
        print(f"  Energy: {energy:.8f} Hartree")
    except Exception as e:
        print(f"  Error: {e}")
    print()
    
    # Test 3: HeH+ (simpler system with different elements)
    print("3. HeH+ with STO-3G:")
    he = Atom('He', np.array([0.0, 0.0, 0.0]))
    h = Atom('H', np.array([0.0, 0.0, 1.0]))
    heh_plus = Molecule([he, h], charge=1, multiplicity=1)
    
    try:
        basis = BasisSet.from_molecule(heh_plus, 'sto-3g')
        print(f"  Basis set: {basis}")
        hf = HartreeFock(heh_plus, basis)
        energy = hf.run()
        print(f"  Energy: {energy:.8f} Hartree")
    except Exception as e:
        print(f"  Error: {e}")
    print()

if __name__ == "__main__":
    test_simple_molecules()
