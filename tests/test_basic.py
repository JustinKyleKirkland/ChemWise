#!/usr/bin/env python3
"""
Basic test script for ChemWise functionality.
"""

import sys
import os
import numpy as np

# Add the project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import chemwise as cw


def test_h2_hf():
    """Test HF calculation on H2."""
    print("Testing H2 HF calculation...")
    
    # Create H2 molecule
    mol = cw.Molecule.from_string("""
    H  0.0  0.0  0.0
    H  0.0  0.0  0.741400
    """)
    
    # Run HF calculation
    calc = cw.Calculator(mol, method="hf", basis="sto-3g")
    result = calc.run(print_results=False)
    
    # Check results
    assert result.converged, "H2 HF calculation did not converge"
    assert abs(result.energy - (-46.4)) < 1.0, f"H2 HF energy unexpected: {result.energy}"
    
    print(f"✓ H2 HF energy: {result.energy:.6f} Hartree")
    return True


def test_h2_dft():
    """Test DFT calculation on H2."""
    print("Testing H2 B3LYP calculation...")
    
    # Create H2 molecule
    mol = cw.Molecule.from_string("""
    H  0.0  0.0  0.0
    H  0.0  0.0  0.741400
    """)
    
    # Run B3LYP calculation
    calc = cw.Calculator(mol, method="b3lyp", basis="sto-3g")
    result = calc.run(print_results=False)
    
    # Check results
    assert result.converged, "H2 B3LYP calculation did not converge"
    assert abs(result.energy - (-47.8)) < 1.0, f"H2 B3LYP energy unexpected: {result.energy}"
    
    print(f"✓ H2 B3LYP energy: {result.energy:.6f} Hartree")
    return True


def test_water_hf():
    """Test HF calculation on water."""
    print("Testing water HF calculation...")
    
    # Create water molecule
    mol = cw.Molecule.from_string("""
    O  0.0000000  0.0000000  0.1173000
    H  0.0000000  0.7572000 -0.4692000
    H  0.0000000 -0.7572000 -0.4692000
    """)
    
    # Run HF calculation with looser convergence for speed
    calc = cw.Calculator(mol, method="hf", basis="sto-3g", 
                        convergence=1e-6, max_iterations=20)
    result = calc.run(print_results=False)
    
    # Check results
    assert result.converged, "Water HF calculation did not converge"
    assert result.energy < 0, f"Water HF energy should be negative: {result.energy}"
    
    print(f"✓ Water HF energy: {result.energy:.6f} Hartree")
    return True


def test_molecule_properties():
    """Test molecular property calculations."""
    print("Testing molecular properties...")
    
    # Create water molecule
    mol = cw.Molecule.from_string("""
    O  0.0000000  0.0000000  0.1173000
    H  0.0000000  0.7572000 -0.4692000
    H  0.0000000 -0.7572000 -0.4692000
    """)
    
    # Check basic properties
    assert mol.n_atoms == 3, f"Expected 3 atoms, got {mol.n_atoms}"
    assert mol.n_electrons == 10, f"Expected 10 electrons, got {mol.n_electrons}"
    assert mol.nuclear_repulsion_energy > 0, "Nuclear repulsion should be positive"
    
    # Check coordinates
    coords = mol.coordinates_angstrom
    assert coords.shape == (3, 3), f"Expected (3,3) coordinates, got {coords.shape}"
    
    print(f"✓ Water molecule: {mol.n_atoms} atoms, {mol.n_electrons} electrons")
    print(f"✓ Nuclear repulsion: {mol.nuclear_repulsion_energy:.6f} Hartree")
    return True


def main():
    """Run all tests."""
    print("Running ChemWise Tests")
    print("=" * 50)
    
    tests = [
        test_molecule_properties,
        test_h2_hf,
        test_h2_dft,
        test_water_hf,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
            print()
        except Exception as e:
            print(f"✗ Test failed: {test.__name__}")
            print(f"  Error: {e}")
            failed += 1
            print()
    
    print("=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("All tests passed! ✓")
        return 0
    else:
        print("Some tests failed! ✗")
        return 1


if __name__ == "__main__":
    sys.exit(main())
