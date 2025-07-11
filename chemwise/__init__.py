"""
ChemWise - A Comprehensive Quantum Chemistry Software Suite

A Python-based quantum chemistry package supporting Hartree-Fock and DFT calculations.
"""

__version__ = "0.1.0"
__author__ = "Justin Kirkland"
__email__ = "your.email@example.com"

# Core imports
try:
    from .core.molecule import Molecule
    from .core.basis import BasisSet
    from .core.calculator import Calculator
    
    # Method imports
    from .methods.hf import HartreeFock
    from .methods.dft import DFT
    
    # Utility imports
    from .utils.constants import ATOMIC_NUMBERS, ATOMIC_SYMBOLS
    from .utils.io import load_molecule, save_results
    
    # Make key classes available at package level
    __all__ = [
        "Molecule",
        "BasisSet", 
        "Calculator",
        "HartreeFock",
        "DFT",
        "load_molecule",
        "save_results",
    ]
except ImportError as e:
    # Handle import errors gracefully during development
    print(f"Warning: Some imports failed during package initialization: {e}")
    __all__ = []
