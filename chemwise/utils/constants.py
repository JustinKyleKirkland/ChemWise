"""
Physical and mathematical constants used in quantum chemistry calculations.
"""

import math
import numpy as np

# Physical constants (in atomic units unless specified)
BOHR_TO_ANGSTROM = 0.52917721067
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_EV = 27.211386245988
EV_TO_HARTREE = 1.0 / HARTREE_TO_EV
HARTREE_TO_KCAL_MOL = 627.5094740631
DEBYE_TO_AU = 0.393430307

# Mathematical constants
PI = math.pi
SQRT_PI = math.sqrt(math.pi)
TWO_PI = 2.0 * math.pi
FOUR_PI = 4.0 * math.pi

# Speed of light (atomic units)
C_LIGHT = 137.035999084

# Fine structure constant
ALPHA = 1.0 / C_LIGHT

# Electron charge (atomic units = 1)
ELECTRON_CHARGE = -1.0

# Electron mass (atomic units = 1)
ELECTRON_MASS = 1.0

# Convergence thresholds
DEFAULT_SCF_CONVERGENCE = 1e-8
DEFAULT_GEOMETRY_CONVERGENCE = 1e-6
DEFAULT_MAX_SCF_ITERATIONS = 100

# Integration grids
DFT_GRID_SIZES = {
    'coarse': (50, 194),    # (radial, angular)
    'medium': (75, 302),
    'fine': (99, 590),
    'ultrafine': (99, 974)
}

# Atomic data
ATOMIC_NUMBERS = {
    'H': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'NE': 10, 'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15,
    'S': 16, 'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'SC': 21, 'TI': 22,
    'V': 23, 'CR': 24, 'MN': 25, 'FE': 26, 'CO': 27, 'NI': 28, 'CU': 29,
    'ZN': 30, 'GA': 31, 'GE': 32, 'AS': 33, 'SE': 34, 'BR': 35, 'KR': 36
}

ATOMIC_SYMBOLS = {v: k for k, v in ATOMIC_NUMBERS.items()}

# Atomic masses (in atomic mass units)
ATOMIC_MASSES = {
    1: 1.00782503223, 2: 4.00260325413, 3: 7.0160034366, 4: 9.0121831,
    5: 11.0093051, 6: 12.0000000, 7: 14.0030740048, 8: 15.9949146196,
    9: 18.99840316273, 10: 19.9924401762, 11: 22.9897692820, 12: 23.985041697,
    13: 26.9815385, 14: 27.9769265350, 15: 30.973761998, 16: 31.972071174,
    17: 34.9688527, 18: 39.9623831225
}

# Covalent radii (in Angstroms)
COVALENT_RADII = {
    1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66,
    9: 0.57, 10: 0.58, 11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07,
    16: 1.05, 17: 1.02, 18: 1.06
}

# Default basis sets for each element
DEFAULT_BASIS_SETS = {
    'sto-3g': list(range(1, 37)),
    'sto-6g': list(range(1, 37)),
    '6-31g': list(range(1, 19)),
    'cc-pvdz': list(range(1, 19)),
    'cc-pvtz': list(range(1, 37)),
    'cc-pvqz': list(range(1, 19)),
    'aug-cc-pvdz': list(range(1, 37)),
    'aug-cc-pvtz': list(range(1, 37)),
    'def2-svp': list(range(1, 87)),
    'def2-tzvp': list(range(1, 87))
}
