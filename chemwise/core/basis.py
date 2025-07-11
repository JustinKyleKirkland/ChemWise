"""
Basis set handling for quantum chemistry calculations.
"""

import numpy as np
from typing import List, Dict, Tuple, Union
from pathlib import Path
import json

from ..utils.constants import ATOMIC_NUMBERS, DEFAULT_BASIS_SETS


class PrimitiveGaussian:
    """A primitive Gaussian function for basis sets."""
    
    def __init__(self, exponent: float, coefficient: float, 
                 angular_momentum: Tuple[int, int, int] = (0, 0, 0)):
        """
        Initialize a primitive Gaussian.
        
        Args:
            exponent: Gaussian exponent (alpha)
            coefficient: Contraction coefficient
            angular_momentum: (l, m, n) angular momentum quantum numbers
        """
        self.exponent = float(exponent)
        self.coefficient = float(coefficient)
        self.angular_momentum = angular_momentum
        self.l, self.m, self.n = angular_momentum
        self.total_angular_momentum = sum(angular_momentum)
    
    def __str__(self) -> str:
        return f"PrimitiveGaussian(α={self.exponent:.6f}, c={self.coefficient:.6f}, l={self.angular_momentum})"


class ContractedGaussian:
    """A contracted Gaussian function (linear combination of primitives)."""
    
    def __init__(self, primitives: List[PrimitiveGaussian], 
                 center: np.ndarray, shell_type: str = 's'):
        """
        Initialize a contracted Gaussian.
        
        Args:
            primitives: List of primitive Gaussians
            center: Nuclear center coordinates
            shell_type: Shell type ('s', 'p', 'd', 'f', etc.)
        """
        self.primitives = primitives
        self.center = np.array(center, dtype=float)
        self.shell_type = shell_type.lower()
        
        # All primitives should have same angular momentum
        if primitives:
            self.angular_momentum = primitives[0].angular_momentum
            self.total_angular_momentum = primitives[0].total_angular_momentum
        else:
            self.angular_momentum = (0, 0, 0)
            self.total_angular_momentum = 0
        
        # Normalize the contracted function
        self._normalize()
    
    def _normalize(self) -> None:
        """Normalize the contracted Gaussian function."""
        # Calculate normalization constant
        overlap = 0.0
        for i, prim_i in enumerate(self.primitives):
            for j, prim_j in enumerate(self.primitives):
                alpha_sum = prim_i.exponent + prim_j.exponent
                overlap += (prim_i.coefficient * prim_j.coefficient * 
                           (np.pi / alpha_sum)**(3/2))
        
        norm_factor = 1.0 / np.sqrt(overlap)
        
        # Apply normalization to coefficients
        for primitive in self.primitives:
            primitive.coefficient *= norm_factor
    
    def evaluate(self, point: np.ndarray) -> float:
        """Evaluate the contracted Gaussian at a given point."""
        r = point - self.center
        r_squared = np.dot(r, r)
        
        value = 0.0
        for primitive in self.primitives:
            # Gaussian part
            gaussian = np.exp(-primitive.exponent * r_squared)
            
            # Angular part (simplified for s-type)
            if self.total_angular_momentum == 0:
                angular = 1.0
            else:
                angular = (r[0]**primitive.l * r[1]**primitive.m * r[2]**primitive.n)
            
            value += primitive.coefficient * gaussian * angular
        
        return value
    
    def __str__(self) -> str:
        return f"ContractedGaussian({self.shell_type}, center={self.center}, n_primitives={len(self.primitives)})"


class BasisFunction:
    """A single basis function with specific angular momentum components."""
    
    def __init__(self, contracted_gaussian: ContractedGaussian, 
                 angular_momentum: Tuple[int, int, int], atom_index: int):
        """
        Initialize a basis function.
        
        Args:
            contracted_gaussian: The contracted Gaussian
            angular_momentum: Specific (l, m, n) for this basis function
            atom_index: Index of the atom this function is centered on
        """
        self.contracted_gaussian = contracted_gaussian
        self.angular_momentum = angular_momentum
        self.atom_index = atom_index
        self.l, self.m, self.n = angular_momentum
        self.total_angular_momentum = sum(angular_momentum)
        self.center = contracted_gaussian.center
    
    def evaluate(self, point: np.ndarray) -> float:
        """Evaluate the basis function at a given point."""
        r = point - self.center
        r_squared = np.dot(r, r)
        
        value = 0.0
        for primitive in self.contracted_gaussian.primitives:
            # Gaussian part
            gaussian = np.exp(-primitive.exponent * r_squared)
            
            # Angular part
            angular = r[0]**self.l * r[1]**self.m * r[2]**self.n
            
            value += primitive.coefficient * gaussian * angular
        
        return value
    
    @property
    def shell_type(self) -> str:
        """Get shell type label."""
        L = self.total_angular_momentum
        if L == 0:
            return 's'
        elif L == 1:
            return 'p'
        elif L == 2:
            return 'd'
        elif L == 3:
            return 'f'
        else:
            return f'L{L}'
    
    def __str__(self) -> str:
        return f"BasisFunction({self.shell_type}, l={self.angular_momentum}, atom={self.atom_index})"


class BasisSet:
    """Collection of basis functions for a molecular system."""
    
    def __init__(self, basis_functions: List[BasisFunction] = None):
        """
        Initialize a basis set.
        
        Args:
            basis_functions: List of basis functions
        """
        self.basis_functions = basis_functions if basis_functions is not None else []
        self.n_basis = len(self.basis_functions)
    
    @classmethod
    def from_molecule(cls, molecule, basis_name: str = 'sto-3g') -> 'BasisSet':
        """
        Create basis set for a molecule.
        
        Args:
            molecule: Molecule object
            basis_name: Name of basis set
            
        Returns:
            BasisSet object
        """
        from .molecule import Molecule
        
        basis_functions = []
        
        for atom_idx, atom in enumerate(molecule.atoms):
            # Load basis set data for this atom
            atom_basis = load_atomic_basis(atom.atomic_number, basis_name)
            
            for shell_data in atom_basis:
                shell_type = shell_data['shell_type']
                exponents = shell_data['exponents']
                coefficients = shell_data['coefficients']
                
                # Create primitives
                primitives = []
                for exp, coeff in zip(exponents, coefficients):
                    primitives.append(PrimitiveGaussian(exp, coeff))
                
                # Create contracted Gaussian
                contracted = ContractedGaussian(primitives, atom.coordinates, shell_type)
                
                # Generate basis functions for this shell
                shell_functions = generate_shell_functions(contracted, shell_type, atom_idx)
                basis_functions.extend(shell_functions)
        
        basis_set = cls(basis_functions)
        basis_set.basis_name = basis_name
        return basis_set
    
    def add_function(self, basis_function: BasisFunction) -> None:
        """Add a basis function to the set."""
        self.basis_functions.append(basis_function)
        self.n_basis = len(self.basis_functions)
    
    def get_functions_on_atom(self, atom_index: int) -> List[BasisFunction]:
        """Get all basis functions centered on a specific atom."""
        return [bf for bf in self.basis_functions if bf.atom_index == atom_index]
    
    def get_functions_by_shell(self, shell_type: str) -> List[BasisFunction]:
        """Get all basis functions of a specific shell type."""
        return [bf for bf in self.basis_functions if bf.shell_type == shell_type]
    
    @property
    def atom_indices(self) -> np.ndarray:
        """Array mapping each basis function to its atom index."""
        return np.array([bf.atom_index for bf in self.basis_functions])
    
    @property
    def shell_types(self) -> List[str]:
        """List of shell types for each basis function."""
        return [bf.shell_type for bf in self.basis_functions]
    
    @property
    def angular_momenta(self) -> List[Tuple[int, int, int]]:
        """List of angular momentum quantum numbers for each basis function."""
        return [bf.angular_momentum for bf in self.basis_functions]
    
    @property
    def centers(self) -> np.ndarray:
        """Array of centers for each basis function."""
        return np.array([bf.center for bf in self.basis_functions])
    
    def __str__(self) -> str:
        shell_counts = {}
        for bf in self.basis_functions:
            shell_type = bf.shell_type
            shell_counts[shell_type] = shell_counts.get(shell_type, 0) + 1
        
        shell_summary = ", ".join([f"{count}{shell}" for shell, count in sorted(shell_counts.items())])
        return f"BasisSet(n_basis={self.n_basis}, shells=[{shell_summary}])"
    
    def __len__(self) -> int:
        return self.n_basis
    
    def __getitem__(self, index: int) -> BasisFunction:
        return self.basis_functions[index]


def generate_shell_functions(contracted: ContractedGaussian, shell_type: str, 
                           atom_index: int) -> List[BasisFunction]:
    """
    Generate basis functions for a shell.
    
    Args:
        contracted: Contracted Gaussian
        shell_type: Shell type ('s', 'p', 'd', etc.)
        atom_index: Index of the atom
        
    Returns:
        List of basis functions
    """
    functions = []
    
    if shell_type == 's':
        # s orbital: (0,0,0)
        functions.append(BasisFunction(contracted, (0, 0, 0), atom_index))
    
    elif shell_type == 'p':
        # p orbitals: px, py, pz
        functions.append(BasisFunction(contracted, (1, 0, 0), atom_index))  # px
        functions.append(BasisFunction(contracted, (0, 1, 0), atom_index))  # py
        functions.append(BasisFunction(contracted, (0, 0, 1), atom_index))  # pz
    
    elif shell_type == 'd':
        # d orbitals: dxx, dyy, dzz, dxy, dxz, dyz
        functions.append(BasisFunction(contracted, (2, 0, 0), atom_index))  # dxx
        functions.append(BasisFunction(contracted, (0, 2, 0), atom_index))  # dyy
        functions.append(BasisFunction(contracted, (0, 0, 2), atom_index))  # dzz
        functions.append(BasisFunction(contracted, (1, 1, 0), atom_index))  # dxy
        functions.append(BasisFunction(contracted, (1, 0, 1), atom_index))  # dxz
        functions.append(BasisFunction(contracted, (0, 1, 1), atom_index))  # dyz
    
    elif shell_type == 'f':
        # f orbitals: (3,0,0), (0,3,0), (0,0,3), (2,1,0), (2,0,1), (1,2,0), (0,2,1), (1,0,2), (0,1,2), (1,1,1)
        for l in range(4):
            for m in range(4-l):
                n = 3 - l - m
                if n >= 0:
                    functions.append(BasisFunction(contracted, (l, m, n), atom_index))
    
    return functions


def load_atomic_basis(atomic_number: int, basis_name: str) -> List[Dict]:
    """
    Load basis set data for a specific atom.
    
    Args:
        atomic_number: Atomic number
        basis_name: Name of basis set
        
    Returns:
        List of shell data dictionaries
    """
    basis_name = basis_name.lower()
    
    # First try to load from JSON files in basis/ directory
    try:
        return load_basis_from_file(atomic_number, basis_name)
    except (FileNotFoundError, KeyError):
        pass
    
    # Fall back to hardcoded basis sets
    if basis_name == 'sto-3g':
        return get_sto_3g_basis(atomic_number)
    elif basis_name == 'sto-6g':
        return get_sto_6g_basis(atomic_number)
    else:
        raise ValueError(f"Basis set '{basis_name}' not implemented")


def load_basis_from_file(atomic_number: int, basis_name: str) -> List[Dict]:
    """
    Load basis set data from JSON files in the basis/ directory.
    
    Args:
        atomic_number: Atomic number
        basis_name: Name of basis set
        
    Returns:
        List of shell data dictionaries
    """
    from ..utils.constants import ATOMIC_SYMBOLS
    
    # Get the path to basis directory
    basis_dir = Path(__file__).parent.parent.parent / "basis"
    basis_file = basis_dir / f"{basis_name}.json"
    
    if not basis_file.exists():
        raise FileNotFoundError(f"Basis set file {basis_file} not found")
    
    # Load basis set data
    with open(basis_file, 'r') as f:
        basis_data = json.load(f)
    
    # Get element symbol
    element_symbol = ATOMIC_SYMBOLS.get(atomic_number)
    if element_symbol is None:
        raise ValueError(f"Unknown atomic number: {atomic_number}")
    
    # Get basis set data for this element
    if str(atomic_number) not in basis_data["elements"]:
        raise KeyError(f"Element {element_symbol} (Z={atomic_number}) not found in basis set {basis_name}")
    
    element_data = basis_data["elements"][str(atomic_number)]
    
    # Convert to the expected format
    shells = []
    for shell in element_data["shells"]:
        shells.append({
            'shell_type': shell['shell_type'],
            'exponents': shell['exponents'],
            'coefficients': shell['coefficients']
        })
    
    return shells


def list_available_basis_sets() -> List[str]:
    """
    List all available basis sets (both hardcoded and from files).
    
    Returns:
        List of available basis set names
    """
    available = ['sto-3g', 'sto-6g']  # Hardcoded basis sets
    
    # Add basis sets from files
    basis_dir = Path(__file__).parent.parent.parent / "basis"
    if basis_dir.exists():
        for basis_file in basis_dir.glob("*.json"):
            basis_name = basis_file.stem
            if basis_name not in available:
                available.append(basis_name)
    
    return sorted(available)


def get_sto_3g_basis(atomic_number: int) -> List[Dict]:
    """Get STO-3G basis set data for an atom."""
    
    # Simplified STO-3G basis sets for first few elements
    sto_3g_data = {
        1: [  # Hydrogen
            {
                'shell_type': 's',
                'exponents': [3.42525091, 0.62391373, 0.16885540],
                'coefficients': [0.15432897, 0.53532814, 0.44463454]
            }
        ],
        2: [  # Helium
            {
                'shell_type': 's',
                'exponents': [6.36242139, 1.15892300, 0.31364979],
                'coefficients': [0.15432897, 0.53532814, 0.44463454]
            }
        ],
        3: [  # Lithium
            {
                'shell_type': 's',
                'exponents': [16.11957475, 2.93620070, 0.79465050],
                'coefficients': [0.15432897, 0.53532814, 0.44463454]
            },
            {
                'shell_type': 's',
                'exponents': [0.63628970, 0.14786010, 0.04808870],
                'coefficients': [-0.09996723, 0.39951283, 0.70011547]
            }
        ],
        6: [  # Carbon
            {
                'shell_type': 's',
                'exponents': [71.6168370, 13.0450960, 3.5305122],
                'coefficients': [0.15432897, 0.53532814, 0.44463454]
            },
            {
                'shell_type': 's',
                'exponents': [2.9412494, 0.6834831, 0.2222899],
                'coefficients': [-0.09996723, 0.39951283, 0.70011547]
            },
            {
                'shell_type': 'p',
                'exponents': [2.9412494, 0.6834831, 0.2222899],
                'coefficients': [0.15591627, 0.60768372, 0.39195739]
            }
        ],
        8: [  # Oxygen
            {
                'shell_type': 's',
                'exponents': [130.7093200, 23.8088610, 6.4436083],
                'coefficients': [0.15432897, 0.53532814, 0.44463454]
            },
            {
                'shell_type': 's',
                'exponents': [5.0331513, 1.1695961, 0.3803890],
                'coefficients': [-0.09996723, 0.39951283, 0.70011547]
            },
            {
                'shell_type': 'p',
                'exponents': [5.0331513, 1.1695961, 0.3803890],
                'coefficients': [0.15591627, 0.60768372, 0.39195739]
            }
        ]
    }
    
    if atomic_number in sto_3g_data:
        return sto_3g_data[atomic_number]
    else:
        raise ValueError(f"STO-3G basis not available for atomic number {atomic_number}")


def get_sto_6g_basis(atomic_number: int) -> List[Dict]:
    """Get STO-6G basis set data for an atom."""
    # Placeholder - would contain actual STO-6G data
    # For now, just use STO-3G
    return get_sto_3g_basis(atomic_number)
