"""
Main calculator interface for ChemWise quantum chemistry calculations.
"""

import numpy as np
from typing import Union, Optional
from pathlib import Path
import time

from .molecule import Molecule
from .basis import BasisSet
from ..methods.hf import HartreeFock
from ..methods.dft import DFT
from ..utils.io import load_molecule, save_results, print_geometry, print_orbital_energies
from ..utils.math_utils import mulliken_charges


class CalculationResult:
    """Container for quantum chemistry calculation results."""
    
    def __init__(self):
        self.method: str = ""
        self.basis_set_name: str = ""
        self.energy: float = 0.0
        self.nuclear_repulsion: float = 0.0
        self.electronic_energy: float = 0.0
        self.orbital_energies: np.ndarray = np.array([])
        self.orbital_coefficients: np.ndarray = np.array([])
        self.density_matrix: np.ndarray = np.array([])
        self.mulliken_charges: np.ndarray = np.array([])
        self.dipole_moment: np.ndarray = np.array([])
        self.converged: bool = False
        self.n_iterations: int = 0
        self.calculation_time: float = 0.0
        self.molecule: Optional[Molecule] = None
        self.basis_set: Optional[BasisSet] = None
    
    def summary(self) -> str:
        """Generate a summary of the calculation results."""
        lines = []
        lines.append("=" * 60)
        lines.append("ChemWise Quantum Chemistry Calculation Results")
        lines.append("=" * 60)
        lines.append(f"Method: {self.method}")
        lines.append(f"Basis Set: {self.basis_set_name}")
        lines.append(f"Convergence: {'Yes' if self.converged else 'No'}")
        if self.converged:
            lines.append(f"Iterations: {self.n_iterations}")
        lines.append(f"Calculation Time: {self.calculation_time:.2f} seconds")
        lines.append("")
        
        # Energy results
        lines.append("Energy Results:")
        lines.append("-" * 30)
        lines.append(f"Nuclear Repulsion Energy: {self.nuclear_repulsion:15.8f} Hartree")
        lines.append(f"Electronic Energy:        {self.electronic_energy:15.8f} Hartree")
        lines.append(f"Total Energy:             {self.energy:15.8f} Hartree")
        lines.append("")
        
        # Molecular properties
        if len(self.mulliken_charges) > 0:
            lines.append("Mulliken Charges:")
            lines.append("-" * 20)
            if self.molecule:
                for i, (symbol, charge) in enumerate(zip(self.molecule.symbols, self.mulliken_charges)):
                    lines.append(f"  {symbol}{i+1:2d}: {charge:8.4f}")
            lines.append("")
        
        if len(self.dipole_moment) > 0:
            dipole_magnitude = np.linalg.norm(self.dipole_moment)
            lines.append(f"Dipole Moment: {dipole_magnitude:.4f} Debye")
            lines.append(f"  X: {self.dipole_moment[0]:8.4f}")
            lines.append(f"  Y: {self.dipole_moment[1]:8.4f}")
            lines.append(f"  Z: {self.dipole_moment[2]:8.4f}")
            lines.append("")
        
        return "\n".join(lines)


class Calculator:
    """Main calculator interface for quantum chemistry calculations."""
    
    def __init__(self, molecule: Union[Molecule, str, Path], 
                 method: str = 'hf', basis: str = 'sto-3g',
                 **kwargs):
        """
        Initialize calculator.
        
        Args:
            molecule: Molecule object or path to molecule file
            method: Calculation method ('hf', 'dft', 'b3lyp', etc.)
            basis: Basis set name
            **kwargs: Additional calculation parameters
        """
        # Load molecule if needed
        if isinstance(molecule, (str, Path)):
            self.molecule = Molecule.from_dict(load_molecule(molecule))
        elif isinstance(molecule, dict):
            self.molecule = Molecule.from_dict(molecule)
        else:
            self.molecule = molecule
        
        self.method = method.lower()
        self.basis_name = basis.lower()
        self.kwargs = kwargs
        
        # Create basis set
        self.basis_set = BasisSet.from_molecule(self.molecule, self.basis_name)
        
        # Initialize calculation engine
        self.engine = None
        self._setup_engine()
    
    def _setup_engine(self) -> None:
        """Setup the appropriate calculation engine."""
        if self.method == 'hf':
            self.engine = HartreeFock(
                self.molecule, 
                self.basis_set,
                **{k: v for k, v in self.kwargs.items() if k in ['convergence', 'max_iterations', 'use_diis']}
            )
        elif self.method in ['dft', 'lda', 'b3lyp', 'b88', 'lyp']:
            functional = self.method if self.method != 'dft' else 'b3lyp'
            self.engine = DFT(
                self.molecule,
                self.basis_set,
                functional=functional,
                **{k: v for k, v in self.kwargs.items() if k in ['convergence', 'max_iterations', 'grid_level']}
            )
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    def run(self, print_results: bool = True) -> CalculationResult:
        """
        Run the quantum chemistry calculation.
        
        Args:
            print_results: Whether to print results to console
            
        Returns:
            Calculation results
        """
        if print_results:
            print(f"Starting {self.method.upper()} calculation")
            print_geometry(self.molecule.symbols, self.molecule.coordinates_angstrom)
            print(f"Charge: {self.molecule.charge}")
            print(f"Multiplicity: {self.molecule.multiplicity}")
            print(f"Number of electrons: {self.molecule.n_electrons}")
            print(f"Basis set: {self.basis_name}")
            print(f"Number of basis functions: {self.basis_set.n_basis}")
            print()
        
        start_time = time.time()
        
        # Run calculation
        scf_result = self.engine.run()
        
        calculation_time = time.time() - start_time
        
        # Create result object
        result = CalculationResult()
        result.method = self.method.upper()
        result.basis_set_name = self.basis_name
        result.energy = scf_result.energy
        result.nuclear_repulsion = scf_result.nuclear_repulsion
        result.electronic_energy = scf_result.electronic_energy
        result.orbital_energies = scf_result.orbital_energies
        result.orbital_coefficients = scf_result.orbital_coefficients
        result.density_matrix = scf_result.density_matrix
        result.converged = scf_result.converged
        result.n_iterations = scf_result.n_iterations
        result.calculation_time = calculation_time
        result.molecule = self.molecule
        result.basis_set = self.basis_set
        
        # Calculate molecular properties
        if result.converged:
            self._calculate_properties(result)
        
        if print_results:
            print(result.summary())
            
            if result.converged:
                # Print orbital energies
                occupations = np.zeros_like(result.orbital_energies)
                occupations[:self.molecule.n_electron_pairs] = 2.0
                print_orbital_energies(result.orbital_energies, occupations)
        
        return result
    
    def _calculate_properties(self, result: CalculationResult) -> None:
        """Calculate molecular properties from SCF results."""
        # Mulliken charges
        result.mulliken_charges = mulliken_charges(
            result.density_matrix,
            self.engine.S,
            self.molecule.nuclear_charges,
            self.basis_set.atom_indices
        )
        
        # Dipole moment (simplified calculation)
        result.dipole_moment = self._calculate_dipole_moment(result)
    
    def _calculate_dipole_moment(self, result: CalculationResult) -> np.ndarray:
        """Calculate electric dipole moment."""
        # Simplified dipole calculation
        # In practice, this would require dipole integral evaluation
        
        # Nuclear contribution
        nuclear_dipole = np.zeros(3)
        for i, atom in enumerate(self.molecule.atoms):
            nuclear_dipole += atom.charge * atom.coordinates
        
        # Electronic contribution (simplified)
        # This is a placeholder - proper implementation requires dipole integrals
        electronic_dipole = np.zeros(3)
        
        # Total dipole (in atomic units)
        total_dipole = nuclear_dipole - electronic_dipole
        
        # Convert to Debye
        from ..utils.constants import DEBYE_TO_AU
        return total_dipole / DEBYE_TO_AU
    
    def save_results(self, filename: Union[str, Path], result: CalculationResult) -> None:
        """
        Save calculation results to file.
        
        Args:
            filename: Output file path
            result: Calculation results to save
        """
        data = {
            'calculation': {
                'method': result.method,
                'basis_set': result.basis_set_name,
                'converged': result.converged,
                'n_iterations': result.n_iterations,
                'calculation_time': result.calculation_time
            },
            'molecule': self.molecule.to_dict(),
            'energies': {
                'total': result.energy,
                'electronic': result.electronic_energy,
                'nuclear_repulsion': result.nuclear_repulsion
            },
            'properties': {
                'mulliken_charges': result.mulliken_charges.tolist() if len(result.mulliken_charges) > 0 else [],
                'dipole_moment': result.dipole_moment.tolist() if len(result.dipole_moment) > 0 else []
            },
            'orbitals': {
                'energies': result.orbital_energies.tolist(),
                'coefficients': result.orbital_coefficients.tolist()
            },
            'matrices': {
                'density': result.density_matrix.tolist(),
                'overlap': self.engine.S.tolist() if self.engine.S is not None else []
            }
        }
        
        save_results(filename, data)


# Convenience functions
def calculate(molecule: Union[Molecule, str, Path], method: str = 'hf', 
              basis: str = 'sto-3g', **kwargs) -> CalculationResult:
    """
    Convenience function to run a quantum chemistry calculation.
    
    Args:
        molecule: Molecule object or path to molecule file
        method: Calculation method
        basis: Basis set name
        **kwargs: Additional calculation parameters
        
    Returns:
        Calculation results
    """
    calc = Calculator(molecule, method, basis, **kwargs)
    return calc.run()


def optimize_geometry(molecule: Union[Molecule, str, Path], method: str = 'hf',
                     basis: str = 'sto-3g', **kwargs) -> CalculationResult:
    """
    Placeholder for geometry optimization.
    
    Args:
        molecule: Molecule object or path to molecule file
        method: Calculation method
        basis: Basis set name
        **kwargs: Additional parameters
        
    Returns:
        Optimized geometry and final energy
    """
    # Placeholder - would implement geometry optimization
    print("Geometry optimization not yet implemented")
    return calculate(molecule, method, basis, **kwargs)
