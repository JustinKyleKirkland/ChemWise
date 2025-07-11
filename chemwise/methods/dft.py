"""
Density Functional Theory (DFT) implementation.
"""

import numpy as np
from typing import Dict, Tuple, Optional
import time

from ..core.basis import BasisSet
from ..core.integrals import IntegralEngine
from ..utils.math_utils import symmetric_orthogonalization, diagonalize_symmetric, density_matrix as build_density_matrix
from ..utils.constants import DEFAULT_SCF_CONVERGENCE, DEFAULT_MAX_SCF_ITERATIONS, DFT_GRID_SIZES
from .functionals import get_functional, XCFunctional
from .hf import SCFResult


class DFTGrid:
    """Numerical integration grid for DFT calculations."""
    
    def __init__(self, molecule, grid_level: str = 'medium'):
        """
        Initialize DFT grid.
        
        Args:
            molecule: Molecule object
            grid_level: Grid quality ('coarse', 'medium', 'fine', 'ultrafine')
        """
        self.molecule = molecule
        self.grid_level = grid_level
        
        if grid_level in DFT_GRID_SIZES:
            self.n_radial, self.n_angular = DFT_GRID_SIZES[grid_level]
        else:
            raise ValueError(f"Unknown grid level: {grid_level}")
        
        self.points = None
        self.weights = None
        self.n_points = 0
        
        self._generate_grid()
    
    def _generate_grid(self) -> None:
        """Generate numerical integration grid points and weights."""
        all_points = []
        all_weights = []
        
        for atom in self.molecule.atoms:
            # Generate atomic grid
            atomic_points, atomic_weights = self._generate_atomic_grid(atom.coordinates)
            all_points.extend(atomic_points)
            all_weights.extend(atomic_weights)
        
        self.points = np.array(all_points)
        self.weights = np.array(all_weights)
        self.n_points = len(self.points)
    
    def _generate_atomic_grid(self, center: np.ndarray) -> Tuple[list, list]:
        """Generate grid points around an atomic center."""
        points = []
        weights = []
        
        # Simplified grid generation - in practice would use more sophisticated methods
        # like Becke partitioning and Lebedev angular grids
        
        # Radial points (simplified Gauss-Chebyshev)
        for i in range(self.n_radial):
            # Simple radial distribution
            r = 0.1 + i * 5.0 / self.n_radial
            r_weight = r**2
            
            # Angular points (simplified spherical grid)
            for j in range(self.n_angular):
                theta = np.pi * j / self.n_angular
                phi = 2 * np.pi * j / self.n_angular
                
                # Convert to Cartesian
                x = center[0] + r * np.sin(theta) * np.cos(phi)
                y = center[1] + r * np.sin(theta) * np.sin(phi)
                z = center[2] + r * np.cos(theta)
                
                points.append([x, y, z])
                weights.append(r_weight * np.sin(theta) * np.pi / self.n_angular)
        
        return points, weights


class DFT:
    """Density Functional Theory implementation."""
    
    def __init__(self, molecule, basis_set: BasisSet, functional: str = 'b3lyp',
                 convergence: float = DEFAULT_SCF_CONVERGENCE,
                 max_iterations: int = DEFAULT_MAX_SCF_ITERATIONS,
                 grid_level: str = 'medium'):
        """
        Initialize DFT calculation.
        
        Args:
            molecule: Molecule object
            basis_set: Basis set for the calculation
            functional: XC functional name
            convergence: SCF convergence threshold
            max_iterations: Maximum number of SCF iterations
            grid_level: Grid quality for numerical integration
        """
        self.molecule = molecule
        self.basis_set = basis_set
        self.functional_name = functional
        self.convergence = convergence
        self.max_iterations = max_iterations
        
        # Initialize functional
        self.functional = get_functional(functional)
        
        # Initialize grid
        self.grid = DFTGrid(molecule, grid_level)
        
        # Initialize integral engine
        self.integral_engine = IntegralEngine(basis_set)
        
        # Calculation properties
        self.n_basis = basis_set.n_basis
        self.n_electrons = molecule.n_electrons
        self.n_occupied = molecule.n_electron_pairs
        
        # Matrices
        self.S = None  # Overlap matrix
        self.T = None  # Kinetic energy matrix
        self.V = None  # Nuclear attraction matrix
        self.H_core = None  # Core Hamiltonian
        self.ERI = None  # Electron repulsion integrals (for hybrids)
        self.X = None  # Orthogonalization matrix
        
        # Check if hybrid functional
        self.is_hybrid = hasattr(self.functional, 'exact_exchange_fraction')
        if self.is_hybrid:
            self.exact_exchange_fraction = self.functional.exact_exchange_fraction
        else:
            self.exact_exchange_fraction = 0.0
    
    def setup_integrals(self) -> None:
        """Calculate all required integrals."""
        print("Calculating integrals...")
        start_time = time.time()
        
        # One-electron integrals
        self.S = self.integral_engine.overlap_matrix()
        self.T = self.integral_engine.kinetic_matrix()
        self.V = self.integral_engine.nuclear_attraction_matrix(
            self.molecule.nuclear_charges,
            self.molecule.coordinates_bohr
        )
        self.H_core = self.T + self.V
        
        # Two-electron integrals (only for hybrid functionals)
        if self.is_hybrid:
            print("Calculating electron repulsion integrals for hybrid functional...")
            self.ERI = self.integral_engine.electron_repulsion_tensor()
        
        # Orthogonalization matrix
        self.X = symmetric_orthogonalization(self.S)
        
        integral_time = time.time() - start_time
        print(f"Integral calculation completed in {integral_time:.2f} seconds")
    
    def evaluate_basis_on_grid(self) -> np.ndarray:
        """Evaluate all basis functions on grid points."""
        basis_values = np.zeros((self.grid.n_points, self.n_basis))
        
        for i, bf in enumerate(self.basis_set.basis_functions):
            for j, point in enumerate(self.grid.points):
                basis_values[j, i] = bf.evaluate(point)
        
        return basis_values
    
    def calculate_density_on_grid(self, density_matrix: np.ndarray,
                                 basis_values: np.ndarray) -> np.ndarray:
        """Calculate electron density on grid points."""
        # ρ(r) = ∑_μν P_μν φ_μ(r) φ_ν(r)
        density = np.zeros(self.grid.n_points)
        
        for i in range(self.grid.n_points):
            phi = basis_values[i, :]
            density[i] = np.sum(density_matrix * np.outer(phi, phi))
        
        return density
    
    def calculate_xc_matrix(self, density_matrix: np.ndarray,
                           basis_values: np.ndarray) -> np.ndarray:
        """Calculate exchange-correlation matrix."""
        # Calculate density on grid
        rho = self.calculate_density_on_grid(density_matrix, basis_values)
        
        # Calculate XC potential on grid
        # For now, assume LDA (no gradient dependence)
        V_xc_grid = self.functional.xc_potential(rho)
        
        # Integrate to get XC matrix elements
        V_xc = np.zeros((self.n_basis, self.n_basis))
        
        for mu in range(self.n_basis):
            for nu in range(self.n_basis):
                integrand = basis_values[:, mu] * V_xc_grid * basis_values[:, nu]
                V_xc[mu, nu] = np.sum(integrand * self.grid.weights)
        
        return V_xc
    
    def calculate_xc_energy(self, density_matrix: np.ndarray,
                           basis_values: np.ndarray) -> float:
        """Calculate exchange-correlation energy."""
        # Calculate density on grid
        rho = self.calculate_density_on_grid(density_matrix, basis_values)
        
        # Calculate XC energy density
        eps_xc = self.functional.xc_energy_density(rho)
        
        # Integrate
        E_xc = np.sum(rho * eps_xc * self.grid.weights)
        
        return E_xc
    
    def build_fock_matrix(self, density_matrix: np.ndarray,
                         basis_values: np.ndarray) -> np.ndarray:
        """Build Kohn-Sham Fock matrix."""
        F = self.H_core.copy()
        
        # Add Coulomb interaction (classical electron repulsion)
        for i in range(self.n_basis):
            for j in range(self.n_basis):
                for k in range(self.n_basis):
                    for ll in range(self.n_basis):
                        # Coulomb term
                        F[i, j] += density_matrix[k, ll] * self.ERI[i, j, k, ll] if self.ERI is not None else 0.0
        
        # Add exchange-correlation potential
        V_xc = self.calculate_xc_matrix(density_matrix, basis_values)
        F += V_xc
        
        # For hybrid functionals, add exact exchange
        if self.is_hybrid and self.ERI is not None:
            for i in range(self.n_basis):
                for j in range(self.n_basis):
                    for k in range(self.n_basis):
                        for ll in range(self.n_basis):
                            # Exact exchange term
                            F[i, j] -= (self.exact_exchange_fraction * 
                                      density_matrix[k, ll] * self.ERI[i, k, j, ll])
        
        return F
    
    def calculate_energy(self, density_matrix: np.ndarray, fock_matrix: np.ndarray,
                        basis_values: np.ndarray) -> float:
        """Calculate total DFT energy."""
        # Electronic energy from core and Coulomb terms
        E_electronic = 0.5 * np.trace(density_matrix @ (self.H_core + fock_matrix))
        
        # Add exchange-correlation energy
        E_xc = self.calculate_xc_energy(density_matrix, basis_values)
        
        # For hybrid functionals, we need to correct for double-counting
        if self.is_hybrid:
            # Remove overcounted XC energy and add correct amount
            E_electronic = E_electronic - 0.5 * np.trace(density_matrix @ self.calculate_xc_matrix(density_matrix, basis_values))
            E_electronic += E_xc
        else:
            E_electronic = E_electronic - 0.5 * np.trace(density_matrix @ self.calculate_xc_matrix(density_matrix, basis_values)) + E_xc
        
        return E_electronic
    
    def guess_density(self, method: str = 'core') -> np.ndarray:
        """Generate initial density matrix guess."""
        if method == 'core':
            # Core Hamiltonian guess
            H_core_ortho = self.X.T @ self.H_core @ self.X
            eigenvals, eigenvecs = diagonalize_symmetric(H_core_ortho)
            C_ortho = eigenvecs[:, :self.n_occupied]
            C = self.X @ C_ortho
            P = build_density_matrix(C, self.n_occupied)
        else:
            # Random density matrix
            P = np.random.random((self.n_basis, self.n_basis))
            P = (P + P.T) / 2
        
        return P
    
    def scf_iteration(self, density_matrix: np.ndarray,
                     basis_values: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Perform one SCF iteration."""
        # Build Fock matrix
        F = self.build_fock_matrix(density_matrix, basis_values)
        
        # Transform to orthogonal basis
        F_ortho = self.X.T @ F @ self.X
        
        # Diagonalize
        orbital_energies, C_ortho = diagonalize_symmetric(F_ortho)
        
        # Transform back to original basis
        C = self.X @ C_ortho
        
        # Build new density matrix
        P_new = build_density_matrix(C, self.n_occupied)
        
        return P_new, F, C, orbital_energies
    
    def run_scf(self, initial_guess: str = 'core') -> SCFResult:
        """Run the SCF procedure for DFT."""
        print(f"Starting DFT calculation with {self.functional_name} functional")
        print(f"Grid level: {self.grid.grid_level} ({self.grid.n_points} points)")
        print(f"Convergence threshold: {self.convergence}")
        print(f"Maximum iterations: {self.max_iterations}")
        
        if self.is_hybrid:
            print(f"Hybrid functional with {self.exact_exchange_fraction:.2%} exact exchange")
        
        print()
        
        # Setup integrals
        self.setup_integrals()
        
        # Evaluate basis functions on grid
        print("Evaluating basis functions on grid...")
        basis_values = self.evaluate_basis_on_grid()
        
        # Initial guess
        P = self.guess_density(initial_guess)
        
        # Initialize result object
        result = SCFResult()
        result.nuclear_repulsion = self.molecule.nuclear_repulsion_energy
        
        print(f"{'Iter':<6} {'E(elec)':<15} {'E(total)':<15} {'ΔE':<12} {'RMS(P)':<12}")
        print("-" * 70)
        
        # SCF iterations
        previous_energy = 0.0
        previous_density = np.zeros_like(P)
        
        for iteration in range(self.max_iterations):
            # Perform SCF iteration
            P_new, F, C, orbital_energies = self.scf_iteration(P, basis_values)
            
            # Calculate energy
            electronic_energy = self.calculate_energy(P_new, F, basis_values)
            total_energy = electronic_energy + result.nuclear_repulsion
            
            # Check convergence
            energy_change = abs(total_energy - previous_energy)
            density_change = np.sqrt(np.mean((P_new - previous_density)**2))
            
            print(f"{iteration+1:<6} {electronic_energy:<15.8f} {total_energy:<15.8f} "
                  f"{energy_change:<12.2e} {density_change:<12.2e}")
            
            # Check for convergence
            if energy_change < self.convergence and density_change < self.convergence:
                result.converged = True
                result.n_iterations = iteration + 1
                result.final_error = max(energy_change, density_change)
                break
            
            # Update for next iteration
            P = P_new
            previous_energy = total_energy
            previous_density = P_new.copy()
        
        if not result.converged:
            print(f"\nSCF did not converge in {self.max_iterations} iterations")
        else:
            print(f"\nSCF converged in {result.n_iterations} iterations")
        
        # Store final results
        result.energy = total_energy
        result.electronic_energy = electronic_energy
        result.orbital_energies = orbital_energies
        result.orbital_coefficients = C
        result.density_matrix = P
        result.fock_matrix = F
        result.overlap_matrix = self.S
        
        return result
    
    def run(self, initial_guess: str = 'core') -> SCFResult:
        """
        Run DFT calculation.
        
        Args:
            initial_guess: Method for initial density guess
            
        Returns:
            SCF results
        """
        return self.run_scf(initial_guess)
