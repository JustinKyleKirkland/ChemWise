"""
Hartree-Fock method implementation.
"""

import numpy as np
from typing import Dict, Tuple, Optional
import time

from ..core.basis import BasisSet
from ..core.integrals import IntegralEngine
from ..utils.math_utils import symmetric_orthogonalization, diagonalize_symmetric, density_matrix as build_density_matrix
from ..utils.constants import DEFAULT_SCF_CONVERGENCE, DEFAULT_MAX_SCF_ITERATIONS


class SCFResult:
    """Container for SCF calculation results."""
    
    def __init__(self):
        self.energy: float = 0.0
        self.nuclear_repulsion: float = 0.0
        self.electronic_energy: float = 0.0
        self.orbital_energies: np.ndarray = np.array([])
        self.orbital_coefficients: np.ndarray = np.array([])
        self.density_matrix: np.ndarray = np.array([])
        self.fock_matrix: np.ndarray = np.array([])
        self.overlap_matrix: np.ndarray = np.array([])
        self.converged: bool = False
        self.n_iterations: int = 0
        self.final_error: float = 0.0
    
    def __format__(self, format_spec: str) -> str:
        """Format the SCF result as the total energy."""
        return format(self.energy, format_spec)
    
    def __float__(self) -> float:
        """Convert SCF result to float (total energy)."""
        return self.energy
    
    def __str__(self) -> str:
        return f"SCFResult(energy={self.energy:.8f}, converged={self.converged}, iterations={self.n_iterations})"


class HartreeFock:
    """Hartree-Fock method implementation."""
    
    def __init__(self, molecule, basis_set: BasisSet, 
                 convergence: float = DEFAULT_SCF_CONVERGENCE,
                 max_iterations: int = DEFAULT_MAX_SCF_ITERATIONS,
                 use_diis: bool = True):
        """
        Initialize Hartree-Fock calculation.
        
        Args:
            molecule: Molecule object
            basis_set: Basis set for the calculation
            convergence: SCF convergence threshold
            max_iterations: Maximum number of SCF iterations
            use_diis: Whether to use DIIS acceleration
        """
        self.molecule = molecule
        self.basis_set = basis_set
        self.convergence = convergence
        self.max_iterations = max_iterations
        self.use_diis = use_diis
        
        # Initialize integral engine
        self.integral_engine = IntegralEngine(basis_set)
        
        # Calculation properties
        self.n_basis = basis_set.n_basis
        self.n_electrons = molecule.n_electrons
        self.n_occupied = molecule.n_electron_pairs  # For RHF
        
        # Matrices
        self.S = None  # Overlap matrix
        self.T = None  # Kinetic energy matrix
        self.V = None  # Nuclear attraction matrix
        self.H_core = None  # Core Hamiltonian
        self.ERI = None  # Electron repulsion integrals
        self.X = None  # Orthogonalization matrix
        
        # DIIS data
        self.diis_errors = []
        self.diis_focks = []
        self.max_diis_vectors = 6
    
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
        
        # Two-electron integrals
        self.ERI = self.integral_engine.electron_repulsion_tensor()
        
        # Orthogonalization matrix
        self.X = symmetric_orthogonalization(self.S)
        
        integral_time = time.time() - start_time
        print(f"Integral calculation completed in {integral_time:.2f} seconds")
    
    def guess_density(self, method: str = 'core') -> np.ndarray:
        """
        Generate initial density matrix guess.
        
        Args:
            method: Method for initial guess ('core', 'sad', 'random')
            
        Returns:
            Initial density matrix
        """
        if method == 'core':
            # Core Hamiltonian guess
            H_core_ortho = self.X.T @ self.H_core @ self.X
            eigenvals, eigenvecs = diagonalize_symmetric(H_core_ortho)
            C_ortho = eigenvecs[:, :self.n_occupied]
            C = self.X @ C_ortho
            P = build_density_matrix(C, self.n_occupied)
            
        elif method == 'sad':
            # Superposition of atomic densities (placeholder)
            P = np.eye(self.n_basis) * 0.1
            
        elif method == 'random':
            # Random density matrix
            P = np.random.random((self.n_basis, self.n_basis))
            P = (P + P.T) / 2  # Make symmetric
            
        else:
            raise ValueError(f"Unknown guess method: {method}")
        
        return P
    
    def build_fock_matrix(self, density_matrix: np.ndarray) -> np.ndarray:
        """
        Build Fock matrix from density matrix.
        
        Args:
            density_matrix: Current density matrix
            
        Returns:
            Fock matrix
        """
        F = self.H_core.copy()
        
        # Add electron repulsion terms
        for i in range(self.n_basis):
            for j in range(self.n_basis):
                for k in range(self.n_basis):
                    for ll in range(self.n_basis):
                        # Coulomb term: J
                        F[i, j] += density_matrix[k, ll] * self.ERI[i, j, k, ll]
                        
                        # Exchange term: K (with factor of -0.5 for RHF)
                        F[i, j] -= 0.5 * density_matrix[k, ll] * self.ERI[i, k, j, ll]
        
        return F
    
    def calculate_energy(self, density_matrix: np.ndarray, fock_matrix: np.ndarray) -> float:
        """
        Calculate electronic energy.
        
        Args:
            density_matrix: Current density matrix
            fock_matrix: Current Fock matrix
            
        Returns:
            Electronic energy
        """
        # Electronic energy = 0.5 * Tr(P * (H_core + F))
        energy = 0.5 * np.trace(density_matrix @ (self.H_core + fock_matrix))
        return energy
    
    def diis_extrapolation(self, fock_matrix: np.ndarray, density_matrix: np.ndarray) -> np.ndarray:
        """
        Perform DIIS extrapolation.
        
        Args:
            fock_matrix: Current Fock matrix
            density_matrix: Current density matrix
            
        Returns:
            Extrapolated Fock matrix
        """
        # Calculate error vector: e = FPS - SPF
        error = fock_matrix @ density_matrix @ self.S - self.S @ density_matrix @ fock_matrix
        
        # Store error and Fock matrix
        self.diis_errors.append(error.flatten())
        self.diis_focks.append(fock_matrix.copy())
        
        # Limit number of stored vectors
        if len(self.diis_errors) > self.max_diis_vectors:
            self.diis_errors.pop(0)
            self.diis_focks.pop(0)
        
        n_vectors = len(self.diis_errors)
        if n_vectors < 2:
            return fock_matrix
        
        # Build B matrix
        B = np.zeros((n_vectors + 1, n_vectors + 1))
        B[-1, :-1] = B[:-1, -1] = -1.0
        
        for i in range(n_vectors):
            for j in range(n_vectors):
                B[i, j] = np.dot(self.diis_errors[i], self.diis_errors[j])
        
        # Solve linear system
        rhs = np.zeros(n_vectors + 1)
        rhs[-1] = -1.0
        
        try:
            coeffs = np.linalg.solve(B, rhs)
            
            # Extrapolate Fock matrix
            F_new = np.zeros_like(fock_matrix)
            for i in range(n_vectors):
                F_new += coeffs[i] * self.diis_focks[i]
            
            return F_new
        
        except np.linalg.LinAlgError:
            # If DIIS fails, return current Fock matrix
            return fock_matrix
    
    def scf_iteration(self, density_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Perform one SCF iteration.
        
        Args:
            density_matrix: Current density matrix
            
        Returns:
            New density matrix, Fock matrix, orbital coefficients, orbital energies
        """
        # Build Fock matrix
        F = self.build_fock_matrix(density_matrix)
        
        # Apply DIIS if requested
        if self.use_diis:
            F = self.diis_extrapolation(F, density_matrix)
        
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
        """
        Run the SCF procedure.
        
        Args:
            initial_guess: Method for initial density guess
            
        Returns:
            SCF results
        """
        print("Starting SCF calculation...")
        print(f"Convergence threshold: {self.convergence}")
        print(f"Maximum iterations: {self.max_iterations}")
        print(f"Number of basis functions: {self.n_basis}")
        print(f"Number of electrons: {self.n_electrons}")
        print(f"Number of occupied orbitals: {self.n_occupied}")
        print()
        
        # Setup integrals
        self.setup_integrals()
        
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
            P_new, F, C, orbital_energies = self.scf_iteration(P)
            
            # Calculate energy
            electronic_energy = self.calculate_energy(P_new, F)
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
        Run Hartree-Fock calculation.
        
        Args:
            initial_guess: Method for initial density guess
            
        Returns:
            SCF results
        """
        return self.run_scf(initial_guess)
