"""
Mathematical utilities for quantum chemistry calculations.
"""

import numpy as np
from typing import Tuple, Union, Optional
import scipy.special
import scipy.linalg
from .constants import PI, SQRT_PI


def factorial2(n: int) -> int:
    """
    Calculate double factorial n!!
    
    Args:
        n: Integer input
        
    Returns:
        Double factorial of n
    """
    if n <= 0:
        return 1
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result


def boys_function(m: int, x: float) -> float:
    """
    Calculate Boys function F_m(x) used in electron repulsion integrals.
    
    Args:
        m: Order of the Boys function
        x: Argument
        
    Returns:
        Value of F_m(x)
    """
    if x < 1e-12:
        return 1.0 / (2 * m + 1)
    
    # Use the incomplete gamma function relation
    # F_m(x) = (1/2) * x^(-m-0.5) * gamma(m+0.5, x)
    return 0.5 * scipy.special.gammainc(m + 0.5, x) * scipy.special.gamma(m + 0.5) / (x ** (m + 0.5))


def hermite_polynomial(n: int, x: float) -> float:
    """
    Calculate physicist's Hermite polynomial H_n(x).
    
    Args:
        n: Order of polynomial
        x: Argument
        
    Returns:
        Value of H_n(x)
    """
    return scipy.special.eval_hermitenorm(n, x) * (2**0.5)**n


def gaussian_product_center(alpha1: float, center1: np.ndarray, 
                          alpha2: float, center2: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Calculate the product of two Gaussian functions parameters.
    
    Args:
        alpha1, alpha2: Exponents of Gaussians
        center1, center2: Centers of Gaussians
        
    Returns:
        Combined exponent and center
    """
    gamma = alpha1 + alpha2
    center = (alpha1 * center1 + alpha2 * center2) / gamma
    return gamma, center


def overlap_1d(alpha1: float, center1: float, alpha2: float, center2: float) -> float:
    """
    Calculate 1D overlap integral between two Gaussian functions.
    
    Args:
        alpha1, alpha2: Exponents
        center1, center2: Centers
        
    Returns:
        Overlap integral value
    """
    gamma = alpha1 + alpha2
    diff = center1 - center2
    prefactor = np.exp(-alpha1 * alpha2 * diff**2 / gamma)
    return prefactor * np.sqrt(PI / gamma)


def kinetic_1d(alpha1: float, center1: float, alpha2: float, center2: float) -> float:
    """
    Calculate 1D kinetic energy integral between two Gaussian functions.
    
    Args:
        alpha1, alpha2: Exponents
        center1, center2: Centers
        
    Returns:
        Kinetic energy integral value
    """
    gamma = alpha1 + alpha2
    diff = center1 - center2
    overlap = overlap_1d(alpha1, center1, alpha2, center2)
    
    kinetic = alpha1 * alpha2 / gamma * (3 - 2 * alpha1 * alpha2 * diff**2 / gamma)
    return kinetic * overlap


def nuclear_attraction_1d(alpha1: float, center1: float, alpha2: float, center2: float, 
                         nuclear_center: float) -> float:
    """
    Calculate 1D nuclear attraction integral.
    
    Args:
        alpha1, alpha2: Exponents
        center1, center2: Centers of basis functions
        nuclear_center: Position of nucleus
        
    Returns:
        Nuclear attraction integral value
    """
    gamma = alpha1 + alpha2
    product_center = (alpha1 * center1 + alpha2 * center2) / gamma
    diff = product_center - nuclear_center
    
    overlap = overlap_1d(alpha1, center1, alpha2, center2)
    boys = boys_function(0, gamma * diff**2)
    
    return -2 * np.sqrt(gamma / PI) * overlap * boys


def normalize_primitive(alpha: float, l: int, m: int, n: int) -> float:
    """
    Calculate normalization constant for a primitive Gaussian function.
    
    Args:
        alpha: Exponent
        l, m, n: Angular momentum quantum numbers
        
    Returns:
        Normalization constant
    """
    L = l + m + n
    numerator = (2 * alpha / PI)**(3/4) * (4 * alpha)**(L/2)
    denominator = np.sqrt(factorial2(2*l - 1) * factorial2(2*m - 1) * factorial2(2*n - 1))
    return numerator / denominator


def diagonalize_symmetric(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonalize a symmetric matrix.
    
    Args:
        matrix: Symmetric matrix to diagonalize
        
    Returns:
        Eigenvalues and eigenvectors
    """
    eigenvalues, eigenvectors = scipy.linalg.eigh(matrix)
    return eigenvalues, eigenvectors


def symmetric_orthogonalization(overlap_matrix: np.ndarray, 
                              threshold: float = 1e-12) -> np.ndarray:
    """
    Perform symmetric orthogonalization.
    
    Args:
        overlap_matrix: Overlap matrix
        threshold: Threshold for linear dependence
        
    Returns:
        Orthogonalization matrix
    """
    eigenvals, eigenvecs = diagonalize_symmetric(overlap_matrix)
    
    # Remove linearly dependent basis functions
    mask = eigenvals > threshold
    eigenvals = eigenvals[mask]
    eigenvecs = eigenvecs[:, mask]
    
    # Form S^(-1/2)
    inv_sqrt_eigenvals = np.diag(1.0 / np.sqrt(eigenvals))
    X = eigenvecs @ inv_sqrt_eigenvals @ eigenvecs.T
    
    return X


def canonical_orthogonalization(overlap_matrix: np.ndarray,
                               threshold: float = 1e-12) -> np.ndarray:
    """
    Perform canonical orthogonalization (Löwdin orthogonalization).
    
    Args:
        overlap_matrix: Overlap matrix
        threshold: Threshold for linear dependence
        
    Returns:
        Orthogonalization matrix
    """
    eigenvals, eigenvecs = diagonalize_symmetric(overlap_matrix)
    
    # Remove linearly dependent basis functions
    mask = eigenvals > threshold
    eigenvals = eigenvals[mask]
    eigenvecs = eigenvecs[:, mask]
    
    # Form X = U * s^(-1/2)
    inv_sqrt_eigenvals = 1.0 / np.sqrt(eigenvals)
    X = eigenvecs * inv_sqrt_eigenvals
    
    return X


def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Calculate commutator [A, B] = AB - BA.
    
    Args:
        A, B: Input matrices
        
    Returns:
        Commutator matrix
    """
    return A @ B - B @ A


def density_matrix(coefficients: np.ndarray, n_occupied: int) -> np.ndarray:
    """
    Build density matrix from molecular orbital coefficients.
    
    Args:
        coefficients: MO coefficient matrix
        n_occupied: Number of occupied orbitals
        
    Returns:
        Density matrix
    """
    C_occ = coefficients[:, :n_occupied]
    return 2.0 * C_occ @ C_occ.T


def mulliken_charges(density_matrix: np.ndarray, overlap_matrix: np.ndarray,
                    nuclear_charges: np.ndarray, basis_on_atom: np.ndarray) -> np.ndarray:
    """
    Calculate Mulliken atomic charges.
    
    Args:
        density_matrix: Electronic density matrix
        overlap_matrix: Overlap matrix
        nuclear_charges: Nuclear charges for each atom
        basis_on_atom: Array mapping basis functions to atoms
        
    Returns:
        Mulliken charges for each atom
    """
    # Calculate Mulliken population matrix
    population_matrix = density_matrix @ overlap_matrix
    
    # Sum populations for each atom
    n_atoms = len(nuclear_charges)
    electron_populations = np.zeros(n_atoms)
    
    for i in range(len(basis_on_atom)):
        atom_idx = basis_on_atom[i]
        electron_populations[atom_idx] += population_matrix[i, i]
    
    # Charges = nuclear charge - electron population
    charges = nuclear_charges - electron_populations
    return charges
