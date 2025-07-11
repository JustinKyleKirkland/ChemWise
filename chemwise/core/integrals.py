"""
Integral evaluation for quantum chemistry calculations.
"""

import numpy as np
from typing import Tuple
from ..utils.math_utils import boys_function, factorial2, overlap_1d, kinetic_1d
from ..utils.constants import PI
from .basis import BasisSet


class IntegralEngine:
    """Engine for evaluating quantum chemistry integrals."""
    
    def __init__(self, basis_set: BasisSet):
        """
        Initialize integral engine.
        
        Args:
            basis_set: Basis set for the calculation
        """
        self.basis_set = basis_set
        self.n_basis = basis_set.n_basis
    
    def overlap_matrix(self) -> np.ndarray:
        """
        Calculate overlap matrix S_ij = <φ_i|φ_j>.
        
        Returns:
            Overlap matrix
        """
        S = np.zeros((self.n_basis, self.n_basis))
        
        for i, bf_i in enumerate(self.basis_set.basis_functions):
            for j, bf_j in enumerate(self.basis_set.basis_functions):
                S[i, j] = self._overlap_integral(bf_i, bf_j)
        
        return S
    
    def kinetic_matrix(self) -> np.ndarray:
        """
        Calculate kinetic energy matrix T_ij = <φ_i|∇²|φ_j>.
        
        Returns:
            Kinetic energy matrix
        """
        T = np.zeros((self.n_basis, self.n_basis))
        
        for i, bf_i in enumerate(self.basis_set.basis_functions):
            for j, bf_j in enumerate(self.basis_set.basis_functions):
                T[i, j] = self._kinetic_integral(bf_i, bf_j)
        
        return T
    
    def nuclear_attraction_matrix(self, nuclear_charges: np.ndarray, 
                                nuclear_coords: np.ndarray) -> np.ndarray:
        """
        Calculate nuclear attraction matrix V_ij = <φ_i|∑_A Z_A/r_iA|φ_j>.
        
        Args:
            nuclear_charges: Array of nuclear charges
            nuclear_coords: Array of nuclear coordinates
            
        Returns:
            Nuclear attraction matrix
        """
        V = np.zeros((self.n_basis, self.n_basis))
        
        for i, bf_i in enumerate(self.basis_set.basis_functions):
            for j, bf_j in enumerate(self.basis_set.basis_functions):
                for k, (charge, coord) in enumerate(zip(nuclear_charges, nuclear_coords)):
                    V[i, j] += charge * self._nuclear_attraction_integral(bf_i, bf_j, coord)
        
        return -V  # Negative because of attractive interaction
    
    def electron_repulsion_tensor(self) -> np.ndarray:
        """
        Calculate electron repulsion integrals (ij|kl).
        
        Returns:
            4D tensor of electron repulsion integrals
        """
        eri = np.zeros((self.n_basis, self.n_basis, self.n_basis, self.n_basis))
        
        for i, bf_i in enumerate(self.basis_set.basis_functions):
            for j, bf_j in enumerate(self.basis_set.basis_functions):
                for k, bf_k in enumerate(self.basis_set.basis_functions):
                    for ll, bf_l in enumerate(self.basis_set.basis_functions):
                        eri[i, j, k, ll] = self._electron_repulsion_integral(bf_i, bf_j, bf_k, bf_l)
        
        return eri
    
    def _overlap_integral(self, bf_i, bf_j) -> float:
        """Calculate overlap integral between two basis functions."""
        total = 0.0
        
        for prim_i in bf_i.contracted_gaussian.primitives:
            for prim_j in bf_j.contracted_gaussian.primitives:
                # Gaussian product theorem
                alpha_i, alpha_j = prim_i.exponent, prim_j.exponent
                center_i, center_j = bf_i.center, bf_j.center
                
                gamma = alpha_i + alpha_j
                diff = center_i - center_j
                
                # Gaussian overlap factor
                gaussian_overlap = np.exp(-alpha_i * alpha_j * np.dot(diff, diff) / gamma)
                gaussian_overlap *= (PI / gamma)**(3/2)
                
                # Angular momentum overlap
                angular_overlap = self._angular_overlap(
                    bf_i.angular_momentum, bf_j.angular_momentum,
                    alpha_i, alpha_j, center_i, center_j
                )
                
                total += (prim_i.coefficient * prim_j.coefficient * 
                         gaussian_overlap * angular_overlap)
        
        return total
    
    def _kinetic_integral(self, bf_i, bf_j) -> float:
        """Calculate kinetic energy integral between two basis functions."""
        total = 0.0
        
        for prim_i in bf_i.contracted_gaussian.primitives:
            for prim_j in bf_j.contracted_gaussian.primitives:
                alpha_i, alpha_j = prim_i.exponent, prim_j.exponent
                center_i, center_j = bf_i.center, bf_j.center
                
                # Use the kinetic energy formula for Gaussian functions
                kinetic = self._gaussian_kinetic_integral(
                    alpha_i, alpha_j, center_i, center_j,
                    bf_i.angular_momentum, bf_j.angular_momentum
                )
                
                total += prim_i.coefficient * prim_j.coefficient * kinetic
        
        return total
    
    def _nuclear_attraction_integral(self, bf_i, bf_j, nuclear_coord: np.ndarray) -> float:
        """Calculate nuclear attraction integral."""
        total = 0.0
        
        for prim_i in bf_i.contracted_gaussian.primitives:
            for prim_j in bf_j.contracted_gaussian.primitives:
                alpha_i, alpha_j = prim_i.exponent, prim_j.exponent
                center_i, center_j = bf_i.center, bf_j.center
                
                nuclear_integral = self._gaussian_nuclear_integral(
                    alpha_i, alpha_j, center_i, center_j,
                    bf_i.angular_momentum, bf_j.angular_momentum,
                    nuclear_coord
                )
                
                total += prim_i.coefficient * prim_j.coefficient * nuclear_integral
        
        return total
    
    def _electron_repulsion_integral(self, bf_i, bf_j, bf_k, bf_l) -> float:
        """Calculate electron repulsion integral (ij|kl)."""
        total = 0.0
        
        for prim_i in bf_i.contracted_gaussian.primitives:
            for prim_j in bf_j.contracted_gaussian.primitives:
                for prim_k in bf_k.contracted_gaussian.primitives:
                    for prim_l in bf_l.contracted_gaussian.primitives:
                        
                        alpha_i, alpha_j = prim_i.exponent, prim_j.exponent
                        alpha_k, alpha_l = prim_k.exponent, prim_l.exponent
                        
                        eri = self._gaussian_eri(
                            alpha_i, alpha_j, alpha_k, alpha_l,
                            bf_i.center, bf_j.center, bf_k.center, bf_l.center,
                            bf_i.angular_momentum, bf_j.angular_momentum,
                            bf_k.angular_momentum, bf_l.angular_momentum
                        )
                        
                        total += (prim_i.coefficient * prim_j.coefficient * 
                                prim_k.coefficient * prim_l.coefficient * eri)
        
        return total
    
    def _angular_overlap(self, l1: Tuple[int, int, int], l2: Tuple[int, int, int],
                        alpha1: float, alpha2: float, 
                        center1: np.ndarray, center2: np.ndarray) -> float:
        """Calculate angular overlap factor."""
        gamma = alpha1 + alpha2
        product_center = (alpha1 * center1 + alpha2 * center2) / gamma
        
        # 1D overlaps for each direction
        overlap_x = self._overlap_1d(l1[0], l2[0], alpha1, alpha2, 
                                    center1[0], center2[0], product_center[0])
        overlap_y = self._overlap_1d(l1[1], l2[1], alpha1, alpha2,
                                    center1[1], center2[1], product_center[1])
        overlap_z = self._overlap_1d(l1[2], l2[2], alpha1, alpha2,
                                    center1[2], center2[2], product_center[2])
        
        return overlap_x * overlap_y * overlap_z
    
    def _overlap_1d(self, l1: int, l2: int, alpha1: float, alpha2: float,
                   center1: float, center2: float, product_center: float) -> float:
        """Calculate 1D overlap integral."""
        if (l1 + l2) % 2 == 1:
            return 0.0  # Odd angular momentum sum gives zero overlap
        
        gamma = alpha1 + alpha2
        diff1 = product_center - center1
        diff2 = product_center - center2
        
        # Use recursion relation for angular momentum integrals
        overlap = 0.0
        for i in range((l1 + l2) // 2 + 1):
            exp1 = l1 - 2*i
            exp2 = l2 - 2*i
            
            # Handle cases where exponents might be negative and diff might be zero
            if exp1 < 0 and abs(diff1) < 1e-12:
                continue  # Skip this term
            if exp2 < 0 and abs(diff2) < 1e-12:
                continue  # Skip this term
            
            factor = self._binomial_coefficient(l1 + l2, 2 * i)
            factor *= diff1**exp1 * diff2**exp2
            factor *= factorial2(2*i - 1) / (2 * gamma)**i
            overlap += factor
        
        return overlap
    
    def _gaussian_kinetic_integral(self, alpha1: float, alpha2: float,
                                  center1: np.ndarray, center2: np.ndarray,
                                  l1: Tuple[int, int, int], l2: Tuple[int, int, int]) -> float:
        """Calculate kinetic energy integral between two Gaussians."""
        # Simplified implementation - in practice this would be more complex
        gamma = alpha1 + alpha2
        diff = center1 - center2
        
        # Gaussian overlap part
        gaussian_overlap = np.exp(-alpha1 * alpha2 * np.dot(diff, diff) / gamma)
        gaussian_overlap *= (PI / gamma)**(3/2)
        
        # Kinetic energy factor (simplified)
        L1, L2 = sum(l1), sum(l2)
        kinetic_factor = alpha1 * alpha2 / gamma * (3 + L1 + L2)
        
        return kinetic_factor * gaussian_overlap
    
    def _gaussian_nuclear_integral(self, alpha1: float, alpha2: float,
                                  center1: np.ndarray, center2: np.ndarray,
                                  l1: Tuple[int, int, int], l2: Tuple[int, int, int],
                                  nuclear_coord: np.ndarray) -> float:
        """Calculate nuclear attraction integral."""
        gamma = alpha1 + alpha2
        product_center = (alpha1 * center1 + alpha2 * center2) / gamma
        diff = center1 - center2
        
        # Gaussian overlap part
        gaussian_overlap = np.exp(-alpha1 * alpha2 * np.dot(diff, diff) / gamma)
        gaussian_overlap *= (PI / gamma)**(3/2)
        
        # Nuclear attraction using Boys function
        pc_diff = product_center - nuclear_coord
        T = gamma * np.dot(pc_diff, pc_diff)
        boys_val = boys_function(sum(l1) + sum(l2), T)
        
        return gaussian_overlap * boys_val * 2 * PI / gamma
    
    def _gaussian_eri(self, alpha1: float, alpha2: float, alpha3: float, alpha4: float,
                     center1: np.ndarray, center2: np.ndarray, 
                     center3: np.ndarray, center4: np.ndarray,
                     l1: Tuple[int, int, int], l2: Tuple[int, int, int],
                     l3: Tuple[int, int, int], l4: Tuple[int, int, int]) -> float:
        """Calculate electron repulsion integral between four Gaussians."""
        # Simplified implementation
        gamma12 = alpha1 + alpha2
        gamma34 = alpha3 + alpha4
        
        product_center12 = (alpha1 * center1 + alpha2 * center2) / gamma12
        product_center34 = (alpha3 * center3 + alpha4 * center4) / gamma34
        
        diff12 = center1 - center2
        diff34 = center3 - center4
        diff_centers = product_center12 - product_center34
        
        # Gaussian overlap factors
        overlap12 = np.exp(-alpha1 * alpha2 * np.dot(diff12, diff12) / gamma12)
        overlap34 = np.exp(-alpha3 * alpha4 * np.dot(diff34, diff34) / gamma34)
        
        # Boys function argument
        delta = gamma12 * gamma34 / (gamma12 + gamma34)
        T = delta * np.dot(diff_centers, diff_centers)
        
        # Total angular momentum
        L_total = sum(l1) + sum(l2) + sum(l3) + sum(l4)
        boys_val = boys_function(L_total, T)
        
        # Prefactor
        prefactor = 2 * PI**(5/2) / (gamma12 * gamma34 * np.sqrt(gamma12 + gamma34))
        
        return prefactor * overlap12 * overlap34 * boys_val
    
    def _binomial_coefficient(self, n: int, k: int) -> int:
        """Calculate binomial coefficient C(n,k)."""
        if k > n or k < 0:
            return 0
        if k == 0 or k == n:
            return 1
        
        # Use multiplicative formula
        result = 1
        for i in range(min(k, n - k)):
            result = result * (n - i) // (i + 1)
        
        return result
