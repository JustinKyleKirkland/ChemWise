"""
Exchange-correlation functionals for DFT calculations.
"""

import numpy as np
from typing import Tuple, Dict, Callable
from abc import ABC, abstractmethod


class XCFunctional(ABC):
    """Abstract base class for exchange-correlation functionals."""
    
    def __init__(self, name: str):
        self.name = name
    
    @abstractmethod
    def exchange_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """Calculate exchange energy density."""
        pass
    
    @abstractmethod
    def correlation_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """Calculate correlation energy density."""
        pass
    
    @abstractmethod
    def exchange_potential(self, rho: np.ndarray) -> np.ndarray:
        """Calculate exchange potential."""
        pass
    
    @abstractmethod
    def correlation_potential(self, rho: np.ndarray) -> np.ndarray:
        """Calculate correlation potential."""
        pass
    
    def xc_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """Calculate total XC energy density."""
        return self.exchange_energy_density(rho) + self.correlation_energy_density(rho)
    
    def xc_potential(self, rho: np.ndarray) -> np.ndarray:
        """Calculate total XC potential."""
        return self.exchange_potential(rho) + self.correlation_potential(rho)


class LDAFunctional(XCFunctional):
    """Local Density Approximation (LDA) functional."""
    
    def __init__(self):
        super().__init__("LDA")
        # Constants for LDA
        self.C_x = -(3.0/4.0) * (3.0/np.pi)**(1.0/3.0)  # Exchange coefficient
        
        # Vosko-Wilk-Nusair correlation parameters
        self.A = 0.0621814
        self.x0 = -0.10498
        self.b = 3.72744
        self.c = 12.9352
        self.Q = np.sqrt(4*self.c - self.b**2)
    
    def exchange_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """LDA exchange energy density."""
        # ε_x = -C_x * ρ^(4/3)
        return self.C_x * rho**(4.0/3.0)
    
    def correlation_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """VWN correlation energy density."""
        rs = (3.0 / (4.0 * np.pi * rho))**(1.0/3.0)
        x = np.sqrt(rs)
        
        X_x = x**2 + self.b*x + self.c
        X_0 = self.x0**2 + self.b*self.x0 + self.c
        
        # Avoid division by zero
        safe_X_x = np.where(X_x != 0, X_x, 1e-16)
        safe_x_minus_x0 = np.where(x != self.x0, x - self.x0, 1e-16)
        
        ln_term = np.log(x**2 / safe_X_x)
        atan_term = np.arctan(self.Q / (2*x + self.b)) - np.arctan(self.Q / (2*self.x0 + self.b))
        ln_term2 = np.log(safe_x_minus_x0**2 / safe_X_x) * self.b*self.x0 / X_0
        
        eps_c = self.A * (ln_term + 2*self.b/self.Q * atan_term - ln_term2)
        
        return eps_c
    
    def exchange_potential(self, rho: np.ndarray) -> np.ndarray:
        """LDA exchange potential."""
        # V_x = d(ρ*ε_x)/dρ = (4/3) * ε_x
        return (4.0/3.0) * self.exchange_energy_density(rho)
    
    def correlation_potential(self, rho: np.ndarray) -> np.ndarray:
        """VWN correlation potential."""
        rs = (3.0 / (4.0 * np.pi * rho))**(1.0/3.0)
        x = np.sqrt(rs)
        
        # This is a simplified version - full VWN potential is more complex
        eps_c = self.correlation_energy_density(rho)
        
        # Approximate derivative (should be analytical for production code)
        drho = 1e-8
        rho_plus = rho + drho
        eps_c_plus = self.correlation_energy_density(rho_plus)
        
        deps_c_drho = (eps_c_plus - eps_c) / drho
        V_c = eps_c + rho * deps_c_drho
        
        return V_c


class B88Exchange(XCFunctional):
    """Becke 1988 exchange functional."""
    
    def __init__(self):
        super().__init__("B88")
        self.beta = 0.0042
        self.C_x = -(3.0/4.0) * (3.0/np.pi)**(1.0/3.0)
    
    def exchange_energy_density(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B88 exchange energy density."""
        if grad_rho is None:
            # Fall back to LDA
            return self.C_x * rho**(4.0/3.0)
        
        # Reduced gradient
        s = np.sqrt(np.sum(grad_rho**2, axis=-1)) / (2 * (3*np.pi**2)**(1.0/3.0) * rho**(4.0/3.0))
        
        # B88 enhancement factor
        x = 6 * self.beta * s * np.arcsinh(s)
        F_x = 1 + x / (1 + 6*self.beta*s*x**(1.0/3.0))
        
        # Exchange energy density
        eps_x_lda = self.C_x * rho**(1.0/3.0)
        return eps_x_lda * F_x
    
    def correlation_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """No correlation in B88."""
        return np.zeros_like(rho)
    
    def exchange_potential(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B88 exchange potential (simplified)."""
        if grad_rho is None:
            return (4.0/3.0) * self.C_x * rho**(1.0/3.0)
        
        # Simplified - should include gradient terms
        s = np.sqrt(np.sum(grad_rho**2, axis=-1)) / (2 * (3*np.pi**2)**(1.0/3.0) * rho**(4.0/3.0))
        F_x = 1 + 6*self.beta*s / (1 + 6*self.beta*s)  # Simplified
        
        V_x_lda = (4.0/3.0) * self.C_x * rho**(1.0/3.0)
        return V_x_lda * F_x
    
    def correlation_potential(self, rho: np.ndarray) -> np.ndarray:
        """No correlation potential in B88."""
        return np.zeros_like(rho)


class LYPCorrelation(XCFunctional):
    """Lee-Yang-Parr correlation functional."""
    
    def __init__(self):
        super().__init__("LYP")
        self.a = 0.04918
        self.b = 0.132
        self.c = 0.2533
        self.d = 0.349
    
    def exchange_energy_density(self, rho: np.ndarray) -> np.ndarray:
        """No exchange in LYP."""
        return np.zeros_like(rho)
    
    def correlation_energy_density(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """LYP correlation energy density (simplified)."""
        if grad_rho is None:
            # Fall back to LDA correlation
            lda = LDAFunctional()
            return lda.correlation_energy_density(rho)
        
        # Simplified LYP - full implementation requires spin densities
        gamma = np.sum(grad_rho**2, axis=-1)
        
        # LYP parameters
        omega = np.exp(-self.c * rho**(-1.0/3.0)) / (1 + self.d * rho**(-1.0/3.0))
        delta = self.c * rho**(-1.0/3.0) + self.d * rho**(-1.0/3.0) / (1 + self.d * rho**(-1.0/3.0))
        
        # Simplified correlation energy
        eps_c = -self.a * omega * (1 + self.b * gamma / (rho**(8.0/3.0) * delta))
        
        return eps_c
    
    def exchange_potential(self, rho: np.ndarray) -> np.ndarray:
        """No exchange potential in LYP."""
        return np.zeros_like(rho)
    
    def correlation_potential(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """LYP correlation potential (simplified)."""
        if grad_rho is None:
            lda = LDAFunctional()
            return lda.correlation_potential(rho)
        
        # Simplified - should be analytical derivative
        eps_c = self.correlation_energy_density(rho, grad_rho)
        
        # Approximate derivative
        drho = 1e-8
        rho_plus = rho + drho
        eps_c_plus = self.correlation_energy_density(rho_plus, grad_rho)
        
        deps_c_drho = (eps_c_plus - eps_c) / drho
        return eps_c + rho * deps_c_drho


class B3LYPFunctional(XCFunctional):
    """B3LYP hybrid functional."""
    
    def __init__(self):
        super().__init__("B3LYP")
        self.a0 = 0.20  # Exact exchange coefficient
        self.ax = 0.72  # B88 exchange coefficient
        self.ac = 0.81  # LYP correlation coefficient
        
        # Component functionals
        self.lda = LDAFunctional()
        self.b88 = B88Exchange()
        self.lyp = LYPCorrelation()
    
    def exchange_energy_density(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B3LYP exchange energy density."""
        eps_x_lda = self.lda.exchange_energy_density(rho)
        
        if grad_rho is not None:
            eps_x_b88 = self.b88.exchange_energy_density(rho, grad_rho)
            # B3LYP: E_x = (1-a0)*E_x^LDA + a0*E_x^exact + ax*(E_x^B88 - E_x^LDA)
            # Here we only compute the DFT part (exact exchange handled separately)
            return (1 - self.a0) * eps_x_lda + self.ax * (eps_x_b88 - eps_x_lda)
        else:
            return (1 - self.a0) * eps_x_lda
    
    def correlation_energy_density(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B3LYP correlation energy density."""
        eps_c_lda = self.lda.correlation_energy_density(rho)
        
        if grad_rho is not None:
            eps_c_lyp = self.lyp.correlation_energy_density(rho, grad_rho)
            return (1 - self.ac) * eps_c_lda + self.ac * eps_c_lyp
        else:
            return eps_c_lda
    
    def exchange_potential(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B3LYP exchange potential."""
        V_x_lda = self.lda.exchange_potential(rho)
        
        if grad_rho is not None:
            V_x_b88 = self.b88.exchange_potential(rho, grad_rho)
            return (1 - self.a0) * V_x_lda + self.ax * (V_x_b88 - V_x_lda)
        else:
            return (1 - self.a0) * V_x_lda
    
    def correlation_potential(self, rho: np.ndarray, grad_rho: np.ndarray = None) -> np.ndarray:
        """B3LYP correlation potential."""
        V_c_lda = self.lda.correlation_potential(rho)
        
        if grad_rho is not None:
            V_c_lyp = self.lyp.correlation_potential(rho, grad_rho)
            return (1 - self.ac) * V_c_lda + self.ac * V_c_lyp
        else:
            return V_c_lda
    
    @property
    def exact_exchange_fraction(self) -> float:
        """Get exact exchange fraction for hybrid functionals."""
        return self.a0


# Registry of available functionals
FUNCTIONALS: Dict[str, Callable[[], XCFunctional]] = {
    'lda': lambda: LDAFunctional(),
    'b88': lambda: B88Exchange(),
    'lyp': lambda: LYPCorrelation(),
    'b3lyp': lambda: B3LYPFunctional(),
}


def get_functional(name: str) -> XCFunctional:
    """
    Get a functional by name.
    
    Args:
        name: Functional name
        
    Returns:
        XC functional object
    """
    name_lower = name.lower()
    if name_lower in FUNCTIONALS:
        return FUNCTIONALS[name_lower]()
    else:
        raise ValueError(f"Unknown functional: {name}")


def list_functionals() -> list:
    """List available functionals."""
    return list(FUNCTIONALS.keys())
