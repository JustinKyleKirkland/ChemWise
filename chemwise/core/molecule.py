"""
Molecular structure and properties for quantum chemistry calculations.
"""

import numpy as np
from typing import List, Dict, Tuple, Union, Optional
from pathlib import Path

from ..utils.constants import (
    ATOMIC_NUMBERS, ATOMIC_SYMBOLS, ATOMIC_MASSES, COVALENT_RADII,
    ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM
)
from ..utils.io import read_xyz_file, parse_geometry_string


class Atom:
    """Represents an individual atom with its properties."""
    
    def __init__(self, symbol: str, coordinates: np.ndarray, charge: int = None):
        """
        Initialize an atom.
        
        Args:
            symbol: Atomic symbol (e.g., 'H', 'C', 'O')
            coordinates: Cartesian coordinates [x, y, z]
            charge: Nuclear charge (defaults to atomic number)
        """
        self.symbol = symbol.upper()
        self.coordinates = np.array(coordinates, dtype=float)
        self.atomic_number = ATOMIC_NUMBERS.get(self.symbol, 0)
        self.charge = charge if charge is not None else self.atomic_number
        self.mass = ATOMIC_MASSES.get(self.atomic_number, 0.0)
        self.covalent_radius = COVALENT_RADII.get(self.atomic_number, 1.0)
    
    def __str__(self) -> str:
        return f"{self.symbol} ({self.coordinates[0]:.6f}, {self.coordinates[1]:.6f}, {self.coordinates[2]:.6f})"
    
    def __repr__(self) -> str:
        return f"Atom('{self.symbol}', {self.coordinates.tolist()})"
    
    def distance_to(self, other: 'Atom') -> float:
        """Calculate distance to another atom."""
        return np.linalg.norm(self.coordinates - other.coordinates)


class Molecule:
    """Represents a molecular system for quantum chemistry calculations."""
    
    def __init__(self, atoms: List[Atom] = None, charge: int = 0, multiplicity: int = 1,
                 units: str = 'angstrom'):
        """
        Initialize a molecule.
        
        Args:
            atoms: List of Atom objects
            charge: Total molecular charge
            multiplicity: Spin multiplicity (2S + 1)
            units: Coordinate units ('angstrom' or 'bohr')
        """
        self.atoms = atoms if atoms is not None else []
        self.charge = charge
        self.multiplicity = multiplicity
        self.units = units.lower()
        
        # Derived properties
        self._nuclear_repulsion = None
        self._center_of_mass = None
        self._moments_of_inertia = None
    
    @classmethod
    def from_xyz(cls, filename: Union[str, Path], charge: int = 0, 
                 multiplicity: int = 1) -> 'Molecule':
        """
        Create molecule from XYZ file.
        
        Args:
            filename: Path to XYZ file
            charge: Molecular charge
            multiplicity: Spin multiplicity
            
        Returns:
            Molecule object
        """
        data = read_xyz_file(filename)
        atoms = []
        
        for symbol, coord in zip(data['symbols'], data['coordinates']):
            atoms.append(Atom(symbol, coord))
        
        return cls(atoms=atoms, charge=charge, multiplicity=multiplicity,
                  units=data.get('units', 'angstrom'))
    
    @classmethod
    def from_string(cls, geometry_str: str, charge: int = 0, 
                   multiplicity: int = 1) -> 'Molecule':
        """
        Create molecule from geometry string.
        
        Args:
            geometry_str: Multi-line string with atomic symbols and coordinates
            charge: Molecular charge
            multiplicity: Spin multiplicity
            
        Returns:
            Molecule object
        """
        data = parse_geometry_string(geometry_str)
        atoms = []
        
        for symbol, coord in zip(data['symbols'], data['coordinates']):
            atoms.append(Atom(symbol, coord))
        
        return cls(atoms=atoms, charge=charge, multiplicity=multiplicity,
                  units=data.get('units', 'angstrom'))
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'Molecule':
        """
        Create molecule from dictionary (e.g., from YAML input).
        
        Args:
            data: Dictionary containing molecular data
            
        Returns:
            Molecule object
        """
        atoms = []
        
        # Handle direct XYZ format (symbols and coordinates arrays)
        if 'symbols' in data and 'coordinates' in data:
            for symbol, coord in zip(data['symbols'], data['coordinates']):
                atoms.append(Atom(symbol, coord))
        
        # Handle nested geometry format
        elif 'geometry' in data:
            # Handle string geometry
            if isinstance(data['geometry'], str):
                geom_data = parse_geometry_string(data['geometry'])
                for symbol, coord in zip(geom_data['symbols'], geom_data['coordinates']):
                    atoms.append(Atom(symbol, coord))
            # Handle list of atoms
            elif isinstance(data['geometry'], list):
                for atom_data in data['geometry']:
                    if isinstance(atom_data, dict):
                        symbol = atom_data['symbol']
                        coord = atom_data['coordinates']
                    else:
                        # Assume [symbol, x, y, z] format
                        symbol, coord = atom_data[0], atom_data[1:4]
                    atoms.append(Atom(symbol, coord))
        
        return cls(
            atoms=atoms,
            charge=data.get('charge', 0),
            multiplicity=data.get('multiplicity', 1),
            units=data.get('units', 'angstrom')
        )
    
    def add_atom(self, symbol: str, coordinates: np.ndarray) -> None:
        """Add an atom to the molecule."""
        self.atoms.append(Atom(symbol, coordinates))
        self._clear_cache()
    
    def remove_atom(self, index: int) -> None:
        """Remove an atom from the molecule."""
        if 0 <= index < len(self.atoms):
            del self.atoms[index]
            self._clear_cache()
    
    def _clear_cache(self) -> None:
        """Clear cached properties."""
        self._nuclear_repulsion = None
        self._center_of_mass = None
        self._moments_of_inertia = None
    
    @property
    def n_atoms(self) -> int:
        """Number of atoms in the molecule."""
        return len(self.atoms)
    
    @property
    def symbols(self) -> List[str]:
        """List of atomic symbols."""
        return [atom.symbol for atom in self.atoms]
    
    @property
    def atomic_numbers(self) -> np.ndarray:
        """Array of atomic numbers."""
        return np.array([atom.atomic_number for atom in self.atoms])
    
    @property
    def nuclear_charges(self) -> np.ndarray:
        """Array of nuclear charges."""
        return np.array([atom.charge for atom in self.atoms])
    
    @property
    def coordinates(self) -> np.ndarray:
        """Atomic coordinates as numpy array."""
        return np.array([atom.coordinates for atom in self.atoms])
    
    @property
    def coordinates_bohr(self) -> np.ndarray:
        """Atomic coordinates in Bohr."""
        coords = self.coordinates
        if self.units == 'angstrom':
            coords = coords * ANGSTROM_TO_BOHR
        return coords
    
    @property
    def coordinates_angstrom(self) -> np.ndarray:
        """Atomic coordinates in Angstrom."""
        coords = self.coordinates
        if self.units == 'bohr':
            coords = coords * BOHR_TO_ANGSTROM
        return coords
    
    @property
    def masses(self) -> np.ndarray:
        """Array of atomic masses."""
        return np.array([atom.mass for atom in self.atoms])
    
    @property
    def n_electrons(self) -> int:
        """Total number of electrons."""
        return int(np.sum(self.atomic_numbers) - self.charge)
    
    @property
    def n_electron_pairs(self) -> int:
        """Number of electron pairs (for closed-shell systems)."""
        return self.n_electrons // 2
    
    @property
    def nuclear_repulsion_energy(self) -> float:
        """Calculate nuclear repulsion energy."""
        if self._nuclear_repulsion is None:
            energy = 0.0
            coords = self.coordinates_bohr  # Use Bohr for atomic units
            
            for i in range(self.n_atoms):
                for j in range(i + 1, self.n_atoms):
                    distance = np.linalg.norm(coords[i] - coords[j])
                    energy += (self.nuclear_charges[i] * self.nuclear_charges[j]) / distance
            
            self._nuclear_repulsion = energy
        
        return self._nuclear_repulsion
    
    @property
    def center_of_mass(self) -> np.ndarray:
        """Calculate center of mass."""
        if self._center_of_mass is None:
            total_mass = np.sum(self.masses)
            weighted_coords = self.coordinates.T * self.masses
            self._center_of_mass = np.sum(weighted_coords, axis=1) / total_mass
        
        return self._center_of_mass
    
    @property
    def moments_of_inertia(self) -> np.ndarray:
        """Calculate moments of inertia tensor."""
        if self._moments_of_inertia is None:
            com = self.center_of_mass
            coords = self.coordinates - com
            masses = self.masses
            
            Ixx = np.sum(masses * (coords[:, 1]**2 + coords[:, 2]**2))
            Iyy = np.sum(masses * (coords[:, 0]**2 + coords[:, 2]**2))
            Izz = np.sum(masses * (coords[:, 0]**2 + coords[:, 1]**2))
            Ixy = -np.sum(masses * coords[:, 0] * coords[:, 1])
            Ixz = -np.sum(masses * coords[:, 0] * coords[:, 2])
            Iyz = -np.sum(masses * coords[:, 1] * coords[:, 2])
            
            self._moments_of_inertia = np.array([
                [Ixx, Ixy, Ixz],
                [Ixy, Iyy, Iyz],
                [Ixz, Iyz, Izz]
            ])
        
        return self._moments_of_inertia
    
    def distance_matrix(self) -> np.ndarray:
        """Calculate distance matrix between all atoms."""
        coords = self.coordinates
        n = self.n_atoms
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                dist = np.linalg.norm(coords[i] - coords[j])
                distances[i, j] = distances[j, i] = dist
        
        return distances
    
    def bond_connectivity(self, scale_factor: float = 1.2) -> np.ndarray:
        """
        Determine bond connectivity based on covalent radii.
        
        Args:
            scale_factor: Factor to scale covalent radii sum
            
        Returns:
            Boolean connectivity matrix
        """
        distances = self.distance_matrix()
        connectivity = np.zeros((self.n_atoms, self.n_atoms), dtype=bool)
        
        for i in range(self.n_atoms):
            for j in range(i + 1, self.n_atoms):
                sum_radii = (self.atoms[i].covalent_radius + 
                           self.atoms[j].covalent_radius) * scale_factor
                
                if self.units == 'bohr':
                    sum_radii *= ANGSTROM_TO_BOHR
                
                if distances[i, j] <= sum_radii:
                    connectivity[i, j] = connectivity[j, i] = True
        
        return connectivity
    
    def translate(self, vector: np.ndarray) -> None:
        """Translate molecule by given vector."""
        for atom in self.atoms:
            atom.coordinates += vector
        self._clear_cache()
    
    def rotate(self, rotation_matrix: np.ndarray) -> None:
        """Rotate molecule using rotation matrix."""
        for atom in self.atoms:
            atom.coordinates = rotation_matrix @ atom.coordinates
        self._clear_cache()
    
    def center_at_origin(self) -> None:
        """Center molecule at origin (center of mass)."""
        com = self.center_of_mass
        self.translate(-com)
    
    def copy(self) -> 'Molecule':
        """Create a deep copy of the molecule."""
        new_atoms = []
        for atom in self.atoms:
            new_atoms.append(Atom(atom.symbol, atom.coordinates.copy(), atom.charge))
        
        return Molecule(
            atoms=new_atoms,
            charge=self.charge,
            multiplicity=self.multiplicity,
            units=self.units
        )
    
    def to_dict(self) -> Dict:
        """Convert molecule to dictionary format."""
        geometry_list = []
        for atom in self.atoms:
            geometry_list.append({
                'symbol': atom.symbol,
                'coordinates': atom.coordinates.tolist()
            })
        
        return {
            'geometry': geometry_list,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'units': self.units
        }
    
    def __str__(self) -> str:
        """String representation of molecule."""
        lines = [f"Molecule (charge={self.charge}, multiplicity={self.multiplicity})"]
        lines.append(f"Units: {self.units}")
        lines.append("-" * 50)
        
        for i, atom in enumerate(self.atoms):
            lines.append(f"{i+1:3d}: {atom}")
        
        lines.append("-" * 50)
        lines.append(f"Nuclear repulsion energy: {self.nuclear_repulsion_energy:.8f} Hartree")
        
        return "\n".join(lines)
    
    def __repr__(self) -> str:
        return f"Molecule(n_atoms={self.n_atoms}, charge={self.charge}, multiplicity={self.multiplicity})"
