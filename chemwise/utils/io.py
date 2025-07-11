"""
Input/output utilities for ChemWise quantum chemistry software.
"""

import yaml
import numpy as np
import h5py
from typing import Dict, Any, List, Optional, Union
from pathlib import Path
import re

from .constants import ATOMIC_NUMBERS, BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR


def read_xyz_file(filename: Union[str, Path]) -> Dict[str, Any]:
    """
    Read molecule geometry from XYZ file.
    
    Args:
        filename: Path to XYZ file
        
    Returns:
        Dictionary containing atomic symbols, coordinates, and metadata
    """
    filename = Path(filename)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse number of atoms
    n_atoms = int(lines[0].strip())
    
    # Parse comment line (optional)
    comment = lines[1].strip() if len(lines) > 1 else ""
    
    # Parse atomic symbols and coordinates
    symbols = []
    coordinates = []
    
    for i in range(2, 2 + n_atoms):
        parts = lines[i].strip().split()
        symbols.append(parts[0].upper())
        coords = [float(x) for x in parts[1:4]]
        coordinates.append(coords)
    
    return {
        'symbols': symbols,
        'coordinates': np.array(coordinates),
        'comment': comment,
        'units': 'angstrom'  # XYZ files typically use Angstroms
    }


def write_xyz_file(filename: Union[str, Path], symbols: List[str], 
                   coordinates: np.ndarray, comment: str = "") -> None:
    """
    Write molecule geometry to XYZ file.
    
    Args:
        filename: Output file path
        symbols: List of atomic symbols
        coordinates: Atomic coordinates (Angstroms)
        comment: Optional comment line
    """
    filename = Path(filename)
    
    with open(filename, 'w') as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        
        for symbol, coord in zip(symbols, coordinates):
            f.write(f"{symbol:<3s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")


def read_input_file(filename: Union[str, Path]) -> Dict[str, Any]:
    """
    Read ChemWise input file in YAML format.
    
    Args:
        filename: Path to input file
        
    Returns:
        Dictionary containing input parameters
    """
    filename = Path(filename)
    
    with open(filename, 'r') as f:
        input_data = yaml.safe_load(f)
    
    # Parse molecule geometry if provided as string
    if 'molecule' in input_data and 'geometry' in input_data['molecule']:
        geometry_str = input_data['molecule']['geometry']
        if isinstance(geometry_str, str):
            input_data['molecule'] = parse_geometry_string(geometry_str)
    
    return input_data


def parse_geometry_string(geometry_str: str) -> Dict[str, Any]:
    """
    Parse geometry from string format.
    
    Args:
        geometry_str: Multi-line string with atomic symbols and coordinates
        
    Returns:
        Dictionary with parsed geometry
    """
    lines = geometry_str.strip().split('\n')
    symbols = []
    coordinates = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        if len(parts) >= 4:
            symbols.append(parts[0].upper())
            coords = [float(x) for x in parts[1:4]]
            coordinates.append(coords)
    
    return {
        'symbols': symbols,
        'coordinates': np.array(coordinates),
        'units': 'angstrom'
    }


def write_input_file(filename: Union[str, Path], input_data: Dict[str, Any]) -> None:
    """
    Write ChemWise input file in YAML format.
    
    Args:
        filename: Output file path
        input_data: Input parameters dictionary
    """
    filename = Path(filename)
    
    with open(filename, 'w') as f:
        yaml.dump(input_data, f, default_flow_style=False, indent=2)


def load_molecule(filename: Union[str, Path]) -> Dict[str, Any]:
    """
    Load molecule from various file formats.
    
    Args:
        filename: Path to molecule file
        
    Returns:
        Dictionary containing molecular data
    """
    filename = Path(filename)
    suffix = filename.suffix.lower()
    
    if suffix == '.xyz':
        return read_xyz_file(filename)
    elif suffix in ['.yml', '.yaml']:
        input_data = read_input_file(filename)
        return input_data.get('molecule', {})
    else:
        raise ValueError(f"Unsupported file format: {suffix}")


def save_results(filename: Union[str, Path], results: Dict[str, Any]) -> None:
    """
    Save calculation results to HDF5 file.
    
    Args:
        filename: Output file path
        results: Results dictionary
    """
    filename = Path(filename)
    
    def write_to_hdf5(group, data):
        """Recursively write data to HDF5 group."""
        for key, value in data.items():
            if isinstance(value, dict):
                subgroup = group.create_group(key)
                write_to_hdf5(subgroup, value)
            elif isinstance(value, np.ndarray):
                group.create_dataset(key, data=value)
            elif isinstance(value, (int, float, str)):
                group.attrs[key] = value
            elif isinstance(value, list):
                # Convert list to numpy array if possible
                try:
                    group.create_dataset(key, data=np.array(value))
                except:
                    # Store as string if conversion fails
                    group.attrs[key] = str(value)
    
    with h5py.File(filename, 'w') as f:
        write_to_hdf5(f, results)


def load_results(filename: Union[str, Path]) -> Dict[str, Any]:
    """
    Load calculation results from HDF5 file.
    
    Args:
        filename: Path to results file
        
    Returns:
        Dictionary containing results
    """
    filename = Path(filename)
    
    def read_from_hdf5(group):
        """Recursively read data from HDF5 group."""
        data = {}
        
        # Read attributes
        for key, value in group.attrs.items():
            data[key] = value
        
        # Read datasets and subgroups
        for key in group.keys():
            item = group[key]
            if isinstance(item, h5py.Group):
                data[key] = read_from_hdf5(item)
            elif isinstance(item, h5py.Dataset):
                data[key] = item[()]
        
        return data
    
    with h5py.File(filename, 'r') as f:
        return read_from_hdf5(f)


def format_energy(energy: float, units: str = 'hartree') -> str:
    """
    Format energy value for output.
    
    Args:
        energy: Energy value in Hartree
        units: Output units ('hartree', 'ev', 'kcal/mol')
        
    Returns:
        Formatted energy string
    """
    if units.lower() == 'hartree':
        return f"{energy:.8f} Hartree"
    elif units.lower() == 'ev':
        from .constants import HARTREE_TO_EV
        return f"{energy * HARTREE_TO_EV:.6f} eV"
    elif units.lower() == 'kcal/mol':
        from .constants import HARTREE_TO_KCAL_MOL
        return f"{energy * HARTREE_TO_KCAL_MOL:.3f} kcal/mol"
    else:
        raise ValueError(f"Unknown units: {units}")


def print_geometry(symbols: List[str], coordinates: np.ndarray, 
                  units: str = 'angstrom') -> None:
    """
    Print molecular geometry in a formatted table.
    
    Args:
        symbols: Atomic symbols
        coordinates: Atomic coordinates
        units: Coordinate units
    """
    print(f"\nMolecular Geometry ({units}):")
    print("-" * 50)
    print(f"{'Atom':<6} {'X':>12} {'Y':>12} {'Z':>12}")
    print("-" * 50)
    
    for symbol, coord in zip(symbols, coordinates):
        print(f"{symbol:<6} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}")
    
    print("-" * 50)


def print_orbital_energies(energies: np.ndarray, occupations: np.ndarray = None) -> None:
    """
    Print molecular orbital energies.
    
    Args:
        energies: Orbital energies in Hartree
        occupations: Orbital occupations (optional)
    """
    print("\nMolecular Orbital Energies:")
    print("-" * 40)
    print(f"{'MO':<6} {'Energy (Hartree)':>15} {'Occ':>6}")
    print("-" * 40)
    
    for i, energy in enumerate(energies):
        occ_str = f"{occupations[i]:.1f}" if occupations is not None else "N/A"
        print(f"{i+1:<6} {energy:15.6f} {occ_str:>6}")
    
    print("-" * 40)
