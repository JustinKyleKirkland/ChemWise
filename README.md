# ChemWise - Quantum Chemistry Software Suite

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

ChemWise is a comprehensive quantum chemistry software suite designed for performing high-quality electronic structure calculations. Currently supports Hartree-Fock (HF) and Density Functional Theory (DFT) methods.

## Features

### Current Implementation
- **Hartree-Fock Theory**
  - Restricted Hartree-Fock (RHF) for closed-shell systems
  - Unrestricted Hartree-Fock (UHF) for open-shell systems
  - Self-Consistent Field (SCF) convergence
  - Various acceleration techniques (DIIS, level shifting)

- **Density Functional Theory**
  - Local Density Approximation (LDA)
  - Generalized Gradient Approximation (GGA)
  - Hybrid functionals (B3LYP, PBE0)
  - Range-separated functionals

- **Basis Sets**
  - Support for Gaussian basis sets
  - STO-nG, cc-pVXZ, aug-cc-pVXZ basis sets
  - Custom basis set input

- **Molecular Properties**
  - Energy calculations
  - Molecular orbitals and orbital energies
  - Mulliken population analysis
  - Dipole moments
  - Vibrational frequencies (planned)

### Planned Features
- Post-Hartree-Fock methods (MP2, CCSD, etc.)
- Excited state calculations (CIS, TD-DFT)
- Solvent effects (PCM)
- Geometry optimization
- Molecular dynamics

## Installation

### From Source
```bash
git clone https://github.com/yourusername/chemwise.git
cd chemwise
pip install -e .
```

### Development Installation
```bash
git clone https://github.com/yourusername/chemwise.git
cd chemwise
pip install -e ".[dev]"
```

## Quick Start

### Command Line Interface
```bash
# Run a simple water molecule calculation
chemwise calculate --input water.yml --method hf --basis sto-3g

# Run DFT calculation with B3LYP
chemwise calculate --input water.yml --method b3lyp --basis cc-pvdz
```

### Python API
```python
import chemwise as cw

# Define molecule
mol = cw.Molecule.from_xyz("water.xyz")

# Set up calculation
calc = cw.Calculator(
    molecule=mol,
    method="b3lyp",
    basis="cc-pvdz"
)

# Run calculation
result = calc.run()

print(f"Total Energy: {result.energy:.8f} Hartree")
print(f"Dipole Moment: {result.dipole} Debye")
```

### Input File Format

ChemWise uses YAML input files for easy configuration:

```yaml
molecule:
  geometry: |
    O  0.0000000  0.0000000  0.1173000
    H  0.0000000  0.7572000 -0.4692000
    H  0.0000000 -0.7572000 -0.4692000
  charge: 0
  multiplicity: 1

calculation:
  method: b3lyp
  basis: cc-pvdz
  scf:
    max_iterations: 100
    convergence: 1e-8
    diis: true

output:
  orbitals: true
  population_analysis: true
  save_wavefunction: result.h5
```

## Project Structure

```
chemwise/
├── core/              # Core quantum chemistry algorithms
│   ├── molecule.py    # Molecular structure and properties
│   ├── basis.py       # Basis set handling
│   ├── integrals.py   # Integral evaluation
│   └── scf.py         # Self-consistent field procedures
├── methods/           # Quantum chemistry methods
│   ├── hf.py          # Hartree-Fock implementation
│   ├── dft.py         # Density functional theory
│   └── functionals.py # Exchange-correlation functionals
├── utils/             # Utility functions
│   ├── constants.py   # Physical constants
│   ├── math_utils.py  # Mathematical utilities
│   └── io.py          # Input/output handling
├── cli.py             # Command-line interface
└── __init__.py        # Package initialization
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use ChemWise in your research, please cite:

```
ChemWise: A Python Quantum Chemistry Software Suite
Author: Justin Kirkland
Version: 0.1.0
```

## Acknowledgments

- The quantum chemistry community for theoretical foundations
- NumPy and SciPy developers for numerical computing tools
- The Python scientific computing ecosystem
