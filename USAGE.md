# ChemWise Usage Guide

## Quick Start

### Installation
```bash
# Clone the repository
git clone https://github.com/yourusername/chemwise.git
cd chemwise

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Running Calculations

#### 1. Hartree-Fock Calculation
```bash
# Simple HF calculation on H2
chemwise calculate -i examples/h2.xyz -m hf -b sto-3g

# With custom parameters
chemwise calculate -i examples/water.xyz -m hf -b sto-3g --convergence 1e-10 --max-iterations 50
```

#### 2. DFT Calculations
```bash
# B3LYP calculation
chemwise calculate -i examples/water.xyz -m b3lyp -b sto-3g

# LDA calculation
chemwise calculate -i examples/methane.xyz -m lda -b sto-3g
```

#### 3. Using YAML Input Files
```bash
# Run calculation from YAML input
chemwise calculate -i examples/water_b3lyp.yml
```

### Viewing Molecular Geometries
```bash
# View geometry and basic properties
chemwise view examples/water.xyz
```

### Generating Input Files
```bash
# Generate YAML input file
chemwise generate-input --molecule examples/h2.xyz --method b3lyp --basis sto-3g -o h2_b3lyp.yml
```

## Example Calculations

### Hydrogen Molecule (H2)
This is the simplest two-electron system, perfect for testing:

```bash
# HF calculation
chemwise calculate -i examples/h2.xyz -m hf -b sto-3g

# Expected energy: ~-1.13 Hartree
```

### Water Molecule (H2O)
Classic benchmark system:

```bash
# HF calculation
chemwise calculate -i examples/water.xyz -m hf -b sto-3g

# B3LYP calculation
chemwise calculate -i examples/water.xyz -m b3lyp -b sto-3g
```

### Using Python API

```python
import chemwise as cw

# Load molecule
mol = cw.Molecule.from_xyz("examples/water.xyz")

# Create calculator
calc = cw.Calculator(mol, method="b3lyp", basis="sto-3g")

# Run calculation
result = calc.run()

print(f"Total Energy: {result.energy:.8f} Hartree")
```

## Available Methods

- **hf**: Hartree-Fock
- **b3lyp**: B3LYP hybrid functional
- **lda**: Local Density Approximation
- **dft**: Alias for B3LYP

## Available Basis Sets

- **sto-3g**: STO-3G minimal basis
- **sto-6g**: STO-6G extended minimal basis

## Output Files

### Text Output
Use `--output` to save results to a text file:
```bash
chemwise calculate -i water.xyz -m hf -b sto-3g --output water_hf.out
```

### HDF5 Wavefunction Files
Use `--save-wavefunction` to save all data in HDF5 format:
```bash
chemwise calculate -i water.xyz -m b3lyp -b sto-3g --save-wavefunction water_b3lyp.h5
```

## Performance Tips

1. **Use smaller molecules** for testing (H2, H2O)
2. **Start with STO-3G basis** for quick calculations
3. **Hybrid functionals are slower** due to exact exchange
4. **Grid quality affects DFT speed**: coarse < medium < fine < ultrafine

## Troubleshooting

### Common Issues

1. **SCF not converging**: Try lower convergence threshold or more iterations
2. **Slow calculations**: Use smaller basis sets or coarse DFT grids
3. **Memory issues**: Reduce system size or basis set

### Debug Mode
```bash
# Add verbose output
chemwise calculate -i molecule.xyz -m hf -b sto-3g --max-iterations 5
```

## Future Features

The following features are planned for future releases:

- More basis sets (cc-pVDZ, cc-pVTZ, def2-SVP, etc.)
- Post-HF methods (MP2, CCSD)
- Geometry optimization
- Vibrational analysis
- Excited state calculations (CIS, TD-DFT)
- Solvent effects (PCM)

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file for details.
