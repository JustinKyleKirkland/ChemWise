"""
Command-line interface for ChemWise quantum chemistry software.
"""

import click
import sys
from pathlib import Path
import yaml

from .core.calculator import Calculator, calculate
from .core.molecule import Molecule
from .utils.io import read_input_file, write_input_file, load_molecule


@click.group()
@click.version_option(version="0.1.0", prog_name="ChemWise")
def main():
    """ChemWise - Quantum Chemistry Software Suite"""
    pass


@main.command()
@click.option('--input', '-i', 'input_file', type=click.Path(exists=True), 
              help='Input file (XYZ or YAML format)')
@click.option('--method', '-m', default='hf', 
              help='Calculation method (hf, dft, b3lyp, etc.)')
@click.option('--basis', '-b', default='sto-3g',
              help='Basis set (sto-3g, cc-pvdz, etc.)')
@click.option('--charge', '-c', default=0, type=int,
              help='Molecular charge')
@click.option('--multiplicity', '-s', default=1, type=int,
              help='Spin multiplicity')
@click.option('--convergence', default=1e-8, type=float,
              help='SCF convergence threshold')
@click.option('--max-iterations', default=100, type=int,
              help='Maximum SCF iterations')
@click.option('--output', '-o', 'output_file', type=click.Path(),
              help='Output file for results')
@click.option('--save-wavefunction', type=click.Path(),
              help='Save wavefunction to HDF5 file')
@click.option('--geometry', type=str,
              help='Molecular geometry as string (alternative to input file)')
def calculate_cmd(input_file, method, basis, charge, multiplicity, convergence,
                  max_iterations, output_file, save_wavefunction, geometry):
    """Run a quantum chemistry calculation."""
    
    try:
        # Load molecule
        if input_file:
            if Path(input_file).suffix.lower() in ['.yml', '.yaml']:
                input_data = read_input_file(input_file)
                molecule = Molecule.from_dict(input_data.get('molecule', {}))
                
                # Override with command line options
                if 'calculation' in input_data:
                    calc_params = input_data['calculation']
                    method = calc_params.get('method', method)
                    basis = calc_params.get('basis', basis)
                    if 'scf' in calc_params:
                        scf_params = calc_params['scf']
                        convergence = scf_params.get('convergence', convergence)
                        max_iterations = scf_params.get('max_iterations', max_iterations)
            else:
                molecule_data = load_molecule(input_file)
                molecule = Molecule.from_dict(molecule_data)
                molecule.charge = charge
                molecule.multiplicity = multiplicity
        
        elif geometry:
            molecule = Molecule.from_string(geometry, charge, multiplicity)
        
        else:
            click.echo("Error: Must provide either --input file or --geometry string")
            sys.exit(1)
        
        # Create calculator
        calc = Calculator(
            molecule=molecule,
            method=method,
            basis=basis,
            convergence=convergence,
            max_iterations=max_iterations
        )
        
        click.echo(f"Running {method.upper()} calculation with {basis} basis set...")
        
        # Run calculation
        result = calc.run(print_results=True)
        
        # Save results if requested
        if output_file:
            click.echo(f"\nSaving results to {output_file}")
            with open(output_file, 'w') as f:
                f.write(result.summary())
        
        if save_wavefunction:
            click.echo(f"Saving wavefunction to {save_wavefunction}")
            calc.save_results(save_wavefunction, result)
        
        if result.converged:
            click.echo(f"\nCalculation completed successfully!")
            click.echo(f"Final energy: {result.energy:.8f} Hartree")
        else:
            click.echo("\nWarning: Calculation did not converge!")
            sys.exit(1)
    
    except Exception as e:
        click.echo(f"Error: {e}")
        sys.exit(1)


@main.command()
@click.option('--molecule', '-m', type=str, required=True,
              help='Molecule geometry (XYZ file or geometry string)')
@click.option('--method', default='hf',
              help='Calculation method')
@click.option('--basis', default='sto-3g',
              help='Basis set')
@click.option('--charge', '-c', default=0, type=int,
              help='Molecular charge')
@click.option('--multiplicity', '-s', default=1, type=int,
              help='Spin multiplicity')
@click.option('--output', '-o', 'output_file', type=click.Path(), required=True,
              help='Output YAML input file')
def generate_input(molecule, method, basis, charge, multiplicity, output_file):
    """Generate a ChemWise input file."""
    
    try:
        # Load or parse molecule
        if Path(molecule).exists():
            molecule_data = load_molecule(molecule)
        else:
            # Treat as geometry string
            mol = Molecule.from_string(molecule, charge, multiplicity)
            molecule_data = mol.to_dict()
        
        # Create input structure
        input_data = {
            'molecule': molecule_data,
            'calculation': {
                'method': method,
                'basis': basis,
                'scf': {
                    'max_iterations': 100,
                    'convergence': 1e-8,
                    'diis': True
                }
            },
            'output': {
                'orbitals': True,
                'population_analysis': True
            }
        }
        
        # Write input file
        write_input_file(output_file, input_data)
        click.echo(f"Generated input file: {output_file}")
    
    except Exception as e:
        click.echo(f"Error: {e}")
        sys.exit(1)


@main.command()
@click.option('--method', default='hf',
              help='Optimization method')
@click.option('--basis', default='sto-3g',
              help='Basis set')
@click.option('--input', '-i', 'input_file', type=click.Path(exists=True), required=True,
              help='Input molecule file')
@click.option('--output', '-o', 'output_file', type=click.Path(),
              help='Output optimized geometry file')
def optimize(method, basis, input_file, output_file):
    """Optimize molecular geometry (placeholder)."""
    
    click.echo("Geometry optimization is not yet implemented.")
    click.echo("This feature will be available in a future version.")


@main.command()
@click.option('--basis-sets', is_flag=True, help='List available basis sets')
def info(basis_sets):
    """Display information about ChemWise."""
    
    if basis_sets:
        from .core.basis import list_available_basis_sets
        click.echo("Available Basis Sets:")
        click.echo("=====================")
        for basis in list_available_basis_sets():
            click.echo(f"  - {basis.upper()}")
        return
    
    info_text = """
ChemWise - Quantum Chemistry Software Suite v0.1.0

Supported Methods:
  - Hartree-Fock (HF)
  - Density Functional Theory (DFT)
    * LDA (Local Density Approximation)
    * B88 exchange functional
    * LYP correlation functional
    * B3LYP hybrid functional

Supported Basis Sets:
  - STO-3G, STO-6G (minimal basis sets)
  - 6-31G (split valence basis set)
  - cc-pVDZ (correlation consistent double zeta)
  - def2-SVP (Ahlrichs split valence)
  - More available via JSON files in basis/ directory

Features:
  - Self-Consistent Field (SCF) calculations
  - DIIS convergence acceleration
  - Mulliken population analysis
  - Molecular orbital analysis
  - HDF5 wavefunction output
  - Extensible basis set framework

Recent Updates:
  - Added support for loading basis sets from JSON files
  - Improved integral calculation stability
  - Enhanced basis set library

Use 'chemwise info --basis-sets' to list all available basis sets.
    """
    
    click.echo(info_text)


@main.command()
@click.argument('molecule_file', type=click.Path(exists=True))
def view(molecule_file):
    """Display molecular geometry."""
    
    try:
        molecule_data = load_molecule(molecule_file)
        molecule = Molecule.from_dict(molecule_data)
        
        from .utils.io import print_geometry
        print_geometry(molecule.symbols, molecule.coordinates_angstrom)
        
        click.echo(f"Charge: {molecule.charge}")
        click.echo(f"Multiplicity: {molecule.multiplicity}")
        click.echo(f"Number of atoms: {molecule.n_atoms}")
        click.echo(f"Number of electrons: {molecule.n_electrons}")
        click.echo(f"Nuclear repulsion energy: {molecule.nuclear_repulsion_energy:.8f} Hartree")
    
    except Exception as e:
        click.echo(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
