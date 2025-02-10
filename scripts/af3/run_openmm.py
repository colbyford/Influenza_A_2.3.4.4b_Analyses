import numpy as np
import matplotlib.pyplot as plt
from sys import stdout
import logging
import argparse

from openff.toolkit.topology import Molecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm import *
from openmm.app import *
from openmm.unit import *

from pdbfixer import PDBFixer
from openmm.app import PDBFile

parser = argparse.ArgumentParser()

parser.add_argument("--input_pdb_path", help="The file path to the input PDB protein file.")
parser.add_argument("--input_sdf_path", help="The file path to the input SDF ligand file.")
parser.add_argument("--output_pdb_path", help="The file path to the output minimized PDB complex file.")

args = parser.parse_args()

input_pdb_path = args.input_pdb_path
input_sdf_path = args.input_sdf_path
output_pdb_path = args.output_pdb_path


def fix_pdb(input_pdb_path, output_pdb_path):
    ## Load the protein PDB file and fix missing atoms
    fixer = PDBFixer(input_pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingHydrogens()
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_path, "w"))
    

def run_simulation(input_pdb_path, input_sdf_path, output_pdb_path, platform_name='CUDA'):
    
    ## Set up platform
    logging.info(f"Setting up {platform_name} platform")
    if platform_name == 'CUDA':
        platform = Platform.getPlatformByName('CUDA')
        platform.setPropertyDefaultValue('Precision', 'mixed')
    elif platform_name == 'CPU':
        platform = Platform.getPlatformByName('CPU')
    
    ##########
    ## PROTEIN SECTION
    ##########

    ## Read in PDB
    logging.info(f"Reading in PDB")
    pdb = PDBFile(input_pdb_path)

    ## Create a Modeller object
    logging.info(f"Creating Modeller object")
    modeller = Modeller(pdb.topology, pdb.positions)

    ##########
    ## LIGAND SECTION
    ##########

    ## Create an OpenFF Molecule object from SDF file
    logging.info(f"Loading SDF file and generating template")
    molecule = Molecule.from_file(input_sdf_path)
    molecule.assign_partial_charges(partial_charge_method="mmff94")

    ## Create the SMIRNOFF template generator
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)

    ## Create an OpenMM ForceField object
    logging.info(f"Creating ForceField object")
    forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml')

    ## Register the SMIRNOFF template generator
    forcefield.registerTemplateGenerator(smirnoff.generator)

    ## Make an OpenFF Topology of the ligand
    molecule_off_topology = offTopology.from_molecules(molecules=[molecule])

    ## Convert it to an OpenMM Topology
    molecule_topology = molecule_off_topology.to_openmm()

    ## Get the positions of the ligand
    molecule_positions = offquantity_to_openmm(molecule.conformers[0])

    ## Add the ligand to the Modeller
    logging.info(f"Combining protein and ligand topologies")
    modeller.add(molecule_topology, molecule_positions)

    ##########
    ## SIMULATION SECTION
    ##########

    ## Solvate
    logging.info(f"Adding solvent")
    # modeller.addSolvent(forcefield, padding=1.0*nanometer, ionicStrength=0.15*molar)
    modeller.addSolvent(forcefield, padding=1.0*nanometer)

    ## Create system
    logging.info(f"Creating system")
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

    ## Setup Simulation
    logging.info(f"Setting up simulation")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    ## Minimize Energy
    logging.info("Minimizing Energy")
    simulation.minimizeEnergy()

    simulation.reporters.append(PDBReporter(output_pdb_path, 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True, volume=True))

    output_log_path = f"{output_pdb_path.split('.')[0]}.log"

    simulation.reporters.append(StateDataReporter(output_log_path, 100, step=True,
            potentialEnergy=True, temperature=True, volume=True))

    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    simulation.context.reinitialize(preserveState=True)

    logging.info("Running NPT")
    simulation.step(10000)

    logging.info("Running NVT")
    simulation.step(10000)

    logging.info("Simulation complete")


##########
## PLOTTING SECTION
##########
def plot_results(output_log_path):

    logging.info("Plotting results")
    ## Read Log Data
    data = np.loadtxt(output_log_path, delimiter=',')

    step = data[:,0]
    potential_energy = data[:,1]
    temperature = data[:,2]
    volume = data[:,3]

    ## Potential Energy
    plt.plot(step, potential_energy)
    plt.xlabel("Step")
    plt.ylabel("Potential energy (kJ/mol)")
    plt.savefig(f"{output_pdb_path.split('.')[0]}_potential_energy.png")

    ## Temperature
    plt.plot(step, temperature)
    plt.xlabel("Step")
    plt.ylabel("Temperature (K)")
    plt.savefig(f"{output_pdb_path.split('.')[0]}_temperature.png")

    ## Volume
    plt.plot(step, volume)
    plt.xlabel("Step")
    plt.ylabel("Volume (nm^3)")
    plt.savefig(f"{output_pdb_path.split('.')[0]}_volume.png")


if __name__ == '__main__':
    ## Set up logging
    logging.basicConfig(level=logging.INFO)
    logging.info("Starting OpenMM Molecular Dynamics...")

    ## Fix the PDB file
    output_fixed_pdb_path = f"{input_pdb_path.split('.')[0]}_fixed.pdb"
    fix_pdb(input_pdb_path, output_fixed_pdb_path)

    ## Run the simulation
    run_simulation(output_fixed_pdb_path,
                   input_sdf_path,
                   output_pdb_path,
                   platform_name='CUDA')

    ## Plot the results
    plot_results(f"{output_pdb_path.split('.')[0]}.log")