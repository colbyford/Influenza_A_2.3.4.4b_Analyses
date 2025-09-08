import torch
import argparse

from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units
from ase.md import MDLogger

from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator


def setup_device():
    """Set up and return the appropriate compute device."""
    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
    print(f"Using device: {device}")
    return device


def run_md_simulation(
    input_file: str = "NaClWater.xyz",
    cell_size: float = 25.25,
    temperature_K: float = 300,
    timestep: float = 0.5 * units.fs,
    friction: float = 0.01 / units.fs,
    total_steps: int = 100,
    traj_interval: int = 20,
    log_interval: int = 1,
    output_file: str = "NaClWaterMD_out.xyz"
):
    """Run molecular dynamics simulation with specified parameters.

    Args:
        input_file: Path to input XYZ file
        cell_size: Size of cubic simulation cell
        temperature_K: Temperature in Kelvin
        timestep: MD timestep
        friction: Langevin friction coefficient
        total_steps: Total number of MD steps
        traj_interval: Interval for trajectory writing
        log_interval: Interval for log writing
    """
    # Set up device
    device = setup_device()

    # Read in the system from file and set the cell size and pbc
    atoms = read(input_file)
    atoms.set_cell([cell_size] * 3)
    atoms.set_pbc([True] * 3)

    # Set the calculator
    atoms.calc = ORBCalculator(pretrained.orb_v3_direct_20_omat(), device=device)#, compile=False)

    # Set the initial velocities
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K)

    # Set the dynamics
    dyn = Langevin(atoms, timestep, temperature_K=temperature_K, friction=friction)

    # Define output functions and attach to dynamics
    dyn.attach(
        lambda: write(output_file, atoms, append=True), interval=traj_interval
    )
    dyn.attach(MDLogger(dyn, atoms, "md_nvt.log"), interval=log_interval)

    # Run the dynamics
    dyn.run(steps=total_steps)


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Run MD simulations with Orb.")
    parser.add_argument("--input_file", type=str,
                        help="Path to XYZ file")
    parser.add_argument("--cell_size", type=float, default=25.25, nargs='?',
                        help="Size of cubic simulation cell (default: 25.25)")
    parser.add_argument("--temperature_K", type=float, default=300, nargs='?',
                        help="Temperature in Kelvin (default: 300)")
    parser.add_argument("--timestep", type=float, default=0.5, nargs='?',
                        help="MD timestep in fs (default: 0.5 fs)")
    parser.add_argument("--friction", type=float, default=0.01, nargs='?',
                        help="Langevin friction coefficient (default: 0.01 / fs)")
    parser.add_argument("--total_steps", type=int, default=100, nargs='?',
                        help="Total number of MD steps (default: 100)")
    parser.add_argument("--traj_interval", type=int, default=20, nargs='?',
                        help="Trajectory writing interval (default: 20 steps)")
    parser.add_argument("--log_interval", type=int, default=1, nargs='?',
                        help="Log writing interval (default: 1 step)")
    parser.add_argument("--output_file", type=str,
                        help="Path to output XYZ file.")

    args = parser.parse_args()

    input_file = args.input_file
    cell_size = args.cell_size
    temperature_K = args.temperature_K
    timestep = args.timestep
    friction = args.friction
    total_steps = args.total_steps
    traj_interval = args.traj_interval
    log_interval = args.log_interval
    output_file = args.output_file

    run_md_simulation(input_file=input_file, cell_size=cell_size, temperature_K=temperature_K,
                      timestep=timestep, friction=friction, total_steps=total_steps,
                      traj_interval=traj_interval, log_interval=log_interval,
                      output_file=output_file)


if __name__ == "__main__":
    main()
