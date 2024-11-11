import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd  # Import the rmsd function directly

def calculate_rmsd(dcd_file, topology_file, selection='backbone'):
    """
    Calculates the RMSD of a molecular dynamics trajectory.

    Args:
        dcd_file: Path to the DCD trajectory file.
        topology_file: Path to the topology file (e.g., PDB or PSF).
        selection: Atom selection string for RMSD calculation (default: 'backbone').

    Returns:
        rmsd_values: A NumPy array of RMSD values for each frame.
    """
    # Load the trajectory and topology
    u = mda.Universe(topology_file, dcd_file)

    # Select the atoms for RMSD calculation
    atoms = u.select_atoms(selection)

    # Create a reference structure from the first frame
    reference_coords = atoms.positions.copy()

    # Initialize an array to store RMSD values
    rmsd_values = np.zeros(len(u.trajectory))

    # Calculate RMSD for each frame
    for ts in u.trajectory:
        rmsd_values[ts.frame] = rmsd(atoms.positions, reference_coords)

    return rmsd_values

# Example usage:
dcd_file = 'md-01.dcd'
topology_file = 'ionized.pdb'
rmsd_values = calculate_rmsd(dcd_file, topology_file)

# Print the RMSD values
print(rmsd_values)
