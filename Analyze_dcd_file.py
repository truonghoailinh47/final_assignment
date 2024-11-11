import MDAnalysis as mda

def save_dcd_data_to_text(dcd_file, topology_file, output_file, selection='all'):
    """
    Opens a DCD trajectory file and saves all atomic position data for each frame to a text file.

    Args:
        dcd_file: Path to the DCD trajectory file.
        topology_file: Path to the topology file (e.g., PDB or PSF).
        output_file: Path to the output text file.
        selection: Atom selection string for saving positions (default: 'all').

    Returns:
        None
    """
    # Load the trajectory and topology
    u = mda.Universe(topology_file, dcd_file)

    # Check the total number of frames
    num_frames = len(u.trajectory)
    print(f"Total frames in the DCD file: {num_frames}")

    # Select the atoms for position extraction
    atoms = u.select_atoms(selection)

    with open(output_file, 'w') as f:
        # Write header for clarity
        f.write("Frame\tAtom Index\tX\tY\tZ\n")

        # Loop over each frame in the trajectory
        for ts in u.trajectory:
            frame_num = ts.frame
            f.write(f"Frame {frame_num}\n")

            # Extract the positions of the selected atoms
            positions = atoms.positions
            for i, pos in enumerate(positions):
                f.write(f"{frame_num}\t{i}\t{pos[0]:.3f}\t{pos[1]:.3f}\t{pos[2]:.3f}\n")
            
            # Optional: Add a blank line between frames for clarity
            f.write("\n")

    print(f"DCD data saved to {output_file}")

# Example usage:
dcd_file = 'md-01.dcd'
topology_file = 'ionized.pdb'
output_file = 'all_dcd_data.txt'
save_dcd_data_to_text(dcd_file, topology_file, output_file)
