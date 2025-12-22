#!/usr/bin/env python

import argparse
import mdtraj as md

def center_in_chunks_h5md(args, chunk_size=20):
    """
    Centers/unwraps + aligns a trajectory in chunks and writes H5MD/HDF5 output.

    Notes:
      - H5MD/HDF5 trajectories in MDTraj typically store coordinates in **nm**
        (so we do NOT convert to Å like in the DCD example).
      - This writes periodic box vectors (lengths/angles) when available.
    """

    # Load reference frame for alignment + anchor detection
    reference = md.load(args.traj, top=args.topo, frame=0)
    alignment_indices = reference.topology.select("backbone")
    protein_anchor = reference.topology.guess_anchor_molecules()

    # Prefer H5MDTrajectoryFile if available; otherwise fall back to HDF5TrajectoryFile
    try:       
        TrajWriter = md.formats.H5MDTrajectoryFile
        print("Excellent: Use H5MDTrajectoryFile")
    except AttributeError:
        TrajWriter = md.formats.HDF5TrajectoryFile  # older MDTraj
        print("All right: Use older HDF5TrajectoryFile")

    # Stream trajectory in chunks and write to H5MD/HDF5
    with TrajWriter(args.traj_center, mode="w", force_overwrite=True) as out:
        for chunk in md.iterload(args.traj, top=args.topo, chunk=chunk_size):

            # 1) ensure molecules are whole first
            chunk.make_molecules_whole(inplace=True)

            # 2) image with molecules kept whole and anchored to the protein
            chunk.image_molecules(
                make_whole=True,
                anchor_molecules=protein_anchor,
                inplace=True,
            )

            # 3) superpose to reference using backbone
            chunk = chunk.superpose(reference, frame=0, atom_indices=alignment_indices)

            # Write (keep MDTraj's native units: xyz in nm; unitcell_lengths in nm)
            # Some formats accept None if unit cell isn't present.
            out.write(
                xyz=chunk.xyz,
                cell_lengths=chunk.unitcell_lengths,
                cell_angles=chunk.unitcell_angles,
            )

    # Save centered/aligned topology (last chunk is fine; topology is the same)
    # If you want the centered coordinates of the final frame, you could save last_chunk[-1].
    chunk[-1].save(args.topo_center)

def center_in_chunks_dcd(args, chunk_size=20):
    # TODO: Calculate chunk_size according to number of trajectory frames

    # Load the first frame from trajectory for reference and topology saving
    reference = md.load(args.traj, top=args.topo, frame=0)
    alignment_indices = reference.topology.select('backbone')
    protein_anchor = reference.topology.guess_anchor_molecules()

    # Prepare DCD writer and center in chunks for better memory efficiency
    with md.formats.DCDTrajectoryFile(args.traj_center, 'w', force_overwrite=True) as dcd_out:
        for chunk in md.iterload(args.traj, top=args.topo, chunk=chunk_size):
            
            # 1) ensure molecules are whole first
            chunk.make_molecules_whole(inplace=True)

            # 2) image with molecules kept whole and anchored to the protein
            chunk.image_molecules(make_whole=True, 
                                  anchor_molecules=protein_anchor,
                                  inplace=True)
           
            # 3) superpose to reference using backbone
            chunk = chunk.superpose(reference, frame=0, atom_indices=alignment_indices)

            # Convert nm → Å for output
            xyz_angstrom = chunk.xyz * 10.0
            cell_lengths = chunk.unitcell_lengths * 10.0

            # Write chunk manually
            dcd_out.write(
                xyz=xyz_angstrom,
                cell_lengths=cell_lengths,
                cell_angles=chunk.unitcell_angles
            )

    chunk[-1].save(args.topo_center)

def transform_h5md_to_dcd(h5md_path, topo_path, dcd_path, selection="all", stride=1):
    """
    Convert an H5MD trajectory to DCD.

    Parameters
    ----------
    h5md_path : str
        Path to input .h5 / .h5md file.
    topo_path : str
        Path to topology (e.g., .pdb, .psf, .prmtop).
    dcd_path : str
        Path to output .dcd file.
    selection : str
        Atom selection in MDAnalysis syntax (default: "all").
    stride : int
        Write every `stride`-th frame (default: 1).
    """
    import MDAnalysis as mda

    u = mda.Universe(topo_path, h5md_path)
    ag = u.select_atoms(selection)

    with mda.Writer(dcd_path, ag.n_atoms) as w:
        for ts in u.trajectory[::stride]:
            w.write(ag)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument('--topo', required=True, help='Input CIF (topology from last frame)')
    parser.add_argument('--traj', required=True, help='Input trajectory (.h5)')

    # Output
    parser.add_argument('--topo_center', required=False, help='Output topology from first frame (.pdb)', default="topo_center.pdb")
    parser.add_argument('--traj_center', required=False, help='Centered output trajectory (.dcd)', default='traj_center.dcd')
    parser.add_argument('--traj_center_h5', required=False, help='Centered output trajectory (.h5md)', default='traj_center.h5')
    args = parser.parse_args()

    # Center protein in middle of water box and remove translation and rotation
    center_in_chunks_h5md(args)

    # Transform h5md to dcd for pymol visulationsion
    transform_h5md_to_dcd(args.traj_center_h5, args.topo_center, args.traj_center)
    