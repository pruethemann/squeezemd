#!/usr/bin/env python

import argparse
import mdtraj as md

def center_in_chunks_h5md(topo, traj,topo_center,traj_center_hdf5, chunk_size=20):
    """
    Centers/unwraps + aligns a trajectory in chunks and writes H5MD/HDF5 output.

    Notes:
      - H5MD/HDF5 trajectories in MDTraj typically store coordinates in **nm**
        (so we do NOT convert to Å like in the DCD example).
      - This writes periodic box vectors (lengths/angles) when available.
    """

    # Load reference frame for alignment + anchor detection
    reference = md.load(traj, top=topo, frame=0)
    alignment_indices = reference.topology.select("backbone")
    protein_anchor = reference.topology.guess_anchor_molecules()

    TrajWriter = md.formats.HDF5TrajectoryFile  # older MDTraj
 
    # Stream trajectory in chunks and write to HDF5
    with TrajWriter(traj_center_hdf5, mode="w", force_overwrite=True) as out:
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

            import numpy as np
            # Ensure numpy arrays (some MDTraj versions are picky)
            xyz = np.asarray(chunk.xyz, dtype=np.float32)  # nm
            cell_lengths = None if chunk.unitcell_lengths is None else np.asarray(chunk.unitcell_lengths, dtype=np.float32)
            cell_angles  = None if chunk.unitcell_angles  is None else np.asarray(chunk.unitcell_angles,  dtype=np.float32)

            # --- Write in a version-tolerant way ---
            # 1) Try the "array + kwargs" signature

            out.write(xyz, time=getattr(chunk, "time", None),cell_lengths=cell_lengths, cell_angles=cell_angles)
     

    # Save centered/aligned topology (last chunk is fine; topology is the same)
    # If you want the centered coordinates of the final frame, you could save last_chunk[-1].
    chunk[-1].save(topo_center)

# deprecated
# TODO: Delete or adjust
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

def convert_hd5f_to_dcd(hd5f_path, topo, dcd_path):

    # Import centered h5md
    traj = md.load(hd5f_path, top=topo)
    traj.save(dcd_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument('--topo', required=True, help='Input CIF (topology from last frame)')
    parser.add_argument('--traj', required=True, help='Input trajectory (.hdf5)')

    # Output
    parser.add_argument('--topo_center', required=False, help='Output topology from first frame (.pdb)', default="topo_center.pdb")
    parser.add_argument('--traj_center', required=False, help='Centered output trajectory (.dcd)', default='traj_center.dcd')
    parser.add_argument('--traj_center_hdf5', required=False, help='Centered output trajectory (.h5md)', default='traj_center.h5')
    args = parser.parse_args()

    # Center protein in middle of water box and remove translation and rotation
    center_in_chunks_h5md(topo=args.topo, traj=args.traj, topo_center=args.topo_center,traj_center_hdf5=args.traj_center_hdf5)

    # Transform h5md to dcd for pymol visulationsion
    convert_hd5f_to_dcd(args.traj_center_hdf5, args.topo_center, args.traj_center)
    