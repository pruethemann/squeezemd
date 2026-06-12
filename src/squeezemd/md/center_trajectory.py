#!/usr/bin/env python

import argparse

import mdtraj as md
import numpy as np


def center_in_chunks_h5(topo, traj, topo_center, traj_center_h5, chunk_size=20):
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
    with TrajWriter(traj_center_h5, mode="w", force_overwrite=True) as out:
        for chunk in md.iterload(traj, top=topo, chunk=chunk_size):
            # 1) ensure molecules are whole first (unwrap)
            chunk.make_molecules_whole(inplace=True)

            # 2) image with molecules kept whole and anchored to the protein
            chunk.image_molecules(
                make_whole=True,
                anchor_molecules=protein_anchor,
                inplace=True,
            )

            # 3) superpose to reference using backbone
            chunk = chunk.superpose(reference, frame=0, atom_indices=alignment_indices)

            # Ensure numpy arrays (some MDTraj versions are picky)
            xyz = np.asarray(chunk.xyz, dtype=np.float32)  # nm
            cell_lengths = (
                None if chunk.unitcell_lengths is None else np.asarray(chunk.unitcell_lengths, dtype=np.float32)
            )
            cell_angles = None if chunk.unitcell_angles is None else np.asarray(chunk.unitcell_angles, dtype=np.float32)

            # --- Write in a version-tolerant way ---
            out.write(xyz, time=getattr(chunk, "time", None), cell_lengths=cell_lengths, cell_angles=cell_angles)

    # Save centered/aligned topology (last chunk is fine; topology is the same)
    chunk[-1].save(topo_center)


def convert_h5_to_dcd(h5_path, topo, dcd_path, stride=10):
    """
    Convert HDF5 trajectory to DCD, saving only every `stride`-th frame.

    Parameters
    ----------
    h5_path : str
        Path to input HDF5 trajectory
    topo : str
        Path to topology file (PDB/PSF/etc.)
    dcd_path : str
        Path to output DCD file
    stride : int, optional
        Save every `stride`-th frame (default: 10)
    """
    traj = md.load(h5_path, top=topo)
    traj_strided = traj[::stride]
    traj_strided.save(dcd_path)


def main():
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument("--topo", required=True, help="Input CIF (topology from last frame)")
    parser.add_argument("--traj", required=True, help="Input trajectory (.hdf5)")

    # Output
    parser.add_argument(
        "--topo_center", required=False, help="Output topology from first frame (.pdb)", default="topo_center.pdb"
    )
    parser.add_argument(
        "--traj_center", required=False, help="Centered output trajectory (.dcd)", default="traj_center.dcd"
    )
    parser.add_argument(
        "--traj_center_h5", required=False, help="Centered output trajectory (.h5)", default="traj_center.h5"
    )
    args = parser.parse_args()

    # Center protein in middle of water box and remove translation and rotation
    center_in_chunks_h5(
        topo=args.topo, traj=args.traj, topo_center=args.topo_center, traj_center_h5=args.traj_center_h5
    )

    # Transform h5 to DCD for PyMOL visualization
    convert_h5_to_dcd(args.traj_center_h5, args.topo_center, args.traj_center)


if __name__ == "__main__":
    main()
