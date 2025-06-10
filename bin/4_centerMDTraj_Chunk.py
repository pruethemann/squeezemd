#!/usr/bin/env python

import argparse
import mdtraj as md

def center_in_chunks(args, chunk_size=100):

    # Load topology separately from .cif file
    top = md.load(args.topo).topology

    # Load the first frame from trajectory for reference and topology saving
    reference = md.load(args.traj, top=top, frame=0)
    alignment_indices = reference.topology.select('backbone')
    

    # Prepare DCD writer
    with md.formats.DCDTrajectoryFile(args.traj_center, 'w', force_overwrite=True) as dcd_out:
        for chunk in md.iterload(args.traj, top=top, chunk=chunk_size):
            chunk.make_molecules_whole()
            chunk.image_molecules(make_whole=False, inplace=True)
            chunk = chunk.superpose(reference, atom_indices=alignment_indices)

            # Write chunk manually
            dcd_out.write(
                xyz=chunk.xyz,
                cell_lengths=chunk.unitcell_lengths,
                cell_angles=chunk.unitcell_angles
            )

    chunk[-1].save(args.topo_center)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--topo', required=True, help='Input CIF (topology from last frame)')
    parser.add_argument('--traj', required=True, help='Input trajectory (.h5)')

    parser.add_argument('--topo_center', required=False, help='Output topology from first frame (.pdb)', default="topo_center.pdb")
    parser.add_argument('--traj_center', required=False, help='Centered output trajectory (.dcd)', default='traj_center.dcd')
    args = parser.parse_args()

    center_in_chunks(args)
    