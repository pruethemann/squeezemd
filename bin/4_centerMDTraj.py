#!/usr/bin/env python

import argparse
import mdtraj


def centerMDTraj(args):

    # Import Trajectory
    traj = mdtraj.load(args.traj, top=args.topo)

    # TODO: this single command works in python 3.10 but not 3.11. Activate as soon bug is fixed in MDTraj
    #traj.image_molecules(inplace=False)

    # Center Trajectory. Image the molecules (deal with periodic boundary conditions)
    traj.make_molecules_whole()
    traj.image_molecules(make_whole=False, inplace=True)

    # Select atoms for alignment (e.g., protein backbone)
    # NOTE: The other way around. Align and then image the traj. doesn't work. There is a bug
    # which leads to water holes
    alignment_indices = traj.topology.select('backbone')

    # Superpose the trajectory to the first frame (or another reference frame)
    traj = traj.superpose(traj, frame=0, atom_indices=alignment_indices)

    # Save centered trajectory as dcd
    traj.save_dcd(args.traj_center)

    traj[-1].save(args.topo_center)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--topo', required=False, help='Cif file of last frame')
    parser.add_argument('--traj', required=False, help='Trajectory file')

    # Output
    parser.add_argument('--topo_center', required=False, help='MD parameter file saved for every MD')
    parser.add_argument('--traj_center', required=False, help='MD parameter file saved for every MD')

    args = parser.parse_args()

    # Parse Arguments
    centerMDTraj(args)
