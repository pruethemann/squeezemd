#!/usr/bin/env python

import argparse
import mdtraj


def centerMDTraj(args):

    # Import Trajectory
    traj = mdtraj.load(args.traj, top=args.topo)

    # TODO: this single command works in python 3.10 but not 3.11. Activate as soon bug is fixed in MDTraj
    #traj.image_molecules(inplace=False)

    # Center Trajectory
    traj.make_molecules_whole()
    traj.image_molecules(make_whole=False)

    # Save centered trajectory as dcd
    traj.save_dcd(args.traj_center)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--topo', required=False, help='Cif file of last frame')
    parser.add_argument('--traj', required=False, help='Trajectory file')

    # Output

    parser.add_argument('--traj_center', required=False, help='MD parameter file saved for every MD')

    args = parser.parse_args()

    # Parse Arguments
    centerMDTraj(args)
