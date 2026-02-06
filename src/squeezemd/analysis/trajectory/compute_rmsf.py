#!/usr/bin/env python

"""Compute RMSF across one or more trajectories and save as parquet."""

import argparse
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import openmm.app as app
from ...helper_functions import remap_MDAnalysis

def calculate_RMSF(u: mda.Universe, i):
    """Calculate Cα RMSF for a single trajectory and label by simulation id."""

    # TODO: separate ligand and receptor. currently all Cα atoms
    c_alphas = u.select_atoms(f'name CA')
    R = rms.RMSF(c_alphas).run()

    # Store RMSF and secondary structure data
    rmsf_df = {'resid':c_alphas.resids, 
               'rmsf': R.results.rmsf, 
               'sim_id': i}

    return pd.DataFrame(rmsf_df)


def parse_arguments():
    """Parse CLI arguments for RMSF computation."""
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', nargs='+', required=False)
    parser.add_argument('--traj', nargs='+', required=False)

    # Output
    parser.add_argument('--output', required=False, default='rmsf.svg', help='')

    return parser.parse_args()

def main():
    args = parse_arguments()

    topos = sorted(args.topo)
    trajs = sorted(args.traj)

    rmsf_data = []

    for i,(topo,traj) in enumerate(zip(topos,trajs)):

        # Import Trajectory
        topo = app.PDBxFile(topo)
        u = mda.Universe(topo, traj, in_memory=False)
        u = remap_MDAnalysis(u, topo)

        # Calculate RMSF
        rmsf = calculate_RMSF(u,i)
        rmsf_data.append(rmsf)

    rmsf = pd.concat(rmsf_data)

    rmsf.to_parquet(args.output)

if __name__ == '__main__':
    main()
    