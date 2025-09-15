#!/usr/bin/env python

import argparse
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import openmm.app as app
from Helper import remap_MDAnalysis

def calculate_RMSF(u: mda.Universe, i):

    c_alphas = u.select_atoms(f'chainID A and name CA')
    R = rms.RMSF(c_alphas).run()

    # Store RMSF and secondary structure data
    rmsf_df = {'resid':c_alphas.resids, 'rmsf': R.results.rmsf, 'sim_id': i}

    return pd.DataFrame(rmsf_df)


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', nargs='+', required=False)
    parser.add_argument('--traj', nargs='+', required=False)

    # Output
    parser.add_argument('--output', required=False, default='rmsf.svg', help='')

    return parser.parse_args()

# Example of running the function
if __name__ == '__main__':

    args = parse_arguments()

    topos = sorted(args.topo)
    trajs = sorted(args.traj)

    rmsf_data = []

    for i,(topo,traj) in enumerate(zip(topos,trajs)):

        # Import Trajectory
        topo = app.PDBxFile(topo)
        u = mda.Universe(topo, traj, in_memory=False)
        u = remap_MDAnalysis(u, topo)

        # calcualte RMSDF
        rmsf = calculate_RMSF(u,i)
        rmsf_data.append(rmsf)

    rmsf = pd.concat(rmsf_data)

    rmsf.to_parquet(args.output)
