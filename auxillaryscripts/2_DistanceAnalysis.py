#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
import argparse
from MDAnalysis.analysis import distances
import numpy as np
import MDAnalysis as mda


def get_distances(u, group_a, group_b, step=2):
    timeseries = []
    for ts in u.trajectory[::step]:
        # calculate distances between group_a and group_b
        distance = distances.distance_array(group_a,group_b, box=u.dimensions)

        timeseries.append([ts.frame, distance[0][0]])
    return np.array(timeseries)


def calculate_distances(args):
    # Just consider all topos and trajectors which have been centered
    topos = glob(f"MASP2_Helo_CM/WT/**/MD/topo_center.pdb", recursive=True)
    trajs = glob(f"MASP2_Helo_CM/WT/**/MD/traj_center.dcd", recursive=True)

    dataset = []

    for topo, traj in zip(topos, trajs):

        print(topo, traj)

        seed = topo.split("/")[-3]
        sim = topo.split("/")[-5]

        lig = "resid 92 and name CA and chainID A"  # Glu-92, Calpha
        rec = "resid 578 and name CZ and chainID B" # Arg-578, Calpha

        u = mda.Universe(topo, traj)

        # N terminus
        lig_grp = u.select_atoms(lig)
        rec_grp = u.select_atoms(rec)


        dists = get_distances(u, lig_grp, rec_grp)
        dists = pd.DataFrame(dists, columns=['time', 'distance'])

        dists['sim'] = sim
        dists['seed'] = seed

        #dists = pd.concat([dists_core, dists_N])
        dataset.append(dists)

    dataset = pd.concat(dataset)

    print(dataset)

    dataset.to_csv(args.distances)


if __name__ == '__main__':

    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Output
    parser.add_argument('--distances', required=False, default='distances.csv')

    args = parser.parse_args()

    calculate_distances(args)

