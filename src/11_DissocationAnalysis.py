#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
import argparse
from MDAnalysis.analysis import distances
import numpy as np
import MDAnalysis as mda

def rmsd_analysis(args):
    # Let's do it the old-fashioned way and just consider all possible RMSD files
    rmsd_files = glob(f"output/{args.job_id}/**/**/**/RMSD.csv")

    print(rmsd_files)

    data = []

    for rmsd in rmsd_files:
        print(rmsd)

        seed = rmsd.split('/')[-3]
        target = rmsd.split('/')[-4]

        rmsd_df = pd.read_csv(rmsd)

        rmsd_df['rec'] = target
        rmsd_df['seed'] = seed

        data.append(rmsd_df)

    data = pd.concat(data)

    del data['Unnamed: 0']

    data['Frame'] = data['Frame'].astype(int)

    print(data.info())

    # Consider only every 25th frame
    data = data[data.Frame %25  == 0]

    data.reset_index(inplace=True)

    sns.lineplot(data=data,
                 x="Time (ns)",
                 y="BD001",
                 hue="rec")

    data.to_csv(args.rmsd)
    plt.savefig(args.rmsd_png)
    plt.show()



def get_distances(u, group_a, group_b):
    timeseries = []
    for ts in u.trajectory[::2]:
        # calculate distances between group_a and group_b
        distance = distances.distance_array(group_a,group_b, box=u.dimensions)

        timeseries.append([ts.frame, distance[0][0]])
    return np.array(timeseries)

# 20 ns simulations

def calculate_distances(args):
    # Just consider all topos and trajectors which have been centered
    topos = glob(f"C1s_BD001/**/**/MD/center/topo_center.pdb", recursive=True)
    trajs = glob(f"C1s_BD001/**/**/MD/center/traj_center.dcd", recursive=True)

    dataset = []

    for topo, traj in zip(topos, trajs):

        print(topo, traj)

        seed = topo.split("/")[-4]
        sim = topo.split("/")[-5]

        trp = "resid 17 and name CA and chainID A"
        gly = "resid 440 and name CZ and chainID B"

        """
        if "C1s" in sim:
            enzyme = 'resid 526 and name CA'
            asp = 'resid 572 and name CG'
        elif "C1r" in sim:
            enzyme = 'resid 526 and name CA'
            asp = 'resid 572 and name CG'
        elif "FactorXa" in sim:
            enzyme = 'resid 526 and name CA'
            asp = 'resid 572 and name CG'
        elif "MASP2" in sim:
            enzyme = 'resid 526 and name CA'
            asp = 'resid 572 and name CG'
        elif "PlasmaKallikrein" in sim:
            enzyme = 'resid 526 and name CA'
            asp = 'resid 572 and name CG'
        else:
            print(sim)
            print("FAIL")

        """

        u = mda.Universe(topo, traj)

        # N terminus
        trp_grp = u.select_atoms(trp)
        enzyme_grp = u.select_atoms(enzyme)


        dists_N = get_distances(u, trp_grp, enzyme_grp)
        dists_N = pd.DataFrame(dists_N, columns=['time', 'distance'])

        dists_N['sim'] = sim
        dists_N['seed'] = seed
        dists_N['dist'] = "N"


        # core distances
        group1 = u.select_atoms(trp_grp)
        group2 = u.select_atoms(enzyme_grp)

        dists_core = get_distances(u, group1, group2)
        dists_core = pd.DataFrame(dists_core, columns=['time', 'distance'])

        dists_core['sim'] = sim
        dists_core['seed'] = seed
        dists_core['dist'] = "core"

        dists = pd.concat([dists_core, dists_N])
        dataset.append(dists)

    dataset = pd.concat(dataset)

    print(dataset)

    dataset.to_csv(args.distances)


if __name__ == '__main__':

    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Output
    parser.add_argument('--rmsd', required=False, default='rmsd.csv')
    parser.add_argument('--rmsd_png', required=False, default='rsmd.png')
    parser.add_argument('--distances', required=False, default='distances.png')

    args = parser.parse_args()

    calculate_distances(args)

    #rmsd_analysis(args)
