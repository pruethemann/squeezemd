#!/usr/bin/env python
import MDAnalysis
from MDAnalysis.analysis import contacts
import os
import MDAnalysis.transformations as trans
import argparse
import subprocess
from pathlib import Path
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def visualize_MDStats(stats_file, output_graph):
    data = pd.read_csv(stats_file, sep='\t')

    data['time (ns)'] = data['Time (ps)'] / 1000

    # Total Energy
    plt.subplot(2,2,1)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Total Energy (kJ/mole)')

    plt.title("Total Energy")

    # Potential Energy
    plt.subplot(2,2,2)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Potential Energy (kJ/mole)')

    plt.title("Potential Energy (kJ/mole)")

    # Temperature
    plt.subplot(2,2,3)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Temperature (K)')

    plt.title("Temperature (K) Mean: " + str(data['Temperature (K)'].mean()))

    # Volume
    plt.subplot(2,2,4)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Box Volume (nm^3)')

    plt.title("Box Volume (nm^3) Mean: " + str(data['Box Volume (nm^3)'].mean()))

    plt.savefig(output_graph)

def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

def calculate_RMSF(u: MDAnalysis.Universe, output):

    print("Init RMSF analysis")

    # RMSF: https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html
    average = align.AverageStructure(u, u, select='protein and name CA',
                                     ref_frame=0).run()
    ref = average.universe

    c_alphas = u.select_atoms('protein and name CA')
    R = mda.analysis.rms.RMSF(c_alphas).run()

    plt.plot(c_alphas.resids, R.results.rmsf)
    plt.xlabel('Residue number')
    plt.ylabel('RMSF ($\AA$)')
    plt.legend()

    plt.savefig(f'{output}analysis/RMSF.svg')

    RMSF_df = pd.DataFrame([c_alphas.resids, R.rmsf])
    RMSF_df.to_csv(f'{output}/analysis/RMSF.csv')

    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms
    for residue, r_value in zip(protein.residues, R.rmsf):
        residue.atoms.tempfactors = r_value

    u.atoms.write(f'{output}/analysis/bfactors.pdb')

def calculate_RMSD(u: MDAnalysis.Universe, output):

    print("Init RMSD analysis")

    print(u)

    BD001 = f"backbone and chainID I"
    target = f"backbone and not chainID I"

    R = mda.analysis.rms.RMSD(u,  # universe to align
                 u,  # reference universe or atomgroup
                 select='backbone',  # group to superimpose and calculate RMSD
                 groupselections=[ BD001, target],  # groups for RMSD
                 ref_frame=0)  # frame index of the reference
    R.run()

    df = pd.DataFrame(R.results.rmsd,
                      columns=['Frame', 'Time (ns)',
                               'Backbone',
                               'BD001', 'target'])



    ax = df.plot(x='Frame', y=['Backbone', 'BD001', 'target'],
                 kind='line')
    ax.set_ylabel(r'RMSD ($\AA$)')

    df.to_csv(f'{output}/analysis/RMSD.csv')
    plt.savefig(f'{output}/analysis/RMSD.svg')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False, help='Topo file', default='output/WT-simple/topo.pdb')
    parser.add_argument('--traj', required=False, help='Trajectory', default='output/WT-simple/traj_150.dcd')
    parser.add_argument('--stats', required=False, help='MD Stats file')
    parser.add_argument('--mapping', required=False, help='working directory')

    # Output
    parser.add_argument('--dir', required=False, help='working directory')
    parser.add_argument('--ana', required=False, help='working directory')
    parser.add_argument('--checkpoint', required=False, help='Checkpoint file path')

    args = parser.parse_args()

    # Import Trajectory
    u = mda.Universe(args.topo, args.traj, in_memory=False)
    traj_length = len(u.trajectory)

    ## Calculate multiple MD trajectery properties
    calculate_RMSF(u, args.dir)                                                     # RMSF
    calculate_RMSD(u, args.dir)                                              # RMSD

    visualize_MDStats(args.stats, os.path.join(args.dir, 'MD', 'MDStats.png'))      # Visualize Energies, T, ...

    # Signal to snakemake that the process completed successfully
    Path(args.checkpoint).touch()
