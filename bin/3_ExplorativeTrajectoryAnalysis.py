#!/usr/bin/env python

"""
    This module analyses typical properties of a MD trajectory:
    - RMSF of receptor and ligand
    - RMSD of recetpor and ligand
    - Statistics of MD properties like T,p,energies

"""

import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

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

    plt.title("Temperature (K) Mean: " + str(data['Temperature (K)'].mean().round(2)))

    # Volume
    plt.subplot(2,2,4)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Box Volume (nm^3)')

    plt.title("Box Volume (nm^3) Mean: " + str(data['Box Volume (nm^3)'].mean().round(2)))

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(output_graph)


def calculate_RMSF(u: mda.Universe, args):
    """
    Calculate the RMSF of the ligand and the receptor

    info:
    https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html
    :param u:
    :param args:
    :return:
    """

    chains = {'ligand': 'I',
              'receptor': 'B'}

    # TODO: Check where the renaming of the chains is coming from!!
    chains = {'ligand': 'C',
              'receptor': 'A'}

    print("Init RMSF analysis")

    for i,(protein, chain) in enumerate(chains.items()):
        c_alphas = u.select_atoms(f'chainID {chain} and name CA')
        R = rms.RMSF(c_alphas).run()

        plt.subplot(2,1,i+1)

        plt.plot(c_alphas.resids, R.results.rmsf)
        plt.xlabel('Residue number')
        plt.ylabel('RMSF ($\AA$)')
        plt.title(f'RMSF {protein}')

        RMSF_df = pd.DataFrame([c_alphas.resids, R.rmsf])
        RMSF_df.to_csv(args.rmsf[-4] + '.csv')

    plt.tight_layout()
    plt.savefig(args.rmsf)
    plt.close()

    # Calculate bfactors
    c_alphas = u.select_atoms('protein and name CA')
    R = mda.analysis.rms.RMSF(c_alphas).run()
    calculate_bfactors(R)


def calculate_bfactors(R):
    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms
    for residue, r_value in zip(protein.residues, R.rmsf):
        residue.atoms.tempfactors = r_value

    u.atoms.write(args.bfactors)

def calculate_RMSD(u: mda.Universe, args):
    """
    Calculate RMSD of receptor and ligand
    :param u:
    :param output:
    :return:
    """

    print("Init RMSD analysis")

    ligand = f"backbone and chainID I"
    target = f"backbone and chainID B"

    R = mda.analysis.rms.RMSD(u,  # universe to align
                 u,  # reference universe or atomgroup
                 select='backbone',  # group to superimpose and calculate RMSD
                 groupselections=[ ligand, target],  # groups for RMSD
                 ref_frame=0)  # frame index of the reference
    R.run()

    df = pd.DataFrame(R.results.rmsd,
                      columns=['Frame', 'Time (ns)',
                               'Backbone',
                               'ligand', 'target'])

    ax = df.plot(x='Frame', y=['Backbone', 'ligand', 'target'],
                 kind='line')
    ax.set_ylabel(r'RMSD ($\AA$)')

    df.to_csv(args.rmsd[-4:] + '.csv')
    plt.savefig(args.rmsd)
    plt.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False, help='Topo file', default='output/WT-simple/topo.pdb')
    parser.add_argument('--traj', required=False, help='Trajectory', default='output/WT-simple/traj_150.dcd')
    parser.add_argument('--stats', required=False, help='MD Stats file')

    # Output
    parser.add_argument('--rmsf', required=False, help='')
    parser.add_argument('--bfactors', required=False, help='')
    parser.add_argument('--rmsd', required=False, help='')
    parser.add_argument('--fig_stats', required=False, help='')

    args = parser.parse_args()

    # Import Trajectory
    u = mda.Universe(args.topo, args.traj, in_memory=False)
    traj_length = len(u.trajectory)

    # Calculate multiple MD trajectery properties
    calculate_RMSF(u, args)
    calculate_RMSD(u, args)

    # Visualize Energies, T, ...
    visualize_MDStats(args.stats, args.fig_stats)
