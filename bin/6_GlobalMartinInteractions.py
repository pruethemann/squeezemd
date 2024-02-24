#!/usr/bin/env python
"""
    Script performs descriptive analysis of the interaction analyizer of the last frames of Molecular dynamics simulations.

    This script can be used independent of snakemake since it detects all folders and seeds.

    1. Imports all interaction energies and merges them
    2. Aggregates over different features:
        - Seed
        - Ligand mutation
        - Residue id
    3. Visualizes the total energy
    4. Visualizes the energy per residue and mutation
    5. Visualizes the energy differences between wildetype and mutation per residue

    Data variable description:
    Group by:
        protein: ligand / receptor
        interaction: inter, intra
        target: receptor (C1s)
        lig: (BD001)
        mutation: WT / Y119E
    Take Mean:
        frame: 1:100
        interaction type: hydrophobic, electrostatic, ..
    Get SD:
        seed: Seed of MD

"""

import pandas as pd
from os import path
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
from glob import glob

sns.set(rc={'figure.figsize':(40,8.27)})

def generate_data(interactions:list, protein='ligand'):
    """
    Imports all the single data files and merges them together to 1 dataframe
    :param sim_df:
    :param frame_number:
    :param data_out:
    :return:
    """

    # Import all interaction analyizer data and combine
    stats = []

    # TODO use metadata from parquet or directly from snakemake
    for interaction_csv in interactions:

        # Extract meta data like an idiot
        metadata = interaction_csv.split('/')

        complex = metadata[-6]
        mutation = metadata[-5]
        target = complex.split('_')[0]
        ligand = complex.split('_')[1]
        frame_id = int(metadata[-1][:-4])
        seed = int(metadata[-4])

        # Import data
        try:
            frame_ana = pd.read_csv(interaction_csv, names=['interaction', 'resname', 'resid', 'energy'])
        except FileNotFoundError:
            print("Error with import from: ", interaction_csv)
            continue

        # Determine metrics lables
        # TODO: Do the same as in 10_interactionsurface and take these values from snakemake
        frame_ana['name'] = complex
        frame_ana['target'] = target
        frame_ana['lig'] = ligand
        frame_ana['mutation'] = mutation
        frame_ana['frame'] = frame_id
        frame_ana['seed'] = seed

        # Save data
        stats.append(frame_ana)


    # Merge all data together
    stats = pd.concat(stats)

    # Determine inter and intramolecular interactions
    stats['interaction type'] = stats.interaction
    stats.loc[stats['interaction type'].str.contains('L-P'), 'interaction'] = 'inter'
    stats.loc[stats['interaction type'].str.contains('L-L'), 'interaction'] = 'intra'

    stats['protein'] = protein

    return stats

def find_analysis_posco_files():
    ligand_files = []
    receptor_files = []
    for dir in args.dirs:
        lig_dir = os.path.join(dir, 'lig/')
        rec_dir = os.path.join(dir, 'rec/')
        ligand_files.extend(glob(lig_dir + '*.csv'))
        receptor_files.extend(glob(rec_dir+ '*.csv'))

    return (ligand_files, receptor_files)

def main(args):

    DEBUG = False
    if DEBUG:
        data = pd.read_feather(args.interactions)
        del data['Unnamed: 0']
    else:
        # Find all paths to Posco analysis files. One file per simutation and frame
        (ligand_files, receptor_files) = find_analysis_posco_files()

        # 3. Create result table including: inter and intramolecular interactions / all simulations / all seeds / ligand and receptor perspective for ligand AND receptor
        data_ligand = generate_data(ligand_files, protein='ligand')
        data_receptor = generate_data(receptor_files, protein='receptor')

        # Merge and export ligand and receptor data
        data = pd.concat([data_ligand, data_receptor])
        data.to_parquet(args.interactions)
        data.to_csv('interactions.csv')

    # Exclude all waters in analysis. TODO perform a separate water analysis
    data = data[data.resname != 'HOH']

    # Aggregate energy over frames
    data_agg = data.groupby(['protein', 'interaction', 'target' , 'lig', 'mutation', 'resid', 'seed']).mean(numeric_only=True)
    del data_agg['frame']

    # TODO Extend for multiple different receptors
    receptor = data.target.unique()[0]
    ligand = data.lig.unique()[0]

    # Per residue interaction: ligand vs receptor
    for i,protein in enumerate(['ligand', 'receptor']):
        plt.subplot(2,1,i+1)
        data_filter = data_agg.loc[(protein,'inter', receptor, ligand)]

        sns.barplot(data=data_filter,
                    x='resid',
                    y='energy',
                    hue='mutation',
                    errorbar='sd',
                    capsize=0.1
                    )

        plt.xticks(rotation=90)
        plt.title(f"{protein}: Total per residue inter molecular interaction energy")

        # Save aggregated data
        data_filter.to_csv(f'{args.residueEnergy[:-4]}_{protein}.csv')

    plt.savefig(args.residueEnergy)
    plt.close()

    # Total Energy
    total_energy = data.groupby(['interaction', 'protein', 'target', 'lig', 'mutation']).sum(numeric_only=True)
    total_energy.drop(columns=['resid','frame','seed'], inplace=True)
    #total_energy.to_csv(path.join(args.output, 'total.csv'))

    sns.barplot(data=total_energy.loc['inter'],
                x='mutation',
                y='energy',
                hue='protein'
                )
    plt.xticks(rotation=90)
    plt.title("Total inter molecular interaction energy")
    plt.savefig(args.totalEnergy)
    plt.close()

def parse_arguments():
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--dirs', nargs='+', required=False)

    # Output
    parser.add_argument('--interactions', required=False, default='dev3/interactions.feather')
    parser.add_argument('--residueEnergy', help="Location for per residue interaction analysis", required=False, default='residues.svg')
    parser.add_argument('--totalEnergy', help="Location for total energy interaction analysis", required=False, default='total_energy.svg')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
