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
import plotly.express as px
import plotly.graph_objects as go
import argparse
import os
from glob import glob

sns.set(rc={'figure.figsize':(40,8.27)})

def generate_data(interactions:list):
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
            frame_ana = pd.read_csv(interaction_csv)
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

    stats['protein'] = protein

    return stats

def find_analysis_posco_files():
    posco_files = []

    for dir in args.dirs:
        posco_dir = os.path.join(dir, 'lig/')

        posco_files.extend(glob(posco_dir + '*.csv'))

    return (posco_files)

def plot_interaction_fingerprint(data, complex_protein:str, output, receptor, ligand):


    data_filter = data.loc[(complex_protein,'inter', receptor, ligand)]

        # Calculate mean and standard deviation for energy
    data_summary = data_filter.groupby(['resid', 'mutation']).agg(
                    energy_mean=('energy', 'mean'),
                    energy_std=('energy', 'std')
                    ).reset_index()

    fig = px.bar(data_summary,
                x='resid',
                y='energy_mean',
                color='mutation',
                error_y='energy_std',  # Use standard deviation as error bars
                title=f"{complex_protein}: Total per residue inter molecular interaction energy"
                )
    
        # Set the bars to be grouped horizontally rather than stacked
    fig.update_layout(xaxis_title="Residue ID",
                    yaxis_title="Energy",
                    xaxis_tickangle=-90,
                    barmode='group')  # Group bars for different mutations side by side

    fig.update_layout(xaxis_title="Residue ID",
                        yaxis_title="Energy",
                        xaxis_tickangle=-90)
    
    # Save figure as an HTML file
    fig.write_html(output)

def main(args):

    DEBUG = False
    if DEBUG:
        data = pd.read_feather(args.interactions)
        del data['Unnamed: 0']
    else:
        # Find all paths to Posco analysis files. One file per simutation and frame
        po_sco_files = find_analysis_posco_files()

        # 3. Create result table including: inter and intramolecular interactions / all simulations / all seeds / ligand and receptor perspective for ligand AND receptor
        data = generate_data(po_sco_files)

        # Merge and export ligand and receptor data
        data.to_parquet(args.interactions)
        data.to_csv('interactions.csv')

    # Exclude all waters in analysis. TODO perform a separate water analysis
    data = data[data.resname != 'HOH']

    # Aggregate energy over frames
    data_agg = data.groupby(['protein', 'interaction', 'target' , 'lig', 'mutation', 'resid', 'seed']).mean(numeric_only=True)
    del data_agg['frame']

    # TODO Extend for multiple different receptors
    # TODO extremly ugly
    receptor = data.target.unique()[0]
    ligand = data.lig.unique()[0]

    plot_interaction_fingerprint(data_agg, 'ligand', args.fingerprint_lig, receptor, ligand)
    plot_interaction_fingerprint(data_agg, 'receptor', args.fingerprint_rec, receptor, ligand)

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
    parser.add_argument('--fingerprint_rec', help="Location for per residue interaction analysis", required=False, default='residues.svg')
    parser.add_argument('--fingerprint_lig', help="Location for per residue interaction analysis", required=False, default='residues.svg')
    parser.add_argument('--totalEnergy', help="Location for total energy interaction analysis", required=False, default='total_energy.svg')


    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
