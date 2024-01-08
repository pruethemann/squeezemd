#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse


def import_data(fingerprints):

    # Import all interaction analyizer data and combine
    data = []

    for fp_feather in fingerprints:

        # Extract meta data like an idiot
        metadata = fp_feather.split('/')

        complex = metadata[-5]
        mutation = metadata[-4]
        target = complex.split('_')[0]
        ligand = complex.split('_')[1]
        seed = int(metadata[-3])

        # Import data
        try:
            fp = pd.read_feather(fp_feather)
            fp = data_engineering(fp)
        except FileNotFoundError:
            print("Error with import from: ", fp_feather)
            continue

        # Determine metrics lables
        fp['name'] = complex
        fp['target'] = target
        fp['lig'] = ligand
        fp['mutation'] = mutation
        fp['seed'] = seed

        fp.to_csv('2_addion.csv')
        # Save data
        data.append(fp)

    # Merge all data together
    data = pd.concat(data)
    data.to_csv('3_total.csv')
    return data

def data_engineering(data):
    n_frames = 100 # TODO data.index.max() + 1

    # Aggregate all interactions
    data_agg = data.sum()
    data_agg = data_agg.reset_index()

    # rename columns
    data_agg.columns = ['ligand', 'target', 'interaction', 'sum']

    data_agg['sum'] /= n_frames

    # Group interactions
    interaction_map = {
        'Cationic': 'salt bridge',
        'Anionic': 'salt bridge',
        'HBAcceptor': 'H bonds',
        'HBDonor': 'H bonds',
        'PiStacking': 'PiStacking',
        "PiCation": 'Cation-Pi',
        "CationPi": 'Cation-Pi',
        "Hydrophobic": "Hydrophobic"
    }


    # interaction_types = ['salt brdige', 'H bonds']

    data_agg['interaction_type'] = data_agg['interaction'].map(interaction_map)

    # extract resids
    data_agg['resid'] = data_agg['ligand'].str.extract('(\d+)').astype(int)
    return data_agg

def visualize(fingerprints, mutation):
    sns.set(rc={'figure.figsize':(25,30)})

    n_residues = 122

    # Group all interactions per residue
    interaction_types = ['salt bridge', 'H bonds', 'PiStacking', 'Hydrophobic', 'Cation-Pi']

    data = fingerprints.reset_index().copy()

    data = data[data.mutation == mutation]
    x = np.arange(1, n_residues + 1)

    for i,interaction_type in enumerate(interaction_types):
        # Add missing interactions to get a even distribution
        data_interaction = data[data.interaction_type == interaction_type]

        data_interaction = data_interaction.groupby(['resid', 'seed']).agg({'interaction_type': 'first', 'sum':'sum'} ).reset_index()

        data_interaction.to_csv('2_seed.csv')

        dummy_df = pd.DataFrame(x, columns=['resid'])
        dummy_df['interaction_type'] = interaction_type
        dummy_df['sum'] = 0

        data_interaction = pd.concat([data_interaction, dummy_df])

        # Get rid of duplicates
        data_interaction.drop_duplicates(subset=['resid','seed'], keep='first', inplace=True)

        plt.subplot(5, 1, i+1)
        sns.barplot(data=data_interaction,
                    x='resid',
                    y='sum'
                    )

        plt.title(interaction_type)
        plt.xticks(rotation=90)

    plt.savefig(args.figure)
    plt.savefig(args.figure[:-4] + mutation + '.svg')
    #plt.show()

if __name__ == '__main__':
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--fingerprints', nargs='+', required=False)

    # Output
    parser.add_argument('--interactions', required=False, default='')
    parser.add_argument('--figure', required=False, default='')

    args = parser.parse_args()
    #job_id = args.output.split('/')[-3]

    # Import all fingerprints data
    fingerprints = import_data(args.fingerprints)

    for mutation in fingerprints.mutation.unique():
        visualize(fingerprints, mutation)

    # Export data
    fingerprints.to_csv(args.interactions)
