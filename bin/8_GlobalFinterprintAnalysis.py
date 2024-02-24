#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse


def import_data(fingerprints):
    """
        Imports and processes fingerprint data from a list of file paths.

        Parameters:
        - fingerprints: A list of file paths to fingerprint data files.

        Returns:
        - A combined DataFrame containing processed fingerprint data from all files.

        Note: This function relies on 'data_engineering' to process individual DataFrames.
        """

    ## Initialize an empty list to store processed data
    combined_data = []

    for fp_path in fingerprints:

        # Import data
        try:
            fp = pd.read_parquet(fp_path)
            fp = data_engineering(fp, args.n_frames)
        except FileNotFoundError:
            print("Error with import from: ", fp_path)
            continue

        # Determine metrics lables
        fp['name'] = fp.attrs['complex']
        fp['target'] = fp.attrs['target']
        fp['lig'] = fp.attrs['ligand']
        fp['mutation'] = fp.attrs['mutation']
        fp['seed'] = fp.attrs['seed']

        # Save data
        combined_data.append(fp)

    # Merge all data together
    data = pd.concat(combined_data)
    data.attrs['n_residues'] = fp.attrs['n_residues_ligand']
    return data

def data_engineering(data, n_frames):
    """
    Processes individual fingerprint data.

    Parameters:
    - fp_data: DataFrame containing fingerprint data.
    - n_frames: An integer indicating the number of frames to process.

    Returns:
    - Processed DataFrame.
    """

    # Aggregate all interactions
    data_agg = data.sum()
    data_agg = data_agg.reset_index()

    # rename columns
    data_agg.columns = ['ligand', 'target', 'interaction', 'sum']
    data_agg['sum']  /= n_frames

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

    data_agg['interaction_type'] = data_agg['interaction'].map(interaction_map)

    # extract resids
    data_agg['resid'] = data_agg['ligand'].str.extract('(\d+)').astype(int)
    return data_agg

def visualize_data(fingerprints, mutation):
    """
    Generates visualizations for the processed fingerprint data.

    Parameters:
    - data: DataFrame containing processed fingerprint data.
    """
    sns.set(rc={'figure.figsize':(25,30)})

    n_residues= fingerprints.attrs['n_residues']

    # Group all interactions per residue
    interaction_types = ['salt bridge', 'H bonds', 'PiStacking', 'Hydrophobic', 'Cation-Pi']

    data = fingerprints.reset_index().copy()

    data = data[data.mutation == mutation]

    for i,interaction_type in enumerate(interaction_types):
        # Add missing interactions to get a even distribution
        data_interaction = data[data.interaction_type == interaction_type]

        data_interaction = data_interaction.groupby(['resid', 'seed']).agg({'interaction_type': 'first', 'sum':'sum'} ).reset_index()

        #Debug
        #data_interaction.to_csv(f'{mutation}_{interaction_type}.seed.csv')

        # Generate a DataFrame with all residue numbers between 1 and 100
        all_resid = pd.DataFrame({'resid': range(1, n_residues + 1),
                                  'interaction_type': interaction_type})

        # Merge the original DataFrame with the complete list of residues
        # Use left join to keep all residues and fill missing energies with -1
        data_interaction = pd.merge(all_resid, data_interaction, on='resid', how='left').fillna(0)

        plt.subplot(5, 1, i+1)
        sns.barplot(data=data_interaction,
                    x='resid',
                    y='sum'
                    )
        plt.title(interaction_type)
        plt.xticks(rotation=90)

    # TODO: Clean up paths
    plt.savefig(f'{args.interactions}_{mutation}.svg')
    plt.close()


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--fingerprints', nargs='+', help="List of fingerprint parquet files", required=True, default=['fp1.parquet', 'fp2.parquet'])
    parser.add_argument('--n_frames', help='How many frames to be analysed. Only the last n frames will be analyzed. Default 100', type=int, required=False, default=100)

    # Output
    parser.add_argument('--interactions', help="Joined fingerprint csv file path", required=False, default='fingerprint.csv')
    parser.add_argument('--figure', help="Fingerprints figure", required=False, default='fingerprint.svg')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    # Import all fingerprints data
    fingerprints = import_data(args.fingerprints)


    for mutation in fingerprints.mutation.unique():
        visualize_data(fingerprints, mutation)

    # Export data
    fingerprints.to_csv(args.interactions)
