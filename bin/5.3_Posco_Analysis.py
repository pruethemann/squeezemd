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
import seaborn as sns
import plotly.express as px
import argparse

sns.set(rc={'figure.figsize':(40,8.27)})

def plot_interaction_fingerprint(data, figure):

    data['ligand_resid'] = data['ligand_resid'].astype(int)

    # TODO: Exclude Water data
    data = data[data.ligand_resid <= 122]

    fig = px.bar(data,
                x='ligand_resid',
                y='energy_mean',
                color='mutation',
                error_y='energy_std'
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
    fig.write_html(figure)

def main(args):

    data = pd.read_parquet(args.input)

    # TODO Exclude all waters in analysis. TODO perform a separate water analysis
    data = data[data.ligand_resname != 'HOH']
    data = data[data.receptor_resname != 'HOH']

    # Sum all energy contributions for every resid
    resid_energy =  data.groupby(['Interaction Type', 'target' , 'lig', 'mutation','seed','frame', 'ligand_resid']).sum(numeric_only=True).reset_index()

    # Take the mean of all residue contributions for all frames
    frame_mean =  resid_energy.groupby(['Interaction Type', 'target' , 'lig', 'mutation','seed', 'ligand_resid']).agg(energy_mean=('Energy (e)', 'mean'),energy_std=('Energy (e)', 'std')).reset_index()

    DEBUG = True
    if DEBUG:
        resid_energy.to_csv('resid_energy.csv')
        frame_mean.to_csv('frame_mean.csv')

    plot_interaction_fingerprint(frame_mean,args.output)


def parse_arguments():
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--input')

    # Output
    parser.add_argument('--output')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    main(args)

