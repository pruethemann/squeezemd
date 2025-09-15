#!/usr/bin/env python

import argparse
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("-i", "--input", required=False, help="Define interaction input file, .parquet or .csv", default="C:/Users/imata/Desktop/peter@gpu_server/H08--MASP2/analysis/posco_interactions.parquet")
    
    
    # Output
    parser.add_argument("-l", "--ligand_interaction", required=False, help="Define ligand analysis output file/directory, .svg", default="lig_heatmap.svg")
    parser.add_argument("-r", "--receptor_interaction", required=False, help="Define receptor analysis output file/directory, .svg", default="rec_heatmap.svg")

    return parser.parse_args()

def interaction_data_aggregation(interaction_partner, interaction_type):
    """"Filter, aggregate and pivot"""

    # filter for each interaction type
    # TODO: remove SB
    if interaction_type == "H-bond":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type) &
                                     (df_filtered['Marked as Salt-Bridge'] == 0)]
    elif interaction_type == "lipophilic":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type)]
    elif interaction_type == "Salt bridge":
        df_interaction = df_filtered[(df_filtered['Marked as Salt-Bridge'] == 1)]
    else:
        df_interaction = df_filtered

    # based on "observed" interaction partner
    resid = f'{interaction_partner}_resid'
    resname = f'{interaction_partner}_resname'

    # number of unique seeds for manual calculation of mean energy over seeds
    n_seeds = len(df_interaction.seed.unique())

    # data wrangling/aggregating for desired values
    # TODO: Generalise for Lig / Rec / mutation
    seed_avg = df_interaction.groupby([resid, resname, 'frame'])['Energy (e)'].sum().reset_index()
    seed_avg["Energy (e)"] = seed_avg["Energy (e)"].div(n_seeds)
    seed_avg["residue_labels"] = seed_avg[resname] + ' ' + seed_avg[resid].astype(str)

    # get maximum binding energy for cbar value limit
    emax = seed_avg["Energy (e)"].min()*-1
    #emax = emax.round(0)

    print(f"Interaction type: {interaction_type}, emax: {emax}")

    # create pivot for sns.heatmap
    heatmap_data = pd.pivot_table(seed_avg, 
                             index=[resid, 'residue_labels'],
                             columns='frame', 
                             values='Energy (e)')

    # sorting by resid, preceeded with resname
    heatmap_data = heatmap_data.sort_index(level=resid)
    return heatmap_data, emax

def plot_interactions(heatmap_data, emax):
    # plotting params based on interaction type
    # TODO: make vmax dynamic based on max interaction energy
    if interaction_type == "total":
        interaction_cmap = "Greys"
        interaction_label = 'Total interaction Energy (kcal/mol)'
        #interaction_min=18
        #interaction_max=0
    elif interaction_type == "H-bond":
        interaction_cmap = "Blues"
        interaction_label = 'H-bond Energy (kcal/mol)'
        #interaction_min=11
        #interaction_max=0.25
    elif interaction_type == "lipophilic":
        interaction_cmap = "Oranges"
        interaction_label = 'Hydrophobic interaction Energy (kcal/mol)'
        #interaction_min=2.5
        #interaction_max=0.25
    elif interaction_type == "Salt bridge":
        interaction_cmap = "Greens"
        interaction_label = 'Salt Bridge interaction Energy (kcal/mol)'
        #interaction_min=9
        #interaction_max=0.25
    else:
        raise Exception("ERROR: Martin introduced a new interaction type")

    sorted_labels = [label for (residue, label) in heatmap_data.index]

    ax=sns.heatmap(heatmap_data*-1, 
                cmap=interaction_cmap, 
                vmin=emax, 
                vmax=0,
                yticklabels=sorted_labels,
                cbar_kws={'label': interaction_label,
                          'pad': 0.012,  #space between heatmap and colorbar
                          })
    ax.set_title(f'Interaction Type: {interaction_type.capitalize()}', fontsize=12)
    ax.set_xlabel('Frame Number', fontsize=12)
    ax.set_ylabel('Residues', fontsize=12)

    # Further adjust colorbar label properties
    cbar = ax.collections[0].colorbar
    cbar.set_label(interaction_label, rotation=90, labelpad=10)

if __name__ == "__main__":
    # Get all arguments
    args = parse_arguments()

    # 1. Data import
    try:
        df = pd.read_parquet(args.input)
        print("File successfully read as .parquet")
    except Exception as parquet_error:
        raise ("Failed to read file")
    
    # 2. Data cleaning 
    # TODO: Perform a water analysis...?
    df_filtered = df[(df['receptor_resname'] != 'HOH') & (df['ligand_resname'] != 'HOH')]


    # 3. Data visualtisation
    figure_files = [args.ligand_interaction, args.receptor_interaction]

    for fig_file, interaction_partner in zip(figure_files, ["ligand", "receptor"]):
        plt.figure(figsize=(15, 30))
        for i,interaction_type in enumerate(["total", 'H-bond', 'lipophilic', "Salt bridge"]): # , 'Marked as Salt-Bridge'
            interaction_pivot, energy_max = interaction_data_aggregation(interaction_partner, interaction_type)
            
            plt.subplot(4, 1, i+1)
            plot_interactions(interaction_pivot, energy_max)

        plt.tight_layout()
        plt.savefig(fig_file)
        
        plt.close()
