#!/usr/bin/env python


import os
import argparse
import pathlib as path
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("-i", "--input", required=False, help="Define interaction input file, .parquet or .csv", default="/home/iman/code/squeeze-toy/demo-H08-MASP2-link/results/posco/posco_interactions.parquet")
    
    
    # Output
    parser.add_argument("-l", "--ligand_interaction", required=False, help="Define ligand analysis output file/directory, .svg", default="lig_interaction_sum3.svg")
    parser.add_argument("-r", "--receptor_interaction", required=False, help="Define receptor analysis output file/directory, .svg", default="rec_interaction_sum3.svg")

    """     
        parser.add_argument("-l", "--ligand", required=True, help="Define interaction partner to be considered (ligand)", default="")
        parser.add_argument("-r", "--receptor", required=True, help="Define interaction partner to be considered (receptor)", default="")
        
        parser.add_argument("-h", "--hbond", required=False, help="Define interaction type to be plotted (H-bonds)", default="")
        parser.add_argument("-n", "--lipophilic", required=True, help="Define interaction type to be plotted (lipophilic/nonpolar)", default="")
        parser.add_argument("-s", "--saltbridge", required=True, help="Define interaction type to be plotted (salt bridges)", default="")
    """

    return parser.parse_args()


def interaction_data_aggregation(interaction_partner, interaction_type):
    """"Filter, aggregate and pivot"""

    # Filter for interaction type

    if interaction_type == "Mar":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type)]
    if interaction_type != "total":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type)]
    else:
        df_interaction = df_filtered


"""     if interaction_type != "total":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type)]
    else:
        df_interaction = df_filtered
 """

    # Data aggregation. Aggreate over 
    resid = f'{interaction_partner}_resid'
    resname = f'{interaction_partner}_resname'

    # TODO: Normalize to numbe of replicates
    n = len(df_interaction.seed.unique())

    seed_avg = df_interaction.groupby([resid, resname, 'frame'])['Energy (e)'].mean().reset_index()
    seed_avg["residue_labels"] = seed_avg[resname] + ' ' + seed_avg[resid].astype(str)

    # Heatmap visualtions
    heatmap_data = pd.pivot_table(seed_avg, 
                             index=[resid, 'residue_labels'],  # Multi-index
                             columns='frame', 
                             values='Energy (e)')

    # Sort the index by residue number
    heatmap_data = heatmap_data.sort_index(level=resid)
    return heatmap_data

def plot_interactions(heatmap_data):

    # Set interaction specific parameters    
    if interaction_type == "total":
        interaction_cmap = "Greys"
        interaction_label = 'Total interaction Energy (kcal/mol)'
        interaction_min=18
        interaction_max=0
    elif interaction_type == "H-bond":
        interaction_cmap = "Blues"
        interaction_label = 'H-bond Energy (kcal/mol)'
        interaction_min=5
        interaction_max=0.25
    elif interaction_type == "lipophilic":
        interaction_cmap = "Oranges"
        interaction_label = 'Hydrophobic interaction Energy (kcal/mol)'
        interaction_min=0.37
        interaction_max=0.25
    elif interaction_type == "Salt bridge":
        interaction_cmap = "Greens"
        interaction_label = 'Salt Bridge interaction Energy (kcal/mol)'
        interaction_min=4.0
        interaction_max=0.25
    else:
        raise Exception("ERROR: Martin introduced a new interaction type")
        print("""Error: Interaction type found. Pls use "hbond", "lipo", or "salt" for H-bonds, hydrophobic interactions, and salt bridges respectively.""")


    sorted_labels = [label for (residue, label) in heatmap_data.index]

    
    ax=sns.heatmap(heatmap_data*-1, 
                cmap=interaction_cmap, 
                vmin=interaction_min, 
                vmax=interaction_max,
                yticklabels=sorted_labels,
                cbar_kws={
            'label': interaction_label,
            'pad': 0.016,  #space between heatmap and colorbar
                }
                )
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
        print("Failed to read file as .parquet")
        raise ValueError("Failed to read file")
    
    # 2. Data cleaning 
    # TODO: Perform a water analysis
    df_filtered = df[(df['receptor_resname'] != 'HOH') & (df['ligand_resname'] != 'HOH')]


    # 3. Data visualtisation

    figure_files = [args.ligand_interaction, args.receptor_interaction]

    for fig_file, interaction_type in zip(figure_files, ["ligand", "receptor"]):
        plt.figure(figsize=(15, 30))
        for i,interaction_type in enumerate(["total", 'H-bond', 'lipophilic', "Salt bridge"]): # , 'Marked as Salt-Bridge'
            interaction_pivot = interaction_data_aggregation('ligand', interaction_type)

            plt.subplot(3, 1, i+1)
            plot_interactions(interaction_pivot)

        plt.tight_layout()
        plt.savefig(fig_file)
        
        plt.close()




