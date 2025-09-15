#!/usr/bin/env python

import os
import argparse
import pathlib as path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from glob import glob

def parse_arguments():
    parser = argparse.ArgumentParser()

    #LINUX PATHS
    # Input
    parser.add_argument("-i", "--input", required=False, help="Define interaction input file, .parquet or .csv", default="/home/iman/caracara/MD/squeeze_MD/S-01_H08_MASP2_30ns/results/posco/posco_interactions.parquet")
    #parser.add_argument("-s", "--sequence", required=False, help="Define sequence range of ligand/receptor for efficient plotting of barplot. Read from sequence.parquet", default="/home/iman/caracara/MD/squeeze_MD/S-01_H08_MASP2_30ns/MASP2_H08/WT/131/frames/sequence.parquet")
    
    # Output
    parser.add_argument("-l", "--ligand_interaction", required=False, help="Define ligand analysis output file/directory, .svg", default="lig_barplot.svg")
    parser.add_argument("-r", "--receptor_interaction", required=False, help="Define receptor analysis output file/directory, .svg", default="rec_barplot.svg")

    return parser.parse_args()

def import_sequence_range(seq_parquet:os.path, protein:str):

    seq_df = pd.read_parquet(seq_parquet).reset_index()
    seq_df = seq_df[(seq_df['protein'] == protein)]
    
    return (seq_df.resid.min(), seq_df.resid.max())

def interaction_data_aggregation(interaction_partner, interaction_type, df_filtered):
    """"Filter, aggregate and pivot"""

    # filter for each interaction type
    if interaction_type == "H-bond":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type) &
                                     (df_filtered['Marked as Salt-Bridge'] == 0)]
    elif interaction_type == "lipophilic":
        df_interaction = df_filtered[(df_filtered['Interaction Type'] == interaction_type)]
    elif interaction_type == "Salt bridge":
        df_interaction = df_filtered[(df_filtered['Marked as Salt-Bridge'] == 1)]
    else:
        df_interaction = df_filtered

    
    # TODO Extremely dirty to get access to sequence
    seq_path = glob(f"**/{mutation}/**/frames/sequence.parquet")

    # based on "observed" interaction partner
    try:
        seq_range = import_sequence_range(seq_path[0], interaction_partner[:3])
        seq_range = range(seq_range[0], seq_range[1])
    except Exception:
        raise Exception("Error: Interaction partner not found.")
    
    resid = f'{interaction_partner}_resid'

    # number of unique seeds for manual calculation of mean energy over seeds
    n_frames = len(df_interaction.frame.unique())
    n_seeds = len(df_interaction.seed.unique())

    # data wrangling/aggregating for desired values, leaving frames
    frame_avg = df_interaction.groupby([resid, "seed"])['Energy (e)'].sum().reset_index()
    frame_avg["Energy (e)"] = frame_avg["Energy (e)"].div(n_frames)

    # data wrangling/aggregating for desired values, leaving seeds
    seed_avg = frame_avg.groupby([resid])['Energy (e)'].sum().reset_index()
    seed_avg["Energy (e)"] = seed_avg["Energy (e)"].div(n_seeds)
    seed_avg.rename(columns={"Energy (e)": "mean"}, inplace=True)

    seed_sd = frame_avg.groupby([resid])['Energy (e)'].std().reset_index()
    seed_sd.rename(columns={"Energy (e)": "sd"}, inplace=True)

    combined = pd.merge(seed_avg, seed_sd, on=resid, how="outer")
    
    all_resid = pd.DataFrame({resid: seq_range})
    final = pd.merge(combined, all_resid, on=resid, how="left")

    # get maximum binding energy for cbar value limit
    emax = seed_avg["mean"].min()
    
    return final, emax

def plot_interactions(plot_data, emax):
    # plotting params based on interaction type
    # TODO: make vmax dynamic based on max interaction energy
    if interaction_type == "total":
        interaction_color = "grey"
        interaction_label = f'Total interaction Energy (kcal/mol) ({mutation})'
    elif interaction_type == "H-bond":
        interaction_color = "dodgerblue"
        interaction_label = f'H-bond Energy (kcal/mol) ({mutation})'
    elif interaction_type == "lipophilic":
        interaction_color = "darkorange"
        interaction_label = 'Hydrophobic interaction Energy (kcal/mol) ({mutation})'
    elif interaction_type == "Salt bridge":
        interaction_color = "seagreen"
        interaction_label = 'Salt Bridge interaction Energy (kcal/mol) ({mutation})'
    else:
        raise Exception("ERROR: Martin introduced a new interaction type")


    # TODO Extremely dirty to get access to sequence
    seq_path = glob(f"**/{mutation}/**/frames/sequence.parquet")

    seq_range = import_sequence_range(seq_path[0], interaction_partner[:3])
    plot_range = range(seq_range[0], seq_range[1], 2)

    resid = f'{interaction_partner}_resid'

    plt.bar(x=plot_data[resid],
            height=plot_data["mean"],
            yerr=plot_data["sd"],
            color=interaction_color,
            )
    
    plt.title(interaction_label)
    #plt.xlabel(plot_data["ligand_resid"])
    plt.xticks(plot_range, fontsize=10)
    
    plt.ylabel('Mean Interaction Energy (kcal/mol)')
    #plt.yticks(np.linspace(0, emax-1, 10), fontsize=10)
    plt.axhline(y=0, color="black", linewidth=0.8)
    #plt.ylim(interaction_min, 0)


if __name__ == "__main__":
    # Get all arguments
    args = parse_arguments()

    # 1. Data import
    try:
        df = pd.read_parquet(args.input)
        print("File successfully read as .parquet")
    except Exception as parquet_error:
        raise ValueError("Failed to read file")
    
    # 2. Data cleaning 
    # TODO: Perform a water analysis...?
    df_filtered = df[(df['receptor_resname'] != 'HOH') & (df['ligand_resname'] != 'HOH')]

    # Iterate over mutations
    mutations = df_filtered.mutation.unique()

    for mutation in mutations:

        data = df_filtered[df_filtered.mutation == mutation]

        # 3. Data visualisation
        figure_files = [args.ligand_interaction, args.receptor_interaction]
        for fig_file, interaction_partner in zip(figure_files, ["ligand", "receptor"]):
            plt.figure(figsize=(15, 30))
            for i,interaction_type in enumerate(["total", 'H-bond', 'lipophilic', "Salt bridge"]): # , 'Marked as Salt-Bridge'
                final, energy_max = interaction_data_aggregation(interaction_partner, interaction_type, data)
                plt.subplot(4, 1, i+1)
                plot_interactions(final, energy_max)

            plt.tight_layout()
            if mutation == 'WT':
                plt.savefig(fig_file)
            else:
                plt.savefig(f'{fig_file[:-4]}_{mutation}.svg')

            plt.close()
