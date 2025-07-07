#!/usr/bin/env python

import os
import argparse
import pathlib as path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def parse_arguments():
    parser = argparse.ArgumentParser()

    #LINUX PATHS
    # Input
    parser.add_argument("-i", "--input", required=False, help="Define interaction input file, .parquet or .csv", default="/home/iman/code/squeeze-toy/notebooks/posco_interactions.parquet")
    # Output
    parser.add_argument("-l", "--ligand_interaction", required=False, help="Define ligand analysis output file/directory, .svg", default="/home/iman/code/squeeze-toy/lig_barplot.svg")
    parser.add_argument("-r", "--receptor_interaction", required=False, help="Define receptor analysis output file/directory, .svg", default="/home/iman/code/squeeze-toy/rec_barplot.svg")

    return parser.parse_args()

def interaction_data_aggregation(interaction_partner, interaction_type):
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

    # based on "observed" interaction partner
    if interaction_partner == "ligand":
        seq_range = range(1, 113)
        plot_range = range(1, 113, 10)
    elif interaction_partner == "receptor":
        seq_range = range(445, 684)
        plot_range = range(445, 684, 20)
    else:
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
        interaction_label = 'Total interaction Energy (kcal/mol)'
    elif interaction_type == "H-bond":
        interaction_color = "dodgerblue"
        interaction_label = 'H-bond Energy (kcal/mol)'
    elif interaction_type == "lipophilic":
        interaction_color = "darkorange"
        interaction_label = 'Hydrophobic interaction Energy (kcal/mol)'
    elif interaction_type == "Salt bridge":
        interaction_color = "seagreen"
        interaction_label = 'Salt Bridge interaction Energy (kcal/mol)'
    else:
        raise Exception("ERROR: Martin introduced a new interaction type")

    if interaction_partner == "ligand":
            seq_range = range(1, 113)
            plot_range = range(1, 113, 10)
            resid="ligand_resid"
    elif interaction_partner == "receptor":
            seq_range = range(445, 684)
            plot_range = range(445, 684, 20)
            resid="receptor_resid"
    else:
            raise Exception("Error: Interaction partner not found.")

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


    # 3. Data visualtisation
    figure_files = [args.ligand_interaction, args.receptor_interaction]
    for fig_file, interaction_partner in zip(figure_files, ["ligand", "receptor"]):
        plt.figure(figsize=(15, 30))
        for i,interaction_type in enumerate(["total", 'H-bond', 'lipophilic', "Salt bridge"]): # , 'Marked as Salt-Bridge'
            final, energy_max = interaction_data_aggregation(interaction_partner, interaction_type)
            plt.subplot(4, 1, i+1)
            plot_interactions(final, energy_max)

        plt.tight_layout()
        plt.savefig(fig_file)
        plt.close()