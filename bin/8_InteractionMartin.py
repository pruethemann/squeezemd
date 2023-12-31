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

    example: python3 8_InteractionMartin.py --simulations output/demo/simulations.csv --number_frames 5 --interactions debug/interactions.csv
"""

import pandas as pd
from os import path
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import ast
from pathlib import Path

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


    for interaction_csv in interactions:

        # Extract meta data like an idiot
        metadata = interaction_csv.split('/')

        complex = metadata[-6]
        mutation = metadata[-5]
        target = complex.split('_')[0]
        ligand = complex.split('_')[1]
        frame_id = int(metadata[-1][:-4])
        seed = int(metadata[-4])

        print(complex, mutation, target, ligand, frame_id)

        #'output/{job_id}/{complex}/{mutation}/{seed}/frames/lig/{frame_id}.csv'

        # Import data
        try:
            frame_ana = pd.read_csv(interaction_csv, names=['interaction', 'resname', 'resid', 'energy'])
        except FileNotFoundError:
            print("Error with import from: ", interaction_csv)
            continue

        # Determine metrics lables
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

def aggregate_data(stats, output_dir):

    # Perform analysis by aggregation #TODO add SD
    # Total energy: Aggregated over frame, resid, frame
    total_energy = stats.groupby(['protein','interaction', 'target', 'mutation']).mean(numeric_only=True)               # Over different simulation
    total_energy.to_csv(path.join(output_dir, 'total.csv'))

    # Per residue energy: Aggregated over frame, resid # todo group for chainID?
    residues = stats.groupby(['protein', 'interaction', 'target','mutation', 'resid', 'seed']).mean(numeric_only=True)  # Over different simulation

    return (total_energy, residues)

def plot_graph(filtered_data, out_path, x, y, title, hue=None):

    filtered_data.to_csv('filter.csv')

    filtered_data.reset_index(inplace=True)

    sns.barplot(data=filtered_data,
                x=x,
                y=y,
                hue=hue,
                errorbar='sd',
                capsize=0.5)

    plt.xticks(rotation=90)

    plt.xlabel('Residue number')
    plt.ylabel('Binding energy kcal/mol')
    #plt.title(title + str(filtered_data.resid.min()) + str(filtered_data.resid.max()))
    plt.tight_layout()

    plt.savefig(out_path)
    #plt.show()
    plt.close()

def generate_total_graphs(total_energy, protein, output_dir):

    for interaction in ['inter', 'intra']:

        # 1. Total interaction energiesj
        x = total_energy.loc[interaction].reset_index()
        out_path = path.join(output_dir, f'total_energy-{interaction}.svg')

        plot_graph(x, out_path, x='mutation', y='energy', title=interaction)

def generate_residues_graphs(data, targets, mutations, output_dir, protein='ligand'):
    # Aggregation degree of data:

    print("Generate Residue Graph")

    for interaction in ['inter', 'intra']:
        # 2. interaction energies per residues
        for target in targets:
            for mutation in mutations:
                out_path = path.join(output_dir, f'figs/{protein}/{interaction}/residues/residues_{target}_{mutation}.svg')

                # TODO is this debug or userfule?
                #data.loc[(interaction, target)].to_csv(f'{protein}_{interaction} {target}.csv')
                data.to_csv('filtered.csv')
                plot_graph(filtered_data=data.loc[(interaction, target, mutation)],
                           out_path=out_path,
                           x='resid',
                           y='energy',
                           hue=None,
                           title=f'Total energies: Simulation: {target}_{mutation}')

def generate_mutation_graphs(residues, targets, mutations, protein, output_dir):

    energy_diff = {}

    for interaction in ['inter', 'intra']:

        # 2. interaction energies per residues
        for target in targets:

            WT = residues.loc[(interaction, target, 'WT')].reset_index()
            WT.set_index(['interaction', 'target', 'seed', 'frame'], inplace=True)

            for mutation in mutations:
                # Wildtype is not a mutation
                if mutation == 'WT': continue

                mutated = residues.loc[(interaction, target, mutation)].reset_index()
                mutated.set_index(['interaction', 'target', 'resid','seed', 'frame'], inplace=True)

                energy_diff[mutation] = mutated.energy - WT.energy

                plt.figure(figsize=(40, 8))

                filtered_data = energy_diff[mutation].reset_index()

                sns.barplot(data=filtered_data,
                            x='resid',
                            y='energy',
                            hue=None,
                            errorbar='sd',
                            capsize=0.5,
                            errwidth=0.5)

                plt.title(f'Energy difference per residue {target}  vs BD001 ({mutation})')
                plt.xlabel('Residue number')
                plt.ylabel('Binding energy kcal/mol')
                plt.xticks(rotation=90)
                plt.tight_layout()
                plt.savefig(path.join(output_dir, f'figs/{protein}/{interaction}/mutations/dG_{target}_{mutation}.svg'))
                plt.close()

    energy_diff_results = path.join(output_dir, 'energy_diff.csv')

    # Check if mutation analysis has been performed
    if len(energy_diff) == 0: # no mutations
        Path(energy_diff_results).touch()
    else:  # Export analysis
        pd.concat(energy_diff).to_csv(energy_diff_results)


def generate_target_graphs(residues, targets, mutations, protein, output_dir):
    """
    Diff between C1s-BD001-WT and TARGET-BD001-WT

    Todo:
    1. Reset index
    2. set index: interaction, mutation, resid, seed

    :param residues:
    :param targets:
    :param mutations:
    :param protein:
    :param output_dir:
    :return:
    """

    energy_diff = {}

    for interaction in ['inter', 'intra']:

        # Reference
        C1s_ene = residues.loc[(interaction, 'C1s', 'WT')].reset_index()
        C1s_ene.set_index(['interaction', 'mutation', 'resid', 'seed'], inplace=True)

        for target in targets:

            if target == "C1s": continue

            target_ene = residues.loc[(interaction, target, 'WT')].reset_index()
            target_ene.set_index(['interaction', 'mutation', 'resid', 'seed'], inplace=True)

            # Calculate difference between WT targets
            energy_diff[target] = target_ene['energy'] - C1s_ene['energy']


            plt.figure(figsize=(40, 8))

            filtered_data = energy_diff[target].reset_index()

            sns.barplot(data=filtered_data,
                        x='resid',
                        y='energy',
                        hue=None,
                        errorbar='sd',
                        capsize=0.5,
                        errwidth=0.5)

            plt.title(f'Energy difference per residue {target}  vs BD001 ({target})')
            plt.xlabel('Residue number')
            plt.ylabel('Binding energy kcal/mol')
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(path.join(output_dir, f'figs/{protein}/{interaction}/targets/dG_{target}.svg'))
            plt.close()

    energy_diff_results = path.join(output_dir, 'energy_diff.csv')

    # Check if mutation analysis has been performed
    if len(energy_diff) == 0: # no mutations
        Path(energy_diff_results).touch()
    else:  # Export analysis
        pd.concat(energy_diff).to_csv(energy_diff_results)


def create_folders(output_dir):
    """
    Create result folders
    :param output_dir: results folder
    :return:
    """
    for protein in ['ligand', 'receptor']:
        for interaction in ['inter', 'intra']:
            for analysis_type in ['targets', 'residues', 'mutations']:
                Path(path.join(output_dir, 'figs', protein, interaction, analysis_type)).mkdir(parents=True, exist_ok=True)

import os
from glob import glob

def init(args):

    # Creates required output folder in "output/job_id/results
    create_folders(args.output)

    ligs_path = []
    recs_path = []
    for dir in args.dirs:
        lig_dir = os.path.join(dir, 'lig/')
        rec_dir = os.path.join(dir, 'rec/')
        ligs_path.extend(glob(lig_dir + '*.csv'))
        recs_path.extend(glob(rec_dir+ '*.csv'))

    print(ligs_path)

    # 3. Create result table including: inter and intramolecular interactions / all simulations / all seeds / ligand and receptor perspective for ligand AND receptor
    data_ligand = generate_data(ligs_path, protein='ligand')
    data_receptor = generate_data(recs_path,  protein='receptor')
    data = pd.concat([data_ligand, data_receptor])          # Join ligand and receptor data

    print(data)

    targets = data_ligand.target.unique()
    mutations = data_ligand.mutation.unique()

    # Aggregate data over energy and per residue
    total_energy, residues = aggregate_data(data, args.output)

    #total_energy.to_csv('total_energy.csv')
    #residues.to_csv('residues.csv')

    # Export data: Raw (data), per residue per MD data
    data.to_csv(args.interactions)

    residues.to_csv(path.join(args.output, 'residues_agg.csv'))

    # TODO: Debug import instead of concat data
    #data = pd.read_csv(args.interactions)
    #total_energy, residues = aggregate_data(data, args.output)


    data.set_index(['protein','interaction', 'target', 'mutation'], inplace=True)
    data.sort_index(inplace=True)
    data.to_csv('delete.csv')

    # 4. Generate residue graphs
    for protein in ['ligand', 'receptor']:

        # Generate total energy graphs. TODO: sep receptor / ligand
        generate_total_graphs(total_energy.loc[protein], protein, args.output)

        # Generate per residue graphs
        generate_residues_graphs(data.loc[protein], targets, mutations, args.output, protein)

        # 5. Calculate differences between mutations
        generate_mutation_graphs(data.loc[protein], targets, mutations, protein, args.output)

        # 6. Calculate differences between targets
        #generate_target_graphs(data.loc[protein], targets, mutations, protein, args.output)



if __name__ == '__main__':
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--dirs', nargs='+', required=False)

    # Output
    parser.add_argument('--interactions', required=False, default='output/1_FactorXa_50ns/results/martin/interactions.csv')
    parser.add_argument('--output', required=False, default='output/1_FactorXa_50ns/martin')

    args = parser.parse_args()
    job_id = args.output.split('/')[-3]

    init(args)
