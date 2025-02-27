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
import argparse

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


def parse_lipophilic(parts):
        
        interaction_info = parts[0].split()
        donor_acceptor = parts[1].strip().split()

        receptor_resname = donor_acceptor[1]
        receptor_resid = int(donor_acceptor[2])

        ligand_resname = donor_acceptor[-2]
        ligand_resid = int(donor_acceptor[-1])

        # TODO: solve the differentiation issue of ligand vs receptor
        if ligand_resid > 122:
             (ligand_resid, receptor_resid) = (receptor_resid, ligand_resid)
             (ligand_resname, receptor_resname) = (receptor_resname,ligand_resname)

        distance = float(interaction_info[2].split("=")[1])
        energy = float(interaction_info[3].split("=")[1])

        interaction = {
            "Interaction Type": 'lipophilic',
            "Distance (r)": distance,
            "Energy (e)": energy,
            'receptor_resname' : receptor_resname,
            'receptor_resid' : receptor_resid,
            'ligand_resname' : ligand_resname,
            'ligand_resid' : ligand_resid,
        }

        return interaction


def parse_hbonds(parts):

        interaction_info = parts[0].split()
        donor_acceptor = parts[1].strip().split()

        ligand_resname = donor_acceptor[1]
        ligand_resid = int(donor_acceptor[2])

        receptor_resname = donor_acceptor[-2]
        receptor_resid = int(donor_acceptor[-1])

        # TODO: solve the differentiation issue of ligand vs receptor
        if ligand_resid > 122:
             (ligand_resid, receptor_resid) = (receptor_resid, ligand_resid)
             (ligand_resname, receptor_resname) = (receptor_resname,ligand_resname)

        distance = float(interaction_info[2].split("=")[1])
        angle = float(interaction_info[3].split("=")[1])
        energy = float(interaction_info[4].split("=")[1])

        # Include salt bridge data
        marked =  "marked as salt-bridge" in parts[2]

        interaction = {
            "Interaction Type": 'H-bond',
            "Distance (r)": distance,
            "Angle (a)": angle,
            "Energy (e)": energy,
            'receptor_resname' : receptor_resname,
            'receptor_resid' : receptor_resid,
            'ligand_resname' : ligand_resname,
            'ligand_resid' : ligand_resid,
            "Marked as Salt-Bridge": marked
        }

        return interaction

# Parse the input data into a pandas DataFrame
def parse(posco_output):
    data = []

    with open(posco_output, 'r') as file:
        for line in file:

            if line.startswith("Lipo_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_lipophilic(parts)
                data.append(interaction)

            if line.startswith("HB_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_hbonds(parts)
                data.append(interaction)

    data = pd.DataFrame(data)
    return extract_metadata(posco_output, data)



def main(args):

    results = []
    for posco_output in args.input:
         posco_result = parse(posco_output)
         results.append(posco_result)

    results = pd.concat(results)

    results.to_parquet(args.output)
    #results.to_csv('posco_interactions.csv')


def extract_metadata(file_path, data):
        # Extract meta data like an idiot
        metadata = file_path.split('/')

        #'{complex}/{mutation}/{seed}/po-sco/{i}.txt'

        complex = metadata[-5]
        mutation = metadata[-4]
        target = complex.split('_')[0]
        ligand = complex.split('_')[1]
        frame_id = int(metadata[-1][:-4])
        seed = int(metadata[-3])

        # Determine metrics lables
        data['name'] = complex
        data['target'] = target
        data['lig'] = ligand
        data['mutation'] = mutation
        data['frame'] = frame_id
        data['seed'] = seed

        return data


def parse_arguments():
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--input', nargs='+', required=False)

    # Output
    parser.add_argument('--output', required=False, default='dev3/interactions.feather')


    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
