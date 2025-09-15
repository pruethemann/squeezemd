#!/usr/bin/env python
"""
This script generates a PyMOL script for labeling mutations in the BD001 complex and calculating interaction surfaces.
It reads interaction data, selects representative frames based on a specified seed, and processes this information
to visualize specific mutations and their interaction energies within the complex.

Terminology:
 - Target: Receptor # TODO rename everywhere
 -

 Example:
    python3 bin/10_InteractionSurface.py --output output/demo/results/interactionSurface   --interactions output/demo/results/martin/interactions.csv     --seed 695   --mutation WT Y117E_Y119E_Y121E --frames output/demo/C1s_BD001/WT/695/MD/frame_end.cif output/demo/C1s_BD001/WT/842/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/695/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/842/MD/frame_end.cif  --receptors C1s

    # TODO extend to multiple targets
    #interactions_agg = interactions[['protein', 'target','mutation', 'resid', 'seed', 'chainID', 'energy']].groupby(['target', 'chainID', 'resid']).mean()
    #interactions_agg = interactions[['protein', 'target', 'mutation', 'resid', 'seed', 'energy']].groupby(['target', 'resid']).mean()

    Data variable description:
    Group by:
        name: same as complex
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
import argparse
import MDAnalysis as mda
import openmm.app as app
from Helper import remap_MDAnalysis
import seaborn as sns
import matplotlib.pyplot as plt
import os

plt.style.use('ggplot')
sns.set_style('ticks')

def create_pml_script(ligand_resids, receptor_resids, pdb, output_file, pymol_script, output_png):
    """
    Generates a PyMOL script from a template, substituting placeholders with actual data.

    Parameters:
    - ligand_resids: Comma-separated string of ligand residue IDs.
    - receptor_resids: Comma-separated string of receptor residue IDs.
    - input_pdb: Path to the input PDB file.
    - output_file: Path where the generated PyMOL script will be saved.
    - target: Name of the target protein.
    """

    # Check if the script is running in a Conda environment and extract path to pymol_template in env
    if 'CONDA_PREFIX' in os.environ:
        conda_env = os.environ['CONDA_PREFIX']
        pymol_template = os.path.join(conda_env, 'bin', 'pymol_template.pml')
    else:
        raise Exception("This script is not running in a Conda environment.")

    with open(pymol_template, 'r') as template_file:
        content = template_file.read().format(input_pdb=pdb,
                                               ligand_resids=ligand_resids,
                                               receptor_resids=receptor_resids,
                                               output=output_file,
                                               output_png=output_png
                                              )
    with open(pymol_script, 'w') as output_pml:
        output_pml.write(content)

def set_residue_interaction_intensity(pdb_path, ligand_resids, receptor_resids, interaction_pdb):
    """
    Sets the interaction intensity per residue and ligand and receptor in a pdb file of the last frame
    of the molecular dynamics simulation. Functions saves the interaction intensities in b factor column.

    Parameters:
    - pdb_path: Path to the PDB file.
    - ligand_resids: List of ligand residue IDs.
    - receptor_resids: List of receptor residue IDs.
    - output_path: Path to save the modified PDB file.
    """
    # This function's implementation will depend on specific requirements for adjusting B-factors.

    # Import trajectory
    u = mda.Universe(pdb_path)

    # probably not necessary
    u.add_TopologyAttr('tempfactors')


    # Select all ligand resids interacting with ligand
    for _,row in ligand_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row['ligand_resid'])} and segid A")
        selected_resid.tempfactors = row['Energy (e)']

    # Select all receptor resids interacting with receptor
    for _,row in receptor_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row['receptor_resid'])} and not segid A")
        selected_resid.tempfactors = row['Energy (e)']

    # Save pdb of protein only
    protein = u.select_atoms("protein")
    protein.write(interaction_pdb)

def data_aggregation (data):
    """
    please provide a pandas dataframe, such as .parquet read by pd.read_parquet
      """

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

    df_interaction = data

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

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
    A namespace object containing the arguments.
    """
    parser = argparse.ArgumentParser(description='Generate PyMOL script for BD001 mutation labeling and interaction surface calculation.')
    parser.add_argument('--interactions', required=True, help='Path to the interactions CSV file.')
    parser.add_argument('--seed', type=int, required=True, help='Seed number for selecting representative frames.')
    parser.add_argument('--mutation', required=True, help='Mutation identifiers (e.g., WT, Y117E_Y119E_Y121E).')
    #parser.add_argument('--frames', nargs='+', required=False, help='Paths to frame files.')
    #parser.add_argument('--receptors', required=False, help='List of all Receptors (e.g. C1s, MASP2, FXa') # TODO: will be used for the future
    parser.add_argument('--complex', required=True, help='')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    # TODO Load aggregated dataframe from 5.5 / 5.4

    # Import interaction data
    interactions = pd.read_parquet(args.interactions)

    # TODO perform a separate water analysis
    n_frames = len(interactions.frame.unique())
    n_seeds = len(interactions.seed.unique())
        
    interactions = interactions[(interactions['receptor_resname'] != 'HOH') & (interactions['ligand_resname'] != 'HOH')]

    interactions.set_index(['name', 'mutation'], inplace=True)
    # Sort index to improve performance
    interactions.sort_index(inplace=True)

    (mutation, complex) = (args.mutation, args.complex)

    pdb = os.path.join(complex, mutation, str(args.seed), 'MD', 'topo_center.pdb')
    
    interactions_filtered = interactions.loc[(complex, mutation)]

    # Extract ligand and receptor interaction data
    #data_ligand = interactions_filtered.groupby(['name', 'mutation', 'ligand_resid']).mean(numeric_only=True).reset_index()
    #data_receptor = interactions_filtered.groupby(['name', 'mutation', 'receptor_resid']).mean(numeric_only=True).reset_index()

    ## data aggregate
    data_ligand = interactions.groupby(['name', 'mutation', 'ligand_resid'])['Energy (e)'].sum().reset_index()
    data_ligand["Energy (e)"] = data_ligand["Energy (e)"].div(n_frames * n_seeds)

    data_receptor = interactions.groupby(['name', 'mutation', 'receptor_resid'])['Energy (e)'].sum().reset_index()
    data_receptor["Energy (e)"] = data_receptor["Energy (e)"].div(n_frames * n_seeds)
    
    # Get all receptor/ligand residues with an interaction energy smaller than -2 and join as string
    # only consider really strong interactions
    ENERGY_THRESHOLD = -0.8
    data_ligand = data_ligand[data_ligand['Energy (e)'] < ENERGY_THRESHOLD]#['ligand_resid']

    # create a string in pymol
    ligand_resids = ','.join(map(str, data_ligand))
    receptor_resids = ','.join(map(str, data_receptor[data_receptor['Energy (e)']  < ENERGY_THRESHOLD]['receptor_resid']))

    # Define output paths. TODO Improve
    dir = os.path.join('results', 'interactionSurface')
    interaction_pdb = os.path.join(dir, f'{complex}.{mutation}.interaction.pdb')
    pymol_out = os.path.join(dir, f'{complex}.{mutation}.final.pse')
    pymol_script = os.path.join(dir, f'{complex}.{mutation}.pml')
    output_png = os.path.join(dir, f'{complex}.{mutation}.png')

    # Set the interaction intensities
    set_residue_interaction_intensity(pdb, data_ligand, data_receptor, interaction_pdb)

    # create a custom pymol script to visualize the relevant interactions
    create_pml_script(ligand_resids, receptor_resids, interaction_pdb, pymol_out, pymol_script, output_png)
