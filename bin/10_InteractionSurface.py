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

def create_pml_script(ligand_resids, receptor_resids, input_pdb, output_file, pymol_script, target):
    """
    Generates a PyMOL script from a template, substituting placeholders with actual data.

    Parameters:
    - ligand_resids: Comma-separated string of ligand residue IDs.
    - receptor_resids: Comma-separated string of receptor residue IDs.
    - input_pdb: Path to the input PDB file.
    - output_file: Path where the generated PyMOL script will be saved.
    - target: Name of the target protein.
    """
    #TODO find this location
    with open('/home/peter/tools/miniconda3/envs/squeezeMD/pymol_template.pml', 'r') as template_file:
        content = template_file.read().format(input_pdb=input_pdb,
                                               ligand_resids=ligand_resids,
                                               receptor_resids=receptor_resids,
                                               output=output_file,
                                               target=target)
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

 # Import pdb
    # Import trajectory
    topo = app.PDBxFile(pdb_path)      # Transform cif to MDAnalysis topology
    u = mda.Universe(topo)
    u = remap_MDAnalysis(u,topo)

    # probably not necessary
    u.add_TopologyAttr('tempfactors')


    for _,row in ligand_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row.resid)} and segid I")
        selected_resid.tempfactors = row.energy

    for _,row in receptor_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row.resid)} and not segid I")
        selected_resid.tempfactors = row.energy

    # Save pdb of protein only
    protein = u.select_atoms("protein")
    print(protein[protein.tempfactors<0])

    protein.write(interaction_pdb)

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
    A namespace object containing the arguments.
    """
    parser = argparse.ArgumentParser(description='Generate PyMOL script for BD001 mutation labeling and interaction surface calculation.')
    parser.add_argument('--interactions', required=True, help='Path to the interactions CSV file.')
    parser.add_argument('--seed', type=int, required=True, help='Seed number for selecting representative frames.')
    parser.add_argument('--output', required=True, help='Output directory for generated files.')
    parser.add_argument('--mutation', nargs='+', required=True, help='Mutation identifiers (e.g., WT, Y117E_Y119E_Y121E).')
    parser.add_argument('--frames', nargs='+', required=True, help='Paths to frame files.')
    parser.add_argument('--receptors', nargs='+', required=True, help='List of all Receptors (e.g. C1s, MASP2, FXa')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    # Go only for one representatitve final position. Choose only the first seed number
    frames = [f for f in args.frames if str(args.seed) in f]

    # TODO: for the moment only WT frames
    frames = [f for f in args.frames if 'WT' in f]

    # Import interaction data
    interactions = pd.read_parquet(args.interactions)
    interactions = interactions[interactions.interaction=='inter']

    # TODO Currently only WT
    interactions = interactions[interactions.mutation == 'WT']

    # Aggregate data for the particular target and residue
    var_names = ['protein', 'target', 'resid', 'energy']        # Relevant var names
    interactions_agg = interactions[var_names].groupby(['protein', 'target', 'resid']).mean()

    ENERGY_THRESHOLD = -2

    for receptor_name, last_frame in zip(args.receptors, frames):
        # Extract ligand and receptor interaction data
        data_ligand = interactions_agg.loc[('ligand',receptor_name)].reset_index()
        data_receptor = interactions_agg.loc['receptor', receptor_name].reset_index()

        # Get all receptor/ligand residues with an interaction energy smaller than -2 and join as string
        ligand_resids = ','.join(map(str, data_ligand[data_ligand.energy < ENERGY_THRESHOLD]['resid']))
        receptor_resids = ','.join(map(str, data_receptor[data_receptor.energy < ENERGY_THRESHOLD]['resid']))

        # TODO: figure out why there is a random wrong resid in each dataset. Probably too many atoms in pdb
        #data_ligand.to_csv('data_ligand.csv')
        #data_receptor.to_csv('data_receptor.csv')

        # Define output paths. TODO Improve
        interactions_pdb = os.path.join(args.output, f'{receptor_name}.interaction.pdb')
        pymol_out = os.path.join(args.output, f'{receptor_name}.final.pse')
        pymol_script = os.path.join(args.output, f'{receptor_name}.pml')

        # Set the interaction intensities
        set_residue_interaction_intensity(last_frame, data_ligand, data_receptor, interactions_pdb)

        # create a custom pymol script to visualize the relevant interactions
        create_pml_script(ligand_resids, receptor_resids, interactions_pdb, pymol_out, pymol_script, receptor_name)

