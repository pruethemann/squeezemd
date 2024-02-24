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

def create_pml_script(ligand_resids, receptor_resids, pdb, output_file, pymol_script):
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
    with open('/home/pixelline/peter/miniconda3/envs/squeezeMD/pymol_template.pml', 'r') as template_file:
        content = template_file.read().format(input_pdb=pdb,
                                               ligand_resids=ligand_resids,
                                               receptor_resids=receptor_resids,
                                               output=output_file
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

    # TODO: fix the issue that the ligand becomes chain B after export to topo_center.pdb
    for _,row in ligand_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row.resid)} and chainID B")
        selected_resid.tempfactors = row.energy

    for _,row in receptor_resids.iterrows():
        selected_resid = u.select_atoms(f"resid {int(row.resid)} and not chainID B")
        selected_resid.tempfactors = row.energy

    # Save pdb of protein only
    protein = u.select_atoms("protein")
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
    parser.add_argument('--mutations', nargs='+', required=True, help='Mutation identifiers (e.g., WT, Y117E_Y119E_Y121E).')
    parser.add_argument('--frames', nargs='+', required=False, help='Paths to frame files.')
    parser.add_argument('--receptors', nargs='+', required=False, help='List of all Receptors (e.g. C1s, MASP2, FXa')
    parser.add_argument('--complexes', nargs='+', required=True, help='')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    ENERGY_THRESHOLD = -2
    # Import interaction data
    interactions = pd.read_parquet(args.interactions)
    interactions.set_index(['interaction', 'name', 'mutation'], inplace=True)
    # Sort index to improve performance
    interactions.sort_index(inplace=True)
    var_names = ['protein', 'target', 'resid', 'energy']  # Relevant var names


    for complex in args.complexes:
        for mutation in args.mutations:
            pdb = os.path.join(complex, mutation, str(args.seed), 'MD', 'topo_center.pdb')

            interactions_filtered = interactions.loc[('inter',complex, mutation)]

            # Aggregate data for the particular target and residue
            interactions_agg = interactions_filtered.groupby(['protein', 'name', 'mutation', 'resid']).mean(numeric_only=True)

            # Extract ligand and receptor interaction data
            data_ligand = interactions_agg.loc[('ligand',complex)].reset_index()
            data_receptor = interactions_agg.loc[('receptor', complex)].reset_index()

            # Get all receptor/ligand residues with an interaction energy smaller than -2 and join as string
            ligand_resids = ','.join(map(str, data_ligand[data_ligand.energy < ENERGY_THRESHOLD]['resid']))
            receptor_resids = ','.join(map(str, data_receptor[data_receptor.energy < ENERGY_THRESHOLD]['resid']))

            DEBUG = True
            if DEBUG:
                data_ligand.to_csv('data_ligand.csv')
                data_receptor.to_csv('data_receptor.csv')

            # Define output paths. TODO Improve
            dir = os.path.join('results', 'interactionSurface')
            interaction_pdb = os.path.join(dir, f'{complex}.{mutation}.interaction.pdb')
            pymol_out = os.path.join(dir, f'{complex}.{mutation}.final.pse')
            pymol_script = os.path.join(dir, f'{complex}.{mutation}.pml')

            # Set the interaction intensities
            set_residue_interaction_intensity(pdb, data_ligand, data_receptor, interaction_pdb)

            # create a custom pymol script to visualize the relevant interactions
            create_pml_script(ligand_resids, receptor_resids, interaction_pdb, pymol_out, pymol_script)
