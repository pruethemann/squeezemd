#!/usr/bin/env python
"""
    This script generates a pymol script which labels all required mutations in BD001
    python3 bin/10_InteractionSurface.py --interactions output/demo/results/martin/interactions.csv   --seed 695 --output output/demo/results/interactionSurface   --complex C1s_BD001 C1s_BD001 --mutation WT Y117E_Y119E_Y121E --frames output/demo/C1s_BD001/WT/695/MD/frame_end.cif output/demo/C1s_BD001/WT/842/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/695/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/842/MD/frame_end.cif output/demo/C1s_BD001/WT/695/MD/frame_end.cif output/demo/C1s_BD001/WT/842/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/695/MD/frame_end.cif output/demo/C1s_BD001/Y117E_Y119E_Y121E/842/MD/frame_end.cif

"""
import pandas as pd
import argparse
import MDAnalysis as mda
import openmm.app as app
from Helper import remap_MDAnalysis
import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use('ggplot')
sns.set_style('ticks')

def create_pml(ligand_resids, rec_resids, input_pdb, output_pdb, output, target):

    # Adapted free energy caluclation file
    # TODO: find config
    with open('config/pymol.pml', 'r') as f:
        content = f.read()

        content = content.replace("INPUT", input_pdb)
        content = content.replace("LIGAND_RESIDS", ligand_resids)
        content = content.replace("REC_RESIDS", rec_resids)
        content = content.replace("OUTPUT", output_pdb)
        content = content.replace("TARGET", target)

    # Save Tleap conataing all file paths
    f = open(output, "w")
    f.write(content)
    f.close()



def set_bfactors(pdb, ligand_resids, rec_resids, output):

    # Import pdb
    # Import trajectory
    topo = app.PDBxFile(pdb)      # Transform cif to MDAnalysis topology
    u = mda.Universe(topo)
    u = remap_MDAnalysis(u,topo)

    # probably not necessary
    u.add_TopologyAttr('tempfactors')

    protein = u.select_atoms("protein")

    for resid in ligand_resids.resid:
        selected_resid = u.select_atoms(f"resid {resid} and segid I")
        selected_resid.tempfactors = float(ligand_resids[(ligand_resids.resid == resid)]['energy'])

    for resid in rec_resids.resid:
        selected_resid = u.select_atoms(f"resid {resid} and not segid I")
        selected_resid.tempfactors = float(rec_resids[(rec_resids.resid == resid)]['energy'])

    protein.write(output)

def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--interactions', required=False, default='/home/pixelline/ownCloud/Institution/code/squeezeMD_run/V4/output/demo1_5ns/results/martin/interactions.csv')
    parser.add_argument('--frames', nargs='+',required=True)     # Simulation overview

    # params
    parser.add_argument('--seed', required=False, type=int, default=0)
    parser.add_argument('--complex', nargs='+', help="", required=True,
                        default=['', ''])
    parser.add_argument('--mutation', nargs='+', help="", required=True,
                        default=['', ''])

    # Output
    parser.add_argument('--output', required=False, help='output folder', default='output/tmp')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    # TODO: Extract from arg.complex
    targets = ['C1s']
    seed = str(args.seed)

    # Go only for one representatitve final position
    frames = [f for f in args.frames if seed in f]

    # Import interaction data
    interactions = pd.read_csv(args.interactions)
    interactions = interactions[interactions.interaction=='inter']

    # Aggregate TODO extend to multiple targets
    #interactions_agg = interactions[['protein', 'target','mutation', 'resid', 'seed', 'chainID', 'energy']].groupby(['target', 'chainID', 'resid']).mean()
    #interactions_agg = interactions[['protein', 'target', 'mutation', 'resid', 'seed', 'energy']].groupby(['target', 'resid']).mean()
    interactions_agg = interactions[['target', 'resid', 'seed', 'energy']].groupby(['target', 'resid']).mean()

    for target, pdb in zip(targets, frames):
        data_ligand = interactions_agg.loc[(target)].reset_index()
        ligand_resids = ','.join(map(str, data_ligand[data_ligand.energy < -2]['resid']))

        data_rec = interactions_agg.loc[target].reset_index()

        # Get all receptor residues with an interaction energy smaller than -2
        rec_resids = ','.join(map(str, data_rec[data_rec.energy < -2]['resid']))

        bfactor_pdbs = f'{args.output}/{target}.interaction.pdb'
        output_pdb = f'{args.output}/{target}.final.pse'

        print(bfactor_pdbs)

        set_bfactors(pdb, data_ligand, data_rec, bfactor_pdbs)

        create_pml(ligand_resids, rec_resids, bfactor_pdbs, output_pdb, f'{args.output}/{target}.pml', target)
