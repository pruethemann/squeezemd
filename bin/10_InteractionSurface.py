#!/usr/bin/env python
"""
    This script generates a pymol script which labels all required mutations in BD001
"""
from pathlib import Path
import pandas as pd
import argparse
import numpy as np
from glob import glob
import MDAnalysis as mda
import ast
import openmm.app as app
from Helper import remap_MDAnylsis

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')
sns.set_style('ticks')

def create_pml(ligand_resids, rec_resids, input_pdb, output_pdb, output, target):

    # Adapted free energy caluclation file

    with open('../config/pymol.pml', 'r') as f:
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
    u = remap_MDAnylsis(u,topo)

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


if __name__ == '__main__':

    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--interactions', required=False, default='/home/pixelline/ownCloud/Institution/code/squeezeMD_run/V4/output/demo1_5ns/results/martin/interactions.csv')
    parser.add_argument('--final_frame', required=False, default='output/19-04-23_proteases/simulations.csv')     # Simulation overview

    # Output
    parser.add_argument('--output', required=False, help='output folder', default='output/tmp')

    args = parser.parse_args()

    Path(args.output).mkdir(parents=True, exist_ok=True)

    # TODO: Do over arguments
    targets = ['C1s']
    frames = ['/home/pixelline/ownCloud/Institution/code/squeezeMD_run/V4/output/demo/C1s_BD001/WT/695/MD/frame_end.cif', ' nur Blödsinn']
    seed = str(695)

    # Go only for one representatitve final position
    frames = [f for f in frames if seed in f]

    # Import interaction data
    interactions = pd.read_csv(args.interactions)
    interactions = interactions[interactions.interaction=='inter']

    # Aggregate
    #interactions_agg = interactions[['protein', 'target','mutation', 'resid', 'seed', 'chainID', 'energy']].groupby(['target', 'chainID', 'resid']).mean()
    print(interactions)
    #interactions_agg = interactions[['protein', 'target', 'mutation', 'resid', 'seed', 'energy']].groupby(['target', 'resid']).mean()
    interactions_agg = interactions[['target', 'resid', 'seed', 'energy']].groupby(['target', 'resid']).mean()

    for target, pdb in zip(targets, frames):

        print(interactions)
        data_ligand = interactions_agg.loc[(target)].reset_index()
        ligand_resids = ','.join(map(str, data_ligand[data_ligand.energy < -2]['resid']))

        data_rec = interactions_agg.loc[target].reset_index()

        # Get all receptor residues with an interaction energy smaller than -2
        rec_resids = ','.join(map(str, data_rec[data_rec.energy < -2]['resid']))

        bfactor_pdbs = f'{args.output}/{target}.interaction.pdb'
        output_pdb = f'{args.output}/{target}.final.pse'

        set_bfactors(pdb, data_ligand, data_rec, bfactor_pdbs)

        create_pml(ligand_resids, rec_resids, bfactor_pdbs, output_pdb, f'{args.output}/{target}.pml', target)
