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

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')
sns.set_style('ticks')

def create_pml(ligand_resids, rec_resids, input_pdb, output_pdb, output, target):

    # Adapted free energy caluclation file

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
    u = mda.Universe(pdb)

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
    parser.add_argument('--interactions', required=False, default='output/19-04-23_proteases/results/martin/interactions.csv')
    parser.add_argument('--simulations', required=False, default='output/19-04-23_proteases/simulations.csv')     # Simulation overview

    # Output
    parser.add_argument('--output', required=False, help='output folder', default='output/tmp')
    parser.add_argument('--checkpoint', required=False, help='output folder')

    args = parser.parse_args()

    Path(args.output).mkdir(parents=True, exist_ok=True)

    #sims = pd.read_csv(args.simulations, index_col='sim_id', converters={"mutations": ast.literal_eval})

    # Add path to last frame
    # TODO replace to last frame
    #sims['frame_end'] = sims['path'] + '/frames/frame_9.pdb'

    # Reduce the amount of end frames to one per replicates
    #end_frames = sims.drop_duplicates(['name'])

    targets = ['C1s']
    frames = ['/home/pixelline/ownCloud/Institution/code/squeezeMD_run/V4/output/demo/C1s_BD001/WT/695/frames/frame_9.pdb',
              '/home/pixelline/ownCloud/Institution/code/squeezeMD_run/V4/output/demo/C1s_BD001/Y117E_Y119E_Y121E/695/frames/frame_9.pdb']

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
        #data_rec = data_rec[data_rec.chainID != 'I']
        rec_resids = ','.join(map(str, data_rec[data_rec.energy < -2]['resid']))

        print(target)
        print(data_ligand)
        #print(rec_resids)
        print(args.output)
        bfactor_pdbs = f'{args.output}/{target}.interaction.pdb'
        output_pdb = f'{args.output}/{target}.final.pse'

        set_bfactors(pdb, data_ligand, data_rec, bfactor_pdbs)

        create_pml(ligand_resids, rec_resids, bfactor_pdbs, output_pdb, f'{args.output}/{target}.pml', target)


    # Calculate total Inter - Interaction energy
    #total_energies = interactions[['target','chainID', 'seed', 'frame', 'energy']].groupby(['target', 'chainID', 'seed', 'frame']).agg(['sum'])
    total_energies = interactions[['target', 'seed', 'frame', 'energy']].groupby(['target', 'seed', 'frame']).agg(['sum'])
    total_energies = total_energies.reset_index()
    #total_energies = total_energies[total_energies.chainID == 'I']
    #print(total_energies)

    #total_energies = interactions[interactions.chainID == 'I']

    sns.barplot(data=total_energies,
                x='target',
                y=('energy', 'sum'))

    plt.show()

    # TODO: uncomment
    Path(args.checkpoint).touch()
