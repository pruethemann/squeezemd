#!/usr/bin/env python
import MDAnalysis
from MDAnalysis.analysis import contacts
import os
import MDAnalysis.transformations as trans
import argparse
import subprocess
from pathlib import Path
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from simtk.openmm.app import pdbxfile



def execute(command):
    """
    Executes commands in console
    :param command:
    :return:
    """

    output_text = subprocess.check_output(command, shell=True)
    return output_text

def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)


def save_frame(u, out_frame, seleAlgebra='not name EPW', frameid=0):
    print(f"Save Frame {frameid}", out_frame)

    selection = u.select_atoms(seleAlgebra)

    with mda.Writer(out_frame, selection.n_atoms) as W:

        for frame_id, ts in enumerate(u.trajectory[frameid:frameid+1]):
            selection = u.select_atoms(seleAlgebra)
            W.write(selection)

def chain2resid_df(file_csv):
    # Find start and end of chain A

    renum = pd.read_csv(file_csv,
                        delim_whitespace=True,
                        names=['resname', 'chainID', 'resid', 'resname amber', 'resid amber'])

    resids = renum[['chainID', 'resid']]
    x = resids.groupby(['chainID']).agg(['min', 'max'])
    renum.set_index(['chainID', 'resid'], inplace=True)

    resid_min = x[('resid', 'min')]
    resid_max = x[('resid', 'max')]

    resid = pd.concat([resid_min, resid_max])

    x = renum.loc[(slice(None),resid), :]
    x = x.reset_index()
    x = x.drop_duplicates('chainID')
    return x


def chain2resid(file_csv):
    # Find start and end of chain A

    renum = pd.read_csv(file_csv,
                        delim_whitespace=True,
                        names=['resname', 'chainID', 'resid', 'resname amber', 'resid amber'])

    chains = {}

    for chainID in renum.chainID.unique():
        chain = renum[renum.chainID == chainID]

        resID_min = chain['resid amber'].min()
        resID_max = chain['resid amber'].max()

        chains[chainID] = (resID_min, resID_max)

    return chains


def export_frames(u, frame_pdb, selection, start_frame:int, frame_count:int):
    steps = 1 # (len(u.trajectory) - start_frame) // frame_count
    print("start frame: ", start_frame, "end frame", len(u.trajectory), steps)

    # aweful counter
    counter = 0
    for i in range(start_frame, len(u.trajectory), steps):

        save_frame(u, frame_pdb, seleAlgebra=selection, frameid=i)
        counter += 1

def interaction_analyzer(frame_pdb, ligand_csv, receptor_csv):
    """
    Execute Martin's interaction analyzer
    :param frames_nr:
    :param args:
    :return:
    """

    print("Analizer ")
    # Analyze interactions of ligand to receptor

    #command = f'./7_interaction_csv.x {frame_pdb} {lig_resname} {lig_resid} > {ligand_csv}'
    command = f'./7_interaction_csv.x {frame_pdb} ALA 329 > {ligand_csv}'
    print(command)
    execute(command)

    # Analyze interaction of receptor to ligand
    #command = f'./7_interaction_csv.x {frame_pdb} {rec_resname} {rec_resid} > {receptor_csv}'

    command = f'./7_interaction_csv.x {frame_pdb} ARG 139 > {receptor_csv}'
    execute(command)
    print(command)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False, help='Topo file', default='output/WT-simple/topo.pdb')
    parser.add_argument('--traj', required=False, help='Trajectory', default='output/WT-simple/traj_150.dcd')
    parser.add_argument('--mapping', required=False, help='working directory')
    parser.add_argument('--frame_id', required=False, help='working directory')
    parser.add_argument('--dir', required=False, help='working directory')

    # Output
    parser.add_argument('--frame_pdb', required=False, help='working directory')
    parser.add_argument('--ligand_csv', required=False, help='working directory')
    parser.add_argument('--receptor_csv', required=False, help='working directory')

    args = parser.parse_args()

    # TODO: Extract target from path
    target = os.path.basename(args.traj).split('_')[0]
    target = 'C1s' # TODO
    number_frames = int(args.frame_id)

    # Import Trajectory
    topo = pdbxfile.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)
    #u = mda.Universe(args.traj, in_memory=False)
    traj_length = len(u.trajectory)

    # Export centered frames for Martin Analyzer
    print("Start extracting last frames")
    selection = f'(protein or byres around 8.0 protein) and not name EPW'          # protein and surface
    #export_frames(u, args.frame_pdb, selection, start_frame=traj_length-frame_id, frame_count=1)

    steps = 1
    frame_id = 0
    for i in range(len(u.trajectory) - number_frames, len(u.trajectory), steps):
        frame_path = os.path.join(args.dir, f'frame_{frame_id}.pdb')
        lig_csv = os.path.join(args.dir, 'lig', f'{frame_id}.csv')
        rec_csv = os.path.join(args.dir, 'rec', f'{frame_id}.csv')

        save_frame(u, frame_path, seleAlgebra=selection, frameid=frame_id)
        interaction_analyzer(frame_path, lig_csv, rec_csv)

        frame_id += 1


"""
    # Get chain information since amber deleted all chain IDs
    chains = chain2resid(args.mapping) # TODo
    chains_df = chain2resid_df(args.mapping) # TODo

    print(chains_df)


    chains_df.set_index('chainID', inplace=True)

    lig_resname = chains_df.loc['I', 'resname']
    lig_resid = chains_df.loc['I', 'resid']

    rec_resname = chains_df.loc['A', 'resname']
    rec_resid = chains_df.loc['A', 'resid']

"""
