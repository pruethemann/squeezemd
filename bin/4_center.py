#!/usr/bin/env python
import MDAnalysis
from MDAnalysis.analysis import contacts
import MDAnalysis.transformations as trans
from tqdm import tqdm
import MDAnalysis as mda
import argparse
import pandas as pd
import subprocess
import openmm.app as app

# python3 4_center.py --topo MD_pdb/C1s_BD001-tleap-water.pdb --traj MD_pdb/trajectory.dcd --mapping output/demo/C1s_BD001/WT/amber/amber_renum.txt --traj_center center.dcd --topo_center center.pdb

def chain2resid(file_csv):
    # Find start and end of chain A

    renum = pd.read_csv(file_csv,
                        delim_whitespace=True,
                        names=['resname', 'chainID', 'resid', 'resname amber', 'resid amber'])

    del renum['resname']
    del renum['resname amber']

    chain_min = renum.groupby('chainID').min().rename(columns={'resid amber': 'amber_start', 'resid': 'start'})
    chain_max = renum.groupby('chainID').max().rename(columns={'resid amber': 'amber_end', 'resid': 'end'})

    chains = pd.concat([chain_min, chain_max], axis=1)

    # TODO: Add one resname of every chain for start resid

    return chains

def align_and_center(u):
    """
    Function centers protein in the middle of the box and wraps water around it.
    All frames are aligned to the first frame.
    :param u: Trajector as universe
    :return: transformed Universe
    """

    # now we have to create atomgroups to represent the protein, non-protein, and c-alpha sections of the system
    protein = u.select_atoms('protein')
    not_protein = u.select_atoms('not protein')

    # next we create a chain of transformations which will be applied each time a frame is loaded from the Universe
    # Note: the transformations are applied linearly, so the order of the transformations matters!

    transforms = [trans.unwrap(protein),
                  trans.center_in_box(protein, wrap=True),
                  trans.wrap(not_protein)
                  ]

    u.trajectory.add_transformations(*transforms)

    return u


def renumber_resids(u, chainID:str, chain_start, chain_end, shift_factor):
    """
    This function restores the original resid numbering destroyed by pdb4amber
    IMPORTANT: tleap names all chains A. Since I am selecting the chain based on the
    start and end resids and directly shift the numbering. It's important to shift
    the orginal chain A only at the end.
    :return:
    """

    # Assign chain ID to amber resid numbering
    chain = u.select_atoms(f"resid {chain_start} to {chain_end}")
    chain.atoms.chainIDs = chainID

    # Shift resid numbering to old numbering
    chain.residues.resids += shift_factor

def renumber_resids_original(u, chainID:str, chain_start, chain_end, shift_factor):
    """
    This function restores the original resid numbering destroyed by pdb4amber
    IMPORTANT: tleap names all chains A. Since I am selecting the chain based on the
    start and end resids and directly shift the numbering. It's important to shift
    the orginal chain A only at the end.
    :return:
    """

    # Assign chain ID to amber resid numbering
    chain = u.select_atoms(f"resid {chain_start} to {chain_end} and chainID A")
    chain.atoms.chainIDs = chainID

    # Shift resid numbering to old numbering
    chain.residues.resids += shift_factor


def remap(args, u):

    print("2. Import reside information")
    # Get chain information since amber deleted all chain IDs
    chains = chain2resid(args.mapping)

    # IMPORTANT: chain ID A needs to be at the end, otherwise it leads to misnumbering
    chains.sort_index(ascending=False, inplace=True)

    # 1. Assign chain Ids to amber resids
    for chainID, r in chains.iterrows():
        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"resid {r.amber_start} to {r.amber_end}")
        chain.atoms.chainIDs = chainID

    # 2. Renumber resids
    for chainID, r in chains.iterrows():

        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"chainID {chainID}")

        # Shift resid numbering to old numbering
        shift_factor = r.start - int(r.amber_start)
        chain.residues.resids += shift_factor

    return u


def center_and_export(u, args):
    # transform traj
    u = align_and_center(u)

    print(u)

    selectionAlgebra = 'all'  # protein and surface

    selection = u.select_atoms(selectionAlgebra)

    with MDAnalysis.Writer(args.topo_center, selection.n_atoms) as Traj:
        for ts in u.trajectory[0:1]:
            Traj.write(selection)

    # https://chemfiles.org/chemfiles/latest/formats.html
    # selection.write("c-alpha.pdb", frames=u.trajectory[0:1])
    # selection.write("c-alpha.cif", frames=u.trajectory[0:1], format='mmCIF')

    # last_frame = u.trajectory[-1]
    # selection.write("last_frame_topology.cif", format="cif")

    print("Print Traj")
    with MDAnalysis.Writer(args.traj_center, selection.n_atoms) as Traj:
        for ts in tqdm(u.trajectory):
            selection = u.select_atoms(selectionAlgebra)
            Traj.write(selection)

    u.atoms.write(args.gro)  #

    return u

def save_frame(u, out_frame, seleAlgebra='not name EPW', frameid=0):
    print(f"Save Frame {frameid}", out_frame)

    selection = u.select_atoms(seleAlgebra)

    with mda.Writer(out_frame, selection.n_atoms) as W:

        for frame_id, ts in enumerate(u.trajectory[frameid:frameid+1]):
            selection = u.select_atoms(seleAlgebra)
            W.write(selection)


def execute(command):
    """
    Executes commands in console
    :param command:
    :return:
    """

    output_text = subprocess.check_output(command, shell=True)
    return output_text

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
    command = f'run.py --pdb {frame_pdb} --resname ALA --resid 1 > {ligand_csv}'
    print(command)
    execute(command)

    # Analyze interaction of receptor to ligand
    #command = f'./7_interaction_csv.x {frame_pdb} {rec_resname} {rec_resid} > {receptor_csv}'

    command = f'run.py --pdb {frame_pdb} --resname SER --resid 632 > {receptor_csv}'
    execute(command)
    print(command)

import os
def export_martin(args, u):
    # TODO: Extract target from path
    #target = os.path.basename(args.traj).split('_')[0]
    target = 'C1s' # TODO
    number_frames = int(args.frame_id)

    traj_length = len(u.trajectory)

    # Export centered frames for Martin Analyzer
    print("Start extracting last frames")
    selection = f'(protein or byres around 8.0 protein) and not name EPW'          # protein and surface
    selection = 'protein'          # TODO: include water
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False,help='', default='last.pdb')
    parser.add_argument('--traj', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--mapping', required=False, help='working directory', default='output/demo/C1s_BD001/WT/amber/amber_renum.txt')

    # Output
    parser.add_argument('--traj_center', required=False,help='', default='center.dcd')
    parser.add_argument('--topo_center', required=False,help='', default='topo.pdb')
    parser.add_argument('--gro', required=False,help='', default='topo.gro')

    # Input
    parser.add_argument('--frame_id', required=False, help='working directory')
    parser.add_argument('--dir', required=False, help='working directory')

    # Output
    parser.add_argument('--frame_pdb', required=False, help='working directory')
    parser.add_argument('--ligand_csv', required=False, help='working directory')
    parser.add_argument('--receptor_csv', required=False, help='working directory')


    args = parser.parse_args()

    print(args)

    print("1. Import")


    # Import Trajectory
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)

    # Remap Residues
    u = remap(args, u)

    # Center protein
    center_and_export(u, args)

    # Export protein
    export_martin(args, u)
