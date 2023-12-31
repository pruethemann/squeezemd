#!/usr/bin/env python

"""
Script takes a Molecular Dynamics trajectory and performs the following tasks:

- Centers the protein in the center of the water box
- export the last n center frames in pdb format
- the last frames contain the protein and all water molecules 8 Angstrom around the binding site
- The Martin interaction-analysis tool is executed on the last n frames

TODO: Implement from MDAnalysis.analysis import contacts
"""
import argparse
import os
import MDAnalysis as mda
import MDAnalysis
import MDAnalysis.transformations as trans
import openmm.app as app
from Helper import remap, execute


def center_trajectory(u, args):
    """
    Function centers protein in the middle of the box and wraps water around it.
    All frames are aligned to the first frame.
    Centered trajectory is saved as dcd. And as pdb as topology.
    However: Due to the fact that the standard pdb format can only save
    a limited amount of atoms. This will lead to the fact that multiple atoms
    have the same id.
    :param u:
    :param args:
    :return: centered Universe
    """
    # Align and center the trajectory, create atomgroups to represent the protein, non-protein, and c-alpha sections of the system
    protein = u.select_atoms('protein')
    not_protein = u.select_atoms('not protein')

    # next we create a chain of transformations which will be applied each time a frame is loaded from the Universe
    # Note: the transformations are applied linearly, so the order of the transformations matters!

    transforms = [trans.unwrap(protein),
                  trans.center_in_box(protein, wrap=True),
                  trans.wrap(not_protein)
                  ]


    # Apply the transformation. Transformation are execute only if trajectory is accessed
    u.trajectory.add_transformations(*transforms)

    # Save centered topology and trajectory
    selection = u.select_atoms('all')
    selection.write(args.topo_center, frames=u.trajectory[0:1])
    selection.write(args.traj_center, frames='all')

    return u


def interaction_analyzer(frame_pdb, ligand_csv, receptor_csv):
    """
    Execute Martin's interaction analyzer.
    In a first step the Analyzer is executed on the pdb file from the ligand perspective
    then from the receptor perspective.
    :param frames_nr:
    :param args:
    :return:
    """

    # Analyze interactions of ligand to receptor
    command = f'interaction-analyzer-csv.x {frame_pdb} ALA 1 > {ligand_csv}'
    execute(command)

    # Analyze interaction of receptor to ligand
    command = f'interaction-analyzer-csv.x {frame_pdb} SER 632 > {receptor_csv}'
    execute(command)


def export_martin_interaction(args, u):
    """
    Exports the last n frames of the trajectory containing the protein and
    all water molecules 8 Angstrom from the binding site.
    The Martin interaction analyzer is executed on the last frames.
    :param args:
    :param u:
    :return:
    """

    number_frames = int(args.n_frames)
    traj_length = len(u.trajectory)

    # Export centered frames for Martin Analyzer
    print(f"Start extracting the frames {traj_length-number_frames}:{traj_length}")

    # Select everything around ligand (chain I) and receptor (chain B)
    # Export all the protein and all the water 8 Angstrom the binding site.
    residues_near_I = u.select_atoms('byres around 8.0 (chainID I)')
    residues_near_B = u.select_atoms('byres around 8.0 (chainID B)')
    binding_site = u.select_atoms('(group residues_near_I and group residues_near_B) or protein',
        residues_near_B=residues_near_B, residues_near_I=residues_near_I)

    frame_id = 0
    for i in range(len(u.trajectory) - number_frames, len(u.trajectory), 1):
        frame_path = os.path.join(args.dir, f'frame_{frame_id}.pdb')
        lig_csv = os.path.join(args.dir, 'lig', f'{frame_id}.csv')
        rec_csv = os.path.join(args.dir, 'rec', f'{frame_id}.csv')
        binding_site.write(frame_path, frames=u.trajectory[i:i + 1])

        interaction_analyzer(frame_path, lig_csv, rec_csv)
        frame_id += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False,help='', default='topo.cif')
    parser.add_argument('--traj', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--mapping', required=False, help='Amber residue mapping', default='amber_renum.txt')
    parser.add_argument('--n_frames', required=False, help='The last number of frames exported from the trajectory')
    parser.add_argument('--dir', required=False, help='The working dir for the analysis')

    # Output
    parser.add_argument('--traj_center', required=False,help='', default='traj_center.dcd')
    parser.add_argument('--topo_center', required=False,help='', default='topo_center.pdb')

    args = parser.parse_args()

    # Import Trajectory
    topo = app.PDBxFile(args.topo)  # Convert cif to internal format
    u = mda.Universe(topo, args.traj, in_memory=False)

    # Remap Residues
    u = remap(args.mapping, u)

    # Center protein
    u_centered = center_trajectory(u, args)

    # Export protein
    export_martin_interaction(args, u_centered)
