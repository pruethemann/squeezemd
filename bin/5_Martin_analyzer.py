#!/usr/bin/env python

"""

"""
import argparse
import os
import MDAnalysis as mda
import MDAnalysis
import MDAnalysis.transformations as trans
import openmm.app as app
from Helper import remap, execute
import mdtraj
import numpy as np



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





def extract_protein_water_shell(traj, cutoff=0.5):

    # Select protein and water
    protein = traj.topology.select('protein')
    water = traj.topology.select('water')

    # Determine all water indices which are in a distance < cutoff
    water_indices = mdtraj.compute_neighbors(traj, cutoff=cutoff, query_indices=protein, haystack_indices=water)
    water_indices = np.array(water_indices)

    # Incomplete water molecules need to be restored
    # Assuming each water molecule consists of 1 oxygen and 2 hydrogens
    complete_water_indices = []
    for frame in water_indices:
        for atom_idx in frame:
            atom = traj.topology.atom(atom_idx)
            if atom.element.symbol == 'O':  # If the atom is an oxygen
                # Get the indices of the water molecule to which this oxygen belongs
                water_molecule = [atom.index for atom in atom.residue.atoms]
                complete_water_indices.append(water_molecule)

    # Deduplicate and flatten the list
    complete_water_indices = np.array(complete_water_indices).flatten()

    # combine water shell and protein
    # ToDO add salts
    combined_indices = np.concatenate([protein, complete_water_indices])

    # Extract the trajectory of complete water molecules
    water_shell = traj.atom_slice(combined_indices)
    return water_shell



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--traj', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--n_frames', required=False, help='The last number of frames exported from the trajectory', default=10)
    parser.add_argument('--dir', required=False, help='The working dir for the analysis', default='tmp')
    parser.add_argument('--final', required=False,help='', default='trajectory.dcd')

    args = parser.parse_args()

    # Import Trajecotry
    traj = mdtraj.load(args.traj, top=args.topo)

    # Slice only last part


    # Export protein
    for frame_id in range(int(args.n_frames)):
        water_sele = extract_protein_water_shell(traj[-frame_id-1], 0.8)
        # Save the new trajectory as a DCD file
        # water_sele.save_dcd('output/water_around_protein.dcd')
        frame_path = os.path.join(args.dir, f'frame_{frame_id}.pdb')
        lig_csv = os.path.join(args.dir, 'lig', f'{frame_id}.csv')
        rec_csv = os.path.join(args.dir, 'rec', f'{frame_id}.csv')

        water_sele.save(frame_path)

        # Execute Martin interaction analyzer
        lig_csv = os.path.join(args.dir, 'lig', f'{frame_id}.csv')
        rec_csv = os.path.join(args.dir, 'rec', f'{frame_id}.csv')

        interaction_analyzer(frame_path, lig_csv, rec_csv)


    # Export last centered frame
    traj.save(args.final)
