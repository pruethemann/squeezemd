#!/usr/bin/env python

"""

"""
import argparse
import os
from Helper import execute, remap_MDAnalysis
import mdtraj
import numpy as np
import multiprocessing
import MDAnalysis as mda


def extract_complex_resids(pdb_ligand: os.path, chainID:str):
    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # Extract ligand at chain I
    protein = u.select_atoms(f'segid {chainID}')

    # Return First resname and position
    return (protein.residues[0].resname, protein.residues[0].resid)

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
    command = f'interaction-analyzer-csv.x {frame_pdb} {resname_lig} {resid_lig} > {ligand_csv}'
    execute(command)

    # Analyze interaction of receptor to ligand
    command = f'interaction-analyzer-csv.x {frame_pdb} {resname_rec} {resid_rec} > {receptor_csv}'
    execute(command)

    # TODO: Error handling check if ligand_csv and receptor_csv are not empty because analysis failed


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

def extract_protein(frame_number:int):
    frame_id = -args.n_frames + frame_number

    water_sele = extract_protein_water_shell(traj[frame_id], 0.8)
    # Save the new trajectory as a DCD file
    frame_path = os.path.join(args.dir, f'frame_{frame_number}.pdb')

    water_sele.save(frame_path)

    # Execute Martin interaction analyzer
    lig_csv = os.path.join(args.dir, 'lig', f'{frame_number}.csv')
    rec_csv = os.path.join(args.dir, 'rec', f'{frame_number}.csv')

    interaction_analyzer(frame_path, lig_csv, rec_csv)

    return "Sucess"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--traj', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--n_frames', required=False, help='The last number of frames exported from the trajectory', default=10, type=int)
    parser.add_argument('--dir', required=False, help='The working dir for the analysis', default='tmp')
    parser.add_argument('--final', required=False,help='', default='trajectory.dcd')
    parser.add_argument('--cpus', required=False, help='', default=1, type=int)
    parser.add_argument('--pdb', required=False, help='')

    args = parser.parse_args()

    # Import Trajecotry
    traj = mdtraj.load(args.traj, top=args.topo)

    (resname_lig, resid_lig) = extract_complex_resids(args.pdb, 'I')
    (resname_rec, resid_rec) = extract_complex_resids(args.pdb, 'A')

    # TODO DEBUG
    print(resname_rec, resid_rec)
    print(resname_lig, resid_lig)

    # Export the last n_frames as pdb files
    # Only the protein and 8 Angstrom around protein is exported
    # frame_number is number from 0:n_frames, frame_id corresponds to number in traj from the end
    # Process i parallized
    with multiprocessing.Pool(args.cpus) as pool:
        for frame_number in pool.map(extract_protein, range(args.n_frames)):
            print(frame_number)

    pool.close()

    # Export last centered frame
    traj[-1].save(args.final)
