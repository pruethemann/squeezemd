#!/usr/bin/env python

"""
This script processes molecular dynamics trajectories and performs interaction analysis between a ligand and a receptor.
"""

import argparse, os
from Helper import execute, remap_MDAnalysis  # Helper functions for execution and MDAnalysis remapping
import MDAnalysis as mda  # MDAnalysis for atom selection and structure manipulation
import openmm.app as app

def extract_complex_resids(pdb_ligand: os.path, chainID: str):
    """
    Extract the residue name and residue ID of the first residue of the specified chain in the pdb file.
    :param pdb_ligand: The PDB file path.
    :param chainID: The chain identifier (e.g., 'A' for ligand, 'B' for receptor).
    :return: Tuple containing the residue name and residue ID of the first residue in the chain.
    """
    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # Extract the ligand or receptor based on the provided chainID
    protein = u.select_atoms(f'segid {chainID}')

    # Return the residue name and ID of the first residue in the chain
    return (protein.residues[0].resname, protein.residues[0].resid)

def extract_binding_surface(u, frame_count, t=8):
    """
    Extracts the protein from the frame plus all complete water molecules t=8 Angstrom from the binding
    surface
    """
    
    print(f"Processing frame {frame_count}: {ts.frame}")

    # Select chain A and chain B
    ligand = u.select_atoms('segid A')
    receptor = u.select_atoms('not segid A and protein')

    # Select water molecules within 5 Å of both chain A and chain B
    water_binding_site = u.select_atoms(f'resname HOH and (around {t} segid A) and (around {t} (not segid A and protein))')

    # Get the residues of selected water molecules
    water_residues = water_binding_site.residues

    # Filter out incomplete water molecules (keep only those with exactly 3 atoms)
    complete_water_residues = water_residues[[len(res.atoms) == 3 for res in water_residues]]

    # Get the atoms of the complete water molecules
    complete_water = complete_water_residues.atoms

    # Combine all selections
    return ligand + receptor + complete_water

def parse_arguments():
    """
    Parse command-line arguments for the script.
    :return: Parsed arguments.
    """
    # Initialize argument parser
    parser = argparse.ArgumentParser()

    # Add arguments for input files, output options, and parallelization settings
    parser.add_argument('--topo', required=False, help='', default='trajectory.dcd')
    parser.add_argument('--traj', required=False, help='', default='trajectory.dcd')
    parser.add_argument('--n_frames', required=False, help='The number of frames exported from the trajectory', default=10, type=int)
    parser.add_argument('--dir', required=False, help='The working directory for analysis', default='tmp')
    parser.add_argument('--final', required=False, help='', default='trajectory.dcd')
    parser.add_argument('--pdb', required=False, help='PDB file for the ligand and receptor')

    return parser.parse_args()

if __name__ == '__main__':
    # Parse command-line arguments
    args = parse_arguments()

    # Import Trajectory #TODO export to helpers
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)
    u = remap_MDAnalysis(u, topo)

    # Extract ligand and receptor residue information from the PDB
    (resname_lig, resid_lig) = extract_complex_resids(args.pdb, 'A')  # For ligand (chain A)
    (resname_rec, resid_rec) = extract_complex_resids(args.pdb, 'B')  # For receptor (chain B)

    frame_count = 0
    for ts in u.trajectory[-args.n_frames:]:
        # Extract protein and water in binding surface
        bind_surface_sele = extract_binding_surface(u, frame_count)

        # Save the selection of the current frame as a PDB file
        pdb_binding_surface = os.path.join(args.dir, f'frame_{frame_count}.pdb')
        bind_surface_sele.write(pdb_binding_surface)

        # Perform interaction analyzer
        lig_csv = os.path.join(args.dir, 'lig', f'{frame_count}.csv')
        rec_csv = os.path.join(args.dir, 'rec', f'{frame_count}.csv')

        # Analyze interactions of the ligand
        command = f'interaction-analyzer-csv.x {pdb_binding_surface} {resname_lig} {resid_lig} > {lig_csv}'
        execute(command)

        # Analyze interactions of the receptor
        command = f'interaction-analyzer-csv.x {pdb_binding_surface} {resname_rec} {resid_rec} > {rec_csv}'
        execute(command)

        frame_count += 1
            

    # Save the last frame
    u.trajectory[-1]

    # Save the last frame as a pdb file
    with mda.Writer(args.final, n_atoms=u.atoms.n_atoms) as W:
        W.write(u)
