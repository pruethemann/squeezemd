#!/usr/bin/env python
"""
This script processes molecular dynamics trajectories and performs interaction analysis between a ligand and a receptor.
"""

import argparse
from Helper import remap_MDAnalysis  # Helper functions for execution and MDAnalysis remapping
import MDAnalysis as mda             # MDAnalysis for atom selection and structure manipulation
import openmm.app as app
import pandas as pd

def extract_sequence(ligand, receptor, sequence_file):
    """
    Extract the amino acid sequence from the structure for the ligand and receptor and saves it as parquet.
    """

    # Extract sequence
    seq_ligand = {"resid": ligand.residues.resids,
                  "resname": ligand.residues.resnames}
    
    seq_receptor = {"resid": receptor.residues.resids,
                    "resname": receptor.residues.resnames}
    
    seq_ligand = pd.DataFrame(seq_ligand)
    seq_receptor = pd.DataFrame(seq_receptor)

    seq_ligand['protein'] = 'lig'
    seq_receptor['protein'] = 'rec'

    seq = pd.concat([seq_ligand, seq_receptor])
    seq = seq.set_index(['resid', 'resname'])
    seq.to_parquet(sequence_file)

def extract_binding_surface(u, t=8):
    """
    Extracts the ligand (segid A or X), receptor, and all complete water molecules within t Angstrom
    from the binding surface.
    A: protein ligand
    X: small molecule ligand
    """

    # Determine the ligand segid
    ligand = u.select_atoms('segid A')
    if len(ligand) == 0:        # For small molecule the chain ID is X
        ligand = u.select_atoms('segid X')
        ligand_segid = 'X'
    else:
        ligand_segid = 'A'

    # Select chain A (must be always ligand) and everything else
    receptor = u.select_atoms(f'not segid {ligand_segid} and protein')

    # Extract and save sequences information for posco
    extract_sequence(ligand, receptor, args.sequence)

    # Select water molecules within 5 Å of both chain A and chain B
    water_binding_site = u.select_atoms(f'resname HOH and (around {t} segid {ligand_segid}) and (around {t} (not segid {ligand_segid} and protein))')

    # Get the residues of selected water molecules
    water_residues = water_binding_site.residues

    # Filter out incomplete water molecules (keep only those with exactly 3 atoms)
    complete_water_residues = water_residues[[len(res.atoms) == 3 for res in water_residues]]

    # Get the atoms of the complete water molecules
    complete_water = complete_water_residues.atoms

    # Combine all selections
    return (ligand, receptor + complete_water)

def parse_arguments():
    """
    Parse command-line arguments for the script.
    :return: Parsed arguments.
    """
    # Initialize argument parser
    parser = argparse.ArgumentParser()

    # Add arguments for input files, output options, and parallelization settings
    parser.add_argument('--topo', required=False, help='', default='frame_end.cif')
    parser.add_argument('--traj', required=False, help='', default='trajectory.dcd')
    parser.add_argument('--frame', type=int,required=False, help='PDB file for the ligand and receptor')

    # Output
    parser.add_argument('--lig_frame', required=False, help='PDB file for the ligand', default='lig.pdb')
    parser.add_argument('--rec_frame', required=False, help='PDB file for the receptor', default='rec.pdb')
    parser.add_argument('--sequence', required=False, help='PDB file for the receptor', default='sequence.parquet')

    return parser.parse_args()

if __name__ == '__main__':
    # Parse command-line arguments
    args = parse_arguments()

    # Import Trajectory
    # I am aware that it is slow to open the whole trajectory in every frame, but snakemake doesn't really like expanding output files
    # TODO Consider exporting all frames in one go.
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)

    # Define residues and chains according to pdb
    u = remap_MDAnalysis(u, topo)

    # Extract frame required
    ts = u.trajectory[-args.frame - 1]

    # Extract protein and water in binding surface
    print(f"Processing frame {args.frame}: {ts.frame}")
    (ligand, receptor) = extract_binding_surface(u)

    # Save ligand and receptor files separatly
    ligand.write(args.lig_frame)
    receptor.write(args.rec_frame)
            