#!/usr/bin/env python

"""
This script processes molecular dynamics trajectories and performs interaction analysis between a ligand and a receptor.
"""

import argparse, os
from Helper import execute, remap_MDAnalysis  # Helper functions for execution and MDAnalysis remapping
import MDAnalysis as mda  # MDAnalysis for atom selection and structure manipulation
import openmm.app as app


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
    return (ligand, receptor + complete_water)

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

    frame_count = 0
    for ts in u.trajectory[-args.n_frames:]:
        # Extract protein and water in binding surface
        (ligand, receptor) = extract_binding_surface(u, frame_count)

        # Save ligand and receptor files separatly
        ligand_pdb = os.path.join(args.dir, f'lig_{frame_count}.pdb')
        ligand.write(ligand_pdb)

        receptor_pdb = os.path.join(args.dir, f'rec_{frame_count}.pdb')
        receptor.write(receptor_pdb)

        frame_count += 1
            

    # Save the last frame
    u.trajectory[-1]

    # Save the last frame as a pdb file
    with mda.Writer(args.final, n_atoms=u.atoms.n_atoms) as W:
        W.write(u)
