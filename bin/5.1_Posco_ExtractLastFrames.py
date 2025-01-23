#!/usr/bin/env python

"""
This script processes molecular dynamics trajectories and performs interaction analysis between a ligand and a receptor.
"""

import argparse
from Helper import remap_MDAnalysis  # Helper functions for execution and MDAnalysis remapping
import MDAnalysis as mda  # MDAnalysis for atom selection and structure manipulation
import openmm.app as app


def extract_binding_surface(u, t=8):
    """
    Extracts the protein from the frame plus all complete water molecules t=8 Angstrom from the binding
    surface
    """

    # Select chain A (must be always ligand) and everything else
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
    parser.add_argument('--frame', type=int,required=False, help='PDB file for the ligand and receptor')
    parser.add_argument('--lig_frame', required=False, help='PDB file for the ligand')
    parser.add_argument('--rec_frame', required=False, help='PDB file for the receptor')

    return parser.parse_args()

if __name__ == '__main__':
    # Parse command-line arguments
    args = parse_arguments()

    # Import Trajectory #TODO export to helpers
    # TODO: I am aware that it is slow to open the whole trajectory in every frame, but snakemake doesn't really like expanding output files
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)
    u = remap_MDAnalysis(u, topo)

    # Extract frame required
    ts = u.trajectory[-args.frame - 1]

    # Extract protein and water in binding surface
    print(f"Processing frame {args.frame}: {ts.frame}")
    (ligand, receptor) = extract_binding_surface(u)

    # Save ligand and receptor files separatly
    ligand.write(args.lig_frame)
    receptor.write(args.rec_frame)
            