#!/usr/bin/env python

import argparse
import prolif as plf
import MDAnalysis as mda

import openmm.app as app
from Helper import remap

def create_interactionFingerprint(args):

    # Import trajectory
    topo = app.PDBxFile(args.topo)      # Transform cif to MDAnalysis topology
    u = mda.Universe(topo, args.traj, in_memory=False)

    # Remap amber residues to original residue numbers
    u = remap(args.mapping, u)

    # Define ligand (gigastasin) and receptor
    ligand = u.select_atoms("chainID I")
    protein = u.select_atoms("chainID A or chainID B")

    # Run interaction fingerprint analysis
    fp = plf.Fingerprint(["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "PiCation", "CationPi", "Anionic", "Cationic"])
    fp.run(u.trajectory[-args.n_frames:], ligand, protein)

    # Export interactions
    interactions_df = fp.to_dataframe()
    interactions_df.to_feather(args.output)
    interactions_df.to_csv(args.output[:-8] + '.csv')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False,help='Topo file as cif', default='topo.cif')
    parser.add_argument('--traj', required=False,help='Trajectory file', default='traj.dcd')
    parser.add_argument('--mapping', required=False, help='Mapping resid to amber resid', default='amber_remap.txt')
    parser.add_argument('--n_frames', required=False, help='Number of last frames to be analyzed', type=int, default=100)

    # Output
    parser.add_argument('--output', required=False, help='Mapping resid to amber resid', default='')

    args = parser.parse_args()
    create_interactionFingerprint(args)
