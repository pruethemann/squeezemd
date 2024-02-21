#!/usr/bin/env python

import argparse
import prolif as plf
import MDAnalysis as mda
import openmm.app as app
from Helper import remap_MDAnylsis


def create_interactionFingerprint(args):

    # Import trajectory
    topo = app.PDBxFile(args.topo)      # Transform cif to MDAnalysis topology
    u = mda.Universe(topo, args.traj, in_memory=False)

    u = remap_MDAnylsis(u, topo)

    # Define ligand (gigastasin) and receptor
    ligand = u.select_atoms("chainID I")            # TODO Fix renaming chainIDs
    protein = u.select_atoms("chainID A or chainID B")

    # Run interaction fingerprint analysis
    fp = plf.Fingerprint(["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "PiCation", "CationPi", "Anionic", "Cationic"], count=True)
    # TODO: introduce multi threads
    fp.run(u.trajectory[-args.n_frames:], ligand, protein,  n_jobs=args.threads)

    # Export interactions
    interactions_df = fp.to_dataframe()

    # Save metadata with parquet
    interactions_df.attrs['description'] = 'Fingerprint metadata'
    interactions_df.attrs['complex'] = args.complex
    interactions_df.attrs['mutation'] = args.mutation
    interactions_df.attrs['target'] = args.complex.split('_')[0]
    interactions_df.attrs['ligand'] = args.complex.split('_')[0]
    interactions_df.attrs['seed'] = args.seed

    interactions_df.to_parquet(args.output)

    interactions_df.to_csv(args.output[:-8] + '.csv')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=True,help='Topo file as cif', default='topo.cif')
    parser.add_argument('--traj', required=True,help='Trajectory file', default='traj.dcd')

    # Input parameters
    parser.add_argument('--n_frames', required=False, help='How many frames to be analysed. Only the last n frames will be analyzed. Default 100', type=int, default=100)
    parser.add_argument('--threads', required=False, help='Optional: Number of threads: Default 4', type=int, default=4)
    parser.add_argument('--complex', required=False, help='Optional: Name of Target and Ligand as Target_Ligand', default='Target_Ligand')
    parser.add_argument('--mutation', required=False, help='Optional: Mutation in ligand', default='Wildtype')
    parser.add_argument('--seed', required=False, help='Optional: Seed used during Molecular Dynamics simulation', type=int, default=-1)

    # Output
    parser.add_argument('--output', required=False, help='Path of parquet datafile including all fingerprint and meta data', default='output.parquet')

    args = parser.parse_args()
    create_interactionFingerprint(args)
