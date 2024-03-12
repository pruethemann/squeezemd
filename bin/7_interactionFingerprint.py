#!/usr/bin/env python

import argparse
import os
import prolif as plf
import MDAnalysis as mda
import openmm.app as app
from Helper import remap_MDAnalysis, extract_ligand_sequence

def create_interaction_fingerprint(topology_file, trajectory_file, output_file, n_frames=100, threads=4, complex_name='Target_Ligand', mutation='Wildtype', seed=-1):
    """
    Generates and exports interaction fingerprints for a given molecular dynamics simulation.

    Args:
        topology_file (str): Path to the topology file in CIF format.
        trajectory_file (str): Path to the molecular dynamics trajectory file.
        output_file (str): Path for the output data file (Parquet and CSV formats).
        n_frames (int, optional): Number of frames from the end of the trajectory to analyze. Defaults to 100.
        threads (int, optional): Number of threads to use for parallel processing. Defaults to 4.
        complex_name (str, optional): Name of the complex, formatted as 'Target_Ligand'. Defaults to 'Target_Ligand'.
        mutation (str, optional): Description of any mutation in the ligand. Defaults to 'Wildtype'.
        seed (int, optional): Seed used during the molecular dynamics simulation. Defaults to -1.

    This function selects the ligand and protein atoms, runs the interaction fingerprint analysis, and exports the results with metadata.
    """
    # Ensure input files exist
    if not os.path.exists(topology_file) or not os.path.exists(trajectory_file):
        raise FileNotFoundError("Topology or trajectory file not found.")

    # Setup analysis environment
    topo = app.PDBxFile(topology_file)
    universe = mda.Universe(topo, trajectory_file, in_memory=False)
    universe = remap_MDAnalysis(universe, topo)

    # Selecting ligand and protein
    # CHAINIDENTIFICAITON
    ligand = universe.select_atoms("chainID A")
    protein = universe.select_atoms("chainID B")

    # Initialize interaction fingerprint analysis
    fingerprint = plf.Fingerprint(["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "PiCation", "CationPi", "Anionic", "Cationic"], count=True)
    fingerprint.run(universe.trajectory[-n_frames:], ligand, protein, n_jobs=threads)

    # Extract number of residues in ligand
    n_residues_ligand = len( str(ligand.residues.sequence().seq))

    # Exporting results
    interactions_df = fingerprint.to_dataframe()
    interactions_df.attrs = {
        'description': 'Fingerprint metadata',
        'complex': complex_name,
        'mutation': mutation,
        'target': complex_name.split('_')[0],
        'ligand': complex_name.split('_')[1],
        'seed': seed,
        'n_residues_ligand': n_residues_ligand
    }



    interactions_df.to_parquet(output_file)
    interactions_df.to_csv(output_file.replace('.parquet', '.csv'))

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate interaction fingerprints from molecular dynamics simulations.")
    # Inputs
    parser.add_argument('--topo', required=True, help='Topology file in CIF format')
    parser.add_argument('--traj', required=True, help='Trajectory file')
    # Parameters
    parser.add_argument('--n_frames', type=int, default=100, help='Number of frames to analyze from the end. Default: 100')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for analysis. Default: 4')
    parser.add_argument('--complex', default='Target_Ligand', help='Name of target and ligand as Target_Ligand. Default: "Target_Ligand"')
    parser.add_argument('--mutation', default='Wildtype', help='Mutation in ligand. Default: "Wildtype"')
    parser.add_argument('--seed', type=int, default=-1, help='Seed used during simulation. Default: -1')
    # Output
    parser.add_argument('--output', default='output.parquet', help='Output file path. Default: "output.parquet"')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    create_interaction_fingerprint(args.topo, args.traj, args.output, args.n_frames, args.threads, args.complex, args.mutation, args.seed)
