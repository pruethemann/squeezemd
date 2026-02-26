#!/usr/bin/env python
"""Run PoSCo on trajectory frames and parse interaction results.

This script extracts ligand/receptor subsets from trajectory frames,
executes PoSCo, and consolidates interactions into a parquet table with
metadata (complex, mutation, seed, frame).
"""

import argparse, os
from ...helper_functions import remap_MDAnalysis, execute # Helper functions for execution and MDAnalysis remapping
import openmm.app as app
import pandas as pd

import warnings

warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"DCDReader currently makes independent timesteps"
)

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Found no information for attr: '.*' Using default value of '.*'"
)

import MDAnalysis as mda

def parse_lipophilic(parts, sequence):
    """Parse a PoSCo lipophilic interaction line into a dict."""
        
    interaction_info = parts[0].split()
    donor_acceptor = parts[1].strip().split()

    ligand_atom = donor_acceptor[0]
    ligand_resname = donor_acceptor[1]
    ligand_resid = int(donor_acceptor[2])

    receptor_atom = donor_acceptor[-3]
    receptor_resname = donor_acceptor[-2]
    receptor_resid = int(donor_acceptor[-1])

    # if any of the conditions holds, swap everything
    should_swap = (
        (ligand_resname == 'HOH' and sequence.loc[(receptor_resid, receptor_resname)]['protein'] == 'lig')  or
        (receptor_resname == 'HOH' and sequence.loc[(ligand_resid, ligand_resname)]['protein'] == 'rec')    or
        (sequence.loc[(ligand_resid, ligand_resname)]['protein'] == 'rec')
    )

    if should_swap:
    # swap ligand ↔ receptor
        (ligand_resid, receptor_resid) = (receptor_resid, ligand_resid)
        (ligand_resname, receptor_resname) = (receptor_resname,ligand_resname)
        (ligand_atom, receptor_atom) = (receptor_atom, ligand_atom)

    distance = float(interaction_info[2].split("=")[1])
    energy = float(interaction_info[3].split("=")[1])

    interaction = {
        "Interaction Type": 'lipophilic',
        "Distance (r)": distance,
        "Energy (e)": energy,
        'receptor_resname' : receptor_resname,
        'receptor_resid' : receptor_resid,
        'ligand_resname' : ligand_resname,
        'ligand_resid' : ligand_resid,
        'receptor_atom' : receptor_atom,
        'ligand_atom' : ligand_atom,
    }

    return interaction


def parse_hbonds(parts, sequence):
    """Parse a PoSCo H‑bond interaction line into a dict."""

    interaction_info = parts[0].split()
    donor_acceptor = parts[1].strip().split()

    ligand_atom = donor_acceptor[0]
    ligand_resname = donor_acceptor[1]
    ligand_resid = int(donor_acceptor[2])

    receptor_atom = donor_acceptor[-3]
    receptor_resname = donor_acceptor[-2]
    receptor_resid = int(donor_acceptor[-1])

    #print(sequence)

    # TODO: That is only necessary because in posco I can't differeniate between ligand and receptors
    # TODO. Do this swap only once
    should_swap = (
        (ligand_resname == 'HOH' and sequence.loc[(receptor_resid, receptor_resname)]['protein'] == 'lig')  or
        (receptor_resname == 'HOH' and sequence.loc[(ligand_resid, ligand_resname)]['protein'] == 'rec')    or
        (sequence.loc[(ligand_resid, ligand_resname)]['protein'] == 'rec')
    )

    if should_swap:
    # swap ligand ↔ receptor
        (ligand_resid, receptor_resid) = (receptor_resid, ligand_resid)
        (ligand_resname, receptor_resname) = (receptor_resname,ligand_resname)
        (ligand_atom, receptor_atom) = (receptor_atom, ligand_atom)

    distance = float(interaction_info[2].split("=")[1])
    angle = float(interaction_info[3].split("=")[1])
    energy = float(interaction_info[4].split("=")[1])

    # Include salt bridge data
    marked =  "marked as salt-bridge" in parts[2]

    interaction = {
        "Interaction Type": 'H-bond',
        "Distance (r)": distance,
        "Angle (a)": angle,
        "Energy (e)": energy,
        'receptor_atom' : receptor_atom,
        'receptor_resname' : receptor_resname,
        'receptor_resid' : receptor_resid,
        'ligand_atom' : ligand_atom,
        'ligand_resname' : ligand_resname,
        'ligand_resid' : ligand_resid,
        "Marked as Salt-Bridge": marked
    }

    return interaction

# Parse the input data into a pandas DataFrame
def parse_posco(posco_output, metadata, frame_id, sequence):
    """
    parse the posco text file and extract relevant interaction data and save
    as parquet.
    
    :param posco_output: Description
    :param metadata: Description
    :param frame_id: Description
    :param sequence_parquet: Description
    """
    data = []

    metadata['target'] = metadata['complex'].split('_')[0]
    metadata['ligand'] = metadata['complex'].split('_')[1]

    # In rare cases the same resname and resid can exist in rec and lig.
    # Keep the first occurrence to avoid ambiguity.
    if not sequence.index.is_unique:
        sequence = sequence[~sequence.index.duplicated(keep='first')]
    
    with open(posco_output, 'r') as file:
        for line in file:

            if line.startswith("Lipo_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_lipophilic(parts, sequence)
                data.append(interaction)

            if line.startswith("HB_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_hbonds(parts, sequence)
                data.append(interaction)

    data = pd.DataFrame(data)

    # Determine metrics lables
    data['name'] = metadata['complex']
    data['target'] = metadata['target']
    data['lig'] = metadata['ligand']
    data['mutation'] = metadata['mutation']
    data['frame'] = frame_id
    data['seed'] = metadata['seed']

    return data

def extract_sequence(ligand, receptor):
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

    return seq
    #seq.to_parquet(sequence_file)

def extract_binding_surface(u, t=8):
    """
    Extracts the ligand (segid A or X), receptor, and all complete water molecules within t Angstrom
    from the binding surface.
    A: protein ligand
    X: small molecule ligand
    """

    # Determine the ligand segid (A for protein ligand, X for small molecule)
    ligand = u.select_atoms('segid A')
    if len(ligand) == 0:        # For small molecule the chain ID is X
        ligand = u.select_atoms('segid X')
        ligand_segid = 'X'
    else:
        ligand_segid = 'A'

    # Select ligand and receptor proteins
    receptor = u.select_atoms(f'not segid {ligand_segid} and protein')

    # Extract and save sequences information for posco
    sequence = extract_sequence(ligand, receptor)

    # Select water molecules within t Å of ligand and receptor
    water_binding_site = u.select_atoms(f'resname HOH and (around {t} segid {ligand_segid}) and (around {t} (not segid {ligand_segid} and protein))')

    # Get the residues of selected water molecules
    water_residues = water_binding_site.residues

    # Filter out incomplete water molecules (keep only those with exactly 3 atoms)
    complete_water_residues = water_residues[[len(res.atoms) == 3 for res in water_residues]]

    # Get the atoms of the complete water molecules
    complete_water = complete_water_residues.atoms

    # Combine all selections
    return (ligand, receptor + complete_water, sequence)

def parse_arguments():
    """
    Parse command-line arguments for the script.
    :return: Parsed arguments.
    """
    # Initialize argument parser
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False, help='Topology file for trajectory in cif format', default='center/structure_end.cif')
    parser.add_argument('--traj', required=False, help='Centered trajectory', default='center/trajectory_centered.dcd')
    
    # Parameters
    parser.add_argument('--number_frames', type=int,required=False, help='Number of last frames to be extracted', default=1)
    parser.add_argument('--complex', required=False, help='', default='C1s_Gigastasin')
    parser.add_argument('--mutation', required=False, help='', default='R65E')
    parser.add_argument('--seed', type=int,required=False, help='', default=222)

    # Output
    parser.add_argument('--posco_interaction', required=False, help='', default='posco.txt')
    parser.add_argument('--posco_parquet', required=False, help='', default='posco.parquet')

    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()

    metadata = {'complex': args.complex,
                'mutation': args.mutation,
                'seed': args.seed
    }

    prefix = f'{args.complex}_{args.mutation}_{args.seed}' # used for tmp file paths

    # Import Trajectory
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)

    # Define residues and chains according to pdb
    u = remap_MDAnalysis(u, topo)

    # Make sure masses and types are correct
    u.guess_TopologyAttrs(to_guess=["masses", "types"])

    posco_interactions = []

    for i in range(args.number_frames):
        # 1. Extract ligand and receptor for this frame
        ts = u.trajectory[-i - 1]

        # Extract protein and water in binding surface
        print(f"Processing frame {i}: {ts.frame}")

        (ligand, receptor, sequence) = extract_binding_surface(u)

        # Save ligand and receptor files separatly
        lig_path = f'.{i}_lig_{prefix}.pdb'
        rec_path = f'.{i}_rec_{prefix}.pdb'

        ligand.write(lig_path)
        receptor.write(rec_path)

        # 2. Perform PoSCo
        posco_result = f'{i}_posco_{prefix}.txt'
        cmd = f"po-sco {rec_path} {lig_path} -b  > {posco_result}"
        execute(cmd)

        # 3. Parse interactions into a single parquet table
        posco_interaction = parse_posco(posco_result, metadata, i, sequence)
        posco_interactions.append(posco_interaction)

        # Only for the last frame perform extensive posco analysis and save
        if i == 0:
            cmd = f"po-sco {rec_path} {lig_path}  > {args.posco_interaction}"
            execute(cmd)

        # Clean up
        os.remove(rec_path)
        os.remove(lig_path)
        os.remove(posco_result)

    posco_interactions = pd.concat(posco_interactions)
    posco_interactions.to_parquet(args.posco_parquet)

if __name__ == '__main__':
    main()
