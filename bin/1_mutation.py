#!/usr/bin/env python

"""

"""

import argparse
import os
import MDAnalysis as mda
from Helper import save_file

def parse_arguments():
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand',required=True, help='PDB file containing the ligand at chain ID I. Required to extract ligand sequence')
    parser.add_argument('--mutation',required=True, help='')
    parser.add_argument('--output', required=True, help='')
    return parser.parse_args()


def extract_ligand_sequence(pdb_ligand: os.path):
    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # Extract ligand at chain I
    ligand = u.select_atoms('segid I')

    # Return sequence
    return str(ligand.residues.sequence().seq)

if __name__ == '__main__':

    args = parse_arguments()

    # Extract ligand sequence and copy for later
    ligand_WT_sequence = extract_ligand_sequence(args.ligand)
    ligand_WT_sequence_original = ligand_WT_sequence

    #TODO DEBUG
    print(ligand_WT_sequence)

    # Get all mutations which are separated by underscore
    mutations = args.mutation.split('_')

    # Check every mutation
    for mutation in mutations:
        resname_WT = mutation[0]
        resname_mutated = mutation[-1]
        resid = int(mutation[1:-1])

        # Checks whether the orginal resname is correct
        if ligand_WT_sequence[resid-1] != resname_WT:
            raise Exception(f"You are mutating the wrong amino acid. AA before: {resname_WT} AA expected: {ligand_WT_sequence[resid-1]} position: {resid}")

        temp = list(ligand_WT_sequence)
        temp[resid-1] = resname_mutated
        mut_seq = "".join(temp)

        print(resname_WT, " Mutate position ", resid, " with amino acid: ", resname_mutated)

        ligand_WT_sequence = mut_seq

    # Save the orginal sequence in line 1 and the mutated sequence in line 2
    mut_seq = ligand_WT_sequence_original + '\n' + mut_seq
    save_file(mut_seq, args.output)
