#!/usr/bin/env python

import argparse
from Helper import save_file, extract_ligand_sequence

def parse_arguments():
    """
    Parse Arguments
    """
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument('--ligand',required=True, help='PDB file containing the ligand at chain ID A. Required to extract ligand sequence')
    parser.add_argument('--mutation',required=True, help='The mutation in the shape of R65E')

    # Output
    parser.add_argument('--output', required=True, help='Mutation file which is required for a foldX mutagenesis process')
    return parser.parse_args()


def main():
    """
    Generate a mutation file required for foldX mutagensis
    """
    args = parse_arguments()

    # Extract ligand sequence and copy for later
    ligand_WT_sequence = extract_ligand_sequence(args.ligand)
    ligand_WT_sequence_original = ligand_WT_sequence

    # Get all mutations which are separated by underscore
    # tools handles single mutations (ex. R65E) and multiple mutations (ex R65E_Y117E)
    mutations = args.mutation.split('_')

    # Check every mutation
    for mutation in mutations:
        resname_WT = mutation[0]            # resname before mutation
        resname_mutated = mutation[-1]      # resname after mutation
        resid = int(mutation[1:-1])         # resid

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

if __name__ == '__main__':
    main()
