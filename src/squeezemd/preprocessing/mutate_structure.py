#!/usr/bin/env python

"""Generate FoldX mutation files from a ligand sequence.

This utility reads a ligand chain from a PDB, validates mutation strings
(e.g., R65E or R65E_Y117E), and writes the WT/mutant sequences in the
format expected by FoldX BuildModel.
"""

import argparse

from ..helper_functions import extract_ligand_sequence, save_file


def parse_arguments():
    """
    Parse CLI arguments for mutation generation.
    """
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument(
        "--ligand",
        required=True,
        help="PDB file containing the ligand at chain ID A. Required to extract ligand sequence",
    )
    parser.add_argument("--mutation", required=True, help="The mutation in the shape of R65E")

    # Output
    parser.add_argument(
        "--output", required=True, help="Mutation file which is required for a foldX mutagenesis process"
    )
    return parser.parse_args()


def main():
    """
    Generate a mutation file required for FoldX mutagenesis.
    """
    args = parse_arguments()

    # Extract ligand sequence and keep the WT sequence for output
    ligand_WT_sequence = extract_ligand_sequence(args.ligand)
    ligand_WT_sequence_original = ligand_WT_sequence

    # Mutations can be single (R65E) or multiple (R65E_Y117E)
    mutations = args.mutation.split("_")

    # Check every mutation
    for mutation in mutations:
        resname_WT = mutation[0]  # resname before mutation
        resname_mutated = mutation[-1]  # resname after mutation
        resid = int(mutation[1:-1])  # resid

        # Validate that the WT residue matches the sequence
        if ligand_WT_sequence[resid - 1] != resname_WT:
            raise Exception(
                f"You are mutating the wrong amino acid. AA before: {resname_WT} AA expected: {ligand_WT_sequence[resid - 1]} position: {resid}"
            )

        temp = list(ligand_WT_sequence)
        temp[resid - 1] = resname_mutated
        mut_seq = "".join(temp)

        print(resname_WT, " Mutate position ", resid, " with amino acid: ", resname_mutated)

        ligand_WT_sequence = mut_seq

    # Save WT sequence in line 1 and mutant sequence in line 2
    mut_seq = ligand_WT_sequence_original + "\n" + mut_seq
    save_file(mut_seq, args.output)


if __name__ == "__main__":
    main()
