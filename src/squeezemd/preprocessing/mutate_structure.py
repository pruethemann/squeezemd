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


def build_mutant_file_content(wt_sequence: str, mutation_string: str) -> str:
    """Build FoldX mutant-file content from a WT sequence and a mutation string.

    The mutation string is one or more single-point mutations joined by ``_``
    (e.g. ``R65E`` or ``R65E_Y117E``), each written as <WT-residue><resid><new
    residue>. Returns two lines: the wild-type sequence and the mutated sequence.

    Raises ``ValueError`` if a mutation cannot be parsed, its position is out of
    range, or its stated wild-type residue does not match the sequence.
    """
    if not mutation_string:
        raise ValueError("Empty mutation string")

    mutant_sequence = wt_sequence
    for mutation in mutation_string.split("_"):
        resname_wt = mutation[0]
        resname_mutated = mutation[-1]
        try:
            resid = int(mutation[1:-1])
        except ValueError as exc:
            raise ValueError(f"Cannot parse mutation '{mutation}' (expected e.g. R65E)") from exc

        if not 1 <= resid <= len(mutant_sequence):
            raise ValueError(f"Mutation '{mutation}' position {resid} is out of range 1..{len(mutant_sequence)}")

        if mutant_sequence[resid - 1] != resname_wt:
            raise ValueError(
                f"Wrong wild-type residue for mutation '{mutation}': sequence has "
                f"'{mutant_sequence[resid - 1]}' at position {resid}, not '{resname_wt}'"
            )

        residues = list(mutant_sequence)
        residues[resid - 1] = resname_mutated
        mutant_sequence = "".join(residues)

    return f"{wt_sequence}\n{mutant_sequence}"


def main():
    """
    Generate a mutation file required for FoldX mutagenesis.
    """
    args = parse_arguments()

    # Extract the WT ligand sequence (chain A) and build the WT/mutant file.
    wt_sequence = extract_ligand_sequence(args.ligand)
    content = build_mutant_file_content(wt_sequence, args.mutation)
    save_file(content, args.output)


if __name__ == "__main__":
    main()
