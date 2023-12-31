#!/usr/bin/env python

"""
    Performs site directed mutatagensis of one or multiple mutations in a protein.
    This script prepares the files for foldX which is doing the local minimsiaotn.

    Further infos:
    https://blog.matteoferla.com/2020/07/filling-missing-loops-proper-way.html

    There used to be a function called extract_ligand_sequence which tried to extract the sequences directly from the pdb
"""

import argparse
from Helper import save_file

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--mutation',required=True, help='Configuration file with all required parameters')
    parser.add_argument('--output', required=True, help='Configuration file with all required parameters')

    args = parser.parse_args()

    # TODO: Decide where to extract that from
    WT_seq = "AKKKLPKCQKQEDCGSWDLKCNNVTKKCECRNQVCGRGCPKERYQRDKYGCRKCLCKGCDGFKCRLGCTYGFKTDKKGCEAFCTCNTKETACVNIWCTDPYKCNPESGRCEDPNEEYEYDYE"
    WT_original = WT_seq

    # Get all mutations
    mutations = args.mutation.split('_')

    # Only execute if mutagensis necessary
    if not 'WT' in mutations:
        WT_seq = WT_original

        for mutation in mutations:
            aa_before = mutation[0]
            aa_after = mutation[-1]
            resid = int(mutation[1:-1])

            # Checks whether the orginal resname is correct
            if WT_seq[resid-1] != aa_before:
                raise Exception(f"You are mutating the wrong amino acid. AA before: {aa_before} AA expected: {WT_seq[resid-1]} position: {resid}")


            temp = list(WT_seq)
            temp[resid-1] = aa_after
            mut_seq = "".join(temp)

            print(aa_before, " Mutate position ", resid, " with amino acid: ", aa_after)

            WT_seq = mut_seq

        mut_seq = WT_original + '\n' + mut_seq
        save_file(mut_seq, args.output)

    else: # Save Wildtype sequence again to make sure Snakemake is happy
        seq = WT_original + '\n' + WT_original
        save_file(seq, args.output)
