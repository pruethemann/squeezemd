#!/usr/bin/env python
"""
Further infos:
https://blog.matteoferla.com/2020/07/filling-missing-loops-proper-way.html

python3 1_Mutate.py --csv config/simulations_tmp.csv --seed 2007 --name C1s_BD001_D18E_G36R

"""



import argparse
import MDAnalysis as mda
from collections import defaultdict



def save_file(content, output_file):
    f = open(output_file, "w")
    f.write(content)
    f.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--mutation',required=True, help='Configuration file with all required parameters')
    parser.add_argument('--output', required=True, help='Configuration file with all required parameters')

    args = parser.parse_args()


    WT_seq = "AKKKLPKCQKQEDCGSWDLKCNNVTKKCECRNQVCGRGCPKERYQRDKYGCRKCLCKGCDGFKCRLGCTYGFKTDKKGCEAFCTCNTKETACVNIWCTDPYKCNPESGRCEDPNEEYEYDYE"
    WT_original = WT_seq

    mutations = args.mutation.split('_')

    print(mutations)


    if not 'WT' in mutations:
        WT_seq = WT_original

        for mutation in mutations:
            aa_before = mutation[0]
            aa_after = mutation[-1]
            resid = int(mutation[1:-1])

            if WT_seq[resid-1] != aa_before:
                raise Exception(f"You are mutating the wrong amino acid. AA before: {aa_before} AA expected: {WT_seq[resid-1]} position: {resid}")


            temp = list(WT_seq)
            temp[resid-1] = aa_after
            mut_seq = "".join(temp)

            print(aa_before, " Mutate position ", resid, " with amino acid: ", aa_after)

            print(mut_seq)
            WT_seq = mut_seq

        print(mut_seq)
        mut_seq = WT_original + '\n' + mut_seq
        save_file(mut_seq, args.output)

    else:
        seq = WT_original + '\n' + WT_original
        save_file(seq, args.output)


def extract_ligand_sequence():
    """
    deprecated
    Extracts the amino acid sequence of the ligand. Chain ID must be I
    :return: amino acid sequence in 1 letter code
    """
    sequence = defaultdict(list)

    # Load PDB file
    u = mda.Universe("input.pdb")

    # Get protein atoms
    protein = u.select_atoms("protein")

    # Extract amino acid sequence and chain IDs
    seq = []
    chains = []
    for residue in protein.residues:

        resname = mda.lib.util.convert_aa_code(residue.resname)
        sequence[residue.segid].append(resname)


    for chain, seq in sequence.items():
        sequence[chain] = ''.join(seq)

    return sequence['I']
