#!/usr/bin/env python

"""
    Currently depracted. Used to split chains for Free Energy calculations

"""
import argparse
import MDAnalysis

def create_tleap_MMGBSA(args):

    # Adapted free energy caluclation file

    with open('config/tleap_MMGBSA.in', 'r') as f:
        content = f.read()

        content = content.replace("COMPLEX", args.complex)
        content = content.replace("LIGAND", args.ligand)
        content = content.replace("RECEPTOR", args.receptor)

        content = content.replace("LIG_PRMTOP", args.lig_prmtop)
        content = content.replace("REC_PRMTOP", args.rec_prmtop)
        content = content.replace("COM_PRMTOP", args.com_prmtop)

    # Save Tleap conataing all file paths
    f = open(args.leap, "w")
    f.write(content)
    f.close()

def split_chains(args):
    """
    Splits a complex into a receptor and a ligand pdb

    :param args:
    :return:
    """

    # Import the pdb file as universe
    u = MDAnalysis.Universe(args.complex)

    # select chain A as receptor
    # TODO: allow multiple chains for receptor AND Ligands
    receptor = u.select_atoms(f"chainID {args.chain_receptor} or chainID B")

    # select chain I as receptor
    ligand = u.select_atoms(f"chainID {args.chain_ligand} or chainID H or chainID I")

    # Save pdb files
    receptor.write(args.receptor)
    ligand.write(args.ligand)


if __name__ == '__main__':
    # Parse Arguments
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--complex', required=True, help='PDB file prepared with pdb4abmer. Ligand + Receptor = Complex')

    # Output
    parser.add_argument('--ligand', required=True, help='PDB file only of Ligand')
    parser.add_argument('--receptor', required=True, help='PDB file only of receptor')

    # Parameters: TODO: Stay flexible
    parser.add_argument('--chain_ligand', required=False, help='', default='I')
    parser.add_argument('--chain_receptor', required=False, help='', default='A')

    # Output MMGBSA
    parser.add_argument('--lig_prmtop', required=False,help='', default='')
    parser.add_argument('--rec_prmtop', required=False,help='', default='')
    parser.add_argument('--com_prmtop', required=False,help='', default='')
    parser.add_argument('--leap', required=False,help='', default='')

    args = parser.parse_args()

    split_chains(args)
    create_tleap_MMGBSA(args)
