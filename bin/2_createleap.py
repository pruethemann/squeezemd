#!/usr/bin/env python
import os
import sys
import argparse
from Helper import save_file


def create_tleap(args):
    """
    Generates a tleap.in file which is required for tleap
    :param args:
    :return:
    """

    # Opens and modifies default tleap.in file
    tpleap_default = os.path.join(sys.prefix, 'tleap_nosolvent.in')

    with open(tpleap_default, 'r') as f:
        content = f.read()

        content = content.replace("AMBERPDB", args.pdb)
        content = content.replace("PRMTOP", args.prmtop)
        content = content.replace("INPCRD", args.inpcrd)
        content = content.replace("PDBLEAP", args.tleappdb)

    # Save Tleap conataing all file paths
    save_file(content, args.leap)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input MD
    parser.add_argument('--pdb', required=False, help='Path the protein pdb file', default='protein.pdb')

    # Output MD
    parser.add_argument('--prmtop', required=False, help='Location of Amber topology file')
    parser.add_argument('--inpcrd', required=False,help='Location of Amber topology file')
    parser.add_argument('--leap', required=False,help='Location of tleap file')
    parser.add_argument('--tleappdb', required=False, help='Amber modified pdb file')

    args = parser.parse_args()

    create_tleap(args)
