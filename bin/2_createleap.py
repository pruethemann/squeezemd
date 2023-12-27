#!/usr/bin/env python
import argparse

def create_tleap(args):

    # Adapted free energy caluclation file

    with open('config/tleap_nosolvent.in', 'r') as f:
        content = f.read()


        content = content.replace("AMBERPDB", args.pdb)
        content = content.replace("PRMTOP", args.prmtop)
        content = content.replace("INPCRD", args.inpcrd)
        content = content.replace("PDBXXX", args.tleappdb)

    # Save Tleap conataing all file paths
    f = open(args.leap, "w")
    f.write(content)
    f.close()





if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input MD
    parser.add_argument('--pdb', required=False,help='', default='')

    # Output MD
    parser.add_argument('--prmtop', required=False,help='', default='')
    parser.add_argument('--inpcrd', required=False,help='', default='')
    parser.add_argument('--leap', required=False,help='', default='')
    parser.add_argument('--tleappdb', required=False,help='', default='')

    args = parser.parse_args()

    create_tleap(args)
