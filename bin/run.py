#!/usr/bin/env python

import subprocess
import os
import argparse

def call_my_binary(arg1, arg2, arg3):
    binary_path = '/home/pixelline/ownCloud/Institution/code/SqueezeMD/squeezemd/bin/7_interaction_csv.x'
    result = subprocess.run([binary_path, str(arg1), str(arg2), str(arg3)], shell=True, capture_output=True)
    return result.communicate()[0].decode('utf-8').strip()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--pdb', required=False, help='Position file (inpcrd)', default='output/demo/C1s_BD001/WT/amber/C1s_BD001.tleap.pdb')
    parser.add_argument('--resname', required=False, help='Configuration file with all required parameters (params.yml')
    parser.add_argument('--resid', required=False, help='Seed for inital velocities', default=23)

    args = parser.parse_args()
    x = call_my_binary(args.pdb, args.resname, args.resid)
    print(x.stdout.decode('utf-8'))
