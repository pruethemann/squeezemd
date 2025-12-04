#!/usr/bin/env python3
"""
analyze_fes.py

Plot a 1D free-energy surface from PLUMED sum_hills output
and optionally perform convergence analysis by truncating the HILLS file.
"""

import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt


def plot_fes(df, label=None):
    plt.plot(df["d1"], df["F"], label=label)
    plt.xlabel("d1 (collective variable)")
    plt.ylabel("Free Energy")
    plt.title("Free Energy Surface")
    plt.grid(True)


# -----------------------------
# Truncate HILLS for convergence
# -----------------------------
def truncate_hills(hills_path, fraction, outfile):
    with open(hills_path) as f:
        lines = f.readlines()

    header = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#")]

    n = max(1, int(len(data) * fraction))
    with open(outfile, "w") as out:
        out.writelines(header)
        out.writelines(data[:n])


def get_fes_from_hills(hills_path, outfile):
    cmd = ["plumed", "sum_hills", "--hills", hills_path, "--outfile", outfile, "--mintozero"]
    subprocess.run(cmd, check=True)

def parse_args():
    parser = argparse.ArgumentParser()

    # input
    parser.add_argument("--fes", help="Path to fes.dat", required=False)

    # output
    parser.add_argument("--freeenergy", help="Path to HILLS file", required=False)

    # parameters
    parser.add_argument("--fractions", nargs="+", type=float,
                        default=[0.25, 0.50, 0.75, 1.00])
    
    return parser.parse_args()

import seaborn as sns

# -----------------------------
# Main CLI
# -----------------------------
def main():
    args = parse_args()

    # Read FES
    df = pd.read_csv(args.fes, delim_whitespace=True, comment="#", header=None, names=['d1', 'F' ,'der_d1'])

    sns.lineplot(data=df,
                 x='d1',
                 y='F')
    
    plt.savefig(args.freeenergy)
    plt.close()
    #plt.show()


if __name__ == "__main__":
    main()
