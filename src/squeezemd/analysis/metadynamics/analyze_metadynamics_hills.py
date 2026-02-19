#!/usr/bin/env python3
"""Plot a 1D free-energy surface from PLUMED sum_hills output.

Reads a `fes.dat` file produced by `plumed sum_hills` and generates
basic FES plots. The helper includes an optional truncation routine
for convergence checks.
"""

import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt


def plot_fes(df, label=None):
    """Plot a single FES curve (d1 vs free energy)."""
    plt.plot(df["d1"], df["F"], label=label)
    plt.xlabel("d1 (collective variable)")
    plt.ylabel("Free Energy")
    plt.title("Free Energy Surface")
    plt.grid(True)


# -----------------------------
# Truncate HILLS for convergence
# -----------------------------
def truncate_hills(hills_path, fraction, outfile):
    """Write a truncated HILLS file for convergence diagnostics."""
    with open(hills_path) as f:
        lines = f.readlines()

    header = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#")]

    n = max(1, int(len(data) * fraction))
    with open(outfile, "w") as out:
        out.writelines(header)
        out.writelines(data[:n])


def get_fes_from_hills(hills_path, outfile):
    """Run `plumed sum_hills` to generate a FES file."""
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

    df = pd.read_csv(args.fes, sep=r'\s+', comment="#", header=None)
    n_cols = df.shape[1]
    if n_cols >= 5:
        df = df.iloc[:, :5]
        df.columns = ['d1', 'c1', 'F', 'der_d1', 'der_c1']
    elif n_cols == 3:
        df.columns = ['d1', 'F', 'der_d1']
    else:
        raise ValueError(f"Unexpected FES format with {n_cols} columns in {args.fes}")

    """
    d1: first CV: COM distance in nm
    c1: second CV: interace contactes (the coordination number, wo unit)
    F: Free energy (F(d1,c1))in k//Mol
    der_d1: Gradient (force) of the free energy along the COM distance
    der_c1: gradient (force) along the contacts CV

    Use gradients to find transition states
    """
    if 'c1' in df.columns:
        heatmap = df.pivot_table(index='c1', columns='d1', values='F', aggfunc='mean')
        sns.heatmap(heatmap.sort_index().sort_index(axis=1), cmap='viridis')
        plt.xlabel('d1 (COM distance)')
        plt.ylabel('c1 (contacts)')
        plt.title('Free Energy Surface (2D)')
    else:
        sns.lineplot(data=df,
                     x='d1',
                     y='F')
        plt.xlabel('d1 (collective variable)')
        plt.ylabel('Free Energy')
        plt.title('Free Energy Surface (1D)')
    
    plt.savefig(args.freeenergy)
    plt.close()


if __name__ == "__main__":
    main()
