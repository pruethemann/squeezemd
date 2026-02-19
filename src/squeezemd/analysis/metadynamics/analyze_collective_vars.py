#!/usr/bin/env python3
"""Analyze PLUMED COLVAR output (time vs collective variable).

Loads a COLVAR file and plots the primary CV (d1) against time.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():
    """Parse CLI arguments for COLVAR plotting."""
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument("--colvar", help="Path to COLVAR file")

    # Output
    parser.add_argument("--colvar_fig", help="Path to COLVAR file")

    return parser.parse_args()

def main():
    args = parse_arguments()

    df = pd.read_csv(args.colvar, sep=r'\s+', comment="#", header=None)
    if df.shape[1] >= 3:
        df = df.iloc[:, :3]
        df.columns = ['time', 'd1', 'c1']
    elif df.shape[1] == 2:
        df.columns = ['time', 'd1']
    else:
        raise ValueError(f"Unexpected COLVAR format with {df.shape[1]} columns in {args.colvar}")
    
    print("Loaded COLVAR with columns:", df.columns.tolist())
    print(df.head())

    if 'c1' in df.columns:
        plt.subplot(2, 1, 1)
        plt.plot(df['time'], df['d1'])
        plt.xlabel("Time (ps)")
        plt.ylabel("distance center of mass")
        plt.title("Collective Variable center of mass vs Time")
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.plot(df['time'], df['c1'])
        plt.xlabel("Time (ps)")
        plt.ylabel("Contacts")
        plt.title("Collective Variable: Number of contacts vs Time")
        plt.grid(True)
    else:
        plt.plot(df['time'], df['d1'])
        plt.xlabel("Time (ps)")
        plt.ylabel("distance center of mass")
        plt.title("Collective Variable center of mass vs Time")
        plt.grid(True)

    plt.tight_layout()
    plt.savefig(args.colvar_fig)
    plt.close()

if __name__ == "__main__":
    main()
