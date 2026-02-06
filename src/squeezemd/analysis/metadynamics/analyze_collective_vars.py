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

    df = pd.read_csv(args.colvar, sep='\s+', comment="#",
                     names=['time', 'd1'])   # ['time', 'd1', 'c1']
    
    print("Loaded COLVAR with columns:", df.columns.tolist())
    print(df.head())

    #plt.subplot(2, 1, 1)
    # Plot CV vs time
    plt.plot(df['time'], df['d1'])
    plt.xlabel("Time (ps)")
    plt.ylabel("distance center of mass")
    plt.title("Collective Variable center of mass vs Time")
    plt.grid(True)

    """
    plt.subplot(2, 1, 2)
    # Plot CV vs time
    plt.plot(df['time'], df['c1'])
    plt.xlabel("Time (ps)")
    plt.ylabel("Contacts")
    plt.title("Collective Variable: Number of contacts vs Time")
    plt.grid(True)
    """

    plt.tight_layout()
    plt.savefig(args.colvar_fig)
    plt.close()

if __name__ == "__main__":
    main()
