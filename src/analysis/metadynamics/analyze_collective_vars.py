#!/usr/bin/env python3
"""
analyze_colvar.py

Analysis and plotting for a PLUMED COLVAR file that contains:
    time   d1

This script:
- Loads the COLVAR file
- Plots d1 vs time

Usage:
    python analyze_colvar.py COLVAR
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser()
    # input
    parser.add_argument("--colvar", help="Path to COLVAR file")

    # output
    parser.add_argument("--colvar_fig", help="Path to COLVAR file")

    return parser.parse_args()

def main():
    args = parse_arguments()

    df = pd.read_csv(args.colvar, sep='\s+', comment="#",
                     names=["time", "d1"])
    
    print("Loaded COLVAR with columns:", df.columns.tolist())
    print(df.head())

    # Plot CV vs time
    plt.figure()
    plt.plot(df["time"], df["d1"])
    plt.xlabel("Time (ps)")
    plt.ylabel("d1 (CV)")
    plt.title("Collective Variable d1 vs Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.colvar_fig)
    plt.close()
    #plt.show()


if __name__ == "__main__":
    main()
