#!/usr/bin/env python

"""Plot RMSF distributions from a parquet table."""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_arguments():
    """Parse CLI arguments for RMSF plotting."""
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("--input", required=False)

    # Output
    parser.add_argument("--output", required=False, default="rmsf.svg", help="")

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Load RMSF data and plot with standard deviation shading
    rmsf_df = pd.read_parquet(args.input)

    sns.lineplot(data=rmsf_df, x="resid", y="rmsf", errorbar="sd")

    plt.savefig(args.output)
    plt.close()


if __name__ == "__main__":
    main()
