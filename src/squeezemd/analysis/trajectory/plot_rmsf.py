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

    # Load RMSF data and plot with standard deviation shading.
    rmsf_df = pd.read_parquet(args.input)

    # Colour by mutation when the metadata is present and distinguishes curves,
    # so different variants are visually separable rather than averaged together.
    hue = "mutation" if "mutation" in rmsf_df.columns and rmsf_df["mutation"].nunique() > 1 else None

    sns.lineplot(data=rmsf_df, x="resid", y="rmsf", hue=hue, errorbar="sd")
    plt.xlabel("Residue")
    plt.ylabel("RMSF (Å)")
    plt.tight_layout()

    plt.savefig(args.output)
    plt.close()


if __name__ == "__main__":
    main()
