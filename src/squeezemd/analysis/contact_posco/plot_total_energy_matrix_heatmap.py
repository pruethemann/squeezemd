#!/usr/bin/env python

"""Plot ligand–receptor total interaction energy matrix heatmaps for PoSCo output."""

import argparse
from pathlib import Path
from os import path

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        help="Define interaction input file, .parquet or .csv",
        default="results/posco/posco_interactions.parquet",
    )
    parser.add_argument(
        "--center",
        required=False,
        type=float,
        default=0.0,
        help="Center value for the diverging colormap.",
    )
    parser.add_argument(
        "--cmap",
        required=False,
        default="vlag",
        help="Seaborn/matplotlib colormap name.",
    )
    parser.add_argument(
        "--annot",
        action="store_true",
        help="Annotate cells with numeric values.",
    )
    return parser.parse_args()


def load_data(input_path: str) -> pd.DataFrame:
    input_file = Path(input_path)
    if input_file.suffix.lower() == ".csv":
        df = pd.read_csv(input_file)
    else:
        df = pd.read_parquet(input_file)

    # Exclude water-mediated interactions
    return df[(df["receptor_resname"] != "HOH") & (df["ligand_resname"] != "HOH")].copy()


def build_energy_matrix(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    frame_count = max(df["frame"].nunique(), 1)

    # Sum interaction energy per pair and seed first, then normalize by frame count
    per_seed_pair = (
        df.groupby(["ligand_resid", "receptor_resid", "seed"], as_index=False)["Energy (e)"]
        .sum()
        .assign(seed_energy=lambda x: x["Energy (e)"] / frame_count)
    )

    # Average pair interaction across seeds / replicates
    mean_pair = (
        per_seed_pair.groupby(["ligand_resid", "receptor_resid"], as_index=False)["seed_energy"]
        .mean()
        .rename(columns={"seed_energy": "mean_energy"})
    )

    matrix = mean_pair.pivot(index="ligand_resid", columns="receptor_resid", values="mean_energy")
    matrix = matrix.sort_index(axis=0).sort_index(axis=1)
    return matrix


def plot_energy_matrix(
    matrix: pd.DataFrame,
    complex_name: str,
    mutation: str,
    output_file: Path,
    center: float,
    cmap: str,
    annot: bool,
) -> None:
    sns.set_theme(style="white", context="talk")

    if matrix.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No interaction data", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    n_rows, n_cols = matrix.shape
    fig_w = max(8, min(0.45 * n_cols + 4, 26))
    fig_h = max(6, min(0.45 * n_rows + 4, 26))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        center=center,
        linewidths=0.3,
        linecolor="#f0f0f0",
        cbar_kws={"label": "Mean interaction energy (kcal/mol)"},
        square=False,
        annot=annot,
        fmt=".2f",
        annot_kws={"fontsize": 7} if annot else None,
    )

    ax.set_title(
        f"Ligand–receptor total interaction energy matrix | {complex_name} | {mutation}",
        fontsize=14,
        pad=14,
        weight="bold",
    )
    ax.set_xlabel("Receptor residue", labelpad=10)
    ax.set_ylabel("Ligand residue", labelpad=10)
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)

    fig.tight_layout(pad=1.2)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print("starting plot_total_energy_matrix_heatmap.py...")
    args = parse_arguments()
    df = load_data(args.input)

    print(df)

    for (complex_name, mutation), group_df in df.groupby(["name", "mutation"], dropna=False):
        matrix = build_energy_matrix(group_df)

        print(matrix)
        output = path.join("results", "posco", f"energy_matrix_heatmap_{complex_name}_{mutation}.svg")
        plot_energy_matrix(
            matrix=matrix,
            complex_name=str(complex_name),
            mutation=str(mutation),
            output_file=output,
            center=args.center,
            cmap=args.cmap,
            annot=args.annot,
        )


if __name__ == "__main__":
    main()
