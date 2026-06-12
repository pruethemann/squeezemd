#!/usr/bin/env python

"""Plot PoSCo interaction-energy heatmaps (residue x frame) per interaction type."""

import argparse

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

INTERACTION_TYPES = ["total", "H-bond", "lipophilic", "Salt bridge"]
INTERACTION_STYLE = {
    "total": ("Greys", "Total interaction energy (kcal/mol)"),
    "H-bond": ("Blues", "H-bond energy (kcal/mol)"),
    "lipophilic": ("Oranges", "Hydrophobic interaction energy (kcal/mol)"),
    "Salt bridge": ("Greens", "Salt-bridge interaction energy (kcal/mol)"),
}


def parse_arguments():
    """Parse CLI arguments for heatmap generation."""
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        help="Interaction input file, .parquet",
        default="results/posco/posco_interactions.parquet",
    )

    # Output
    parser.add_argument(
        "-l",
        "--ligand_interaction",
        required=False,
        help="Ligand heatmap output file, .svg",
        default="lig_heatmap.svg",
    )
    parser.add_argument(
        "-r",
        "--receptor_interaction",
        required=False,
        help="Receptor heatmap output file, .svg",
        default="rec_heatmap.svg",
    )

    return parser.parse_args()


def select_interaction_type(df: pd.DataFrame, interaction_type: str) -> pd.DataFrame:
    """Filter the interaction table down to a single interaction category."""
    if interaction_type == "H-bond":
        return df[(df["Interaction Type"] == "H-bond") & (df["Marked as Salt-Bridge"] == 0)]
    if interaction_type == "lipophilic":
        return df[df["Interaction Type"] == "lipophilic"]
    if interaction_type == "Salt bridge":
        return df[df["Marked as Salt-Bridge"] == 1]
    return df  # "total"


def interaction_data_aggregation(df_filtered, interaction_partner, interaction_type):
    """Filter, aggregate (seed-averaged per frame) and pivot to residue x frame."""
    df_interaction = select_interaction_type(df_filtered, interaction_type)

    resid = f"{interaction_partner}_resid"
    resname = f"{interaction_partner}_resname"

    n_seeds = max(df_interaction["seed"].nunique(), 1)

    seed_avg = df_interaction.groupby([resid, resname, "frame"])["Energy (e)"].sum().reset_index()
    seed_avg["Energy (e)"] = seed_avg["Energy (e)"].div(n_seeds)
    seed_avg["residue_labels"] = seed_avg[resname] + " " + seed_avg[resid].astype(str)

    # Colorbar limit from the strongest (most negative) energy.
    emax = seed_avg["Energy (e)"].min() * -1 if not seed_avg.empty else 1.0

    heatmap_data = pd.pivot_table(seed_avg, index=[resid, "residue_labels"], columns="frame", values="Energy (e)")
    heatmap_data = heatmap_data.sort_index(level=resid)
    return heatmap_data, emax


def plot_interactions(heatmap_data, emax, interaction_type):
    """Render a single interaction-type heatmap panel."""
    interaction_cmap, interaction_label = INTERACTION_STYLE[interaction_type]

    if heatmap_data.empty:
        plt.text(0.5, 0.5, f"No {interaction_type} data", ha="center", va="center")
        return

    sorted_labels = [label for (_residue, label) in heatmap_data.index]

    ax = sns.heatmap(
        heatmap_data * -1,
        cmap=interaction_cmap,
        vmin=0,
        vmax=emax,
        yticklabels=sorted_labels,
        cbar_kws={"label": interaction_label, "pad": 0.012},
    )
    ax.set_title(f"Interaction type: {interaction_type.capitalize()}", fontsize=12)
    ax.set_xlabel("Frame number", fontsize=12)
    ax.set_ylabel("Residues", fontsize=12)

    cbar = ax.collections[0].colorbar
    cbar.set_label(interaction_label, rotation=90, labelpad=10)


def main():
    args = parse_arguments()

    # 1. Data import
    try:
        df = pd.read_parquet(args.input)
    except Exception as exc:
        raise RuntimeError(f"Failed to read interaction file '{args.input}'") from exc

    # 2. Exclude water-mediated interactions.
    df_filtered = df[(df["receptor_resname"] != "HOH") & (df["ligand_resname"] != "HOH")]

    # 3. One figure per interaction partner, one subplot per interaction type.
    figure_files = [args.ligand_interaction, args.receptor_interaction]
    for fig_file, interaction_partner in zip(figure_files, ["ligand", "receptor"], strict=True):
        fig = plt.figure(figsize=(15, 30))
        for i, interaction_type in enumerate(INTERACTION_TYPES):
            heatmap_data, energy_max = interaction_data_aggregation(df_filtered, interaction_partner, interaction_type)
            plt.subplot(len(INTERACTION_TYPES), 1, i + 1)
            plot_interactions(heatmap_data, energy_max, interaction_type)

        plt.tight_layout()
        fig.savefig(fig_file)
        plt.close(fig)


if __name__ == "__main__":
    main()
