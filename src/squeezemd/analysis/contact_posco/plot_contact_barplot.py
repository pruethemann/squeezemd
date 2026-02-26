#!/usr/bin/env python

"""Plot PoSCo interaction energy barplots by residue."""

import argparse
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
from os import path


INTERACTION_TYPES = ["total", "H-bond", "lipophilic", "Salt bridge"]
INTERACTION_STYLES = {
    "total": ("grey", "Total interaction energy"),
    "H-bond": ("dodgerblue", "H-bond interaction energy"),
    "lipophilic": ("darkorange", "Lipophilic interaction energy"),
    "Salt bridge": ("seagreen", "Salt bridge interaction energy"),
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",required=False,help="Define interaction input file, .parquet or .csv",default="results/posco/posco_interactions.parquet",)
    return parser.parse_args()


def load_data(input_path: str) -> pd.DataFrame:
    input_file = Path(input_path)
    df = pd.read_parquet(input_file)
    # Exclude water-mediated interactions # TODO: consider keeping these and adding a separate category for them
    return df[(df["receptor_resname"] != "HOH") & (df["ligand_resname"] != "HOH")].copy()


def select_interaction_type(df: pd.DataFrame, interaction_type: str) -> pd.DataFrame:
    if interaction_type == "total":
        return df
    if interaction_type == "H-bond":
        return df[(df["Interaction Type"] == "H-bond") & (df["Marked as Salt-Bridge"] == 0)]
    if interaction_type == "lipophilic":
        return df[df["Interaction Type"] == "lipophilic"]
    if interaction_type == "Salt bridge":
        return df[df["Marked as Salt-Bridge"] == 1]
    raise ValueError(f"Unsupported interaction type: {interaction_type}")


def aggregate_partner_energy(df: pd.DataFrame, interaction_partner: str, interaction_type: str) -> pd.DataFrame:
    resid_col = f"{interaction_partner}_resid"
    filtered = select_interaction_type(df, interaction_type)
    if filtered.empty:
        return pd.DataFrame(columns=[resid_col, "mean", "sd"])

    frame_count = max(filtered["frame"].nunique(), 1)
    # Note: don't group by mean beacuse we want to sum energies across residues for each seed before averaging across seeds, to avoid underestimating the energy of residues that have multiple interactions. This is a manual implementation of a groupby with nested aggregation to achieve this.
    per_seed = (filtered.groupby([resid_col, "seed"], as_index=False)["Energy (e)"].sum().assign(seed_energy=lambda x: x["Energy (e)"] / frame_count))
    # Now we can group by residue to get the mean and standard deviation across seeds / replicates
    out = (per_seed.groupby(resid_col, as_index=False)["seed_energy"].agg(mean="mean", sd="std").fillna({"sd": 0.0}).sort_values(resid_col))
    return out


def plot_partner(df: pd.DataFrame, interaction_partner: str, complex_name: str, mutation: str, output_file: Path) -> None:
    resid_col = f"{interaction_partner}_resid"
    fig, axes = plt.subplots(4, 1, figsize=(14, 22), sharex=True)

    for axis, interaction_type in zip(axes, INTERACTION_TYPES):
        data = aggregate_partner_energy(df, interaction_partner, interaction_type)

        # DEBUG
        #data.to_csv(f"debug_{complex_name}_{mutation}_{interaction_type}_{interaction_partner}.csv", index=False)

        color, label = INTERACTION_STYLES[interaction_type]
        if data.empty:
            axis.text(0.5, 0.5, "No data", ha="center", va="center", transform=axis.transAxes)
        else:
            axis.bar(data[resid_col], data["mean"], yerr=data["sd"], color=color)
        axis.axhline(y=0, color="black", linewidth=0.8)
        axis.set_ylabel("Energy (kcal/mol)")
        axis.set_title(f"{label} | {complex_name} | {mutation}")

    axes[-1].set_xlabel(f"{interaction_partner.capitalize()} residue")
    fig.tight_layout()
    #output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file)
    plt.close(fig)


def main() -> None:
    args = parse_arguments()
    df = load_data(args.input)

    # Loop through each complex/mutation group and generate plots for ligand and receptor interactions
    for (complex_name, mutation), group_df in df.groupby(["name", "mutation"], dropna=False):
        lig_output = path.join('results', 'posco', f'lig_barplot_{complex_name}_{mutation}.svg')
        rec_output = path.join('results', 'posco', f'rec_barplot_{complex_name}_{mutation}.svg')
        plot_partner(group_df, "ligand", str(complex_name), str(mutation), lig_output)
        plot_partner(group_df, "receptor", str(complex_name), str(mutation), rec_output)


if __name__ == "__main__":
    main()