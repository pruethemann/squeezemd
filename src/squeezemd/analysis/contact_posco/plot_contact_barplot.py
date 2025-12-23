#!/usr/bin/env python

import argparse
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
from squeezemd import io

INTERACTION_TYPES = ["total", "H-bond", "lipophilic", "salt-bridge"]

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Per-residue interaction barplots (mean ± SD) for ligand and receptor using the squeezeMD unified tables.")
    parser.add_argument("--complex",required=True,help="Name of the complex (must match 'complex' column in tables).")
    parser.add_argument("--condition",default=None,help="Optional condition filter (e.g. 'NPT_300K', 'WT_meta') must match 'condition' in runs table if provided.")
    parser.add_argument("--output-dir",default="figures/posco_barplots",help="Directory to store output SVGs.")
    parser.add_argument("--dpi",type=int,default=300,help="Figure DPI.")
    return parser.parse_args()

def get_runs(complex_name: str, condition: str | None) -> pd.DataFrame:
    filters = {"complex": complex_name}
    runs = io.load_runs(filters)

    if condition is not None:
        runs = runs[runs["condition"] == condition]

    if runs.empty:
        raise ValueError(f"No runs found for complex='{complex_name}' "
                         f"and condition='{condition}'.")
    return runs


def load_residue_interactions_for_runs(runs: pd.DataFrame) -> pd.DataFrame:
    """
    Load residue_interactions only for the desired runs.
    """
    run_ids = runs["run_id"].unique().tolist()
    ri = io.load_residue_interactions({"run_id": run_ids})
    # Optional safety: keep only 'energy' metrics
    ri = ri[ri["metric"] == "energy"]
    return ri


def load_residues_for_runs(runs: pd.DataFrame) -> pd.DataFrame:
    """
    Load residues table and collapse across runs so we get a single
    sequence definition per (complex, mutation, partner, resid).
    """
    run_ids = runs["run_id"].unique().tolist()
    res = io.load_residues({"run_id": run_ids})

    # If residues are identical for all seeds of a mutation, this collapse is safe
    # Group by fields that define "a residue" and take the first row.
    res = (
        res.groupby(
            [
                "run_id",
                "partner",
                "chain_id",
                "resid",
                "resname",
                "protein_label",
                "sequence_index",
            ],
            as_index=False,
        )
        .first()
    )

    # We will later filter by mutation by joining via run_id
    return res


def aggregate_interactions(
    ri: pd.DataFrame,
    runs: pd.DataFrame,
    partner: str,
    interaction_type: str,
) -> pd.DataFrame:
    """
    Aggregate interactions for one partner and one interaction_type.

    Returns a DataFrame with columns:
        complex, mutation, partner, chain_id, resid, resname,
        mean, sd
    """

    # 1) Filter to partner + interaction_type
    df = ri[
        (ri["partner"] == partner)
        & (ri["interaction_type"] == interaction_type)
        & (ri["metric"] == "energy")
    ].copy()
    if df.empty:
        raise ValueError(
            f"No residue_interactions rows for partner='{partner}', "
            f"interaction_type='{interaction_type}'."
        )

    # 2) Map run_id -> mutation from runs table
    df = df.merge(
        runs[["run_id", "complex", "mutation"]],
        on="run_id",
        how="left",
        validate="many_to_one",
    )

    # 3) First average over frames per run
    per_run = (
        df.groupby(
            [
                "complex",
                "mutation",
                "run_id",
                "partner",
                "chain_id",
                "resid",
                "resname",
            ],
            as_index=False,
        )["value"]
        .mean()
        .rename(columns={"value": "mean_per_run"})
    )

    # 4) Then average over runs (seeds) per mutation
    per_residue = (
        per_run.groupby(
            [
                "complex",
                "mutation",
                "partner",
                "chain_id",
                "resid",
                "resname",
            ],
            as_index=False,
        )["mean_per_run"]
        .agg(["mean", "std"])
        .reset_index()
    )

    per_residue = per_residue.rename(
        columns={"mean": "mean", "std": "sd"}
    )

    return per_residue


def merge_with_sequence(
    per_residue: pd.DataFrame,
    residues: pd.DataFrame,
    runs: pd.DataFrame,
    partner: str,
    mutation: str,
) -> pd.DataFrame:
    """
    For a given mutation and partner, merge aggregated per-residue energies
    with the residue metadata (sequence_index etc.), and fill missing residues
    with 0 mean (no interaction) and 0 SD.
    """
    # Subset residues to runs with this mutation & partner
    mut_run_ids = runs.loc[runs["mutation"] == mutation, "run_id"].unique()
    res_sub = residues[
        (residues["run_id"].isin(mut_run_ids))
        & (residues["partner"] == partner)
    ].copy()

    # Collapse across run_id for that mutation: we just want one row per residue
    res_sub = (
        res_sub.groupby(
            ["partner", "chain_id", "resid", "resname", "protein_label", "sequence_index"],
            as_index=False,
        )
        .first()
    )

    # Merge energies
    merged = res_sub.merge(
        per_residue[
            [
                "complex",
                "mutation",
                "partner",
                "chain_id",
                "resid",
                "resname",
                "mean",
                "sd",
            ]
        ],
        on=["partner", "chain_id", "resid", "resname"],
        how="left",
    )

    # For residues with no data, set 0 and 0
    merged["mean"] = merged["mean"].fillna(0.0)
    merged["sd"] = merged["sd"].fillna(0.0)

    # Sort by sequence_index for plotting
    merged = merged.sort_values("sequence_index").reset_index(drop=True)
    return merged


def plot_partner_for_mutation(
    complex_name: str,
    mutation: str,
    partner: str,
    per_type_data: dict[str, pd.DataFrame],
    out_path: Path,
    dpi: int,
) -> None:
    """
    per_type_data: mapping interaction_type -> merged dataframe (with sequence_index, mean, sd)
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        nrows=len(INTERACTION_TYPES), ncols=1, figsize=(15, 3 * len(INTERACTION_TYPES)), sharex=True
    )

    if len(INTERACTION_TYPES) == 1:
        axes = [axes]

    for ax, interaction_type in zip(axes, INTERACTION_TYPES):
        df = per_type_data[interaction_type]
        x = df["sequence_index"]
        y = df["mean"]
        yerr = df["sd"]

        ax.bar(x, y, yerr=yerr)
        ax.set_ylabel(f"{interaction_type} energy\n(kcal/mol)")
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_xlim(x.min() - 1, x.max() + 1)

    axes[-1].set_xlabel("Residue index (sequence_index)")

    fig.suptitle(
        f"{complex_name} – {mutation} – {partner}",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_arguments()

    complex_name = args.complex
    condition = args.condition
    out_dir = Path(args.output_dir)

    runs = get_runs(complex_name, condition)
    ri = load_residue_interactions_for_runs(runs)
    residues = load_residues_for_runs(runs)

    mutations = sorted(runs["mutation"].unique())

    for mutation in mutations:
        for partner in ["ligand", "receptor"]:
            per_type_data: dict[str, pd.DataFrame] = {}

            for interaction_type in INTERACTION_TYPES:
                per_residue = aggregate_interactions(
                    ri=ri,
                    runs=runs,
                    partner=partner,
                    interaction_type=interaction_type,
                )
                merged = merge_with_sequence(
                    per_residue=per_residue,
                    residues=residues,
                    runs=runs,
                    partner=partner,
                    mutation=mutation,
                )
                per_type_data[interaction_type] = merged

            out_path = (
                out_dir
                / f"{complex_name}__{mutation}__{partner}_residue_interactions.svg"
            )
            plot_partner_for_mutation(
                complex_name=complex_name,
                mutation=mutation,
                partner=partner,
                per_type_data=per_type_data,
                out_path=out_path,
                dpi=args.dpi,
            )


if __name__ == "__main__":
    main()
