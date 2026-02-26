from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


FES_COLUMN_NAMES = {
    3: ["d1", "F", "der_d1"],
    5: ["d1", "c1", "F", "der_d1", "der_c1"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate metadynamics free-energy/colvar files and create summary plots."
    )
    parser.add_argument(
        "root",
        type=Path,
        help="Project root containing metadynamics/C1s_Gigastasin/...",
    )
    parser.add_argument(
        "--project-id",
        default="P_22",
        help="Identifier used in output table (default: P_22).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for parquet and plot outputs.",
    )
    parser.add_argument(
        "--base-pattern",
        default="metadynamics/C1s_Gigastasin/**/**/metadynamics",
        help="Glob pattern (relative to root) where metadynamics files live.",
    )
    parser.add_argument(
        "--colvar-downsample",
        type=int,
        default=10,
        help="Read every N-th row from Colvar.dat (default: 10).",
    )
    parser.add_argument(
        "--skip-colvar",
        action="store_true",
        help="Skip loading and plotting Colvar.dat.",
    )
    return parser.parse_args()


def discover_files(root: Path, base_pattern: str, filename: str) -> list[Path]:
    pattern = Path(base_pattern) / filename
    files = sorted(root.glob(str(pattern)))
    return [file for file in files if file.is_file()]


def extract_metadata(file_path: Path) -> tuple[str, str]:
    parts = file_path.parts
    mutation = parts[-4]
    seed = parts[-3]
    return mutation, seed


def infer_fes_columns(file_path: Path) -> list[str]:
    with file_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            n_columns = len(stripped.split())
            if n_columns in FES_COLUMN_NAMES:
                return FES_COLUMN_NAMES[n_columns]
            return [f"col_{index + 1}" for index in range(n_columns)]
    raise ValueError(f"No data rows found in {file_path}")


def add_run_metadata(df: pd.DataFrame, mutation: str, seed: str, project_id: str) -> pd.DataFrame:
    result = df.copy()
    result["seed"] = seed
    result["mutation"] = mutation
    result["id"] = project_id
    result["sim"] = result["id"] + "_" + result["seed"] + "_" + result["mutation"]
    return result


def load_fes_data(fes_files: Iterable[Path], project_id: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for fes_file in fes_files:
        columns = infer_fes_columns(fes_file)
        df = pd.read_csv(fes_file, sep=r"\s+", comment="#", header=None, names=columns)
        mutation, seed = extract_metadata(fes_file)
        frames.append(add_run_metadata(df, mutation, seed, project_id))

    if not frames:
        raise FileNotFoundError("No fes.dat files found for the provided root/pattern.")
    return pd.concat(frames, ignore_index=True)


def load_colvar_data(colvar_files: Iterable[Path], project_id: str, downsample: int) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for colvar_file in colvar_files:
        df = pd.read_csv(
            colvar_file,
            sep=r"\s+",
            comment="#",
            header=None,
            names=["time", "d1", "c1"],
            skiprows=lambda idx: idx % downsample != 0,
        )
        mutation, seed = extract_metadata(colvar_file)
        frames.append(add_run_metadata(df, mutation, seed, project_id))

    if not frames:
        raise FileNotFoundError("No Colvar.dat files found for the provided root/pattern.")
    return pd.concat(frames, ignore_index=True)


def plot_energy_by_mutation(energy_df: pd.DataFrame, output_dir: Path) -> None:
    for mutation in sorted(energy_df["mutation"].dropna().unique()):
        subset = energy_df[energy_df["mutation"] == mutation]
        if subset.empty:
            continue

        plt.figure(figsize=(10, 6))
        for sim in sorted(subset["sim"].unique()):
            sns.lineplot(
                data=subset[subset["sim"] == sim],
                x="d1",
                y="F",
                legend=False,
            )
        plt.title(f"Free Energy vs CVs ({mutation})")
        plt.tight_layout()
        plt.savefig(output_dir / f"{mutation}.png", dpi=200)
        plt.close()


def plot_colvar_summary(colvar_df: pd.DataFrame, output_dir: Path) -> None:
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=colvar_df, x="time", y="c1", hue="mutation")
    plt.title("Colvar c1 over time")
    plt.tight_layout()
    plt.savefig(output_dir / "colvar_c1_vs_time.png", dpi=200)
    plt.close()


def main() -> None:
    args = parse_args()
    root = args.root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    fes_files = discover_files(root, args.base_pattern, "fes.dat")
    energy_df = load_fes_data(fes_files, args.project_id)
    energy_df.to_parquet(output_dir / "energy.parquet", index=False)
    plot_energy_by_mutation(energy_df, output_dir)

    print(f"Loaded {len(fes_files)} fes.dat files")
    print(f"Saved {output_dir / 'energy.parquet'}")

    if args.skip_colvar:
        return

    colvar_files = discover_files(root, args.base_pattern, "Colvar.dat")
    colvar_df = load_colvar_data(colvar_files, args.project_id, args.colvar_downsample)
    colvar_df.to_parquet(output_dir / "colvar.parquet", index=False)
    plot_colvar_summary(colvar_df, output_dir)

    print(f"Loaded {len(colvar_files)} Colvar.dat files")
    print(f"Saved {output_dir / 'colvar.parquet'}")


if __name__ == "__main__":
    main()

