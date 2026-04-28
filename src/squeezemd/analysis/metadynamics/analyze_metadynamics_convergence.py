#!/usr/bin/env python3
"""Metadynamics convergence analysis.

Checks convergence by computing the free-energy surface (FES) at multiple
time fractions using `plumed sum_hills --stride` and plotting how ΔF (the
free-energy difference between the bound minimum and the unbound plateau)
evolves over simulation time.

A converged metadynamics run should show ΔF stabilising as time increases.
"""

import argparse
import os
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_fes(path: str) -> pd.DataFrame:
    """Parse a PLUMED fes.dat file into a DataFrame."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
    n = df.shape[1]
    if n >= 5:
        df = df.iloc[:, :5]
        df.columns = ["d1", "c1", "F", "der_d1", "der_c1"]
    elif n == 3:
        df.columns = ["d1", "F", "der_d1"]
    else:
        raise ValueError(f"Unexpected FES format with {n} columns in {path}")
    # Replace PLUMED sentinel (1e308) with NaN so they don't distort stats
    df["F"] = df["F"].replace(1e308, np.nan)
    df["F"] = df["F"].replace(-1e308, np.nan)
    return df


def _delta_f(df: pd.DataFrame, bound_cutoff_nm: float = 0.5) -> float:
    """Return ΔF = F(unbound plateau) − F(bound minimum) in kJ/mol.

    The bound state is defined as d1 < *bound_cutoff_nm*; the unbound state
    is d1 > max(d1) − bound_cutoff_nm.  Both are taken as the minimum F in
    their respective windows to be robust against noise.
    """
    d = df.dropna(subset=["F"])
    if d.empty:
        return np.nan

    d_max = d["d1"].max()
    bound = d[d["d1"] < bound_cutoff_nm]
    unbound = d[d["d1"] > (d_max - bound_cutoff_nm)]

    if bound.empty or unbound.empty:
        # Fall back to global range as proxy
        return float(d["F"].max() - d["F"].min())

    return float(unbound["F"].min() - bound["F"].min())


def _count_hills(hills_path: str) -> int:
    """Count deposited hills (non-comment lines) in a HILLS file."""
    count = 0
    with open(hills_path) as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                count += 1
    return count


def compute_convergence(hills_path: str, fractions: list[float], bound_cutoff_nm: float) -> pd.DataFrame:
    """Run sum_hills at each fraction and return a DataFrame with ΔF vs fraction."""
    n_hills = _count_hills(hills_path)
    records = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for frac in fractions:
            stride = max(1, int(n_hills * frac))
            fes_path = os.path.join(tmpdir, f"fes_{frac:.2f}.dat")
            cmd = [
                "plumed", "sum_hills",
                "--hills", hills_path,
                "--outfile", fes_path,
                "--stride", str(stride),
                "--mintozero",
                "--kt", "2.479"
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"WARNING: sum_hills failed for fraction {frac}: {result.stderr}")
                continue

            # sum_hills with --stride writes multiple fes files: fes_<n>.dat
            # Pick the last one produced (= fullest coverage up to stride*n)
            produced = sorted(
                [f for f in os.listdir(tmpdir) if f.startswith("fes_") and f.endswith(".dat")],
                key=lambda x: int(x.replace("fes_", "").replace(".dat", "")) if x.replace("fes_", "").replace(".dat", "").isdigit() else 0,
            )
            # The file we requested directly (non-strided run) or the last strided one
            target = fes_path if os.path.exists(fes_path) else (
                os.path.join(tmpdir, produced[-1]) if produced else None
            )
            if target is None or not os.path.exists(target):
                continue

            try:
                df = _read_fes(target)
                dF = _delta_f(df, bound_cutoff_nm=bound_cutoff_nm)
                records.append({"fraction": frac, "hills": stride, "delta_F_kJmol": dF})
            except Exception as e:
                print(f"WARNING: could not parse FES for fraction {frac}: {e}")

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_convergence(conv_df: pd.DataFrame, output_path: str) -> None:
    """Plot ΔF vs fraction of simulation time."""
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(conv_df["fraction"] * 100, conv_df["delta_F_kJmol"],
            marker="o", linewidth=2, color="steelblue")
    ax.axhline(conv_df["delta_F_kJmol"].iloc[-1], linestyle="--",
               color="gray", linewidth=1, label="Final ΔF")

    ax.set_xlabel("Simulation progress (%)")
    ax.set_ylabel("ΔF (kJ/mol)")
    ax.set_title("Metadynamics convergence: ΔF vs time")
    ax.legend()
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assess metadynamics convergence by computing FES at multiple time fractions."
    )
    parser.add_argument("--hills", required=True, help="Path to PLUMED HILLS file.")
    parser.add_argument("--convergence", required=True, help="Output PNG path for convergence plot.")
    parser.add_argument(
        "--fractions", nargs="+", type=float,
        default=[0.25, 0.50, 0.75, 1.00],
        help="Fractions of simulation to evaluate (default: 0.25 0.50 0.75 1.00).",
    )
    parser.add_argument(
        "--bound_cutoff_nm", type=float, default=0.5,
        help="d1 threshold (nm) separating bound from unbound state (default: 0.5 nm).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print(f"Computing metadynamics convergence from {args.hills}")
    conv_df = compute_convergence(args.hills, args.fractions, args.bound_cutoff_nm)

    if conv_df.empty:
        raise RuntimeError("No FES could be computed. Check that PLUMED is available and HILLS file is valid.")

    print(conv_df.to_string(index=False))
    plot_convergence(conv_df, args.convergence)
    print(f"Convergence plot saved to {args.convergence}")


if __name__ == "__main__":
    main()
