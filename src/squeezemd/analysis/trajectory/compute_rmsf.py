#!/usr/bin/env python

"""Compute RMSF across one or more trajectories and save as parquet."""

import argparse

import MDAnalysis as mda
import openmm.app as app
import pandas as pd
from MDAnalysis.analysis import rms

from ...helper_functions import parse_run_metadata, remap_MDAnalysis


def calculate_RMSF(u: mda.Universe, i, metadata=None):
    """Calculate Cα RMSF for a single trajectory and tag it with run metadata.

    ``metadata`` (complex/mutation/seed, from :func:`parse_run_metadata`) is added
    as columns so the aggregated table can be traced back to and grouped by the
    exact simulation that produced each curve.
    """

    # TODO: separate ligand and receptor; currently all Cα atoms are pooled.
    c_alphas = u.select_atoms("name CA")
    R = rms.RMSF(c_alphas).run()

    rmsf_df = pd.DataFrame({"resid": c_alphas.resids, "rmsf": R.results.rmsf, "sim_id": i})
    for key, value in (metadata or {}).items():
        rmsf_df[key] = value

    return rmsf_df


def parse_arguments():
    """Parse CLI arguments for RMSF computation."""
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("--topo", nargs="+", required=False)
    parser.add_argument("--traj", nargs="+", required=False)

    # Output
    parser.add_argument("--output", required=False, default="rmsf.svg", help="")

    return parser.parse_args()


def main():
    args = parse_arguments()

    topos = sorted(args.topo)
    trajs = sorted(args.traj)

    rmsf_data = []

    for i, (topo_path, traj) in enumerate(zip(topos, trajs, strict=True)):
        # Recover the run identity from the topology path before opening the file.
        metadata = parse_run_metadata(topo_path)

        # Import Trajectory
        topo = app.PDBxFile(topo_path)
        u = mda.Universe(topo, traj, in_memory=False)
        u = remap_MDAnalysis(u, topo)

        # Calculate RMSF
        rmsf = calculate_RMSF(u, i, metadata)
        rmsf_data.append(rmsf)

    rmsf = pd.concat(rmsf_data)

    rmsf.to_parquet(args.output)


if __name__ == "__main__":
    main()
