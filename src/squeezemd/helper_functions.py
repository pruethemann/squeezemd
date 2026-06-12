#!/usr/bin/env python

"""Shared helper utilities used across the package.

This module groups small utilities used by multiple scripts, including
YAML handling, command execution, and residue/chain remapping helpers
for MDAnalysis/OpenMM interoperability.

The heavy MD dependencies (MDAnalysis, OpenMM) are imported lazily inside the
functions that need them, so the pure config/path helpers can be imported and
unit-tested without the full native stack installed.
"""

from __future__ import annotations

import os
import subprocess
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import yaml

if TYPE_CHECKING:
    import MDAnalysis as mda
    import openmm.app as app


def parse_run_metadata(path) -> dict:
    """Recover the run identity from a squeezemd result path.

    Result paths follow ``<mode>/<receptor>_<ligand>/<mutation>/<seed>/MD/...``.
    Using the ``MD`` directory as an anchor, return the ``mode``, ``complex``
    (``receptor_ligand``), ``mutation`` and ``seed`` so per-run analysis outputs
    can be traced back to the exact simulation that produced them.
    """
    parts = Path(path).parts
    if "MD" not in parts:
        raise ValueError(f"Cannot parse run metadata from path without an 'MD' component: {path!r}")
    md_idx = len(parts) - 1 - parts[::-1].index("MD")  # last 'MD'
    if md_idx < 4:
        raise ValueError(f"Path is too shallow to contain <mode>/<complex>/<mutation>/<seed>/MD/...: {path!r}")
    mode, complex_name, mutation, seed = parts[md_idx - 4 : md_idx]
    return {"mode": mode, "complex": complex_name, "mutation": mutation, "seed": seed}


def update_md_overview(config):
    exp = {}

    # Find a solution for that. Save it in conda env?
    md_overview = pd.read_parquet("/home/peter/caracara/Squeeze/md_overview.parquet")

    exp["ID"] = config["ID"]
    exp["Time"] = config["simulation"]["time_ns"]
    exp["Replicates"] = config["simulation"]["replicates"]

    exp["Receptor"] = ",".join(config["receptors"].keys())
    exp["Ligands"] = ",".join(config["ligands"])

    exp["Mode"] = config["mode"]
    exp["Name"] = config["name"]

    exp_id = exp.pop("ID")
    exp_df = pd.DataFrame([exp], index=pd.Index([exp_id], name="ID"))

    if exp_id in md_overview.index:
        md_overview.loc[exp_id] = exp_df.loc[exp_id]
    else:
        md_overview = pd.concat([md_overview, exp_df])

    md_overview.to_parquet("/home/peter/caracara/Squeeze/md_overview.parquet")

    print(md_overview)


def setup_testrun(config):
    if "debug" in config and config["debug"]:
        print("ATTENTION: This is a testrun")
        test_md_config = files("squeezemd").joinpath("resources", "md_test_config.yaml")
        md_test = import_yaml(test_md_config)
        config = config_deep_update(config, md_test)

        # Keep only the first few ligands to speed up the test run
        config["ligands"] = (config.get("ligands") or [])[0:3]

        save_yaml(config, "config/md_test_config.yaml")

        return ("config/md_test_config.yaml", config)

    print("Production run")
    return ("config/md_config.yaml", config)


def config_deep_update(base: dict, override: dict) -> dict:
    """Recursively merge nested dicts, updating ``base`` in place.

    Parameters
    ----------
    base
        Original configuration dictionary.
    override
        New values to merge into ``base``.
    """
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            config_deep_update(base[key], value)
        else:
            base[key] = value
    return base


def save_file(content, output_file):
    """
    Saves a string (content) to a text file and closes it.
    :param content:
    :param output_file:
    :return:
    """

    with open(output_file, "w") as file:
        file.write(content)


def execute(command):
    """
    Executes commands in console
    :param command:
    :return:
    """

    output_text = subprocess.check_output(command, shell=True)
    return output_text


def import_yaml(yaml_path: os.path):
    """
    Open a YAML configuration file and return it as a dict.

    Raises a clear error if the file is missing or malformed instead of
    silently returning ``None`` (which used to surface much later as a
    confusing ``NoneType`` error).

    :param yaml_path: File path to yaml
    :return: dictionary with parameters
    """
    try:
        with open(yaml_path) as stream:
            return yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise ValueError(f"Could not parse YAML file '{yaml_path}': {exc}") from exc


# Keys that must be present in the merged sim+md configuration.
REQUIRED_CONFIG_KEYS = ("mode", "receptors", "ligands", "mutations", "simulation")
REQUIRED_SIMULATION_KEYS = ("replicates", "time_ns", "number_frames")


def validate_config(config: dict) -> None:
    """Fail fast with a clear message if the merged config is incomplete.

    Called from the Snakefile after sim_config.yaml and md_config.yaml have been
    merged. Catches the common mistakes (missing ``mode``, empty ``ligands``,
    missing ``simulation`` block) before any expensive rule runs.
    """
    missing = [key for key in REQUIRED_CONFIG_KEYS if key not in config or config[key] is None]
    if missing:
        raise ValueError(
            f"Invalid squeezemd configuration: missing required key(s) {missing}. "
            "Define them in config/sim_config.yaml (mode, receptors, ligands) "
            "and config/md_config.yaml (simulation)."
        )

    if not config["ligands"]:
        raise ValueError(
            "Invalid squeezemd configuration: 'ligands' is empty. "
            "List at least one ligand (for protein-protein runs, the name of the "
            "binding partner present in the PDB, e.g. 'Gigastasin')."
        )

    if not config["receptors"]:
        raise ValueError("Invalid squeezemd configuration: 'receptors' is empty.")

    sim = config["simulation"]
    missing_sim = [key for key in REQUIRED_SIMULATION_KEYS if key not in sim]
    if missing_sim:
        raise ValueError(
            f"Invalid squeezemd configuration: config['simulation'] is missing {missing_sim}. "
            "Check config/md_config.yaml."
        )


def extract_ligand_sequence(pdb_ligand: os.path):
    """Extract a single-chain ligand sequence from a PDB file.

    Assumes the ligand is on chain A and normalizes CYX -> CYS for
    compatibility with sequence extraction.
    """
    import MDAnalysis as mda

    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # --- Normalize CYX → CYS ---
    for res in u.residues:
        if res.resname == "CYX":
            res.resname = "CYS"

    # Extract ligand at chain A
    ligand = u.select_atoms("chainID A")

    # Return sequence
    sequence = ligand.residues.sequence().seq

    return str(sequence)


def save_yaml(d, filepath):
    """
    Save a yml file
    :param d:
    :param filepath:
    :return:
    """
    with open(filepath, "w") as file:
        yaml.dump(d, file)


def chain2resid(file_csv):
    """
    Amber preparation removes chain ids and starts renumbering residues from 1.
    This function remaps the numbering based on the amber mapping file.
    :param file_csv:
    :return:
    """
    # Find start and end of chain A

    renum = pd.read_csv(file_csv, sep=r"\s+", names=["resname", "chainID", "resid", "resname amber", "resid amber"])

    del renum["resname"]
    del renum["resname amber"]

    chain_min = renum.groupby("chainID").min().rename(columns={"resid amber": "amber_start", "resid": "start"})
    chain_max = renum.groupby("chainID").max().rename(columns={"resid amber": "amber_end", "resid": "end"})

    chains = pd.concat([chain_min, chain_max], axis=1)
    return chains


def is_numeric(character):
    """
    This function checks if a given character is numeric.

    :param character: A single character (string) to check.
    :return: True if the character is numeric, False otherwise.
    """
    if len(character) != 1:
        raise ValueError("Input must be a single character.")
    return character.isdigit()


def remap_MDAnalysis(u: mda.Universe, topo: app.PDBxFile):
    """
    Remaps the correct residue and chain IDs from the OpenMM PDBxFile topology
    to an MDAnalysis universe.

    :param u: MDAnalysis Universe
    :param topo: openmm.app.PDBxFile object
    :return: updated MDAnalysis Universe
    """
    chains = list(topo.topology.chains())
    residues = list(topo.topology.residues())

    if len(u.segments) != len(chains):
        raise ValueError("Mismatch in number of segments and chains")

    for mda_seg, omm_chain in zip(u.segments, chains, strict=True):
        mda_seg.segid = omm_chain.id  # Safe: assigns chainID

    if len(u.residues) != len(residues):
        raise ValueError("Mismatch in number of residues between MDAnalysis and OpenMM topology")

    for mda_res, omm_res in zip(u.residues, residues, strict=True):
        # Optional: Only remap if different
        mda_res.resid = int(omm_res.id)
        mda_res.resname = omm_res.name

    return u


def remap_amber(mapping_file, u):
    """
    Amber preparation removes chain ids and starts renumbering residues from 1.
    This function remaps the numbering based on the amber mapping file.
    :param args:
    :param u:
    :return:
    """

    print("2. Import reside information")
    # Get chain information since amber deleted all chain IDs
    chains = chain2resid(mapping_file)

    # IMPORTANT: chain ID A needs to be at the end, otherwise it leads to misnumbering
    chains.sort_index(ascending=False, inplace=True)

    # 1. Assign chain Ids to amber resids
    for chainID, r in chains.iterrows():
        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"resid {r.amber_start} to {r.amber_end}")
        chain.atoms.chainIDs = chainID

    # 2. Renumber resids
    for chainID, r in chains.iterrows():
        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"chainID {chainID}")

        # Shift resid numbering to old numbering
        shift_factor = r.start - int(r.amber_start)
        chain.residues.resids += shift_factor

    return u
