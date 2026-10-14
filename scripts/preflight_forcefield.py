#!/usr/bin/env python

"""Check that structures can be templated by the force field, before burning GPU time.

Runs exactly the preparation path squeezeMD's ``run-md`` uses
(``PDBFile`` -> ``Modeller`` -> ``addHydrogens`` -> ``addExtraParticles``) and
reports per file. Takes seconds, needs no GPU, and catches missing residue
templates and malformed termini that would otherwise fail every MD job.

Must run in the ``squeeze`` conda env (needs openmm).

Usage:
    python scripts/preflight_forcefield.py demo/S-1-C4B/pdb_md/*.pdb
    python scripts/preflight_forcefield.py --md-config demo/S-1-C4B/config/md_config.yaml <pdbs>
"""

import argparse
import sys
from pathlib import Path

import openmm.app as app
import yaml


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdbs", nargs="+", type=Path)
    parser.add_argument(
        "--md-config",
        type=Path,
        help="Read forcefield/pH from this md_config.yaml instead of the amber19 defaults.",
    )
    parser.add_argument("--ph", type=float, default=7.4)
    return parser.parse_args()


def main():
    args = parse_arguments()

    protein_ff, water_ff, ph = "amber19-all.xml", "amber19/tip4pew.xml", args.ph
    if args.md_config:
        sim = yaml.safe_load(args.md_config.read_text())["simulation"]
        protein_ff = sim["forcefield"]["protein"]
        water_ff = sim["forcefield"]["water"]
        ph = sim.get("system", {}).get("ph", args.ph)

    print(f"forcefield: {protein_ff} + {water_ff}   pH {ph}\n")
    forcefield = app.ForceField(protein_ff, water_ff)

    failures = 0
    for pdb_path in args.pdbs:
        try:
            pdb = app.PDBFile(str(pdb_path))
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(forcefield, pH=ph)
            modeller.addExtraParticles(forcefield)
            print(f"  PASS  {pdb_path.name:12s} {modeller.topology.getNumAtoms()} atoms")
        except Exception as exc:  # noqa: BLE001 - we want to report any failure, not crash
            failures += 1
            print(f"  FAIL  {pdb_path.name:12s} {str(exc).splitlines()[0]}")

    print(f"\n{len(args.pdbs) - failures} pass / {failures} fail")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
