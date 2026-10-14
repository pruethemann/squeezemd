#!/usr/bin/env python

"""Generate MD-ready copies of the S-1-C4B complexes.

The PrepWizard exports in ``demo/S-1-C4B/pdb/`` cannot be parameterised by Amber:
chains B and D end in a Maestro ``HXT`` pseudo-aldehyde instead of a carboxylate.
This writes repaired copies to ``demo/S-1-C4B/pdb_md/``, which is what the
squeezeMD configs reference. ``pdb/`` is left untouched as the provenance record.

Usage:
    python scripts/prepare_c4b_for_md.py
    python scripts/prepare_c4b_for_md.py --in-dir ... --out-dir ... --expected-repairs 2
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from squeezemd.preprocessing.repair_termini import repair_cterm_file  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-dir", type=Path, default=REPO / "demo" / "S-1-C4B" / "pdb")
    parser.add_argument("--out-dir", type=Path, default=REPO / "demo" / "S-1-C4B" / "pdb_md")
    parser.add_argument(
        "--expected-repairs",
        type=int,
        default=2,
        help="C-termini expected per file (2 = chains B and D). Use -1 to disable the check.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    expected = None if args.expected_repairs < 0 else args.expected_repairs

    pdbs = sorted(args.in_dir.glob("*.pdb"))
    if not pdbs:
        raise SystemExit(f"No PDB files found in {args.in_dir}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for pdb in pdbs:
        out = args.out_dir / pdb.name
        repaired = repair_cterm_file(pdb, out, expected_repairs=expected)
        chains = ",".join(sorted(key[0] for key in repaired))
        print(f"{pdb.name:12s} -> {out.relative_to(REPO)}  (repaired chains: {chains})")

    print(f"\n{len(pdbs)} structures written to {args.out_dir.relative_to(REPO)}")


if __name__ == "__main__":
    main()
