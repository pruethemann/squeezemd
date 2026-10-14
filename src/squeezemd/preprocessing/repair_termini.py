#!/usr/bin/env python

"""Repair Maestro-style neutral C-termini so Amber force fields can template them.

Schrodinger/Maestro writes a truncated protein C-terminus as an ``HXT`` hydrogen
bonded directly to the backbone carbonyl carbon, with no ``OXT``. That is an
aldehyde-like terminus which matches no Amber template, so
``Modeller.addHydrogens`` fails with::

    No template found for residue N (ILE). The set of heavy atoms matches NILE,
    but the residue is missing 1 H atom.

Deleting ``HXT`` is not sufficient - the residue then looks like a mid-chain
residue missing its successor. The fix is to turn it into a charged carboxylate
by replacing ``HXT`` with an ``OXT`` oxygen at proper sp2 geometry.
"""

import numpy as np

# PDB fixed-column slices (0-based, end-exclusive).
_NAME = slice(12, 16)
_RESNAME = slice(17, 20)
_CHAIN = slice(21, 22)
_RESSEQ = slice(22, 26)
_ICODE = slice(26, 27)
_XYZ = (slice(30, 38), slice(38, 46), slice(46, 54))
_ELEMENT = slice(76, 78)

_C_O_BOND_A = 1.25  # carboxylate C-O bond length in Angstrom


def _is_atom(line):
    return line.startswith(("ATOM  ", "HETATM"))


def _residue_key(line):
    """Chain + residue number + insertion code, i.e. what identifies one residue."""
    return (line[_CHAIN], line[_RESSEQ], line[_ICODE])


def _coords(line):
    return np.array([float(line[s]) for s in _XYZ])


def _set_coords(line, xyz):
    line = line.ljust(80)
    for s, value in zip(_XYZ, xyz, strict=True):
        line = line[: s.start] + f"{value:8.3f}" + line[s.stop :]
    return line


def build_oxt_position(ca, c, o):
    """Place OXT as the third trigonal-planar substituent on the backbone carbon.

    The two existing substituents are CA and O. The missing one points opposite
    their bisector, which puts OXT in the CA-C-O plane at ~120 degrees to both.
    """
    ca, c, o = (np.asarray(v, dtype=float) for v in (ca, c, o))
    u_ca = (ca - c) / np.linalg.norm(ca - c)
    u_o = (o - c) / np.linalg.norm(o - c)
    bisector = u_ca + u_o
    norm = np.linalg.norm(bisector)
    if norm < 1e-6:  # degenerate: CA, C, O collinear. Should not happen in real structures.
        raise ValueError("Cannot place OXT: CA-C-O geometry is degenerate")
    return c - _C_O_BOND_A * (bisector / norm)


def repair_cterm(lines, expected_repairs=None, drop_anisou=True):
    """Convert every ``HXT``-without-``OXT`` C-terminus into a charged carboxylate.

    :param lines: PDB file contents as a list of lines (newline-terminated or not).
    :param expected_repairs: if given, raise unless exactly this many residues were repaired.
    :param drop_anisou: discard ANISOU records (they roughly double file size and
                        nothing in the pipeline reads them).
    :return: ``(repaired_lines, repaired_residue_keys)``
    """
    lines = [line.rstrip("\n") for line in lines]

    # Index the atoms of every residue so we can find HXT and its backbone partners.
    residues = {}
    for i, line in enumerate(lines):
        if _is_atom(line):
            residues.setdefault(_residue_key(line), {})[line[_NAME].strip()] = i

    targets = [key for key, atoms in residues.items() if "HXT" in atoms and "OXT" not in atoms]

    repaired = []
    for key in targets:
        atoms = residues[key]
        missing = [name for name in ("CA", "C", "O") if name not in atoms]
        if missing:
            raise ValueError(f"Residue {key} has HXT but is missing backbone atom(s) {missing}")

        oxt = build_oxt_position(_coords(lines[atoms["CA"]]), _coords(lines[atoms["C"]]), _coords(lines[atoms["O"]]))

        # Reuse the HXT record: same serial (so CONECT stays valid), new name/element/position.
        idx = atoms["HXT"]
        line = _set_coords(lines[idx], oxt)
        line = line[: _NAME.start] + " OXT" + line[_NAME.stop :]
        line = line[: _ELEMENT.start] + " O" + line[_ELEMENT.stop :]
        lines[idx] = line.rstrip()
        repaired.append(key)

    if expected_repairs is not None and len(repaired) != expected_repairs:
        raise ValueError(f"Expected {expected_repairs} C-termini to repair, found {len(repaired)}: {repaired}")

    if drop_anisou:
        lines = [line for line in lines if not line.startswith("ANISOU")]

    return [line + "\n" for line in lines], repaired


def repair_cterm_file(in_pdb, out_pdb, expected_repairs=None, drop_anisou=True):
    """Read a PDB, repair its C-termini, and write the result. Returns repaired keys."""
    with open(in_pdb) as handle:
        lines = handle.readlines()

    repaired_lines, repaired = repair_cterm(lines, expected_repairs=expected_repairs, drop_anisou=drop_anisou)

    with open(out_pdb, "w") as handle:
        handle.writelines(repaired_lines)

    return repaired
