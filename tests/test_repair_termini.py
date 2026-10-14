"""Unit tests for the Maestro C-terminus repair (no MD stack required)."""

import math

import pytest

from squeezemd.preprocessing.repair_termini import build_oxt_position, repair_cterm

# A minimal C-terminal residue: N, CA, C, O and the offending HXT, plus an ANISOU
# record that the repair is expected to drop.
_CTERM = [
    "ATOM      1  N   ILE B 124      74.385  45.573  79.586  1.00 87.80           N  ",
    "ATOM      2  CA  ILE B 124      75.481  45.467  80.575  1.00 87.25           C  ",
    "ANISOU    2  CA  ILE B 124     8241   7589   5838   1197   1817   3243       C  ",
    "ATOM      3  C   ILE B 124      74.835  45.168  81.926  1.00 92.69           C  ",
    "ATOM      4  O   ILE B 124      74.499  44.025  82.250  1.00 88.49           O  ",
    "ATOM      5  HXT ILE B 124      74.500  46.000  82.500  1.00 92.69           H  ",
]


def _atoms(lines):
    return {line[12:16].strip(): line for line in lines if line.startswith(("ATOM  ", "HETATM"))}


def _xyz(line):
    return tuple(float(line[s]) for s in (slice(30, 38), slice(38, 46), slice(46, 54)))


def test_hxt_becomes_oxt_with_carboxylate_geometry():
    out, repaired = repair_cterm(_CTERM, expected_repairs=1)
    atoms = _atoms(out)

    assert repaired == [("B", " 124", " ")]
    assert "HXT" not in atoms
    assert atoms["OXT"][76:78] == " O"

    c, o, ca, oxt = (_xyz(atoms[n]) for n in ("C", "O", "CA", "OXT"))
    bond = math.dist(c, oxt)
    assert bond == pytest.approx(1.25, abs=1e-3)

    def angle(p, q):
        u = [a - b for a, b in zip(p, c, strict=True)]
        v = [a - b for a, b in zip(q, c, strict=True)]
        dot = sum(a * b for a, b in zip(u, v, strict=True))
        nu = math.sqrt(sum(a * a for a in u))
        nv = math.sqrt(sum(a * a for a in v))
        return math.degrees(math.acos(dot / (nu * nv)))

    # sp2 carboxylate: both angles near 120 degrees
    assert angle(ca, oxt) == pytest.approx(120, abs=3)
    assert angle(o, oxt) == pytest.approx(120, abs=3)


def test_serial_is_reused_so_conect_stays_valid():
    out, _ = repair_cterm(_CTERM, expected_repairs=1)
    assert _atoms(out)["OXT"][6:11] == "    5"


def test_anisou_records_are_dropped():
    out, _ = repair_cterm(_CTERM, expected_repairs=1)
    assert not any(line.startswith("ANISOU") for line in out)
    assert any(line.startswith("ANISOU") for line in repair_cterm(_CTERM, drop_anisou=False)[0])


def test_residue_with_existing_oxt_is_left_alone():
    already_fine = [line for line in _CTERM if "HXT" not in line]
    out, repaired = repair_cterm(already_fine, expected_repairs=0)
    assert repaired == []
    assert "OXT" not in _atoms(out)


def test_wrong_repair_count_raises():
    with pytest.raises(ValueError, match="Expected 2 C-termini"):
        repair_cterm(_CTERM, expected_repairs=2)


def test_missing_backbone_atom_raises():
    without_o = [line for line in _CTERM if line[12:16].strip() != "O"]
    with pytest.raises(ValueError, match="missing backbone atom"):
        repair_cterm(without_o)


def test_degenerate_geometry_raises():
    # CA, C, O collinear -> the bisector vanishes and OXT cannot be placed.
    with pytest.raises(ValueError, match="degenerate"):
        build_oxt_position((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))
