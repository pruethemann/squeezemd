"""Unit tests for FoldX mutation-string parsing/validation (no MDAnalysis needed)."""

import pytest

from squeezemd.preprocessing.mutate_structure import build_mutant_file_content

WT = "ARNDCEQGHILKMFPSTWYV"  # 20-residue toy sequence, one of each


def test_single_point_mutation():
    # Position 1 is 'A' -> mutate to 'G'
    content = build_mutant_file_content(WT, "A1G")
    wt_line, mut_line = content.split("\n")
    assert wt_line == WT
    assert mut_line == "G" + WT[1:]


def test_multiple_mutations_applied_in_order():
    # A1G and N3E (position 3 is 'N')
    content = build_mutant_file_content(WT, "A1G_N3E")
    _wt_line, mut_line = content.split("\n")
    assert mut_line[0] == "G"
    assert mut_line[2] == "E"
    # untouched positions preserved
    assert mut_line[1] == WT[1]


def test_wrong_wildtype_residue_raises():
    # Position 1 is 'A', so claiming 'C1G' must fail
    with pytest.raises(ValueError, match="Wrong wild-type residue"):
        build_mutant_file_content(WT, "C1G")


def test_out_of_range_position_raises():
    with pytest.raises(ValueError, match="out of range"):
        build_mutant_file_content(WT, "A99G")


def test_unparseable_mutation_raises():
    with pytest.raises(ValueError, match="Cannot parse"):
        build_mutant_file_content(WT, "AXG")


def test_empty_mutation_raises():
    with pytest.raises(ValueError):
        build_mutant_file_content(WT, "")
