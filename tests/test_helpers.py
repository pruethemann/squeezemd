"""Unit tests for the pure helper utilities (no MD run / native stack required)."""

import pytest

from squeezemd.helper_functions import (
    chain2resid,
    config_deep_update,
    import_yaml,
    is_numeric,
    parse_run_metadata,
    save_yaml,
    validate_config,
)


def _valid_config():
    return {
        "mode": "protein_protein",
        "receptors": {"C1s": {"pdb": "pdb/C1s.pdb"}},
        "ligands": ["Gigastasin"],
        "mutations": ["WT"],
        "simulation": {"replicates": 1, "time_ns": 0.1, "number_frames": 2},
    }


# --- config_deep_update ---------------------------------------------------


def test_config_deep_update_merges_nested_without_clobbering():
    base = {"simulation": {"time_ns": 1.0, "system": {"temperature_K": 300}}, "mode": "x"}
    override = {"simulation": {"system": {"temperature_K": 310, "ph": 7.0}}}
    merged = config_deep_update(base, override)
    assert merged["mode"] == "x"
    assert merged["simulation"]["time_ns"] == 1.0  # untouched
    assert merged["simulation"]["system"]["temperature_K"] == 310  # overridden
    assert merged["simulation"]["system"]["ph"] == 7.0  # added


# --- is_numeric -----------------------------------------------------------


@pytest.mark.parametrize(("char", "expected"), [("5", True), ("0", True), ("a", False), ("-", False)])
def test_is_numeric(char, expected):
    assert is_numeric(char) is expected


def test_is_numeric_rejects_multichar():
    with pytest.raises(ValueError):
        is_numeric("12")


# --- validate_config ------------------------------------------------------


def test_validate_config_accepts_complete_config():
    validate_config(_valid_config())  # should not raise


@pytest.mark.parametrize(
    "mutate",
    [
        lambda c: c.pop("mode"),
        lambda c: c.update(mode=None),
        lambda c: c.update(ligands=[]),
        lambda c: c.update(receptors={}),
        lambda c: c.update(simulation={"replicates": 1}),  # missing time_ns/number_frames
    ],
)
def test_validate_config_rejects_incomplete(mutate):
    cfg = _valid_config()
    mutate(cfg)
    with pytest.raises(ValueError):
        validate_config(cfg)


# --- import_yaml / save_yaml ---------------------------------------------


def test_yaml_round_trip(tmp_path):
    data = {"a": 1, "nested": {"b": [1, 2, 3]}}
    path = tmp_path / "cfg.yaml"
    save_yaml(data, path)
    assert import_yaml(path) == data


def test_import_yaml_raises_on_malformed(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("a: [1, 2\nb: : :")
    with pytest.raises(ValueError):
        import_yaml(path)


# --- parse_run_metadata ---------------------------------------------------


def test_parse_run_metadata_structure_path():
    meta = parse_run_metadata("protein_protein/C1s_Gigastasin/WT/513/MD/structure_end.cif")
    assert meta == {"mode": "protein_protein", "complex": "C1s_Gigastasin", "mutation": "WT", "seed": "513"}


def test_parse_run_metadata_centered_path_with_underscored_ligand():
    meta = parse_run_metadata("protein_molecule/MASP2_mol_21/R65E/777/MD/center/trajectory_centered.dcd")
    assert meta["complex"] == "MASP2_mol_21"
    assert meta["mutation"] == "R65E"
    assert meta["seed"] == "777"


def test_parse_run_metadata_rejects_path_without_md():
    with pytest.raises(ValueError):
        parse_run_metadata("results/posco/posco_interactions.parquet")


# --- chain2resid ----------------------------------------------------------


def test_chain2resid_maps_amber_numbering(tmp_path):
    mapping = tmp_path / "renum.txt"
    mapping.write_text("ALA A 10 ALA 1\nGLY A 11 GLY 2\nSER B 20 SER 3\nLEU B 21 LEU 4\n")
    chains = chain2resid(mapping)
    assert chains.loc["A", "start"] == 10
    assert chains.loc["A", "end"] == 11
    assert chains.loc["A", "amber_start"] == 1
    assert chains.loc["A", "amber_end"] == 2
    assert chains.loc["B", "start"] == 20
    assert chains.loc["B", "amber_end"] == 4
