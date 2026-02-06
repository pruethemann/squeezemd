"""Legacy Streamlit config builder (v1).

This is an older configuration UI retained for reference.
"""

# gui/app.py
from __future__ import annotations

import io
from pathlib import Path
from typing import Literal, Optional, Dict, Any

import streamlit as st
import yaml
from pydantic import BaseModel, Field, ValidationError, model_validator


# -----------------------------
# 1) Pydantic schema (sanity checks)
# -----------------------------
Mode = Literal["single_protein", "protein_small_molecule", "protein_protein"]
Solvent = Literal["implicit", "explicit"]
Nonbonded = Literal["PME", "CutoffNonPeriodic"]


class EquilibrationConfig(BaseModel):
    NVT_heating: int = Field(default=5000, ge=0)
    NPT_equilibration: int = Field(default=50000, ge=0)
    NPT_unrestrained: int = Field(default=100000, ge=0)


class MetadynamicsConfig(BaseModel):
    enabled: bool = False
    cv: Literal["COM_distance", "contacts"] = "COM_distance"
    pace: int = Field(default=500, ge=1)
    hill_height: float = Field(default=1.2, gt=0)
    sigma: float = Field(default=0.1, gt=0)


class SimulationConfig(BaseModel):
    replicates: int = Field(default=1, ge=1, le=50)
    time_ns: float = Field(default=1.0, gt=0, le=5000)
    recordingInterval_ps: float = Field(default=10.0, gt=0)
    seed: int = Field(default=2025, ge=0)

    solvent: Solvent = "explicit"
    nonbondedMethod: Nonbonded = "PME"
    cutoff_nm: float = Field(default=1.0, gt=0)

    equilibration: EquilibrationConfig = EquilibrationConfig()
    metadynamics: MetadynamicsConfig = MetadynamicsConfig()

    @model_validator(mode="after")
    def sanity_checks(self):
        # Rough frame count check
        n_frames = (self.time_ns * 1e6) / self.recordingInterval_ps  # ns -> ps
        if n_frames > 2_000_000:
            raise ValueError(
                f"Too many frames ({int(n_frames):,}). "
                "Increase recordingInterval_ps or reduce time_ns."
            )

        if self.solvent == "implicit" and self.nonbondedMethod == "PME":
            raise ValueError("PME is not compatible with implicit solvent.")

        if self.cutoff_nm > 2.0:
            st.warning("cutoff_nm > 2.0 nm is unusual; make sure you want this.")

        return self


class SystemConfig(BaseModel):
    mode: Mode = "protein_small_molecule"

    # These paths are examples; adapt to your pipeline conventions
    protein_path: str = "inputs/protein.pdb"
    ligand_path: Optional[str] = "inputs/ligand.sdf"
    partner_path: Optional[str] = None  # for protein_protein

    chain_A: str = "A"
    chain_B: Optional[str] = "B"  # used in protein_protein

    @model_validator(mode="after")
    def mode_checks(self):
        if self.mode == "protein_small_molecule":
            if not self.ligand_path:
                raise ValueError("ligand_path is required for protein_small_molecule mode.")
        if self.mode == "protein_protein":
            if not self.partner_path:
                raise ValueError("partner_path is required for protein_protein mode.")
            if not self.chain_B:
                raise ValueError("chain_B is required for protein_protein mode.")
        if self.mode == "single_protein":
            # no ligand, no partner
            return self
        return self


class SqueezeMDConfig(BaseModel):
    system: SystemConfig = SystemConfig()
    simulation: SimulationConfig = SimulationConfig()

    @model_validator(mode="after")
    def cross_checks(self):
        # MetaD allowed only for certain modes (customize as you like)
        if self.simulation.metadynamics.enabled and self.system.mode == "single_protein":
            raise ValueError("Metadynamics enabled, but mode is single_protein. Disable MetaD or switch mode.")

        return self


# -----------------------------
# 2) YAML helpers
# -----------------------------
def load_yaml(path: Path) -> Dict[str, Any]:
    data = yaml.safe_load(path.read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML must contain a mapping (top-level dict).")
    return data


def dump_yaml(data: Dict[str, Any]) -> str:
    return yaml.safe_dump(data, sort_keys=False, default_flow_style=False)


def validate_config(raw: Dict[str, Any]) -> SqueezeMDConfig:
    # Pydantic v2: use model_validate
    return SqueezeMDConfig.model_validate(raw)


# -----------------------------
# 3) Streamlit UI
# -----------------------------
st.set_page_config(page_title="squeezeMD Config Builder", layout="wide")
st.title("squeezeMD Config Builder")

# Sidebar: load / save
st.sidebar.header("Config I/O")

uploaded = st.sidebar.file_uploader("Load YAML", type=["yml", "yaml"])
default_path_str = st.sidebar.text_input("Or load from path", value="config/md_config.yml")
load_from_path = st.sidebar.button("Load from path")

if "raw_config" not in st.session_state:
    # default config dict
    st.session_state.raw_config = {
        "system": SystemConfig().model_dump(),
        "simulation": SimulationConfig().model_dump(),
    }

# Load config (upload has priority)
if uploaded is not None:
    try:
        st.session_state.raw_config = yaml.safe_load(uploaded.getvalue()) or {}
        st.sidebar.success("Loaded uploaded YAML.")
    except Exception as e:
        st.sidebar.error(f"Failed to load uploaded YAML: {e}")

elif load_from_path:
    try:
        p = Path(default_path_str)
        st.session_state.raw_config = load_yaml(p)
        st.sidebar.success(f"Loaded from {p}.")
    except Exception as e:
        st.sidebar.error(f"Failed to load from path: {e}")


raw = st.session_state.raw_config

# Validate current raw config (for pre-fill errors)
cfg: Optional[SqueezeMDConfig] = None
validation_error: Optional[str] = None
try:
    cfg = validate_config(raw)
except ValidationError as e:
    validation_error = e.__str__()
except Exception as e:
    validation_error = str(e)

# Layout columns
col_form, col_yaml = st.columns([1.1, 0.9], gap="large")

with col_form:
    st.subheader("Parameters")

    # ---- SYSTEM ----
    st.markdown("### System")
    system = raw.get("system", {}) if isinstance(raw.get("system", {}), dict) else {}

    mode = st.selectbox(
        "Mode",
        options=["single_protein", "protein_small_molecule", "protein_protein"],
        index=["single_protein", "protein_small_molecule", "protein_protein"].index(system.get("mode", "protein_small_molecule")),
    )
    protein_path = st.text_input("Protein path", value=system.get("protein_path", "inputs/protein.pdb"))
    chain_A = st.text_input("Chain A", value=system.get("chain_A", "A"))

    ligand_path = system.get("ligand_path", "inputs/ligand.sdf")
    partner_path = system.get("partner_path", None)
    chain_B = system.get("chain_B", "B")

    if mode == "protein_small_molecule":
        ligand_path = st.text_input("Ligand path", value=ligand_path or "inputs/ligand.sdf")
        partner_path = None
        chain_B = None

    elif mode == "protein_protein":
        partner_path = st.text_input("Partner path", value=partner_path or "inputs/partner.pdb")
        chain_B = st.text_input("Chain B", value=chain_B or "B")
        ligand_path = None

    else:  # single_protein
        ligand_path = None
        partner_path = None
        chain_B = None

    # write back
    raw["system"] = {
        "mode": mode,
        "protein_path": protein_path,
        "ligand_path": ligand_path,
        "partner_path": partner_path,
        "chain_A": chain_A,
        "chain_B": chain_B,
    }

    # ---- SIMULATION ----
    st.markdown("### Simulation")
    sim = raw.get("simulation", {}) if isinstance(raw.get("simulation", {}), dict) else {}

    replicates = st.number_input("Replicates", min_value=1, max_value=50, value=int(sim.get("replicates", 1)))
    time_ns = st.number_input("Time (ns)", min_value=0.0001, max_value=5000.0, value=float(sim.get("time_ns", 1.0)))
    recording_ps = st.number_input("Recording interval (ps)", min_value=0.001, value=float(sim.get("recordingInterval_ps", 10.0)))
    seed = st.number_input("Seed", min_value=0, value=int(sim.get("seed", 2025)))

    solvent = st.selectbox("Solvent", options=["explicit", "implicit"], index=["explicit", "implicit"].index(sim.get("solvent", "explicit")))
    nonbonded = st.selectbox("Nonbonded method", options=["PME", "CutoffNonPeriodic"], index=["PME", "CutoffNonPeriodic"].index(sim.get("nonbondedMethod", "PME")))
    cutoff_nm = st.number_input("Nonbonded cutoff (nm)", min_value=0.1, max_value=5.0, value=float(sim.get("cutoff_nm", 1.0)))

    st.markdown("#### Equilibration")
    eq = sim.get("equilibration", {}) if isinstance(sim.get("equilibration", {}), dict) else {}
    NVT_heating = st.number_input("NVT heating steps", min_value=0, value=int(eq.get("NVT_heating", 5000)))
    NPT_equil = st.number_input("NPT equilibration steps", min_value=0, value=int(eq.get("NPT_equilibration", 50000)))
    NPT_unrest = st.number_input("NPT unrestrained steps", min_value=0, value=int(eq.get("NPT_unrestrained", 100000)))

    st.markdown("#### Metadynamics")
    meta = sim.get("metadynamics", {}) if isinstance(sim.get("metadynamics", {}), dict) else {}
    meta_enabled = st.checkbox("Enable metadynamics", value=bool(meta.get("enabled", False)))

    meta_cv = meta.get("cv", "COM_distance")
    meta_pace = int(meta.get("pace", 500))
    meta_hill = float(meta.get("hill_height", 1.2))
    meta_sigma = float(meta.get("sigma", 0.1))

    if meta_enabled:
        meta_cv = st.selectbox("CV", options=["COM_distance", "contacts"], index=["COM_distance", "contacts"].index(meta_cv))
        meta_pace = st.number_input("Pace (steps)", min_value=1, value=meta_pace)
        meta_hill = st.number_input("Hill height", min_value=1e-6, value=meta_hill)
        meta_sigma = st.number_input("Sigma", min_value=1e-6, value=meta_sigma)

    raw["simulation"] = {
        "replicates": int(replicates),
        "time_ns": float(time_ns),
        "recordingInterval_ps": float(recording_ps),
        "seed": int(seed),
        "solvent": solvent,
        "nonbondedMethod": nonbonded,
        "cutoff_nm": float(cutoff_nm),
        "equilibration": {
            "NVT_heating": int(NVT_heating),
            "NPT_equilibration": int(NPT_equil),
            "NPT_unrestrained": int(NPT_unrest),
        },
        "metadynamics": {
            "enabled": bool(meta_enabled),
            "cv": meta_cv,
            "pace": int(meta_pace),
            "hill_height": float(meta_hill),
            "sigma": float(meta_sigma),
        },
    }

    st.divider()

    # Validate button
    if st.button("Validate config", type="primary"):
        try:
            cfg = validate_config(raw)
            st.success("Config is valid ✔")
        except Exception as e:
            st.error(f"Config invalid: {e}")

    # Quick “effective summary”
    if cfg:
        st.info(
            f"**Summary:** mode={cfg.system.mode}, time_ns={cfg.simulation.time_ns}, "
            f"recordingInterval_ps={cfg.simulation.recordingInterval_ps}, "
            f"replicates={cfg.simulation.replicates}, MetaD={cfg.simulation.metadynamics.enabled}"
        )
    elif validation_error:
        st.warning("Current config has validation errors (see YAML panel).")


with col_yaml:
    st.subheader("YAML Preview")

    yaml_text = dump_yaml(raw)
    st.code(yaml_text, language="yaml")

    if validation_error:
        st.error("Validation error")
        st.code(validation_error)

    st.download_button(
        "Download YAML",
        data=yaml_text.encode("utf-8"),
        file_name="md_config.yml",
        mime="text/yaml",
    )

    save_path = st.text_input("Save to path", value="config/md_config.yml", key="save_path")
    if st.button("Save to path"):
        try:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            Path(save_path).write_text(yaml_text)
            st.success(f"Saved to {save_path}")
        except Exception as e:
            st.error(f"Failed to save: {e}")