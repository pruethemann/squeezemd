# gui/app.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, Optional, List

import streamlit as st
import yaml
import pandas as pd
from pydantic import BaseModel, Field, ValidationError, model_validator


# -----------------------------
# YAML helpers
# -----------------------------
def load_yaml(path: Path) -> Dict[str, Any]:
    data = yaml.safe_load(path.read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML must contain a mapping (top-level dict).")
    return data


def dump_yaml(data: Dict[str, Any]) -> str:
    return yaml.safe_dump(data, sort_keys=False, default_flow_style=False)


# -----------------------------
# Conversions: steps <-> time
# -----------------------------
def steps_to_ps(steps: int, dt_fs: float) -> float:
    # dt_fs: fs/step. 1 ps = 1000 fs
    return steps * dt_fs / 1000.0


def ps_to_steps(time_ps: float, dt_fs: float) -> int:
    # round to nearest integer step
    return int(round(time_ps * 1000.0 / dt_fs))


# -----------------------------
# MD config schema (based on your md_config.yaml)
# -----------------------------
class EquilibrationConfig(BaseModel):
    protein_k: float = Field(default=10.0, ge=0)

    NVT_heating: int = Field(default=5000, ge=0)
    NPT_equilibration: int = Field(default=50000, ge=0)
    NPT_unrestrained: int = Field(default=1000, ge=0)


class ForcefieldConfig(BaseModel):
    protein: str = "amber19-all.xml"
    water: str = "amber19/tip4pew.xml"
    watermodel: str = "tip4pew"


class SystemConditionsConfig(BaseModel):
    salt_molar: float = Field(default=0.0, ge=0.0)
    temperature_K: float = Field(default=300.0, gt=0.0)
    pressure_atm: float = Field(default=1.0, gt=0.0)
    barostat_interval_steps: int = Field(default=25, ge=1)


class ConstraintsConfig(BaseModel):
    dt_fs: float = Field(default=2.0, gt=0.0)
    bonds: str = "HBonds"  # None, HBonds, AllBonds
    rigid_water: bool = True
    cutoff_nm: float = Field(default=1.0, gt=0.0)
    ewald_error_tolerance: float = Field(default=0.0005, gt=0.0)
    constraint_tolerance: float = Field(default=0.0001, gt=0.0)
    friction_ps: float = Field(default=1.0, gt=0.0)


class FlexibleBindingPocketConfig(BaseModel):
    protein_k: float = Field(default=5000.0, ge=0.0)
    flexible_resids: Dict[int, str] = Field(default_factory=dict)


class MetadynamicsConfig(BaseModel):
    CV: str = "COM"
    SIGMA: float = Field(default=0.2, gt=0.0)
    HEIGHT: float = Field(default=0.03, gt=0.0)
    PACE: int = Field(default=1000, ge=1)
    STRIDE: int = Field(default=50, ge=1)


from pydantic import model_validator

class SimulationConfig(BaseModel):
    replicates: int = Field(default=2, ge=1, le=100)
    time_ns: float = Field(default=3.0, gt=0.0, le=50000.0)

    recording_interval_ps: float = Field(default=1.0, gt=0.0)
    number_frames: int = Field(default=2000, ge=1)  # <-- give a sensible default
    seed: int = Field(default=2025, ge=0)

    # ... other fields ...

    @model_validator(mode="after")
    def sanity(self):
        total_ps = self.time_ns * 1e6

        # If number_frames is set, treat that as the controlling knob.
        # Derive the effective recording interval for info/sanity.
        if self.number_frames is not None and self.number_frames > 0:
            effective_interval_ps = total_ps / self.number_frames
            # keep your YAML's recording_interval_ps in sync (optional but nice)
            self.recording_interval_ps = effective_interval_ps

            # sanity check on number_frames (disk)
            if self.number_frames > 5_000_000:
                raise ValueError(
                    f"Too many frames ({self.number_frames:,}). Reduce number_frames."
                )

            # sanity check on effective interval (performance)
            if effective_interval_ps < 0.1:
                # don't hard-fail; this is often intentional but expensive
                # (replace with st.warning in UI layer if you prefer)
                pass

        else:
            # fall back to interval-based recording
            n_frames_est = total_ps / self.recording_interval_ps
            if n_frames_est > 5_000_000:
                raise ValueError(
                    f"Too many frames (~{int(n_frames_est):,}). "
                    "Increase recording_interval_ps or reduce time_ns."
                )

        return self


class MDConfig(BaseModel):
    simulation: SimulationConfig = SimulationConfig()


# -----------------------------
# SIM config schema (based on your sim_config.yaml)
# -----------------------------
class ComplexEntry(BaseModel):
    receptor: str
    ligand: str
    pdb: str
    sdf: str


class SimConfig(BaseModel):
    mutations: List[str] = Field(default_factory=lambda: ["WT"])
    complexes: Dict[str, ComplexEntry] = Field(default_factory=dict)

    @model_validator(mode="after")
    def sanity(self):
        if len(self.mutations) == 0:
            raise ValueError("mutations must contain at least one entry (e.g., WT).")
        if len(self.complexes) == 0:
            st.warning("No complexes defined yet.")
        return self


# -----------------------------
# App setup
# -----------------------------
st.set_page_config(page_title="squeezeMD Config Builder", layout="wide")
st.title("squeezeMD Config Builder (md_config + sim_config)")


# -----------------------------
# Sidebar: Load / save
# -----------------------------
st.sidebar.header("Config I/O")

# Defaults: point to your local repo paths if you want; keeping generic here
md_default_path = st.sidebar.text_input("MD config path", value="config/md_config.yaml")
sim_default_path = st.sidebar.text_input("SIM config path", value="config/sim_config.yaml")

md_upload = st.sidebar.file_uploader("Upload md_config.yaml", type=["yml", "yaml"], key="md_upload")
sim_upload = st.sidebar.file_uploader("Upload sim_config.yaml", type=["yml", "yaml"], key="sim_upload")

col_load1, col_load2 = st.sidebar.columns(2)
load_md_from_path = col_load1.button("Load MD from path")
load_sim_from_path = col_load2.button("Load SIM from path")


def init_defaults():
    # Minimal defaults aligned to your uploaded examples
    return (
        {"simulation": SimulationConfig().model_dump()},  # md
        {"mutations": ["WT"], "complexes": {}},           # sim
    )


if "raw_md" not in st.session_state or "raw_sim" not in st.session_state:
    st.session_state.raw_md, st.session_state.raw_sim = init_defaults()

# Load MD
if md_upload is not None:
    try:
        st.session_state.raw_md = yaml.safe_load(md_upload.getvalue()) or {}
        st.sidebar.success("Loaded uploaded md_config.yaml")
    except Exception as e:
        st.sidebar.error(f"Failed to load md_config upload: {e}")
elif load_md_from_path:
    try:
        st.session_state.raw_md = load_yaml(Path(md_default_path))
        st.sidebar.success(f"Loaded MD from {md_default_path}")
    except Exception as e:
        st.sidebar.error(f"Failed to load MD from path: {e}")

# Load SIM
if sim_upload is not None:
    try:
        st.session_state.raw_sim = yaml.safe_load(sim_upload.getvalue()) or {}
        st.sidebar.success("Loaded uploaded sim_config.yaml")
    except Exception as e:
        st.sidebar.error(f"Failed to load sim_config upload: {e}")
elif load_sim_from_path:
    try:
        st.session_state.raw_sim = load_yaml(Path(sim_default_path))
        st.sidebar.success(f"Loaded SIM from {sim_default_path}")
    except Exception as e:
        st.sidebar.error(f"Failed to load SIM from path: {e}")


raw_md: Dict[str, Any] = st.session_state.raw_md
raw_sim: Dict[str, Any] = st.session_state.raw_sim


# -----------------------------
# Validate current configs
# -----------------------------
md_cfg: Optional[MDConfig] = None
sim_cfg: Optional[SimConfig] = None
md_error: Optional[str] = None
sim_error: Optional[str] = None

try:
    md_cfg = MDConfig.model_validate(raw_md)
except Exception as e:
    md_error = str(e)

try:
    sim_cfg = SimConfig.model_validate(raw_sim)
except Exception as e:
    sim_error = str(e)


# -----------------------------
# Main UI
# -----------------------------
tab_md, tab_sim, tab_export = st.tabs(["MD config", "SIM config", "Export / Save"])

# ---------- MD TAB ----------
with tab_md:
    st.subheader("Molecular Dynamics config (md_config.yaml)")

    sim = raw_md.get("simulation", {}) if isinstance(raw_md.get("simulation", {}), dict) else {}

    # Simulation basics
    c0, c1, c2, c3, c4 = st.columns(5)
    #sim["mode"] = int(c0.number_input("test", min_value=1, max_value=100, value=int(sim.get("test", 2))))
    sim["replicates"] = int(c1.number_input("replicates", min_value=1, max_value=100, value=int(sim.get("replicates", 2))))
    sim["time_ns"] = float(c2.number_input("time_ns", min_value=0.0001, value=float(sim.get("time_ns", 3.0))))
    sim["recording_interval_ps"] = float(c3.number_input("recording_interval_ps", min_value=0.001, value=float(sim.get("recording_interval_ps", 1.0))))
    sim["seed"] = int(c4.number_input("seed", min_value=0, value=int(sim.get("seed", 2025))))

    sim["mode"] = c0.selectbox(
    "Mode",
    options=["single_protein", "protein_small_molecule", "protein_protein"],
    index=["single_protein", "protein_small_molecule", "protein_protein"].index(sim.get("mode", "protein_small_molecule")),
)


    # Constraints (dt_fs needed for conversions)
    st.markdown("### constraints")
    constraints = sim.get("constraints", {}) if isinstance(sim.get("constraints", {}), dict) else {}
    cc1, cc2, cc3 = st.columns(3)
    constraints["dt_fs"] = float(cc1.number_input("dt_fs", min_value=0.1, value=float(constraints.get("dt_fs", 2.0))))
    constraints["bonds"] = cc2.selectbox("bonds", options=["None", "HBonds", "AllBonds"], index=["None", "HBonds", "AllBonds"].index(constraints.get("bonds", "HBonds")))
    constraints["rigid_water"] = cc3.checkbox("rigid_water", value=bool(constraints.get("rigid_water", True)))

    cc4, cc5, cc6 = st.columns(3)
    constraints["cutoff_nm"] = float(cc4.number_input("cutoff_nm", min_value=0.1, value=float(constraints.get("cutoff_nm", 1.0))))
    constraints["ewald_error_tolerance"] = float(cc5.number_input("ewald_error_tolerance", min_value=1e-8, value=float(constraints.get("ewald_error_tolerance", 0.0005))))
    constraints["constraint_tolerance"] = float(cc6.number_input("constraint_tolerance", min_value=1e-8, value=float(constraints.get("constraint_tolerance", 0.0001))))

    constraints["friction_ps"] = float(st.number_input("friction_ps", min_value=1e-6, value=float(constraints.get("friction_ps", 1.0))))
    sim["constraints"] = constraints

    dt_fs = float(constraints["dt_fs"])

    # Equilibration with steps<->time conversion
    st.markdown("### equilibration")
    eq = sim.get("equilibration", {}) if isinstance(sim.get("equilibration", {}), dict) else {}

    eq["protein_k"] = float(st.number_input("protein_k (harmonic)", min_value=0.0, value=float(eq.get("protein_k", 10.0))))

    entry_mode = st.radio(
        "Equilibration input mode",
        ["steps", "time (ps)"],
        horizontal=True,
        help="Enter equilibration lengths as steps OR time; values will be converted and saved as steps in YAML.",
    )

    def equil_row(label: str, key: str, default_steps: int) -> int:
        steps_val = int(eq.get(key, default_steps))
        time_ps_val = steps_to_ps(steps_val, dt_fs)

        cA, cB, cC = st.columns([1.2, 1, 1])
        cA.markdown(f"**{label}**")

        if entry_mode == "steps":
            steps_new = int(cB.number_input(f"{key} (steps)", min_value=0, value=int(steps_val), key=f"{key}_steps"))
            time_new = steps_to_ps(steps_new, dt_fs)
            cC.write(f"{time_new:.3f} ps")
            return steps_new
        else:
            time_new = float(cB.number_input(f"{key} (ps)", min_value=0.0, value=float(time_ps_val), key=f"{key}_ps"))
            steps_new = ps_to_steps(time_new, dt_fs)
            cC.write(f"{steps_new:,} steps")
            return steps_new

    eq["NVT_heating"] = equil_row("NVT heating", "NVT_heating", 5000)
    eq["NPT_equilibration"] = equil_row("NPT equilibration", "NPT_equilibration", 50000)
    eq["NPT_unrestrained"] = equil_row("NPT unrestrained", "NPT_unrestrained", 1000)

    sim["equilibration"] = eq

    # Forcefield
    st.markdown("### forcefield")
    ff = sim.get("forcefield", {}) if isinstance(sim.get("forcefield", {}), dict) else {}
    f1, f2, f3 = st.columns(3)
    ff["protein"] = f1.text_input("protein", value=str(ff.get("protein", "amber19-all.xml")))
    ff["water"] = f2.text_input("water", value=str(ff.get("water", "amber19/tip4pew.xml")))
    ff["watermodel"] = f3.text_input("watermodel", value=str(ff.get("watermodel", "tip4pew")))
    sim["forcefield"] = ff

    # System
    st.markdown("### system")
    sysc = sim.get("system", {}) if isinstance(sim.get("system", {}), dict) else {}
    s1, s2, s3, s4 = st.columns(4)
    sysc["salt_molar"] = float(s1.number_input("salt_molar", min_value=0.0, value=float(sysc.get("salt_molar", 0.0))))
    sysc["temperature_K"] = float(s2.number_input("temperature_K", min_value=1.0, value=float(sysc.get("temperature_K", 300.0))))
    sysc["pressure_atm"] = float(s3.number_input("pressure_atm", min_value=0.1, value=float(sysc.get("pressure_atm", 1.0))))
    sysc["barostat_interval_steps"] = int(s4.number_input("barostat_interval_steps", min_value=1, value=int(sysc.get("barostat_interval_steps", 25))))
    sim["system"] = sysc

    # Flexible binding pocket (simple editor)
    st.markdown("### flexible_binding_pocket")
    fbp = sim.get("flexible_binding_pocket", {}) if isinstance(sim.get("flexible_binding_pocket", {}), dict) else {}
    fbp["protein_k"] = float(st.number_input("flexible_binding_pocket.protein_k", min_value=0.0, value=float(fbp.get("protein_k", 5000.0))))

    # flexible_resids: show as table (resid:int, aa:str)
    flex = fbp.get("flexible_resids", {}) if isinstance(fbp.get("flexible_resids", {}), dict) else {}
    flex_df = pd.DataFrame([{"resid": int(k), "aa": str(v)} for k, v in flex.items()]).sort_values("resid") if flex else pd.DataFrame(columns=["resid", "aa"])
    flex_df = st.data_editor(
        flex_df,
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "resid": st.column_config.NumberColumn("resid", min_value=1, step=1),
            "aa": st.column_config.TextColumn("aa"),
        },
        key="flex_resids_editor",
    )
    # write back
    flex_clean = {}
    for _, row in flex_df.dropna().iterrows():
        try:
            flex_clean[int(row["resid"])] = str(row["aa"])
        except Exception:
            pass
    fbp["flexible_resids"] = flex_clean
    sim["flexible_binding_pocket"] = fbp

    # Metadynamics
    st.markdown("### metadynamics")
    meta = sim.get("metadynamics", {}) if isinstance(sim.get("metadynamics", {}), dict) else {}
    m1, m2, m3, m4, m5 = st.columns(5)
    meta["CV"] = m1.selectbox("CV", options=["COM", "contacts"], index=["COM", "contacts"].index(meta.get("CV", "COM")))
    meta["SIGMA"] = float(m2.number_input("SIGMA", min_value=1e-6, value=float(meta.get("SIGMA", 0.2))))
    meta["HEIGHT"] = float(m3.number_input("HEIGHT", min_value=1e-6, value=float(meta.get("HEIGHT", 0.03))))
    meta["PACE"] = int(m4.number_input("PACE", min_value=1, value=int(meta.get("PACE", 1000))))
    meta["STRIDE"] = int(m5.number_input("STRIDE", min_value=1, value=int(meta.get("STRIDE", 50))))
    sim["metadynamics"] = meta

    # Write back to raw_md
    raw_md["simulation"] = sim
    st.session_state.raw_md = raw_md

    # Validate + summary
    st.divider()
    if st.button("Validate MD config", type="primary"):
        try:
            md_cfg = MDConfig.model_validate(raw_md)
            st.success("MD config is valid ✔")
        except Exception as e:
            st.error(f"MD config invalid: {e}")

    if md_error:
        st.warning("Current MD config has validation issues.")
        st.code(md_error)


# ---------- SIM TAB ----------
with tab_sim:
    st.subheader("Simulation Info config (sim_config.yaml)")

    # mutations
    muts = raw_sim.get("mutations", ["WT"])
    if not isinstance(muts, list):
        muts = ["WT"]

    st.markdown("### mutations")
    muts_text = st.text_input("Mutations (comma-separated)", value=",".join([str(m) for m in muts]))
    mutations = [m.strip() for m in muts_text.split(",") if m.strip()]
    raw_sim["mutations"] = mutations

    # complexes
    st.markdown("### complexes")
    complexes = raw_sim.get("complexes", {})
    if not isinstance(complexes, dict):
        complexes = {}

    rows = []
    for name, entry in complexes.items():
        if not isinstance(entry, dict):
            continue
        rows.append(
            {
                "complex_name": name,
                "receptor": entry.get("receptor", ""),
                "ligand": entry.get("ligand", ""),
                "pdb": entry.get("pdb", ""),
                "sdf": entry.get("sdf", ""),
            }
        )

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["complex_name", "receptor", "ligand", "pdb", "sdf"])
    df = st.data_editor(
        df,
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "complex_name": st.column_config.TextColumn("complex_name", required=True),
            "receptor": st.column_config.TextColumn("receptor"),
            "ligand": st.column_config.TextColumn("ligand"),
            "pdb": st.column_config.TextColumn("pdb"),
            "sdf": st.column_config.TextColumn("sdf"),
        },
        key="complexes_editor",
    )

    # write back complexes dict
    new_complexes: Dict[str, Dict[str, str]] = {}
    for _, row in df.dropna(subset=["complex_name"]).iterrows():
        name = str(row["complex_name"]).strip()
        if not name:
            continue
        new_complexes[name] = {
            "receptor": str(row.get("receptor", "")).strip(),
            "ligand": str(row.get("ligand", "")).strip(),
            "pdb": str(row.get("pdb", "")).strip(),
            "sdf": str(row.get("sdf", "")).strip(),
        }

    raw_sim["complexes"] = new_complexes
    st.session_state.raw_sim = raw_sim

    st.divider()
    if st.button("Validate SIM config", type="primary"):
        try:
            sim_cfg = SimConfig.model_validate(raw_sim)
            st.success("SIM config is valid ✔")
        except Exception as e:
            st.error(f"SIM config invalid: {e}")

    if sim_error:
        st.warning("Current SIM config has validation issues.")
        st.code(sim_error)


# ---------- EXPORT TAB ----------
with tab_export:
    st.subheader("Export / Download / Save both YAMLs")

    md_yaml = dump_yaml(st.session_state.raw_md)
    sim_yaml = dump_yaml(st.session_state.raw_sim)

    cA, cB = st.columns(2, gap="large")

    with cA:
        st.markdown("### md_config.yaml")
        st.code(md_yaml, language="yaml")
        st.download_button("Download md_config.yaml", data=md_yaml.encode("utf-8"), file_name="md_config.yaml", mime="text/yaml")

        md_save_path = st.text_input("Save md_config.yaml to path", value=md_default_path, key="md_save_path")
        if st.button("Save md_config.yaml to path"):
            try:
                Path(md_save_path).parent.mkdir(parents=True, exist_ok=True)
                Path(md_save_path).write_text(md_yaml)
                st.success(f"Saved md_config.yaml to {md_save_path}")
            except Exception as e:
                st.error(f"Failed to save md_config.yaml: {e}")

    with cB:
        st.markdown("### sim_config.yaml")
        st.code(sim_yaml, language="yaml")
        st.download_button("Download sim_config.yaml", data=sim_yaml.encode("utf-8"), file_name="sim_config.yaml", mime="text/yaml")

        sim_save_path = st.text_input("Save sim_config.yaml to path", value=sim_default_path, key="sim_save_path")
        if st.button("Save sim_config.yaml to path"):
            try:
                Path(sim_save_path).parent.mkdir(parents=True, exist_ok=True)
                Path(sim_save_path).write_text(sim_yaml)
                st.success(f"Saved sim_config.yaml to {sim_save_path}")
            except Exception as e:
                st.error(f"Failed to save sim_config.yaml: {e}")