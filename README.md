# squeezeMD – Molecular Dynamics Workflow + Analysis

squeezeMD is an automated molecular dynamics (MD) workflow built around Snakemake. It performs optional mutagenesis, runs OpenMM simulations, and generates standardized analysis outputs (trajectory statistics, interaction contacts, fingerprints, and visualization artifacts).

This repository contains the Snakemake pipeline, the MD/analysis scripts, and helper tools for configuration and visualization.

---

## Install

Follow the setup guide in [install/INSTALL.md](install/INSTALL.md).

---

## What the package does

### 1) Preprocessing / Mutagenesis
- Generates FoldX mutation files from a ligand chain and applies FoldX BuildModel for non‑WT variants.
- See: `mutate` command and the Snakemake `prep_mutagenesis` rule.

### 2) Molecular Dynamics (OpenMM)
- Protein–protein, protein–small‑molecule, and metadynamics modes.
- Includes restrained minimization, NVT heating, NPT equilibration, and production.
- Optional flexible binding pocket restraints and metadynamics (PLUMED).

### 3) Trajectory preparation
- Centers and aligns trajectories (H5MD → DCD) for analysis and visualization.

### 4) Analysis
- RMSF/RMSD, MD statistics plots, and B‑factor outputs.
- Posco contact analysis and summary plots (heatmaps + barplots).
- ProLIF interaction fingerprints across trajectories.
- Interaction surface visualization (PyMOL session + PNG + PDB with B‑factors).

### 5) Visualization
- Aligns final structures into a single PyMOL session.
- Creates interaction-surface sessions and images per complex/mutation.

---

## Workflow overview

The core pipeline is defined in the Snakemake [src/squeezemd/Snakefile](src/squeezemd/Snakefile). You run it via the `squeeze` wrapper, which locates the packaged Snakefile and forwards arguments to Snakemake.

Pipeline modes:
- `protein_protein` (alias: `PPi`)
- `protein_molecule` (alias: `molecule`)
- `metadynamics`
- `protein` (apo protein only)

Outputs are written per mode/complex/mutation/seed, plus summarized results in `results/`.

---

## Quick start

1) Create a `config/` directory in your working folder and add:
- `config/sim_config.yaml`
- `config/md_config.yaml`

2) Run the pipeline:

```bash
squeeze PPi --resources gpu=1 -j4
```

Dry‑run first if needed:

```bash
squeeze PPi --resources gpu=1 -j4 -n
```

The wrapper writes an `execute.sh` with the full Snakemake command for reproducibility.

---

## Configuration

### sim_config.yaml
Defines complexes and mutations for the Snakemake pipeline.

Expected structure (example):

```yaml
mutations:
  - WT
  - R65E
complexes:
  C1s_Gigastasin:
    receptor: C1s
    ligand: Gigastasin
    pdb: /abs/path/to/complex.pdb
    sdf: /abs/path/to/ligand.sdf   # required for protein_molecule mode
```
# modify pdb
pdb4amber -i input.pdb input.amber.pdb

rename resids in pymol
alter (chain A), resv += 432


### md_config.yaml
Defines MD protocol and system parameters (equilibration, forcefield, salt, temperature, recording interval, etc.).

You can generate both configs with the Streamlit app:

```bash
streamlit run src/squeezemd/streamlit/app.py
```

---

## Key CLI tools (installed entry points)

### Preprocessing
- `mutate` – create FoldX mutation files from ligand chain and mutation string.

### MD
- `run-md` – run OpenMM MD for a prepared structure.
- `center-traj` – center/align H5MD trajectory and export DCD.

### Trajectory analysis
- `explore-trajectory` – RMSF, RMSD, B‑factors, and MD statistics plots.
- `compute-rmsf` / `plot-rmsf` – aggregate RMSF across replicas.
- `analyze-potential-energy` – ligand potential energy breakdown (small molecule).

### Contact analysis
- `compute-posco-contacts` – run PoSCo on selected frames.
- `plot-contact-heatmap` / `plot-contact-barplot` – summarize contacts.

### ProLIF fingerprints
- `compute-proflif-fingerprints` – generate interaction fingerprints.
- `analyze-proflif-fingerprints` – aggregate and plot fingerprints.

### Visualization
- `visualize-interaction-surface` – generate PyMOL script, session, and PNG.

Most users should call these via the Snakemake workflow (`squeeze`) rather than running each script manually.

---

## Output structure (high level)

```
<mode>/<complex>/<mutation>/<seed>/
  MD/                     # raw and centered trajectories
  po-sco/                 # PoSCo contact data
  fingerprint/            # ProLIF outputs
  analysis/               # RMSF/RMSD/Stats
results/
  alignment/
  posco/
  rmsf/
  fingerprints/
  interactionSurface/
```

---

## External tools used

squeezeMD integrates several external tools. Most are installed via the conda environment in [install/INSTALL.md](install/INSTALL.md):

- OpenMM + OpenMMForceFields + OpenFF Toolkit
- MDAnalysis + MDTraj
- PoSCo (po‑sco)
- FoldX (mutagenesis)
- PLUMED (metadynamics)
- PyMOL (visualization)
- ProLIF (fingerprints)
- Aquaduct (optional channel analysis)

---

## Demo

There are demo folders under `demo/` for protein–protein workflows. From the repo root:

```bash
cd demo
squeeze PPi --resources gpu=1 -j4 -n
squeeze PPi --resources gpu=1 -j4
```

---

## Notes

- The Snakemake pipeline expects `config/sim_config.yaml` and `config/md_config.yaml` in the working directory.
- GPU execution is preferred; CPU fallback is supported but slower.
- Some analysis scripts assume chain A is the ligand and chain B/C are receptor chains.

---

## License

MIT