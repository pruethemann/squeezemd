# squeezeMD — Project Summary

A reference description of the squeezeMD molecular-dynamics workflow: what it does,
how to install and run it, the input/output structure, the configuration schema,
and current limitations.

---

## 1. Project purpose

squeezeMD is an automated, Snakemake-driven molecular-dynamics (MD) pipeline for
studying protein interactions. From a prepared structure it can:

- introduce in-silico point mutations (FoldX),
- run solvated OpenMM MD for protein–protein, protein–small-molecule, and
  metadynamics systems,
- and produce a standardized set of post-MD analyses (contacts, fingerprints,
  RMSF/RMSD, free energy, and PyMOL visualizations)

so that many complexes / mutations / replicates can be processed reproducibly and
compared.

---

## 2. Main features

- **Snakemake orchestration** — one workflow, multiple "modes" (targets); jobs fan
  out over receptor × ligand × mutation × seed.
- **Mutagenesis** — FoldX `BuildModel` driven from a parsed/validated mutation
  string (e.g. `R65E`, `R65E_Y117E`); `WT` is passed through unchanged.
- **OpenMM MD** — restrained minimization → staged NVT heating → NPT
  equilibration with tapering restraints → unrestrained NPT → production, with an
  optional well-tempered metadynamics stage (PLUMED).
- **Reproducibility controls** — a single configured seed drives the integrator,
  initial velocities, and barostat; force field, water model, temperature, salt,
  cutoffs, and (optionally) charge method / pH / box padding are read from config.
- **Analyses** — PoSCo contacts (parquet + heatmaps + barplots), ProLIF
  fingerprints, RMSF/RMSD/B-factors and MD statistics, metadynamics free-energy
  and convergence, interaction-surface PyMOL sessions, and structure alignment.
- **Traceability** — analysis outputs carry the originating complex / mutation /
  seed so results map back to the exact simulation.

---

## 3. Installation / setup

Linux (Ubuntu recommended) + an NVIDIA GPU for production MD. See
[install/INSTALL.md](install/INSTALL.md) for the full guide. In short:

```sh
git clone https://github.com/pruethemann/squeezemd.git
cd squeezemd
conda env create -f install/environment.yml      # creates the 'squeeze' env
conda activate squeeze
./install/install_bins_linux.sh                   # installs FoldX + PoSCo binaries
```

Verify:

```sh
python -m openmm.testInstallation
foldx_20261231 --version
po-sco --version
```

Developer install (editable + dev tools):

```sh
pip install -e ".[test]"     # ruff, pytest, coverage
```

---

## 4. Required dependencies

The conda environment (`install/environment.yml`) is the **source of truth** for
the full stack. Key components:

| Area | Packages |
|---|---|
| MD engine | `openmm`, `openmm-plumed`, `openmmforcefields`, `openff-toolkit` |
| Trajectory / analysis | `mdanalysis`, `mdtraj`, `prolif` |
| Workflow | `snakemake-minimal` |
| Data / plotting | `pandas`, `pyarrow`, `numpy`, `seaborn`, `plotly`, `matplotlib` |
| Visualization | `pymol-open-source` |
| External binaries | FoldX (mutagenesis), PoSCo (contacts) — installed via script |
| Optional | Aquaduct (water-channel analysis), packaged env `resources/aquaduct_env.yaml` |

`pyproject.toml` pins only the lightweight pip-installable core
(`typer`, `pyyaml`, `numpy`, `pandas`, `pyarrow`) — enough to import the pure
helpers and run the unit tests, but **not** a substitute for the conda env.

---

## 5. Input file structure

Run squeezeMD from a working directory that contains a `config/` folder and the
input structures it references:

```
<job-dir>/
  config/
    sim_config.yaml      # complexes, mutations, mode, resources
    md_config.yaml       # MD protocol + system parameters
  pdb/
    <receptor>.pdb       # prepared (e.g. pdb4amber) structure; chain A = ligand
  ligands/               # only for protein_molecule / metadynamics_molecule modes
    <ligand>.sdf
```

Conventions: **chain A is the ligand**; other protein chains are the receptor.
Small-molecule ligands are matched by residue name `UNK` (OpenFF).

See the worked examples under [demo/](demo/) (`T-1` protein–protein,
`T-2` metadynamics small-molecule).

---

## 6. Configuration options

### `config/sim_config.yaml`

```yaml
ID: T-1                       # short experiment id
name: Test_C1s_Gigastasin     # human-readable name
comment: ...                  # free text
mode: protein_protein         # protein_protein | protein_molecule |
                              #   metadynamics_ppi | metadynamics_molecule | protein
cpu: 4                        # cores passed to snakemake
gpu: 1                        # gpu resource budget
debug: True                   # use the fast test MD profile (md_test_config.yaml)

mutations:                    # one entry per variant; WT is passed through
  - WT
  - R65E

receptors:
  C1s:
    pdb: pdb/C1s-BD001.pdb
    # sdf_dir: ligands/        # optional override for small-molecule sdf location

ligands:                      # PPI: name of the chain-A partner; molecule: sdf basenames
  - Gigastasin

# Optional FoldX overrides (defaults shown):
# foldx:
#   binary: foldx_20261231
#   rotabase: ~/tools/foldX/foldx_Linux/rotabase.txt
```

### `config/md_config.yaml`

```yaml
simulation:
  replicates: 1
  time_ns: 0.1
  recording_interval_ps: 1
  number_frames: 2            # last-N frames used by contact/fingerprint analyses
  seed: 23                    # master seed (integrator + velocities + barostat)
  analysis_stride: 10         # (optional) keep every Nth frame in the analysis DCD

  equilibration:
    protein_k: 10             # restraint force constant (kJ/mol/nm²)
    NVT_heating: 10           # steps per heating step
    NPT_equilibration: 10
    NPT_unrestrained: 10

  forcefield:
    protein: amber19-all.xml
    water:   amber19/tip4pew.xml
    watermodel: tip4pew
    # small_molecule: openff-2.2.0          # (optional) OpenFF small-molecule FF
    # ligand_charge_method: am1bcc          # (optional) partial-charge method

  system:
    salt_molar: 0.15
    temperature_K: 300
    pressure_atm: 1
    barostat_interval_steps: 25
    # ph: 7.4                                # (optional) addHydrogens pH
    # box_padding_nm: 1.2                     # (optional) solvent padding

  constraints:
    dt_fs: 2
    bonds: HBonds
    rigid_water: True
    cutoff_nm: 1.0
    ewald_error_tolerance: 0.0005
    constraint_tolerance: 0.0001
    friction_ps: 1.0

  metadynamics:               # only used by the metadynamics_* modes
    SIGMA_COM: 0.08
    SIGMA_CONTACTS: 8
    CONTACT_R0: 0.45
    CONTACT_NN: 6
    CONTACT_ATOM_MODE: ca     # ca | heavy | backbone | all
    CONTACT_MAX_ATOMS_PER_PARTNER: 120
    HEIGHT: 0.08
    PACE: 3000
    STRIDE: 200
    # BIASFACTOR: 10          # (optional) well-tempered bias factor
```

The Snakefile calls `validate_config()` after merging the two files and exits with
a clear message if a required key (`mode`, `receptors`, `ligands`, `mutations`,
`simulation.*`) is missing or empty.

---

## 7. Example usage

From a job directory containing `config/`:

```sh
squeeze -n              # dry-run: show the DAG (target/mode read from config)
squeeze -j4             # run with 4 jobs (gpu/cpu budgets come from sim_config)
squeeze continue        # re-run the exact last command (saved in execute.sh)
squeeze upgrade         # reinstall the package from the git checkout
```

The wrapper reads `mode`/`gpu`/`cpu` from `config/sim_config.yaml`, locates the
packaged Snakefile, always passes `--rerun-incomplete`, and writes the exact
command to `execute.sh` for reproducibility. Extra arguments (`-n`, `-j`,
`--keep-going`, …) are forwarded to Snakemake.

Demo:

```sh
cd demo/T-1_Gigastasin
squeeze -n              # non-empty DAG for the C1s_Gigastasin protein–protein run
squeeze -j2             # run (debug: True uses the fast test MD profile)
```

---

## 8. Generated outputs

Per-run artifacts live under `<mode>/<receptor>_<ligand>/<mutation>/<seed>/`:

```
MD/
  structure_equilibrated.pdb        # equilibrated, solvated system
  structure_end.cif                 # final-frame topology (PDBx/mmCIF)
  trajectory_raw.h5                 # production trajectory (temporary)
  mdstats.csv                       # energies / temperature / volume time series
  center/                           # centered + aligned trajectory (DCD/H5) + topo
  energy/potential_energy.parquet   # ligand internal energy (molecule mode)
po-sco/posco_interaction.parquet    # per-frame PoSCo contacts (+ .txt dump)
fingerprint/fingerprint.parquet     # ProLIF interaction fingerprints
analysis/RMSF.html, RMSD.svg, RMSD.csv, bfactors.pdb, Stats.svg
metadynamics/                       # hills, Colvar, fes.dat, free-energy + convergence png
aquaduct/                           # optional water-channel analysis + PyMOL session
```

Aggregated results under `results/`:

```
results/
  alignment/align.pse               # all final structures aligned in one PyMOL session
  posco/                            # merged contacts parquet, heatmaps, barplots,
                                    #   per-complex interaction-surface PDB/PML/PSE/PNG
  rmsf/rmsf.parquet, rmsf.svg       # RMSF across replicates (coloured by mutation)
  fingerprints/                     # aggregated ProLIF interactions + plot
```

A consolidated `report.html` is generated on success (`squeeze` `onsuccess` hook).

---

## 9. Mutation workflow

1. The `prep_mutagenesis` rule runs per (mode, complex, mutation).
2. For `WT`, the receptor PDB is copied to `mutation.pdb` unchanged.
3. Otherwise, `mutate` extracts the chain-A ligand sequence
   (`extract_ligand_sequence`), validates the mutation string and builds the
   FoldX mutant file (`build_mutant_file_content` — raises on an unparseable
   mutation, out-of-range position, or wrong wild-type residue), and FoldX
   `BuildModel` produces the mutant structure. The rule verifies the FoldX output
   exists and fails loudly (printing the FoldX log) if it does not.

The FoldX binary name and rotabase path are configurable (`config['foldx']`).

---

## 10. MD workflow

`run-md` (one of the `md_*` rules per mode) builds and runs the system:

1. Solvate + ionize (`create_model_ppi` or `create_model_smallmolecule`); OpenFF
   parameterizes small molecules.
2. Add heavy-atom positional restraints; energy-minimize.
3. NVT heating ramp (50 → production temperature).
4. NPT equilibration with tapering restraints, then unrestrained NPT.
5. Optional well-tempered metadynamics (PLUMED) on COM distance + interface
   contacts.
6. Production run → HDF5 trajectory + `mdstats.csv` + final-frame CIF.

The integrator, initial velocities, and barostat are all seeded from the
configured master seed. `center-traj` then unwraps/centers/aligns the trajectory
and exports a (strided) DCD for analysis and visualization.

---

## 11. Post-MD analysis workflow

- **Contacts (PoSCo)** — `compute-posco-contacts` runs PoSCo on the last
  `number_frames` frames and writes a per-frame parquet; `merge_posco_contact`
  concatenates them; `plot-contact-barplot` / `plot-contact-heatmap` summarize
  per-residue interaction energies by type.
- **Fingerprints (ProLIF)** — `compute-prolif-fingerprints` then
  `analyze-prolif-fingerprints` aggregate interaction fingerprints across runs.
- **Trajectory** — `explore-trajectory` produces RMSF (per chain, with secondary
  structure), RMSD (ligand vs. all receptor chains), B-factors, and MD-statistics
  plots; `compute-rmsf` / `plot-rmsf` aggregate RMSF across replicates, coloured
  by mutation.
- **Metadynamics** — `plumed sum_hills` builds the FES; `analyze-collective-vars`,
  `analyze-metadynamics-hills`, and `analyze-metadynamics-convergence` plot the
  collective variables, free energy, and convergence.
- **Visualization** — `visualize-interaction-surface` writes per-residue energies
  into the B-factor column and renders a PyMOL session/PNG;
  `visualize_aligned_structures` aligns all final structures into one session.

---

## 12. Known limitations

- **GPU-centric, Linux only.** Production MD expects CUDA; CPU is a slow fallback
  for MD and for `analyze-potential-energy`. Only tested on Ubuntu.
- **Machine-specific bookkeeping.** `update_md_overview()` writes, and the
  `onsuccess` hook rsyncs to, a hard-coded `~/caracara/Squeeze/` path. The
  workflow will not run unmodified on a machine without that path. (Left as-is by
  request; see "future improvements".)
- **Metadynamics double-bias (open item).** The well-tempered PLUMED script in
  `md/metadynamics_auxillary.py` currently emits **two** `METAD` actions writing
  the same HILLS file (one plain, one well-tempered). This is flagged with a
  `FIXME` and **left unchanged pending author confirmation** — metadynamics
  free-energy results should be treated with care until it is resolved.
- **Charge-method mismatch (open item).** `analyze-potential-energy` parameterizes
  the ligand with `gasteiger` charges while production MD uses `am1bcc`; flagged
  with a `FIXME`, not silently changed.
- **Chain-A ligand assumption.** Most analyses assume chain A is the ligand and
  that residue numbering is consistent across runs.
- **Experimental code.** `src/squeezemd/toy/` and `auxillaryscripts/` are
  unmaintained scratch scripts, excluded from linting/tests (see their READMEs).
- **CI runs unit tests only** — it does not exercise a full MD run (no GPU).

---

## 13. Recommended future improvements

- Resolve the two flagged scientific items (single well-tempered `METAD`;
  consistent ligand charge method) after confirming intent.
- Make the `~/caracara` overview/rsync bookkeeping config-gated (off by default)
  so the pipeline is portable.
- Add `log:` / `benchmark:` directives to the long-running MD/analysis rules and
  ship a Snakemake profile for SLURM/HPC execution.
- Tighten wildcard constraints (`seed="\d+"`, mode alternation) as defensive
  hygiene, and add a JSON-schema validation of the config.
- Generalize the chain-A/ligand assumption and expose the remaining hard-coded MD
  parameters (already partly done) so a run is fully reconstructable from config.
- Record effective force-field / package / FoldX / PoSCo versions into each run's
  outputs for provenance.
