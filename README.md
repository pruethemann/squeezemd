# squeezeMD – A Comprehensive Molecular Dynamics Analysis Workflow

[![PyPI version](https://img.shields.io/pypi/v/squeezemd.svg)](https://pypi.python.org/pypi/squeezemd)
[![Documentation Status](https://readthedocs.org/projects/squeezemd/badge/?version=latest)](https://squeezemd.readthedocs.io/en/latest/?version=latest)

* MIT license
* Documentation: [squeezemd.readthedocs.io](https://squeezemd.readthedocs.io)

---

## Install

Please follow the instructions in [install/INSTALL.md](install/INSTALL.md).

---

## Summary

*squeezeMD* provides an integrated solution for comprehensive molecular dynamics (MD) analysis.
It includes functionality for:

* in-silico mutagensis
* MD simulations
* Explorative trajectory analysis
* Interaction fingerprinting analysis
* Visualisation of the interaction surface


The workflow streamlines the analysis of complex MD simulations, enabling detailed examination of molecular interactions, stability, and conformational changes.

---

## Detailed Summary

### Snakefile

Serves as the backbone of the workflow, orchestrating the execution of analysis scripts.
Ensures pipeline steps are executed efficiently and in the correct order.

### Mutagenesis

Performs mutation analysis, characterizing structural and functional effects of amino acid changes.

### MD

Runs molecular dynamics simulations, revealing conformational changes and dynamics of molecules.

### Explorative Trajectory Analysis

Provides tools for trajectory inspection, helping identify key molecular events and interactions.

### Centering Trajectories

Centers and aligns trajectories, enabling accurate comparisons across simulations.

### Interaction Analysis

Analyzes molecular interactions at atom/residue level and extracts global interaction patterns.

### Interaction Fingerprints

Generates and compares interaction fingerprints, summarizing interaction motifs across simulations.

### Interaction Surface

Analyzes binding and interaction surfaces, critical for studying molecular recognition.

---

## Demo Workflow

```bash
cd demo

# Perform a dry run
squeeze PPi --resources gpu=1 -j4 -n

# Perform the demo production run
squeeze PPi --resources gpu=1 -j4
```


## Infos

* [Python Packaging & CLI](https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html)
* [GitHub workflow for PyPI](https://github.com/pypa/packaging.python.org/blob/main/source/guides/github-actions-ci-cd-sample/publish-to-test-pypi.yml)

---

## Upload to pyPi

```bash
python -m build && pip install --upgrade .
twine upload dist/*
# Username: __token__
# Password: PyPI token
```

---

## Prepare AMBER PDBs

*in progress*

1. Prepare protein in Maestro
2. Convert with `pdb4amber`:

   ```bash
   pdb4amber -i input.pdb -o input.amber.pdb
   ```
3. Example PyMOL commands:

   ```python
   alter (chain A), chain='B'
   alter (chain B), resi=str(int(resi)+315)
   alter (chain C), resi=str(int(resi)+51)
   ```

---

