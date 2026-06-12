#!/usr/bin/env python3
"""
Simple one-file workflow for protein-ligand metadynamics in OpenMM.

Workflow summary
1) Import a prepared PDB structure (`structure_equilibrated.pdb` by default).
2) Select protein atoms from chain A and ligand atoms from chain X.
3) Build an OpenMM System.
   - Preferred: load a pre-parameterized system XML (`--system-xml`).
   - Fallback: parameterize directly from force field files (protein/water only).
4) Define one collective variable (distance between chain A and chain X centroids).
5) Run well-tempered metadynamics with default settings via `metadynamics.py`.
6) Save trajectory, log, checkpoint, and free-energy grid.

Example
python metadynamics_workflow.py \
  --pdb structure_equilibrated.pdb \
  --protein-chain A \
  --ligand-chain X \
  --output-dir meta_run \
  --system-xml ../NPT/system.xml

Notes
- If your ligand is non-standard and you do not provide `--system-xml`, force-field
  assignment may fail. In that case, use the system XML from a prior equilibration step.
- This script imports `Metadynamics` and `BiasVariable` from
  `/home/peter/tools/Funnel-Metadynamics/source/metadynamics.py` by default.
"""

import argparse
import os
import sys

import numpy as np
from openmm import CustomCentroidBondForce, LangevinMiddleIntegrator, Platform, XmlSerializer, unit
from openmm.app import (
    DCDReporter,
    ForceField,
    PDBFile,
    Simulation,
    StateDataReporter,
)


def _import_metadynamics_module(source_dir: str):
    """Import Metadynamics classes from the provided source directory."""
    source_dir = os.path.abspath(source_dir)
    if source_dir not in sys.path:
        sys.path.insert(0, source_dir)

    try:
        from metadynamics import BiasVariable, Metadynamics
    except ImportError as exc:
        raise ImportError(
            f"Could not import metadynamics.py from: {source_dir}. "
            "Set --meta-source to the folder containing metadynamics.py."
        ) from exc

    return Metadynamics, BiasVariable


def _select_chain_atom_indices(topology, chain_id: str):
    """Return atom indices belonging to one chain ID."""
    indices = []
    for chain in topology.chains():
        if chain.id == chain_id:
            for residue in chain.residues():
                for atom in residue.atoms():
                    indices.append(atom.index)
    return indices


def _load_system(args, pdb):
    """Load pre-parameterized system XML or attempt direct force-field build."""
    if args.system_xml:
        with open(args.system_xml, encoding="utf-8") as handle:
            return XmlSerializer.deserialize(handle.read())

    forcefield = ForceField(*args.forcefield)
    return forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=args.nonbonded_method,
        nonbondedCutoff=args.nonbonded_cutoff * unit.nanometer,
        constraints=args.constraints,
        rigidWater=True,
    )


def build_parser():
    parser = argparse.ArgumentParser(description="Simple metadynamics workflow for a protein-ligand complex.")
    parser.add_argument(
        "--pdb",
        default="structure_equilibrated.pdb",
        help="Input equilibrated PDB file.",
    )
    parser.add_argument(
        "--protein-chain",
        default="A",
        help="Protein chain ID.",
    )
    parser.add_argument(
        "--ligand-chain",
        default="X",
        help="Ligand chain ID.",
    )
    parser.add_argument(
        "--output-dir",
        default="metadynamics_output",
        help="Directory for outputs.",
    )
    parser.add_argument(
        "--meta-source",
        default="/home/peter/tools/Funnel-Metadynamics/source",
        help="Directory that contains metadynamics.py.",
    )

    parser.add_argument(
        "--system-xml",
        default=None,
        help="Optional pre-parameterized OpenMM system.xml.",
    )
    parser.add_argument(
        "--forcefield",
        nargs="+",
        default=["amber14-all.xml", "amber14/tip3pfb.xml"],
        help="OpenMM force field XML files used when --system-xml is not provided.",
    )
    parser.add_argument(
        "--platform",
        default="CUDA",
        choices=["CUDA", "CPU", "OpenCL", "Reference"],
        help="OpenMM platform.",
    )

    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep", type=float, default=0.002, help="ps")
    parser.add_argument("--steps", type=int, default=200000)

    # Default metadynamics settings mirror common values used in the project scripts.
    parser.add_argument("--bias-factor", type=float, default=10.0)
    parser.add_argument("--hill-height", type=float, default=1.5, help="kJ/mol")
    parser.add_argument("--hill-frequency", type=int, default=1000)
    parser.add_argument("--save-frequency", type=int, default=1000)

    # CV settings for protein-ligand COM distance.
    parser.add_argument("--cv-min", type=float, default=0.0, help="nm")
    parser.add_argument("--cv-max", type=float, default=3.0, help="nm")
    parser.add_argument("--cv-width", type=float, default=0.02, help="nm")
    parser.add_argument("--cv-grid", type=int, default=200)

    parser.set_defaults(
        nonbonded_method="PME",
        nonbonded_cutoff=1.0,
        constraints="HBonds",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    metadyn_cls, biasvar_cls = _import_metadynamics_module(args.meta_source)

    os.makedirs(args.output_dir, exist_ok=True)
    bias_dir = os.path.join(args.output_dir, "bias")
    os.makedirs(bias_dir, exist_ok=True)

    pdb = PDBFile(args.pdb)

    protein_atoms = _select_chain_atom_indices(pdb.topology, args.protein_chain)
    ligand_atoms = _select_chain_atom_indices(pdb.topology, args.ligand_chain)

    if not protein_atoms:
        raise ValueError(f"No atoms found in protein chain '{args.protein_chain}'.")
    if not ligand_atoms:
        raise ValueError(f"No atoms found in ligand chain '{args.ligand_chain}'.")

    system = _load_system(args, pdb)

    # 1D CV: centroid distance between protein chain A and ligand chain X.
    cv_force = CustomCentroidBondForce(2, "distance(g1,g2)")
    cv_force.addGroup(ligand_atoms)
    cv_force.addGroup(protein_atoms)
    cv_force.addBond([0, 1])
    cv_force.setUsesPeriodicBoundaryConditions(True)

    distance_cv = biasvar_cls(
        cv_force,
        0.0
        * unit.nanometer,  # args.cv_min *  the lower bound of the CV range, in nm. In your script it defaults to 0.0, so the bias starts at zero separation.
        3.0
        * unit.nanometer,  # args.cv_max the upper bound of the CV range, in nm. It defaults to 3.0, so the bias is tabulated up to 3 nm.
        0.02
        * unit.nanometer,  # args.cv_widththe Gaussian width used when adding metadynamics hills, in nm. It defaults to 0.02, so each deposited hill is fairly narrow.
        periodic=False,
        gridWidth=args.cv_grid,
    )

    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.timestep * unit.picoseconds,
    )

    platform = Platform.getPlatformByName(args.platform)
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    metad = metadyn_cls(
        system=system,
        variables=[distance_cv],
        temperature=args.temperature * unit.kelvin,
        biasFactor=args.bias_factor,
        height=args.hill_height * unit.kilojoules_per_mole,
        frequency=args.hill_frequency,
        saveFrequency=args.save_frequency,
        biasDir=bias_dir,
    )

    simulation.reporters.append(DCDReporter(os.path.join(args.output_dir, "trajectory.dcd"), 5000))
    simulation.reporters.append(
        StateDataReporter(
            os.path.join(args.output_dir, "state.log"),
            5000,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            speed=True,
            separator=",",
        )
    )

    metad.step(simulation, args.steps)

    simulation.saveCheckpoint(os.path.join(args.output_dir, "checkpoint.chk"))

    free_energy = metad.getFreeEnergy().value_in_unit(unit.kilojoules_per_mole)
    np.save(os.path.join(args.output_dir, "free_energy_kjmol.npy"), free_energy)

    with open(os.path.join(args.output_dir, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("Metadynamics run completed.\n")
        f.write(f"PDB: {args.pdb}\n")
        f.write(f"Protein chain: {args.protein_chain}\n")
        f.write(f"Ligand chain: {args.ligand_chain}\n")
        f.write(f"Steps: {args.steps}\n")
        f.write(f"Bias factor: {args.bias_factor}\n")
        f.write(f"Hill height (kJ/mol): {args.hill_height}\n")
        f.write(f"Hill frequency (steps): {args.hill_frequency}\n")


if __name__ == "__main__":
    main()
