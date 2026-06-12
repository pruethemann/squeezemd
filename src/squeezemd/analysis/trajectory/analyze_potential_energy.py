#!/usr/bin/env python

"""Compute ligand potential energy components from an MD trajectory."""

import argparse

import mdtraj as md
import numpy as np
import pandas as pd
from openff.toolkit.topology import Molecule
from openmm import Context, OpenMMException, Platform, VerletIntegrator, unit
from openmm.openmm import System
from openmmforcefields.generators import SystemGenerator

from ...helper_functions import parse_run_metadata


def generate_ligand_system(ligand_path):
    """Build an OpenMM system for the ligand only (OpenFF parameters)."""
    ligand = Molecule.from_file(ligand_path)
    ligand_topology = ligand.to_topology().to_openmm()
    # FIXME(review): this uses 'gasteiger' charges, but the production MD
    # (run_md.create_model_smallmolecule) parameterizes the ligand with 'am1bcc'.
    # The ligand internal energies computed here are therefore not on the same
    # charge model as the trajectory. Left unchanged pending author confirmation.
    ligand.assign_partial_charges("gasteiger")

    protein_forcefield = "amber19-all.xml"  # params['simulation']['forcefield']['protein']
    water_model = "amber19/tip4pew.xml"  # params['simulation']['forcefield']['water']

    # 3. Use SystemGenerator to combine force fields
    generator = SystemGenerator(
        forcefields=[protein_forcefield, water_model],
        small_molecule_forcefield="openff-2.2.0",  # TODO: make sure to update to 3.0 if released soon
        molecules=[ligand],
        cache=None,
    )

    ligand_system = generator.create_system(ligand_topology)
    return ligand_system


def assign_force_groups(system):
    """
    Assign each Force in the system to its own force group (0..31).
    Returns a dict: name -> group_index.
    """
    group_map = {}
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
        # Ignore CMmotion
        if force.__class__.__name__.startswith("CM"):
            continue
        group_map[force.__class__.__name__] = i
    return group_map


def compute_potential_energy(
    traj_h5: str, top: str, ligand_selection: str, ligand_system: System, group_map
) -> np.ndarray:
    """
    Compute ligand internal potential energy per frame using a provided ligand-only OpenMM System.

    Parameters
    ----------
    traj_h5 : str
        Path to H5MD trajectory.
    top : str
        Path to topology file (PDB/mmCIF) matching traj.
    ligand_selection : str
        MDTraj selection string (e.g. "resname UNK" or "chainid 0").
    ligand_system : openmm.System
        System containing ONLY ligand atoms, in the SAME order as ligand_selection returns.
        Must match the ligand-only topology used below.

    Returns
    -------
    energies_kj_mol : np.ndarray
        Potential energies in kJ/mol for each frame.
    """
    # Load trajectory and select ligand
    traj = md.load(traj_h5, top=top)
    lig_idx = traj.topology.select(ligand_selection)

    if lig_idx.size == 0:
        raise ValueError(f"No atoms matched ligand_selection='{ligand_selection}'")

    lig_traj = traj.atom_slice(lig_idx)
    lig_positions_nm = lig_traj.xyz  # shape (n_frames, n_atoms, 3) in nm

    # Integrator is not used, but OpenMM requires an integrator instance
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)

    # Prefer the GPU but fall back to CPU so this analysis also runs on CPU-only hosts.
    try:
        platform = Platform.getPlatformByName("CUDA")
    except OpenMMException:
        print("ATTENTION: No CUDA GPU detected. Computing potential energy on CPU.")
        platform = Platform.getPlatformByName("CPU")
    context = Context(ligand_system, integrator, platform)

    # Allocate one array per force-group term (plus the total). Building this from
    # group_map avoids a duplicate-key bug and silently dropping force types such
    # as UreyBradley / CMAP that the ligand force field may contain.
    energies = {name: np.empty(lig_traj.n_frames, dtype=float) for name in (*group_map.keys(), "potential")}

    # Calculate the potential energy for every frame
    for i in range(lig_traj.n_frames):
        # Determine atom positions of seleciton in paricular frame
        context.setPositions(lig_positions_nm[i] * unit.nanometer)

        # Extract all energy terms individually
        for ene_name, g in group_map.items():
            state = context.getState(getEnergy=True, groups=(1 << g))
            energies[ene_name][i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

        # Extract total energy: sum of all components
        state = context.getState(getEnergy=True)
        energies["potential"][i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    return energies


def parse_arguments():
    """Parse CLI arguments for potential energy analysis."""
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("--topo")
    parser.add_argument("--traj")
    parser.add_argument("--sdf")
    parser.add_argument("--selection", default="resname UNK", help="MDTraj selection for the ligand")
    parser.add_argument("--config", default=None, help="MD config (reserved; currently unused)")

    # Output
    parser.add_argument("--energy")
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Build the ligand-only OpenMM system and assign per-term force groups.
    ligand_system = generate_ligand_system(args.sdf)
    group_map = assign_force_groups(ligand_system)

    # Compute the energies from the trajectory
    energy = compute_potential_energy(
        traj_h5=args.traj,
        top=args.topo,
        ligand_selection=args.selection,
        ligand_system=ligand_system,
        group_map=group_map,
    )

    data_df = pd.DataFrame(energy)
    data_df["frame"] = data_df.index

    # Tag with the run identity so the energy table is traceable to its simulation.
    for key, value in parse_run_metadata(args.topo).items():
        data_df[key] = value

    data_df.to_parquet(args.energy)


if __name__ == "__main__":
    main()
