#!/usr/bin/env python
import numpy as np
import mdtraj as md
from openmm import unit, Platform, Context, app, VerletIntegrator
from openmm.openmm import System
import argparse
import pandas as pd
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule

def assign_force_groups(system):
    """
    Assign each Force in the system to its own force group (0..31).
    Returns a dict: name -> group_index.
    """
    group_map = {}
    for i, force in enumerate(system.getForces()):
        print(force)
        force.setForceGroup(i)
        # Ignore CMmotion
        if force.__class__.__name__.startswith('CM'):
            continue
        group_map[force.__class__.__name__] = i
    return group_map

def compute_potential_energy(
    traj_h5: str,
    top: str,
    ligand_selection: str,
    ligand_system: System,
    group_map) -> np.ndarray:
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
    # Load Traj and select ligand
    traj = md.load(traj_h5, top=top)
    lig_idx = traj.topology.select(ligand_selection)

    if lig_idx.size == 0: raise ValueError(f"No atoms matched ligand_selection='{ligand_selection}'")

    lig_traj = traj.atom_slice(lig_idx)
    lig_positions_nm = lig_traj.xyz  # shape (n_frames, n_atoms, 3) in nm

    # Integrator is NOT used. But OpenMM requires integrator for API
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)

    #platform = Platform.getPlatformByName('CUDA')
    platform = Platform.getPlatformByName('CUDA')
    context = Context(ligand_system, integrator, platform)

    energies = {'potential':np.empty(lig_traj.n_frames, dtype=float),
                'NonbondedForce':np.empty(lig_traj.n_frames, dtype=float),
                'HarmonicBondForce':np.empty(lig_traj.n_frames, dtype=float),
                'PeriodicTorsionForce':np.empty(lig_traj.n_frames, dtype=float),
                'PeriodicTorsionForce':np.empty(lig_traj.n_frames, dtype=float),
                'HarmonicAngleForce':np.empty(lig_traj.n_frames, dtype=float),
    }

    # Calculate the potential energy for every frame
    for i in range(lig_traj.n_frames):
        # Determine atom positions of seleciton in paricular frame
        context.setPositions(lig_positions_nm[i] * unit.nanometer)

        # Extract all energy terms individually
        for ene_name, g in group_map.items():
            state = context.getState(getEnergy=True, groups=(1 << g))
            energies[ene_name][i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)       

        # Extract total energy: Sum of the above
        state = context.getState(getEnergy=True)
        energies['potential'][i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    return energies

def generate_ligand_system(ligand_path):
    ligand = Molecule.from_file(ligand_path)

    ligand_topology = ligand.to_topology().to_openmm()
    ligand_positions = ligand.conformers[0].to_openmm()

    ligand.assign_partial_charges('gasteiger')   

    ff_kwargs = {
        'constraints':app.HBonds,
        'rigidWater': True,# TODO standardize with yaml
        'ewaldErrorTolerance':0.0001 # params['simulation']['constraints']['ewald_error_tolerance']
    }
    periodic_forcefield_kwargs = {
        'nonbondedMethod': app.PME,
        'nonbondedCutoff': 1.0 #params['simulation']['constraints']['cutoff_nm'] * nanometers
    }

    protein_forcefield = "amber19-all.xml" # params['simulation']['forcefield']['protein']
    water_model = "amber19/tip4pew.xml"#params['simulation']['forcefield']['water']

    # 3. Use SystemGenerator to combine force fields
    generator = SystemGenerator(
        forcefields=[protein_forcefield, water_model],
        small_molecule_forcefield="openff-2.2.0",           # TODO: make sure to update to 3.0 if released soon
        molecules=[ligand],
        cache=None,
        forcefield_kwargs=ff_kwargs,
        periodic_forcefield_kwargs=periodic_forcefield_kwargs
    )

    ligand_system = generator.create_system(ligand_topology)
    return ligand_system

def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo')
    parser.add_argument('--traj')
    parser.add_argument('--sdf')
    parser.add_argument('--selection', default='resname UNK')
    parser.add_argument('--config', default='resname UNK')

    # Output
    parser.add_argument('--energy')
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    # define the system if energies
    ligand_system = generate_ligand_system(args.sdf)

    # Define bonded and non-bonded energy terms
    group_map = assign_force_groups(ligand_system)

    # Compute the energies from the trajectory
    energy = compute_potential_energy(args.traj, args.topo,args.selection,ligand_system,group_map)

    data_df = pd.DataFrame(energy)
    data_df['frame'] = data_df.index

    print(data_df)
    data_df.to_parquet(args.energy)
