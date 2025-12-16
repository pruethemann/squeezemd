#!/usr/bin/env python
"""
    Molecular Dynamics Workflow using OpenMM
    -------------------------------------------------
    Updated protocol (fixed heating instability):
    1. Add positional restraints before minimization.
    2. Energy minimization with heavy atoms restrained.
    3. Gradual heating under NVT using smooth ramp (50→300 K).
    4. Switch to NPT, gradually reduce restraints (10→5→1 kcal/mol/Å²).
    5. Unrestrained NPT equilibration.
    6. Production run.

    Designed for solvated protein(-protein) complexes.
"""

import argparse, os
from openmm.unit import *
from openmm import app, OpenMMException, Platform, LangevinMiddleIntegrator, MonteCarloBarostat, CustomExternalForce
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import mdtraj
import mdtraj.reporters
from openmmplumed import PlumedForce
import MDAnalysis as mda
import numpy as np
from Helper import import_yaml
from openmm.unit import kilojoule_per_mole,  nanometer


def add_positional_restraints(system, topology, positions, k=10.0):
    """
    Add harmonic restraints to heavy atoms (kcal/mol/Å²).
    Applied to all non-solvent heavy atoms.
    """
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    
    system.addForce(restraint)

    restraint.addGlobalParameter('k', k*kilojoule_per_mole/nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in topology.atoms():
        resname = atom.residue.name
        # Exclude water and ions from the restraints
        if resname not in ('HOH', 'Na+', 'Cl-') and atom.element.symbol != 'H':
            restraint.addParticle(atom.index, positions[atom.index])

    return system, restraint


def define_platform():
    """
    Detect NVIDIA GPU (CUDA) or fallback to CPU.
    """
    try:
        return Platform.getPlatformByName('CUDA')
    except OpenMMException:
        print("ATTENTION: No CUDA GPU detected. Running on CPU.")
        return Platform.getPlatformByName('CPU')


def energy_minimisation(simulation):
    """Run energy minimization and print energy difference."""
    e_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy()
    e_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print('Energy difference (minimization):', e_before - e_after)


def create_model_smallmolecule(modeller, salt_concentration, params, sdf):
    """
    Build solvated system with ions.
    This does only work for protein protein interaction. See legacy MD for small molecules
    """
    protein_forcefield = params['simulation']['forcefield']['protein']
    water_model = params['simulation']['forcefield']['water']

    ligand = Molecule.from_file(sdf)
    if not ligand.conformers:
        raise ValueError("Ligand SDF has no 3D conformers – please provide a 3D SDF.")
    
    #is this necessary? ligand.assign_partial_charges('gasteiger')   

    ligand_topology = ligand.to_topology().to_openmm()
    ligand_positions = ligand.conformers[0].to_openmm()

    ff_kwargs = {
        'constraints':app.HBonds,
        'rigidWater': True,                     # Allows to increase step size to 4 fs
        'ewaldErrorTolerance':params['ewaldErrorTolerance']
    }
    periodic_forcefield_kwargs = {
        'nonbondedMethod':app.PME,
        'nonbondedCutoff':params['nonbondedCutoff'] * nanometers
    }

    #'removeCMMotion': False

    # 3. Use SystemGenerator to combine force fields
    generator = SystemGenerator(
        forcefields=[protein_forcefield, water_model],
        small_molecule_forcefield='openff-2.0.0',
        molecules=[ligand],
        cache=None,
        forcefield_kwargs=ff_kwargs,
        periodic_forcefield_kwargs=periodic_forcefield_kwargs
    )

    # Add ligand to modell
    modeller.add(ligand_topology, ligand_positions)

    modeller.addHydrogens(generator.forcefield)       # TODO: Check whether His protonation states are changed
    #modeller.addExtraParticles(forcefield_generated.forcefield)          # Required for tip4p (orbital)

    # Add solvent
    modeller.addSolvent(generator.forcefield,
                        model=params['simulation']['forcefield']['watermodel'],                
                        boxShape='cube',
                        ionicStrength=salt_concentration,
                        positiveIon='Na+',
                        negativeIon='Cl-',
                        neutralize=True,
                        padding=1.2 * nanometers)
    
    # Create the MD system
    # TODO: add constraints before
    system = generator.create_system(modeller.topology)
    return system

def create_model_ppi(modeller, salt_concentration, params):
    """
    Build solvated system with ions.
    This does only work for protein protein interaction. See legacy MD for small molecules
    """

    protein_forcefield = params['simulation']['forcefield']['protein']
    water_model = params['simulation']['forcefield']['water']

    print(f'Initializing ForceField: {protein_forcefield} + {water_model}')
    forcefield = app.ForceField(protein_forcefield, water_model)

    modeller.addHydrogens(forcefield)       # TODO: Check whether His protonation states are changed
    modeller.addExtraParticles(forcefield)  # Required for tip4p (orbital)

    # Add solvent
    modeller.addSolvent(forcefield,
                        model=params['simulation']['forcefield']['watermodel'],                
                        boxShape='cube',
                        ionicStrength=salt_concentration,
                        positiveIon='Na+',
                        negativeIon='Cl-',
                        neutralize=True,
                        padding=1.2 * nanometers)
    
    # Create the MD system
    system = forcefield.createSystem(modeller.topology,
                             nonbondedMethod=app.PME,
                             nonbondedCutoff=params['simulation']['constraints']['cutoff_nm'] * nanometers,
                             constraints=app.HBonds,
                             rigidWater=params['simulation']['constraints']['rigid_water'],
                             ewaldErrorTolerance=params['simulation']['constraints']['ewald_error_tolerance'])
    return system

def save_cif(simulation, cif_file: os.path):
    positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    with open(cif_file, "w") as f:
        app.PDBxFile.writeFile(simulation.topology, positions, f, keepIds=True)

def save_pdb(simulation, pdb_file:os.path):
    positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    with open(pdb_file, "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions,f, keepIds=True)

def extract_atom_indices(pdf_file: os.path, cutoff = 5.0):

    u = mda.Universe(pdf_file)

    # get all heavy atoms of lig and rec
    lig_heavy = u.select_atoms("chainID A and not name H*")
    rec_heavy = u.select_atoms("(chainID B or chainID C) and not name H*")

    # get all heavy atoms of lig and rec
    lig = u.select_atoms("chainID A")
    rec = u.select_atoms("chainID B or chainID C")

    # Compute distance matrix between all atoms of the two chains
    dist = mda.lib.distances.distance_array(lig_heavy.positions, rec_heavy.positions)

    # Boolean masks of interface atoms
    lig_interface_mask = np.any(dist < cutoff, axis=1)
    rec_interface_mask = np.any(dist < cutoff, axis=0)

    # Interface atoms selections
    lig_interface = lig_heavy[lig_interface_mask]
    rec_interface = rec_heavy[rec_interface_mask]

    # Print in PLUMED-friendly format. Add +1 because plumed starts at atom id 1 and not 0
    lig_plumed = ",".join(map(str, lig.indices + 1))
    rec_plumed = ",".join(map(str, rec.indices + 1))
    lig_interface_plumed = ",".join(map(str, lig_interface.indices + 1))
    rec_interface_plumed = ",".join(map(str, rec_interface.indices + 1))

    atom_indices = {'lig_index': lig_plumed,
                    'rec_index': rec_plumed,
                    'lig_interface_index':lig_interface_plumed,
                    'rec_interface_index':rec_interface_plumed,
                    'lig_min':lig.indices.min() + 1,
                    'lig_max':lig.indices.max() + 1,
                    'rec_min':rec.indices.min() + 1,
                    'rec_max':rec.indices.max() + 1,
    }

    return atom_indices

def add_metadynamics_forces_centerofmass(params, system):
    # Metadynamics params
    sigma = params['simulation']['metadynamics']['SIGMA']
    height = params['simulation']['metadynamics']['HEIGHT']
    pace = params['simulation']['metadynamics']['PACE']
    stride = params['simulation']['metadynamics']['STRIDE']

    # Get absolute paths for outputs
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # get relevant atom indexes
    id = extract_atom_indices(args.equilibrated)

    script = f"""
            # get residue and chainID information
            MOLINFO STRUCTURE={args.equilibrated}

            # Define two groups (ligand:Entity0 and receptor:Entity1)
            WHOLEMOLECULES ENTITY0={id['lig_min']}-{id['lig_max']} ENTITY1={id['rec_min']}-{id['rec_max']}

            # Group heavy atoms for contact
            grp_lig: GROUP ATOMS={id['lig_min']}-{id['lig_max']}
            grp_rec: GROUP ATOMS={id['rec_min']}-{id['rec_max']}

            # Define center of mass of the two partners
            lig: COM ATOMS=grp_lig
            rec: COM ATOMS=grp_rec

            # Distance between the two COMs (in nm)
            d1: DISTANCE ATOMS=lig,rec

            METAD ARG=d1 SIGMA={sigma} HEIGHT={height} PACE={pace} FILE={hills_path}
            PRINT ARG=d1 STRIDE={stride} FILE={colvar_path}
            """
    
    print(script)

    plumed = PlumedForce(script)
    plumed.setTemperature(T*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system


def add_metadynamics_forces_singledistance(metadynamics_params, T:int, system):
    """
    in progress
    """
    # Example: add a distance-based collective variable
    atomId1 = metadynamics_params[0]['d1'][0]['atomId1'] + 1
    atomId2 = metadynamics_params[0]['d1'][1]['atomId2'] + 1

    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    script = f"""
            d1: DISTANCE ATOMS={atomId1},{atomId2}
            METAD ARG=d1 SIGMA=0.1 HEIGHT=0.3 PACE=50 FILE={hills_path}
            PRINT ARG=d1 STRIDE=50 FILE={colvar_path}
            """
    plumed = PlumedForce(script)
    plumed.setTemperature(T*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system

# ---------------------------
# Simulation procedure
# ---------------------------

def simulate(args, params):
    """
    Set up and start the simulation
    """
    # Detect GPU
    platform = define_platform()

    # Load structure
    protein = app.PDBFile(args.pdb)
    modeller = app.Modeller(protein.topology, protein.positions)

    salt_concentration = params['simulation']['system']['salt_molar'] * molar

    # Create solvated system depending on whether ligand is small molecule or protein
    if args.sdf is None: # ligand is protein
        system = create_model_ppi(modeller, salt_concentration, params)
    else: # ligand is small molecule
        system = create_model_smallmolecule(modeller, salt_concentration, params, args.sdf)
        
    # Add restraints BEFORE minimization
    system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=10.0)

    # Integrator setup
    dt_fs = params['simulation']['constraints']['dt_fs']
    temperature = params['simulation']['system']['temperature_K'] * kelvin
    friction = params['simulation']['constraints']['friction_ps'] / picoseconds
    
    integrator = LangevinMiddleIntegrator(temperature, friction, dt_fs * femtoseconds)
    integrator.setConstraintTolerance(params['simulation']['constraints']['constraint_tolerance'])
    integrator.setRandomNumberSeed(args.seed)

    # Set up the simulation. Add the integrator and the position
    properties = {"Precision": "mixed", "DeterministicForces": "true"}
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # ---------------------
    # Stage 0: Minimization
    # ---------------------
    print('\n=== Stage 0: Energy minimization with restraints ===')
    energy_minimisation(simulation)

    # ---------------------
    # Stage 1: NVT Heating with restrains
    # ---------------------
    print('\n=== Stage 1: Smooth NVT heating ===')
    simulation.context.setVelocitiesToTemperature(50 * kelvin)

    # Smooth temperature ramp 50 → 300 K
    temp_steps = [50, 100, 150, 200, 250, 300, params['simulation']['system']['temperature_K']]
    steps_per_temp = params['simulation']['equilibration']['NVT_heating']  # e.g. 5000 = 10 ps
    for T in temp_steps:
        print(f" → Heating to {T} K ...")
        simulation.integrator.setTemperature(T * kelvin)
        simulation.step(steps_per_temp)

    # ---------------------
    # Stage 2: NPT equilibration with tapering restraints
    # ---------------------
    print('\n=== Stage 2: NPT equilibration with tapering restraints ===')
    
    pressure = params['simulation']['system']['pressure_atm'] * atmospheres
    barostat_interval_steps = params['simulation']['system']['barostat_interval_steps']

    barostat = MonteCarloBarostat(pressure, temperature, barostat_interval_steps)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)

    # Define tapering schedule for restraints (kcal/mol/Å²)
    for k in [5.0, 1.0]:
        print(f"Tapering restraints to {k} kcal/mol/Å²")
        # Update global k parameter (not per particle!)
        restraint_force.setGlobalParameterDefaultValue(0, k * kilojoule_per_mole / nanometer**2)
        restraint_force.updateParametersInContext(simulation.context)

        # Run equilibration for each step
        simulation.step(params['simulation']['equilibration']['NPT_equilibration'])

    # ---------------------
    # Stage 3: Unrestrained NPT at simulation T
    # ---------------------
    print('\n=== Stage 3: Unrestrained NPT equilibration ===')
    # remove the CustomForce which restrained the system
    for i, force in enumerate(system.getForces()):
        if force.__class__.__name__ == restraint_force.__class__.__name__:
            system.removeForce(i)
            break

    simulation.context.reinitialize(preserveState=True)
    simulation.step(params['simulation']['equilibration']['NPT_unrestrained'])

    # Save equilibrated pdb
    save_pdb(simulation, args.equilibrated)

    # ---------------------
    # Stage 4: Metadynamics (optional)
    # ---------------------

    if params['simulation']['metadynamics']['enabled']:
        print(f'\n=== Stage 4: Initiate Metadynamics')
        simulation.system = add_metadynamics_forces_centerofmass(params, simulation.system)
        simulation.context.reinitialize(preserveState=True)  # keep positions/velocities

    # ---------------------
    # Stage 5: Production
    # ---------------------
    print(f'\n=== Stage 5: Production run ({params['simulation']['time_ns']} ns) ===')
    recordInterval = int(params['simulation']['recording_interval_ps'] * 1000 / params['simulation']['time_ns'])
    total_steps=params['simulation']['time_ns'] * 1e6 / dt_fs
    
    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj, recordInterval)
    dataReporter = app.StateDataReporter(
        args.stats, recordInterval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True,
        volume=True, density=True, progress=True, remainingTime=True,
        speed=True, totalSteps=total_steps)
    
    simulation.reporters.append(HDF5Reporter)
    simulation.reporters.append(dataReporter)

    simulation.step(total_steps)

    # Save the end file as cif
    save_cif(simulation, args.topo_cif)


# ---------------------------
# Argument parsing
# ---------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    # Input
    parser.add_argument('--pdb', default='input/fix1.pdb')
    parser.add_argument('--md_settings', default='input/params.yml')
    parser.add_argument('--seed', type=int, default=12)

    # Output
    parser.add_argument('--equilibrated', default='output/equilibrated.pdb', help="Equilibrated system in water box")
    parser.add_argument('--topo_cif', default='output/top.cif', help="Last uncentered frame of MD")
    parser.add_argument('--traj', default='output/traj.h5', help="h5md trajectory")
    parser.add_argument('--stats', default='output/stats.txt', help="Molecular dynamics statistics and progress")
    parser.add_argument('--metadynamics_hills', default="output/metadynamics_hills.txt", help='Metadynamics hill output file.')
    parser.add_argument('--metadynamics_colvar', default="output/metadynamics_colvar.txt", help='Metadynamics output file.')
    parser.add_argument('--sdf', required=False,help='Small molecule sdf file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    params = import_yaml(args.md_settings)
    simulate(args, params)
