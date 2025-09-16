#!/usr/bin/env python
"""
    Module performs Molecular Dynamics using OpenMM
    - Import of amber prepared pdb structure
    - Adding water box with 150 mM NaCl
    - Performs MD

NOTES:
    A seed is set for the integrator and the initial velocities

Updated protocol:
1. Energy minimization with heavy atoms fixed.
2. Gradual heating under NVT with strong restraints.
3. Switch to NPT, gradually reduce restraints (10 → 5 → 1 kcal/mol/Å²).
4. Unrestrained NPT equilibration.
5. Production run.

This follows modern AMBER/CHARMM best practices for solvated proteins.
"""

import argparse, os
from openmm.unit import nanometers, kelvin, femtoseconds, picoseconds, atmospheres, molar, kilojoule_per_mole
from openmm import app, OpenMMException, Platform, LangevinMiddleIntegrator, MonteCarloBarostat, CustomExternalForce
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import mdtraj
import mdtraj.reporters
from openmmplumed import PlumedForce
from Helper import import_yaml, save_yaml

# ---------------------------
# Restraint helper
# ---------------------------
def add_positional_restraints(system, topology, positions, k=10.0):
    """
    Add harmonic restraints to heavy atoms (kcal/mol/Å²).
    Required for NVT equilibration.
    """
    force = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for atom in topology.atoms():
        if atom.element.symbol != 'H':  # heavy atoms only
            pos = positions[atom.index]
            force.addParticle(atom.index, [k, pos.x, pos.y, pos.z])

    system.addForce(force)
    return system, force


# ---------------------------
# GPU Platform detection
# ---------------------------
def define_platform():
    """
    Detect NVIDIA GPU with CUDA, fallback to CPU if not available.
    """
    try:
        return Platform.getPlatformByName('CUDA')
    except OpenMMException:
        print("ATTENTION: no CUDA driver or GPU detected. Simulation runs on CPU")
        return Platform.getPlatformByName('CPU')

# ---------------------------
# Parameter handling
# ---------------------------
def set_parameters(params):
    global nonbondedCutoff, ewaldErrorTolerance, constraintTolerance, temperature
    global dt, recordInterval, friction, pressure, constraint, barostatInterval, platform
    global ff_kwargs

    # Physical parameters
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = 0.00001
    temperature = params['temperature'] * kelvin  # physiological temperature

    # Time parameters
    args.steps = int(params['time'] * 1e6 / params['dt'])
    args.time = params['time']
    args.recorded_steps = int(params['time'] * 1000 / params['recordingInterval'])
    dt = params['dt'] * femtoseconds
    recordInterval = args.steps * params['recordingInterval'] // (params['time'] * 1000)

    # Constraints
    friction = 1.0 / picoseconds
    pressure = 1.0 * atmospheres
    constraints = {'HBonds': app.HBonds, 'AllBonds': app.AllBonds, 'None': None}
    constraint = constraints[params['constraints']]
    barostatInterval = 25

    # Force field kwargs
    ff_kwargs = {
        'constraints': constraint,
        'rigidWater': True,    # Allows time step up to 4 fs with HMR
        'removeCMMotion': False
    }

    platform = define_platform()
    save_yaml(params, args.params)

# ---------------------------
# Energy minimization
# ---------------------------
def energy_minimisation(simulation):
    """
    Minimize the system to relieve bad contacts.
    """
    energy_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=1000)
    energy_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print('Energy difference during minimization:', energy_before - energy_after)

# ---------------------------
# Create protein system
# ---------------------------
def create_model_ppi(modeller, salt_concentration, params):
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield)

    print('Adding solvent...')
    modeller.addSolvent(forcefield,
                        boxShape='cube',
                        ionicStrength=salt_concentration * molar,
                        positiveIon='Na+',
                        negativeIon='Cl-',
                        model='tip3p',
                        neutralize=True,
                        padding=1 * nanometers)

    print('Creating force field system...')
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     constraints=constraint,
                                     rigidWater=params['rigidWater'],
                                     ewaldErrorTolerance=ewaldErrorTolerance)
    return system

# ---------------------------
# Metadynamics setup
# ---------------------------
def compute_metadynamics(metadynamics_params, system):
    # Example: add a distance-based collective variable
    atomId1 = metadynamics_params[0]['d1'][0]['atomId1'] + 1
    atomId2 = metadynamics_params[0]['d1'][1]['atomId2'] + 1

    hills_path = os.path.abspath(args.metadynamics)
    script = f"""
            d1: DISTANCE ATOMS={atomId1},{atomId2}
            METAD ARG=d1 SIGMA=0.1 HEIGHT=0.3 PACE=50 FILE={hills_path}
            PRINT ARG=d1 STRIDE=50 FILE=COLVAR
            """
    plumed = PlumedForce(script)
    plumed.setTemperature(temperature*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system

def simulate(args, params, salt_concentration=0.15):
    set_parameters(params)

    # Load initial structure
    protein = app.PDBFile(args.pdb)
    modeller = app.Modeller(protein.topology, protein.positions)

    # Create solvated protein system
    system = create_model_ppi(modeller, salt_concentration, params)

    # MetaDynamics (optional)
    # TODO move
    if params['metadynamics'] is not None:
        system = compute_metadynamics(params['metadynamics'], system)
    else:
        with open(args.metadynamics, 'w') as f:
            pass  # create dummy file

    # Precision and determinism
    properties = {
        "Precision": "mixed",
        "DeterministicForces": "true"
    }

    # ---------------------
    # Stage 0: Energy minimization
    # ---------------------
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(args.seed)
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print('STAGE 0: Running energy minimization...')
    energy_minimisation(simulation)


    # ---------------------
    # Stage 1: NVT heating with restraints
    # ---------------------
    print('STAGE 1: NVT heating with heavy atom restraints...')
    system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=10.0)

    integrator = LangevinMiddleIntegrator(100*kelvin, friction, dt)
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    for T in [100, 150, 200, 250, 300, temperature]:  # temperature ramp
        integrator.setTemperature(T*kelvin)
        simulation.step(5000)  # ~10 ps per increment
       

    # ---------------------
    # Stage 2: NPT equilibration with tapering restraints
    # ---------------------
    print('STAGE 2: Switching to NPT for density equilibration...')
    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
    simulation.context.reinitialize(preserveState=True)

    # TODO: I have no clue what is happening here
    for k in [5.0, 1.0]:
        print(f"Tapering restraints to {k} kcal/mol/Å²")
        for i in range(restraint_force.getNumParticles()):
            pos = restraint_force.getParticleParameters(i)
            particle_index = pos[0]
            (_, x0, y0, z0) = pos[1]
            restraint_force.setParticleParameters(i, particle_index, [k, x0, y0, z0])
        restraint_force.updateParametersInContext(simulation.context)
        simulation.step(50000)  # ~100 ps at each stage


    print("CHECK for forces before")
    forces = system.getForces()
    for i, force in enumerate(forces):
        print(f"Force {i}: {force}")


    # ---------------------
    # Stage 3: Unrestrained NPT equilibration
    # ---------------------
    print('STAGE 3: Removing restraints...')
    system.removeForce(system.getNumForces() - 1)  # assumes restraint is last added
    simulation.context.reinitialize(preserveState=True)
    simulation.step(100000)  # 200 ps unrestrained NPT

    print("CHECK for forces after")
    forces = system.getForces()
    for i, force in enumerate(forces):
        print(f"Force {i}: {force}")

    # ---------------------
    # Stage 4: Production run
    # ---------------------
    print(f"STAGE 4: Starting production run for {args.time} ns...")

    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj, recordInterval)
    dataReporter = app.StateDataReporter(
        args.stats, recordInterval, totalSteps=args.steps,progress=True,
        step=True, time=True, speed=True, elapsedTime=True,
        remainingTime=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True,
        separator='\t'
    )
    simulation.reporters.append(HDF5Reporter)
    simulation.reporters.append(dataReporter)

    simulation.currentStep = 0
    simulation.step(args.steps)

    # Save final frame
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(args.topo, mode="w") as file:
        app.PDBxFile.writeFile(simulation.topology, state.getPositions(), file, keepIds=True)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    parser.add_argument('--pdb', default='input/fix1.pdb', help='Path to protein PDB file.')
    parser.add_argument('--sdf', nargs='?', const='', help='Optional small molecule SDF file.')
    parser.add_argument('--md_settings', default='input/params.yml', help='MD configuration YAML.')
    parser.add_argument('--seed', type=int, default=12, help='Random seed for initial velocities.')
    parser.add_argument('--topo', default="output/top.cif", help='Output CIF of last frame.')
    parser.add_argument('--traj', default="output/traj.h5", help='Output trajectory file.')
    parser.add_argument('--stats', default="output/stats.txt", help='Output energy/statistics file.')
    parser.add_argument('--params', default="output/params.txt", help='Copy of used parameters.')
    parser.add_argument('--metadynamics', default="output/metadynamics.txt", help='Metadynamics output file.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    yaml_params = import_yaml(args.md_settings)
    simulate(args, yaml_params)
