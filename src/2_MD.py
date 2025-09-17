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
        resname = atom.residue.name
        if resname not in ('HOH', 'Na+', 'Cl-') and atom.element.symbol != 'H':
            pos = positions[atom.index]
            force.addParticle(atom.index, [k, pos.x, pos.y, pos.z])

    system.addForce(force)
    return system, force


def define_platform():
    """
    Detect NVIDIA GPU with CUDA, fallback to CPU if not available.
    """
    try:
        return Platform.getPlatformByName('CUDA')
    except OpenMMException:
        print("ATTENTION: no CUDA driver or GPU detected. Simulation runs on CPU")
        return Platform.getPlatformByName('CPU')

def set_parameters(params):
    global nonbondedCutoff, ewaldErrorTolerance, constraintTolerance, temperature
    global dt, recordInterval, friction, pressure, constraint, barostatInterval, platform
    global ff_kwargs

    # Physical parameters
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = params['constraintTolerance']
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
    # TODO currently not used. Required for small molecules!
    ff_kwargs = {
        'constraints': constraint,
        'rigidWater': True,    # Allows time step up to 4 fs with HMR
        'removeCMMotion': False
    }

    platform = define_platform()
    save_yaml(params, args.params)

def energy_minimisation(simulation):
    """
    Minimize the system to relieve bad contacts.
    """
    energy_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=1000)
    energy_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print('Energy difference during minimization:', energy_before - energy_after)

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


def save_cif(simulation, cif_path:os.path):
    # Save final frame
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(cif_path, mode="w") as file:
        app.PDBxFile.writeFile(simulation.topology, state.getPositions(), file, keepIds=True)

def compute_metadynamics(metadynamics_params, system):
    """
    in progress
    """
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

def get_force_paramters(system, stage):
    forces = system.getForces()
    print("STAGE ", stage)
    for i, f in enumerate(forces):
        print(f"Force {i}: {f.__class__.__name__}")

    print()

def get_force_paramters_extended(system, stage):
    import openmm
    print(f"CHECK for forces STAGE {stage}")
    forces = system.getForces()
    for i, force in enumerate(forces):
        print(f"\nForce {i}: {force.__class__.__name__}")

        # HarmonicBondForce
        if isinstance(force, openmm.HarmonicBondForce):
            print(f"  Number of bonds: {force.getNumBonds()}")
            for j in range(force.getNumBonds()):
                p1, p2, length, k = force.getBondParameters(j)
                print(f"    Bond {j}: particles=({p1},{p2}), length={length}, k={k}")

        # HarmonicAngleForce
        elif isinstance(force, openmm.HarmonicAngleForce):
            print(f"  Number of angles: {force.getNumAngles()}")
            for j in range(force.getNumAngles()):
                p1, p2, p3, angle, k = force.getAngleParameters(j)
                print(f"    Angle {j}: particles=({p1},{p2},{p3}), angle={angle}, k={k}")

        # PeriodicTorsionForce
        elif isinstance(force, openmm.PeriodicTorsionForce):
            print(f"  Number of torsions: {force.getNumTorsions()}")
            for j in range(force.getNumTorsions()):
                p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(j)
                print(f"    Torsion {j}: particles=({p1},{p2},{p3},{p4}), "
                    f"periodicity={periodicity}, phase={phase}, k={k}")

        # NonbondedForce
        elif isinstance(force, openmm.NonbondedForce):
            print(f"  Number of particles: {force.getNumParticles()}")
            for j in range(min(5, force.getNumParticles())):  # only show first few
                charge, sigma, epsilon = force.getParticleParameters(j)
                print(f"    Particle {j}: charge={charge}, sigma={sigma}, epsilon={epsilon}")

        # CustomExternalForce
        elif isinstance(force, openmm.CustomExternalForce):
            print(f"  Number of particles: {force.getNumParticles()}")
            for j in range(force.getNumParticles()):
                p, params = force.getParticleParameters(j)
                print(f"    Particle {p}: params={params}")

        # CustomBondForce
        elif isinstance(force, openmm.CustomBondForce):
            print(f"  Number of bonds: {force.getNumBonds()}")
            for j in range(force.getNumBonds()):
                p1, p2, params = force.getBondParameters(j)
                print(f"    Bond {j}: particles=({p1},{p2}), params={params}")

        # Catch-all
        else:
            print("  Parameters not implemented for this force type.")

def debug_traj(simulation, traj_path):
    from openmm.app import DCDReporter
    dcd = DCDReporter(traj_path, 200)
    simulation.reporters.append(dcd)

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
    debug_traj(simulation, 'minimize.dcd')
    energy_minimisation(simulation)

    get_force_paramters(system, 0)

    save_cif(simulation, 'minimize.cif')

    # ---------------------
    # Stage 1: NVT heating with restraints
    # ---------------------
    print('STAGE 1: NVT heating with heavy atom restraints...')
    system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=10.0)

    integrator = LangevinMiddleIntegrator(100*kelvin, friction, dt)
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print(simulation.context.getPlatform().getName())      # e.g. 'CUDA'

    debug_traj(simulation, 'heat.dcd')

    for T in [100, 150, 200, 250, 300]:  # temperature ramp
        integrator.setTemperature(T*kelvin)
        simulation.step(params['NVT_heating'])  # ~10 ps per increment
       
    get_force_paramters(system, 1)
    save_cif(simulation, "stage_1.cif")

    # ---------------------
    # Stage 2: NPT equilibration with tapering restraints
    # ---------------------
    print('STAGE 2: Switching to NPT for density equilibration...')
    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
    simulation.context.reinitialize(preserveState=True)

    debug_traj(simulation, 'npt.dcd')

    # Reduce protein restrain
    for k in [5.0, 1.0]:
        print(f"Tapering restraints to {k} kcal/mol/Å²")
        for i in range(restraint_force.getNumParticles()):
            (particle_index, parameters) = restraint_force.getParticleParameters(i)
            (_, x0, y0, z0) = parameters
            restraint_force.setParticleParameters(i, particle_index, [k, x0, y0, z0])
        restraint_force.updateParametersInContext(simulation.context)
        simulation.step(params['NPT_equilibration'])  # ~100 ps at each stage

    get_force_paramters(system, 2)
    save_cif(simulation, "stage_2.cif")

    # ---------------------
    # Stage 3: Unrestrained NPT equilibration
    # ---------------------
    print('STAGE 3: Removing restraints...')
   
    # remove the custom force
    for i, force in enumerate(system.getForces()):
        if force.__class__.__name__ == restraint_force.__class__.__name__:
            system.removeForce(i)
            break

    # system.removeForce(system.getNumForces() - 1)  # only removes barostat assumes restraint is last added
    simulation.context.reinitialize(preserveState=True)
    simulation.step(params['NPT_unrestrained'])  # 200 ps unrestrained NPT

    get_force_paramters(system, 3)

    save_cif(simulation, "stage_3.cif")

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
    save_cif(simulation, args.topo)

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
