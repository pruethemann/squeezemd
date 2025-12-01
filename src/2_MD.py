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
from Helper import import_yaml, save_yaml


# ---------------------------
# Helper functions
# ---------------------------
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
    protein_forcefield = params['forcefield']['protein']
    water_model = params['forcefield']['water']

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
                        model=params['forcefield']['watermodel'],                
                        boxShape='cube',
                        ionicStrength=salt_concentration * molar,
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

    protein_forcefield = params['forcefield']['protein']
    water_model = params['forcefield']['water']

    print(f'Initializing ForceField: {protein_forcefield} + {water_model}')
    forcefield = app.ForceField(protein_forcefield, water_model)

    modeller.addHydrogens(forcefield)       # TODO: Check whether His protonation states are changed
    modeller.addExtraParticles(forcefield)          # Required for tip4p (orbital)

    # Add solvent
    modeller.addSolvent(forcefield,
                        model=params['forcefield']['watermodel'],                
                        boxShape='cube',
                        ionicStrength=salt_concentration * molar,
                        positiveIon='Na+',
                        negativeIon='Cl-',
                        neutralize=True,
                        padding=1.2 * nanometers)
    
    # Create the MD system
    system = forcefield.createSystem(modeller.topology,
                             nonbondedMethod=app.PME,
                             nonbondedCutoff=params['nonbondedCutoff'] * nanometers,
                             constraints=app.HBonds,
                             rigidWater=True,
                             ewaldErrorTolerance=params['ewaldErrorTolerance'])
    return system

def save_cif(simulation, cif_file: os.path):
    positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    with open(cif_file, "w") as f:
        app.PDBxFile.writeFile(simulation.topology, positions, f, keepIds=True)

def save_pdb(simulation, pdb_file:os.path):
    positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    with open(pdb_file, "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions,f, keepIds=True)

# ---------------------------
# Simulation procedure
# ---------------------------

def simulate(args, params, salt_concentration=0.15):
    """
    Set up and start the simulation
    """
    # Detect GPU
    platform = define_platform()

    # Load structure
    protein = app.PDBFile(args.pdb)
    modeller = app.Modeller(protein.topology, protein.positions)

    # Create solvated system
    if args.sdf == "-1":
        system = create_model_ppi(modeller, salt_concentration, params)
    else:
        system = create_model_smallmolecule(modeller, salt_concentration, params, args.sdf)
        

    # MetaDynamics (optional)
    # TODO move
    if params['metadynamics'] is not None:
        # Only import if required. Currently doesn't work with openmm 8.3.1
        from openmmplumed import PlumedForce
        system = compute_metadynamics(params['metadynamics'], system)
    else:
        with open(args.metadynamics, 'w') as f:
            pass  # create dummy file

    # Add restraints BEFORE minimization
    system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=10.0)

    # Integrator setup
    dt = params['dt'] * femtoseconds
    temperature = params['temperature'] * kelvin
    friction = 1.0 / picoseconds
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(params['constraintTolerance'])
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
    temp_steps = [50, 100, 150, 200, 250, 300, params['temperature']]
    steps_per_temp = params['NVT_heating']  # e.g. 5000 = 10 ps
    for T in temp_steps:
        print(f" → Heating to {T} K ...")
        simulation.integrator.setTemperature(T * kelvin)
        simulation.step(steps_per_temp)

    # ---------------------
    # Stage 2: NPT equilibration with tapering restraints
    # ---------------------
    print('\n=== Stage 2: NPT equilibration with tapering restraints ===')
    barostat = MonteCarloBarostat(1.0 * atmospheres, params['temperature'] * kelvin, 25)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)

    # Define tapering schedule for restraints (kcal/mol/Å²)
    for k in [5.0, 1.0]:
        print(f"Tapering restraints to {k} kcal/mol/Å²")
        # Update global k parameter (not per particle!)
        restraint_force.setGlobalParameterDefaultValue(0, k * kilojoule_per_mole / nanometer**2)
        restraint_force.updateParametersInContext(simulation.context)

        # Run equilibration for each step
        simulation.step(params['NPT_equilibration'])

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
    simulation.step(params['NPT_unrestrained'])

    # ---------------------
    # Stage 4: Production
    # ---------------------
    print(f'\n=== Stage 4: Production run ({params["time"]} ns) ===')
    recordInterval = int(params['recordingInterval'] * 1000 / params['dt'])
    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj, recordInterval)
    dataReporter = app.StateDataReporter(
        args.stats, recordInterval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True,
        volume=True, density=True, progress=True, remainingTime=True,
        speed=True, totalSteps=int(params['time'] * 1e6 / params['dt'])
    )
    simulation.reporters.append(HDF5Reporter)
    simulation.reporters.append(dataReporter)

    simulation.step(int(params['time'] * 1e6 / params['dt']))

    # Save the end file as cif
    save_cif(simulation, args.topo_cif)


# ---------------------------
# Argument parsing
# ---------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    parser.add_argument('--pdb', default='input/fix1.pdb')
    parser.add_argument('--md_settings', default='input/params.yml')
    parser.add_argument('--seed', type=int, default=12)
    parser.add_argument('--topo_cif', default='output/top.cif')
    parser.add_argument('--traj', default='output/traj.h5')
    parser.add_argument('--stats', default='output/stats.txt')
    parser.add_argument('--metadynamics', default="output/metadynamics.txt", help='Metadynamics output file.')
    parser.add_argument('--sdf', required=False,help='Small molecule sdf file', default='0')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    params = import_yaml(args.md_settings)
    simulate(args, params, salt_concentration=params['salt'])
