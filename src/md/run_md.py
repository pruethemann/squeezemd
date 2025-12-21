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
from openmm import app, OpenMMException, Platform, LangevinMiddleIntegrator, MonteCarloBarostat, CustomExternalForce
from openmm.unit import kilojoule_per_mole, nanometers, femtoseconds, kelvin, molar, picoseconds, atmospheres
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import mdtraj
from Helper import import_yaml
from metadynamics_auxillary import add_metadynamics_forces_centerofmass, save_active_forces

def add_positional_restraints(system, topology, positions, k=10.0, flexible_resids={}):
    """
    Add harmonic restraints to heavy atoms (kcal/mol/Å²).
    Applied to all non-solvent heavy atoms.
    """
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    
    system.addForce(restraint)

    restraint.addGlobalParameter('k', k*kilojoule_per_mole/nanometers**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    if len(flexible_resids) > 0:
        # Do not restrain water, ions and ligands during equilibration
        unrestrained_residues = ('HOH', 'Na+', 'Cl-', 'CL', 'UNK')
    else: # keep flexible binding pocket -> ligand (UNK) can move
        unrestrained_residues = ('HOH', 'Na+', 'Cl-', 'CL')

    for atom in topology.atoms():
        res = atom.residue      # information about residue
        resname = res.name      # Either amino acid name or ion name or ligand name
        resid = int(res.id)     # PDB residue number (string!)
        resindex = res.index    # 0-based OpenMM index (int)

        # Exclude water and ions and ligand from the restraints
        if resname in unrestrained_residues or atom.element.symbol == 'H':
            continue

        # Exclude flexible amino acids from restraining
        if resid in flexible_resids:
            print(resid, flexible_resids[resid])
            continue

        # Restrain the rest
        if args.verbose:
            print(f"Residue {resname} {resid} (index {resindex}), "f"Atom {atom.name}, element {atom.element.symbol}, "f"atom index {atom.index}")
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
    # TODO remove steps
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
    
    # Assign partial charges
    ligand.assign_partial_charges('gasteiger')   

    ligand_topology = ligand.to_topology().to_openmm()
    ligand_positions = ligand.conformers[0].to_openmm()

    ff_kwargs = {
        'constraints':app.HBonds,
        'rigidWater': True,# TODO standardize with yaml
        'ewaldErrorTolerance':params['simulation']['constraints']['ewald_error_tolerance']
    }
    periodic_forcefield_kwargs = {
        'nonbondedMethod':app.PME,
        'nonbondedCutoff':params['simulation']['constraints']['cutoff_nm'] * nanometers
    }

    # 3. Use SystemGenerator to combine force fields
    generator = SystemGenerator(
        forcefields=[protein_forcefield, water_model],
        small_molecule_forcefield="openff-2.2.0",           # TODO: make sure to update to 3.0 if released soon
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
    if args.mode == 'protein_molecule':
        system = create_model_smallmolecule(modeller, salt_concentration, params, args.sdf)
    else: # PPi, metadynamics, protein
        system = create_model_ppi(modeller, salt_concentration, params)

    # ---------------------
    # Stage 0: Minimization
    # ---------------------
    print('\n=== Stage 0: Energy minimization with restraints ===')

    # Add restraints for equilibration BEFORE minimization
    k = params['simulation']['equilibration']['protein_k']
    system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=k)

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
  
    for k in [k/2, k/10]: # usually k = 10 -> 5 -> 1
        print(f"Tapering restraints to {k} kcal/mol/Å²")
        # Update global k parameter (not per particle!)
        simulation.context.setParameter('k', k * kilojoule_per_mole / nanometers**2)
        print("k in context:", simulation.context.getParameter('k'))

        # Run equilibration for each step
        simulation.step(params['simulation']['equilibration']['NPT_equilibration'])

    # ---------------------
    # Stage 3: Unrestrained NPT at simulation T
    # ---------------------
    print('\n=== Stage 3: Unrestrained NPT equilibration ===')

    if args.verbose:
        save_active_forces(system, simulation.context, logfile='before.txt')

    print('\n=== Removing positional restraints ===')
    for i, force in enumerate(system.getForces()):
        # Remove 06 CustomExternalForce
        if force.__class__.__name__ == restraint_force.__class__.__name__:
            system.removeForce(i)
            break

    if args.verbose:
        save_active_forces(system, simulation.context, logfile='after.txt')

    simulation.context.reinitialize(preserveState=True)
    simulation.step(params['simulation']['equilibration']['NPT_unrestrained'])

    # Save equilibrated pdb
    save_pdb(simulation, args.equilibrated)

    # ---------------------
    # Stage 4: Metadynamics (optional)
    # ---------------------

    if args.mode == 'metadynamics':
        print(f'\n=== Stage 4: Initiate Metadynamics')
        simulation.system = add_metadynamics_forces_centerofmass(params, simulation.system, args, T)
        simulation.context.reinitialize(preserveState=True)  # keep positions/velocities

    # ---------------------
    # Stage 5: Rigidify receptor except flexible binding pocket (optional)
    # ---------------------

    # Define tapering schedule for restraints (kcal/mol/Å²)
    # Rigify everything except 
    if 'flexible_binding_pocket' in params['simulation']:
        flexible_resids = params['simulation']['flexible_binding_pocket']['flexible_resids']
        k = params['simulation']['flexible_binding_pocket']['protein_k']
        system, restraint_force = add_positional_restraints(system, modeller.topology, modeller.positions, k=k, flexible_resids=flexible_resids)
        
        # is this required?
        simulation.context.setParameter('k', k * kilojoule_per_mole / nanometers**2)
        print("k in context:", simulation.context.getParameter('k'))
        simulation.context.reinitialize(preserveState=True)
        save_active_forces(system, simulation.context, logfile='flexible_binding_pocket.txt')

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


def parse_arguments():
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    # Required Input
    parser.add_argument('--pdb', default='input/protein.pdb')
    parser.add_argument('--md_settings', default='input/params.yml')  
    parser.add_argument('--mode', default='PPi')

    # Optional Input
    parser.add_argument('--seed', type=int, default=12)
    parser.add_argument('--verbose', default=False)

    # Output
    parser.add_argument('--equilibrated', default='output/equilibrated.pdb', help="Equilibrated system in water box")
    parser.add_argument('--topo_cif', default='output/top.cif', help="Last uncentered frame of MD")
    parser.add_argument('--traj', default='output/traj.h5', help="h5md trajectory")
    parser.add_argument('--stats', default='output/stats.txt', help="Molecular dynamics statistics and progress")
    parser.add_argument('--metadynamics_hills', help='Metadynamics hill output file.')
    parser.add_argument('--metadynamics_colvar', default="output/metadynamics_colvar.txt", help='Metadynamics output file.')
    parser.add_argument('--sdf', required=False, default="input/ligand.sdf",help='Small molecule sdf file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    params = import_yaml(args.md_settings)

    simulate(args, params)
