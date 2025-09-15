#!/usr/bin/env python
"""
    Module performs Molecular Dynamics using OpenMM
    - Import of amber prepared pdb structure
    - Adding water box with 150 mM NaCl
    - Performs MD

NOTES:
    A seed is set for the integrator and the initial velocities

Export of pmrtop:
    https://github.com/openforcefield/smarty/pull/187#issuecomment-262381974
"""

import argparse, os
from openmm.unit import nanometers, kelvin,femtoseconds, picoseconds, atmospheres, molar
from openmm import app, OpenMMException, Platform, LangevinMiddleIntegrator, MonteCarloBarostat
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import mdtraj
import mdtraj.reporters
from openmmplumed import PlumedForce
from Helper import import_yaml, save_yaml


def define_platform():
    """
    Functions tries to detect if nvidia gpu and driver is available. 
    Otherwise falls back to CPU
    """
    try:
        return Platform.getPlatformByName('CUDA')
    except OpenMMException:
        # TODO: Write a log
        print("ATTENTION: no CUDA driver or GPU detected. Simulation runs on CPU")
        return Platform.getPlatformByName('CPU')

def set_parameters(params):

    global nonbondedCutoff, ewaldErrorTolerance, constraintTolerance, temperature
    global dt, recordInterval, friction, pressure, constraint, barostatInterval, platform
    global ff_kwargs

    # Physical parameters
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = 0.00001
    temperature = 310 * kelvin                    # TODO get from param file Simulation temperature

    # Time parameter
    args.steps = int(params['time'] * 1e6 / params['dt'])
    args.time = params['time']
    args.recorded_steps = int(params['time'] * 1000 / params['recordingInterval'])
    dt = params['dt'] * femtoseconds     # Simulation time steps
    # TODO deleteargs.equilibrationSteps = params['equilibrationSteps']
    recordInterval  = args.steps * params['recordingInterval'] // (params['time'] * 1000)

    # Constraints
    friction = 1.0 / picoseconds
    pressure = 1.0 * atmospheres        # Simulation pressure
    constraints = {'HBonds': app.HBonds, 'AllBonds': app.AllBonds, 'None': None}
    constraint = constraints[params['constraints']]
    barostatInterval = 25               # Fix Barostat every 25 simulations steps

    # Force field parameters
    ff_kwargs = {
        'constraints': constraint,
        'rigidWater': True,                     # Allows to increase step size to 4 fs
        'removeCMMotion': False                # System should not drift
    }

    platform = define_platform()

    # Save parameters to simulation folderTODO: combine args and md_settings
    #save_yaml(args, args.params)
    save_yaml(params, args.params)

def energy_minimisation(simulation):
    # Minimize and Equilibrate
    
    energy_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()  

    simulation.minimizeEnergy()

    energy_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()  

    print('Energy before minimization: ', energy_before, energy_after)


def create_model_ppi(modeller, salt_concentration, params):
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

    print('Adding hydrogens..')
    modeller.addHydrogens(forcefield)

    print('Adding solvent..')
    modeller.addSolvent(forcefield,
                        boxShape='cube', # 'dodecahedron'
                        ionicStrength=salt_concentration * molar,
                        positiveIon = 'Na+',
                        negativeIon = 'Cl-',
                        model='tip3p',
                        neutralize=True,
                        padding=1 * nanometers
                        )

    print('Create Forcefield..')
    system = forcefield.createSystem(modeller.topology,
                                        nonbondedMethod=app.PME,
                                        nonbondedCutoff=nonbondedCutoff,
                                        constraints=constraint,
                                        rigidWater=params['rigidWater'],
                                        ewaldErrorTolerance=ewaldErrorTolerance
                                        )
    
    return system
    


def create_model_smallmolecule(modeller, salt_concentration):
    ligand = Molecule.from_file(args.sdf)
    #ligand.assign_partial_charges('gasteiger')   

    ligand_topology = ligand.to_topology().to_openmm()
    ligand_positions = ligand.conformers[0].to_openmm()

    # 3. Use SystemGenerator to combine force fields
    system_generator = SystemGenerator(
        forcefields=['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml'],
        small_molecule_forcefield='openff-2.0.0',
        molecules=[ligand],
        cache=None,
        forcefield_kwargs=ff_kwargs
    )

    # Add ligand to modell
    modeller.add(ligand_topology, ligand_positions)

    # Adding hydrogens, fixes small issues like protonation states 
    # at N/C terminus but may change HiS protonation
    modeller.addHydrogens(system_generator.forcefield)

    print(f'Adding solvent and {salt_concentration} M NaCl ..')
    modeller.addSolvent(system_generator.forcefield,
                        boxShape='cube', # 'dodecahedron'
                        ionicStrength=salt_concentration * molar,
                        positiveIon = 'Na+',
                        negativeIon = 'Cl-',
                        model='tip3p',
                        neutralize=True,
                        padding=1 * nanometers
                        )
    

    # 5. Create system using the generator
    system = system_generator.create_system(modeller.topology)

    return system

def compute_metadynamics(metadynamics_params, system):
    # Generate collective variables
    # here we have 2 distances together using the PLUMED syntax
    # Remember: Add +1 to each atom index when adding to PLUMED
    # C-terminus GLu-92: 863
    # Receptor: Arg-522:2177

    # Distance 1
    print(metadynamics_params[0]['d1'][1]['atomId2'])
    atomId1 = metadynamics_params[0]['d1'][0]['atomId1'] + 1
    atomId2 = metadynamics_params[0]['d1'][1]['atomId2'] + 1

    # distance 2
    #atomId3 = metadynamics_params[1]['d2'][0]['atomId1'] + 1
    #atomId4 = metadynamics_params[1]['d2'][1]['atomId2'] + 1
    
    hills_path = os.path.abspath(args.metadynamics)
    print(args.metadynamics)
    print(hills_path)

    script = f"""
            d1: DISTANCE ATOMS={atomId1},{atomId2}
            METAD ARG=d1 SIGMA=0.1 HEIGHT=0.3 PACE=50 FILE={hills_path}
            PRINT ARG=d1 STRIDE=50 FILE=COLVAR
            """
    
    """ For two distances
            d1: DISTANCE ATOMS={atomId1},{atomId2}
            d2: DISTANCE ATOMS={atomId3},{atomId4}
            METAD ARG=d1,d2 SIGMA=0.1,0.1 HEIGHT=0.3 PACE=50 FILE={args.metadynamics}
    """
    plumed = PlumedForce(script)
    plumed.setTemperature(310*kelvin)  # sets kBT internally
    system.addForce(plumed)
    print("metdadynics variable added")
    return system

def simulate(args, params, salt_concentration=0.15):
    """
    Function that handles a molecular dynamics simulation.
    Most MD parameters are saved in a job specific params.yml

    Barostat:   Monte Carlo Barostat
    Integrator: Langevin Middle Integrator


    Examples:
        https://notebooks.githubusercontent.com/view/ipynb?browser=chrome&color_mode=auto&commit=1bc6a022bb6c07b1389a2f18749d8fcc01304ea3&device=unknown&enc_url=68747470733a2f2f7261772e67697468756275736572636f6e74656e742e636f6d2f63686f646572616c61622f6f70656e6d6d2d7475746f7269616c732f316263366130323262623663303762313338396132663138373439643866636330313330346561332f30322532302d253230496e7465677261746f7273253230616e6425323073616d706c696e672e6970796e62&logged_in=false&nwo=choderalab%2Fopenmm-tutorials&path=02+-+Integrators+and+sampling.ipynb&platform=android&repository_id=100135600&repository_type=Repository&version=99
    Further information:
        http://docs.openmm.org/latest/userguide/application/02_running_sims.html

    Attributes:
        dt (Quantity): The time step for the simulation.
        temperature (Quantity): The temperature of the simulation.
        friction (Quantity): The friction coefficient for the simulation.
        pressure (Quantity): The pressure of the simulation.
        barostatInterval (int): The interval at which to apply the barostat.
        equilibrationSteps (int): The number of equilibration steps to run.
        recordInterval (int): The number of steps between saving frames.
        simulation (Simulation): The OpenMM Simulation object.
        traj (str): The name of the trajectory file.
        dataReporter (StateDataReporter): The StateDataReporter object for recording statistics.
        forcefield (ForceField): The OpenMM ForceField object.
        modeller (Modeller): The OpenMM Modeller object.
        pdb (PDBFile): The PDB file object.

    """

    # Add protein (receptor) to model
    set_parameters(params)

    # Define integrator for standard MD
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(args.seed)

    # Import protein pdb file. Prepared and checked for amber import
    protein = app.PDBFile(args.pdb)

    # Add protein (receptor) to model
    modeller = app.Modeller(protein.topology, protein.positions)

    # add small molecule to system


    system = create_model_ppi(modeller, salt_concentration, params)


    print('Add MonteCarloBarostat')
    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

    # MetaDynamics option
    if params['metadynamics'] != None:
        print(params['metadynamics'])
        system = compute_metadynamics(params['metadynamics'], system)
    else:
        # create dummy file to satisfy snakemake
        with open(args.metadynamics, 'w') as f:
            pass  # Ensure no content is written


    simulation = app.Simulation(modeller.topology,
                            system,
                            integrator,
                            platform
                            )

    # Set coordinates
    simulation.context.setPositions(modeller.positions)

    # Minimise the energy of the system
    energy_minimisation(simulation)

    # TODO probably not necessary
    print('Checking energy components before minimization:')
    """
    for i in range(system.getNumForces()):
        energy = simulation.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
        print(f'Force {i}: {system.getForce(i).__class__.__name__} = {energy}')
    """

    print('Equilibrating..')
    simulation.context.setVelocitiesToTemperature(temperature, args.seed)
    simulation.step(params['equilibrationSteps'])

    # Set up log file and trajectory dcd
    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj, recordInterval)

    dataReporter = app.StateDataReporter(args.stats,
                                        recordInterval,
                                        totalSteps=args.steps,
                                        step=True,
                                        time=True,
                                        speed=True,
                                        progress=True, elapsedTime=True,
                                        remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                        temperature=True, volume=True, density=True, separator='\t')

    simulation.reporters.append(HDF5Reporter)
    simulation.reporters.append(dataReporter)

    print(f"Running simulation for {args.time} ns.")
    simulation.currentStep = 0
    simulation.step(args.steps)

    # Save final frame as topology.cif
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)

    with open(args.topo, mode="w") as file:
        app.PDBxFile.writeFile(simulation.topology,
                           state.getPositions(),
                           file,
                           keepIds=True)
        

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    
    # Input files
    parser.add_argument('--pdb', required=False, help='Paenergy_afterth to single protein or protein protein complex.', default='input/fix1.pdb')
    parser.add_argument('--sdf', required=False, help='Path to small molecule sdf file', nargs='?', const='')
    parser.add_argument('--md_settings', required=False, help='Configuration file with all required parameters (params.yml', default='input/params.yml')
    parser.add_argument('--seed', required=False, help='Seed for inital velocities', type=int, default=12)
    
    # Output
    parser.add_argument('--topo', required=False, help='Cif file of last frame', default="output/top.cif")
    parser.add_argument('--traj', required=False, help='Trajectory file', default="output/traj.h5")
    parser.add_argument('--stats', required=False, help='Energy saves for every 1000 frames', default="output/stats.txt")
    parser.add_argument('--params', required=False, help='MD parameter file saved for every MD', default="output/params.txt")
    parser.add_argument('--metadynamics', required=False, help='MD parameter file saved for every MD', default="output/metadynamics.txt")
    return parser.parse_args()


if __name__ == '__main__':

    # Import Argparse and MDsettings from yaml
    args = parse_arguments()
    yaml_params = import_yaml(args.md_settings)

    simulate(args, yaml_params)