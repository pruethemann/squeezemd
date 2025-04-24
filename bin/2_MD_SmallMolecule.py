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

TODO
- Catch this error for wrong cuda version: Error setting up simulation: Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)
"""

import argparse
from openmm import *
from openmm.app import *
from openmm.unit import *
#from Helper import import_yaml, save_yaml
import warnings
from openmm import app
from openmm.amd import AMDIntegrator

from openff.toolkit.topology import Molecule
from simtk import unit
from openmm import app, Platform, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller
from openmmforcefields.generators import SystemGenerator

warnings.filterwarnings('ignore')
import mdtraj
import mdtraj.reporters

import yaml
def save_yaml(d, filepath):
    """
    Save a yml file
    :param d:
    :param filepath:
    :return:
    """
    with open(filepath, 'w') as file:
        documents = yaml.dump(d, file)


def import_yaml(yaml_path: os.path):
    """
    Opens yaml file containing hyper parameters.

    :param yaml_path: File path to yaml
    :return: dictionary with parameters
    """
    try:
        with open(yaml_path, 'r') as stream:
            return yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

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

    global nonbondedCutoff, ewaldErrorTolerance, constraintTolerance, hydrogenMass, temperature
    global dt, recordInterval, friction, pressure, constraint, barostatInterval, platform
    
    # Physical parameters
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = 0.00001
    hydrogenMass = 1.5 * amu            # check optimal weight
    temperature = 310 * kelvin                    # TODO get from param file Simulation temperature

    # Time parameter
    args.steps = int(params['time'] * 1000 / params['dt'])
    args.time = params['time']
    args.recorded_steps = int(params['time'] * 1000 / params['recordingInterval'])
    dt = params['dt'] * picoseconds     # Simulation time steps
    # TODO deleteargs.equilibrationSteps = params['equilibrationSteps']
    recordInterval  = args.steps * params['recordingInterval'] // (params['time'] * 1000)

    # Constraints
    friction = 1.0 / picosecond
    pressure = 1.0 * atmospheres        # Simulation pressure
    constraints = {'HBonds': HBonds, 'AllBonds': AllBonds, 'None': None}
    constraint = constraints[params['constraints']]
    barostatInterval = 25               # Fix Barostat every 25 simulations steps

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
    set_parameters(params)

    # Define integrator for standard MD
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(args.seed)

    # Import protein pdb file. Prepared and checked for amber import
    protein = app.PDBFile(args.pdb)

    # Import small molecule
    ligand = Molecule.from_file(args.sdf)
    #ligand.assign_partial_charges('gasteiger')   

    ligand_topology = ligand.to_topology().to_openmm()
    ligand_positions = ligand.conformers[0].to_openmm()

    ff_kwargs = {
        'constraints': app.HBonds,
        'rigidWater': True,
        'removeCMMotion': False,
        'hydrogenMass': 4 * unit.amu
    }

    # 3. Use SystemGenerator to combine force fields
    system_generator = SystemGenerator(
        forcefields=['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml'],
        small_molecule_forcefield='openff-2.0.0',
        molecules=[ligand],
        cache=None,
        forcefield_kwargs=ff_kwargs
    )


    # 4. Combine protein and ligand
    modeller = Modeller(protein.topology, protein.positions)

    modeller.add(ligand_topology, ligand_positions)

    #print('Adding hydrogens..')
    #modeller.addHydrogens(system_generator.forcefield)

    print('Adding solvent..')
    modeller.addSolvent(system_generator.forcefield,
                        boxShape='cube', # 'dodecahedron'
                        ionicStrength=salt_concentration * molar,
                        positiveIon = 'Na+',
                        negativeIon = 'Cl-',
                        model='tip3p',
                        neutralize=True,
                        padding=1 * nanometer
                        )

    # 5. Create system using the generator
    system = system_generator.create_system(modeller.topology)

    # Set up integrator and simulation
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(args.seed)

    #print('Add MonteCarloBarostat')
    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

    simulation = Simulation(modeller.topology,
                            system,
                            integrator,
                            platform
                            )

    # Set coordinates
    simulation.context.setPositions(modeller.positions)

    print('Performing energy minimization..')
    energy_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()  
    simulation.minimizeEnergy()
    energy_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()  
    print('Energy before minimization: ', energy_before, energy_after)

    #######
    print('Checking energy components before minimization:')
    for i in range(system.getNumForces()):
        energy = simulation.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
        print(f'Force {i}: {system.getForce(i).__class__.__name__} = {energy}')


    print('Equilibrating..')
    simulation.context.setVelocitiesToTemperature(temperature, args.seed)
    simulation.step(params['equilibrationSteps'])

    # Set up log file and trajectory dcd
    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj, recordInterval)

    dataReporter = StateDataReporter(args.stats,
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
        PDBxFile.writeFile(simulation.topology,
                           state.getPositions(),
                           file,
                           keepIds=True)
        


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run Molecular Dynamics simulations.')
    
    # Input files
    parser.add_argument('--pdb', required=False, help='Path to single protein or protein protein complex.', default='input/fix1.pdb')
    parser.add_argument('--sdf', required=False, help='Path to small molecule sdf file', default='input/1.sdf')
    parser.add_argument('--md_settings', required=False, help='Configuration file with all required parameters (params.yml', default='input/params.yml')
    parser.add_argument('--seed', required=False, help='Seed for inital velocities', type=int, default=12)
    
    # Output
    parser.add_argument('--topo', required=False, help='Cif file of last frame', default="output/top.cif")
    parser.add_argument('--traj', required=False, help='Trajectory file', default="output/traj.h5")
    parser.add_argument('--traj_dcd', required=False, help='Trajectory file', default="output/traj.dcd")
    parser.add_argument('--stats', required=False, help='Energy saves for every 1000 frames', default="output/stats.txt")
    parser.add_argument('--params', required=False, help='MD parameter file saved for every MD', default="output/params.txt")
    return parser.parse_args()


if __name__ == '__main__':

    # Import Argparse and MDsettings from yaml
    args = parse_arguments()
    yaml_params = import_yaml(args.md_settings)

    # Set up the simulation
    simulate(args, yaml_params)

