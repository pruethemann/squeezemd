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

import argparse
from openmm import *
from openmm.app import *
from openmm.unit import *
from Helper import import_yaml, save_yaml
import warnings
from openmm import app
# suppress some MDAnalysis warnings when writing PDB files
warnings.filterwarnings('ignore')
import mdtraj
import mdtraj.reporters
import shutil

def simulate(args, params):
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

    # System Configuration
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = 0.000001
    hydrogenMass = 1.5 * amu
    constraints = {'HBonds': HBonds, 'AllBonds': AllBonds, 'None': None}

    dt = params['dt'] * picoseconds     # Simulation time steps
    temperature = 310 * kelvin          # Simulation temperature
    friction = 1.0 / picosecond
    pressure = 1.0 * atmospheres        # Simulation pressure
    barostatInterval = 25               # Fix Barostat every 25 simulations steps
    steps = int(args.steps)
    equilibrationSteps = params['equilibrationSteps']

    recordInterval = int(steps * params['recordingInterval'] / (params['time'] * 1000))

    platform = Platform.getPlatformByName('CUDA')

    # Define integrator
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(args.seed)

    # Init MD model
    pdb = PDBFile(args.input_pdb)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    print('Adding hydrogens..')
    modeller.addHydrogens(forcefield)

    print('Adding solvent..')
    # Define amber forcefield
    modeller.addSolvent(forcefield,
                        ionicStrength=0.15 * molar,
                        model='tip3p',
                        padding=1 * nanometer)


    print('Create Forcefield..')
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     constraints=constraints[params['constraints']],
                                     rigidWater=params['rigidWater'],
                                     ewaldErrorTolerance=ewaldErrorTolerance,
                                     hydrogenMass=hydrogenMass
    )

    print('Add MonteCarloBarostat')
    system.addForce(MonteCarloBarostat(pressure,temperature,barostatInterval))

    simulation = Simulation(modeller.topology,
                            system,
                            integrator,
                            platform
                            )


    # Minimize and Equilibrate
    print('Performing energy minimization..')
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()

    print('Equilibrating..')
    simulation.context.setVelocitiesToTemperature(temperature, args.seed)
    simulation.step(equilibrationSteps)


    # Set up log file and trajectory dcd
    HDF5Reporter = mdtraj.reporters.HDF5Reporter(args.traj,recordInterval)

    dataReporter = StateDataReporter(args.stats,
                                     recordInterval,
                                     totalSteps=steps,
                                     step=True,
                                     time=True,
                                     speed=True,
                                     progress=True, elapsedTime=True,
                                     remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                     temperature=True, volume=True, density=True, separator='\t')

    simulation.reporters.append(HDF5Reporter)
    simulation.reporters.append(DCDReporter('test.dcd', recordInterval))
    simulation.reporters.append(dataReporter)



    print('Simulating..')
    simulation.currentStep = 0
    simulation.step(steps)

    # Save final frame as topology.cif
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())

    with open(args.topo, mode="w") as file:
        PDBxFile.writeFile(simulation.topology,
                           state.getPositions(),
                           file,
                           keepIds=True)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--input_pdb', required=False, help='Amber and Tleap prepated pdb file from complex or single protein')
    parser.add_argument('--md_settings', required=False, help='Configuration file with all required parameters (params.yml')
    parser.add_argument('--seed', required=False, help='Seed for inital velocities', default=23, type=int)

    # Output
    parser.add_argument('--topo', required=False, help='Cif file of last frame')
    parser.add_argument('--traj', required=False, help='Trajectory file')
    parser.add_argument('--traj_center', required=False, help='MD parameter file saved for every MD')
    parser.add_argument('--stats', required=False, help='Energy saves for every 1000 frames')
    parser.add_argument('--params', required=False, help='MD parameter file saved for every MD')

    args = parser.parse_args()

    # Import standard MD parameters
    params = import_yaml(args.md_settings)

    args.steps = int(params['time'] * 1000 / params['dt'])
    args.recorded_steps = int(params['time'] * 1000 / params['recordingInterval'])

    # Save parameters to simulation folder
    save_yaml(args, args.params)

    # Parse Arguments
    simulate(args, params)
