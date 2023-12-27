#!/usr/bin/env python
import argparse
from openmm import *
from openmm.app import *
from openmm.unit import *
#from Helper import import_yaml, save_yaml
import yaml

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


def save_yaml(d, filepath):
    with open(filepath, 'w') as file:
        documents = yaml.dump(d, file)



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

    # Input Files


    # Define forcefield
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # System Configuration
    nonbondedCutoff = params['nonbondedCutoff'] * nanometers
    ewaldErrorTolerance = params['ewaldErrorTolerance']
    constraintTolerance = 0.000001
    hydrogenMass = 1.5 * amu
    constraints = {'HBonds': HBonds, 'AllBonds': AllBonds, 'None': None}

    dt = params['dt'] * picoseconds
    temperature = 310 * kelvin
    friction = 1.0 / picosecond
    pressure = 1.0 * atmospheres # TODO
    barostatInterval = 25  # TODO
    steps = int(args.steps)
    equilibrationSteps = params['equilibrationSteps']

    recordInterval = int(steps * params['recordingInterval'] / (params['time'] * 1000))

    platform = Platform.getPlatformByName('CUDA')

    # Define integrator
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    integrator.setRandomNumberSeed(int(args.seed))

    # Init MD model
    pdb = PDBFile(args.tleap)
    modeller = Modeller(pdb.topology, pdb.positions)

    print('Adding hydrogens..')
    modeller.addHydrogens(forcefield)

    print('Adding solvent..')
    modeller.addSolvent(forcefield,
                        ionicStrength=0.15 * molar,
                        model='tip3p',
                        padding=1 * nanometer)

    print('Create Forcefield..')
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=PME,
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


    # Set up log file and trajectory dcd
    dcdReporter = DCDReporter(args.traj,
                              recordInterval)

    #import mdtraj.reporters

    #h5Reporter = mdtraj.reporters.HDF5Reporter('h5md.h5', recordInterval)

    dataReporter = StateDataReporter(args.stats,
                                     recordInterval,
                                     totalSteps=steps,
                                     step=True,
                                     time=True,
                                     speed=True,
                                     progress=True, elapsedTime=True,
                                     remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                     temperature=True, volume=True, density=True, separator='\t')

    # Minimize and Equilibrate
    print('Performing energy minimization..')
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()

    print('Equilibrating..')
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(equilibrationSteps)

    print('Simulating..')
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    #simulation.reporters.append(h5Reporter)
    simulation.currentStep = 0
    simulation.step(steps)

    print('Save final frame')
    # 1. Standard export
    state = simulation.context.getState(getPositions=True,
                                        enforcePeriodicBox=system.usesPeriodicBoundaryConditions())

    # 2. Try
    with open('end1.pdb', mode="w") as file:
        PDBFile.writeFile(simulation.topology,
                          state.getPositions(),
                          file)

    with open(args.pdb_last, mode="w") as file:
        PDBFile.writeFile(simulation.topology,
                          state.getPositions(),
                          file,
                          keepIds=False)

    with open('end1.cif', mode="w") as file:
        PDBxFile.writeFile(simulation.topology,
                          state.getPositions(),
                          file,
                          keepIds=False)

    with open('end2.cif', mode="w") as file:
        PDBxFile.writeModel(simulation.topology,
                          state.getPositions(),
                          file,
                          keepIds=False)

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())
    with open(args.pdbx, mode="w") as file:
        PDBxFile.writeFile(simulation.topology,
                           state.getPositions(),
                           file,
                           keepIds=False)

    with open("final_state.cif", mode="w") as file:
        PDBxFile.writeFile(simulation.topology, state.getPositions(), file)

    with open(args.last, mode="w") as file:
        PDBxFile.writeFile(simulation.topology,
                           state.getPositions(),
                           file,
                           keepIds=False)

    import MDAnalysis as mda
    # Assuming 'simulation' is your OpenMM Simulation object
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)

    import numpy as np

    """

    # https://docs.mdanalysis.org/1.1.0/documentation_pages/core/topology.html
    # Create an empty MDAnalysis topology
    n_atoms = simulation.topology.getNumAtoms()
    n_residues = simulation.topology.getNumResidues()
    n_chains = simulation.topology.getNumChains()
    mda_topology = mda.core.topology.Topology(n_atoms, n_residues, n_chains,
                                              attrs=atom_attrs,
                                              residue_attrs=residue_attrs,
                                              segment_attrs=segment_attrs)

    # Create an MDAnalysis Universe using the topology and positions
    u = mda.Universe(mda_topology, positions)
    u.atoms.write('output.gro')  #
    """

    simulation.saveState(args.xml)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--tleap', required=False, help='Position file (inpcrd)', default='output/demo/C1s_BD001/WT/amber/C1s_BD001.tleap.pdb')

    # Parameters
    parser.add_argument('--config', required=False, help='Configuration file with all required parameters (params.yml')
    parser.add_argument('--seed', required=False, help='Seed for inital velocities', default=23)
    parser.add_argument('--threads', required=False, help='Threads', default=4)
    parser.add_argument('--simulation_name', required=False, help='Simulation name', default='test')
    parser.add_argument('--job_id', required=False, help='Output directory', default='demo')

    # Output
    parser.add_argument('--traj', required=False, help='Trajectory file')
    parser.add_argument('--stats', required=False, help='Stats File', default='test')
    parser.add_argument('--last', required=False, help='Stats File', default='last.pdb')

    parser.add_argument('--topo', required=False, help='Topology file (prmtop)', default='test')
    parser.add_argument('--crd', required=False, help='Position file (inpcrd)', default='test')
    parser.add_argument('--amber', required=False, help='Position file (inpcrd)', default='test')

    parser.add_argument('--xml', required=False, help='Position file (inpcrd)', default='test')
    parser.add_argument('--pdb_last', required=False, help='Position file (inpcrd)', default='test')
    parser.add_argument('--pdbx', required=False, help='Position file (inpcrd)', default='test')

    parser.add_argument('--directory', required=False, help='Output directory', default='test')

    args = parser.parse_args()

    # Import standard MD parameters
    params = import_yaml(args.config)

    args.steps = int(params['time'] * 1000 / params['dt'])
    args.recorded_steps = int(params['time'] * 1000 / params['recordingInterval'])

    # Save parameters to simulation folder
    save_yaml(args, os.path.join(args.directory, 'params.yml'))

    # Parse Arguments
    simulate(args, params)



# Deprecated

"""
    parser.add_argument('--topo', required=False, help='Topology file (prmtop)', default='test')
    parser.add_argument('--crd', required=False, help='Position file (inpcrd)', default='test')
    parser.add_argument('--amber', required=False, help='Position file (inpcrd)', default='test')

    parser.add_argument('--directory', required=False, help='Output directory', default='test')
"""
