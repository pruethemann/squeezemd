from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Load your PDB file
pdb = PDBFile('13_C1s-BD001.pdb')

# Select a force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create a modeller object
modeller = Modeller(pdb.topology, pdb.positions)

# Add hydrogens to the modeller object if not already present
modeller.addHydrogens(forcefield)

# Add a water box around the protein with padding. The padding is the distance from the protein to the box edge.
padding = 1.0 * nanometers
modeller.addSolvent(forcefield, model='tip3p', padding=padding, neutralize=True)


# Create a system from the force field and pdb topology
system = forcefield.createSystem(pdb.topology)

# Specify the target point in 3D space (x0, y0, z0) you want to pull the atom towards
x0, y0, z0 = 68,  50,  25  # Example coordinates in nanometers

# Identify the atom index for the atom you want to pull on in amino acid 122, chain A
atom_index = None
for atom in pdb.topology.atoms():
    if atom.residue.id == '122' and atom.residue.chain.id == 'A':
        atom_index = atom.index
        break

if atom_index is None:
    raise ValueError("Atom to pull on not found. Please check the residue and chain.")

print(atom_index)

# Define a custom external force that points towards the specified 3D point
force_expression = '0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)'
force = CustomExternalForce(force_expression)
force.addGlobalParameter('k', 5.0 *10**-2 * kilocalories_per_mole/angstroms**2)  # Spring constant
force.addGlobalParameter('x0', x0 * nanometers)
force.addGlobalParameter('y0', y0 * nanometers)
force.addGlobalParameter('z0', z0 * nanometers)
force.addParticle(atom_index, [])

system.addForce(force)

# Simulation setup
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Minimize energy and equilibrate
simulation.minimizeEnergy()
#simulation.step(10000)  # Number of steps for equilibration

# Add reporters for output (optional)
simulation.reporters.append(StateDataReporter('traj.csv', 1000, step=True, potentialEnergy=True, temperature=True))

# Set up log file and trajectory dcd
dcdReporter = app.dcdreporter.DCDReporter('traj.dcd', 5)
simulation.reporters.append(dcdReporter)

# Run the simulation
simulation.step(1000)  # Adjust the number of steps as needed

print("Steered MD simulation towards a 3D point completed.")
