import MDAnalysis as mda

# Load the PDB file
u = mda.Universe('13_C1s-BD001.pdb')

# Select atoms in chain A
chain_a = u.select_atoms('chainID A')

# Save the atom ids in a list
atom_ids_chain_a = [atom.id for atom in chain_a.atoms]

# Print the list of atom ids
print(atom_ids_chain_a)
