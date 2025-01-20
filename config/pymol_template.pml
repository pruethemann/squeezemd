# This pyMol command file is used as a template to visualize the binding surface


# Load the pdb file
load {input_pdb};

# Clean up the strucutre
remove solvent;
show cartoon;
color grey;

# Show sticks of interacting reisudes: ligand / Receptor
show sticks, (resid {ligand_resids}) and chain A;
show sticks, (resid {receptor_resids}) and not chain A;

# The interaction "entalphies" are encoded in the b factors of the pdb file
# color coding for interactions for ligand
spectrum b, red blue grey;
spectrum b, red blue grey, chain A;

# Color coding for interactions for receptor
spectrum b, red blue grey, not chain A;

# ChainID
show surface, not chain A;
set transparency, 0.2, not chain A;
extract ligand, chain A;   # To avoid pymol bug of broken surface

# Save all the output as pymol file and png
save {output};
png {output_png};

dele all;
