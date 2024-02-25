
# Load the pdb file
load {input_pdb};
remove solvent
show cartoon;
color grey;
spectrum b, red blue grey;


# Show sticks of interacting reisudes
# TODO switch back to I
#color dblue, (resid {ligand_resids}) and chain B;
show sticks, (resid {ligand_resids}) and chain B;

#color red, (resid {receptor_resids}) and not chain B;
show sticks, (resid {receptor_resids}) and not chain B;

# Gradient for Inhibitor
spectrum b, red blue grey, chain B

# Gradient for receptor
spectrum b, red blue grey, not chain B

# ChainID depends on remap
show surface, chain B
set transparency, 0.1, chain B
extract ligand, chain B
save {output};
png {output_png};
dele all;
