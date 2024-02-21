
# Load the pdb file
load {input_pdb};
remove solvent
show cartoon;
color grey;
spectrum b, red blue grey;


# Show sticks of interacting reisudes
#color dblue, (resid {ligand_resids}) and chain I;
show sticks, (resid {ligand_resids}) and chain I;

#color red, (resid {receptor_resids}) and not chain I;
show sticks, (resid {receptor_resids}) and not chain I;

# Gradient for Inhibitor
spectrum b, red blue grey, chain I

# Gradient for Inhibitor
spectrum b, red blue grey, not chain I

# ChainID depends on remap
show surface, chain B
set transparency, 0.1, chain B
extract {target}, chain I
save {output};
